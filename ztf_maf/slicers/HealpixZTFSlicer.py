import numpy as np
from functools import wraps
import matplotlib.path as mplPath

import lsst.sims.maf.slicers as slicers
from lsst.sims.maf.plots import HealpixSkyMap, HealpixHistogram
from lsst.sims.maf.utils.mafUtils import gnomonic_project_toxy


# create a class for pixelizing the data
class HealpixZTFSlicer(slicers.HealpixSDSSSlicer):
    """For use with ZTF simulated """

    def __init__(self, nside=128, lonCol='FieldRA', latCol='FieldDec', verbose=True,
                 useCache=True, radius=3.5 * np.sqrt(2), leafsize=100, **kwargs):
        """Using one corner of the chip as the spatial key and the diagonal as the radius.  """
        super(HealpixZTFSlicer, self).__init__(verbose=verbose,
                                               lonCol=lonCol, latCol=latCol,
                                               radius=radius, leafsize=leafsize,
                                               useCache=useCache, nside=nside)
        self.plotFuncs = [HealpixSkyMap, ]
        # TODO: HealpixHistogram raises ValueError: HealpixHistogram is for use with healpix slicer.
        #self.plotFuncs = [HealpixSkyMap, HealpixHistogram]
        # radius needs to be chip size

    def setupSlicer(self, simData, maps=None):
        """
        Use simData[self.lonCol] and simData[self.latCol]
        (in radians) to set up KDTree.
        """
        self._runMaps(maps)
        self._buildTree(simData[self.lonCol], simData[
                        self.latCol], self.leafsize)
        self._setRad(self.radius)

        @wraps(self._sliceSimData)
        def _sliceSimData(islice):
            """Return indexes for relevant opsim data at slicepoint
            (slicepoint=lonCol/latCol value .. usually ra/dec)."""
            sx, sy, sz = self._treexyz(self.slicePoints['ra'][
                                       islice], self.slicePoints['dec'][islice])
            # Query against tree.
            initIndices = self.opsimtree.query_ball_point(
                (sx, sy, sz), self.rad)
            # Loop through all the images and check if the slicepoint is inside the corners of the chip
            # XXX--should check if there's a better/faster way to do this.
            # Maybe in the setupSlicer loop through each image, and use the contains_points method to test all the
            # healpixels simultaneously?  Then just have a dict with keys = healpix id and values = list of indices?
            # That way _sliceSimData is just doing a dict look-up and we can
            # get rid of the spatialkey kwargs.

            indices = []

            # add positions of chip corners as ra1-ra4.  1 in upper left, 4 in
            # lower left
            ZTF_RA_HWID = np.radians(3.50)
            ZTF_DEC_HWID = np.radians(3.50)
            ra1 = simData['fieldRA'][initIndices] - ZTF_RA_HWID
            ra2 = simData['fieldRA'][initIndices] + ZTF_RA_HWID
            ra3 = simData['fieldRA'][initIndices] + ZTF_RA_HWID
            ra4 = simData['fieldRA'][initIndices] - ZTF_RA_HWID

            dec1 = simData['fieldDec'][initIndices] + ZTF_DEC_HWID
            dec2 = simData['fieldDec'][initIndices] + ZTF_DEC_HWID
            dec3 = simData['fieldDec'][initIndices] - ZTF_DEC_HWID
            dec4 = simData['fieldDec'][initIndices] - ZTF_DEC_HWID

            # Gnomic project all the corners that are near the slice point,
            # centered on slice point
            x1, y1 = gnomonic_project_toxy(ra1, dec1,
                                           self.slicePoints['ra'][islice], self.slicePoints['dec'][islice])
            x2, y2 = gnomonic_project_toxy(ra2, dec2,
                                           self.slicePoints['ra'][islice], self.slicePoints['dec'][islice])
            x3, y3 = gnomonic_project_toxy(ra3, dec3,
                                           self.slicePoints['ra'][islice], self.slicePoints['dec'][islice])
            x4, y4 = gnomonic_project_toxy(ra4, dec4,
                                           self.slicePoints['ra'][islice], self.slicePoints['dec'][islice])

            for i, ind in enumerate(initIndices):
                # Use matplotlib to make a polygon on
                bbPath = mplPath.Path(np.array([[x1[i], y1[i]],
                                                [x2[i], y2[i]],
                                                [x3[i], y3[i]],
                                                [x4[i], y4[i]],
                                                [x1[i], y1[i]]]))
                # Check if the slicepoint is inside the image corners and
                # append to list if it is
                if bbPath.contains_point((0., 0.)) == 1:
                    indices.append(ind)

            return {'idxs': indices,
                    'slicePoint': {'sid': self.slicePoints['sid'][islice],
                                   'ra': self.slicePoints['ra'][islice],
                                   'dec': self.slicePoints['dec'][islice]}}
        setattr(self, '_sliceSimData', _sliceSimData)
