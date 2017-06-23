
import numpy as np

from lsst.sims.utils import _equatorialFromGalactic
from lsst.sims.maf.plots import BaseSkyMap
from matplotlib.patches import Polygon


class ZTFBaseSkyMap(BaseSkyMap):

    def __init__(self, aspect_ratio = 1.0):
        super(ZTFBaseSkyMap, self).__init__()
        # width / length
        self.aspect_ratio = aspect_ratio

    def _plot_tissot_ellipse(self, lon, lat, radius, ax=None, **kwargs):
        """Hack to plot PTF/ZTF rectangular/square FOV
        Parameters
        ----------
        lon : float or array_like
        longitude-like of field centers (radians)
        lat : float or array_like
        latitude-like of field centers (radians)
        radius : float or array_like
        radius of *circumscribed* circle (radians)
        ax : Axes object (optional)
        matplotlib axes instance on which to draw polygon patches
        Other Parameters
        ----------------
        other keyword arguments will be passed to matplotlib.patches.Polygon.
        # The code in this method adapted from astroML, which is BSD-licensed.
        # See http: //github.com/astroML/astroML for details.
        """
        # Code adapted from astroML, which is BSD-licensed.
        # See http: //github.com/astroML/astroML for details.
        polygons = []
        if ax is None:
            ax = plt.gca()
        for l, b, r in np.broadcast(lon, lat, radius):
            theta = np.arctan(1/self.aspect_ratio)
            dy = r * np.sin(theta)
            dx = r * np.cos(theta)
            bp = b + dy 
            bm = b - dy
            xy = np.array([[l - dx / np.cos(bm), bm],
                           [l + dx / np.cos(bm), bm],
                           [l + dx / np.cos(bp), bp],
                           [l - dx / np.cos(bp), bp]])
            poly = Polygon(xy, closed=True, **kwargs)
            polygons.append(poly)
        return polygons

    def _plot_mwZone(self, raCen=0, peakWidth=np.radians(10.), taperLength=np.radians(80.), ax=None):
        """
        Plot blue lines to mark the milky way galactic exclusion zone.
        """
        if ax is None:
            ax = plt.gca()
        # Calculate galactic coordinates for mw location.
        step = 0.02
        galL = np.arange(-np.pi, np.pi + step / 2., step)
        galB = np.zeros(len(galL))
        galBm20 = np.zeros(len(galL)) - np.radians(20.)
        galBp20 = np.zeros(len(galL)) + np.radians(20.)
        #val = peakWidth * np.cos(galL / taperLength * np.pi / 2.)
        #galB1 = np.where(np.abs(galL) <= taperLength, val, 0)
        #galB2 = np.where(np.abs(galL) <= taperLength, -val, 0)
        # Convert to ra/dec.
        # Convert to lon/lat and plot.
        ra, dec = _equatorialFromGalactic(galL, galB)
        lon = -(ra - raCen - np.pi) % (np.pi * 2) - np.pi
        ax.plot(lon, dec, 'b.', markersize=1.8, alpha=0.4)

        ra, dec = _equatorialFromGalactic(galL, galBm20)
        lon = -(ra - raCen - np.pi) % (np.pi * 2) - np.pi
        ax.plot(lon, dec, 'b.', markersize=1.8, alpha=0.2)

        ra, dec = _equatorialFromGalactic(galL, galBp20)
        lon = -(ra - raCen - np.pi) % (np.pi * 2) - np.pi
        ax.plot(lon, dec, 'b.', markersize=1.8, alpha=0.2)
