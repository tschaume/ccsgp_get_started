import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
from scipy.optimize import curve_fit
import numpy as np

def linear(x, a, b):
    return a*x+b

def exp(x, c, d):
    return np.exp(-c*x+d)

def fitfunc(x, a, b, c, d):
    return linear(x, a, b) - exp(x, c, d)

def gp_ccX():
    """fit experimental data"""
    inDir, outDir = getWorkDirs()
    data, alldata = OrderedDict(), None
    for infile in os.listdir(inDir):
        # get key and import data
        key = os.path.splitext(infile)[0].replace('_', '/')
        data_import = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
        # convert to log10 vs log10 plot, z=log10(y) => dz=0.434*dy/y
        data_import[:,3] = 0.434*(data_import[:,3]/data_import[:,1])
        data_import[:,(0,1)] = np.log10(data_import[:,(0,1)])
        data_import[:,2] = 0
        # fill dictionary
        data[key] = data_import
        alldata = data[key] if alldata is None else np.vstack((alldata, data[key]))
    # fit linear part first
    lindata = alldata[alldata[:,0]>2.5]
    m = (lindata[-1,1]-lindata[0,1])/(lindata[-1,0]-lindata[0,0])
    t = lindata[0,1] - m * lindata[0,0]
    popt1, pcov = curve_fit(
        linear, lindata[:,0], lindata[:,1], p0=[m, t],
        sigma=lindata[:,3], absolute_sigma=True
    )
    # fit full range
    popt2, pcov = curve_fit(
        lambda x, c, d: fitfunc(x, popt1[0], popt1[1], c, d),
        alldata[:,0], alldata[:,1], sigma=alldata[:,3], absolute_sigma=True,
    )
    fitdata = np.array([[
        x, fitfunc(x, popt1[0], popt1[1], popt2[0], popt2[1]), 0, 0, 0
    ] for x in np.linspace(1,4)])
    make_plot(
        data = data.values() + [fitdata],
        properties = [
            'lc %s lw 4 lt 1 pt 18 ps 1.5' % (default_colors[i])
            for i in xrange(len(data))
        ] + ['with lines lc 0 lw 4 lt 2'],
        titles = data.keys() + ['parameterization'],
        xlabel = 'log_{10}[{/Symbol \326}s_{NN} (GeV)]',
        ylabel = 'log_{10}[{/Symbol \163}@_{c@^{/=18-}c}^{NN} ({/Symbol \155}b)]',
        name = os.path.join(outDir, 'ccX'),
        size = '8.5in,8in', xr = [1, 4], yr = [0.5,4.3],
        key = ['bottom right', 'nobox', 'width -3'],
        tmargin = 0.99, rmargin = 0.97, bmargin = 0.13
    )

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  gp_ccX()
