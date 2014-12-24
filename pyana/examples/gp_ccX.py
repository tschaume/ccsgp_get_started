import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
from scipy.optimize import curve_fit
import numpy as np
from uncertainties import ufloat

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
    popt = np.hstack((popt1, popt2))
    model = lambda x: fitfunc(x, *popt)
    # calculate mean standard deviation of data from parameterization
    yfit = np.array([model(x) for x in alldata[:,0]])
    stddev = 1.5*np.sqrt( # multiple of "sigma"
        np.average((alldata[:,1]-yfit)**2, weights=1./alldata[:,3])
    )
    print 'stddev = %.2g' % stddev
    errorband = np.array([[x, model(x), 0, 0, stddev] for x in np.linspace(1,4)])
    # make plot
    fitdata = np.array([[x, model(x), 0, 0, 0] for x in np.linspace(1,4)])
    par_names = ['a', 'b', 'c', 'd']
    energies = [19.6, 27, 39, 62.4, 200]
    labels = dict(
        ('%s = %.3g' % (par_name, popt[i]), [3.3, 3-i*0.2, True])
        for i,par_name in enumerate(par_names)
    )
    ccX_vals = [10**model(np.log10(energy)) for energy in energies]
    ccX = [' '.join([
        '%g GeV:' % energy,
        '({})'.format(ufloat(ccX_vals[i], stddev/0.434*ccX_vals[i])),
        '{/Symbol \155}b'
    ]) for i,energy in enumerate(energies)]
    print ccX
    #labels.update(dict(
    #    (cc, [1+i*0.5, 4.5+(i%2+1)*0.2, True]) for i,cc in enumerate(ccX)
    #))
    make_plot(
        data = [errorband] + data.values() + [fitdata],
        properties = [
            'with filledcurves lt 1 lw 5 pt 0 lc %s' % default_colors[8]
        ] + [
            'lc %s lw 4 lt 1 pt 18 ps 1.5' % (default_colors[i])
            for i in xrange(len(data))
        ] + ['with lines lc 0 lw 4 lt 1'],
        titles = [''] + data.keys() + ['y = ax+b - e^{-cx+d}'],
        xlabel = 'x = log_{10}[{/Symbol \326}s_{NN} (GeV)]',
        ylabel = 'y = log_{10}[{/Symbol \163}@_{c@^{/=18-}c}^{NN} ({/Symbol \155}b)]',
        name = os.path.join(outDir, 'ccX'),
        size = '11.4in,8.3in', xr = [1, 4], yr = [0.5,4.5],
        key = ['bottom right', 'nobox', 'width -5'],
        tmargin = 0.98, rmargin = 0.99, bmargin = 0.13,
        lines = dict(
            ('y=%f' % (np.log10(energy)), 'lc 0 lw 2 lt 2') for energy in energies
        ), labels = labels,
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
