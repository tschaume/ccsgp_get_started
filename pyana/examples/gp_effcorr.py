import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.config import default_colors
import numpy as np

def gp_syserr():
    energies = ['19.6', '27', '39', '62.4']
    inDir, outDir = getWorkDirs()
    inDir = os.path.join(inDir, 'syserr')
    data = OrderedDict()
    shift, factor = 0.7, -1.5
    ymin = OrderedDict((energy, 100) for energy in energies)
    ymax = OrderedDict((energy, -100) for energy in energies)
    for particle in ['e', 'pi']:
        for charge in ['plus', 'minus']:
            fname = particle + charge + '.dat'
            key = '{/Symbol \160}' if particle == 'pi' else 'e'
            key += '^{+}' if charge == 'plus' else '^{-}'
            data[key] = np.loadtxt(open(os.path.join(inDir, fname), 'rb'))
            # get y-values for error box
            for dp in data[key]:
                energy = '%g' % dp[0]
                curymin, curymax = dp[1]-dp[3], dp[1]+dp[3]
                if curymin < ymin[energy]: ymin[energy] = curymin
                if curymax > ymax[energy]: ymax[energy] = curymax
            # shift data for visibility
            data[key][:,0] += factor*shift
            factor += 1
    print ymax
    make_plot(
        data = data.values(),
        properties = [
            'lc %s lw 5 ps 2 pt %d' % (default_colors[10+i/2],18+i%2)
            for i in xrange(len(data))
        ],
        titles = data.keys(),
        tmargin = 0.98, rmargin = 0.99, bmargin = 0.13,
        yr = [0,22], xr = [15,68],
        xlabel = '{/Symbol \326}s_{NN} (GeV)',
        ylabel = 'Relative Systematic Uncertainty on {/Symbol \145}_{qual}{/Symbol \327}{/Symbol \145}_{glDCA} (%)',
        name = os.path.join(outDir, 'syserr'),
        size = '11in,8.5in', key = ['nobox', 'at graph 0.7,0.2', 'maxrows 2'],
        lines = dict(('y=%s' % energy,'lw 4 lt 2 lc "black"') for energy in energies),
        gpcalls = [
            'object %d rectangle back from %f,%f to %f,%f fc rgb "#FF9999" lw 2 fs border lc rgb "#FF6666"' % (
                50+i, float(energy)-1.5*shift, ymin[energy], float(energy)+1.5*shift, ymax[energy]
            ) for i,energy in enumerate(energies)
        ]
    )

def gp_tpc_select_eff():
    inDir, outDir = getWorkDirs()
    infile = os.path.join(inDir, 'tpc_select_eff_electrons_39GeV.dat')
    data = np.loadtxt(open(infile, 'rb'))
    nrows = len(data)
    data[:,1:] *= 100. # convert to %
    data = np.c_[data[:,:2], np.zeros(nrows), np.zeros(nrows), data[:,-1] ]
    make_plot(
        data = [data], titles = [''],
        properties = [ 'with filledcurves lt 1 lc %s lw 5 pt 0' % default_colors[1] ],
        tmargin = 0.98, rmargin = 0.99, yr = [80,100], xr = [0.2,2.052],
        gpcalls = ['nokey', 'xtics 0.5', 'mxtics 5', 'mytics 2'],
        xlabel = 'momentum, p (GeV/c)', ylabel = 'TPC Selection Efficiency (%)',
        name = os.path.join(outDir, 'tpc_select_eff'), size = '8.8in,6.8in',
        labels = {'39 GeV Electrons': [1.0,93,True]}
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
  #gp_syserr()
  gp_tpc_select_eff()
