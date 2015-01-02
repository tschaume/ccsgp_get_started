import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.config import default_colors
import numpy as np

def gp_syserr():
    inDir, outDir = getWorkDirs()
    inDir = os.path.join(inDir, 'syserr')
    data = {}
    for particle in ['e', 'pi']:
        for charge in ['plus', 'minus']:
            fname = particle + charge + '.dat'
            data[fname] = np.loadtxt(open(os.path.join(inDir, fname), 'rb'))
    make_plot(
        data = data.values(),
        properties = [ 'lc %s lw 4 pt 18' % default_colors[i] for i in xrange(len(data)) ],
        titles = data.keys(),
        xlabel = 'sqrt(s)', ylabel = 'syserr (%)',
        name = os.path.join(outDir, 'syserr'),
        size = '10in,8in',
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
  gp_syserr()
