import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np


def gp_ccX():
    """fit experimental data"""
    inDir, outDir = getWorkDirs()
    data = OrderedDict()
    for infile in os.listdir(inDir):
        key = os.path.splitext(infile)[0].replace('_', '/')
        data[key] = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
    print data
    make_plot(
        data = data.values(),
        properties = [
            'lc %s lw 4 lt 1 pt 18 ps 1.5' % (default_colors[i])
            for i in xrange(len(data))
        ],
        titles = data.keys(),
        xlabel = '{/Symbol \326}s_{NN} (GeV)',
        ylabel = '{/Symbol \163}@_{c@^{/=18-}c}^{NN} ({/Symbol \155}b)',
        name = os.path.join(outDir, 'ccX'),
        ylog = True, xlog = True, size = '8.5in,8in',
        yr = [1,2e4], xr = [10, 1e4], key = ['bottom right', 'nobox'],
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
