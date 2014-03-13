import logging, argparse, os, re
import numpy as np
from ..ccsgp.config import default_colors
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs, checkSymLink, getMassRangesSums

def gp_rapp():
  """rho in-medium ratios by Rapp (based on protected data)"""
  inDir, outDir = getWorkDirs()
  # prepare data
  yields = {}
  for infile in os.listdir(inDir):
    energy = re.compile('\d+').search(infile).group()
    medium = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
    getMassRangesSums(energy, medium, yields)
  data = dict( # sort by energy
    (k, np.array(sorted(v)))
    for k, v in yields.iteritems()
  )
  for k in data: data[k][:,1] /= data[k][-1,1] # divide by 200
  # make plot
  nSets = len(data)
  make_plot(
    data = [ data[k][:-1] for k in data ],
    properties = [
      'with linespoints lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
      for i in xrange(nSets)
    ],
    titles = data.keys(), # TODO: titles order correct?
    name = os.path.join(outDir, 'rapp'),
    xlabel = '{/Symbol \326}s_{NN} (GeV)', ylabel = 'Rapp Ratio to 200 GeV',
    lmargin = 0.1, key = ['left'], yr = [0.1, 0.8]
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_rapp()
