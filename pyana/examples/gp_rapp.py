import logging, argparse, os, re
import numpy as np
from .utils import getWorkDirs, checkSymLink
from .gp_rdiff import getUArray, getEdges, getCocktailSum, enumzipEdges
from decimal import Decimal

ranges = [0.3, 0.7]

def gp_rapp():
  inDir, outDir = getWorkDirs()
  eRanges = np.array([Decimal(str(e)) for e in ranges])
  for infile in os.listdir(inDir):
    energy = re.compile('\d+').search(infile).group()
    medium = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
    uMedium = getUArray(medium)
    eMedium = getEdges(medium)
    for i, (e0, e1) in enumzipEdges(eRanges):
      uSum = getCocktailSum(e0, e1, eMedium, uMedium)
      logging.debug('%s> %g - %g: %r' % (energy, e0, e1, uSum))
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
