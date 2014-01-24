import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors
from uncertainties import unumpy as unp
from uncertainties.umath import fsum

def gp_rdiff(version):
  """example for ratio or difference plots using QM12 data (see gp_panel)

  - uses uncertainties package for easier error propagation and rebinning

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data, cocktail = OrderedDict(), OrderedDict()
  for file in os.listdir(inDir):
    energy = re.compile('\d+').search(file).group()
    data_type = re.sub('%s\.dat' % energy, '', file)
    file_url = os.path.join(inDir, file)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if fnmatch(file, 'data*'): data[energy] = data_import
    else: cocktail[energy] = data_import
  dataOrdered = OrderedDict()
  for energy in sorted(data, key=int):
    # TODO: 39 GeV has a datapoint with negative error!
    # TODO: how about statistical error bar on data?
    uData = unp.uarray(data[energy][:,1], data[energy][:,4])
    uData *= data[energy][:,2] # multiply data by binwidth
    # stat. error for cocktail = 0!
    uCocktail = unp.uarray(cocktail[energy][:,1], cocktail[energy][:,4])
    uCocktail *= cocktail[energy][:,2] # multiply cocktail by binwidth
    edges = np.concatenate(([0], data[energy][:,0] + data[energy][:,2] * 0.5))
    xc = cocktail[energy][:,0]
    for i, (e0, e1) in enumerate(zip(edges[:-1], edges[1:])):
      # TODO: check whether cocktail & data edges coincide!
      # TODO: implement ratio!
      uDiff = uData[i] - fsum(uCocktail[(xc > e0) & (xc < e1)])
      dp = [
        data[energy][i,0], uDiff.nominal_value,
        # TODO: adjust statistical error on data for ratio!
        0., data[energy][i,3], uDiff.std_dev
      ]
      key = ' '.join([energy, 'GeV'])
      print dp
      if key in dataOrdered:
        print 'appending ', dp
        np.append(dataOrdered[key], [dp], axis=0)
        dataOrdered[key]
      else: dataOrdered[key] = np.array([ dp ])
      print '%r' % dataOrdered[key]
    break
  #nSets = len(dataOrdered)
  #make_plot(
  #  data = dataOrdered.values(),
  #  properties = [
  #    'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i] for i in xrange(nSets)
  #  ],
  #  titles = dataOrdered.keys(),
  #  # TODO: adjust name and ylabel for ratio
  #  name = os.path.join(outDir, 'diff%s' % version), ylabel = 'diff',
  #  xlabel = 'dielectron mass (GeV/c^{2})',
  #)
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_rdiff(args.version)
