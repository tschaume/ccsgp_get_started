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

xshift = 0.01
yunit = 1.0e-3

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
    data_import = data_import[data_import[:,0] < 1.5]
    if data_type == 'data': data[energy] = data_import
    else: cocktail[energy] = data_import
  dataOrdered = OrderedDict()
  for energy in sorted(data, key=int):
    # data uncertainty array multiplied by binwidth (col2 = dx)
    # TODO: 39 GeV has a datapoint with negative error!
    uData = unp.uarray(data[energy][:,1], data[energy][:,4])
    uData *= data[energy][:,2] * 2
    # cocktail uncertainty array multiplied by binwidth
    # stat. error for cocktail ~ 0!
    uCocktail = unp.uarray(cocktail[energy][:,1], cocktail[energy][:,4])
    uCocktail *= cocktail[energy][:,2] * 2
    # data bin edges and cocktail bin centers
    edges = np.concatenate(([0], data[energy][:,0] + data[energy][:,2]))
    xc = cocktail[energy][:,0]
    # loop data bins
    # TODO: check whether cocktail & data edges coincide!
    for i, (e0, e1) in enumerate(zip(edges[:-1], edges[1:])):
      # calc. difference and divide by data binwidth again
      # TODO: implement ratio!
      uDiff = uData[i] - fsum(uCocktail[(xc > e0) & (xc < e1)])
      uDiff /= data[energy][i,2] * 2 * yunit
      # set data point
      xs = xshift if energy == '39' else 0.
      dp = [
        data[energy][i,0] + xs, uDiff.nominal_value,
        # statistical error bar on data stays the same for diff
        # TODO: adjust statistical error on data for ratio!
        0., data[energy][i,3] / yunit, uDiff.std_dev
      ]
      # build list of data points
      key = ' '.join([energy, 'GeV'])
      if key in dataOrdered: dataOrdered[key].append(dp)
      else: dataOrdered[key] = [ dp ]
  # make plot
  nSets = len(dataOrdered)
  make_plot(
    data = [ np.array(d) for d in dataOrdered.values()],
    properties = [
      'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i] for i in xrange(nSets)
    ],
    titles = dataOrdered.keys(),
    # TODO: adjust name and ylabel for ratio
    name = os.path.join(outDir, 'diffAbs%s' % version),
    xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
    ylabel = 'data - (cocktail w/o {/Symbol \162}) ({/Symbol \264} 10^{-3})',
    xr = [0.2,0.77], yr = [-1,9],
    labels = {
      '{/Symbol \104}M_{ee}(39GeV) = +%g GeV/c^{2}' % xshift: [0.1, 0.9, False]
    },
    key = ['at graph 1.,1.1', 'maxrows 1'],
    lines = { 'x=0': 'lc 0 lw 4 lt 2' }
  )
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
