import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors

def gp_rdiff(version):
  """example for ratio or difference plots using QM12 data (see gp_panel)

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
    data_import = data_import[data_import[:,0] < 1.3]
    if fnmatch(file, 'data*'): data[energy] = data_import
    else: cocktail[energy] = data_import
  # TODO: rebin cocktail to data
  # TODO: implement difference
  # TODO: implement ratio!
  dataOrdered = OrderedDict(
    (' '.join([k, 'GeV']), '''data[k] - cocktailRebin[k]''')
    for k in sorted(data, key=int)
  )
  nSets = len(dataOrdered)
  make_plot(
    data = dataOrdered.values(),
    properties = [
      'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i] for i in xrange(nSets)
    ],
    titles = dataOrdered.keys(),
    # TODO: adjust name and ylabel for ratio
    name = os.path.join(outDir, 'diff%s' % version), ylabel = 'diff',
    xlabel = 'dielectron mass (GeV/c^{2})',
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
  print gp_stack(args.version)
