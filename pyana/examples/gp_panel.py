import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.utils import getOpts

def gp_panel(version):
  """example for a panel plot using QM12 data (see gp_xfac)"""
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = {}
  for file in os.listdir(inDir):
    energy = re.compile('\d+').search(file).group()
    file_url = os.path.join(inDir, file)
    data[energy] = np.loadtxt(open(file_url, 'rb')).reshape((-1,5))
  make_panel(
    dpt_dict = OrderedDict(
      (' '.join([k, 'GeV']), [[data[k]], [getOpts(0)], ['data'] ])
      for k in sorted(data, key=int)
    ),
    name = os.path.join(outDir, 'panel%s' % version),
    ylabel = 'invariant yield',
    xlabel = 'invariant mass (GeV/c^{2})',
    ylog = True, xr = [0, 1.29], yr = [1e-4, 20],
    lmargin = 0.08, bmargin = 0.15
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
  print gp_panel(args.version)