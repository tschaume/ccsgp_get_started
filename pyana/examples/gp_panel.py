import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors

def gp_panel(version):
  """example for a panel plot using QM12 data (see gp_xfac)

  .. image:: pics/panelQM12.png
     :width: 700 px

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = {}
  for file in os.listdir(inDir):
    energy = re.compile('\d+').search(file).group()
    data_type = re.sub('%s\.dat' % energy, '', file)
    file_url = os.path.join(inDir, file)
    data_import = np.loadtxt(open(file_url, 'rb'))
    data_import = data_import[data_import[:,0] < 1.1]
    if data_type == 'cocktail': data_import[:,(2,4)] = 0.
    elif data_type == '+medium': data_import[:,2] = 0.
    if energy not in data: data[energy] = [ data_import ]
    else: data[energy].append(data_import)
  make_panel(
    dpt_dict = OrderedDict(
      (' '.join([k, 'GeV']), [
        data[k], [
          'with filledcurves lt 1 lw 4 pt 0 lc %s' % default_colors[8],
          'with lines lc 0 lw 5 lt 1',
          'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
        ], [ '+medium', 'cocktail', 'data' ]
      ]) for k in sorted(data, key=int)
    ), # 'lc 0' works here because no error plotting necessary
    name = os.path.join(outDir, 'panel%s' % version),
    ylabel = 'dielectron pair production rate',
    xlabel = 'dielectron mass (GeV/c^{2})',
    ylog = True, xr = [0, 1.1], yr = [1e-4, 20],
    lmargin = 0.08, bmargin = 0.15,
    arrow_length = 0.4, arrow_bar = 0.002,
    gpcalls = ['mxtics 2']
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
