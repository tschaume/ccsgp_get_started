import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors

def gp_panel(version, skip):
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
    if skip is not None and energy == skip: continue
    data_type = re.sub('%s\.dat' % energy, '', file)
    file_url = os.path.join(inDir, file)
    data_import = np.loadtxt(open(file_url, 'rb'))
    data_import = data_import[data_import[:,0] < 1.1]
    if data_type == 'cocktail': data_import[:,2:] = 0.
    elif data_type == '+medium': data_import[:,2] = 0.
    key = getEnergy4Key(energy)
    if key not in data: data[key] = {}
    data[key][data_type] = data_import
  plot_order = ['+medium', 'cocktail', 'data']
  plot_opts = {
    '+medium': 'with filledcurves lt 1 lw 4 pt 0 lc %s' % default_colors[8],
    'cocktail': 'with lines lc %s lw 5 lt 1' % default_colors[-2],
    'data': 'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
  }
  make_panel(
    dpt_dict = OrderedDict(
      (' '.join([k, 'GeV']), [
        [ data[k][dt] for dt in plot_order if dt in data[k] ],
        [ plot_opts[dt] for dt in plot_order if dt in data[k] ],
        [ dt for dt in plot_order if dt in data[k] ]
      ]) for k in sorted(data, key=float)
    ),
    name = os.path.join(
      outDir, 'panel%s%s' % (version, 'No'+skip if skip is not None else '')
    ),
    ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]',
    xlabel = 'invariant dielectron mass, M_{ee} (GeV/c^{2})',
    ylog = True, xr = [0, 1.1], yr = [1e-4, 20],
    lmargin = 0.1, bmargin = 0.15,
    arrow_length = 0.4, arrow_bar = 0.002,
    gpcalls = ['mxtics 2'] + (['label %d "" at graph 0.4,0.7' % (
      8 if skip is None else 6
    )] if version == 'QM12Latest200' else []),
    labels = {'STAR Preliminary': [0.4,0.7,False]}
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--skip", help="skip an energy", metavar="energy")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_panel(args.version, args.skip)
