import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_panel
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors
from fnmatch import fnmatch

def gp_panel(version, skip):
  """example for a panel plot using QM12 data (see gp_xfac)

  .. image:: pics/panelQM12.png
     :width: 700 px

  :param version: plot version / input subdir name
  :type version: str
  """
  #scale = { # LatestPatrickJieYi
  #    '19': 0.4274654744079354, '200': 1.0, '39': 0.4362451929487654,
  #    '27': 0.47464918475541873, '62': 0.5800852553921563
  #}
  scale = { # QM14
    '19': 0.47144704299427165, '200': 1.0, '39': 0.7776170174498098,
    '27': 0.47464918475541873, '62': 0.9173998879333009
  }
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = {}
  for infile in os.listdir(inDir):
    if infile == "cocktail_contribs": continue
    energy = re.compile('\d+').search(infile).group()
    if skip is not None and energy == skip: continue
    data_type = re.sub('%s\.dat' % energy, '', infile)
    file_url = os.path.join(inDir, infile)
    data_import = np.loadtxt(open(file_url, 'rb'))
    data_import = data_import[data_import[:,0] < 1.1]
    if data_type == 'data' and (
        version == 'LatestPatrickJieYi' or (
            version == 'QM14' and energy != '19'
            and energy != '27' # skip when using Joey's dataset
        )
    ):
        data_import[:,(1,3,4)] *= scale[energy]
    if data_type == 'cocktail': data_import[:,2:] = 0.
    elif fnmatch(data_type, '*medium*'): data_import[:,2] = 0.
    key = getEnergy4Key(energy)
    if key not in data: data[key] = {}
    data_type_mod = data_type
    if data_type_mod == 'mediumMedOnly': data_type_mod = 'in-medium'
    elif data_type_mod == 'mediumQgpOnly': data_type_mod = 'QGP'
    elif data_type_mod == '+medium': data_type_mod = 'cocktail + model'
    data[key][data_type_mod] = data_import
  plot_order = ['in-medium', 'QGP', 'cocktail + model', 'cocktail', 'data']
  plot_opts = {
    'QGP': 'with lines lt 2 lw 5 lc %s' % default_colors[1],
    'in-medium': 'with lines lt 2 lw 5 lc %s' % default_colors[2],
    'cocktail + model': 'with filledcurves lt 1 lw 5 pt 0 lc %s' % default_colors[16],
    'cocktail': 'with lines lc %s lw 5 lt 1' % default_colors[8],
    'data': 'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
  }
  panel2D_versions = (version == 'LatestPatrickJieYi' or version == 'QM14')
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
    lmargin = 0.12 if panel2D_versions else 0.1,
    bmargin = 0.11 if panel2D_versions else 0.15,
    arrow_length = 0.4, arrow_bar = 0.002,
    gpcalls = ['mxtics 2'] + (['label %d "" at graph 0.4,0.7' % (
      8 if skip is None else 6
    )] if version == 'QM12Latest200' else []),
    labels = {'STAR Preliminary': [0.4,0.5,False]},
    layout = '3x2' if panel2D_versions else ('%dx1' % len(data)),
    key = ['width -4', 'at graph 0.95,0.85']
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
