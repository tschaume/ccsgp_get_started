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
  #scale = { # QM14 (19 GeV skip later, factor here only informational)
  #  '19': 1.0340571932983775, '200': 1.0, '39': 0.7776679085382481,
  #  '27': 0.6412140408244136, '62': 0.9174700031778402
  #}
  scale = {
    '19': 1.0812324298238396, '200': 1.1051002240771077,
    '39': 1.12093451186094, '27': 1.206453891072013,
    '62': 1.4992194614005152
  }
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = {}
  vacRhoTitle = '{/Symbol \162}/{/Symbol \167} {/Symbol \256} e^{+}e^{-} (vac.)'
  for infile in os.listdir(inDir):
    if infile == "cocktail_contribs": continue
    if infile == 'mediumDmOnly200.dat': continue
    energy = re.compile('\d+').search(infile).group()
    if skip is not None and energy == skip: continue
    data_type = re.sub('%s\.dat' % energy, '', infile)
    file_url = os.path.join(inDir, infile)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if (data_type == 'cocktail' or fnmatch(data_type, '*medium*')) \
       and (version == 'QM14' and energy != '19'):
       data_import[:,(1,3,4)] /= scale[energy]
    elif data_type == 'data' and version == 'LatestPatrickJieYi':
       data_import[:,(1,3,4)] *= scale[energy]
    if data_type == 'cocktail': data_import[:,2:] = 0.
    elif fnmatch(data_type, '*medium*') or data_type == 'vacRho':
       data_import = data_import[data_import[:,0] < 0.9] \
               if energy == '200' and data_type == '+medium' \
               else data_import
       data_import[:,2] = 0.
    key = getEnergy4Key(energy)
    if key not in data: data[key] = {}
    data_type_mod = data_type
    if data_type_mod == 'mediumMedOnly': data_type_mod = 'HMBT'
    elif data_type_mod == 'mediumQgpOnly': data_type_mod = 'QGP'
    elif data_type_mod == '+medium': data_type_mod = 'cocktail + HMBT'
    elif data_type_mod == 'vacRho': data_type_mod = vacRhoTitle
    data[key][data_type_mod] = data_import
  plot_order = ['cocktail', vacRhoTitle, 'QGP', 'HMBT', 'cocktail + HMBT', 'data']
  plot_opts = {
    vacRhoTitle: 'with lines lt 2 lw 5 lc %s' % default_colors[6],
    'QGP': 'with lines lt 2 lw 5 lc %s' % default_colors[1],
    'HMBT': 'with lines lt 2 lw 5 lc %s' % default_colors[2],
    'cocktail + HMBT': 'with filledcurves lt 1 lw 5 pt 0 lc %s' % default_colors[16],
    'cocktail': 'with lines lc %s lw 5 lt 1' % default_colors[8],
    'data': 'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
  }
  panel2D_versions = (version == 'LatestPatrickJieYi' or version == 'QM14')
  make_panel(
    dpt_dict = OrderedDict(
      (' '.join([k, 'GeV %s' % (
          '{/=18 PRL 113 022301}'
          if k == '200' and (
              version == 'QM12Latest200' or version == 'QM14' or
              version == 'LatestPatrickJieYi'
          ) else '{/=18 STAR Preliminary}'
      )]), [
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
    ylog = True, xr = [0.05, 1.1], yr = [1e-4, 0.5],
    lmargin = 0.12 if panel2D_versions else 0.1,
    bmargin = 0.11 if panel2D_versions else 0.15,
    arrow_length = 0.4, arrow_bar = 0.002,
    gpcalls = [
        'mxtics 2',
        'object 50 rectangle back fc rgb "#C6E2FF" from 0.4,1e-4 to 0.74,2e-2'
    ],
    layout = '3x2' if panel2D_versions else ('%dx1' % len(data)),
    key = ['width -3', 'at graph 0.95,0.85'],
    key_subplot_id = 5, size = '8in,8in'
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
