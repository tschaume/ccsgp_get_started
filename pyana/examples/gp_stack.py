import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors

shift = { '200': 200., '62': 20., '39': 1., '19': 0.05 }
cocktail_style = 'with filledcurves pt 0 lc %s lw 5 lt 1' % default_colors[8]
pseudo_point = np.array([ [-1,1e-7,0,0,1] ])

def gp_stack(version):
  """example for a plot w/ stacked graphs using QM12 data (see gp_panel)

  * how to omit keys from the legend
  * manually add legend entries
  * automatically plot arrows for error bars larger than data point value

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
    # following scaling is wrong for y < 0 && shift != 1
    data_import[:, 1:] *= shift[energy]
    if fnmatch(file, 'data*'):
      data[energy] = data_import
    elif energy == '19': # cut of cocktail above 1.1 GeV/c^2
      cocktail[energy] = data_import[data_import[:,0] < 1.3]
    else:
      cocktail[energy] = data_import
  dataOrdered = OrderedDict(
    (' '.join([k, 'GeV', '{/Symbol \264} %g' % shift[k]]), data[k])
    for k in sorted(data, key=int)
  )
  nSetsData, nSetsCocktail = len(dataOrdered), len(cocktail)
  make_plot(
    data = cocktail.values() + [ pseudo_point ] + dataOrdered.values(),
    properties = [ cocktail_style ] * (nSetsCocktail+1) + [
      'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
      for i in xrange(nSetsData)
    ],
    titles = [''] * nSetsCocktail + ['Cocktail w/o {/Symbol \162}'] + dataOrdered.keys(),
    name = os.path.join(outDir, 'stack%s' % version),
    ylabel = 'dielectron pair production rate',
    xlabel = 'dielectron mass (GeV/c^{2})',
    ylog = True, xr = [0, 3.5], yr = [1e-6, 2e3],
    lmargin = 0.07, key = ['width -3'],
    #arrows = [ # example arrow
    #  [ [2.4, 5e-5], [2.3, 1e-5], 'head filled lc 1 lw 5 lt 1 front' ],
    #],
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
