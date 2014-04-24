import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, eRanges, getEnergy4Key
from .utils import getUArray, getEdges, getCocktailSum, enumzipEdges, getMassRangesSums
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts, zip_flat
from ..ccsgp.config import default_colors

xshift = 0.01
yunit = 1.0e-3

def gp_rdiff(version, nomed, noxerr):
  """example for ratio or difference plots using QM12 data (see gp_panel)

  - uses uncertainties package for easier error propagation and rebinning
  - stat. error for medium = 0!
  - stat. error for cocktail ~ 0!
  - statistical error bar on data stays the same for diff
  - TODO: implement ratio!
  - TODO: adjust statistical error on data for ratio!
  - TODO: adjust name and ylabel for ratio

  .. image:: pics/diffAbsQM12.png
     :width: 450 px

  :param version: plot version
  :type version: str
  :param nomed: don't plot medium
  :type nomed: bool
  :param noxerr: don't plot x-errors
  :type noxerr: bool
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data, cocktail, medium = OrderedDict(), OrderedDict(), OrderedDict()
  for file in os.listdir(inDir):
    energy = re.compile('\d+').search(file).group()
    data_type = re.sub('%s\.dat' % energy, '', file)
    energy = getEnergy4Key(energy)
    file_url = os.path.join(inDir, file)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if data_type == 'data': data[energy] = data_import[data_import[:,0] < 0.8]
    elif data_type == 'cocktail': cocktail[energy] = data_import
    elif not nomed: medium[energy] = data_import
  nSetsData = len(data)

  dataOrdered = OrderedDict()
  for energy in sorted(data, key=float):
    # data & bin edges
    uData = getUArray(data[energy])
    eData = getEdges(data[energy])
    uCocktail = getUArray(cocktail[energy])
    eCocktail = getEdges(cocktail[energy])
    loop = [eData]
    if energy in medium:
      uMedium = getUArray(medium[energy])
      eMedium = getEdges(medium[energy])
      loop.append(eMedium)
    # loop data/medium bins
    for l, eArr in enumerate(loop):
      for i, (e0, e1) in enumzipEdges(eArr):
        logging.debug('%s/%d> %g - %g:' % (energy, l, e0, e1))
        # get cocktail sum in data bin range
        uCocktailSum = getCocktailSum(e0, e1, eCocktail, uCocktail)
        # calc. difference and divide by data binwidth again
        # + set data point
        #xs = xshift if energy == '39' else 0.
        xs = 0. # TODO: cannot use shift when calculating excess yield based on dataOrdered
        if not l:
          uDiff = uData[i] - uCocktailSum
          uDiff /= data[energy][i,2] * 2 * yunit
          dp = [
            data[energy][i,0] + xs, uDiff.nominal_value,
            data[energy][i,2] if not noxerr else 0.,
            data[energy][i,3] / yunit, uDiff.std_dev
          ]
          key = ' '.join([energy, 'GeV'])
        else:
          uDiff = uMedium[i] - uCocktailSum
          uDiff /= medium[energy][i,2] * 2 * yunit
          # cut off medium/cocktail at omega
          if medium[energy][i,0] > 0.74:
            continue
          dp = [
            medium[energy][i,0] + xs, uDiff.nominal_value,
            medium[energy][i,2] if not noxerr else 0.,
            0., 0. # uDiff.std_dev
          ]
          key = ' '.join([energy, 'GeV (Med.)'])
        # build list of data points
        if key in dataOrdered: dataOrdered[key].append(dp)
        else: dataOrdered[key] = [ dp ]

  # make plot
  nSets = len(dataOrdered)
  nSetsPlot = nSets/2 if nSets > nSetsData else nSets
  ylabel = 'data/medium' if nSets > nSetsData and not nomed else 'data'
  props = [
    'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i] for i in xrange(nSetsPlot)
  ]
  titles = dataOrdered.keys()
  if nSets > nSetsData:
    props = zip_flat(props, [
      'with filledcurves pt 0 lt 1 lw 4 lc %s' % default_colors[i]
      for i in xrange(nSetsPlot)
    ])
    titles = zip_flat(dataOrdered.keys()[::2], [''] * nSetsPlot)
  labels = {
    #'{/Symbol \104}M_{ee}(39GeV) = +%g GeV/c^{2}' % xshift: [0.1, 0.9, False],
    'BES: STAR Preliminary' if version == 'QM12Latest200'
    else 'STAR Preliminary': [0.5,0.05,False],
    '200 GeV: [arXiv:1312.7397]' if version == 'QM12Latest200'
    else '': [0.25,0.93,False]
  }
  make_plot(
    data = [ np.array(d) for d in dataOrdered.values()],
    properties = props, titles = titles,
    name = os.path.join(outDir, 'diffAbs%s%s%s' % (
      version, 'NoMed' if nomed else '',
      'NoXErr' if noxerr else ''
    )),
    xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
    ylabel = '%s - (cocktail w/o {/Symbol \162}) ({/Symbol \264} 10^{-3})' % ylabel,
    xr = [0.2,0.76], yr = [-1,10.3], labels = labels,
    key = ['at graph 1.,1.1', 'maxrows 1', 'width -1.5'],
    lines = { 'x=0': 'lc 0 lw 4 lt 2' }
  )

  # integrated excess yield in mass ranges
  if nomed or noxerr or version == 'QM12': return 'done'
  excess = {}
  for k, v in dataOrdered.iteritems():
    suffix = ''
    if fnmatch(k, '*Med.*'): suffix = '_Med'
    energy = getEnergy4Key(re.compile('\d+').search(k).group())
    getMassRangesSums(energy, np.array(v), excess, onlyLMR = True, suffix = suffix)
  logging.debug(excess)
  make_plot(
    data = [ np.array(excess['LMR']), np.array(excess['LMR_Med']) ],
    properties = [
        'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0],
        'with linespoints lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[1],
        ],
    titles = [ 'data - cocktail', 'medium' ],
    name = os.path.join(outDir, 'excess%s' % version),
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    ylabel = 'LMR Excess Yield for %.2f < M_{ee} < %.2f ({/Symbol \264} 10^{-3})' % (
      eRanges[1], eRanges[2]
    ),
    lmargin = 0.08, xlog = True, #xr = [0.2,0.76],
    key = ['width -4'],
    yr = [0,3], gpcalls = [
      'format x "%g"',
      'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
    ], labels = labels
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--nomed", help="don't plot medium", action="store_true")
  parser.add_argument("--noxerr", help="no dx errors", action="store_true")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_rdiff(version, nomed, noxerr)
