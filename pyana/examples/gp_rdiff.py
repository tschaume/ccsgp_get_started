import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts, zip_flat
from ..ccsgp.config import default_colors
from uncertainties import unumpy as unp
from uncertainties.umath import fsum
from decimal import Decimal

xshift = 0.01
yunit = 1.0e-3

def getUArray(npArr):
  """uncertainty array multiplied by binwidth (col2 = dx)"""
  return unp.uarray(npArr[:,1], npArr[:,4]) * npArr[:,2] * 2

def getEdges(npArr):
  """get np array of bin edges"""
  edges = np.concatenate(([0], npArr[:,0] + npArr[:,2]))
  return np.array([Decimal(str(i)) for i in edges])

def getMaskIndices(mask):
  """get lower and upper index of mask"""
  return [
    list(mask).index(True), len(mask) - 1 - list(mask)[::-1].index(True)
  ]

def gp_rdiff(version):
  """example for ratio or difference plots using QM12 data (see gp_panel)

  - uses uncertainties package for easier error propagation and rebinning
  - stat. error for medium = 0!
  - stat. error for cocktail ~ 0!
  - statistical error bar on data stays the same for diff
  - TODO: implement ratio!
  - TODO: adjust statistical error on data for ratio!
  - TODO: adjust name and ylabel for ratio

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data, cocktail, medium = OrderedDict(), OrderedDict(), OrderedDict()
  for file in os.listdir(inDir):
    energy = re.compile('\d+').search(file).group()
    data_type = re.sub('%s\.dat' % energy, '', file)
    file_url = os.path.join(inDir, file)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if data_type == 'data': data[energy] = data_import[data_import[:,0] < 0.8]
    elif data_type == 'cocktail': cocktail[energy] = data_import
    else: medium[energy] = data_import

  dataOrdered = OrderedDict()
  for energy in sorted(data, key=int):
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
      for i, (e0, e1) in enumerate(zip(eArr[:-1], eArr[1:])):
        logging.debug('%s/%d> %g - %g:' % (energy, l, e0, e1))
        # get mask and according indices
        mask = (eCocktail >= e0) & (eCocktail <= e1)
        # data/medium bin range larger than single cocktail bin
        if np.any(mask):
          idx = getMaskIndices(mask)
          # determine coinciding flags
          eCl, eCu = eCocktail[idx[0]], eCocktail[idx[1]]
          not_coinc_low, not_coinc_upp = (eCl != e0), (eCu != e1)
          # get cocktail sum in data bin (always w/o last bin)
          uCocktailSum = fsum(uCocktail[mask[:-1]][:-1])
          logging.debug('    sum: {}'.format(uCocktailSum))
          # get correction for non-coinciding edges
          if not_coinc_low:
            eCl_bw = eCl - eCocktail[idx[0]-1]
            corr_low = (eCl - e0) / eCl_bw
            abs_corr_low = float(corr_low) * uCocktail[idx[0]-1]
            uCocktailSum += abs_corr_low
            logging.debug(('    low: %g == %g -> %g (%g) -> %g -> {} -> {}' % (
              e0, eCl, eCl - e0, eCl_bw, corr_low
            )).format(abs_corr_low, uCocktailSum))
          if not_coinc_upp:
            eCu_bw = eCocktail[idx[1]+1] - eCu
            corr_upp = (e1 - eCu) / eCu_bw
            abs_corr_upp = float(corr_upp) * uCocktail[idx[1]]
            uCocktailSum += abs_corr_upp
            logging.debug(('    upp: %g == %g -> %g (%g) -> %g -> {} -> {}' % (
              e1, eCu, e1 - eCu, eCu_bw, corr_upp
            )).format(abs_corr_upp, uCocktailSum))
        else:
          mask = (eCocktail >= e0)
          idx = getMaskIndices(mask) # only use first index
          corr = (e1 - e0) / (eCocktail[idx[0]+1] - eCocktail[idx[0]])
          uCocktailSum = float(corr) * uCocktail[idx[0]]
          logging.debug('    sum: {}'.format(uCocktailSum))
        # calc. difference and divide by data binwidth again
        # + set data point
        xs = xshift if energy == '39' else 0.
        if not l:
          uDiff = uData[i] - uCocktailSum
          uDiff /= data[energy][i,2] * 2 * yunit
          dp = [
            data[energy][i,0] + xs, uDiff.nominal_value,
            0., data[energy][i,3] / yunit, uDiff.std_dev
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
            0., 0., 0. # uDiff.std_dev
          ]
          key = ' '.join([energy, 'GeV (Med.)'])
        # build list of data points
        if key in dataOrdered: dataOrdered[key].append(dp)
        else: dataOrdered[key] = [ dp ]

  # make plot
  nSets = len(dataOrdered)
  nSetsPlot = nSets/2 if nSets > 4 else nSets
  ylabel = 'data/medium' if nSets > 4 else 'data'
  props = [
    'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i] for i in xrange(nSetsPlot)
  ]
  titles = dataOrdered.keys()
  if nSets > 4:
    props = zip_flat(props, [
      'with filledcurves pt 0 lt 1 lw 4 lc %s' % default_colors[i]
      for i in xrange(nSetsPlot)
    ])
    titles = zip_flat(dataOrdered.keys()[::2], [''] * nSetsPlot)
  make_plot(
    data = [ np.array(d) for d in dataOrdered.values()],
    properties = props, titles = titles,
    name = os.path.join(outDir, 'diffAbs%s' % version),
    xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
    ylabel = '%s - (cocktail w/o {/Symbol \162}) ({/Symbol \264} 10^{-3})' % ylabel,
    xr = [0.2,0.76], yr = [-1,9],
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
