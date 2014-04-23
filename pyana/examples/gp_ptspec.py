import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_panel, make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors
from decimal import Decimal
import uncertainties.umath as umath
import uncertainties.unumpy as unp

mee_keys = ['pi0', 'LMR', 'omega', 'phi', 'IMR', 'jpsi']

def getMeeLabel(s):
  if s == 'pi0': return '{/Symbol \160}^0'
  if s == 'omega': return '{/Symbol \167}'
  if s == 'phi': return '{/Symbol \152}'
  if s == 'jpsi': return 'J/{/Symbol \171}'
  return s

def splitFileName(fn):
  # energy, mee_name, mee_range, data_type
  split_arr = fn.split('_')
  return (
    re.compile('\d+').search(split_arr[0]).group(),
    split_arr[1], split_arr[2],
    re.compile('[a-z]+').search(split_arr[0]).group()
  )

def gp_ptspec():
  """example for a 2D-panel plot"""
  # import data and calculate average pT's
  inDir, outDir = getWorkDirs()
  data, data_avpt = OrderedDict(), OrderedDict()
  yvalsPt = []
  for filename in os.listdir(inDir):
    file_url = os.path.join(inDir, filename)
    filebase = os.path.splitext(filename)[0] # unique
    data[filebase] = np.loadtxt(open(file_url, 'rb'))
    probs = unp.uarray(data[filebase][:,1], data[filebase][:,3]) # dN/pTdpT
    probs *= data[filebase][:,0] # = dN/pTdpT * pT
    probs /= umath.fsum(probs) # probabilities
    pTs = unp.uarray(data[filebase][:,0], data[filebase][:,2]) # pT
    avpt = umath.fsum(pTs*probs)
    weights = data[filebase][:,1] * data[filebase][:,0]
    logging.info(('%s: {} %g' % (filebase, np.average(data[filebase][:,0], weights=weights))).format(avpt))
    energy, mee_name, mee_range, data_type = splitFileName(filebase)
    dp = [ float(energy), avpt.nominal_value, 0., avpt.std_dev, 0. ] # TODO: syst. uncertainties
    avpt_key = mee_name if data_type == 'data' else mee_name + '_c'
    if avpt_key in data_avpt: data_avpt[avpt_key].append(dp)
    else: data_avpt[avpt_key] = [ dp ]
    yvalsPt.append(avpt)
  # generate input dpt_dict
  dpt_dict = OrderedDict(), OrderedDict()
  yvals = []
  yscale = { '62': '1e6', '39': '1e4', '27': '1e2', '19': '1.' }
  for k in sorted(data.keys()):
    energy, mee_name, mee_range, data_type = splitFileName(k)
    if energy == '200': continue
    data[k][:,(1,3,4)] *= float(yscale[energy])
    if data_type == 'cocktail': data[k][:,2:] = 0.
    yvals += [v for v in data[k][:,1] if v > 0]
    if mee_name not in dpt_dict: dpt_dict[mee_name] = [ [], [], [] ]
    col = len(dpt_dict[mee_name][1]) # TODO: fix colors to match between data/cocktail
    dpt_dict[mee_name][0].append(data[k])
    dpt_dict[mee_name][1].append(
        'lt 1 lw 4 ps 1.5 lc %s pt 18' % (default_colors[col])
        if data_type == 'data' else
        'with lines lt 1 lw 4 lc %s' % (default_colors[col])
        )
    dpt_dict[mee_name][2].append(
        ' '.join([
            getEnergy4Key(energy), 'GeV', '{/Symbol \264} 10^{%d}' % (
                Decimal(yscale[energy]).as_tuple().exponent
                )
            ]) if data_type == 'data' else ''
    )
  # order dict for panel plot
  dpt_dict_sort = OrderedDict(
    (' '.join([getMeeLabel(k), ':', mee_range[k], ' GeV/c^{2}']), dpt_dict[k])
    for k in mee_keys
  )
  yMin, yMax = 0.5*min(yvals), 3*max(yvals)
  # make panel plot
  make_panel(
    dpt_dict = dpt_dict_sort,
    name = os.path.join(outDir, 'ptspec'),
    ylabel = '1/N@_{mb}^{evt} d^{2}N@_{ee}^{acc.}/p_{T}dp_{T}dM_{ee} (c^4/GeV^3)',
    xlabel = 'dielectron transverse momentum, p_{T} (GeV/c)',
    ylog = True, xr = [0, 1.05], yr = [yMin, yMax],
    lmargin = 0.045, bmargin = 0.15, rmargin = 0.998,
    key = ['bottom left', 'samplen 0.5', 'width -1'],
    arrow_bar = 0.002,
  )
  # make mean pt plot
  yMinPt, yMaxPt = 0.95*min(yvalsPt), 1.25 #1.05*max(yvalsPt)
  make_plot(
    data = [ # cocktail
      np.array(data_avpt[k+'_c']) for k in mee_keys
    ] + [ # data
      np.array(data_avpt[k]) for k in mee_keys
    ],
    properties = [
     'with lines lt 2 lw 4 lc %s' % default_colors[i if i < 5 else i+1]
      for i in xrange(len(mee_keys))
    ] + [
     'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i if i < 5 else i+1]
      for i in xrange(len(mee_keys))
    ],
    titles = [ getMeeLabel(k) for k in mee_keys ] + ['']*len(mee_keys),
    name = os.path.join(outDir, 'meanPt'),
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    ylabel = '{/Symbol \341}p_{T}{/Symbol \361} (GeV/c)',
    lmargin = 0.08, xlog = True,
    xr = [17,220], yr = [yMinPt, yMaxPt],
    key = [ 'maxrows 1', 'at graph 1, 1.1' ],
    gpcalls = [
      'format x "%g"',
      'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
    ],
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_ptspec()
