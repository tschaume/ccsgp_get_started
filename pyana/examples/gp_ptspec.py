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

def getSubplotTitle(mn, mr):
  return ' '.join([getMeeLabel(mn), ':', mr, ' GeV/c^{2}'])

def gp_ptspec():
  """example for a 2D-panel plot (TODO)"""
  fenergies = ['19', '27', '39', '62']
  nen = len(fenergies)
  mee_keys = ['pi0', 'LMR', 'omega', 'phi', 'IMR', 'jpsi']
  mee_dict = OrderedDict((k,'') for k in mee_keys)
  yscale = { '62': '1e6', '39': '1e4', '27': '1e2', '19': '1.' }
  inDir, outDir = getWorkDirs()
  data, data_avpt, dpt_dict = {}, {}, {}
  yvals, yvalsPt = [], []
  for filename in os.listdir(inDir):
    # import data
    file_url = os.path.join(inDir, filename)
    filebase = os.path.splitext(filename)[0] # unique
    energy, mee_name, mee_range, data_type = splitFileName(filebase)
    if energy == '200': continue
    if mee_name not in mee_keys: continue
    mee_dict[mee_name] = mee_range
    data[filebase] = np.loadtxt(open(file_url, 'rb'))
    # calculate average pT first
    probs = unp.uarray(data[filebase][:,1], data[filebase][:,3]) # dN/pTdpT
    probs *= data[filebase][:,0] # = dN/pTdpT * pT
    probs /= umath.fsum(probs) # probabilities
    pTs = unp.uarray(data[filebase][:,0], data[filebase][:,2]) # pT
    avpt = umath.fsum(pTs*probs)
    weights = data[filebase][:,1] * data[filebase][:,0]
    logging.info(('%s: {} %g' % (
      filebase, np.average(data[filebase][:,0], weights=weights)
    )).format(avpt)) # TODO: syst. uncertainties
    # save datapoint for average pT and append to yvalsPt for yaxis range
    dp = [ float(getEnergy4Key(energy)), avpt.nominal_value, 0., avpt.std_dev, 0. ]
    avpt_key = mee_name if data_type == 'data' else mee_name + '_c'
    if avpt_key in data_avpt: data_avpt[avpt_key].append(dp)
    else: data_avpt[avpt_key] = [ dp ]
    yvalsPt.append(avpt.nominal_value)
    # now adjust data for panel plot and append to yvals
    data[filebase][:,(1,3,4)] *= float(yscale[energy])
    if data_type == 'cocktail': data[filebase][:,2:] = 0.
    yvals += [v for v in data[filebase][:,1] if v > 0]
    # prepare dict for panel plot
    dpt_dict_key = getSubplotTitle(mee_name, mee_range)
    if dpt_dict_key not in dpt_dict:
        dpt_dict[dpt_dict_key] = [ [None]*(nen*2), [None]*(nen*2), [None]*(nen*2) ]
    enidx = fenergies.index(energy)
    dsidx = (data_type == 'data') * nen + enidx
    dpt_dict[dpt_dict_key][0][dsidx] = data[filebase] # data
    if data_type == 'data': # properties
      dpt_dict[dpt_dict_key][1][dsidx] = 'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[enidx]
    else:
      dpt_dict[dpt_dict_key][1][dsidx] = 'with lines lt 1 lw 4 lc %s' % default_colors[enidx]
    dpt_dict[dpt_dict_key][2][dsidx] = ' '.join([ # legend titles
        getEnergy4Key(energy), 'GeV', '{/Symbol \264} 10^{%d}' % (
          Decimal(yscale[energy]).as_tuple().exponent
        )
      ]) if data_type == 'data' else ''
  # use mass range in dict key to sort dpt_dict with increasing mass
  plot_key_order = dpt_dict.keys()
  plot_key_order.sort(key=lambda x: float(x.split(':')[1].split('-')[0]))
  # sort data_avpt by energy and apply x-shift for better visibility
  for k in data_avpt: data_avpt[k].sort(key=lambda x: x[0])
  energies = [ dp[0] for dp in data_avpt[mee_keys[0]] ]
  energies.append(75.) # TODO: think of better upper limit
  linsp = {}
  for start,stop in zip(energies[:-1],energies[1:]):
    linsp[start] = np.linspace(start, stop, num = 4*len(mee_keys))
  for k in data_avpt:
    key = k.split('_')[0]
    for i in xrange(len(data_avpt[k])):
      data_avpt[k][i][0] = linsp[energies[i]][mee_keys.index(key)]
  # make panel plot
  yMin, yMax = 0.5*min(yvals), 3*max(yvals)
  make_panel(
    dpt_dict = OrderedDict((k,dpt_dict[k]) for k in plot_key_order),
    name = os.path.join(outDir, 'ptspec'),
    ylabel = '1/N@_{mb}^{evt} d^{2}N@_{ee}^{acc.}/p_{T}dp_{T}dM_{ee} (c^4/GeV^3)',
    xlabel = 'dielectron transverse momentum, p_{T} (GeV/c)',
    ylog = True, xr = [0, 1.05], yr = [yMin, yMax],
    lmargin = 0.045, bmargin = 0.15, rmargin = 0.998,
    key = ['bottom left', 'samplen 0.5', 'width -1'],
    arrow_bar = 0.002,
  )
  # make mean pt plot
  yMinPt, yMaxPt = 0.95*min(yvalsPt), 1.05*max(yvalsPt)
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
    lmargin = 0.08, xlog = True, #xr = [17,220],
    yr = [0,1.4], #yr = [yMinPt, yMaxPt],
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
