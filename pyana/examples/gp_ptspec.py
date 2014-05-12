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
  fenergies = ['19', '27', '39', '62', '200']
  nen = len(fenergies)
  #mee_keys = ['pi0', 'LMR', 'omega', 'phi', 'IMR', 'jpsi']
  mee_keys = ['LMR', ]
  mee_dict = OrderedDict((k,'') for k in mee_keys)
  yscale = { '200': '300', '62': '50', '39': '3', '27': '0.15', '19': '0.01' }
  inDir, outDir = getWorkDirs()
  data, data_avpt, dpt_dict = {}, {}, {}
  yvals, yvalsPt = [], []
  scale = {
    '19': 0.360245401469014, '200': 1.0, '39': 0.37353401814177617,
    '27': 0.39374082474968564, '62': 0.47675708579800646
  }
  for filename in os.listdir(inDir):
    # import data
    file_url = os.path.join(inDir, filename)
    filebase = os.path.splitext(filename)[0] # unique
    energy, mee_name, mee_range, data_type = splitFileName(filebase)
    if mee_name not in mee_keys: continue
    mee_dict[mee_name] = mee_range
    data[filebase] = np.loadtxt(open(file_url, 'rb'))
    if energy == '200': data[filebase][:,(1,3,4)] /= 0.5
    if data_type == 'data': data[filebase][:,(1,3,4)] *= scale[energy]
    # calculate average pT first
    probs = unp.uarray(data[filebase][:,1], data[filebase][:,3]) # dN/pT
    probs /= umath.fsum(probs) # probabilities
    pTs = unp.uarray(data[filebase][:,0], data[filebase][:,2]) # pT
    avpt = umath.fsum(pTs*probs)
    weights = data[filebase][:,1] * data[filebase][:,0]
    logging.info(('%s: {} %g' % (
      filebase, np.average(data[filebase][:,0], weights=weights)
    )).format(avpt)) # TODO: syst. uncertainties
    # save datapoint for average pT and append to yvalsPt for yaxis range
    dp = [ float(getEnergy4Key(energy)), avpt.nominal_value, 0., avpt.std_dev, 0. ]
    avpt_key = mee_name
    if data_type == 'cocktail': avpt_key += '_c'
    if data_type == 'medium': avpt_key += '_m'
    if avpt_key in data_avpt: data_avpt[avpt_key].append(dp)
    else: data_avpt[avpt_key] = [ dp ]
    yvalsPt.append(avpt.nominal_value)
    # now adjust data for panel plot and append to yvals
    data[filebase][:,(1,3,4)] *= float(yscale[energy])
    if data_type == 'cocktail' or data_type == 'medium':
        data[filebase][:,2:] = 0.
    yvals += [v for v in data[filebase][:,1] if v > 0]
    # prepare dict for panel plot
    dpt_dict_key = getSubplotTitle(mee_name, mee_range)
    if dpt_dict_key not in dpt_dict:
        ndsets = nen*2
        if mee_name == 'LMR': ndsets += 3 # TODO: currently only 19/39/62 medium avail.
        dpt_dict[dpt_dict_key] = [ [None]*ndsets, [None]*ndsets, [None]*ndsets ]
    enidx = fenergies.index(energy)
    dsidx = enidx
    if data_type == 'medium':
      dsidx = (energy=='19')*0 + (energy=='39')*1 + (energy=='62')*2
    else:
      dsidx += int(mee_name == 'LMR') * 3 # number of medium calc avail.
    dsidx += int(data_type == 'data') * len(fenergies)
    dpt_dict[dpt_dict_key][0][dsidx] = data[filebase] # data
    if data_type == 'data': # properties
      dpt_dict[dpt_dict_key][1][dsidx] = 'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[enidx]
    elif data_type == 'cocktail':
      dpt_dict[dpt_dict_key][1][dsidx] = 'with lines lt 1 lw 4 lc %s' % default_colors[enidx]
    else:
      dpt_dict[dpt_dict_key][1][dsidx] = 'with lines lt 2 lw 4 lc %s' % default_colors[enidx]
    dpt_dict[dpt_dict_key][2][dsidx] = ' '.join([ # legend titles
        getEnergy4Key(energy), 'GeV', '{/Symbol \264} %g' % (
          Decimal(yscale[energy])#.as_tuple().exponent
        )
      ]) if data_type == 'data' else ''
  # use mass range in dict key to sort dpt_dict with increasing mass
  plot_key_order = dpt_dict.keys()
  plot_key_order.sort(key=lambda x: float(x.split(':')[1].split('-')[0]))
  # sort data_avpt by energy and apply x-shift for better visibility
  for k in data_avpt: data_avpt[k].sort(key=lambda x: x[0])
  energies = [ dp[0] for dp in data_avpt[mee_keys[0]] ]
  energies.append(215.) # TODO: think of better upper limit
  #linsp = {}
  #for start,stop in zip(energies[:-1],energies[1:]):
  #  linsp[start] = np.linspace(start, stop, num = 4*len(mee_keys))
  #for k in data_avpt:
  #  key = k.split('_')[0]
  #  for i in xrange(len(data_avpt[k])):
  #    data_avpt[k][i][0] = linsp[energies[i]][mee_keys.index(key)]
  # make panel plot
  yMin, yMax = 0.5*min(yvals), 3*max(yvals)
  #make_panel(
  #  dpt_dict = OrderedDict((k,dpt_dict[k]) for k in plot_key_order),
  #  name = os.path.join(outDir, 'ptspec'),
  #  ylabel = '1/N@_{mb}^{evt} d^{2}N@_{ee}^{acc.}/dp_{T}dM_{ee} (c^3/GeV^2)',
  #  xlabel = 'dielectron transverse momentum, p_{T} (GeV/c)',
  #  ylog = True, xr = [0, 2.0], yr = [1e-6, 5e9],
  #  lmargin = 0.12, bmargin = 0.11, rmargin = 0.998,
  #  key = ['bottom left', 'samplen 0.5', 'width -1', 'opaque'],
  #  arrow_bar = 0.002, layout = '3x2'
  #)
  #make plot for LMR spectra only
  lmr_key = getSubplotTitle('LMR', '0.3-0.76')
  make_plot(
    data = dpt_dict[lmr_key][0],
    properties = dpt_dict[lmr_key][1],
    titles = dpt_dict[lmr_key][2],
    name = os.path.join(outDir, 'ptspecLMR'),
    ylabel = '1/N@_{mb}^{evt} d^{2}N@_{ee}^{acc.}/dp_{T}dM_{ee} (c^3/GeV^2)',
    xlabel = 'dielectron transverse momentum, p_{T} (GeV/c)',
    ylog = True, xr = [0, 2.0], yr = [2e-7, 5],
    lmargin = 0.14, bmargin = 0.08, rmargin = 0.98, tmargin = 0.9,
    key = ['maxrows 2', 'samplen 0.5', 'width -3', 'at graph 1.,1.1'],
    arrow_bar = 0.005, size = '10in,13in'
  )
  # make mean pt plot
  yMinPt, yMaxPt = 0.95*min(yvalsPt), 1.05*max(yvalsPt)
  #make_plot(
  #  data = [ # cocktail
  #    np.array(data_avpt[k+'_c']) for k in mee_keys
  #  ] + [ # medium
  #    np.array(data_avpt['LMR_m'])
  #  ] + [ # data
  #    np.array(data_avpt[k]) for k in mee_keys
  #  ],
  #  properties = [
  #   'with lines lt 1 lw 4 lc %s' % default_colors[i if i < 5 else i+1]
  #    for i in xrange(len(mee_keys))
  #  ] + [
  #   'with lines lt 2 lw 4 lc %s' % default_colors[mee_keys.index('LMR')]
  #  ] + [
  #   'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i if i < 5 else i+1]
  #    for i in xrange(len(mee_keys))
  #  ],
  #  titles = [ getMeeLabel(k) for k in mee_keys ] + ['']*(len(mee_keys)+1),
  #  name = os.path.join(outDir, 'meanPt'),
  #  xlabel = '{/Symbol \326}s_{NN} (GeV)',
  #  ylabel = '{/Symbol \341}p_{T}{/Symbol \361} in STAR Acceptance (GeV/c)',
  #  lmargin = 0.08, xlog = True, xr = [17,220],
  #  yr = [0,1.4], #yr = [yMinPt, yMaxPt],
  #  key = [ 'maxrows 1', 'at graph 1, 1.1' ],
  #  gpcalls = [
  #    'format x "%g"',
  #    'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
  #  ]
  #)
  # make mean pt plot for LMR only
  make_plot(
    data = [
      np.array(data_avpt['LMR_c']),
      np.array(data_avpt['LMR_m']),
      np.array(data_avpt['LMR'])
    ],
    properties = [
     'with lines lt 1 lw 4 lc %s' % default_colors[0],
     'with lines lt 2 lw 4 lc %s' % default_colors[0],
     'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
    ],
    titles = [ 'cocktail', '+medium', getMeeLabel('data') ],
    name = os.path.join(outDir, 'meanPtLMR'),
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    ylabel = 'LMR {/Symbol \341}p_{T}{/Symbol \361} in STAR Acceptance (GeV/c)',
    lmargin = 0.15, bmargin = 0.15, xlog = True, xr = [17,220],
    yr = [0.4,1.], #yr = [yMinPt, yMaxPt],
    key = [ 'bottom right' ],
    gpcalls = [
      'format x "%g"',
      'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
    ]
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
