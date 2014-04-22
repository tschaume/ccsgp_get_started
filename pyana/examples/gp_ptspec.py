import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_panel, make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors
from decimal import Decimal

mee_keys = ['pi0', 'LMR', 'omega', 'phi', 'IMR', 'jpsi']

def getMeeLabel(s):
  if s == 'pi0': return '{/Symbol \160}^0'
  if s == 'omega': return '{/Symbol \167}'
  if s == 'phi': return '{/Symbol \152}'
  if s == 'jpsi': return 'J/{/Symbol \171}'
  return s

def gp_ptspec():
  """example for a 2D-panel plot"""
  # import data and calculate average pT's
  inDir, outDir = getWorkDirs()
  data, avpt = OrderedDict(), OrderedDict()
  for file in os.listdir(inDir):
    file_url = os.path.join(inDir, file)
    key = os.path.splitext(file)[0] # unique
    data[key] = np.loadtxt(open(file_url, 'rb'))
    avpt[key] = np.average(data[key][:,0], weights=data[key][:,0]*data[key][:,1])
    # TODO: statistical errors on avpt!
  # generate input dpt_dict
  dpt_dict, avpt_data = OrderedDict(), OrderedDict()
  mee_range = {}
  yvals, yvalsPt = [], []
  yscale = { '62': '1e6', '39': '1e4', '27': '1e2', '19': '1.' }
  for k in sorted(data.keys()):
    # energy, mee_name, mass range, data type
    base, mee_name, rnge = k.split('_')
    if mee_name not in mee_range: mee_range[mee_name] = rnge
    energy = re.compile('\d+').search(base).group()
    if energy == '200': continue
    data_type = re.sub(str(energy), '', base)
    # pt spectra
    data[k][:,(1,3,4)] *= float(yscale[energy])
    if data_type == 'cocktail': data[k][:,2:] = 0.
    yvals += [v for v in data[k][:,1] if v > 0]
    if mee_name not in dpt_dict: dpt_dict[mee_name] = [ [], [], [] ]
    dpt_dict[mee_name][0].append(data[k])
    dpt_dict[mee_name][1].append(
        'lt 1 lw 4 ps 1.5 lc %s pt 18' % (default_colors[len(dpt_dict[mee_name][1])])
        if data_type == 'data' else
        'with lines lt 1 lw 4 lc %s' % (default_colors[len(dpt_dict[mee_name][1])])
        )
    dpt_dict[mee_name][2].append(
        ' '.join([
            getEnergy4Key(energy), 'GeV', '{/Symbol \264} 10^{%d}' % (
                Decimal(yscale[energy]).as_tuple().exponent
                )
            ]) if data_type == 'data' else ''
    )
    # mean pt
    if data_type == 'cocktail':
        yvalsPt.append(avpt[k])
        dp = [ float(energy), avpt[k], 0., 0., 0. ] # TODO: stat./syst. uncertainties
        if mee_name in avpt_data: avpt_data[mee_name].append(dp)
        else: avpt_data[mee_name] = [ dp ]
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
    data = [ np.array(avpt_data[k]) for k in mee_keys ],
    properties = [
        'with linespoints lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
        for i in xrange(len(mee_keys))
        ],
    titles = [ getMeeLabel(k) for k in mee_keys ],
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
