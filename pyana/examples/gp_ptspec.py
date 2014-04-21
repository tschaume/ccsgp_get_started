import logging, argparse, os, sys, re
import numpy as np
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_panel
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
  # import data
  inDir, outDir = getWorkDirs()
  data = {}
  for file in os.listdir(inDir):
    file_url = os.path.join(inDir, file)
    key = os.path.splitext(file)[0] # unique
    data[key] = np.loadtxt(open(file_url, 'rb'))
  # generate input dpt_dict
  dpt_dict = OrderedDict()
  mee_range = {}
  yvals = []
  yscale = {
    '62': '1e6', '39': '1e4', '27': '1e2', '19': '1.'
  }
  for k in sorted(data.keys()):
    base, mee_name, rnge = k.split('_')
    if mee_name not in mee_range: mee_range[mee_name] = rnge
    energy = re.compile('\d+').search(base).group()
    data[k][:,(1,3,4)] *= float(yscale[energy])
    yvals += [v for v in data[k][:,1] if v > 0]
    if mee_name not in dpt_dict: dpt_dict[mee_name] = [ [], [], [] ]
    dpt_dict[mee_name][0].append(data[k])
    dpt_dict[mee_name][1].append('lt 1 lw 4 ps 1.5 lc %s pt 18' % (
      default_colors[len(dpt_dict[mee_name][1])]
    ))
    dpt_dict[mee_name][2].append(' '.join([
      getEnergy4Key(energy), 'GeV', '{/Symbol \264} 10^{%d}' % (
        Decimal(yscale[energy]).as_tuple().exponent
      )
    ]))
  dpt_dict_sort = OrderedDict(
    (' '.join([getMeeLabel(k), ':', mee_range[k]]), dpt_dict[k])
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
  return 'done'
  #data_type = re.sub(str(energy), '', base)

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
