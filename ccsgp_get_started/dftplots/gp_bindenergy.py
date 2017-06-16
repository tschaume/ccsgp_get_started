import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_panel
from ..examples.utils import getWorkDirs
from ..ccsgp.utils import getOpts, colorscale
from ..ccsgp.config import default_colors

def gp_bindenergy(guest):
  """example for plotting from a text file via numpy.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. call ccsgp.make_panel with data from 2

  Also see::

    $ python -m ccsgp_get_started.dftplots.gp_bindenergy -h

  for help on the command line options.

  :param guest: guest molecule
  :type guest: str
  :ivar inDir: input directory according to package structure and guest
  :ivar outDir: output directory according to package structure
  :ivar data: OrderedDict with datasets to plot as separate keys
  :ivar file: data input file for specific guest, format: [x y] OR [x y dx dy]
  :ivar funct: functional, filename stem of input file
  :ivar file_url: absolute url to input file
  :ivar nfiles: number of files
  :ivar dx: bar width estimated for specific amount of datasets
  """
  # prepare color arrays for dataset
  my_color_set = [default_colors[3]] + [default_colors[i] for i in range(0, 3)]
  my_color_array = []
  for icol,color in enumerate(my_color_set):
      my_color_array.append(colorscale(color[-7:-1], 0.8))
      if icol == 0: continue
      my_color_array.append(colorscale(color[-7:-1], 1.4))
  # prepare input/output directories
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, guest)
  if not os.path.exists(inDir): # catch missing guest
    return "guest %s doesn't exist" % guest
  # prepare data
  dpt_dict = OrderedDict()
  if 'CO2' in guest: key = '(a)'
  elif 'CH4' in guest: key = '(b)'
  else: key = '(a)'
  dpt_dict[key] = [[], [], []] #dpt
  # exp-data
  exp_filename = os.path.join(inDir,'experiments.dat')
  use_exp_data = os.path.exists(exp_filename)
  isM = guest.startswith('M-')
  if use_exp_data:
    exp_data = np.loadtxt(open(exp_filename, 'rb')) # load exp_data
    if len(exp_data.shape) < 2:
        exp_data = np.array([exp_data])
    if isM: exp_data[:,1] /= 10.
  yfake = -5 if use_exp_data else -40
  if isM: yfake = 1
  dpt_dict[key][0].append(np.array([[1,yfake]]))
  dpt_dict[key][1].append(
    'with boxes lt -1 lc {}'.format(default_colors[-10]) if use_exp_data else \
    'with lines lc rgb "white"'
  )
  dpt_dict[key][2].append('experiments' if use_exp_data else ' ')
  # sim-data
  files = [
      os.path.join(inDir, funct+ext)
      for funct,ext in map(os.path.splitext, os.listdir(inDir))
      if funct != 'experiments' and ext == '.dat'
  ]
  nfiles = len(files)
  dx = 0.7/nfiles
  for idx,file_url in enumerate(files):
      funct = os.path.splitext(os.path.basename(file_url))[0]
      data_import = np.loadtxt(open(file_url, 'rb')) # load data
      data_import[:,0] += (idx - (nfiles-1)/2.) * dx
      if isM: data_import[:,1] /= 10.
      dpt_dict[key][0].append(data_import)
      dpt_dict[key][1].append('with boxes lt -1 lc %s' % my_color_array[idx])
      dpt_dict[key][2].append(funct)
  # relative differences
  if guest != 'M-O_H2O_relax' and use_exp_data:
    dpt_dict[' '] = [[], [], []] #dpt
    for idx,d in enumerate(dpt_dict[key][0][1:]):
        if not isM and idx == 4: continue
        try:
          exp_data_mod = exp_data if idx != 0 else exp_data[(0,6),:]
        except:
          exp_data_mod = np.array([exp_data[0,:]])
        diff = np.copy(d)
        if len(diff) != len(exp_data_mod):
          for i,dp in enumerate(diff):
            if i >= len(exp_data_mod) or round(dp[0]) != exp_data_mod[i,0]:
              exp_data_mod = np.insert(exp_data_mod, i, [round(dp[0]), dp[1]], 0)
        diff[:,1] /= exp_data_mod[:,1]
        diff[:,1] -= 1.
        dpt_dict[' '][0].append(diff)
        dpt_dict[' '][1].append('with boxes lt -1 lc %s' % my_color_array[idx])
        dpt_dict[' '][2].append(dpt_dict[key][2][idx+1])
  logging.debug(dpt_dict) # shown if --log flag given on command line
  # generate plot using ccsgp.make_panel
  exp_data_gpcalls = [
    'arrow {} from {},{} to {},{} lw 4 lc {} nohead front'.format(
      i+5, dp[0]-nfiles*dx/2., dp[1], dp[0]+nfiles*dx/2., dp[1], default_colors[-10]
    ) for i,dp in enumerate(exp_data)
  ] if use_exp_data else []
  guest_split = guest.split('_')
  ylabel = '{} distance for {} binding (nm)'.format(guest_split[0], guest_split[1]) if isM \
    else 'electronic binding energy (kJ/mol)'
  make_panel(
      dpt_dict = dpt_dict,
      gpcalls = [
          'boxwidth {} absolute'.format(dx),
          'style fill solid 1.0 border lt -1',
          'xtics ("Mg"1, "Mn"2, "Fe"3, "Co"4, "Ni"5, "Cu"6, "Zn"7)',
      ] + exp_data_gpcalls,
      name = os.path.join(outDir, guest), yreverse = (not isM),
      key = [ 'at graph 1.02, 1.17', 'maxrows 2', 'width -1.1', 'nobox' ],
      ylabel = ylabel, xlabel = 'metals', rmargin = 0.99,
      tmargin = 0.93 if guest != 'M-O_H2O_relax' and use_exp_data else 0.88,
      size='7in,5.5in' if guest != 'M-O_H2O_relax' and use_exp_data else '4in,5.5in',
      debug = True, layout = '1x2' if guest != 'M-O_H2O_relax' and use_exp_data else '1x1'
  )
  return 'done'

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("guest", help="guest = input subdir with txt files")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_bindenergy(args.guest)
