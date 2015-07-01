import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..examples.utils import getWorkDirs
from ..ccsgp.utils import getOpts, colorscale
from ..ccsgp.config import default_colors

def gp_bindenergy(guest):
  """example for plotting from a text file via numpy.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. call ccsgp.make_plot with data from 2

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
  :ivar nSets: number of datasets
  :ivar nfiles: number of files
  :ivar dx: bar width estimated for specific amount of datasets
  :ivar gap: space between 2 bars
  """
  # prepare color arrays for dataset
  my_color_set = [default_colors[i] for i in range(0, 4)]
  my_color_array = []
  for color in my_color_set:
      my_color_array.append(colorscale(color[-7:-1], 0.8))
      my_color_array.append(colorscale(color[-7:-1], 1.4))
  # prepare input/output directories
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, guest)
  if not os.path.exists(inDir): # catch missing guest
    return "guest %s doesn't exist" % guest
  # prepare data
  data = OrderedDict()
  files = os.listdir(inDir)
  nfiles = len(files)-1
  dx = 0.7/nfiles
  gap = (1. - dx*nfiles)/2
  for idx,file in enumerate(files):
    funct, ext = os.path.splitext(file)
    if funct == 'experiments' or ext != '.dat': continue
    file_url = os.path.join(inDir, file)
    data[funct] = np.loadtxt(open(file_url, 'rb')) # load data
    data[funct][:,0] += (idx - (nfiles+1)/2. + 0.5*(not nfiles%2)) * dx
  data['M06-L'] = data.pop('M06-L')
  exp_data = np.loadtxt(open(os.path.join(inDir,'experiments.dat'), 'rb')) # load exp_data
  logging.debug(data) # shown if --log flag given on command line
  # generate plot using ccsgp.make_plot
  nSets = len(data)
  make_plot(
    data = [np.array([[1,-5]])] + data.values(),
    properties = ['with boxes lt -1 lc %s' % default_colors[-10]] + [
        'with boxes lt -1 lc %s' % (my_color_array[i]) for i in range(nSets)
    ], gpcalls = [
        'boxwidth {} absolute'.format(dx),
        'style fill solid 1.0 border lt -1',
        'xtics ("Mg"1, "Mn"2, "Fe"3, "Co"4, "Ni"5, "Cu"6, "Zn"7)',
    ] + [
        'arrow {} from {},{} to {},{} lw 4 lc {} nohead front'.format(
            i+5, dp[0]-nfiles*dx/2., dp[1], dp[0]+nfiles*dx/2., dp[1], default_colors[-10]
        ) for i,dp in enumerate(exp_data)
    ], titles = ['experiments'] + data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, guest), yreverse = True,
    key = [ 'at graph 1.07, 1.25', 'maxrows 3', 'width -1.1', 'nobox' ],
    ylabel = 'electronic binding energy (kJ/mol)',
    xlabel = 'metals', rmargin = 0.99, tmargin = 0.83, size='8.5in,8in',
    debug = True
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
