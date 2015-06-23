import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..examples.utils import getWorkDirs
from ..ccsgp.utils import getOpts
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
  """
  # prepare input/output directories
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, guest)
  if not os.path.exists(inDir): # catch missing guest
    return "guest %s doesn't exist" % guest
  # prepare data
  data = OrderedDict()
  for file in os.listdir(inDir):
    funct = os.path.splitext(file)[0]
    file_url = os.path.join(inDir, file)
    data[funct] = np.loadtxt(open(file_url, 'rb')) # load data
  logging.debug(data) # shown if --log flag given on command line
  # generate plot using ccsgp.make_plot
  nSets = len(data)
  make_plot(
    data = data.values(),
    properties = [ 'with boxes lc %s' % (default_colors[i]) for i in
                                         range(nSets)],
    gpcalls = [ 'boxwidth 0.2 absolute', 'style fill solid 1.0 border lt 0',
               'xtics ("Mg"1, "Mn"2, "Fe"3, "Co"4, "Ni"5, "Cu"6, "Zn"7)'],
    titles = data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, guest),
    key = [ 'at graph 1., 1.2', 'maxrows 2' ],
    ylabel = 'binding energy (kJ/mol)',
    xlabel = 'metals', rmargin = 0.99, tmargin = 0.85, size='8.5in,8in',
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
