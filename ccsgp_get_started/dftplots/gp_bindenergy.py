import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..examples.utils import getWorkDirs
from ..ccsgp.utils import getOpts

def gp_bindenergy(guest):
  """example for plotting from a text file via numpy.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. sort countries from highest to lowest population
  4. select the <topN> most populated countries
  5. call ccsgp.make_plot with data from 4

  Below is an output image for country initial T and the 4 most populated
  countries for this initial (click to enlarge). Also see::

    $ python -m ccsgp_get_started.examples.gp_datdir -h

  for help on the command line options.

  .. image:: pics/T.png
     :width: 450 px

  .. image:: pics/U.png
     :width: 450 px

  :param initial: country initial
  :type initial: str
  :param topN: number of most populated countries to plot
  :type topN: int
  :ivar inDir: input directory according to package structure and initial
  :ivar outDir: output directory according to package structure
  :ivar data: OrderedDict with datasets to plot as separate keys
  :ivar file: data input file for specific country, format: [x y] OR [x y dx dy]
  :ivar country: country, filename stem of input file
  :ivar file_url: absolute url to input file
  :ivar nSets: number of datasets
  """
  # prepare input/output directories
  inDir, outDir = getWorkDirs()
  guest = guest.capitalize()
  inDir = os.path.join(inDir, guest)
  if not os.path.exists(inDir): # catch missing initial
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
    properties = [ getOpts(i) for i in xrange(nSets) ],
    titles = data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, guest),
    key = [ 'at graph 1., 1.2', 'maxrows 2' ],
    ylabel = 'binding energy (kJ/mol)',
    xlabel = 'metals', rmargin = 0.99, tmargin = 0.85, size='8.5in,8in'
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
