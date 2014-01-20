import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.utils import getOpts

def gp_datdir(initial, topN):
  """example for plotting from a text file via numpy.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. sort countries from highest to lowest population
  4. select the <topN> most populated countries
  5. call ccsgp.make_plot with data from 4

  Below is an output image for country initial T and the 4 most populated
  countries for this initial (click to enlarge). Also see::

    $ python -m pyana.examples.gp_datdir -h

  for help on the command line options.

  .. image:: ../ccsgp_get_started_data/examples/gp_datdir/T.png
     :width: 450 px

  .. image:: ../ccsgp_get_started_data/examples/gp_datdir/U.png
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
  initial = initial.capitalize()
  inDir = os.path.join(inDir, initial)
  if not os.path.exists(inDir): # catch missing initial
    return "initial %s doesn't exist" % initial
  # prepare data
  data = OrderedDict()
  for file in os.listdir(inDir):
    country = os.path.splitext(file)[0]
    file_url = os.path.join(inDir, file)
    data[country] = np.loadtxt(open(file_url, 'rb')) # load data
    # set y-axis unit to 1M
    data[country][:, 1] /= 1e6
    if data[country].shape[1] > 2: data[country][:, 3:] /= 1e6
  logging.debug(data) # shown if --log flag given on command line
  # sort countries according to mean population (highest -> lowest)
  sorted_data = OrderedDict(sorted(
    data.items(), key = lambda t: np.mean(t[1][:,1]), reverse = True
  ))
  # "pop" (select) N most populated countries
  top_data = OrderedDict(
    sorted_data.popitem(last = False) for i in xrange(topN)
    if sorted_data
  )
  # generate plot using ccsgp.make_plot
  nSets = len(top_data)
  make_plot(
    data = top_data.values(),
    properties = [ getOpts(i) for i in xrange(nSets) ],
    titles = top_data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, initial),
    key = [ 'at graph 1., 1.2', 'maxrows 2' ],
    ylabel = 'total population ({/Symbol \664} 10^{6})',
    xlabel = 'year', lmargin = 0.08, tmargin = 0.9
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("initial", help="country initial = input subdir with txt files")
  parser.add_argument("topN", help="number of most populated countries to plot")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_datdir(args.initial, int(args.topN))
