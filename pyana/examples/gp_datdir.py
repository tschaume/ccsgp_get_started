import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.utils import getOpts

def gp_datdir(initial, topN):
  """example for plotting from a text file via np.loadtxt

  1. prepare input/output directories
  2. load the data into an OrderedDict() [adjust axes units]
  3. sort countries from highest to lowest population
  4. select the <topN> most populated countries
  5. call ccsgp.make_plot with data from 4

  also see ``$ python -m pyana.examples.gp_datdir -h``.

  :param initial: country initial
  :type initial: str
  :param topN: number of most populated countries to plot
  :type topN: int
  :ivar inDir: input directory according to package structure and initial
  :ivar outDir: output directory according to package structure
  :ivar data: OrderedDict with datasets to plot as separate keys
  :ivar file: data input file for specific country, format: x y y_err bin_width
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
    data[country][:, 1:3] /= 1e6 # set unit to 1M
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
    styles = ['points'] * nSets,
    properties = [ getOpts(i) for i in xrange(nSets) ],
    titles = top_data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, initial),
    key = [ 'at graph 1., 1.2', 'maxrows 2' ],
    ylabel = 'total population ({/Symbol \664} 10^{6})',
    xlabel = 'year', tmargin = 2.0
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
