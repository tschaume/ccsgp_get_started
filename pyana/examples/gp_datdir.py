import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..aux.utils import getWorkDirs, getOpts

def gp_datdir(initial):
  """example for plotting from a text file using ccsgp

  - this is an easy way to use ccsgp by accepting its defaults
  - see ccsgp documentation for more on make_plot

  :param initial: country initial
  :type initial: str
  :ivar inDir: input directory according to package structure and initial
  :ivar outDir: output directory according to package structure
  :ivar data: OrderedDict with datasets to plot as separate keys
  :ivar file: data input file for specific country, format: x y y_err bin_width
  :ivar country: country, filename stem of input file
  :ivar file_url: absolute url to input file
  :ivar nSets: number of datasets
  """
  inDir, outDir = getWorkDirs()
  initial = initial.capitalize()
  inDir = os.path.join(inDir, initial)
  if not os.path.exists(inDir): # catch missing initial
    return "initial %s doesn't exist" % initial
  data = OrderedDict()
  for file in os.listdir(inDir):
    country = os.path.splitext(file)[0]
    file_url = os.path.join(inDir, file)
    data[country] = np.loadtxt(open(file_url, 'rb')) # load data
    data[country][:, 1:3] /= 1e6 # set unit to 1M
    if len(data) > 10: break # don't plot more than 10 countries
  logging.debug(data) # shown if --log flag given on command line
  nSets = len(data)
  make_plot(
    data = data.values(),
    styles = ['points'] * nSets,
    properties = [ getOpts(i) for i in xrange(nSets) ],
    titles = data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, initial),
    key = [ 'width -1' ],
    xlabel = 'year',
    ylabel = 'total population ({/Symbol \664} 10^{6})'
  )
  return 'done'

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("initial", help="country initial = input subdir with txt files")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_datdir(args.initial)
