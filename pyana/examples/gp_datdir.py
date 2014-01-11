import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..aux.utils import getWorkDirs, zip_flat, getOpts
from ..aux.config import defaultkey

# subname = country initial in examples case
# TODO: write pydoc's
def gp_datdir(subname):
  inDir, outDir = getWorkDirs()
  initial = subname.capitalize() # catch small case initials
  inDir = os.path.join(inDir, initial) # inDir for country initial
  if not os.path.exists(inDir): # catch missing initial
    return "initial %s doesn't exist" % initial
  data = OrderedDict() # OrderedDict with datasets to plot as separate keys
  for file in os.listdir(inDir): # loop all files for country initial
    country = os.path.splitext(file)[0] # filename stem is country
    file_url = os.path.join(inDir, file) # absolute file url
    data[country] = np.loadtxt(open(file_url, 'rb')) # load data
    data[country][:, 1:3] /= 1e6 # set unit to 1M
    if len(data) > 10: break # don't plot more than 10 countries
  yVals = [ n for v in data.values() for n in v[:, 1] ] # all y values
  yMin, yMax = min(yVals), max(yVals) # yMin, yMax for axis range
  logging.debug(data) # shown if --log flag given on command line
  nSets = len(data) # number of datasets
  make_plot(
    name = os.path.join(outDir, initial), # extension-less name of output file
    log = [False, False], # logarithmic axes? [x, y]
    xr = [1999, 2020], yr = [0.9*yMin, 1.1*yMax],
    data = data.values(), # list of numpy arrays
    using = ['1:2:3'] * nSets, # datapoint format: x y yerr bw
    main_opts = ['points'] * nSets, # other options: yerrorbars, boxerrorbars
    extra_opts = [ # len(extra_opts) = nSets*2! pts/err alternately
      getOpts(i) for i in xrange(nSets*2) # see below for more control
    ], # using getOpts helper to generate gnuplot option string
    xlabel = 'year', # for symbol numbers see http://bit.ly/1erBgIk
    ylabel = 'total population ({/Symbol \664} 10^{6})',
    titles = data.keys(), # using data keys as legend titles
    key = defaultkey + [ 'font ",22"', 'width -1' ], # add extra options
    #write = True # write data to hdf5 file (optional), install h5py if you like to use it
  )
  return 'done'

#    MORE CONTROL for extra_opts:
#    zip_flat from aux/utils.py takes to lists and zips them flat alternately
#    extra_opts = zip_flat(
#      [ 'lt 1 lw 4 ps 2 lc %d pt %d' % (i/2, 18+i%2) for i in xrange(nSets) ],  # datapoints
#      [ "lt 1 lw 4 lc rgb 'grey'" ] * nSets # errors
#    ),

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("subname", help="input subdir with txt files")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_datdir(args.subname)
