import logging, os
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..aux.utils import getWorkDir, zip_flat, getOpts
from ..aux.config import defaultkey

def gp_datdir(subname): # subname = country initial in examples case
  inDir = os.path.join(getWorkDir(__name__, True), subname) # inDir for country initial
  data = OrderedDict() # OrderedDict with datasets to plot as separate keys
  for file in os.listdir(inDir): # loop all files for country initial
    country = os.path.splitext(file)[0] # filename stem is country
    data[country] = np.loadtxt( open(os.path.join(inDir, file), 'rb') ) # load data
    data[country][:, 1:3] /= 1e6 # set unit to 1M
    if len(data) > 10: break
  yVals = [ n for v in data.values() for n in v[:, 1] ] # all y values
  yMin = min(yVals)
  yMax = max(yVals)
  logging.debug(data) # shown if --log flag given on command line
  outname = os.path.join(getWorkDir(__name__), subname) # output filename for make_plot
  nSets = len(data) # number of datasets
  make_plot(
    name = outname, # extension-less name of output file
    log = [False, False], # logarithmic axes? [x, y]
    xr = [1999, 2020], yr = [0.9*yMin, 1.1*yMax],
    data = data.values(),
    using = ['1:2:3'] * nSets, # datapoint format: x y yerr bw
    main_opts = ['points'] * nSets, # other options: yerrorbars, boxerrorbars
    extra_opts = [ getOpts(i) for i in xrange(nSets*2) ], # len(extra_opts) = nSets*2! pts/err alternately
    xlabel = 'year', ylabel = 'total population ({/Symbol \664} 10^{6})',
    titles = data.keys(), # using data keys as legend titles
    key = defaultkey + [ 'font ",22"', 'width -1' ] # add extra options
  )
  return 'done'

#    MORE CONTROL for extra_opts:
#    zip_flat from aux/utils.py takes to lists and zips them flat alternately
#    extra_opts = zip_flat(
#      [ 'lt 1 lw 4 ps 2 lc %d pt %d' % (i/2, 18+i%2) for i in xrange(nSets) ],  # datapoints
#      [ "lt 1 lw 4 lc rgb 'grey'" ] * nSets # errors
#    ),
