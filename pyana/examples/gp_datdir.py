import logging, os
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..aux.utils import getWorkDir, zip_flat, getOpts
from ..aux.config import defaultkey

def gp_datdir(subname):
  inDir = os.path.join(getWorkDir(__name__, True), subname)
  return 'inDir'
  # TODO: prepare data = OrderedDict
  outDir = getWorkDir(__name__)
  logging.debug(data)
  nSets = len(data)
  make_plot(
    name = os.path.join(outDir, 'gp_from_txt'), # basename of output file
    log = [False, True], # logarithmic axes? [x, y]
    xr = [xMin, xMax], yr = [yMin, yMax],
    data = data,
    using = ['1:2:3'] * nSets,
    main_opts = ['points'] * nSets, # yerrorbars, boxerrorbars
    extra_opts = zip_flat(
      [ 'lt 1 lw 4 ps 2 lc %d pt %d' % (i/2, 18+i%2) for i in xrange(nSets) ],
      [ "lt 1 lw 4 lc rgb 'grey'" ] * nSets
    ),
    xlabel = 'xlabel', ylabel = 'ylabel',
    titles = data.keys(),
    key = defaultkey + [
      'maxrows 3', 'top', 'outside', 'opaque', 'samplen 0.8',
      'font ",22"', 'width -2', 'height 0.2'
    ]
  )
  return 'done'
