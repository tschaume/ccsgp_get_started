import os
import numpy as np
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs
from ..ccsgp.config import default_colors

def gp_lcltpt():
  """ example plot to display linecolors, linetypes and pointtypes

  .. image:: pics/gp_lcltpt.png
     :width: 450 px
  """
  inDir, outDir = getWorkDirs()
  nSets = len(default_colors)
  make_plot(
    data = [
      np.array([ [0,i,0,0,0], [1,i,0,0,0] ])
      for i in xrange(nSets)
    ],
    properties = [
      'with linespoints lw 4 lc %s lt %d pt %d' % (col, i, i)
      for i, col in enumerate(default_colors)
    ],
    titles = [''] * nSets, yr = [-1, 51],
    name = os.path.join(outDir, 'gp_lcltpt'),
    ylabel = 'linecolor / linetype / pointtype', xlabel = '',
  )
  return 'done'

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser()
  args = parser.parse_args()
  print gp_lcltpt()
