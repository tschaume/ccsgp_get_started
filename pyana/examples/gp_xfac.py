import logging, argparse, os, sys
import numpy as np
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs, checkSymLink
from ..ccsgp.utils import getOpts

shift = {
  'STAR Au+Au : 0.3 - 0.75 GeV': 1.02, 'STAR Au+Au : 0.2 - 0.6 GeV': 1.04
}

def gp_xfac():
  """example using QM12 enhancement factors

  - uses `gpcalls` kwarg to reset xtics
  - numpy.loadtxt needs reshaping for input files w/ only one datapoint
  - according poster presentations see QM12_ & NSD_ review

  .. _QM12: http://indico.cern.ch/getFile.py/access?contribId=268&sessionId=10&resId=0&materialId=slides&confId=181055
  .. _NSD: http://rnc.lbl.gov/~xdong/RNC/DirectorReview2012/posters/Huck.pdf

  .. image:: ../ccsgp_get_started_data/examples/gp_xfac/xfac.png
     :width: 450 px

  :ivar key: translates filename into legend/key label
  :ivar shift: slightly shift selected data points
  """
  # prepare data
  inDir, outDir = getWorkDirs()
  data = OrderedDict()
  for file in os.listdir(inDir):
    info = os.path.splitext(file)[0].split('_')
    key = ' '.join(info[:2] + [':',
      ' - '.join([
        str(float(s)/1e3) for s in info[-1][:7].split('-')
      ]) + ' GeV'
    ])
    file_url = os.path.join(inDir, file)
    data[key] = np.loadtxt(open(file_url, 'rb')).reshape((-1,5))
    data[key][:, 0] *= shift.get(key, 1)
  logging.debug(data) # shown if --log flag given on command line
  # generate plot
  nSets = len(data)
  make_plot(
    data = data.values(),
    properties = [ getOpts(i) for i in xrange(nSets) ],
    titles = data.keys(), # use data keys as legend titles
    name = os.path.join(outDir, 'xfac'),
    key = [ 'top center', 'maxcols 2', 'width -7', 'font ",20"' ],
    ylabel = 'LMR Enhancement Factor',
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    xlog = True, lmargin = 4.3, gpcalls = [ 'format x "%g"',
      'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)'
    ],
    labels = { 'STAR Preliminary': [0.5, 0.5, False] }
  )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  #parser.add_argument("initial", help="country initial = input subdir with txt files")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_xfac()
