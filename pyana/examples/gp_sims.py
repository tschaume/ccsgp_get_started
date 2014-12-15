import os, argparse, logging
from .utils import getWorkDirs, checkSymLink
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.config import default_colors
import numpy as np

def gp_sims(version):
  """example for a batch generating simple plots (cocktail contributions)

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version, 'cocktail_contribs')
  xmax = {
      'pion': 0.13, 'eta': 0.55, 'etap': 0.95, 'rho': 1.1,
      'phi': 1.3, 'jpsi': 3.5
  }
  for particles in [
      'pion', 'eta', 'etap', ['rho', 'omega'], 'phi', 'jpsi'
  ]:
      if isinstance(particles, str):
          particles = [particles]
      contribs = OrderedDict()
      for energy in [19, 27, 39, 62]:
          for particle in particles:
              fstem = particle+str(energy)
              fname = os.path.join(inDir, fstem+'.dat')
              contribs[fstem] = np.loadtxt(open(fname, 'rb'))
              contribs[fstem][:,2:] = 0
      make_plot(
          data = contribs.values(),
          properties = [
              'with lines lc %s lw 4 lt 2' % default_colors[i]
              for i in xrange(len(contribs))
          ],
          titles = contribs.keys(),
          name = os.path.join(outDir, '_'.join(['sims']+particles)),
          ylog = True, lmargin = 0.1, bmargin = 0.15,
          gpcalls = [ 'nokey' ] if particles[0] != 'pion' else [],
          xr = [0,xmax[particles[0]]]
      )

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  gp_sims(args.version)
