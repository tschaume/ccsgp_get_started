import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np

def gp_sims(version):
  """example for a batch generating simple plots (cocktail contributions)

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version, 'cocktail_contribs')
  energies = [19, 27, 39, 62]
  xmax = {
      'pion': 0.125, 'eta': 0.52, 'etap': 0.92, 'omega': 1.22,
      'phi': 1.22, 'jpsi': 3.52
  }
  for particles in [
      'pion', 'eta', 'etap', ['omega', 'rho'], 'phi', 'jpsi'
  ]:
      if isinstance(particles, str):
          particles = [particles]
      contribs = OrderedDict()
      for particle in particles:
          for energy in energies:
              fstem = particle+str(energy)
              fname = os.path.join(inDir, fstem+'.dat')
              contribs[fstem] = np.loadtxt(open(fname, 'rb'))
              contribs[fstem][:,2:] = 0
      print contribs.keys()
      titles = [
          ' '.join([getEnergy4Key(str(energy)), 'GeV'])
          for energy in energies
      ] + [ '' for k in range(len(contribs)-len(energies)) ]
      make_plot(
          data = contribs.values(),
          properties = [
              'with lines lc %s lw 4 lt %d' % (
                  default_colors[i%len(energies)], i/len(energies)+1
              ) for i in xrange(len(contribs))
          ],
          titles = titles,
          xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})' if particles[0] == 'phi' else '',
          ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]' \
          if particles[0] == 'pion' or particles[0] == 'omega' else '',
          name = os.path.join(outDir, '_'.join(['sims']+particles)),
          ylog = True, lmargin = 0.18, bmargin = 0.13, tmargin = 0.96,
          gpcalls = [ 'nokey' ] if particles[0] != 'pion' else [],
          xr = [1. if particles[0] == 'jpsi' else 0,xmax[particles[0]]],
          size = '8.5in,8in',
          labels = {
              particleLabel4Key(particles[0]): [0.15,0.9,False],
              particleLabel4Key(particles[1]) if len(particles) > 1 else '': [0.15,0.1,False],
          }
      )

def gp_sims_panel(version):
  """panel plot of cocktail simulations at all energies, includ. total

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  energies = [19, 27, 39, 62]
  mesons = ['pion', 'eta', 'etap', 'rho', 'omega', 'phi', 'jpsi']
  fstems = ['cocktail_contribs/'+m for m in mesons] + ['cocktail']
  data = OrderedDict((energy, [
      np.loadtxt(open(os.path.join(inDir, fstem+str(energy)+'.dat'), 'rb'))
      for fstem in fstems
  ]) for energy in energies)
  for v in data.values():
      for d in v:
          d[:,2:] = 0
  make_panel(
      dpt_dict = OrderedDict((
          ' '.join([getEnergy4Key(str(energy)), 'GeV']),
          [ data[energy], [
              'with lines lc %s lw 4 lt %d' % (
                  default_colors[(-2*i-2) if i!=len(fstems)-1 else 0], int(i==3)+1
              ) for i in xrange(len(fstems))
          ], [particleLabel4Key(m) for m in mesons] + ['Cocktail (w/o {/Symbol \162})'] ]
      ) for energy in energies),
      name = os.path.join(outDir, 'sims_panel'),
      ylog = True, xr = [0.,3.2], yr = [1e-6,9],
      ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]',
      xlabel = 'invariant dielectron mass, M_{ee} (GeV/c^{2})',
      layout = '2x2', size = '7in,9in', key = ['width -4']
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
  #gp_sims(args.version)
  gp_sims_panel(args.version)
