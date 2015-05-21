import os, argparse, logging
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np

energies = [19, 27, 39, 62]
xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})'
ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]'

def gp_sims(version):
  """example for a batch generating simple plots (cocktail contributions)

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version, 'cocktail_contribs')
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
          xlabel = xlabel, # if particles[0] == 'phi' else '',
          ylabel = ylabel, # if particles[0] == 'pion' or particles[0] == 'omega' else '',
          name = os.path.join(outDir, '_'.join(['sims']+particles)),
          ylog = True, lmargin = 0.15, bmargin = 0.17, tmargin = 0.96, rmargin = 0.98,
          gpcalls = [ 'nokey' ] if particles[0] != 'pion' else [],
          xr = [1. if particles[0] == 'jpsi' else 0,xmax[particles[0]]],
          size = '10in,6.5in',
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
  mesons = ['pion', 'eta', 'etap', 'rho', 'omega', 'phi', 'jpsi']
  fstems = ['cocktail_contribs/'+m for m in mesons] + ['cocktail', 'cocktail_contribs/ccbar']
  data = OrderedDict((energy, [
      np.loadtxt(open(os.path.join(inDir, fstem+str(energy)+'.dat'), 'rb'))
      for fstem in fstems
  ]) for energy in energies)
  for v in data.values():
      # keep syserrs for total cocktail
      v[-2][:,(2,3)] = 0
      # all errors zero for cocktail contribs
      v[-1][:,2:] = 0
      for d in v[:-2]: d[:,2:] = 0
  make_panel(
      dpt_dict = OrderedDict((
          ' '.join([getEnergy4Key(str(energy)), 'GeV']),
          [ data[energy], [
              'with %s lc %s lw 5 lt %d' % (
                  'lines' if i!=len(fstems)-2 else 'filledcurves pt 0',
                  default_colors[(-2*i-2) if i!=len(fstems)-2 else 0] \
                  if i!=len(fstems)-1 else default_colors[1], int(i==3)+1
              ) for i in xrange(len(fstems))
          ], [particleLabel4Key(m) for m in mesons] + [
              'Cocktail (w/o {/Symbol \162})', particleLabel4Key('ccbar')
          ]]
      ) for energy in energies),
      name = os.path.join(outDir, 'sims_panel'),
      ylog = True, xr = [0.,3.2], yr = [1e-6,9],
      xlabel = xlabel, ylabel = ylabel,
      layout = '2x2', size = '7in,9in',
      key = ['width -4', 'spacing 1.5', 'nobox', 'at graph 0.9,0.95']
  )

def gp_sims_total_overlay(version):
  """single plot comparing total cocktails at all energies

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = OrderedDict()
  for energy in energies:
      fname = os.path.join(inDir, 'cocktail'+str(energy)+'.dat')
      data[energy] = np.loadtxt(open(fname, 'rb'))
      data[energy][:,2:] = 0
  make_plot(
      data = data.values(),
      properties = [
          'with lines lc %s lw 4 lt 1' % (default_colors[i])
          for i in xrange(len(energies))
      ],
      titles = [
          ' '.join([getEnergy4Key(str(energy)), 'GeV'])
          for energy in energies
      ],
      xlabel = xlabel, ylabel = ylabel,
      name = os.path.join(outDir, 'sims_total_overlay'),
      ylog = True, xr = [0.,3.2], yr = [1e-6,9],
      tmargin = 0.98, rmargin = 0.99, bmargin = 0.14,
      size = '8.5in,8in',
  )

def gp_sims_totalerrors_overlay(version):
  """single plot comparing syst. uncertainties on total cocktails at all energies

  :param version: plot version / input subdir name
  :type version: str
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data = OrderedDict()
  for energy in energies:
      fname = os.path.join(inDir, 'cocktail'+str(energy)+'.dat')
      data[energy] = np.loadtxt(open(fname, 'rb'))
      data[energy][:,1] = data[energy][:,4]/data[energy][:,1]
      data[energy][:,2:] = 0
  make_plot(
      data = data.values(),
      properties = [
          'with lines lc %s lw 4 lt 1' % (default_colors[i])
          for i in xrange(len(energies))
      ],
      titles = [
          ' '.join([getEnergy4Key(str(energy)), 'GeV'])
          for energy in energies
      ],
      xlabel = xlabel, ylabel = 'total relative systematic uncertainty',
      name = os.path.join(outDir, 'sims_totalerrors_overlay'),
      xr = [0.,3.2], yr = [0.15,0.65], key = ['at graph 0.7,0.3'],
      tmargin = 0.98, rmargin = 0.99, bmargin = 0.14,
      size = '8.5in,8in',
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
  #gp_sims_panel(args.version)
  #gp_sims_total_overlay(args.version)
  gp_sims_totalerrors_overlay(args.version)
