import logging, argparse, os, re
import numpy as np
from ..ccsgp.config import default_colors
from ..ccsgp.ccsgp import make_plot
from .utils import getWorkDirs, checkSymLink, getMassRangesSums
from math import pi, log

def gp_rapp():
  """rho in-medium ratios by Rapp (based on protected data)"""
  inDir, outDir = getWorkDirs()
  # prepare data
  yields = {}
  for infile in os.listdir(inDir):
    energy = re.compile('\d+').search(infile).group()
    medium = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
    getMassRangesSums(energy, medium, yields)
  data = dict( # sort by energy
    (k, np.array(sorted(v)))
    for k, v in yields.iteritems()
  )
  for k in data: data[k][:,1] /= data[k][-1,1] # divide by 200
  # make plot
  nSets = len(data)
  make_plot(
    data = [ data[k][:-1] for k in data ],
    properties = [
      'with linespoints lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
      for i in xrange(nSets)
    ],
    titles = data.keys(), # TODO: titles order correct?
    name = os.path.join(outDir, 'rapp'),
    xlabel = '{/Symbol \326}s_{NN} (GeV)', ylabel = 'Rapp Ratio to 200 GeV',
    lmargin = 0.1, key = ['left'], yr = [0.1, 0.8]
  )
  return 'done'

def gp_ee_hadrons_xsec():
  inDir, outDir = getWorkDirs()
  infile = os.path.join(inDir, 'ee_hadrons_xsec.dat')
  data = np.loadtxt(open(infile, 'rb'))
  data[:,-1] = 0 # ignore systematic uncertainties
  # constants for pQCD curve
  nf, L = 3, 0.217
  beta0 = (11-2./3.*nf)/(4*pi)
  eq = np.array([2./3., -1./3., -1./3.]) # u/d/s
  REW = 3*np.sum(eq*eq)
  eta = np.sum(eq)**2/REW
  cn = [
      1,
      1.9857-0.1152*nf,
      -6.63694-1.20013*nf-0.00518*nf**2-1.240*eta,
      -156.61+18.775*nf-0.7974*nf**2+0.0215*nf**3+(17.828-0.575*nf)*eta
  ]

  # center-of-mass energy Q
  def alpha_s(Q):
      return 1./(beta0*log((Q/L)**2))

  def delta_QCD(Q):
      delta = 0.
      for n, c in enumerate(cn):
          delta += c * (alpha_s(Q)/pi)**(n+1)
          return delta

  pQCD = np.array([
      [roots, REW*(1.+delta_QCD(roots)), 0, 0, 0]
      for roots in np.linspace(1.5, 3)
  ])
  make_plot(
      data = [data, np.array([[1.5, REW, 0, 0, 0],[3, REW, 0, 0, 0]]), pQCD],
      properties = [
          'lc %s lw 2 lt 1 pt 18 ps 0.8' % (default_colors[0]),
          'with lines lt 2 lw 3 lc %s' % (default_colors[1]),
          'with lines lt 2 lw 3 lc %s' % (default_colors[2])
      ],
      titles = ['world data', 'naive quark-parton model', 'perturbative QCD'],
      name = os.path.join(outDir, 'ee_hadrons_xsec'),
      xlabel = '{/Symbol \326}s = M_{ee} (GeV)',
      ylabel = 'R = {/Symbol \163}(e^{+}e^{-}{/Symbol \256}hadrons) / {/Symbol \163}(e^{+}e^{-}{/Symbol \256}{/Symbol \155}^{+}{/Symbol \155}^{-})',
      size = '8in,8in', xr = [0.35, 3], yr = [0.1,60], ylog=True,
      bmargin = 0.13, rmargin = 0.98, tmargin = 0.99,
      gpcalls = ['bars small'], key = ['width -6', 'bottom right']
  )

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  #print gp_rapp()
  gp_ee_hadrons_xsec()
