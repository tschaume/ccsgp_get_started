import logging, argparse, os, re
import numpy as np
from ..ccsgp.config import default_colors
from ..ccsgp.ccsgp import make_plot, make_panel
from .utils import getWorkDirs, checkSymLink, getMassRangesSums, getEnergy4Key
from math import pi, log
from collections import OrderedDict

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

def calc_REW_eta(nf):
  charges = [2./3., -1./3., -1./3.] # u/d/s
  eq = np.array(charges[:nf])
  REW = 3*np.sum(eq*eq)
  eta = np.sum(eq)**2/REW
  return REW, eta

def calc_pQCD(nf): # nf = 2,3
  REW, eta = calc_REW_eta(nf)
  limits = [0.5, 1.02, 3]
  L = 0.217
  beta0 = (11-2./3.*nf)/(4*pi)
  cn = [
      1, 1.9857-0.1152*nf,
      -6.63694-1.20013*nf-0.00518*nf**2-1.240*eta,
      -156.61+18.775*nf-0.7974*nf**2+0.0215*nf**3+(17.828-0.575*nf)*eta
  ]
  print 'beta0 = ', beta0, ', eta = ', eta
  print 'cn = ', cn
  alpha_s = lambda Q: 1./(beta0*log((Q/L)**2)) # center-of-mass energy Q

  def delta_QCD(Q):
      delta = 0.
      for n, c in enumerate(cn):
          delta += c * (alpha_s(Q)/pi)**(n+1)
          return delta

  return np.array([
      [roots, REW*(1.+delta_QCD(roots)), 0, 0, 0]
      for roots in np.linspace(limits[nf-2], limits[nf-1])
  ])

def gp_ee_hadrons_xsec():
  inDir, outDir = getWorkDirs()
  infile = os.path.join(inDir, 'ee_hadrons_xsec.dat')
  data = np.loadtxt(open(infile, 'rb'))
  data[:,-1] = 0 # ignore systematic uncertainties
  pQCD = calc_pQCD(2)
  pQCD = np.vstack((pQCD, calc_pQCD(3)))
  REW = [ calc_REW_eta(nf)[0] for nf in [2,3] ]
  print 'REW = ', REW
  rew = np.array([
      [0.5, REW[0], 0, 0, 0], [1.02, REW[0], 0, 0, 0],
      [1.02, REW[1], 0, 0, 0], [3, REW[1], 0, 0, 0],
  ])
  make_plot(
      data = [data, rew, pQCD],
      properties = [
          'lc %s lw 2 lt 1 pt 18 ps 0.8' % (default_colors[0]),
          'with lines lt 2 lw 3 lc %s' % (default_colors[1]),
          'with lines lt 2 lw 3 lc %s' % (default_colors[2])
      ],
      titles = ['world data', 'naive quark-parton model', 'perturbative QCD'],
      name = os.path.join(outDir, 'ee_hadrons_xsec'),
      xlabel = '{/Symbol \326}s = M_{ee} (GeV)',
      ylabel = 'R = {/Symbol \163}(e^{+}e^{-}{/Symbol \256}hadrons) / {/Symbol \163}(e^{+}e^{-}{/Symbol \256}{/Symbol \155}^{+}{/Symbol \155}^{-})',
      size = '10in,7.5in', xr = [0.5, 3], yr = [0.7,60], ylog=True,
      bmargin = 0.14, rmargin = 0.99, tmargin = 0.99,
      gpcalls = ['bars small', 'format y "%g"'], key = ['width -6', 'nobox'],
      labels = {
          '{/Symbol \162}': [0.7,1.1,True], '{/Symbol \167}': [0.8,25,True],
          '{/Symbol \146}': [1.05,40,True], "{/Symbol \162}'": [1.7,4,True],
      }
  )

def gp_rapp_overview_panel():
  inDir, outDir = getWorkDirs()
  energies = ['19', '27', '39', '62']
  subkeys = ['Energy Dependence', '27 GeV Medium Effects']
  dpt_dict = OrderedDict((subkey, [[], [], []]) for subkey in subkeys)
  pseudo_data = np.array([[0,1,0,0,0]])
  for i,title in enumerate([
      'HMBT+QGP', 'HMBT', 'QGP', 'VacSF+FB+FO', 'Cocktail (w/o {/Symbol \162})'
  ]):
      dpt_dict[subkeys[0]][0].append(pseudo_data)
      dpt_dict[subkeys[0]][1].append(
          'with lines lt %d lc rgb "black" lw 5' % (i+1)
      )
      dpt_dict[subkeys[0]][2].append(title)
  for i,modeltype in enumerate(['MedOnly', 'QgpOnly']):
      infile = os.path.join(inDir, 'medium'+modeltype+'19.dat')
      data = np.loadtxt(open(infile, 'rb'))
      data[:,2:] = 0
      dpt_dict[subkeys[0]][0].append(data)
      dpt_dict[subkeys[0]][1].append(
          'with lines lt %d lc %s lw 5' % (i+2, default_colors[0])
      )
      dpt_dict[subkeys[0]][2].append('')
  for i,energy in enumerate(energies):
      infile = os.path.join(inDir, 'medium'+energy+'.dat')
      data = np.loadtxt(open(infile, 'rb'))
      data[:,2:] = 0
      dpt_dict[subkeys[0]][0].append(data)
      dpt_dict[subkeys[0]][1].append('with lines lt 1 lc %s lw 5' % default_colors[i])
      dpt_dict[subkeys[0]][2].append(' '.join([getEnergy4Key(energy), 'GeV']))
  linetypes = [5,4,1]
  for i,infile in enumerate([
      '../../gp_panel/input/LatestPatrickJieYi/cocktail27.dat',
      'vacRho27.dat', 'medium27.dat'
  ]):
      data = np.loadtxt(open(os.path.join(inDir, infile), 'rb'))
      data[:,(2,3)] = 0
      if i != 2: data[:,4] = 0
      dpt_dict[subkeys[1]][0].append(data)
      dpt_dict[subkeys[1]][1].append('with %s lt %d lc %s lw 5' % (
          'filledcurves pt 0' if i == 2 else 'lines', linetypes[i], default_colors[1]
      ))
      dpt_dict[subkeys[1]][2].append('')
  yr = [2e-5, 0.07]
  make_panel(
    dpt_dict = dpt_dict,
    name = os.path.join(outDir, 'rapp_overview_panel'),
    ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]',
    xlabel = 'invariant dielectron mass, M_{ee} (GeV/c^{2})',
    ylog = True, xr = [0.3, 1.45], yr = yr,
    layout = '2x1', size = '4in,9in', key = ['opaque', 'width -5'],
    gpcalls = [
        'object 51 rectangle back fc rgb "grey" from 0.75,%f to 0.825,%f' % (yr[0]*2, yr[1]/4),
        'object 52 rectangle back fc rgb "grey" from 0.95,%f to 1.08,%f' % (yr[0]*2, yr[1]/4),
        'object 53 rectangle back fc rgb "#C6E2FF" from 0.4,%f to 0.75,%f' % (yr[0]*2, yr[1]/4),
    ]
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
  #gp_ee_hadrons_xsec()
  gp_rapp_overview_panel()
