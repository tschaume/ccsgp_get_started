import logging, argparse, os, sys, re, math, random
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, getMassRangesSums
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts
from ..ccsgp.config import default_colors
from pymodelfit import LinearModel
from uncertainties import ufloat
from decimal import Decimal

dataIMRfit_style = 'with lines lc %s lw 4 lt 1' % default_colors[-2]
cocktailIMRfit_style = 'with lines lc %s lw 4 lt 2' % default_colors[-2]
pseudo_point = np.array([ [-1,1e-7,0,0,1] ])

def truncated_gaus(r, mu, sig):
    while 1:
        x = r.gauss(mu, sig)
        if x > 0: return x

def particleLabel4Key(k):
    if k == 'pion': return '{/Symbol \160}^0 {/Symbol \256} e^{+}e^{-}{/Symbol \147}'
    if k == 'eta': return '{/Symbol \150} {/Symbol \256} e^{+}e^{-}{/Symbol \147}'
    if k == 'etap': return '{/Symbol \150}\' {/Symbol \256} e^{+}e^{-}{/Symbol \147}'
    if k == 'omega': return '{/Symbol \167} {/Symbol \256} e^{+}e^{-}({/Symbol \160})'
    if k == 'phi': return '{/Symbol \146} {/Symbol \256} e^{+}e^{-}({/Symbol \150})'
    if k == 'jpsi': return 'J/{/Symbol \171} {/Symbol \256} e^{+}e^{-}'
    if k == 'ccbar':
        return 'c@^{/=18-}c {/Symbol \256} D/{/Symbol \514} {/Symbol \256} e^{+}e^{-}'
    return k

def gp_stack(version, energies, inclMed, inclFits):
  """example for a plot w/ stacked graphs using QM12 data (see gp_panel)

  * how to omit keys from the legend
  * manually add legend entries
  * automatically plot arrows for error bars larger than data point value

  .. image:: pics/stackQM12.png
     :width: 550px

  :param version: plot version / input subdir name
  :type version: str
  """
  inclMed = (inclMed and version != 'QM12')
  inclFits = (inclFits and version == 'LatestPatrickJieYi')
  cocktail_style = 'with filledcurves pt 0 lc %s lw 4 lt 1' % default_colors[8]
  medium_style = 'with lines lc %s lw 4 lt 2' % default_colors[4]
  if inclMed:
    medium_style = 'with filledcurves pt 0 lc %s lw 4 lt 2' % default_colors[4]
  shift = {
    '200': 200., '62': 15., '39': 0.5, '27': 0.01, '19': 2e-4
  } if (
    version != 'QM12' and version != 'Latest19200_PatrickQM12' and version != 'QM12Latest200'
  ) else {
    '200': 200., '62': 20., '39': 1., '19': 0.05
  }
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data, cocktail, medium = OrderedDict(), OrderedDict(), OrderedDict()
  dataIMRfit, cocktailIMRfit, dataTvsS = OrderedDict(), OrderedDict(), OrderedDict()
  cocktailContribs = OrderedDict()
  linmod = LinearModel()
  rangeIMR = [1.15, 2.5]
  nPtsMC = 1000 # number of MC points per data point
  cRanges = map(Decimal, ['0.', '0.1'])
  pi0yld = {}
  for filename in os.listdir(inDir):
    file_url = os.path.join(inDir, filename)
    # take care of cocktail contributions first
    if os.path.isdir(file_url):
        for fn in os.listdir(file_url):
            energy = re.compile('\d+').search(fn).group()
            particle = re.sub('%s\.dat' % energy, '', fn)
            cocktailContribs[particle] = np.loadtxt(open(
                os.path.join(file_url, fn), 'rb'
            ))
            if particle == 'omega' or particle == 'phi' or particle == 'ccbar':
                thr = 0.95 if particle == 'omega' else 1.4
                if particle == 'ccbar': thr = 2.6
                mask = cocktailContribs[particle][:,0] < thr
                cocktailContribs[particle] = cocktailContribs[particle][mask]
            cocktailContribs[particle][:,(1,3,4)] *= shift[energy]
            cocktailContribs[particle][:,2:] = 0
        continue
    # normal input files
    energy = re.compile('\d+').search(filename).group()
    data_type = re.sub('%s\.dat' % energy, '', filename)
    data_import = np.array([[-1, 1, 0, 0, 0]]) if (
      energies is not None and energy not in energies
    ) else np.loadtxt(open(file_url, 'rb'))
    # fit IMR region with exp(-M/kT+C)
    if (
      inclFits and energies is None and
      (data_type == 'data' or data_type == 'cocktail')
    ):
      # data in IMR
      mask = (data_import[:,0] > rangeIMR[0]) & (data_import[:,0] < rangeIMR[1])
      dataIMR = data_import[mask]
      # exp fit in IMR region -> slope parameter
      mIMR, bIMR = linmod.fitErrxy(
          dataIMR[:,0], np.log10(dataIMR[:,1]), dataIMR[:,2],
          np.log10(dataIMR[:,3]) # TODO: include syst. uncertainties
      )
      slope_par = -1./mIMR
      logging.info('%s: m = %g , b = %g => T = %g' % (filename, mIMR, bIMR, slope_par))
      # Monte-Carlo the datapoints within dx/dy -> parameter mean & std.dev.
      slope_par_err = 0.
      if data_type == 'data':
        # one random generator per x,y for each data point in IMR
        rndm = OrderedDict((ax,[]) for ax in ['x','y'])
        rndm_jump = 0
        dataMC = OrderedDict((n, []) for n in xrange(nPtsMC))
        for i,dp in enumerate(dataIMR): # for each datapoint
          logging.info(('MC %d: x = {}, y = {}' % i).format(
            ufloat(dp[0], dp[2]), ufloat(dp[1], dp[3]) # TODO: syst. uncertainties
          ))
          for ax in rndm: # for each axis
            rndm[ax].append(random.Random())
            if ax == 'y' and i == 0: # jumpahead y-axis
              rndm[ax][i].setstate(rndm['x'][i].getstate())
              rndm[ax][i].jumpahead(nPtsMC)
            if i > 0: # jumpahead within axes
              rndm[ax][i].setstate(rndm[ax][i-1].getstate())
              rndm[ax][i].jumpahead(nPtsMC)
          for n in xrange(nPtsMC): # generate nPtsMC new points for current datapoint
            dataMC[n].append([
              rndm['x'][i].uniform(dp[0] - dp[2], dp[0] + dp[2]),
              truncated_gaus(rndm['y'][i], dp[1], dp[3])
            ])
            #logging.info(' %d: x = %g, y = %g' % (n, dataMC[n][i][0], dataMC[n][i][1]))
        mIMRMC = []
        for dp in dataMC.itervalues():
          dp = np.array(dp)
          mIMRMC.append([
            linmod.fitBasic(dp[:,0], np.log10(dp[:,1]))[0][0],
            linmod.stdData(x=dp[:,0], y=np.log10(dp[:,1]))
          ])
        mIMRMC = np.array(mIMRMC)
        weights = 1./mIMRMC[:,1]
        mIMRMC_avg = np.average(mIMRMC[:,0], weights=weights)
        mIMRMC_var = np.average((mIMRMC[:,0]-mIMRMC_avg)**2, weights=weights)
        umIMR = ufloat(mIMRMC_avg, math.sqrt(mIMRMC_var))
        slope_par_err = abs(math.sqrt(mIMRMC_var)/mIMRMC_avg * slope_par)
        logging.info(('=> %g == {} => err = %g' % (mIMR, slope_par_err)).format(umIMR))
      # set IMR slope datapoint
      IMRfit = np.array([ [x, math.pow(10.,mIMR*x+bIMR), 0., 0., 0.] for x in rangeIMR ])
      IMRfit[:,(1,3,4)] *= shift[energy]
      if data_type == 'data': dataIMRfit[energy] = IMRfit
      else: cocktailIMRfit[energy] = IMRfit
      # fill array for T vs sqrt(s) plot
      dp = [ float(getEnergy4Key(energy)), slope_par, 0., slope_par_err, 0. ]
      if data_type in dataTvsS: dataTvsS[data_type].append(dp)
      else: dataTvsS[data_type] = [ dp ]
    if data_type != '+medium':
      pi0yld['_'.join([energy,data_type])] = getMassRangesSums(
        np.copy(data_import), customRanges = cRanges , singleRange = True
      )
    # function changes syst. uncertainties of input numpy array
    # following scaling is wrong for y < 0 && shift != 1
    data_import[:,(1,3,4)] *= shift[energy]
    if fnmatch(filename, 'data*'):
      data[energy] = data_import
    elif fnmatch(filename, 'cocktail*'):
      data_import[:,(2,3)] = 0 # don't plot dx,dy for cocktail
      if inclMed:
          for di in data_import:
              if di[0] < 1.1: di[4] = 0 # don't plot dy2 for cocktail
      if energy == '19' and version == 'QM12':
        # cut off cocktail above 1.1 GeV/c^2
        cocktail[energy] = data_import[data_import[:,0] < 1.3]
      else:
        cocktail[energy] = data_import
    elif inclMed and fnmatch(filename, '+medium*'):
      data_import[:,(2,3)] = 0 # don't plot dx, dy1 for medium
      medium[energy] = data_import
  # calculate data-to-cocktail scaling factors in pi0 region < 0.1 GeV/c2
  # cocktail/data
  scale = {}
  for e in ['19', '27', '39', '62' ]:
    scale[e] = (pi0yld[e+'_cocktail'] / pi0yld[e+'_data']).nominal_value
  scale['200'] = 1.
  if version == 'LatestPatrickJieYi' or version == 'QM14':
    for k in cocktail:
        if k != '19': cocktail[k][:,(1,3,4)] /= scale[k]
    for k in medium:
        if k != '19': medium[k][:,(1,3,4)] /= scale[k]
  print scale
  # ordered
  dataOrdered = OrderedDict(
    (' '.join([
      getEnergy4Key(k), 'GeV', '{/Symbol \264} %g' % shift[k],
      '            STAR Preliminary' if version == 'QM12Latest200' and k == '39' else '',
      '    [arXiv:1312.7397]' if version == 'QM12Latest200' and k == '200' else ''
    ]), data[k]) for k in sorted(data, key=int)
  )
  cocktailOrdered = OrderedDict((k, cocktail[k]) for k in sorted(cocktail, key=int))
  mediumOrdered = OrderedDict((k, medium[k]) for k in sorted(medium, key=int))
  dataIMRfitOrdered = OrderedDict((k, dataIMRfit[k]) for k in sorted(dataIMRfit, key=int))
  cocktailIMRfitOrdered = OrderedDict((k, cocktailIMRfit[k]) for k in sorted(cocktailIMRfit, key=int))
  nSetsData, nSetsCocktail, nSetsMedium = len(dataOrdered), len(cocktail), len(medium)
  nSetsDataIMRfit, nSetsCocktailIMRfit = len(dataIMRfitOrdered), len(cocktailIMRfitOrdered)
  nSetsCocktailContribs = len(cocktailContribs)
  yr_low = 3e-7 if version == 'QM12' else 1e-10
  if version == 'Latest19200_PatrickQM12': yr_low = 1e-7
  if version == 'QM12Latest200': yr_low = 2e-6
  make_plot(
    data = cocktailContribs.values()
    + cocktailOrdered.values() + ([ pseudo_point ] if inclMed else [])
    + mediumOrdered.values() + [ pseudo_point ] + dataOrdered.values()
    + dataIMRfitOrdered.values() + ([ pseudo_point ] if inclFits else [])
    + cocktailIMRfitOrdered.values() + ([ pseudo_point ] if inclFits else []),
    properties = [
      'with lines lc %s lw 4 lt 1' % default_colors[i if i < 5 else i+1]
      for i in xrange(nSetsCocktailContribs)
    ] + [ cocktail_style ] * (nSetsCocktail+1)
    + [ medium_style ] * (nSetsMedium+bool(nSetsMedium)) + [
      'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
      for i in xrange(nSetsData)
    ] + [ dataIMRfit_style ] * (nSetsDataIMRfit+1)
    + [ cocktailIMRfit_style ] * (nSetsCocktailIMRfit+1),
    titles = [ particleLabel4Key(k) for k in cocktailContribs.keys() ]
    + [''] * nSetsCocktail + ['Cocktail w/o {/Symbol \162}']
    + [''] * nSetsMedium + ['+ Medium'] * bool(nSetsMedium) + dataOrdered.keys()
    + [''] * nSetsDataIMRfit + [''] * inclFits
    + [''] * nSetsCocktailIMRfit + [''] * inclFits,
    name = os.path.join(outDir, 'stack%s%s%s%s' % (
      version, 'InclMed' if inclMed else '', 'InclFits' if inclFits else '',
      '_' + '-'.join(energies) if energies is not None else ''
    )),
    ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]',
    xlabel = 'invariant dielectron mass, M_{ee} (GeV/c^{2})',
    ylog = True, xr = [0, 3.5], yr = [yr_low, 3e3],
    lmargin = 0.17, rmargin = 0.97, bmargin = 0.08, arrow_offset = 0.8,
    #tmargin = 0.9 if version != 'QM12Latest200' else 0.99,
    key = [
      'width -7.5', 'maxrows 7', 'font ",19"', 'samplen 0.5', 'spacing 0.9'
    ] if version != 'QM12Latest200' else [
      'width -14', 'maxcols 1'
    ],
    labels = {
      'BES: STAR Preliminary': [0.38,0.9,False],
      '200 GeV: [arXiv:1312.7397]': [0.38,0.85,False],
      '{/Symbol=50 \775}': [0.64,0.81 if not inclMed else 0.75,False]
    } if version == 'QM12Latest200' else {}, size = '10in,13in',
    #arrows = [ # example arrow
    #  [ [2.4, 5e-5], [2.3, 1e-5], 'head filled lc 1 lw 4 lt 1 front' ],
    #],
  )
  if inclFits:
    for t in dataTvsS: dataTvsS[t].sort(key=lambda x: x[0])
    make_plot(
      data = [ np.array(dataTvsS['cocktail']), np.array(dataTvsS['data']) ],
      properties = [
       'with lines lt 2 lw 4 lc %s' % default_colors[0],
       'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[0]
      ],
      titles = [ 'cocktail', 'data' ],
      name = os.path.join(outDir, 'IMRslope%s' % version),
      xlabel = '{/Symbol \326}s_{NN} (GeV)',
      ylabel = 'Slope Parameter T [log(dN/dM) = -M/T+C] (GeV/c^{2})',
      lmargin = 0.1, xlog = True, yr = [0,3.5], gpcalls = [
        'format x "%g"',
        'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
      ]
    )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--energies", nargs='*', help="list of energies to plot (for animation)")
  parser.add_argument("--med", help="include medium calculations", action="store_true")
  parser.add_argument("--fit", help="include IMR fits", action="store_true")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_stack(args.version, args.energies, args.med, args.fit)
