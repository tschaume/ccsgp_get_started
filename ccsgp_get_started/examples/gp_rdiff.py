import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, eRanges, getEnergy4Key
from .utils import getUArray, getEdges, getCocktailSum, enumzipEdges, getMassRangesSums
from .utils import getErrorComponent
from ..ccsgp.ccsgp import make_plot
from ..ccsgp.utils import getOpts, zip_flat
from ..ccsgp.config import default_colors
from uncertainties import ufloat

labels = None

dNdyPi0 = { '19.6': 52.8, '27': 57.6, '39': 60.8, '62.4': 77.2, '200': 105 }

def gp_rdiff(version, nomed, noxerr, diffRel, divdNdy):
  """example for ratio or difference plots using QM12 data (see gp_panel)

  - uses uncertainties package for easier error propagation and rebinning
  - stat. error for medium = 0!
  - stat. error for cocktail ~ 0!
  - statistical error bar on data stays the same for diff
  - TODO: implement ratio!
  - TODO: adjust statistical error on data for ratio!
  - TODO: adjust name and ylabel for ratio

  .. image:: pics/diffAbsQM12.png
     :width: 450 px

  :param version: plot version
  :type version: str
  :param nomed: don't plot medium
  :type nomed: bool
  :param noxerr: don't plot x-errors
  :type noxerr: bool
  """
  inDir, outDir = getWorkDirs()
  inDir = os.path.join(inDir, version)
  data, cocktail, medium, rhofo, vacrho = \
          OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
  #scale = { # QM14 (19 GeV skip later, factor here only informational)
  #  '19.6': 1.0340571932983775, '200': 1.0, '39': 0.7776679085382481,
  #  '27': 0.6412140408244136, '62.4': 0.9174700031778402
  #}
  scale = {
      '19.6': 1.3410566491548412, '200': 1.1051002240771077,
      '39': 1.2719203877292842, '27': 1.350873678084769,
      '62.4': 1.2664666321635087
  }
  yunit = 1.0e-3 if not diffRel else 1.
  for infile in os.listdir(inDir):
    if infile == "cocktail_contribs": continue
    energy = re.compile('\d+').search(infile).group()
    data_type = re.sub('%s\.dat' % energy, '', infile)
    energy = getEnergy4Key(energy)
    file_url = os.path.join(inDir, infile)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if data_type != 'data' and (
        (version == 'QM14' and energy != '19.6') or version == 'LatestPatrickJieYi'
    ):
       data_import[:,(1,3,4)] /= scale[energy]
    if version == 'LatestPatrickJieYi':
        if data_type == 'data':
            data_import = data_import[(data_import[:,0] > 0.14) & (data_import[:,0] < 1.0)]
        else:
            data_import = data_import[data_import[:,0] < 1.0]
    if data_type == 'data': data[energy] = data_import
    elif data_type == 'cocktail': cocktail[energy] = data_import
    elif data_type == 'rho' or data_type == 'vacRho' or data_type == 'medium':
        if noxerr and not diffRel: data_import[:,2:] = 0.
        data_import[:,1] /= yunit
        if data_type == 'rho':
            mask = data_import[:,1] > 0.1
            rhofo[energy] = data_import if diffRel else data_import[mask]
        elif data_type == 'vacRho':
            mask = (data_import[:,0] > 0.35) & (data_import[:,1] > 0.01)
            vacrho[energy] = data_import if diffRel else data_import[mask]
        elif not nomed and data_type == 'medium':
            medium[energy] = data_import
  nSetsData = len(data)

  shift = { '19.6': '1e0', '27': '1e1', '39': '1e2', '62.4': '1e3', '200': '1e4'
  } if not diffRel else {
      '19.6': '1', '27': '8', '39': '50', '62.4': '200', '200': '900'
  }
  dataOrdered = OrderedDict()
  for energy in sorted(data, key=float, reverse=True):
    # data & bin edges
    # getUArray propagates stat/syst errors separately internally but
    # errors need to be doubled to retrieve correct errors
    uData = getUArray(data[energy])
    eData = getEdges(data[energy])
    uCocktail = getUArray(cocktail[energy])
    eCocktail = getEdges(cocktail[energy])
    loop = [eData]
    if energy in medium and diffRel:
      uMedium = getUArray(medium[energy])
      eMedium = getEdges(medium[energy])
      loop.append(eMedium)
    if energy in rhofo and diffRel:
      uRho = getUArray(rhofo[energy])
      eRho = getEdges(rhofo[energy])
      loop.append(eRho)
    if energy in vacrho and diffRel:
      uVacRho = getUArray(vacrho[energy])
      eVacRho = getEdges(vacrho[energy])
      loop.append(eVacRho)
    # loop data/medium bins
    for l, eArr in enumerate(loop):
      for i, (e0, e1) in enumzipEdges(eArr):
        logging.debug('%s/%d> %g - %g:' % (energy, l, e0, e1))
        # get cocktail sum in data bin range
        # value+/-0.5*tot.uncert.
        uCocktailSum = getCocktailSum(e0, e1, eCocktail, uCocktail)
        if uCocktailSum == 0.: continue
        # calc. difference and divide by data binwidth again
        # + set data point
        if l == 0:
          uDiff = uData[i] # value+/-0.5*tot.uncert.
          if diffRel:
            uDiff /= uCocktailSum # value+/-0.5*tot.uncert.
          else:
            uDiff -= uCocktailSum
            uDiff /= data[energy][i,2] * 2 * yunit
          dp = [
            data[energy][i,0], uDiff.nominal_value,
            data[energy][i,2] if not noxerr else 0.,
            getErrorComponent(uDiff, 'stat'),
            getErrorComponent(uDiff, 'syst')
          ]
          key = ' '.join([energy, 'GeV'])
          if noxerr:
              if diffRel:
                  key += ' {/Symbol \264} %s' % shift[energy]
              else:
                  expon = shift[energy].split('e')[1]
                  key += ' {/Symbol \264} 10^{%s}' % expon
        elif l == 1:
          # only done if diffRel
          uDiff = uMedium[i]
          uDiff /= uCocktailSum
          dp = [
            medium[energy][i,0], uDiff.nominal_value+1,
            medium[energy][i,2] if not noxerr else 0.,
            0., 0. # both errors included in data points
          ]
          key = ' '.join([energy, 'GeV (Med.)'])
        elif l == 2:
          # only done if diffRel
          uDiff = uRho[i]
          uDiff /= uCocktailSum
          dp = [
            rhofo[energy][i,0], uDiff.nominal_value+1.,
            rhofo[energy][i,2] if not noxerr else 0.,
            0., 0. # both errors included in data points
          ]
          key = ' '.join([energy, 'GeV (RhoFO.)'])
        elif l == 3:
          # only done if diffRel
          uDiff = uVacRho[i]
          uDiff /= uCocktailSum
          dp = [
            vacrho[energy][i,0], uDiff.nominal_value+1.,
            vacrho[energy][i,2] if not noxerr else 0.,
            0., 0. # both errors included in data points
          ]
          key = ' '.join([energy, 'GeV (VacRho.)'])
        # build list of data points
        if diffRel or l == 0:
          if dp[0] > 0.7425 and dp[0] < 0.825: continue # mask out omega region
          if dp[0] > 0.97 and dp[0] < 1.0495: continue # mask out phi region
        if key in dataOrdered:
            dataOrdered[key] = np.vstack([dataOrdered[key], dp])
        else:
            dataOrdered[key] = np.array([ dp ])
    if not diffRel:
      if energy in medium:
        dataOrdered[' '.join([energy, 'GeV (Med.)'])] = medium[energy]
      if energy in rhofo:
        dataOrdered[' '.join([energy, 'GeV (RhoFO.)'])] = rhofo[energy]
      if energy in vacrho:
        dataOrdered[' '.join([energy, 'GeV (VacRho.)'])] = vacrho[energy]

  # make plot
  nSets = len(dataOrdered)
  nCats = 4
  nSetsPlot = nSets/nCats if nSets > nSetsData else nSets
  props = [
    'lt 1 lw 4 ps 1.5 lc %s pt 18' % default_colors[i]
      for i in reversed(range(nSetsPlot))
  ]
  titles = dataOrdered.keys()
  if nSets > nSetsData:
    props = zip_flat(props, *[
        [
            'with lines lt %d lw 4 lc %s' % (j+1, default_colors[i])
            for i in reversed(range(nSetsPlot))
        ]
        for j in xrange(nCats-1)
    ])
    titles = zip_flat(dataOrdered.keys()[::nCats], *[ [''] * nSetsPlot for j in xrange(nCats-1) ])
  global labels
  labels = {
    '{BES: STAR Preliminary}' if version == 'QM12Latest200' or \
      version == 'QM14' or version == 'LatestPatrickJieYi'
    else 'STAR Preliminary': [
        0.4 if diffRel else 0.2,0.09 if not diffRel and noxerr else 0.75,False
    ],
    '{200 GeV: PRL 113 022301' if version == 'QM12Latest200' \
      or version == 'QM14' or version == 'LatestPatrickJieYi'
    else '': [0.4 if diffRel else 0.2,0.04 if not diffRel and noxerr else 0.7,False],
  }
  yr = [.6,2.5e3] if diffRel else [0.05,1.5e5]
  if noxerr:
      for k,d in dataOrdered.iteritems():
          energy = getEnergy4Key(re.compile('\d+').search(k).group())
          d[:,(1,3,4)] *= float(shift[energy])
  gpcalls = [
      'object 1 rectangle back fc rgb "grey" from 0.75,%f to 0.825,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'object 2 rectangle back fc rgb "grey" from 0.96,%f to 1.0495,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'object 3 rectangle back fc rgb "#C6E2FF" from 0.4,%f to 0.75,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'boxwidth 0.01 absolute',
  ]
  hline = 1. if diffRel else .5
  lines = dict(
      (('x=%g' % (hline*float(shift[energy]))), 'lc rgb "black" lw 4 lt 4')
      for energy in shift
  )
  pseudo_point = np.array([[-1,1,0,0,0]])
  make_plot(
    data = dataOrdered.values() + [
        pseudo_point, pseudo_point, pseudo_point, pseudo_point
    ],
    properties = props + [
        'with lines lt %d lw 4 lc rgb "black"' % (lt+1)
        for lt in xrange(nCats)
    ],
    titles = titles + [
        'HMBT + QGP', 'BW/FO-{/Symbol \162}', '{/Symbol \162}/{/Symbol \167} VacSF+FB+FO',
        'baseline', #'%g%s' % (hline, ' {/Symbol \264} 10^{-3}' if not diffRel else '')
    ],
    name = os.path.join(outDir, 'diff%s%s%s%s' % (
      'Rel' if diffRel else 'Abs', version,
      'NoMed' if nomed else '', 'NoXErr' if noxerr else ''
    )),
    xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
    ylabel = 'Enhancement Ratio' if diffRel else 'Excess Yield / dM_{ee} ({/Symbol \264} 10^{-3} (GeV/c^2)^{-1})',
    #labels = labels,
    xr = [0.18,0.97], yr = yr, ylog = True,
    key = ['at graph 0.96,1.17', 'maxrows 3', 'width -4', 'nobox', 'samplen 0.9'],
    lines = lines if noxerr else {},
    gpcalls = gpcalls,
    lmargin = 0.17, bmargin = 0.1, tmargin = 0.86, rmargin = 0.98,
    size = '9in,11in', arrow_offset = 0.9, #arrow_length = 0.4,
  )

  if nomed or noxerr or version == 'QM12': return 'done'

  # integrated enhancement factor
  if diffRel:
    enhance = {}
    data_enhance, medium_enhance, rhofo_enhance, vacrho_enhance = None, None, None, None
    for energy in sorted(data, key=float):
      for systLMR in [False, True]:
        suffix = str(energy)
        uEnhanceData = getMassRangesSums(
          data[energy], onlyLMR = True,
          systLMR = systLMR, suffix = suffix
        )
        uEnhanceCocktail = getMassRangesSums(
          cocktail[energy], onlyLMR = True,
          systLMR = systLMR, suffix = suffix
        )
        if energy in medium:
          uEnhanceMed = getMassRangesSums(
            medium[energy], onlyLMR = True,
            systLMR = systLMR, suffix = suffix
          )
        if energy in rhofo:
          uEnhanceRhoFO = getMassRangesSums(
            rhofo[energy], onlyLMR = True,
            systLMR = systLMR, suffix = suffix
          )
        if energy in vacrho:
          uEnhanceVacRho = getMassRangesSums(
            vacrho[energy], onlyLMR = True,
            systLMR = systLMR, suffix = suffix
          )
        if not systLMR: # uEnhance's are ufloats
          uEnhanceData /= uEnhanceCocktail
          dp = [
              float(energy), uEnhanceData.nominal_value, 0,
              getErrorComponent(uEnhanceData, 'stat'),
              getErrorComponent(uEnhanceData, 'syst')
          ]
          if data_enhance is None: data_enhance = [ dp ]
          else: data_enhance.append(dp)
          if energy in medium:
            uEnhanceMed /= uEnhanceCocktail
            dpM = [ float(energy), uEnhanceMed.nominal_value+1., 0, 0, 0 ]
            if medium_enhance is None: medium_enhance = [ dpM ]
            else: medium_enhance.append(dpM)
          if energy in rhofo:
            uEnhanceRhoFO /= uEnhanceCocktail
            dpM = [ float(energy), uEnhanceRhoFO.nominal_value+1., 0, 0, 0 ]
            if rhofo_enhance is None: rhofo_enhance = [ dpM ]
            else: rhofo_enhance.append(dpM)
          if energy in vacrho:
            uEnhanceVacRho /= uEnhanceCocktail
            dpM = [ float(energy), uEnhanceVacRho.nominal_value+1., 0, 0, 0 ]
            if vacrho_enhance is None: vacrho_enhance = [ dpM ]
            else: vacrho_enhance.append(dpM)
        else: # uEnhance's are dicts of ufloats
          for k in uEnhanceData:
            uEnhanceData[k] /= uEnhanceCocktail[k]
            dp = [
              float(energy), uEnhanceData[k].nominal_value, 0,
              getErrorComponent(uEnhanceData[k], 'stat'),
              getErrorComponent(uEnhanceData[k], 'syst')
            ]
            rngstr = k.split('_')[-1]
            data_key = 'data_' + rngstr
            if data_key not in enhance: enhance[data_key] = [ dp ]
            else: enhance[data_key].append(dp)
            if k in uEnhanceMed:
              uEnhanceMed[k] /= uEnhanceCocktail[k]
              dpM = [ float(energy), uEnhanceMed[k].nominal_value ]
              med_key = 'model_' + rngstr
              if med_key not in enhance: enhance[med_key] = [ dpM ]
              else: enhance[med_key].append(dpM)
            if k in uEnhanceRhoFO:
              uEnhanceRhoFO[k] /= uEnhanceCocktail[k]
              dpM = [ float(energy), uEnhanceRhoFO[k].nominal_value+1. ]
              rhofo_key = 'rhofo_' + rngstr
              if rhofo_key not in enhance: enhance[rhofo_key] = [ dpM ]
              else: enhance[rhofo_key].append(dpM)
            if k in uEnhanceVacRho:
              uEnhanceVacRho[k] /= uEnhanceCocktail[k]
              dpM = [ float(energy), uEnhanceVacRho[k].nominal_value+1. ]
              vacrho_key = 'vacrho_' + rngstr
              if vacrho_key not in enhance: enhance[vacrho_key] = [ dpM ]
              else: enhance[vacrho_key].append(dpM)
    xfacs = os.path.join(outDir, 'xfacs%s.dat' % version)
    if os.path.exists(xfacs): os.remove(xfacs)
    fSystLMR = open(xfacs, 'ab')
    for k in sorted(enhance.keys()):
      np.savetxt(fSystLMR, enhance[k], fmt = '%g', header = k, comments = '\n\n')
    fSystLMR.close()
    yr_upp = 4 if version == 'QM12Latest200' or version == 'QM14' else 7
    if version == 'LatestPatrickJieYi': yr_upp = 5.5
    #labels.update({
    #    '{LMR: %.2f < M_{ee} < %.2f GeV/c^{2}}' % (eRanges[1], eRanges[2]): [0.4,0.15,False]
    #})
    make_plot(
      data = [
          pseudo_point, pseudo_point, pseudo_point,
          np.array([[17.3,2.73,0,0.25,1.47]]),
          np.array([[200,4.7,0,0.4,1.5]]),
          np.array(enhance['data_0.15-0.75']),
          np.array(enhance['data_0.4-0.75']),
          np.array(medium_enhance),
          np.array(rhofo_enhance), np.array(vacrho_enhance)
      ],
      properties = [
          'lt 1 lw 4 ps 2 lc rgb "white" pt 19',
          'lt 1 lw 4 ps 2 lc rgb "white" pt 20',
          'lt 1 lw 4 ps 2 lc rgb "white" pt 18',
          'lt 1 lw 4 ps 2 lc %s pt 19' % default_colors[1],
          'lt 1 lw 4 ps 2 lc %s pt 20' % default_colors[3],
          'lt 1 lw 4 ps 2 lc %s pt 18' % default_colors[4],
          'lt 1 lw 4 ps 2 lc %s pt 18' % default_colors[0],
          'with lines lt 2 lw 4 lc %s' % default_colors[-1],
          'with lines lt 3 lw 4 lc %s' % default_colors[-1],
          'with lines lt 4 lw 4 lc %s' % default_colors[-1],
      ],
      titles = [
          'CERES Pb+Au', 'PHENIX Au+Au', 'STAR Au+Au',
          '', '', '', '',
          'HMBT + QGP', 'BW/FO-{/Symbol \162}', '{/Symbol \162}/{/Symbol \167} VacSF+FB',
      ],
      name = os.path.join(outDir, 'enhance%s' % version),
      xlabel = '{/Symbol \326}s_{NN} (GeV)',
      ylabel = 'LMR Enhancement Factor',
      xlog = True, key = [ 'at graph 0.9,0.98', 'nobox', 'maxrows 4' ],
      size = '10in,8in', bmargin = 0.13, tmargin = 0.92, rmargin = 0.99,
      yr = [1.,yr_upp], xr = [14,220], gpcalls = [
        'format x "%g"',
        'xtics (10, 20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        'boxwidth 0.025 absolute',
        'label 50 "{/=18 0.2 < M_{ee} < 0.6 GeV/c^{2}}" at 15.5,3 tc %s rotate center' % default_colors[1],
        'label 51 "{/=18 0.15 < M_{ee} < 0.75 GeV/c^{2}}" at 180,4.2 tc %s rotate center' % default_colors[3],
        'label 52 "{/=18 0.4 < M_{ee} < 0.75 GeV/c^{2}}" at 75,3.1 tc %s rotate by -20' % default_colors[0],
        'label 53 "{/=18 0.15 < M_{ee} < 0.75 GeV/c^{2}}" at 50,1.2 tc %s' % default_colors[4]
      ], #labels = labels
    )
    return 'done'

  # integrated excess yield in mass ranges
  excess = {}
  for k, v in dataOrdered.iteritems():
    suffix = ''
    energy = getEnergy4Key(re.compile('\d+').search(k).group())
    if fnmatch(k, '*Med.*'):
      suffix = '_Med'
      if version != 'LatestPatrickJieYi' and energy == '27': continue # TODO
    if fnmatch(k, '*RhoFO.*'): suffix = '_RhoFO'
    if fnmatch(k, '*VacRho.*'): suffix = '_VacRho'
    exc = getMassRangesSums(np.array(v), onlyLMR = True)
    if divdNdy: exc /= dNdyPi0[energy] * 1e-2
    dp = [
        float(energy), exc.nominal_value, 0,
        getErrorComponent(exc, 'stat'), getErrorComponent(exc, 'syst')
    ]
    if suffix == '_Med' and not diffRel and not divdNdy:
        print dp
    key = 'LMR' + suffix
    if key not in excess: excess[key] = [ dp ]
    else: excess[key].append(dp)
  logging.debug(excess)
  avdata = np.array(excess['LMR'])
  avg = np.average(avdata[:,1], weights = avdata[:,3])
  graph_data = [
      np.array([
          [ 7.7, avg, 0, 0, avdata[-1][-1]],
          [ 19.6, avg, 0, 0, avdata[-1][-1]]
      ]),
      np.array([
          [ 19.6, avg, 0, 0, 0], [ 200., avg, 0, 0, 0]
      ]),
      np.array([
          [ 7.7, 2*avg, 0, 0, 0], [ 19.6, avg, 0, 0, 0],
      ]),
      np.array(excess['LMR']),
      np.array(excess['LMR_Med']),
      np.array(excess['LMR_VacRho']),
      np.array(excess['LMR_RhoFO']),
  ]
  props = [
      'with filledcurves pt 0 lc %s lw 4 lt 2' % default_colors[8],
      'with lines lc %s lw 4 lt 2' % default_colors[8],
      'with lines lc %s lw 8 lt 2' % default_colors[1],
      'lt 1 lw 4 ps 2 lc %s pt 18' % default_colors[0],
      'with lines lt 2 lw 4 lc %s' % default_colors[-1],
      'with lines lt 3 lw 4 lc %s' % default_colors[-1],
      'with lines lt 4 lw 4 lc %s' % default_colors[-1],
  ]
  tits = [
      'BES-I extrapolation', '', 'model expectation', 'STAR Au+Au',
      'HMBT + QGP', '{/Symbol \162}/{/Symbol \167} VacSF+FB', 'BW/FO-{/Symbol \162}',
  ]
  yr_upp = 4.5 if version == 'QM12Latest200' or version == 'QM14' else 7
  if version == 'LatestPatrickJieYi': yr_upp = 2 if divdNdy else 2.
  labels = {}
  if divdNdy:
    labels.update(dict((str(v), [float(k)*0.9,yr_upp*1.05,True]) for k,v in dNdyPi0.items()))
    labels.update({ 'dN/dy|_{/Symbol \\160}': [100,yr_upp*1.05,True]})
  #labels.update({
  #    '{LMR: %.2f < M_{ee} < %.2f GeV/c^{2}}' % (eRanges[1], eRanges[2]): [0.4,0.15,False],
  #})
  make_plot(
    data = graph_data, properties = props, titles = tits,
    name = os.path.join(outDir, 'excess%s%s' % (version,'DivdNdy' if divdNdy else '')),
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    ylabel = 'LMR Excess Yield %s({/Symbol \264} 10^{-%d})' % (
        '/ dN/dy|_{/Symbol \\160}  ' if divdNdy else '', 5 if divdNdy else 3
    ),
    xlog = True, xr = [7,220], size = '10in,8in',
    key = ['at graph 1.05,0.98', 'width -3', 'nobox', 'maxrows 3'],
    bmargin = 0.13, tmargin = 0.92, rmargin = 0.99,
    yr = [0,yr_upp], gpcalls = [
      'format x "%g"',
      'xtics (7,10,20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
      'boxwidth 0.025 absolute',
      'label 52 "{/=18 0.4 < M_{ee} < 0.75 GeV/c^{2}}" at 60,0.4 tc %s' % default_colors[0],
    ], labels = labels,
  )
  return 'done'

def gp_rdiff_merged(version, divdNdy):
  # merge plots for excess yields and enhancement factors if output files exist
  inDir, outDir = getWorkDirs() # inDir not used
  enhance_datdir = os.path.join(outDir, 'enhance%s' % version)
  excess_datdir = os.path.join(outDir, 'excess%s%s' % (version,'DivdNdy' if divdNdy else ''))
  print enhance_datdir, excess_datdir
  weird_key = 'LMR Excess Yield %s({/Symbol \264} 10^{-%d} (GeV/c^2)^{-1}))' % (
      '/ dN/dy|_{/Symbol \\160}  ' if divdNdy else '', 5 if divdNdy else 3
  )
  if os.path.exists(enhance_datdir) and os.path.exists(excess_datdir):
      excess_data = np.loadtxt(
          open(os.path.join(
              excess_datdir,
              'LMR_Excess_Yield_%s_Symbol_10_%d_.dat' % (
                  'dN_dy___Symbol_160' if divdNdy else '', 5 if divdNdy else 3
              )
          ), 'rb')
      )
      avdata = np.array(excess_data)
      avg = np.average(avdata[:,1], weights = avdata[:,4])
      data = OrderedDict()
      data['BES-II Extrapolation'] = np.array([
          [ 7.7, avg, 0, 0, excess_data[-1][-1]],
          [ 9.1, avg, 0, 0, excess_data[-1][-1]],
          [ 11.5, avg, 0, 0, excess_data[-1][-1]],
          [ 14.6, avg, 0, 0, excess_data[-1][-1]],
          [ 19.6, avg, 0, 0, excess_data[-1][-1]]
      ])
      data['Model Expectation at BES-II'] = np.array([
          [ 7.7, 2*avg, 0, 0, 0], [ 19.6, avg, 0, 0, 0],
      ])
      data[weird_key] = excess_data
      #data['LMR Enhancement Factor'] = np.loadtxt(
      #    open(os.path.join(enhance_datdir, 'LMR_Enhancement_Factor.dat'), 'rb')
      #)
      data['Model for Excess'] = np.loadtxt(
          open(os.path.join(excess_datdir, 'Model.dat'), 'rb')
      )
      #data['Model for Enhancement'] = np.loadtxt(
      #    open(os.path.join(enhance_datdir, 'Model.dat'), 'rb')
      #)
      #xshift = 1.05
      #data['LMR Enhancement Factor'][:,0] *= xshift
      #data['Model for Enhancement'][:,0] *= xshift
      if divdNdy:
          labels.update(dict((str(v), [float(k)*0.9,5.2,True]) for k,v in dNdyPi0.items()))
          labels.update({ 'dN/dy|_{/Symbol \\160}': [0.73,1.043,False]})
      make_plot(
        data = data.values(),
        properties = [
            'with filledcurves pt 0 lc %s lw 4 lt 2' % default_colors[8]
        ] + [
            'with lines lc %s lw 10 lt 2' % default_colors[3]
        ] + [
          'lt 1 lw 4 ps 1.5 lc %s pt %d' % (default_colors[i], 18+i) for i in xrange(1) #2
        ] + [
          'with lines lt %d lw 4 lc %s' % (i+2, default_colors[i]) for i in xrange(1) #2
        ],
        titles = data.keys(),
        name = os.path.join(outDir, 'enhanceexcess%s' % version),
        xlabel = '{/Symbol \326}s_{NN} (GeV)', ylabel = '',
        lmargin = 0.02, rmargin = 0.99, xlog = True,
        xr = [7,220], key = ['width -10'],#, 'font ",18"', 'spacing 0.9'],
        yr = [0.,5 if version == 'QM12Latest200' or version == 'QM14' else 7],
        #labels = labels,
        gpcalls = [
          'format x "%g"',
          'xtics (7,10,20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        ],
      )
  return 'done'

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("version", help="version = subdir name of input dir")
  parser.add_argument("--nomed", help="don't plot medium", action="store_true")
  parser.add_argument("--noxerr", help="no dx errors", action="store_true")
  parser.add_argument("--diffRel", help="plot relative difference (ratio)", action="store_true")
  parser.add_argument("--divdNdy", help="divide excess plot by dNdy_pi0", action="store_true")
  parser.add_argument("--log", help="show log output", action="store_true")
  args = parser.parse_args()
  loglevel = 'DEBUG' if args.log else 'WARNING'
  logging.basicConfig(
    format='%(message)s', level=getattr(logging, loglevel)
  )
  print gp_rdiff(args.version, args.nomed, args.noxerr, args.diffRel, args.divdNdy)
  #print gp_rdiff_merged(args.version,args.divdNdy)
