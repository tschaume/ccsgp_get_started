import logging, argparse, os, sys, re
import numpy as np
from fnmatch import fnmatch
from collections import OrderedDict
from .utils import getWorkDirs, checkSymLink, eRanges, getEnergy4Key
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
  data, cocktail, medium, qgpvac = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
  #scale = { # QM14 (19 GeV skip later, factor here only informational)
  #  '19.6': 1.0340571932983775, '200': 1.0, '39': 0.7776679085382481,
  #  '27': 0.6412140408244136, '62.4': 0.9174700031778402
  #}
  scale = {
    '19.6': 1.0812324298238396, '200': 1.1051002240771077,
    '39': 1.12093451186094, '27': 1.206453891072013,
    '62.4': 1.4992194614005152
  }
  yunit = 1.0e-3 if not diffRel else 1.
  for infile in os.listdir(inDir):
    if infile == "cocktail_contribs": continue
    if fnmatch(infile, '*vacRho*'): continue
    energy = re.compile('\d+').search(infile).group()
    data_type = re.sub('%s\.dat' % energy, '', infile)
    if diffRel and data_type == 'qgp+vac': continue
    if fnmatch(data_type, 'medium*Only'): continue
    energy = getEnergy4Key(energy)
    file_url = os.path.join(inDir, infile)
    data_import = np.loadtxt(open(file_url, 'rb'))
    if (data_type == 'cocktail' or fnmatch(data_type, '*medium*')) \
       and (version == 'QM14' and energy != '19.6'):
       data_import[:,(1,3,4)] /= scale[energy]
    elif data_type == 'data' and version == 'LatestPatrickJieYi':
       data_import[:,(1,3,4)] *= scale[energy]
    if version == 'LatestPatrickJieYi':
       data_import = data_import[(data_import[:,0] > 0.15) & (data_import[:,0] < 1.0)]
    if data_type == 'data': data[energy] = data_import
    elif data_type == 'cocktail': cocktail[energy] = data_import
    elif not diffRel and data_type == 'qgp+vac':
        if noxerr: data_import[:,2] = 0.
        data_import[:,1] /= yunit
        qgpvac[energy] = data_import[
            (data_import[:,0] < 0.7425) | (data_import[:,0] > 0.825)
        ]
    elif not nomed: medium[energy] = data_import
  nSetsData = len(data)

  dataOrdered = OrderedDict()
  for energy in sorted(data, key=float):
    # data & bin edges
    # getUArray propagates stat/syst errors separately internally but
    # errors need to be doubled to retrieve correct errors
    uData = getUArray(data[energy])
    eData = getEdges(data[energy])
    uCocktail = getUArray(cocktail[energy])
    eCocktail = getEdges(cocktail[energy])
    loop = [eData]
    if energy in medium:
      uMedium = getUArray(medium[energy])
      eMedium = getEdges(medium[energy])
      loop.append(eMedium)
    # loop data/medium bins
    for l, eArr in enumerate(loop):
      for i, (e0, e1) in enumzipEdges(eArr):
        logging.debug('%s/%d> %g - %g:' % (energy, l, e0, e1))
        # get cocktail sum in data bin range
        # value+/-0.5*tot.uncert.
        uCocktailSum = getCocktailSum(e0, e1, eCocktail, uCocktail)
        # calc. difference and divide by data binwidth again
        # + set data point
        if not l:
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
        else:
          uDiff = uMedium[i]
          if diffRel:
            uDiff /= uCocktailSum
          else:
            uDiff -= uCocktailSum
            uDiff /= medium[energy][i,2] * 2 * yunit
          dp = [
            medium[energy][i,0], uDiff.nominal_value,
            medium[energy][i,2] if not noxerr else 0.,
            0., 0. # both errors included in data points
          ]
          key = ' '.join([energy, 'GeV (Med.)'])
        # build list of data points
        if dp[0] > 0.7425 and dp[0] < 0.825: continue # mask out omega region
        if dp[0] > 0.97 and dp[0] < 1.0495: continue # mask out phi region
        if key in dataOrdered:
            dataOrdered[key] = np.vstack([dataOrdered[key], dp])
        else:
            dataOrdered[key] = np.array([ dp ])
    if energy not in medium: # TODO add pseudo-data for missing medium (27GeV)
      key = ' '.join([energy, 'GeV (Med.)'])
      dataOrdered[key] = np.array([ [0.1,-1,0,0,0], [1,-1,0,0,0] ])
    if energy in qgpvac:
      dataOrdered[' '.join([energy, 'GeV (QgpVac.)'])] = qgpvac[energy]

  # make plot
  nSets = len(dataOrdered)
  nCats = 2 if diffRel else 3
  nSetsPlot = nSets/nCats if nSets > nSetsData else nSets
  props = [
    'lt 1 lw 6 ps 3 lc %s pt 18' % default_colors[i] for i in xrange(nSetsPlot)
  ]
  titles = dataOrdered.keys()
  if nSets > nSetsData:
    props = zip_flat(props, *[
        [
            'with lines lt %d lw 6 lc %s' % (j+1, default_colors[i])
            for i in xrange(nSetsPlot)
        ]
        for j in xrange(nCats-1)
    ])
    titles = zip_flat(dataOrdered.keys()[::nCats], *[ [''] * nSetsPlot for j in xrange(nCats-1) ])
  global labels
  labels = {
    '{BES: STAR Preliminary}' if version == 'QM12Latest200' or \
      version == 'QM14' or version == 'LatestPatrickJieYi'
    else 'STAR Preliminary': [0.4,0.10,False],
    '{200 GeV: PRL 113 022301' if version == 'QM12Latest200' \
      or version == 'QM14' or version == 'LatestPatrickJieYi'
    else '': [0.4,0.05,False],
  }
  yr = [1.,40] if diffRel else [0.15,1e3]
  if noxerr:
      shift = {
          '19': 1., '27': 3., '39': 10., '62': 30., '200': 100.
      } if not diffRel else {
          '19': 1., '27': 2., '39': 4., '62': 8., '200': 15.
      }
      for k,d in dataOrdered.iteritems():
          energy = re.compile('\d+').search(k).group()
          d[:,(1,3,4)] *= shift[energy]
  gpcalls = [
      'object 1 rectangle back fc rgb "grey" from 0.7425,%f to 0.825,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'object 2 rectangle back fc rgb "grey" from 0.96,%f to 1.0495,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'object 3 rectangle back fc rgb "#C6E2FF" from 0.4,%f to 0.74,%f' % \
      (1.7 if diffRel else 0.5, yr[1]),
      'boxwidth 0.01 absolute',
  ]
  if diffRel:
      gpcalls += [
          'format y "%g"',
          'ytics (1, 2,"" 3,"" 4,5,"" 6,"" 7,"" 8,"" 9, 10, 20, 40)',
      ]
  make_plot(
    data = dataOrdered.values(),
    properties = props, titles = titles,
    name = os.path.join(outDir, 'diff%s%s%s%s' % (
      'Rel' if diffRel else 'Abs', version,
      'NoMed' if nomed else '', 'NoXErr' if noxerr else ''
    )),
    xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
    ylabel = 'Enhancement Ratio' if diffRel else 'Excess Yield ({/Symbol \264} 10^{-3})',
    labels = labels, ylog = True,
    xr = [0.27,0.97], yr = yr,
    key = ['at graph 1.,1.1', 'maxrows 1', 'width -1.5'],
    #lines = { ('x=1' if diffRel else 'x=0'): 'lc 0 lw 4 lt 2' },
    gpcalls = gpcalls,
    lmargin = 0.12, bmargin = 0.12, tmargin = 0.9, rmargin = 0.98,
    size = '12in,9in', arrow_length = 0.4, arrow_offset = 0.9,
  )

  if nomed or noxerr or version == 'QM12': return 'done'

  # integrated enhancement factor
  if diffRel:
    enhance = {}
    data_enhance, medium_enhance = None, None
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
            dpM = [ float(energy), uEnhanceMed.nominal_value, 0, 0, 0 ]
            if medium_enhance is None: medium_enhance = [ dpM ]
            else: medium_enhance.append(dpM)
        else: # uEnhance's are dicts of ufloats
          for k in uEnhanceData:
            uEnhanceData[k] /= uEnhanceCocktail[k]
            dp = [
              float(energy), uEnhanceData[k].nominal_value,
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
    xfacs = os.path.join(outDir, 'xfacs%s.dat' % version)
    if os.path.exists(xfacs): os.remove(xfacs)
    fSystLMR = open(xfacs, 'ab')
    for k in sorted(enhance.keys()):
      np.savetxt(fSystLMR, enhance[k], fmt = '%g', header = k, comments = '\n\n')
    fSystLMR.close()
    yr_upp = 4 if version == 'QM12Latest200' or version == 'QM14' else 7
    if version == 'LatestPatrickJieYi': yr_upp = 5.2
    labels.update({
        '{LMR: %.2f < M_{ee} < %.2f GeV/c^{2}}' % (eRanges[1], eRanges[2]): [0.4,0.15,False]
    })
    make_plot(
      data = [
          np.array([[17.3,2.73,0,0.25,1.47]]),
          np.array([[200,4.7,0,0.4,1.5]]),
          np.array(data_enhance), np.array(medium_enhance)
      ],
      properties = [
          'lt 1 lw 6 ps 3 lc %s pt 18' % default_colors[1],
          'lt 1 lw 6 ps 3 lc %s pt 18' % default_colors[3],
          'lt 1 lw 6 ps 3 lc %s pt 18' % default_colors[0],
          'with lines lt 1 lw 6 lc %s' % default_colors[4],
          ],
      titles = [ 'CERES Pb+Au', 'PHENIX Au+Au', 'STAR Au+Au', 'HMBT + QGP' ],
      name = os.path.join(outDir, 'enhance%s' % version),
      xlabel = '{/Symbol \326}s_{NN} (GeV)',
      ylabel = 'LMR Enhancement Factor',
      xlog = True, key = [ 'at graph 0.7,0.98' ],
      lmargin = 0.15, bmargin = 0.12, tmargin = 0.89, size = '10in,10in',
      yr = [1.,yr_upp], xr = [14,220], gpcalls = [
        'format x "%g"',
        'xtics (10, 20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        'boxwidth 0.03 absolute',
        'label 50 "{/=18 0.2 < M_{ee} < 0.6 GeV/c^{2}}" at 15.5,3 tc %s rotate center' % default_colors[1],
        'label 51 "{/=18 0.15 < M_{ee} < 0.75 GeV/c^{2}}" at 180,4.2 tc %s rotate center' % default_colors[3]
      ], labels = labels
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
    if fnmatch(k, '*QgpVac.*'): suffix = '_QgpVac'
    exc = getMassRangesSums(np.array(v), onlyLMR = True)
    if divdNdy: exc /= dNdyPi0[energy] * 1e-2
    dp = [
        float(energy), exc.nominal_value, 0,
        getErrorComponent(exc, 'stat'), getErrorComponent(exc, 'syst')
    ]
    key = 'LMR' + suffix
    if key not in excess: excess[key] = [ dp ]
    else: excess[key].append(dp)
  logging.debug(excess)
  avdata = np.array(excess['LMR'])
  avg = np.average(avdata[:,1], weights = avdata[:,3])
  graph_data = [
      np.array([
          [ 7.7, avg, 0, 0, avdata[-1][-1]],
          [ 200, avg, 0, 0, avdata[-1][-1]]
      ]),
      np.array([
          [ 7.7, 2*avg, 0, 0, 0], [ 19.6, avg, 0, 0, 0],
      ]),
      np.array(excess['LMR']),
      np.array(excess['LMR_Med']),
      np.array(excess['LMR_QgpVac']),
  ]
  props = [
      'with filledcurves pt 0 lc %s lw 6 lt 2' % default_colors[8],
      'with lines lc %s lw 8 lt 2' % default_colors[1],
      'lt 1 lw 6 ps 3 lc %s pt 18' % default_colors[0],
      'with lines lt 1 lw 6 lc %s' % default_colors[4],
      'with lines lt 1 lw 6 lc %s' % default_colors[6],
  ]
  tits = [
      'BES-I extrapolation for BES-II',
      'model expectation at BES-II',
      'STAR Au+Au',
      'HMBT + QGP',
      'Vacuum + QGP',
  ]
  yr_upp = 4.5 if version == 'QM12Latest200' or version == 'QM14' else 7
  if version == 'LatestPatrickJieYi': yr_upp = 2.5 if divdNdy else 2.
  if divdNdy:
    labels.update(dict((str(v), [float(k)*0.9,yr_upp*1.05,True]) for k,v in dNdyPi0.items()))
    labels.update({ 'dN/dy|_{/Symbol \\160}': [100,yr_upp*1.05,True]})
  labels.update({
      '{LMR: %.2f < M_{ee} < %.2f GeV/c^{2}}' % (eRanges[1], eRanges[2]): [0.4,0.15,False],
  })
  make_plot(
    data = graph_data, properties = props, titles = tits,
    name = os.path.join(outDir, 'excess%s%s' % (version,'DivdNdy' if divdNdy else '')),
    xlabel = '{/Symbol \326}s_{NN} (GeV)',
    ylabel = 'LMR Excess Yield %s({/Symbol \264} 10^{-%d})' % (
        '/ dN/dy|_{/Symbol \\160}  ' if divdNdy else '', 5 if divdNdy else 3
    ),
    xlog = True, xr = [7,220], key = ['at graph 0.9,0.98', 'width -7'],
    lmargin = 0.16, bmargin = 0.12, tmargin = 0.89, size = '10in,10in',
    yr = [0.3,yr_upp], gpcalls = [
      'format x "%g"',
      'xtics (7,10,20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
      'boxwidth 0.03 absolute',
    ], labels = labels,
  )
  return 'done'

def gp_rdiff_merged(version, divdNdy):
  # merge plots for excess yields and enhancement factors if output files exist
  inDir, outDir = getWorkDirs() # inDir not used
  enhance_datdir = os.path.join(outDir, 'enhance%s' % version)
  excess_datdir = os.path.join(outDir, 'excess%s%s' % (version,'DivdNdy' if divdNdy else ''))
  print enhance_datdir, excess_datdir
  weird_key = 'LMR Excess Yield %s({/Symbol \264} 10^{-%d})' % (
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
            'with filledcurves pt 0 lc %s lw 6 lt 2' % default_colors[8]
        ] + [
            'with lines lc %s lw 10 lt 2' % default_colors[3]
        ] + [
          'lt 1 lw 6 ps 1.5 lc %s pt %d' % (default_colors[i], 18+i) for i in xrange(1) #2
        ] + [
          'with lines lt %d lw 6 lc %s' % (i+2, default_colors[i]) for i in xrange(1) #2
        ],
        titles = data.keys(),
        name = os.path.join(outDir, 'enhanceexcess%s' % version),
        xlabel = '{/Symbol \326}s_{NN} (GeV)', ylabel = '',
        lmargin = 0.02, rmargin = 0.99, xlog = True,
        xr = [7,220], key = ['width -10'],#, 'font ",18"', 'spacing 0.9'],
        yr = [0.,5 if version == 'QM12Latest200' or version == 'QM14' else 7],
        labels = labels,
        gpcalls = [
          'format x "%g"',
          'xtics (7,10,20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        ],
      )
  return 'done'

if __name__ == '__main__':
  checkSymLink()
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
