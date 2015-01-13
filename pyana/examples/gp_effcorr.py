import os, argparse, logging, math
from .utils import getWorkDirs, checkSymLink, getEnergy4Key, particleLabel4Key
from collections import OrderedDict
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
import numpy as np
from itertools import groupby
from operator import itemgetter

def gp_syserr():
    energies = ['19.6', '27', '39', '62.4']
    inDir, outDir = getWorkDirs()
    inDir = os.path.join(inDir, 'syserr')
    data = OrderedDict()
    for charge in ['minus', 'plus']:
        subkey = 'positive particles' if charge == 'plus' else 'negative particles'
        data[subkey] = [[], [], []]
        shift, factor = 2., -0.5
        for p,particle in enumerate(['e', 'pi']):
            key = '{/Symbol \160}' if particle == 'pi' else 'e'
            key += '^{+}' if charge == 'plus' else '^{-}'
            shiftPt, factorPt = 0.3, -3
            for ipt,pt in enumerate([0.35+0.3*i for i in range(6)]):
                ptstr = 'Pt%.2f' % pt
                fname = particle + charge + ptstr + '.dat'
                data[subkey][0].append(
                    np.loadtxt(open(os.path.join(inDir, fname), 'rb'))
                )
                data[subkey][0][-1][:,0] += factor*shift + factorPt*shiftPt # shift data for visibility
                data[subkey][1].append(
                    'lc %s lw 4 ps 2 pt %d' % (default_colors[10+p],18+ipt)
                )
                pttit = '%g - %g GeV/c' % (pt-0.15, pt+0.15)
                data[subkey][2].append(pttit if particle == 'e' else '')
                factorPt += 1
            factor += 1
    fake_point = np.array([[-1,-1,0,0,0]])
    for p in range(2):
        data['negative particles'][0].append(fake_point)
        data['negative particles'][1].append('with lines lc %s lw 4 lt 1' % default_colors[10+p])
        data['negative particles'][2].append('{/Symbol \160}^{+/-}' if p else 'e^{+/-}')
    dx = 0.5*3+4*0.3
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'syserr'),
        yr = [0,35], xr = [16,66],
        xlabel = '{/Symbol \326}s_{NN} (GeV)',
        ylabel = 'syst. {/Symbol \504}{/Symbol \145}/{/Symbol \145} with {/Symbol \145} = {/Symbol \145}_{qual} {/Symbol \327} {/Symbol \145}_{glDCA} (%)',
        layout = '2x1', size = '4.5in,8in', tmargin = 0.9,
        key = ['nobox', 'width -3', 'at screen 1.0,1.0', 'maxrows 2'],
        key_subplot_id = 0, gpcalls = [
            'object %d rectangle back fc rgb "#C6E2FF" from %f,%f to %f,%f' % (
                50+e, float(energy)-dx, 1, float(energy)+dx, 32
            ) for e,energy in enumerate(energies)
        ] + [
            'label %d "{/Helvetica-Bold=18 %s GeV}" at %f,%f rotate by 90 center' % (
                55+e, energy, float(energy)+0.5*dx, 28
            ) for e,energy in enumerate(energies)
        ],
    )

def gp_tpc_select_eff():
    data = OrderedDict()
    inDir, outDir = getWorkDirs()
    for energy in ['19', '27', '39', '62']:
        infile = os.path.join(inDir, 'tpc_select_eff', 'electrons_%sGeV.dat' % energy)
        data_import = np.loadtxt(open(infile, 'rb'))
        nrows = len(data_import)
        data_import[:,1:] *= 100. # convert to %
        #if energy != '19': data_import[:,2:] = 0
        key = '%s GeV' % getEnergy4Key(energy)
        data[key] = np.c_[
            data_import[:,:2], np.zeros(nrows),
            np.zeros(nrows), data_import[:,-1]
        ]
    make_plot(
        data = data.values(), titles = data.keys(),
        properties = [
            'with filledcurves lt 1 lc %s lw 5 pt 0' % default_colors[i]
            for i in range(len(data))
        ],
        tmargin = 0.98, rmargin = 0.99, yr = [77,100], xr = [0.2,2.052],
        gpcalls = ['xtics 0.5', 'mxtics 5', 'mytics 2'], key = ['nobox'],
        xlabel = 'momentum, p (GeV/c)', ylabel = 'TPC Selection Efficiency (%)',
        name = os.path.join(outDir, 'tpc_select_eff'), size = '8.8in,6.8in',
        labels = {'{/Helvetica-Bold=22 Electrons}': [1.0,93,True]}
    )

def gp_tof_match():
    data = OrderedDict()
    inDir, outDir = getWorkDirs()
    subtitles = [
        '-1 < {/Symbol \150} < 1, -180 < {/Symbol \146} < 180',
        '0 < {/Symbol \150} < 0.25, -180 < {/Symbol \146} < 180',
        '-1 < {/Symbol \150} < 1, 45 < {/Symbol \146} < 60',
        '0 < {/Symbol \150} < 0.25, 45 < {/Symbol \146} < 60',
    ]
    for isuf,suffix in enumerate([
        'Eta1_Phi1', 'Eta8_Phi1', 'Eta1_Phi24', 'Eta8_Phi24'
    ]):
        subkey = subtitles[isuf]
        data[subkey] = [[], [], []]
        d = OrderedDict()
        for ip,particle in enumerate(['pi', 'e']):
            infile = os.path.join(inDir, 'tof_match', '%sminus_39_%s.dat' % (particle, suffix))
            data_import = np.loadtxt(open(infile, 'rb'))
            data_import[:,3:] *= 100. # convert to %
            data_import = data_import[data_import[:,0]<2.]
            if particle == 'pi': data_import[:,4] = 0
            if suffix != 'Eta8_Phi24': data_import[:,5] = 0
            nrows = len(data_import)
            d[particle] = np.c_[ data_import[:,0], data_import[:,3], np.zeros(nrows), data_import[:,4:] ]
            data[subkey][0].append(d[particle])
            data[subkey][1].append('with %s lc %s lw 5 lt 1' % (
                'points pt 18 ps 1.5' if particle == 'e' else 'filledcurves pt 0', default_colors[ip]
            ))
            data[subkey][2].append('e^{-} (39 GeV)' if particle == 'e' else 'scaled {/Symbol \160}^{-} (39 GeV)')
        chi2i = [
            ((o-e)/s)**2
            for o,e,s in np.c_[d['e'][:,1], d['pi'][:,1], d['e'][:,3]]
            if o > 0 and e > 0 and s > 0
        ]
        newsubkey = subkey + ', {/Symbol \543}@^{2}_{red} = %.2g' % (sum(chi2i)/len(chi2i))
        data[newsubkey] = data[subkey]
        del data[subkey]
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'tof_match'),
        xlog = True, xr = [0.18, 2.1], yr = [35, 97],
        xlabel = 'transverse momentum, p_{T} (GeV/c)',
        ylabel = 'TOF Matching Efficiency (%)',
        layout = '2x2', size = '5in,7.5in',
        key = ['nobox', 'bottom right'],
        gpcalls = [ 'xtics add (.2,.5,1,2)' ],
    )

def gp_tof_match_extra():
    inDir, outDir = getWorkDirs()
    data_import = np.loadtxt(open(os.path.join(inDir, 'tof_match', 'extra.dat'), 'rb'))
    energies = ['19.6 GeV', '27 GeV', '39 GeV', '62.4 GeV']
    particles = ['e^{-}', 'e^{+}']
    columns = ['{/Symbol \104}{/Symbol \145}', 'F', '{/Symbol \143}^{2}']
    NE, NP, NC = len(energies), len(particles), len(columns)
    data = OrderedDict()
    for particle in particles:
        for column in columns:
            key = 'x = %s_{%s}' % (column, particle)
            data[key] = [[], [], []]
    for col,column in enumerate(data_import.T[2:]):
        eidx, cidx = col/NP/NC, col%NC
        key = 'x = %s_{%s}' % (columns[cidx], particles[(col/NC)%NP])
        mean, std = np.mean(column), np.std(column)
        xmin, xmax = mean-3*std, mean+3*std
        trunc_column = np.array([x for x in column if x > xmin and x < xmax])
        bins = int(math.ceil(math.sqrt(len(trunc_column))))
        histo = list(np.histogram(column, bins=bins, range=(xmin,xmax)))
        binwidths = np.diff(histo[1])
        data[key][0].append(np.c_[
            0.5*(histo[1][:-1]+histo[1][1:]), histo[0],
            0.5*binwidths, np.zeros(bins), np.zeros(bins)
        ])
        data[key][1].append('with histeps lc %s lt 1 lw 3' % default_colors[eidx])
        data[key][2].append(energies[eidx])
    keys = data.keys()
    for key in keys:
        Ns, means, sigmas = [], [], []
        for d in data[key][0]:
            d[:,1] /= sum(d[:,1])
            Ns.append(len(d))
            mean = np.average(d[:,0], weights=d[:,1])
            sigma = np.sqrt(np.average((d[:,0]-mean)**2, weights=d[:,1]))
            means.append(mean)
            sigmas.append(sigma)
        Ns, means, sigmas = np.array(Ns), np.array(means), np.array(sigmas)
        glmean = sum(Ns*means)/sum(Ns)
        glsigma = sum(Ns*sigmas)/sum(Ns)
        for d in data[key][0]:
            d[:,0] -= glmean
            d[:,(0,2)] /= glsigma
            d[:,1] /= 2*d[:,2]
        newkey = ', '.join([
            key, '{/Symbol \155} = %.3g' % glmean, '{/Symbol \163} = %.2g' % glsigma
        ])
        data[newkey] = data[key]
        del data[key]
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'tof_match_extra'),
        xr = [-3.5,6.5], yr = [0,.98],
        xlabel = '(x - {/Symbol \155}) / {/Symbol \163}', ylabel = 'dN / N dx',
        layout = '3x2', size = '5in,8in',
        key = ['nobox', 'at graph 0.99,0.8'],
        gpcalls = ['bars small', 'ytics 0.2'],
        key_subplot_id = 2
    )

def gp_total():
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    particles = ['eplus', 'eminus']
    data = OrderedDict()
    for energy in energies:
        ekey = ' '.join([getEnergy4Key(energy), 'GeV'])
        data[ekey] = [[], [], []]
        for p,particle in enumerate(particles):
            data_import = np.loadtxt(
                open(os.path.join(inDir, 'total', particle+energy+'.dat'), 'rb')
            )
            data_import[:,3:] *= 100.
            if p == 1: data_import[:,0] *= 1.03
            data[ekey][0].append(np.c_[
                data_import[:,0], data_import[:,3],
                np.zeros(len(data_import)), data_import[:,4:]
            ])
            data[ekey][1].append('lt 1 lc %s lw 3 pt 18 ps 1.3' % default_colors[p])
            data[ekey][2].append('e^{-}' if particle == 'eminus' else 'e^{+}')
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'total_overview'),
        xlog = True, xr = [0.2, 3.6], yr = [10, 75],
        xlabel = 'transverse momentum, p_{T} (GeV/c)',
        ylabel = 'Total Single Track Efficiency (%)',
        layout = '2x2', size = '5in,7.5in',
        gpcalls = [ 'xtics add (.2,.5,1,2,3,5)'],
        key = ['bottom left', 'width 1.5'],
        #labels = { '0 < {/Symbol \550} < 0.25, 45 < {/Symbol \556} < 60': [0.7, 1.0, False]}
    )

def gp_pair():
    inDir, outDir = getWorkDirs()
    energies = ['19', '27', '39', '62']
    data = OrderedDict()
    data['19.6 GeV'] = [[], [], []]
    data['19.6 GeV'][0].append(np.array([[0]*5]))
    data['19.6 GeV'][1].append('with filledcurves lt 1 lc rgb "black" lw 3 pt 0')
    data['19.6 GeV'][2].append('{/Symbol \545\661\104\545}_{stat}')
    data['19.6 GeV'][0].append(np.array([[0]*5]))
    data['19.6 GeV'][1].append('with lines lt 2 lc rgb "black" lw 3')
    data['19.6 GeV'][2].append('{/Symbol \104\545}_{syst} / {/Symbol \545}')
    for energy in energies:
        ekey = ' '.join([getEnergy4Key(energy), 'GeV'])
        if ekey != '19.6 GeV': data[ekey] = [[], [], []]
        data_import = np.loadtxt(
            open(os.path.join(inDir, 'pair', 'pair'+energy+'.dat'), 'rb')
        )
        centers = np.array(sorted(set(data_import[:,1])))
        edges = [0.]
        for i,center in enumerate(centers):
            edges.append(edges[-1]+2*(center-edges[-1]))
        for i,(pt,g) in enumerate(groupby(data_import, itemgetter(1))):
            d = np.array(list(g))
            d[:,2:] *= 100.
            data[ekey][0].append(np.c_[
                d[:,(0,2)], np.zeros(len(d)), np.zeros(len(d)), d[:,3]
            ])
            data[ekey][1].append('with filledcurves lt 1 lc %s lw 3 pt 0' % default_colors[i])
            data[ekey][2].append('%g - %g GeV/c' % (edges[i], edges[i+1]))
            data[ekey][0].append(np.c_[
                d[:,0], d[:,4]/d[:,2]*100., np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
            ])
            data[ekey][1].append('with lines lt 2 lc %s lw 3' % default_colors[i])
            data[ekey][2].append('')
    make_panel(
        dpt_dict = data,
        name = os.path.join(outDir, 'pair'),
        xlog = True, xr = [0.002, 3.35], yr = [0, 33],
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = 'pair efficiency and systematic uncertainty (%)',
        layout = '2x2', size = '5in,7.5in', tmargin = 0.9,
        key = ['nobox', 'at screen 1.0,1.0', 'maxrows 2', 'width -3'],
    )

def gp_phiv():
    inDir, outDir = getWorkDirs()
    infile = os.path.join(inDir, 'phiVeff.dat')
    data_import = np.loadtxt(open(infile, 'rb'))
    nrows = len(data_import)
    data_import[:,1] *= 100. # convert to %
    data = np.c_[ data_import, np.zeros(nrows), np.zeros(nrows), np.zeros(nrows) ]
    make_plot(
        data = [data], titles = [''],
        properties = [ 'with linespoints lt 1 lc %s lw 3 pt 18 ps 1.7' % default_colors[0] ],
        tmargin = 0.98, rmargin = 0.99, bmargin = 0.16,
        yr = [70,103], xr = [0,0.31],
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = '{/Symbol \146}_{V} Efficiency (%)', gpcalls = ['nokey'],
        name = os.path.join(outDir, 'phiv'), size = '8.8in,6.8in',
        lines = {'y=0.2': 'lw 4 lt 2 lc 0'}
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
  #gp_syserr()
  #gp_tpc_select_eff()
  #gp_tof_match()
  #gp_tof_match_extra()
  #gp_total()
  #gp_pair()
  gp_phiv()
