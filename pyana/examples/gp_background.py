import logging, argparse, re, os, glob
import numpy as np
from .utils import getWorkDirs, checkSymLink, getEnergy4Key
from ..ccsgp.ccsgp import make_plot, make_panel
from ..ccsgp.config import default_colors
from collections import OrderedDict

def _clamp(val, minimum = 0, maximum = 255):
    """convenience function to clamp number into min..max range"""
    if val < minimum: return minimum
    if val > maximum: return maximum
    return val

def _colorscale(gpcol, scalefactor = 1.4):
    hexstr = gpcol.split('#')[-1][:-1]
    if scalefactor < 0 or len(hexstr) != 6: return hexstr
    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)
    r = _clamp(r * scalefactor)
    g = _clamp(g * scalefactor)
    b = _clamp(b * scalefactor)
    return 'rgb "#%02x%02x%02x"' % (r, g, b)

def gp_background():
    """ plot background methods and S/B vs energy """
    inDir, outDir = getWorkDirs()
    data, REBIN = OrderedDict(), None
    titles = [ 'SE_{+-}', 'SE@^{corr}_{/Symbol \\261\\261}', 'ME@^{N}_{+-}' ]
    Apm = OrderedDict([
        ('19', 0.026668), ('27', 0.026554), ('39', 0.026816), ('62', 0.026726)
    ])
    fake = np.array([[-1, 1, 0, 0, 0]])
    lines = {'y=0.9': 'lc {} lt 2 lw 3'.format(default_colors[-4])}
    for energy in ['19', '27', '39', '62']:
        ekey = ' '.join([getEnergy4Key(energy), 'GeV'])
        data[ekey] = [[], [], []]
        for didx,dtype in enumerate(['epsPt', 'ngmPt_corr', 'epmPt']):
            for idx,infile in enumerate(glob.glob(os.path.realpath(os.path.join(
                inDir, 'rawdata', energy, 'pt-differential', '%s_*.dat' % dtype
            )))):
                file_url = os.path.realpath(os.path.join(inDir, infile))
                data_import = np.loadtxt(open(file_url, 'rb'))
                if REBIN is None: REBIN = int(data_import[-1][2]*2*1000) # MeV
                data_import[:,4] = data_import[:,3]
                data_import[:,(2,3)] = 0
                if dtype == 'ngmPt_corr':
                    data_import = data_import[data_import[:,0] <= 0.9]
                if dtype == 'epmPt':
                    data_import = data_import[data_import[:,0] > 0.9]
                    data_import[:,(1,4)] *= Apm[energy]
                col = _colorscale(default_colors[didx], 1.+idx*0.2)
                momrange =  os.path.basename(infile).split('_')[-1][:-4]
                if idx < 1:
                    data[ekey][0].append(fake)
                    data[ekey][1].append('with lines lt 1 lw 5 lc %s' % col)
                    data[ekey][2].append(titles[didx])
                data[ekey][0].append(data_import)
                data[ekey][1].append('lt 1 lw 3 lc %s pt 0' % col)
                data[ekey][2].append('')
    # unsubtracted background
    make_panel(
        name = '%s/methods' % outDir, dpt_dict = data,
        xr = [0,3.35], yr = [0.9,2e5], ylog = True, lines = lines,
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = 'counts / %d MeV/c^{2}' % REBIN, layout = '2x2',
        key = ['spacing 1.6', 'nobox'], gpcalls = ['boxwidth 0.002'],
        labels = {
            '{/=20 the lighter the color, the higher the p_{T}}': (1.5, 1e5)
        }, key_subplot_id = 1,
    )
    return 'done'
    # background ratio and acc.corr.
    make_plot(
        name = '%s/ratios%s' % (outDir, energy),
        xr = [0,1.6], yr = [0.95,1.2],
        data = graph_data[3:],
        properties = [
            'with filledcurves lt 1 lw 3 lc %s pt 0' % default_colors[i]
            for i in xrange(2)
        ],
        titles = [
            'SE_{/Symbol \\261\\261} / ME_{/Symbol \\261\\261}',
            'f_{acc} = ME_{+-} / ME_{/Symbol \\261\\261}'
        ],
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = '', key = [ 'width -2' ],
        labels = { '%s GeV' % energy: (0.4, 0.97) }
    )
    # signal-to-background ratio in rho/omega region vs. energy
    graph_data_sn = []
    for infile in os.listdir(os.path.join(inDir, 'sn')):
        energy = re.compile('\d+').search(infile).group()
        file_url = os.path.join(inDir, 'sn', infile)
        data_import = np.loadtxt(open(file_url, 'rb'))
        mask = (data_import[:,0] > 0.3) & (data_import[:,0] < 0.75)
        data_import = data_import[mask]
        weights = 1./data_import[:,3]
        sn = np.average(data_import[:,1], weights = weights)
        sn_err = np.average((data_import[:,1]-sn)**2, weights = weights)
        graph_data_sn.append([float(getEnergy4Key(energy)), sn, 0, sn_err, 0])
    graph_data_sn = np.array(graph_data_sn)
    make_plot(
        name = '%s/SNvsEnergy' % (outDir), xr = [15,210],
        yr = [1e-3, .11], xlog = True, ylog = True,
        data = [ np.array([[15,0.1,0,0,0],[210,0.1,0,0,0]]), graph_data_sn ],
        properties = [
            'with lines lt 2 lw 4 lc 0',
            'lt 1 lw 3 lc %s pt 18 ps 2' % default_colors[0]
        ],
        titles = ['']*2,
        xlabel = '{/Symbol \326}s_{NN} (GeV)',
        ylabel = 'S/B for 0.3 < M_{ee} < 0.75 GeV/c^{2}',
        lmargin = 0.1, gpcalls = [
            'nokey', 'format x "%g"',
            'xtics (20,"" 30, 40,"" 50, 60,"" 70,"" 80,"" 90, 100, 200)',
        ],
        labels = { 'p+p': (100, 0.09) }
    )
    return 'done'

def gp_rebin():
    MIL = 1e6
    NEVTS = { '19': 32.2307*MIL, '27': 63.6828*MIL, '39': 122.390*MIL, '62': 59.4631*MIL}
    inDir, outDir = getWorkDirs()
    data = OrderedDict()
    titles = [ 'raw signal', 'rebinned raw signal']
    lines = {'y=0.9': 'lc {} lt 2 lw 3'.format(default_colors[-4])}
    min_content, max_content = 1e20, -1e20
    colors = [default_colors[-9], default_colors[0]]
    points = [1,6]
    for eidx,energy in enumerate(['19', '27', '19', '27', '39', '62', '39', '62']):
        sgn_idx = (eidx%4)/2
        ekey = ' '.join([getEnergy4Key(energy), 'GeV'])
        if sgn_idx == 0: ekey += ' {}'.format(sgn_idx)
        data[ekey] = [[], [], []]
        for didx,dtype in enumerate(['sig', 'sigRb']):
            for idx,infile in enumerate(glob.glob(os.path.realpath(os.path.join(
                inDir, 'rawdata', energy, 'pt-integrated', '%s.dat' % dtype
            )))):
                if sgn_idx == 1 and dtype == 'sigRb': continue
                file_url = os.path.realpath(os.path.join(inDir, infile))
                data_import = np.loadtxt(open(file_url, 'rb'))
                data_import = data_import[data_import[:,0]>0.1]
                for i in [1,3,4]:
                    data_import[:,i] /= NEVTS[energy]
                    if dtype == 'sig': data_import[:,i] /= 2*data_import[:,2]
                if dtype == 'sig':
                    data_import[:,4] = data_import[:,3]
                    data_import[:,(2,3)] = 0
                if sgn_idx == 0: data_import = data_import[data_import[:,1]>0]
                else: data_import = np.abs(data_import[data_import[:,1]<0])
                cur_min, cur_max = min(data_import[:,1]), max(data_import[:,1]) 
                if ( cur_min < min_content ): min_content = cur_min
                if ( cur_max > max_content ): max_content = cur_max
                data[ekey][0].append(data_import)
                data[ekey][1].append('lt 1 lw %d lc %s pt %d' % (
                    3+didx, colors[didx], points[didx]))
                data[ekey][2].append(titles[didx])
    make_panel(
        name = '%s/rebin' % outDir, dpt_dict = data, ylog = True,
        xr = [0.1,3.2], yr = [2e-6, 7e-3], lines = lines, size = '7in,8in',
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        ylabel = '1/N@_{mb}^{evt} dN@_{ee}^{acc.}/dM_{ee} [ (GeV/c^2)^{-1} ]',
        layout = '2x4', key = ['width -5', 'nobox'], rmargin = 0.999, tmargin = 0.999,
        gpcalls = ['boxwidth 0.002', 'bars small', 'xtics (0.1,0.5,1,1.5,2,2.5,3)'],
    )
    return 'done'

def gp_norm(infile):
    """indentify normalization region"""
    inDir, outDir = getWorkDirs()
    data, titles = [], []
    for eidx,energy in enumerate(['19', '27', '39', '62']):
        file_url = os.path.realpath(os.path.join(
            inDir, 'rawdata', energy, 'pt-integrated', infile+'.dat'
        ))
        data_import = np.loadtxt(open(file_url, 'rb'))
        data_import[:,1] += eidx * 0.2
        data_import[:,4] = data_import[:,3]
        data_import[:,(2,3)] = 0
        data.append(data_import)
        titles.append(' '.join([getEnergy4Key(energy), 'GeV']))
    nData = len(data)
    lines = dict(
        ('x={}'.format(1+i*0.2), 'lc {} lt 2 lw 4'.format(default_colors[-2]))
        for i in range(nData)
    )
    lines.update(dict(
        ('x={}'.format(1+i*0.2+0.02), 'lc {} lt 3 lw 4'.format(default_colors[-5]))
        for i in range(nData)
    ))
    lines.update(dict(
        ('x={}'.format(1+i*0.2-0.02), 'lc {} lt 3 lw 4'.format(default_colors[-5]))
        for i in range(nData)
    ))
    lines.update({'y=0.9': 'lc {} lt 1 lw 4'.format(default_colors[-2])})
    make_plot(
        name = '%s/norm_range_%s' % (outDir,infile), xr = [0,2], yr = [0.9,1.7],
        data = data, properties = [
            'lt 1 lw 3 lc %s pt 1' % (default_colors[i]) # (i/2)%4
            for i in range(nData)
        ], titles = titles, size = '8in,8in',
        lmargin = 0.05, rmargin = 0.99, tmargin = 0.93, bmargin = 0.14,
        xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
        lines = lines, key = [
            'maxrows 1', 'nobox', 'samplen 0.1', 'width -1', 'at graph 1,1.1'
        ], labels = { 'SE_{--} / ME@_{--}^N': (0.3, 1.3) }, gpcalls = [
            'ytics (1,"1" 1.2, "1" 1.4, "1" 1.6)', 'boxwidth 0.002',
        ],
    )

def gp_acc():
    """acceptance correction"""
    inDir, outDir = getWorkDirs()
    for energy in ['19', '27', '39', '62']:
        data, titles = [], []
        for idx,infile in enumerate(glob.glob(os.path.realpath(os.path.join(
            inDir, 'rawdata', energy, 'pt-differential', 'acPt_*.dat'
        )))):
            data_import = np.loadtxt(open(infile, 'rb'))
            data_import[:,1] += idx * 0.2
            data_import[:,4] = data_import[:,3]
            data_import[:,(2,3)] = 0
            data.append(data_import)
            titles.append(os.path.basename(infile).split('_')[-1][:-4])
        nData = len(data)
        lines = dict(
            ('x={}'.format(1+i*0.2), 'lc {} lt 2 lw 4'.format(default_colors[-2]))
            for i in range(nData)
        )
        lines.update(dict(
            ('x={}'.format(1+i*0.2+0.02), 'lc {} lt 3 lw 4'.format(default_colors[-5]))
            for i in range(nData)
        ))
        lines.update(dict(
            ('x={}'.format(1+i*0.2-0.02), 'lc {} lt 3 lw 4'.format(default_colors[-5]))
            for i in range(nData)
        ))
        make_plot(
            name = '%s/accfac%s' % (outDir,energy), xr = [0,2], yr = [0.8,2],
            data = data, properties = [
                'lt 1 lw 3 lc %s pt 1' % (default_colors[i])
                for i in range(nData)
            ], titles = titles, size = '8in,8in',
            lmargin = 0.05, rmargin = 0.98, tmargin = 0.93, bmargin = 0.14,
            xlabel = 'dielectron invariant mass, M_{ee} (GeV/c^{2})',
            lines = lines, key = [
                'maxrows 1', 'nobox', 'samplen 0.1', 'width -2', 'at graph 1,1.1'
            ], labels = {
                'ME_{+-} / 2{/Symbol \326}ME_{++}ME_{--}': (0.3, 0.85),
                ' '.join([getEnergy4Key(energy), 'GeV']): (1.3, 0.85)
            },
            gpcalls = [ 'ytics (1,"1" 1.2, "1" 1.4, "1" 1.6, "1" 1.8)', 'boxwidth 0.002', ],
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
    #gp_background()
    #gp_norm('rmm')
    #gp_norm('rpp')
    #gp_acc()
    gp_rebin()
