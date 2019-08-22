import numpy as np
import os
import subprocess as subp
import time
import glob
from astropy.io import fits as pyfits
import scipy.optimize as opt
import scipy.stats as st
from itertools import groupby
import matplotlib.pylab as plt
from multiprocessing import Pool
import re

import plag
import aztools as az

from javelin.zylc import get_data
from javelin.lcmodel import Cont_Model, Rmap_Model

def fit_xspec_model(fit_func, spec_ids, base_dir, suff='', **kwargs):
    '''Call fit_func from fit.tcl to model the spectra in spec_ids
        and read the fit parameters. If the fitting has already been done,
        just read it.
    
    Parameters:
        fit_func: tcl function in fit.tcl to all; 
            It should be of the form proc fit_2a {sfile {suff ""}} {...}
        spec_ids: list of spec ids so the files are: spec_$id.grp
        base_dir: directory containing the fit.tcl file
        suff: any extra suffix for the saved fit files. 
            The saved files will be: fits/{fit_func}{suff}__{ispec}
     
     Keywords:
        read_fit: read fit result? use False when errors are not needed, so
            no log files exist. Default: True
        spec_root: root name for the spectra. Default: spec_%d.grp
        ext_check: file extention to use when checking whether the fit has already run
            or not. Default: xcm 
        extra_args: extra arguments to fit_func; default nothing
            
    Returns: an array of the fit parameters of shape: (nspec, npar, 4(par, err, err+, err-))
    '''    
    read_fit = kwargs.get('read_fit', True)
    spec_root = kwargs.get('spec_root', 'spec_%d.grp')
    ext_check = kwargs.get('ext_check', 'xcm')
    extra_args = kwargs.get('extra_args', '')
           
    
    procs = []
    for ispec in spec_ids:
        # if fit already done, skip
        if os.path.exists('fits/%s%s__%d.%s'%(fit_func, suff, ispec, ext_check)): continue
        
        tcl  = 'source %s/fit.tcl\n'%base_dir
        tcl += '%s %s %s__%d %s\nexit\n'%(fit_func, spec_root%ispec, suff, ispec, extra_args) 
        xcm = 'tmp_%d.xcm'%ispec
        with open(xcm, 'w') as fp: fp.write(tcl)
        cmd = 'xspec - %s > tmp_%d.log 2>&1'%(xcm, ispec)
        time.sleep(0.1)
        p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(p)
        if len(procs)==30:
            for p in procs: p.wait()
            procs = []
    # wait for the tasks to end
    for p in procs: p.wait()
    _ = os.system('rm tmp_*.??? >/dev/null 2>&1')
    
    # read the fit #
    if not read_fit: return
    fit = []
    for ispec in spec_ids:
        # exception for missing mos-1 data
        if not os.path.exists('fits/%s%s__%d.log'%(fit_func, suff, ispec)):
            print('missing fits/%s%s__%d.log'%(fit_func, suff, ispec))
            fit.append(fit[-1]*np.nan)
            continue
        fit.append(np.loadtxt('fits/%s%s__%d.log'%(fit_func, suff, ispec), usecols=[0,1,2,3]))
    return np.array(fit)

def ftest(c2, d2, c1, d1):
    """Do F-test"""
    fstat = ((c1-c2)/(d1-d2)) / (c2/d2)
    fprob = st.f.cdf(fstat, d1-d2, d2)
    fpval = 1 - fprob
    return fstat, fprob, fpval


def write_resid(base_dir, spec_ids, suff, extra_cmds='', avg_iref=-1, avg_bin=True,
                     outdir='results', z=3.319e-3, ispec_suff=''):
    """Plot the residuals from fits of the form fit_{suff}__{ispec}
    
    spec_ids: [1,2,3,...]
    suff: e.g. indiv_1l, 2a etc
    extra_cmds: any extra xspec commands between loading the data and plotting.
        e.g. removing cflux and renormalizing.
    avg_iref: reference ispec when averaging; -1, select the first even array
    ispec_suff: what to add to spec_ids (e.g. a, b for nustar)
    outdir: output directory
    """
    os.system('mkdir -p %s'%outdir)
    outfile = '%s/fit_%s.plot'%(outdir,suff)
    # individual fits #
    for ispec in spec_ids:
        tcl  = 'source %s/fit.tcl\nsetpl ener\nsetpl redshift %g\n'%(base_dir, z)
        tcl += '@fits/fit_%s__%d%s.xcm\n%s\n'%(suff, ispec, ispec_suff, extra_cmds)
        tcl += 'fit 1000\naz_plot_unfold u tmp_%d%s %s__%d%s 1 1\n'%(
                    ispec, ispec_suff, suff, ispec, ispec_suff)
        with open('tmp.xcm', 'w') as fp: fp.write(tcl)
        cmd = 'xspec - tmp.xcm > tmp.log 2>&1'
        p = subp.call(['/bin/bash', '-i', '-c', cmd])
    os.system("ls -1 tmp_*.plot|sort -t'_' -n -k2| xargs cat > %s/fit_%s.plot"%
              (outdir, suff))
    _ = os.system('rm tmp.??? tmp_*plot')
    
    # average residuals #
    # read and group the spectral data #
    lines = open(outfile).readlines()
    grp = [list(v) for k,v in groupby(lines, lambda l: (len(l)==0 or l[0]=='d') )]
    grp = [np.array([x.split() for x in g if x!='\n'], np.double) for g in grp if len(g)>4]

    dat = grp[::4]
    mod = grp[1::4]
    mod_spec = []
    for m,d in zip(mod, dat):
        mod_spec.append(np.array(d))
        mod_spec[-1][:,3] = 0
        mod_spec[-1][:,2] = m[:,0]

    # choose  some reference grid #
    iref = avg_iref
    if iref == -1:
        ilen = [i for i,d in enumerate(dat) if len(d)%2==0]
        iref = ilen[0]
    egrid_ref = np.array([dat[iref][:,0]-dat[iref][:,1], dat[iref][:,0]+dat[iref][:,1]]).T
    if avg_bin:
        egrid_ref = np.array([egrid_ref[::2,0], egrid_ref[1::2,1]]).T

    # map of spectra and models to the same reference grid #
    en_new, spec_new = spec_common_grid(dat, egrid_ref)
    _, mod_new = spec_common_grid(mod_spec, egrid_ref)

    # calculate residuals #
    dat_tot = np.array([ np.mean(spec_new[:,0], 0), 
                        (np.sum(spec_new[:,1]**2, 0)**0.5)/len(spec_new)])
    mod_tot = np.mean(mod_new[:,0], 0)
    #return en_new, spec_new, mod_new, dat_tot, mod_tot
    del_tot = (dat_tot[0]-mod_tot) / dat_tot[1]
    rat_tot = [dat_tot[0]/mod_tot, dat_tot[1]/mod_tot]

    # update the results file #results/fit_indiv_1.plot
    text = '\n\ndescriptor en_%s__tot,+- del_%s__tot,+- rat_%s__tot,+-\n'%tuple([suff]*3)
    text += '\n'.join(['{:.5} {:.5} {:.5} 1.0 {:.5} {:.5}'.format(*z) for z in 
                       zip(en_new[0], en_new[1], del_tot, rat_tot[0], rat_tot[1])])
    
    # add binned individual spectra #
    txt1 = ' '.join(['del_%s__g%d%s,+- rat_%s__g%d%s,+-'%(
                    suff, ispec, ispec_suff, suff, ispec, ispec_suff) 
                        for ispec in spec_ids])
    txt2 = '\n'.join([' '.join(['{:.5} 1.0 {:.5} {:.5}'.format(
            (spec_new[ispec,0,ie]-mod_new[ispec,0,ie])/spec_new[ispec,1,ie],
            spec_new[ispec,0,ie]/mod_new[ispec,0,ie], spec_new[ispec,1,ie]/mod_new[ispec,0,ie]
        )
        for ispec in range(len(spec_ids))]) for ie in range(len(en_new[0]))])
    text += '\n\ndescriptor ' + txt1 + '\n' + txt2
    with open(outfile, 'a') as fp: fp.write(text)



def spec_common_grid(spectra, egrid):
    """Map a list of spectra into a single energy grid
    spectra: a list of (x,xe,y,ye) spectra
    egrid: array(nen, 2) giving low/high bin boundaries
    
    Returns: en_new(2,nen), spectra_new(nspec,2,nen)
    """
    nspec = len(spectra)
    spectra_new = np.zeros((nspec, 3, len(egrid)))
    for ispec in range(nspec):
        # the energy grid to be aligned #
        ebins = np.array([spectra[ispec][:,0]-spectra[ispec][:,1],
                          spectra[ispec][:,0]+spectra[ispec][:,1]]).T
        # loop through elements of the reference grid #
        # for each, sum the bins from the other grids that
        # fall into the current bin of the reference
        for iref, eref in enumerate(egrid):
            for ibin, ebin in enumerate(ebins):
                if ebin[1] < eref[0] or ebin[0] > eref[1]:
                    continue
                erange = (np.min([eref[1], ebin[1]]) - np.max([eref[0], ebin[0]]))
                efrac  = erange / (ebin[1] - ebin[0])
                spectra_new[ispec, 0, iref] += efrac * spectra[ispec][ibin,2]
                spectra_new[ispec, 1, iref] += (efrac * spectra[ispec][ibin,3])**2
                spectra_new[ispec, 2, iref] += efrac
    # spectra_new[:,2] is the fractional exposure in each bin
    spectra_new[:,0] /= spectra_new[:,2]
    spectra_new[:,1] = (spectra_new[:,1]/spectra_new[:,2])**0.5
    spectra_new = spectra_new[:,:2,]
    en_new = np.array([np.sum(egrid, 1)/2, np.diff(egrid, 1)[:,0]/2])
    return en_new, spectra_new



def fit_linear_model(X, Y, name, spec_ids=None, suff=''):
    '''Fit a linear model to two variables. The uncertainty in the variables
    is accounted for through mcmc
    
    Parameters:
        X,Y: the two parameters each has shape: (nchain, npoints)
        name: name of parameters being studies for printing: e.g. pf_gf
        spec_ids: a list of spec id's for printing corresponding to npoints.
            if None, use range(len(npoints))
        suff: suffix for parameter names when printing. _ will be added
    '''
    if spec_ids is None: spec_ids = [i for i in range(len(X[0]))]
    if suff != '': suff = '__%s'%suff
    
    # correlation test
    npoints = len(X[0])
    sr = np.median([st.spearmanr(x,y)[0] for x,y in zip(X, Y)])
    ttest = sr*np.sqrt((npoints-2)/(1-sr*sr))
    tprob = 1-st.t.cdf(ttest, df=npoints-2)
    
    # fit each set in the chain separately #
    def fit_f(x, *p): return p[0]*x + p[1]
    fit = np.array([opt.curve_fit(fit_f, x, y, [0.4,-7])[0]
                    for x,y in zip(X, Y)])
    
    # model for plotting #
    Xm = np.median(X, 0)
    xmod = np.linspace(np.min(Xm), np.max(Xm),40)
    ymod = np.array([fit_f(xmod, *_x) for _x in fit])
    ymod = [np.median(ymod, 0), ymod.std(0)]
    xdat = [np.median(X, 0), X.std(0)]
    ydat = [np.median(Y, 0), Y.std(0)]
    fit = [np.median(fit, 0), fit.std(0)]
    
    # output text #
    text = '\n# {}: spearman r,pvalue: {:8.4} {:8.4}\n'.format(name, sr, tprob)
    print(text)
    text += '# fit pars: {:8.4} +- {:8.4}, {:8.4} +- {:8.4}\n'.format(
                fit[0][0], fit[1][0], fit[0][1], fit[1][1])
    
    text += '\ndescriptor {0}_x{1} {0}_y{1},+-\n'.format(name, suff)
    text += '\n'.join(['{:8.6} {:8.6} {:8.6}'.format(
        xmod[i], ymod[0][i], ymod[1][i]) for i in range(len(xmod))])
    text += '\ndescriptor id{1} {0}_xd{1},+- {0}_yd{1},+-\n'.format(name, suff)
    text += '\n'.join(['{:4} {:8.6} {:8.6} {:8.6} {:8.6}'.format(spec_ids[i],
        xdat[0][i], xdat[1][i], ydat[0][i], ydat[1][i]) for i in range(len(xdat[0]))])
    
    return xdat, ydat, xmod, ymod, text

def model_correlations(froot, spec_ids, ipars, names, ilog=[], plot=True, read_only=False):
    """read mcmc of parameters and model correlations by calling fit_linear_model on each
    correlation
    
    Parameters:
        froot: root for the fit. e.g. 'fit_2'
        spec_ids: list of spectral id's
        ipars: a list of pairs of parameter indices to model
        names: name of the parameter ids for plotting, corresponds to ipars
        ilog: a list of parameters to log. e.g sigma
        plot?
        read_only: read and return the chain data only. Default False.
            
    """
    # read mcmc data #
    fdata = []
    for ispec in spec_ids:
        d = pyfits.getdata('fits/%s_mc__%d.fits'%(froot, ispec))
        fdata.append([d.field(i) for i in d.names])
    # nspec, npar, nmcmc
    imax = np.max(ipars)+1
    fdata = np.array(fdata)[:,:imax,:]
    
    if read_only:
        return fdata
    
    nspec, npar, nsim = fdata.shape
    
    # log parameters #
    for i in ilog:
        fdata[:,i] = np.log10(fdata[:,i])
    
    if plot: fig = plt.figure(figsize=(12,3))
    
    text = ''
    for ip in range(len(ipars)):
        print('doing %s'%names[ip])
        X, Y = fdata[:, ipars[ip][0], :].T, fdata[:, ipars[ip][1], :].T
        xdat, ydat, xmod, ymod, txt = fit_linear_model(X, Y, names[ip], spec_ids, froot)
        text += txt

        # plot #
        if not plot: continue
        ax = plt.subplot(1, len(ipars), ip+1)
        plt.errorbar(xdat[0], ydat[0], ydat[1], xerr=xdat[1], fmt='o', ms=8, lw=0.5)
        plt.fill_between(xmod, ymod[0]-ymod[1], ymod[0]+ymod[1], alpha=0.5)
        xl,yl = names[ip].split('_')
        ax.set_xlabel(xl); ax.set_ylabel(yl)
    if plot: plt.tight_layout(pad=0)
    return text

def continuum_line_split_text(text, suff):
    """Split the veusz text for the continuum line into groups like old/new 
    and those separated by less than two days
    
    Parameters:
        text: veusz text resulting from running correlation modeling.   
        suff: the suffix used in test
    """

    grp_ids = [[1,2,3], [4,5,6], [16,17,18], [19,20], [22,23], 
              [1,2,3,4,5,6,7,8,10,11,12,13,15], [16,17,18,19,20,21,22,23,24]]
    grp_lab = ['s1', 's2', 's3', 's4', 's5', 'O', 'N']
    lines = text.split('\n')
    do_seg, data, desc = False, [], []
    for l in lines:
        if len(l) == 0 or l[0] == '#': continue
        if l[0] == 'd':  
            do_seg = 'xd' in l
            if not do_seg: data.append([])
            if do_seg: desc.append(l)
            continue
        if not do_seg: continue
        data[-1].append(l.split())
    text_sep = '\n\n'
    for dsc,dat in zip(desc,data):
        for ig, gl in zip(grp_ids, grp_lab):
            g = [d for d in dat if np.int(d[0]) in ig]
            text_sep += '\n%s\n'%dsc.replace('__%s'%suff, '_%s__%s'%(gl, suff))
            text_sep += '\n'.join([' '.join(x) for x in g])
    return text_sep






def javelin_modeling(time, cflux, lflux, suff='', nburn=1000, nchain=500, mods=None):
    """Do lag calculations with Javelin
    
    Parameters:
        time: [n] array of time axis
        cflux: [2,n] continuum flux and error
        lflux: [2,n] line flux and error
        suff: suffix for printing
        mods: javelin models if already calculated [cont_mod, cont_line_mod]. 
            In that case, this function just does the plotting. 
    """
    
    if mods is None:
        # write light curves to file so javelin can read them #
        irand = np.random.randint(100000)
        text = '\n'.join(['{} {} {}'.format(*x) for x in zip(time, cflux[0], cflux[1])])
        with open('tmp_c%d.dat'%irand, 'w') as fp: fp.write(text)
        text = '\n'.join(['{} {} {}'.format(*x) for x in zip(time, lflux[0], lflux[1])])
        with open('tmp_l%d.dat'%irand, 'w') as fp: fp.write(text)

        # continuum first #
        cont_lc = get_data(['tmp_c%d.dat'%irand])
        cont_mod = Cont_Model(cont_lc)
        cont_mod.do_mcmc(set_verbose=False)
        #cont_mod.show_hist();return 0,cont_mod,0

        # continuum and line #
        cont_line_lc = get_data(['tmp_c%d.dat'%irand, 'tmp_l%d.dat'%irand]) 
        cont_line_mod = Rmap_Model(cont_line_lc)



        llimit = [-100, 100]
        pmap,_ = cont_line_mod.do_map([-0.6, 1.5, 1.0, 0.1, .2], 
                                      fixed=[1,1,1,1,1], set_verbose=False)
        cont_line_mod.do_mcmc(conthpd=cont_mod.hpd, laglimit=[llimit], nburn=nburn, nchain=nchain,
                         #fixed=[1,1,1,1,1], p_fix=pmap, threads=20,set_verbose=False)
                         threads=30,set_verbose=False)

        os.system('rm tmp_c%d.dat tmp_l%d.dat'%(irand, irand))
    else:
        cont_mod, cont_line_mod = mods
    
    # plot the result  of javelin fit #
    chains = cont_line_mod.flatchain
    lag_javelin = plt.histogram(chains[:,2], 400, density=1)

    bins_cent = (lag_javelin[1][1:] + lag_javelin[1][:-1])/2
    bins_err  = (lag_javelin[1][1:] - lag_javelin[1][:-1])/2

    percentile = lambda l: '[{:.4}, {:.4}]'.format(
        *np.percentile(chains[:,2], [(100-l)/2, l+(100-l)/2]))
    text = '# percentiles: 68, 90, 99: {}, {}, {}'.format(*[percentile(x) for x in [68, 90, 99]])
    text += '\n# mean lag: {:.4} +{:.4} {:.4}'.format(
        np.percentile(chains[:,2],50),
        np.percentile(chains[:,2],68+16)-np.percentile(chains[:,2],50),
        np.percentile(chains[:,2],16)-np.percentile(chains[:,2],50))
    print(text)
    text += '\ndescriptor lag_javelin{0},+- lag_javelin_prob{0}\n'.format(suff)
    text += '\n'.join(['{:.4} {:.4} {:.4}'.format(*x) 
                       for x in zip(bins_cent, bins_err, lag_javelin[0])])
    return text, cont_mod, cont_line_mod, lag_javelin



def estimate_psd(xtime, xarr, xerr, mod, p0=None, **kwargs):
    """Fit a light curve with plag
    
    xtime,xarr,xerr: time, flux, flux_err
    mod: 1 for powerlaw, 2 for bending powerlaw
    """
    isort = np.argsort(xtime)
    xtime, xarr, xerr = xtime[isort], xarr[isort], xerr[isort]
    
    fqL = [0.1/(xtime[-1]-xtime[0]), 2./np.diff(xtime).min()]
    #print(fqL)
    pm = plag.PLag('psdf', [xtime], [xarr], [xerr], 0.1, fqL, 'rms', 0, mod)
    if p0 is None:
        p0 = [-9, -2] if mod==1 else [-9.4,-2, -3.6]
    p,pe = plag.optimize(pm, p0, **kwargs)[:2]
    return p,pe,pm

def calc_scatter(x, xe, y, ye, nsim=100):
    """Fit y~x with a linear model and calculate the scatter
    Uncertainties are handled with mcmc
    """
    def fit_func(x, *p): return x*p[0] + p[1]
    nx = len(x)
    if nsim==0:
        X,Y = [x], [y]
    else:
        X = np.random.randn(nsim, nx) * xe + x
        Y = np.random.randn(nsim, nx) * ye + y
    fit = np.array([opt.curve_fit(fit_func, _x, _y, [0.4, -7])[0] for _x,_y in zip(X,Y)])
    scat = np.array([np.sum((y - fit_func(x, *f))**2) for f in fit])
    scat_e = np.mean(0.5 * (scat**2 * np.sqrt(2/(nx-1))) / scat)
    # here we include both systematic (stochastic) and statistical (measurement) uncertainties.
    return scat.mean(), (scat.std()**2 + scat_e**2)**0.5


def simulate_scatter(lag, xtime, cflux, lflux, psdpar, dt, lag_model='const', 
                     nsim=200, tf_width_frac=0.5):
    """Simulate the scatter in cflux-lflux relation
    
    Parameters:
        lag: lag value to simulate
        xtime, cflux,lflux: light curves of dims [2,nval]
        psdpar: input continuum psd parameters for az.SimLC
            shape of (2,2) for poweralw and (2,3) for bending powerlaw.
            the first dim is for value,err
        dt: time sampling
        lag_model: 'const' for a delta function otherwise use top-hat with tf_width_frac
            controlling the fractional half width
    
    """
        
    xe, ye = cflux[1], lflux[1]
    mu_c, mu_l = cflux[0].mean(), lflux[0].mean()

    xtime = np.array(xtime/dt, np.int) * dt
    gtime = np.arange(np.int((xtime[-1]-xtime[0])/dt)+10) * dt + xtime[0]
    itime = np.in1d(gtime, xtime)
    assert(len(itime[itime]) == len(xtime))
    ngtime = len(gtime)
    len_sim = 2**(1+np.int(np.log(ngtime)/np.log(2)))
    
    
    scat = []
    for isim in range(nsim):
        sim = az.SimLC()
        psdpar_sim = np.random.randn(len(psdpar[0])) * psdpar[1] + psdpar[0]
        psdpar_sim[0] = np.exp(psdpar_sim[0])
        if len(psdpar_sim) == 2:
            # powerlaw psd fit #
            sim.add_model('broken_powerlaw', [psdpar_sim[0], -1, psdpar_sim[1], 1e-7])
        else:
            # bending power psd fit #
            psdpar_sim[2] = np.exp(psdpar_sim[2])
            sim.add_model('bending_powerlaw', [psdpar_sim[0], psdpar_sim[1], psdpar_sim[2]])
        sim.simulate(len_sim, dt, mu_c, 'rms')
        
        if lag_model == 'const':
            # simple delta function #
            sim.add_model('constant', lag, clear=True, lag=True)
            sim.apply_lag(phase=False)
            y = sim.y
        else:
            # top hat 
            tf = np.zeros_like(sim.x)
            if tf_width_frac > 0:
                # tf_width_frac is used as a fraction of the lag #
                tf[(sim.t>=np.abs(lag)*(1-tf_width_frac)) & 
                   (sim.t<=np.abs(lag)*(1+tf_width_frac))] = 1
            else:
                if tf_width_frac==0: # randomize to marginalize
                    w = np.random.rand() * 5
                    #w = np.clip(w, dt/2, np.abs(lag)/2)
                else:
                    w = np.abs(tf_width_frac)/2 # simple fixed width
                    w = np.clip(w, 0, np.abs(lag)/2)
                tf[(sim.t>=np.abs(lag) - w) & (sim.t<=np.abs(lag) + w)] = 1
                
            #tf = st.norm.pdf(sim.t, loc=lag, scale=lag*tf_width_frac)
            if not np.all(tf==0):
                tf /= tf.sum()
            # tf2lag #
            n = len(sim.t)
            fq = np.arange(1,n/2+1)/(n*dt)
            ff = np.fft.fft(tf)[1:(n//2+1)]
            denom = 1+np.real(ff)
            denom[denom==0] = 1e-7
            tan_phi = np.imag(ff) / denom
            lag_ = -np.sign(lag)*np.arctan(tan_phi) / (np.pi*fq)
            lag_ = np.concatenate([[lag_[0]], lag_])
            fq = np.concatenate([[0], fq])
            y = sim.lag_array(sim.x, lag_, False, fq)
            
        X  = sim.x[:ngtime][itime]
        Y = (y[:ngtime][itime] - mu_c) * (mu_l/mu_c) + mu_l   
        X += np.random.randn(len(xe)) * xe
        Y += np.random.randn(len(ye)) * ye
        scat.append(calc_scatter(X, xe, Y, ye, nsim=0)[0])
    return np.array(scat)

def run_parallel_simulate_scatter(npz_file, args):
    """A wrapper to call simulate_scatter in parallel for multiple lag values
       
    Parameters:
        npz_file: to save data, and if present read it from there.
        args: a list of parameters to be passed to simulate_scatter
    
    """
    if not os.path.exists(npz_file):
        pool = Pool(35)
        l1 = np.linspace(0, 10, 31)[1:-1]
        l2 = np.linspace(10, 30, 21)
        lag = np.concatenate([-l2[::-1], -l1[::-1], [0.0], l1, l2])
        #lag = np.linspace(-30,30,101)
        dlag = lag[1:] - lag[:-1]
        lag = (lag[1:] + lag[:-1])/2
        args = [[l]+args for l in lag]
        sim_scat = np.array(pool.starmap(simulate_scatter, args))
        pool.close()
        np.savez(npz_file, data=[sim_scat, lag, args, dlag, []])
    else:
        sdata = np.load(npz_file)['data']
        sim_scat, lag, args, dlag, _ = sdata
    return sim_scat, lag, dlag