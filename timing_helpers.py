import numpy as np
import os
import subprocess as subp
import time
import glob
from astropy.io import fits as pyfits
import scipy.stats as st
import matplotlib.pylab as plt

import aztools as az


def read_pn_lc(obs, dt, enL, lcdir, data_dir, interp=False, **kwargs):
    """Read PN light curves
    
    obs: a list of observation folders
    dt: time bin
    enL: a string of space-separated bin boundaries
    lcdir: directory name indicating the bin type. e.g 8l for 8 bins in log space 
    data_dir: folder containing the obsids
    interp: if true, interpolate small gaps
    kwargs for read_pn_lcurve
    """
    enL       = np.array(enL.split(), np.double)
    en, ene   = (enL[1:] + enL[:-1])/2, (enL[1:] - enL[:-1])/2
    Lc        = [[az.LCurve.read_pn_lcurve('{}/{}/pn/lc_r0.5/{}/lc_{:03g}__{}.fits'.format(
                        data_dir, o, lcdir, dt, ie+1), **kwargs)
                        for o in obs] for ie in range(len(en))]
    Lc = [[ll.make_even() for ll in l] for l in Lc]
    if interp:
        for l in Lc:
            for ll in l: ll.interp_small_gaps(np.int(1e3/dt), noise='norm')
    return Lc, en, ene


def lc_to_segments(Lc, seglen=10000, strict=False, **kwargs):
    """Split a Lc: a list of lists of LCurve to segments.
    No energy sum is done
    
    Lc: list of lists of LCurve: dims: nen, nobs; this is the output of read_pn_lc
    seglen, strict and kwargs: to be passed to az.misc.split_array
    
    Returns: rate_all, rerr_all, time_all, seg_indx. The first three are of dims: (nen, nseg_total),
        and seg_indx is a list of rate_all indices to which individual segments belong. 
    """

    nen, nobs = len(Lc), len(Lc[0])
    Lc = [[ll.make_even() for ll in l] for l in Lc]

    # arrays to be split, including rerr
    # arrs: (nobs, 2*nen(rate, rerr)); splt: (nobs, 2*nen, nseg)
    arrs = [[Lc[ie][io].rate for ie in range(nen)] + 
            [Lc[ie][io].rerr for ie in range(nen)] + 
            [Lc[ie][io].time for ie in range(nen)] for io in range(nobs)]
    splt = [az.misc.split_array(arr[0], seglen, strict, *arr, **kwargs)[2:] for arr in arrs]

    # segments: dim: (nseg_total, nen)
    rate_all = [s[:nen] for s in splt]
    rerr_all = [s[nen:(2*nen)] for s in splt]
    time_all = [s[(2*nen):] for s in splt]
    
    # indices for seprating the observations if needed #
    ic, seg_indx = 0, []
    for iar in splt:
        j = len(iar[0])
        seg_indx.append(list(range(ic, ic+j)))
        ic += j
        
    # flatten the lists so they have dims: (nseg_total, nen)
    rate_all = [[k for i in range(nobs) for k in rate_all[i][ie]] for ie in range(nen)]
    rerr_all = [[k for i in range(nobs) for k in rerr_all[i][ie]] for ie in range(nen)]
    time_all = [[k for i in range(nobs) for k in time_all[i][ie]] for ie in range(nen)]
    seg_indx = np.concatenate([[i]*len(s) for i,s in enumerate(seg_indx)])

    return rate_all, rerr_all, time_all, seg_indx


def simulate_like(x, dt, nsim=1, seed=None):
    """simulate nsim light curve arrays like x
    powerlaw parameters are the average
    Use seed to create the same underlying light curve
    that is sampled through a different poisson noise.
    This allows zero lag light curves with perfect coherence
    to be simulated
    """
    # 2.2018e+00 2.2170e-01 2.4064e-01 -2.0275e-01 "PhoIndex "
    # 3.9876e-09 4.3754e-09 5.2287e-09 -3.5220e-09 "norm "
    psdpar = [1.36e-10, -2.548]
    nx  = len(x)
    mux = np.nanmean(x)
    inan = np.isnan(x)

    sim = az.SimLC(seed)
    # use the break from Papadakis+95
    sim.add_model('broken_powerlaw', [psdpar[0], -1, psdpar[1], 1e-7])
    ns  = max([2, nsim])
    y   = []
    for ii in range(ns//2+1):
        sim.simulate(nx*2, dt, mux, 'rms')
        y += np.split(sim.x, 2)
    y   = sim.add_noise(np.array(y), seed=None, dt=dt)[:nsim]
    y[:, inan] = np.nan
    return y


def combine_en_segments(rate_all, rerr_all, ibin, iref=None):
    """Combined segments from different energies to produce 1 (iref=None) or 2 bins.
    We assume the arrays from different energies have the same shape
    
    rate_all, rerr_all: the output of lc_to_segments with dim: (nen, nseg_total)
    ibin: bin of interest
    iref: reference band. -1 for all (excluding ibin of course)
    
    return: rate, rerr, Rate, Rerr
    """

    nen = len(rate_all)

    # make sure we are dealing with lists #
    if not isinstance(ibin, list): ibin = [ibin]

    # get the rate at the bins of interest #    
    rate = np.sum(np.array(rate_all)[ibin], 0)
    rerr = np.sum(np.square(rerr_all)[ibin], 0)**0.5
    
    # and the reference if needed #
    Rate, Rerr = None, None
    if not iref is None:
        if not isinstance(iref, list):
            iref = list(range(nen)) if iref==-1 else [iref]
        iref = [i for i in iref if not i in ibin]

        Rate = np.sum(np.array(rate_all)[iref], 0)
        Rerr = np.sum(np.square(rerr_all)[iref], 0)**0.5

    return rate, rerr, Rate, Rerr


def calc_lag_en(rate_all, rerr_all, dt, fqbin, indv=None, iref=-1, overlap=None, **kwargs):
    """Calculate lag vs energy for the total and individual observations
    
    rate_all, rerr_all: the output of lc_to_segments. dim: (nen, nseg)
    dt: time bin
    fqbin: binning dict
    indv: a list of obs indices to group. e.g. [[0,1,2], [3,4,5]] etc.
    iref: reference band. -1 for all (excluding ibin of course)
    overlap: number of bins to overlap e.g. 1, 2, or 3. default None
    **kwargs: any parameters to pass to az.LCurve.calculate_lag
    
    Return: lag (nen, 3, nfq), ilag (nen, nindv, 3, nfq)
    
    """

    # these are arrays of objects because the segments may have different lengths #
    rate_all = np.array(rate_all)
    rerr_all = np.array(rerr_all)
    
    nen, nseg = rate_all.shape

    ibins = list(range(nen))
    if not overlap is None:
        ibins = [ibins[i:i+overlap] for i in range(nen-overlap+1)]
    
    Lag, iLag, LagE, iLagE = [], [], [], []
    for ibin in ibins:
        rate, rerr, Rate, Rerr = combine_en_segments(rate_all, rerr_all, ibin, iref)
        lag  = az.LCurve.calculate_lag(rate, Rate, dt, fqbin, rerr=rerr, Rerr=Rerr, **kwargs)
        Lag.append(lag[:3])
        LagE.append(lag[3])

        # indiviudal obs #
        ilag, ilagE = [], []
        if indv is None: continue
        for ii in indv:
            rate, rerr, Rate, Rerr = combine_en_segments(rate_all[:,ii], rerr_all[:,ii], ibin, iref)
            lag  = az.LCurve.calculate_lag(rate, Rate, dt, fqbin, rerr=rerr, Rerr=Rerr, **kwargs)
            ilag.append(lag[:3])
            ilagE.append(lag[3])
        iLag.append(ilag)
        iLagE.append(ilagE)

    # nen, 3, nfq
    lag    = np.array(Lag)
    # nen, nindv, 3, nfq
    ilag   = np.array(iLag)
    return lag, ilag, LagE, iLagE

def lag_en_null_test(en, l, le, verbosity=True):
    """Test lag-enery against a null hypothesis of a const and log-linear model"""
    import scipy.optimize as opt
    import scipy.stats as st

    # const #
    def fit_1(x, *p): return x*0+p[0]
    f_1  = opt.curve_fit(fit_1, en, l, [l.mean()], sigma=le)
    c2_1 = np.sum(((l - fit_1(en, *f_1[0]))/le)**2)
    p_1  = 1 - st.chi2.cdf(c2_1, df=len(en)-1)
    nsigma = st.norm.ppf(p_1/2)
    text = '\n- fit 1: {} {:6.3} {:6.3} {:6.3}'.format(f_1[0], c2_1, p_1, nsigma)

    # log-linear
    def fit_2(x, *p): return p[0] + p[1] * np.log10(x)
    f_2  = opt.curve_fit(fit_2, en, l, [l.mean(), 0.1], sigma=le)
    c2_2 = np.sum(((l - fit_2(en, *f_2[0]))/le)**2)
    p_2  = 1 - st.chi2.cdf(c2_2, df=len(en)-2)
    nsigma = st.norm.ppf(p_2/2)
    text += '\n- fit 2: {} {:6.3} {:6.3} {:6.3}'.format(f_2[0], c2_2, p_2, nsigma)
    if verbosity:
        print(text)
    return [f_1, c2_1, p_1], [f_2, c2_2, p_2], text


def lag_en_pipeline(lcdata, fqbin=-1, indv=None, iref=-1, nsim=50, overlap=None):
    """read lc, calc lag, run simulations and plot it"""

    # some common input #
    if fqbin == -1:
        fqbin  = {'bins': [3e-5, 5e-4]}
    if indv == -1:
        indv   = [[0,1,2], [3,4,5], [6,7], [8,9,10,11,12,13,14,15,16], 
                  [8,9,10,11] ,[12,13,14,15,16],
                  [0,1,2,6,7,8,9,10,11], [3,4,5,12,13,14,15,16]]
    kwargs = dict(taper=True)
    

    # read light curve files: Lc has dims: (nen, nobs) #
    print('prepare data ...')
    Lc, en, ene = lcdata 
    nen, nobs   = len(Lc), len(Lc[0])
    dt = Lc[0][0].dt
    #nrebin = 256//np.int(dt)
    #Lc = [[l.rebin(nrebin) for l in ll] for ll in Lc]; dt *= nrebin
    # update en overlap is requested #
    if not overlap is None:
        ibins = list(range(nen))
        ibins = [ibins[i:i+overlap] for i in range(nen-overlap+1)]
        ene = np.array([( (en[i[-1]]+ene[i[-1]])-(en[i[0]]-ene[i[0]]) )/2 for i in ibins])
        en  = np.array([(en[i[0]]+en[i[-1]])/2 for i in ibins])
    
    
    min_length = np.int(2e3/dt)    
    rate_all, rerr_all, time_all, seg_idx = lc_to_segments(Lc, min_seg_length=min_length)
    
    # map indv which is for the obs number to indv_seg which is for the segment number #
    indv_seg = None
    if not indv is None:
        indv_seg = [[j for j,jj in enumerate(seg_idx) if jj in i] for i in indv]

    # calculate lag #
    print('calculate lag ...')
    lag,ilag, lagE, ilagE = calc_lag_en(rate_all, rerr_all, dt, fqbin, indv_seg, iref, overlap, **kwargs)
    
    # simulations #
    print('run simulations ...')
    rate_sim = []
    # use same seed for all energies; seperate for every segment, to simulate 0-lag
    seeds = np.random.randint(10000, size=len(seg_idx))
    for R in rate_all:
        rate_sim.append([simulate_like(r, dt, nsim, s) for r,s in zip(R, seeds)])
    # make rate_sim dims: (nsim, nseg, nen, ..)
    rate_sim = [[np.array([rate_sim[ie][iseg][isim] for iseg in range(len(seg_idx))])
             for ie in range(nen)] for isim in range(nsim)]
    rerr_sim = [[(r2/dt)**0.5 for r2 in r1] for r1 in rate_sim]
    
    LagS = [calc_lag_en(rs, rse, dt, fqbin, indv_seg, iref, overlap, **kwargs)[:2]
               for rs,rse in zip(rate_sim, rerr_sim)]
    # unzip so we have dims: (nsim, nen, 3, nfq) and (nsim, nen, nindv, 3, nfq)
    lagS, ilagS = zip(*LagS)
    lagS, ilagS = np.array(lagS), np.array(ilagS)
    

    # null tests for lag-vs-energy #
    print('calculating null tests ...')
    nfq = lag.shape[-1]
    fit_data  = [lag_en_null_test(en, lag[:,1,ifq], lag[:,2,ifq], 0) for ifq in range(nfq)]
    fit_sim  = [lag_en_null_test(en, lag[:,1,ifq], lagS[:,:,1,ifq].std(0), 0)
                        for ifq in range(nfq)]
    if not indv is None:
        ifit_data = [[lag_en_null_test(en, ilag[:,ii,1,ifq], ilag[:,ii,2,ifq], 0)
                        for ii in range(len(indv))] for ifq in range(nfq)]
        ifit_sim = [[lag_en_null_test(en, ilag[:,ii,1,ifq], ilagS[:,:,ii,1,ifq].std(0), 0)
                        for ii in range(len(indv))] for ifq in range(nfq)]

    text  = ''
    for ifq in range(nfq):
        text += fit_data[ifq][2].replace('\n','\n#') + fit_sim[ifq][2].replace('\n','\n#') + '\n'
        if not indv is None:
            text += ('\n'.join(['#- {} -# {}'.format(i+1, idv) + 
                   ifit_data[ifq][i][2].replace('\n','\n#') +
                   ifit_sim[ifq][i][2].replace('\n','\n#') for i,idv in enumerate(indv)]))
        text += '\n'
    
    extra = [text, lagE, ilagE]
    return en, ene, lag, lagS, ilag, ilagS, extra


def plot_ilag(en, lag, ilag, lagS=None, ilagS=None, figsize=None):
    """plot the results of lag_en_pipeline, with simulations if given"""
    
    doindv = False
    try:
        _,nindv,_,nfq = ilag.shape
        doindv = True 
    except:
        nfq, nindv = lag.shape[-1], 1
    if figsize is None: figsize=(14,3*nfq)
    fig = plt.figure(figsize=figsize)
    for ifq in range(nfq):
        if not lagS is None:
            lm, ls = lagS[:,:,1,ifq].mean(0), lagS[:,:,1,ifq].std(0)
        for il in range(nindv):
            if nindv < 10:
                ax = plt.subplot(nfq, nindv, il+ifq*nfq+1)
            else:
                ax = plt.subplot(nfq*2, nindv//2+1, il+ifq*nfq+1)
            ax.set_xscale('log')
            ax.set_ylim([-2,2])
            ax.set_xticks([3,6])
            ax.xaxis.set_major_formatter(plt.ScalarFormatter())
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.errorbar(en, lag[:,1,ifq]/1e3, lag[:,2,ifq]/1e3, fmt='o-', alpha=0.5)
            if doindv:
                ax.errorbar(en, ilag[:,il,1,ifq]/1e3, ilag[:,il,2,ifq]/1e3, fmt='s-')

            # simulations #
            if not lagS is None:
                ax.fill_between(en, (lm-ls)/1e3, (lm+ls)/1e3, alpha=0.3)
                if doindv:
                    ilm, ils = ilagS[:,:,il,1,ifq].mean(0), ilagS[:,:,il,1,ifq].std(0)
                    ax.fill_between(en, (ilm-ils)/1e3, (ilm+ils)/1e3, alpha=0.3)
    plt.tight_layout(pad=0)

    
def write_ilag(en, ene, lag, ilag, lagS=None, ilagS=None, suff=''):
    """Write the lag results to a string to be written to a veusz file"""

    doindv = False
    try:
        nen,nindv,_,nfq = ilag.shape
        doindv = True
    except:
        nfq, nindv,nen = lag.shape[-1], 1, lag.shape[0]
        
    if suff != '' and suff[0] != '_': suff = '_%s'%suff

    text = ''
    for ifq in range(nfq):
        
        text += '\ndescriptor en_f{0}{1},+- lag_f{0}{1},+-\n'.format(ifq+1, suff)
        text += '\n'.join(['{} {} {} {}'.format(en[ie], ene[ie], lag[ie,1,ifq]/1e3,
                lag[ie,2,ifq]/1e3) for ie in range(nen)])
        if not lagS is None:
            lm, ls = lagS[:,:,1,ifq].mean(0), lagS[:,:,1,ifq].std(0)
            text += '\ndescriptor slag_f{0}{1},+- sLag_f{0}{1},+-\n'.format(ifq+1, suff)
            text += '\n'.join(['{0} {1} {2} {1}'.format(lm[ie]/1e3, ls[ie]/1e3, lag[ie,1,ifq]/1e3)
                          for ie in range(nen)])

        if not doindv: continue
        for il in range(nindv):
            text += '\ndescriptor lag_f{}_i{}{},+-\n'.format(ifq+1, il+1, suff)
            text += '\n'.join(['{} {}'.format(ilag[ie,il,1,ifq]/1e3, ilag[ie,il,2,ifq]/1e3)
                              for ie in range(nen)])
            # simulations #
            if not lagS is None:
                ilm, ils = ilagS[:,:,il,1,ifq].mean(0), ilagS[:,:,il,1,ifq].std(0)
                text += '\ndescriptor slag_f{0}_i{1}{2},+- sLag_f{0}_i{1}{2},+-\n'.format(
                            ifq+1, il+1, suff)
                text += '\n'.join(['{0} {1} {2} {1}'.format(ilm[ie]/1e3, ils[ie]/1e3,
                            ilag[ie,il,1,ifq]/1e3) for ie in range(nen)])
    return text
