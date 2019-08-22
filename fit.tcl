source ~/codes/xspec/az.tcl
load /home/abzoghbi/soft/xspec_local/relxill/librelxill.so
lmod relxill

# ---------------------------------------------- #
# fit po+zga to 6-7 keV to get enegy of the line #
proc fit_1 {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-6., 7.-**
    mo po+zga & 1 & 1e-2& 6.4 & 0.02 & 3.319e-3& 1e-4
    fit 1000
    az_calc_errors [az_free_params] fits/fit_1$suff 1.0
}


# ----------------------------------------------------------- #
# fit po+zga to 0.3-1 keV to get enegy of the strongest lines #
proc fit_1a {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-.3, 1.-**
    mo po+zga+zga & 2.5 & 1e-2& .55 .1 .5 .5 .6 .6 & 0 -1 & 3.319e-3& 1e-4 &\
        .9 .1 .85 .85 1 1 & 0 -1 & 3.319e-3& 1e-4
    fit 1000
    az_calc_errors [az_free_params] fits/fit_1a$suff 1.0
}

# ----------------------------- #
# Fit rgs with sum of gaussians #
proc fit_gain_1 {sfile {suff ""}} {
    para st 5; para lev 1; para err 5; query yes
    statis chi
    da $sfile
    ign 0.0-.35,1.-**
    mo po & 2 & 1e-3; fit 100
    set p [tcloutr null]
    
    while { $p < 0.05 } {
        scan [tcloutr peakr] "%f %f %f %f" e1 n1 e2 n2
        if { $n1 > [expr -$n2] } {
            set e $e1; set n "$n1 .1 1e-8 1e-8"; set s "1e-3 .1 1e-4 1e-4 0.01 0.01"
        } else {
            set e $e2; set n "$n2 .1 -1 -1 -1e-8 -1e-8"; set s "0 -1"
        }
        set e [expr $e/(1-3.319e-3)]
        
        add 1 zga & $e 0.1 0.35 .35 1 1& $s & 3.319e-3 & $n; fit
        stepp best 1 0.35 1. 1000
        err max 200 1. 1 4
        
        set p [tcloutr null]
        freez 1 2 4
    }
    rm fits/fit_gain_1$suff.xcm > /dev/null 2>&1
    save all fits/fit_gain_1$suff
}

# ----------------------------------------- #
# apply rgs fit in gain_1 to pn to get gain #
proc fit_gain_2 {ispec} {
    para lev 1; para st 40; query yes
    set suff __$ispec
    rm fits/fit_gain_2$suff*
    @fits/fit_gain_1$suff
    da spec_$ispec.grp
    statis chi
    ign 0.0-.35, 1.-6.,7.-**
    
    set ng [expr ([tcloutr modpar]-2)/4]
    set txt "zga"; for {set i 1} { $i<$ng } {incr i} {append txt "+zga"}
    eval "editm const*($txt)+po & 1"
    add 1 po & 0 & 1e-2
    # cut the hard powerlaw so it doesn't affect the soft stuff
    add 1 tbabs & 10 -1 
    #add 1 zga & 6.4 -1 & 0.01 & 3.319e-3& 1e-4
    add 2 xillver & 1.6 -1 & 1 -1 & 1000 -1& 0 -1 & 3.319e-3& 30 -1 & -1 -1 & 1e-3
    fit 1000
    set c0 [tcloutr stat]; set d0 [lindex [tcloutr dof] 0]
    save all fits/fit_gain_2$suff
    gain fit 1 & 1 & 0; fit
    rm tmp.xcm; save all tmp.xcm
    @tmp.xcm; fit
    stepp best r1 0.98 1.02 50 r2 -0.02 0.02 50; fit
    set c1 [tcloutr stat]; set d1 [lindex [tcloutr dof] 0]
    save all fits/fit_gain_2g$suff
    ftest $c1 $d1 $c0 $d0
    set fprob [tcloutr ftest]
    set slope [lindex [tcloutr gain 1 slope] 0]
    set offset [lindex [tcloutr gain 1 offset] 0]
    set fp [open fits/fit_gain_2g$suff.log w]
    puts $fp "$c0 $d0 $c1 $d1 $fprob $slope $offset"
    close $fp
}



############################################
############################################
############## EXPLORE SPECTRA #############
############################################
############################################

# ---------------------------------- #
# K edge energy, using a simple edge #
proc fit_2 {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-6., 10.-**
    
    
    mo zedge*cflux*po+cflux*zga+cflux*zga & 7 .1 7 7 7.5 7.5 & 0.1 & 3.319e-3& \
        7 & 10 & -10& 1.6 & 1 -1 & \
        2 & 10 & -11& 6.4 -1 & 0.01 & 3.319e-3 & 1 -1 & \
        2 & 10 & =11-0.920& 7.06 -1 & =13 & 3.319e-3& 1 -1
    fit 1000
    az_calc_errors [az_free_params] fits/fit_2$suff 1.0
}

# -------------------------------- #
# neutral vs ionized edge at 7 keV #
proc fit_2ab {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    da $sfile
    ign 0.0-6., 10.-**
    
    # neutral absorber #
    mo ztbabs*cflux*po + cflux*zga+cflux*zga & 5 & 3.319e-3 & \
        7 & 10 & -10 & 1.6 & 1 -1 & \
        2 & 10 & -11& 6.4 -1 & 0.01 & 3.319e-3 & 1 -1 & \
        2 & 10 & =10-0.920& 7.06 -1 & =12 & 3.319e-3& 1 -1
    fit 1000
    az_calc_errors [az_free_params] fits/fit_2a$suff 1.0
    
    # ionized absorber #
    set nh [lindex [tcloutr par 1] 0]
    add 1 zxipc & $nh & 0 & 1 -1& 3.319e-3
    del 2; fit
    stepp best 2 -2 3 200
    az_calc_errors [az_free_params] fits/fit_2b$suff 1.0
}


# ---------------------------- #
# redo fit_2a fitting 4-10 keV #
proc fit_2c {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2a$suff
    no **; ign 0.0-4.,10.-**; fit 1000
    az_calc_errors [az_free_params] fits/fit_2c$suff 1.0
}


# ----------------------------- #
# add hi-xi to fit_2c; 4-10 keV #
proc fit_2d {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2c$suff
    add 1 zxipc & 2 & 3.7 .1 3.2 3.2 & 1 -1 & 3.319e-3; fit 1000
    az_calc_errors [az_free_params] fits/fit_2d$suff 1.0
}


# ----------------------------- #
# add lo-xi to fit_2d; 4-10 keV #
proc fit_2e {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2d$suff
    add 1 zxipc & 2 & 1.5 & 1 -1 & 3.319e-3; fit 1000
    stepp best 2 0 4 200; fit
    az_calc_errors [az_free_params] fits/fit_2e$suff 1.0
}

# -------------------------------- #
# add broad ga to fit_2d; 4-10 keV #
proc fit_2f {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2d$suff
    add 3 zga & 6. .1 4 4 7.5 7.5 & 0.4 .1 .1 .1 1.2 1.2 & 3.319e-3 & 1e-5; fit 1000
    az_calc_errors [az_free_params] fits/fit_2f$suff 1.0
}



# ------------------------------------------- #
# extend fit_2d (neutral+hi-xi WA) to 2.5 keV #
proc fit_2g {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2d$suff
    no **; ign 0.0-2.5,10.-**;fit 1000
    az_calc_errors [az_free_params] fits/fit_2g$suff 1.0
}


# ------------------------------------ #
# Add an extra unabsorbed po to fit_2g #
proc fit_2h {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2g$suff
    editm tbabs(zxipcf*zTBabs(cflux*po) + cflux*zga + cflux*zga) & 0.0198 -1
    editm tbabs(zxipcf*zTBabs(cflux*po) + po + cflux*zga + cflux*zga) & =11& 1 -1
    editm tbabs(zxipcf*zTBabs(cflux*po) + cflux*po + cflux*zga + cflux*zga) & 7 & 10 & -11
    fit 1000
    az_calc_errors [az_free_params] fits/fit_2h$suff 1.0
}


# --------------------------- #
# test WA instead of extra PL #
proc fit_2h1 {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2h$suff
    del 6; del 6
    add 2 zxipc & 2 .1 0 0 50 50 & 1.5 & 1 -1 & 3.319e-3; fit 100
    stepp best 3 -1 4 200; fit
    az_calc_errors [az_free_params] fits/fit_2h1$suff 1.0
}


# ---------------------------------------- #
# use xillver for lines, keep the extra PL #
proc fit_2i {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2h$suff
    add 8 xillver & =11 & 1 -1& 1000 -1& 0 -1& 3.319e-3 & 30 -1& -1 -1 & 1 -1
    add 8 cflux & 6 & 6.8 & -11
    del 10; del 10; del 10; del 10
    fit 1000
    az_calc_errors [az_free_params] fits/fit_2i$suff 1.0
}

# ----------------------------------------- #
# test WA instead of extra PL; xillver case #
proc fit_2i1 {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2i$suff
    del 6; del 6
    add 2 zxipc & 2 & 1.5 & 1 -1 & 3.319e-3; fit 100
    new 10 ,,.01 .01
    stepp best 3 -1 4 200; fit
    stepp best 7 3 4.5 80; fit
    az_calc_errors [az_free_params] fits/fit_2i1$suff 1.0
}

# --------------------------------------------- #
# rm extra PL from fit_2i, make xillver ionized #
proc fit_2j {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_2i$suff
    del 6; del 6; new 19 0 .1; fit 1000
    az_calc_errors [az_free_params] fits/fit_2j$suff 1.0
}



# ------------------------------------ #
# fit descriptive models to 0.3-10 keV #
proc fit_3 {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile
    ign 0.0-0.3, 10.-**
    
    mo tbabs*(ztbabs*cflux*po+cflux*zga + cflux*po+cflux*zga+cflux*zga + cflux*ga) & \
        0.0198 -1 & 10 & 3.319e-3 & \
        2 & 10 & -10 & 1.7 .1 1.4 1.4 1.9 1.9& 1 -1& \
        2 & 10 & -11 & 6.4 -1& 0.02 0.1 0 0 0.1 0.1 & 3.319e-3 & 1 -1 & \
        0.3 & 2 & -10 & 2.7& 1 -1 & \
        0.3 & 2 & -11& 0.56 .1 .5 .5 .6 .6 & 0 -1& 3.319e-3 & 1 -1 & \
        0.3 & 3 & -11 & 0.9 .1 .8 .8 1. 1. & 0 -1 & 3.319e-3& 1 -1 & \
        1 & 10 & -10 & 2.7 -.1 1 1 4 4& 0.8 -.1 & 1 -1
    fit 1000
    thaw 38 39
    az_calc_errors [az_free_params] fits/fit_3$suff 1.0
}

# ---------------------------------------- #
# Use PC in the descriptive model of fit_3 #
proc fit_3a {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3$suff
    add 3 tbpcf & =2 & 1. .1 & 3.319e-3
    untie 4; new 2 2
    del 14; del 14; fit 1000
    az_calc_errors [az_free_params] fits/fit_3a$suff 1.0
}

# --------------------------------------- #
# Use PC similar 3a, but constrain cf>0.9 #
proc fit_3a1 {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3a$suff
    new 5 .95 .1 .9 .9 1 1; fit 1000
    az_calc_errors [az_free_params] fits/fit_3a1$suff 1.0
}




# ------------------------------------ #
# Make full absorber in fit_3a ionized #
proc fit_3b {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3a$suff
    add 3 zxipc & =2 & 1 .1 -1 -1 2 2 & 1 -1 & 3.319e-3
    del 2; 
    new 6,,,,,40 40
    fit 1000
    stepp best 3 -1 2 40
    az_calc_errors [az_free_params] fits/fit_3b$suff 1.0
}


# -------------------------------------- #
# fit_3b but fix full nh at fit_3a value #
proc fit_3c {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_3a$suff
    add 3 zxipc & =2 & 1 .1 -1 -1 2 2 & 1 -1 & 3.319e-3
    del 2; 
    freez 6 7;fit 1000
    stepp best 3 -1 2 40
    az_calc_errors [az_free_params] fits/fit_3c$suff 1.0
}


# --------------------------------------------- #
# leaked xillver; similar to 2i1 but 0.3-10 keV #
proc fit_3d {sfile {suff ""}} {
    para err 5; para lev 5
    da $sfile; ign 0.0-0.3,10.-**
    mo TBabs(zxipcf*ztbabs*(cflux*powerlaw) + cflux*xillver +cflux*zga + cflux*zga) &\
        0.0198 -1 & 3 & 1 .1 -2 -2 2.5 2.5 & 1 -1 & 3.319e-3 & 10 & 3.319e-3 & \
        2 & 10 & -10 & 1.7 .1 1.4 1.4 1.9 1.9& 1 -1& \
        6 & 6.8 & -11 & =11 & 1 -1 & 1000 & 0 -1 & 3.319e-3& 30 -1 & -1 -1 & 1 -1 &\
        0.3 & 2 & -11& 0.56 .1 .5 .5 .6 .6 & 0 -1& 3.319e-3 & 1 -1 & \
        0.3 & 3 & -11 & 0.9 .1 .8 .8 1. 1. & 0 -1 & 3.319e-3& 1 -1    
    fit 1000
    thaw 19; fit
    az_calc_errors [az_free_params] fits/fit_3d$suff 1.0
}




# -------------------------------------------- #
# difference spectra to check soft variability #
proc soft_variability_diff {ihigh ilow} {
    # ihigh,ilow are the spectra ids of the high/low fluxes to consider
    
    setpl ener; setpl area
    setpl del all
    
    # remove tmp files if present 
    rm tmp_?.dat >/dev/null 2>&1
    
    # low #
    da spec_$ilow.grp
    ign 0.0-0.3,10.-**
    setpl co wd tmp_l.dat; pl d; setpl del 1
    
    # high #
    da spec_$ihigh.grp
    ign 0.0-0.3,10.-**
    setpl co wd tmp_h.dat; pl d; setpl del 1
    
    # difference #
    back spec_$ilow.grp
    
    
    # simple WA model #
    mo TBabs(TBpcf*zxipcf*powerlaw + powerlaw) & 0.0198 -1 & \
        2 .1 0 0 50 50 & .98 .1 .7 .7 & 3.319e-3 & \
        3 & .5 & 1 -1 & 3.319e-3 & 1.7 -1 & 1e-2 & 2. -.1 & 5e-5
    fit 1000; thaw 11; fit
    rm fits/fit_soft_var_diff_${ihigh}_${ilow}.* >/dev/null 2>&1
    save all fits/fit_soft_var_diff_${ihigh}_${ilow}
    
    setpl add
    setpl co wd tmp_d.dat; pl d; setpl del 1
}

# ---------------------------------------- #
# illustrate the need for partial covering #
proc fit_hard_pc {{sfile spec_24.grp}} {
    
    load /u/Reverb1/abzoghbi/soft/xspec_local/tbfeopc/libtbfeopc.so
    lmod tbfeopc
    setpl ener; setpl area
    setpl del all
    setpl add
    rm tmp_?.dat >/dev/null 2>&1
    
    da $sfile
    ign 0.0-2.,10.-**
    
    # full covering
    mo tbpc*po + zga & 12 & 1 -1 & 3.319e-3& 1.3& 1e-2 & 6.4 -1 & 0.02 -1& 3.319e-3& 1e-4
    fit 1000
    setpl co wd tmp_1.dat; pl d; setpl del 1
    
    # parital covering
    thaw 2; fit
    setpl co wd tmp_2.dat; pl d; setpl del 1
    
    # variable Fe abundance
    add 1 tbfeopc & 12 & 1 -1 & 1 .1& 1 -1 & 3.319e-3
    del 2; fit
    setpl co wd tmp_3.dat; pl d; setpl del 1
    
    # warm absorber #
    add 1 zxipcf & 10 & 2 & 1 -1 & 3.319e-3
    del 2; fit
    para st 20; query yes
    stepp best 2 -2 2 200
    setpl co wd tmp_4.dat; pl d; setpl del 1
    
    # extra gaussian #
    add 1 tbpcf & 12 & 1 -1& 3.319e-3
    del 2; fit
    add 4 zga & 3 & 1 & 3.319e-3& 1e-4; fit
    setpl co wd tmp_5.dat; pl d; setpl del 1
    
    # extra powerlaw #
    add 4 po & =4& 1e-3
    del 5; fit
    setpl co wd tmp_6.dat; pl d; setpl del 1
}

# ----------------------------------------------------------- #
# similar to fit_3a; replace lines by xillver, add soft lines #
proc fit_4 {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    da $sfile; ign 0.0-.3,10.-**
    mo tbabs*(zxipc(zxipc*tbpcf*cflux*po) + cflux*xillver + brems+zga+zga) & 0.0198 -1 & \
        2 -.1 & 1 -.1 -1 -1 2 2 & 1 -1 & 3.319e-3 & \
        2 & 3.5 -.1 3.2 3.2 5 5 & 1 -1 & 3.319e-3 & \
        3 .1 0.01 0.01 & .95 .1 .8 .8& 3.319e-3 & \
        7 & 10 & -10 & 1.6 -.1 & 1 -1 & \
        6 & 6.8 & -11 & =16& 1 -1 & 1000 -1 & 0 -1 & 3.319e-3& 30 -1 & -1 -1 & 1 -1 & \
        0.5 -.1 & 1e-3 & .56 & 0 -1 & 3.319e-3& 1e-4 & 0.9 & 0 -1 & 3.319e-3 & 1e-4
    fit 1000
    thaw 2 3 7 16 29; fit
    stepp best 3 -1 1.7 200
    az_calc_errors [az_free_params] fits/fit_4$suff 1.0
}


# ---------------------- #
# add a single soft line #
proc _add_soft_line {iadd ien {abs 0} {extra_cmd ""}} {
    para lev 5; para st 30; fit 100

    set c0 [tcloutr stat]; set d0 [lindex [tcloutr dof] 0]
    query no
    if {$abs == 0} {
        add $iadd zga & 0.3 .1 .3 .3 4. 4. & 0 -1 & 3.319e-3& 1e-5 .1 1e-8 1e-8 & /*; 
    } else {
        add $iadd zga & 2 .1 .2 .2 5.5 5.5 & 0 -1 & 3.319e-3& -1e-7 .1 -1 -1 -1e-8 -1e-8 & /*; 
    }
    fit 1
    query yes
    if {$abs == 0} {
        stepp best $ien 0.3 4. 200
    } else {
        stepp best $ien 2 5.5 200
    }
    
    eval $extra_cmd
    
    # check improvement #
    set c1 [tcloutr stat]; set d1 [lindex [tcloutr dof] 0]
    set fprob 1
    if { $c1 < $c0 } {
        ftest $c1 $d1 $c0 $d0
        set fprob [tcloutr ftest]
    }
    return $fprob
}


# ------------------------------ #
# base function for running mcmc #
proc _run_mcmc {fit suff} {
    # fit: e.g. fit_2; suff: e.g. __1
    para walk 5
    @fits/${fit}${suff}
    fit 1000
    chain len 10000
    chain burn 100000
    chain walk 200
    rm fits/${fit}_mc${suff}.fits >/dev/null 2>&1
    chain run fits/${fit}_mc${suff}.fits
}


# -------------------------------- #
# add soft emission lines to fit_4 #
proc fit_4a {sfile {suff ""}} {
    # stop once 3 lines are added and don't improve the fit
    para err 5; para lev 5; query yes
    @fits/fit_4$suff
    set ifreez {6 7} 
    freez $ifreez
    set ic 0; rm tmp${suff}_${ic}.xcm >/dev/null 2>&1; save all tmp${suff}_${ic}
    set stop 0
    for {set ic 1} {$ic<=20} {incr ic} {
        set fprob [_add_soft_line 10 31 0 "stepp best 3 -1 1.7 100"]
        rm tmp${suff}_${ic}.xcm >/dev/null 2>&1; save all tmp${suff}_${ic}
        if { $fprob > 0.1 } { 
            if {$stop == 1} { break }
            incr stop
        } else {
            set stop 0
        }
    }
    @tmp${suff}_[expr $ic-2].xcm; fit
    thaw $ifreez; fit
    stepp best 3 -1 1.7 100
    stepp best 7 3.2 5 100
    az_calc_errors [lrange [az_free_params] 0 10] fits/fit_4a$suff 1.0
    rm tmp${suff}_*.xcm >/dev/null 2>&1
}


# ------------------------------------------------ #
# add absorption lines between 2-5.5 keV to fit_4c #
proc fit_4b {sfile {suff ""}} {
    # stop once two lines are added and don't improve the fit
    para err 5; para lev 5
    @fits/fit_4a$suff
    set ic 0; rm tmp${suff}_${ic}.xcm >/dev/null 2>&1; save all tmp${suff}_${ic}
    set stop 0
    for {set ic 1} {$ic<=4} {incr ic} {
        set fprob [_add_soft_line 10 31 1 "stepp best 3 -1 1.7 100"]
        rm tmp${suff}_${ic}.xcm >/dev/null 2>&1; save all tmp${suff}_${ic}
        if { $fprob > 0.1 } { 
            if {$stop == 1} { break }
            incr stop
        } else {
            set stop 0
        }
    }
    @tmp${suff}_[expr $ic-2].xcm; fit
    az_calc_errors [lrange [az_free_params] 0 10] fits/fit_4b$suff 1.0
    rm tmp${suff}_*.xcm >/dev/null 2>&1
}


# -------------------------------- #
# add mission and absorption lines #
proc fit_4c {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_4b$suff
    #freez 31-**
    add 10 zga & 8.04 -1 & 0 -1 & 3.319e-3& 1e-5 .01 1e-8 1e-8;fit 1000
    new 31 ,, .1 7.9 7.9 8.5 8.5; fit
    add 10 zga & 8.7 -1 & 0 -1 & 3.319e-3& -1e-5 .01 -1 -1 -1e-8 -1e-8;fit 1000
    new 31 ,, .1 8.5 8.5 9.5 9.5; fit
    az_calc_errors [lrange [az_free_params] 0 10] fits/fit_4c$suff 1.0
}

# -------------------------------- #
# freeze emission/absorption lines #
proc fit_4d {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_4c$suff
    freez 31-**
    az_calc_errors [az_free_params] fits/fit_4d$suff 1.0
}
proc fit_4d_mc {sfile suff} { _run_mcmc fit_4d $suff }


# --------------------------------------- #
# similar to 4d, replace pc with extra po #
proc fit_4d1 {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_4d$suff
    add 3 po & & 1 -1
    add 3 cflux & 7 & 10 & -11
    new 16 1 -1
    new 9=21
    set ifree [az_free]
    freez $ifree
    thaw 20 8; fit
    thaw $ifree
    az_calc_errors [az_free_params] fits/fit_4d1$suff 1.0
}

# ---------------------------------------- #
# fit_4d simplified, fitted to nustar data #
proc fit_4dnu {sfile {suff ""}} {
    para err 5; para lev 5; para st 10; query yes
    da $sfile
    ign 0.0-3.,79.-**
    mo TBabs(zxipcf(zxipcf*TBpcf*cutoffpl) + xillver) & 0.0198 -1 & \
        6 & 1 & 1 -1 & 3.319e-3 & 20 & 4 & 1 -1 & 3.319e-3 & 30 & 0.5 .1 .1 .1 & 3.319e-3 & \
        1.8 & 300 .1 30 30 1000 1000 &1e-2 & \
        =13 & 1 -1 & =14 & 0 -1 & 3.319e-3 & 30 -1 & -1 -1 & 1e-4
    renorm; fit 1000
    stepp best 3 -1 1.5 30
    stepp best 7 3.2 4.5 30
    stepp best 14 50 400 50
    puts "az_calc_errors [az_free_params] fits/fit_4dnu$suff 1.0"
}



# --------------------------------------------- #
# similar to 4d, allow more gaussian to be free #
proc fit_4e {sfile {suff ""}} {
    para err 5; para lev 5
    @fits/fit_4d$suff
    set np [tcloutr modpar]
    thaw [expr $np-3] [expr $np]
    thaw [expr $np-7] [expr $np-4]
    az_calc_errors [az_free_params] fits/fit_4e$suff 1.0
}



################ SUBSPEC ###############

proc _identify_group {suff} {
    # identify the group (1 of 22) to which the suff belongs
    set ispec [string trim $suff _]
    if {$ispec <= 4} { set igrp 1
    } elseif {$ispec <=14} { set igrp 2
    } elseif {$ispec <= 17} {set igrp 3 
    } elseif {$ispec <= 19} {set igrp 4
    } elseif {$ispec <= 21} {set igrp 5
    } elseif {$ispec <= 23} {set igrp 6
    } elseif {$ispec <= 28} {set igrp 7
    } elseif {$ispec <= 32} {set igrp 8
    } elseif {$ispec <= 33} {set igrp 10
    } elseif {$ispec <= 34} {set igrp 11
    } elseif {$ispec <= 35} {set igrp 12
    } elseif {$ispec <= 36} {set igrp 13
    } elseif {$ispec <= 37} {set igrp 15
    } elseif {$ispec <= 41} {set igrp 16
    } elseif {$ispec <= 45} {set igrp 17
    } elseif {$ispec <= 51} {set igrp 18
    } elseif {$ispec <= 55} {set igrp 19
    } elseif {$ispec <= 60} {set igrp 20
    } elseif {$ispec <= 65} {set igrp 21
    } elseif {$ispec <= 70} {set igrp 22
    } elseif {$ispec <= 75} {set igrp 23
    } elseif {$ispec <= 81} {set igrp 24}
    return $igrp
}

# ----------------------------- #
# fit_4d
proc fit_sub_4a {sfile suff} {
    para err 5; para lev 5
    set igrp [_identify_group $suff]
    @fits/fit_4a__$igrp
    da $sfile; ign 0.0-0.3,10.-**
    fit
    freez 6 7 31-**  
    az_calc_errors [az_free_params] fits/fit_sub_4a$suff 1.0
}

# ----------------------------- #
# fit_4a for suzaku 
proc fit_sub_suz_4a {sfile suff} {
    para err 5; para lev 5
    @fits/fit_4a_base.xcm
    da $sfile; ign 0.0-0.35,10.-**
    set ifree [az_free]
    foreach i {31 35 39 43 47 51 55 59 63} {
        new $i,,,.3 .3 4.5 4.5
    }
    thaw 6 7 31 34 35 38 39 42 43 46 47 50 51 54 55 58 59 62 63 66; fit 100
    freez 1-**; thaw $ifree; fit
    az_calc_errors [az_free_params] fits/fit_sub_suz_4a$suff 1.0
}

########################



###############################
### tests for a broad lines ###
###############################


# ------------------------- #
# similar to published work #
# po + xillver 
proc fit_re_1 {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    da $sfile; ign 0.0-2.,10.-**
    mo tbabs*(tbpcf*(po+xillver)) & 0.0198 -1 & \
        2 & 1 -1 & 3.319e-3 & 1.5 & 1e-3 & \
        =5 & 1 -1& 1000 -1& 0 -1& 3.319e-3& ,-1& -1 -1& 1e-3
    fit 1000
    az_calc_errors [az_free_params] fits/fit_re_1$suff 1.0
}

# ------------------------- #
# make absorber partial #
proc fit_re_1a {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_re_1$suff
    thaw 3; fit 1000
    az_calc_errors [az_free_params] fits/fit_re_1a$suff 1.0
}

# ------------------------------- #
# Add warm absorbers; if addition
# doesn't imporve fit; set cf=1
proc fit_re_1b {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_re_1a$suff
    
    # add low-xi WA #
    set c0 [tcloutr stat]; set d0 [lindex [tcloutr dof] 0]
    add 3 zxipc & 1 & 1 .1 -1 -1 1.7 1.7 & 1 -1 & 3.319e-3
    fit 1000
    
    # check improvement #
    set c1 [tcloutr stat]; set d1 [lindex [tcloutr dof] 0]
    set fprob 1
    if { $c1 < $c0 } {
        ftest $c1 $d1 $c0 $d0
        set fprob [tcloutr ftest]
    }
    if { $fprob > 0.1 } { new 7 0 -1; freez 5 6 ; fit}

    ##------------##
    # add hi-xi WA #
    set c0 [tcloutr stat]; set d0 [lindex [tcloutr dof] 0]
    add 3 zxipc & 1 & 3.5 .1 3.2 3.2 5 5 & 1 -1 & 3.319e-3
    fit 1000
    
    # check improvement #
    set c1 [tcloutr stat]; set d1 [lindex [tcloutr dof] 0]
    set fprob 1
    if { $c1 < $c0 } {
        ftest $c1 $d1 $c0 $d0
        set fprob [tcloutr ftest]
    }
    if { $fprob > 0.1 } { new 7 0 -1; freez 5 6 }
    
    ##-------------##
    az_calc_errors [az_free_params] fits/fit_re_1b$suff 1.0
}

# ------------------------- #
# use relxill instead of PC #
proc fit_re_1c {sfile {suff ""}} {
    para err 5; para lev 5; query yes
    @fits/fit_re_1$suff
    add 4 relxill & 7 .1&=7&300&,-1&,.1&,.1&&3.319e-3&=5&2.7&,.1&1000&-1 -1& 1e-4
    new 22=17; new 26=11
    new 6 ,,,1e-8 1e-8
    fit 1000
    stepp best 16 0 3 60; fit
    rm fits/fit_re_1c${suff}.xcm >/dev/null 2>&1; save all fits/fit_re_1c$suff
        
    # hi-xi WA
    editm TBabs(TBpcf(zxipc(powerlaw + relxill) + xillver)) &2& 3.5 .1 3.2 3.2 5 5 &1 -1& 3.319e-3
    fit
    rm fits/fit_re_1c1${suff}.xcm >/dev/null 2>&1; save all fits/fit_re_1c1$suff
    
    # lo-xi WA
    editm TBabs(TBpcf(zxipc*zxipc(powerlaw+relxill)+xillver)) &0.05& 1,,-1,-1,1.7 1.7 &1 -1& 3.319e-3
    fit
    stepp best 6 -1 1.7 60
    stepp best 10 2.3 5 60
    stepp best 24 0 3 60; fit
    rm fits/fit_re_1c2${suff}.xcm >/dev/null 2>&1; save all fits/fit_re_1c2$suff
}



#########################
#########################
#########################
####### OLD #############




# ------------------------------------- #
# Towrads a full model; subtset of data #
proc _fit_sub_1 {} {
    # tbabs(tbpc*zxipcf*zxipcf*po + zga+zga + po+zga+zga)
    # to 2, 6, 7, 16, 23
    query yes; para lev 10; para st 20
    
    set obs {2 6 7 16 23}
    set nobs [llength $obs]
    
    for {set i 0} {$i<$nobs} {incr i} {
        da [expr $i+1]:[expr $i+1] spec_[lindex $obs $i].grp
    }
    ign 1-5:0.0-0.3,10.-**
    
    mo tbabs(tbpc*zxipcf*zxipcf*po + zga+zga + po+zga+zga) & \
        0.0198 -1& 15& 0.9& 3.319e-3 & \
        6 & 1.3 & 1 -1& 3.319e-3& 4 & 3.5 -1 & 1 -1& 3.319e-3 & \
        1.6 & 1e-2 & \
        6.4 -1 & 0.02 -.1 & 3.319e-3& 1e-4 & 7.06 -1 & 0 -1 & 3.319e-3& 1e-5 & \
        2.5 & 1e-3& 0.56 & 0 -1& 3.319e-3& 4e-4 & 0.89 & 0 -1 & 3.319e-3& 1e-4
        
    az_untie_spec {2 3 4 5} {2 3 14}
    foreach i [az_find 2] {new $i 15 .1 0.01 0.01 40 40}
    foreach i [az_find 3] {new $i .9 .1 0.1 0.1}
    fit 1000
    return
    stepp best 6 1 2.5 35; fit
    set suff 1a; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 3894/863
    
    # untie norm of zga@6.4 #
    untie [az_find 18]; fit
    set suff 1b; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 3700/859
    
    
    # untie nh of lo-xi zxipc
    untie [az_find 5]; fit
    az_rand_pars; fit
    set suff 1c; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 3028/855
    
    # untie gamma of powerlaw and randomize #
    untie [az_find 13]; fit
    az_rand_pars; fit
    set suff 1d; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2988/851
    
    
    # untie xi of lo-xi zxipc
    untie [az_find 6]; fit
    foreach i [az_find 6] { stepp best $i 1 2.5 35 }; fit
    set suff 1e; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2931/847
    
    
    # untie nh of hi-xi zxipc
    untie [az_find 9]; fit
    set suff 1f; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2825/843
      
    
    
    
    # add a full covering ztbabs
    editm tbabs(ztbabs*tbpcf*zxipcf*zxipcf*po + zga + zga + po + zga + zga) & \
            0.05 0.01 0.001 0.001 & 3.319e-3
    fit
    set suff 1g; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2657/842
    
    
    # thaw xi of h-xi zxipc; we froze it to guide the fit only
    thaw 12; fit; stepp best 12 3 4 50
    set suff 1h; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2647/841
    
    ## thaw cf of lo-xi zxipc and allow it to vary #
    # limit xi of lo-xi zxipc to >1; limit cf of tbpcf of > 0.8
    foreach i [az_find 9] {new $i 1 0.01 .8 .8}; fit
    foreach i [az_find 5] {new $i,,, .8 .8}; fit
    foreach i [az_find 8] { new $i ,,0.01 1 1}; fit
    foreach i [az_find 8] { stepp best $i 1 2.5 45}
    set suff 1i; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2495/836
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_1i 1.0
    
    # ------------------------------------ #
    # ------------------------------------ #
    
    #<< untie kbeta norm; no significant improvement; leave it tied #
    #<< tested blueshift of lo-xi zxipc; nothing
    
    
    for {set ic 1} {$ic<=20} {incr ic} {
        set fprob [_sub_add_soft_line 10 27 8]
        if { $fprob < 0.01 } {
            set suff 1j$ic; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
        } else { break }
    }
    # 2141/834, 2005/832, 1829/830, 1756/828, 1622/826, 1584/824, 1550/622, 1501/820,
    # 1464/818
    
    
    
    # add variable soft lines #
    set fprob [_sub_add_soft_line 10 27 8 1]
    set suff 1k1; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1336/812
    set fprob [_sub_add_soft_line 10 27 8 1]
    set suff 1k2; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1301/806
    set fprob [_sub_add_soft_line 10 27 8]
    set suff 1k3; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1271/804
    set fprob [_sub_add_soft_line 10 27 8 1]
    set suff 1k4; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1229/798
    
    
    
    # thaw/untie sigma of 6.4 line #
    foreach i [az_find 18] {new $i 0.02 0.01}
    fit
    set suff 1l; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 1195/793
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_1l 1.0
   
    
    ###########################################
    
    
    # split the 6.4 line into a constant and variable part #
    add 9 zga & 6.4 -1 & 0 -1 & 3.319e-3 & 1e-5; new 24=0.12*28; fit
    set suff 1m; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 1167/793
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_1m 1.0 1
    
    
    
    # split the 6.4 line into two gaussians, one at 6.4 and the other at ~6.2 #
    @fits/fit_sub_1l
    foreach i [az_find 18] {new $i 0 -1}
    add 8 zga & 6.2 & 0 -1 & 3.319e-3 & 1e-4
    untie [az_find 24]; fit
    set suff 1n; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 1195/792
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_1n 1.0
    
}

# -------------------------------------- #
# sub_1 but use xillver instead of lines #
proc _fit_sub_2 {} {
    query yes; para lev 10; para st 20
    # start from sub_1i; 2495/836
    @fits/fit_sub_1i
    add 7 xillver & & 1 -1& 1000 -1 & 0 -1& 3.319e-3& 30 -1& & 1e-4
    foreach i [az_find 17] {new $i=[expr $i-2]}
    untie [az_find 24]
    del 8; del 8; fit 100
    foreach i [az_find 8] { stepp best $i 1 2.5 35}
    set suff 2a; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 2290/837
    
    # add soft lines #
    for {set ic 1} {$ic<=20} {incr ic} {
        set fprob [_sub_add_soft_line 9 27 8]
        if { $fprob < 0.01 } {
            set suff 2b$ic; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
        } else { break }
    }
    # 2101/835, 1877/833, 1700/831, 1622/829, 1590/827, 1564/825, 1537/823, 1437/821
    
    # add variable soft lines #
    set fprob [_sub_add_soft_line 9 27 8 1]
    set suff 2c1; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1281/815
    set fprob [_sub_add_soft_line 9 27 8 1]
    set suff 2c2; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1245/809
    set fprob [_sub_add_soft_line 9 27 8 1]
    set suff 2c3; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1219/803
    
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_2c 1.0
    # stat: 1214/803
    
    # add a variable absorption line at 9 keV #
    add 8 zga & 9.2 & 0 -1 & 3.319e-3& -1e-6 .1 -1 -1 0 0
    foreach i [az_find 28] { new $i -1e-6 0.01 -1 -1 0 0}
    stepp best 25 7 10 50
    set suff 2d; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # 1187/797
    para err 30; para lev 5; az_calc_errors [az_free] fits/fit_sub_2d 1.0 1
    
}

# --------------------------------------- #
# Towrads a full model; add suzaku sub_2c #
proc fit_sub_3 {} {
    # suz_1 is similar to xmm_2; 
    # suz_2-4 are somewhat similar, but not very similar; use only suz_2, which overlaps
    # with nu_2
    para lev 1; para st 30
    @fits/fit_sub_2c
    da 6:6 spec_fi_2.grp
    ign 6:0.0-0.5,1.5-2.5,10.-**
    da 7:7 spec_2_a.grp
    da 8:8 spec_2_b.grp
    ign 7-8:0.0-5.,79.-**
    return
    #return
    az_untie_spec {6}  {4 5 7 8 9 11 15 16 24 30 34 38}
    
    add 6 cutoffp & 1.7 & 1000 -.1 10 10 1000 1000 & 1e-3 & /*
    foreach i [az_find 15] {new $i=[expr $i+3]}
    foreach i [az_find 17] {new $i=[expr $i+2]}
    new 22=16
    del 7
    foreach i [az_find 18] {new $i=[expr $i-3]}
    set npar 79
    for {set i 1} {$i<=$npar} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    add 1 const & 1 -1 & /*
    foreach i {481 561} {new $i 1. .1}
    
    set ifree [az_free]
    freez 1-400; fit 100
    thaw $ifree; fit 100

    set suff sub_3; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 2697/1586
    
    
    # thaw Ecut
    new 17 200 1.0; fit
    foreach i {9 89 169 249 329 409} {stepp best $i 1. 2.5 30}
    set suff sub_3a; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 2570/1585
    
    
    # thaw inclination #
    new 24 30 .1; fit
    stepp best 24 3 75 30
    foreach j {32 36 40} { foreach i [lrange [az_find $j] 0 5] {new $i,,,1e-8 1e-8} }
    # xi
    foreach i {9 89 169 249 329 409} {stepp best $i 1 2.5 25}
    stepp best 13 3.1 4 20
    set suff sub_3b; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 2518/1584
    
    
    
    # add another blueshifted high-xi zxipc to mainly the model suzaku/nustar data  
    add 5 zxipc & 0.05 0.1 & 3.9 0.1 & 1 -1& -.1 0.01 -.2 -.2 3.319e-3 3.319e-3
    set npar 84
    for {set i 8} {$i<=11} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    untie 8 92 176 260 344 428; fit
    stepp best 9 3.1 4 30
    foreach i {36 40 44 48 52 56 60 64 68 72 76 80 84} {new $i,,,1e-8 1e-8}
    set suff sub_3c; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 2356/1576
    

    
    # add relxill to model the extra suzaku emission #
    editm const*TBabs(zTBabs*TBpcf*zxipcf*zxipcf*zxipcf*(cutoffpl+relxill) + xillver + po \
        + zga+zga+zga+zga+zga+zga+zga+zga+zga+zga+zga+zga+zga) & 3 -.1 & =23& 300 & ,-1& \
        & -10 .1 & &3.319e-3 & & 2.7 & ,-1& & -1 -1& 1e-5
    new 27 [tcloutr par 42]; new 42=27
    foreach i [az_find 31] {new $i=[expr $i-11]}
    new 34=21
    untie [az_find 36]
    set npar 98
    for {set i 23} {$i<=36} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    foreach i {36 134 232 330 428 526} {new $i 1e-5,,1e-8 1e-8}
    fit
    new 526 5e-4; fit
    stepp best 27 3 80 20
    set suff sub_3d; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # 2341/1568
    
    
    # force low incl
    new 27 10 -1; fit
    foreach i {13 111 209 307 405 503 } {stepp best $i 1 2.5 30}
    set suff sub_3e; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # 2418/1569
    
    
    # add another relxill #
    add 9 relxill & 10 .1& =23& 300 & ,-1& & -1.5 .1 & & 3.319e-3& & 0 & 1 -1& &-1 -1& 1e-4
    new 27 [tcloutr par 41]; new 41=27
    foreach i [az_find 31] {new $i=[expr $i-11]}
    new 34=21
    untie [az_find 36]
    set npar 112
    for {set i 23} {$i<=36} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    foreach i {36 148 260 372 484 596} {new $i,,,1e-8 1e-8}\
    stepp best 9 3 4 20
    stepp best 13 3 4 20
    foreach i {13 111 209 307 405 503 } {stepp best $i 1 2.5 25}
    thaw 27; fit; thaw 37 ; fit
    new 27 ,,,3 3 80 80
    stepp best 27 3 80 30
    set suff sub_3f; rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    stepp best 23 3 10 10
    stepp best 32 0 4 30
    # 2253/1558
    para err 30; para lev 1; az_calc_errors [az_free] fits/fit_sub_3f 1.0 1
    
}


# ------------------------------------------------------------ #
# Towrads a full model; Simple model and try to force rel_refl #
proc fit_sub_4 {} {
    para lev 5; para st 30
    
    set obs {2 6 7}
    set nobs [llength $obs]
    
    for {set i 0} {$i<$nobs} {incr i} {
        da [expr $i+1]:[expr $i+1] spec_[lindex $obs $i].grp
    }
    ign 1-$nobs:0.0-2.5,10.-**
    
    mo tbabs*(tbpcf*po + xillver) & 0.0198 -1 & 2 & 1 -1& 3.319e-3 & \
        1.5 .1 1.4 1.4 & 1e-3 & \
        =5& 1 -.1& 1000 -1& 0 -1& 3.319e-3& 30 -1& -1 -1& 1e-3 
    
    foreach i {2 6 14} {untie [az_find $i] }
    foreach i [az_find 7] {new $i=[expr $i-2]}
    return
    fit 100
    set suff 4a; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 839/357
    
    # untie gamma #
    untie [az_find 5]
    set suff 4b; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 712/355
    
    # add relxill #
    editm TBabs(TBpcf*(powerlaw+relxill) + xillver) & 4 .1 & =7 & 300 -1& ,-1& 30 -1 & \
        -1.,.1 & & 3.319e-3 & & 2.7 & 1 -1& 1000 -1& -1 -1 & 1e-4
    foreach i [az_find 15] {new $i=[expr $i-10]}
    new 22=17; new 23=18; new 26=11; unite [az_find 20]; fit
    set suff 4c; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 553/349
    
    thaw 11; fit
    stepp best 16 0 4 30
    set suff 4d; rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 503/348
    para err 20; para lev 5; az_calc_errors [az_free] fits/fit_sub_4d 1.0
    ################
    
    # add datasets 16 and 23
    da 4:4 spec_16.grp; da 5:5 spec_23.grp; ign 4-5:0.0-2.5,10.-**
    az_untie_spec {4 5} {2 5 6 20 28}
    new 127=117; new 133=117; new 105=89; new 99=89
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    set ifree [az_free]; freez 1-84; fit 100; thaw $ifree; fit
    set suff 4e;rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 1299/581
    para err 20; para lev 5; az_calc_errors [az_free] fits/fit_sub_4e 1.0
    
    
    # add hi-zxipc
    editm TBabs(TBpcf*zxipc(powerlaw + relxill) + xillver) & 0.05 & 3.5 & 1 -1 & 3.319e-3
    untie [az_find 5]
    foreach i [az_find 9] {new $i,,,1.4 1.4 1.9 1.9}
    stepp best 20 0 4 30
    stepp best 15 5 80 30
    set suff 4f;rm fits/fit_sub_$suff.xcm >/dev/null 2>&1; save all fits/fit_sub_$suff
    # stat: 1235/575; there is another better solution but the gamma's are >2.3
    
    
    # replace tbpcf with zxipc #
    add 3 zxipc & 0.05 & 1.5 & 1 -1 & 3.319e-3; del 2; untie [az_find 2]; fit
    set suff sub_4g;rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    az_rand_pars;fit
    # stat: 966/574
    
    # add another hi-zxipc with blueshift #
    add 2 zxipc & 0.05 & 3.8 & 1 -1 & -.1 .1 -.2 -.2 3.319e-3 3.319e-3;fit; untie [az_find 2]; fit
    set suff sub_4h;rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 898/567
    
    
    # 4f (using tbpcf) and add a hi-zxipc
    add 3 zxipc & 0.05 & 3.8 .1 3.1 3.1 4 4 & 1 -1 & -.1 .1 -.2 -.2 3.319e-3 3.319e-3
    fit; untie [az_find 5]; fit
    set suff sub_4i;rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 1179/568;
    para err 20; para lev 5; az_calc_errors [az_free] fits/fit_sub_4i 1.0
    
    
    # make nuetral absorber partial #
    thaw 3; unite [az_find 3]; fit
    stepp best 24 0 3 20
    set suff sub_4j;rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 749/563;
    
    # add gaussian at 3.7 keV
    editm TBabs(TBpcf*zxipcf*zxipcf(powerlaw + relxill) + xillver+zga) & \
        3.7 & 0 -1 & 3.319e-3 & 1e-5
    fit 100
    set suff sub_4k;rm fits/fit_$suff.xcm >/dev/null 2>&1; save all fits/fit_$suff
    # stat: 721/561;
    para err 20; para lev 5; az_calc_errors [az_free] fits/fit_sub_4k 1.0
}


# ---------------------------------------- #
# sub_4 (2.5-10 keV) and add suzaku/nustar #
proc fit_sub_5 {} {
    # start with sub_4f
    para lev 1; para st 30
    @fits/fit_sub_4k
    da 6:6 spec_fi_2.grp
    ign 6:0.0-2.5,10.-**
    da 7:7 spec_2_a.grp; da 8:8 spec_2_b.grp
    ign 7-8:0.0-5.,79.-**
    return
    az_untie_spec {6}  {2 3 5 9 13 14 28 36}
    
    # pl -> cutoff
    add 5 cutoffp & 1.7 & 1000 -.1 10 10 1000 1000 & 1e-3 & /*
    foreach i [az_find 13] {new $i=[expr $i+3] }
    foreach i [az_find 15] {new $i=[expr $i+2]}
    foreach i [az_find 14] {new [expr $i+15]=$i; new [expr $i+20]=$i }
    foreach i [az_find 13] {new [expr $i+13]=$i; new [expr $i+19]=$i }
    del 6
    
    set npar 41
    for {set i 2} {$i<=$npar} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    add 1 const & 1 -1 & /*
    untie 253 295; thaw 253 295
    
    
    set ifree [az_free]
    freez 1-210; fit 100
    thaw $ifree; fit 100
    foreach i {16 58 100 142 184 226 30 72 114 156 198 240 } {new $i,,,1e-8 1e-8}
    new 17,,,2.5 2.5
    new 21 20 -.1; fit; thaw 21;fit
    stepp best 26 0 3 30
    set suff fit_sub_5a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1915/1311
    
    
    # add wa with fixed xiy #
    add 4 zxipc & 0.05 & 1.5 .1 1. 1. 2.5 2.5 & 1 -1 & 3.319e-3
    untie [az_find 6]
    set npar 46
    foreach i {6 7 8 9} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    fit
    stepp best 30 0 3 20; fit
    set suff fit_sub_5b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1803/1304


    # add another relxill to model suzaku 
    add 8 relxill & 4. .1 & & 300 -1 & =38 & =39& -10 .1 & 400 & 3.319e-3 & & 2.5 &=45 \
            &&-1 -1& 1e-4 
    foreach i [az_find 22] {new $i=[expr $i-1] }
    foreach i [az_find 29] {new $i=[expr $i-11] }
    foreach i [az_find 32] {new $i=[expr $i-13] }
    foreach i [az_find 34] {new $i 1e-8 -.1 1e-8 1e-8}
    set npar 60
    for {set i 21} {$i<=34} {incr i} {
        new [expr $npar*6+$i]=[expr $npar*5+$i]
        new [expr $npar*7+$i]=[expr $npar*5+$i]
    }
    
    thaw 334
    stepp best 30 0 3 20
    set suff fit_sub_5c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1744/1300
    
    
    
    # fitting is slow; reparameterize  #
    new 24 .998 -1; new 25 3. .1 3 3 80 80; new 31 1 -1
    new 38=24; new 39=25; new 45=31
    new 50=31; new 54=25
    fit
    set suff fit_sub_5d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1702/1300
    
    # thaw noms of 2nd relxill in spectra other than suzaku
    foreach i {34 94 154 214 274} {new $i 1e-4 .1}
    renorm; fit
    stepp best 7 0 2.5 30
    set suff fit_sub_5e; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1684/1295
    
    # thaw Ec and error #
    thaw 19
    para err 30; para lev 1; az_calc_errors [az_free] fits/fit_sub_5e 1.0
    # stat: 1680/1294
    
    # start noticing the soft data #
    no 1-5:.3-.4; ign 1-5:0.0-.3,.4-1.; fit
    set suff fit_sub_5f1; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1726/1301
    no 1-6:2.-3.;fit
    set suff fit_sub_5f2; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 1988/1374
    no 1-6:1.5-3.;fit
    add 12 zga & 1.9 & 0 -1 & 3.319e-3& 1e-5; fit
    add 12 zga & 2.3 & 0 -1 & 3.319e-3& 1e-5; fit
    set suff fit_sub_5f3; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 2275/1447
    no 1-6:1.-3.; fit
    add 12 zga & 1.2 & 0 -1 & 3.319e-3& 1e-5; fit
    set suff fit_sub_5f4; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    no 1-6:.8-2.
    add 12 zga & .89 & 0 -1 & 3.319e-3& 1e-4; fit
    set suff fit_sub_5f5; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    no 1-6:.5-1.
    add 12 zga & .56 & 0 -1 & 3.319e-3& 4e-4; fit
    untie [az_find 62]; new 554=472; new 636=472; fit
    add 12 zga & .7 & 0 -1 & 3.319e-3& 1e-5; fit
    untie [az_find 62]; new 578=492; new 664=492; fit
    set suff fit_sub_5f6; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    no 1-5:.35-1.; no 6:.5-1.
    untie [az_find 70]; new 586=500; new 672=500
    foreach i {62 148 234 320 406 492} {new $i,,,1e-8 1e-8}
    foreach i {66 152 238 324 410 496} {new $i,,,1e-8 1e-8}
    foreach i {70 156 242 328 414 500} {new $i,,,1e-8 1e-8}
    set suff fit_sub_5f7; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 3081/1611
    
    stepp best 7 1 2.5 40
    stepp best 11 3.1 4 30
    stepp best 15 3.1 4 30
    stepp best 30 0 4 40
    stepp best 44 0 4 40
    stepp best 19 100 1000 30
    set suff fit_sub_5f; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 3064/1611
    
    # untie xi of main WA
    untie 7 93 179 265 351 437; fit
    foreach i {7 93 179 265 351 437} {stepp best $i 1 2.5 30}
    set suff fit_sub_5g; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 3011/1606
    para err 40; para lev 1; az_calc_errors [az_free] fits/fit_sub_5g 1.0
    
}


# ---------------------------------------------- #
# Simple tests for the presenece of a broad line #
proc fit_sub_6 {} {
    para lev 1; para st 30
    
    # tbabs*(tbpc*po+xillver); untied gamma
    @fits/fit_sub_4b
    
    # add datasets 16 and 23
    da 4:4 spec_16.grp; da 5:5 spec_23.grp; ign 4-5:0.0-2.5,10.-**
    az_untie_spec {4 5} {2 5 6 14}
    new 49=47; new 63=61; fit 1000
    set suff fit_sub_6a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 3132/590
    
    
    # add a constant powerlaw for the ~3 keV excess
    editm TBabs(TBpcf*powerlaw + xillver + po) & 2 & 1e-3; fit
    set suff fit_sub_6b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 1137/588
    
    
    # add a gaussian emission and absorption
    add 5 zga & 6.2 & 0.3 & 3.319e-3 & 1e-5
    add 5 zga & 6.7 & 0 -1 & 3.319e-3 & 1e-5
    foreach i [az_find 22] {new $i 1e-5 0.01 1e-8 1e-8 1 1}
    foreach i [az_find 18] {new $i -1e-5 0.01 -1 -1 -1e-8 -1e-8}
    fit
    set suff fit_sub_6c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # 905/575
    
}


# ---------------------------------- #
# Descriptive model; subtset of data #
proc fit_sub_7 {ispec} {
    da spec_$ispec.grp
    ign 0.0-6.8,10.-**
    mo ztbabs*po & 3 & 3.319e-3 & 1.6 & 1e-3; fit
    data none
    exec grppha spec_$ispec.pha !tmp_$ispec.grp "chkey respfile spec_$ispec.rmf& \
            chkey ancrfile spec_$ispec.arf&chkey backfile spec_$ispec.bgd&exit"
    da tmp_$ispec.grp
    setpl ch; setpl noarea
    
    # absorbed po #
    new 3 0; new 4 1
    set x [tcloutr plot d x 1]
    set y1 [tcloutr plot d mod 1]
    # unabsorbed po #
    del 1
    set y2 [tcloutr plot d mod 1]
    set txt ""
    for {set i 0} {$i<[llength $x]} {incr i} {
        set a1 [lindex $y1 $i]
        set a2 [lindex $y2 $i]
        if {$a2 == 0} {
            set r 1
        } else {
            set r [expr $a1/$a2]
        }
        append txt "[lindex $x $i] $r\n"
    }
    set fp [open tmp_sub_7.dat w]
    puts $fp $txt
    close $fp
}



# ----------------------------------------------------- #
# Relativistic ref. in published paper to sub spectra   #
# start with the old spectra only, similar to published #
# we use residuals to illustrate the point              #
proc fit_sub_re_1 {} {
    query yes; para lev 1; para st 30
    
    set obs {1 4}
    set nobs [llength $obs]
    
    for {set i 0} {$i<$nobs} {incr i} {
        da [expr $i+1]:[expr $i+1] spec_[lindex $obs $i].grp
    }
    ign 1-$nobs:0.0-2.5,10.-**
    mo tbabs*(tbpcf*po + zga +zga) & 0.0198 -1 & 2 & 1 -1 & 3.319e-3 & \
        1.6 & 7e-3 & 6.4 -1 & 0.02 -1 & 3.319e-3& 1e-4 & 7.06 -1 & 0 -1 & 3.319e-3& 1e-5
    
    return
    az_untie_spec {2} {2 6 10}
    fit
    set suff fit_sub_re_1a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 524/235
    
    # add 16 23 
    da 3:3 spec_16.grp; da 4:4 spec_23.grp; ign 3-4:0.0-2.5,10.-**
    az_untie_spec {3 4} {2 6 10}
    fit
    set suff fit_sub_re_1b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 4582/472
    
    # add a const po for the extra component
    editm TBabs(TBpcf*po +zga+zga+po) & .9 & 1e-3; fit
    set suff fit_sub_re_1c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1128/470
    
    # make model pc instead of adding a po
    @fits/fit_sub_re_1b
    new 3 .99 .1; untie [az_find 3]; fit
    set suff fit_sub_re_1d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 961/468
    
    # add zxipc #
    add 3 zxipc & 1 & 3.5 .1 3 3 & 1 -1 & 3.319e-3
    untie [az_find 5]
    set suff fit_sub_re_1e; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 813/463
}

# ----------------------------------------------------- #
# similar to _re_1, but using xillver instead of guass  #
proc fit_sub_re_2 {} {
    query yes; para lev 1; para st 30
    
    set obs {1 4}
    set nobs [llength $obs]
    
    for {set i 0} {$i<$nobs} {incr i} {
        da [expr $i+1]:[expr $i+1] spec_[lindex $obs $i].grp
    }
    ign 1-$nobs:0.0-2.5,10.-**
    mo tbabs*(tbpcf*po+xillver) & 0.0198 -1 & 2 & 1 -1 & 3.319e-3 & \
        1.6 & 7e-3 & \
        =5 & 1 -1 & 1000 -1 & 0 -1 & 3.319e-3& 30 -1 & -1 -1 & 1e-3
    
    return
    az_untie_spec {2} {2 6 14}
    fit
    set suff fit_sub_re_2a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 359/236
    
    # add 16 23 
    da 3:3 spec_16.grp; da 4:4 spec_23.grp; ign 3-4:0.0-2.5,10.-**
    az_untie_spec {3 4} {2 6 14}
    fit
    set suff fit_sub_re_2b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2904/473
    
    # add a const po for the extra component
    editm TBabs(TBpcf*po +xillver+po) & .9 & 1e-3; fit
    set suff fit_sub_re_2c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 959/471
    
    # make model pc instead of adding a po
    @fits/fit_sub_re_2b
    new 3 .99 .1; untie [az_find 3]; fit
    set suff fit_sub_re_2d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 961/469
    
    # add zxipc #
    add 3 zxipc & 1 & 3.5 .1 3 3 & 1 -1 & 3.319e-3
    untie [az_find 5]
    set suff fit_sub_re_2e; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 705/464
    
    # untie gamma #
    untie [az_find 9]
    foreach i [az_find 9] {new [expr $i+2]=$i}
    fit
    set suff fit_sub_re_2f; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 632/461
}




# ----------------------------------------------------- #
# Relativistic ref. in published paper to sub spectra   #
# start with the old spectra only, similar to published #
# we use residuals to illustrate the point              #
proc _fit_sub_re_1 {} {
    query yes; para lev 1; para st 30
    
    set obs {1 2 3 4 5 6 7 8}
    set nobs [llength $obs]
    
    for {set i 0} {$i<$nobs} {incr i} {
        da [expr $i+1]:[expr $i+1] spec_[lindex $obs $i].grp
    }
    ign 1-$nobs:0.0-2.5,10.-**
    mo tbabs*(ztbabs*(cutoffp+relxill+xillver)) & 0.0198 -1 & 2 & 3.319e-3 & \
        1.6 & 1000 -1 20 20 1000 1000 & 7e-3 & \
        5. -1&=7 & 300 & ,-1& 40 -1 & -2 .1 & & 3.319e-3& =4 & 2.8 & 1 -1& =5& -1 -1 & 6e-4 & \
        =4& =17& =5& 0 -1& 3.319e-3& =11& -1 -1& 1e-3 & /* 
    return
    az_untie_spec {2 3 4 5 6 7 8} {2 6 20 28}
    fit 100
    thaw 11; fit
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    foreach i [az_find 20] {new $i,,,1e-8 1e-8}
    set suff fit_sub_re_1a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1230/937
    
    
    # fix xi at 2.8; so things are closer to the Z12 model
    new 16 2.8 -1; fit; thaw 16; fit
    set suff fit_sub_re_1b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1244/937
    
    # untie xi
    untie [az_find 16]; fit
    for i [az_find 16] {stepp best $i 2 3.5 10}
    set suff fit_sub_re_1c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1150/930
    
    ########################
    # add datasets 16 and 23
    da 9:9 spec_16.grp; da 10:10 spec_23.grp; ign 9-10:0.0-2.5,10.-**
    az_untie_spec {9 10} {2 6 16 20 28}
    set ifree [az_free]; freez 1-224; fit 100
    set suff fit_sub_re_1d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2655/1206
    
    # a const powerlaw
    editm TBabs(zTBabs(cutoffpl + relxill + xillver) + po) & 0.8 -1 & 0 -1 & /*
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    new 299=269; new 300=270
    untie 269 270; thaw 269 270
    fit
    foreach i {256 286} {stepp best $i 2 3.5 10}
    set suff fit_sub_re_1e; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1724/1204
    
    # replace extra po with patrtial covering
        editm TBabs(zTBabs*tbpcf(cutoffpl + relxill + xillver) + po) & 3 & 1 -1 & 3.319e-3
        foreach i [az_find 4] {new $i [tcloutr par [expr $i-2]]}
        del 2; del 6
        new 235 .98 .1 .8 .8; new 264 .98 .1 .8 .8
        fit
        set suff fit_sub_re_1e1; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
        # stat: 1681/1204
        
        # add zxipc 
        add 2 zxipc & 0.05 -1 & 3.5 -1.1 3 3 & 0 -1 & 3.319e-3
        untie 266 299; thaw 266 299; thaw 3
        new 268 1 -1; new 301=268; fit
        set suff fit_sub_re_1e2; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
        # stat: 1589/1201
    
    
    # continue from sub_re_1e
    # add 21 and 22
    da 11:11 spec_21.grp; da 12:12 spec_22.grp; ign 11-12:0.0-2.5,10.-**
    az_untie_spec {11 12} {2 6 16 20 28}
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    foreach i [az_find 20] {new $i,,,1e-8 1e-8}
    thaw 302 306 316 320 328 332 336 346 350 358
    new 359=269; new 360=270
    new 329=269; new 330=270
    fit
    foreach i {256 286 316 346} {stepp best $i 2 3.5 10}
    set suff fit_sub_re_1f; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2139/1440
    
    # add zxipc 
    add 2 zxipc & 0.05 -1 & 3.5 -1.1 3 3 & 0 -1 & 3.319e-3
    untie 274 308 342 376; thaw 274 308 342 376
    foreach i {309 343 377} {new $i=275}
    new 275 3.5 .1 3 3
    foreach i {310 344 378} {new $i=276}
    new 276 1 -1
    set suff fit_sub_re_1g; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1996/1437
    
    
    # notice more soft data
    no 1-12:2.-3.; ign 1-12:0.0-2.; fit
    set suff fit_sub_re_1h; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 5294/1570
}


# --------------------------------------------- #
# similar to sub_re_1, unfix the old parameters #
proc _fit_sub_re_2 {} {
    # start from fit_sub_re_1h
    @fits/fit_sub_re_1h
    
    # thaw the variable parameters
    untie 6 10 20 24 32 40 44 54 58 66 74 78 88 92 100 108 112 122 126 134 142 146 156 160 168 176 180 190 194 202 210 214 224 228 236 244 248 258 262 270
    
    new 33 0.62 .1; new 34 7e-4 .1
    foreach i { 305 339 373 407} {new $i=33}
    foreach i { 306 340 374 408} {new $i=34}
    
    new 3 3.35 .1 3 3 ; new 4 1 -1; 
    foreach i {275 309 343 377} {new $i=3}
    foreach i {276 310 344 378} {new $i=4}
    untie [az_find 2]; thaw [az_find 2]
    foreach i [az_find 20] {new $i,,,.1 .1}
    thaw 8 15 16 
    set suff fit_sub_re_2a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 3186/1517
    
    
    # reposition xillver outside absorbers #
    editm TBabs(zxipcf*zTBabs(cutoffpl + relxill + xillver) + xillver + powerlaw) & \
        =8 & =21& =9& 0 -1 & 3.319e-3& =15& -1 -1 & 1e-3 
    foreach i [az_find 40] {new $i=[expr $i-8]}
    del 6; fit
    foreach i [az_find 20] {stepp best $i 2 3.5 10}
    set suff fit_sub_re_2b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2556/1517
    
    
    # replace extra po with partial covering #
    add 3 tbpc & 3 & 1 -1 & 3.319e-3
    foreach i [az_find 6] {new $i [tcloutr par [expr $i+3]]}
    del 4
    del 7 
    new 7 .98;
    foreach i [az_find 7] {new $i .98 .1 .8 .8}
    fit
    foreach i [az_find 21] {stepp best $i 2 3.5 10}
    set suff fit_sub_re_2c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2267/1507
    
    
    # add WA
    add 3 zxipc & 1 & 1.5 & 1 -1 & 3.319e-3
    untie [az_find 6]; fit
    set suff fit_sub_re_2d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1955/1494
}


# ---------------------------------------------------- #
# similar to sub_re_1, fix the parameters of old data  #
proc _fit_sub_re_3 {} {
    # start from fit_sub_re_1b
    @fits/fit_sub_re_1c
    
    # add datasets 16 and 23
    da 9:9 spec_16.grp; da 10:10 spec_23.grp; ign 9-10:0.0-2.5,10.-**
    az_untie_spec {9 10} {2 6 16 20 28}
    set ifree [az_free]; freez 1-224; fit 100
    set suff fit_sub_re_3a; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 3657/1208
    
    # a const powerlaw
    editm TBabs(zTBabs(cutoffpl + relxill + xillver) + po) & 0.8 -1 & 0 -1 & /*
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    new 299=269; new 300=270
    untie 269 270; thaw 269 270
    fit
    set suff fit_sub_re_3b; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 1853/1206
    
    
    # add 21 and 22
    da 11:11 spec_21.grp; da 12:12 spec_22.grp; ign 11-12:0.0-2.5,10.-**
    az_untie_spec {11 12} {2 6 20 28}
    foreach i [az_find 6] {new $i,,,1e-8 1e-8}
    foreach i [az_find 20] {new $i,,,1e-8 1e-8}
    thaw 302 306 320 328 332 336 350 358
    new 359=269; new 360=270
    new 329=269; new 330=270
    fit
    set suff fit_sub_re_3c; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2328/1444
    
    # add zxipc 
    add 2 zxipc & 0.05 -1 & 3.5 -1.1 3 3 & 0 -1 & 3.319e-3
    untie 274 308 342 376; thaw 274 308 342 376
    foreach i {309 343 377} {new $i=275}
    new 275 3.5 .1 3 3
    foreach i {310 344 378} {new $i=276}
    new 276 1 -1
    set suff fit_sub_re_3d; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 2135/1439
    
    # notice more soft data
    no 1-12:2.-3.; ign 1-12:0.0-2.; fit
    set suff fit_sub_re_3e; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 5629/1572
    
    # notice more soft data
    no 1-12:1.5-3.; ign 1-12:0.0-1.5; fit
    set suff fit_sub_re_3f; rm fits/${suff}.* >/dev/null 2>&1; save all fits/$suff
    # stat: 5629/1572
    
}




# ---------------------------------------------------------- #
# setup fit_indiv_?; distribute the ispec to the closest sub #
proc _fit_indiv_base {sfile {sub 0}} {
    para err 8; para st 8; para lev 5; query yes;
    # which starting parameters to use for each ispec #
    if {$sub == 0} {
        set obs  {1 2 3 4 5 6 7 8 10 11 12 13 15 16 17 18 19 20 21 22 23 24}
        set igrp {0 0 0 1 1 1 2 2 1 1 2 2 1 3 3 3 3 3 3 4 4 4}
    } else {
        set obs {}; set igrp {}
        for {set i 1} {$i<=81} {incr i} {
            lappend obs $i
            if {$i<=17} { lappend igrp 0
            } elseif {($i<=23) || ($i==33) || ($i==34) || ($i==37) } { lappend igrp 1
            } elseif {$i<=36} { lappend igrp 2
            } elseif {$i<=65} { lappend igrp 3
            } else {lappend igrp 4}
        } 
    }
    set ispec [lindex [split [string trimright $sfile .grp] _] 1]
    set idx [lindex $igrp [lsearch $obs $ispec]]
    
    set pars_list [split [az_save_pars] &]
    set pars_list [lreplace $pars_list  0 0]
    set nmod [expr [tcloutr modpar]/[tcloutr datagrp]]
    
    set a0 [lrange $pars_list 0 [expr $nmod-1]]
    set pars_grp {}
    lappend pars_grp [join $a0 &]
    for {set i 1} {$i<5} {incr i} {
        set a1 [lrange $pars_list [expr $i*$nmod] [expr ($i+1)*$nmod-1]]
        for {set j 0} {$j<[llength $a0]} {incr j} {
            if { [string range [lindex $a1 $j] 0 0] eq "=" } {
                lset a1 $j [lindex $a0 $j]
            }
        }
        lappend pars_grp [join $a1 &]
    }
    da $sfile; ign 1:0.0-0.3,10.-**
    newpar 1-$nmod [lindex $pars_grp $idx]

}



# --------------------------------------------------- #
# Fit the best in sub_1l to all spectra, individually #
proc fit_indiv_1l {sfile {suff ""} {sub 0}} {

    # load the relevant sub fit #
    @fits/fit_sub_1l
    
    # setup the model base #
    _fit_indiv_base $sfile $sub
    
    # do other tweaking #
    freez 2 12 24-**; thaw 30 38 42
    foreach i {30 38 42} {
        set n [lindex [tcloutr par $i] 0]
        if {$n<1e-8} {set n 1e-8}
        new $i $n 0.01 1e-8 1e-8
    }
    #new 24 1e-5 .01 1e-8 1e-8
    new 18 ,,,1e-4 1e-4
    fit 1000
    stepp best 8 1.0 2.5 35
    
    # switch to fluxes #
    add 6 cflux & 7 & 10 & -10
    new 19 1 -1; fit 1000
    add 8 cflux & 2 & 10 & -11.5
    new 26 1 -1; fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_1l$suff 1.0
}
proc fit_indiv_1l_mc {sfile suff} { run_mcmc fit_indiv_1l $suff }


# --------------------------------------------------- #
# Fit the best in sub_1n to all spectra, individually #
proc fit_indiv_1n {sfile {suff ""} {sub 0}} {

    # load the relevant sub fit #
    @fits/fit_sub_1n
    
    # setup the model base #
    _fit_indiv_base $sfile $sub
    
    # do other tweaking #
    freez 2 12 21 28-**; thaw 34 42 46 
    foreach i {34 42 46} {
        set n [lindex [tcloutr par $i] 0]
        if {$n<1e-8} {set n 1e-8}
        new $i $n 0.01 1e-8 1e-8
    }
    new 24 ,,,1e-8 1e-8
    fit 1000
    stepp best 8 1.0 2.5 35
     
    # switch to fluxes #
    add 6 cflux & 7 & 10 & -10; new 19 1 -1; fit 1000
    add 8 cflux & 2 & 10 & -11.5; new 26 1 -1; fit
    add 10 cflux & 2 & 10 & -11.5; new 33 1 -1; new 29 -11.5 .1 -14 -14;fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_1n$suff 1.0
}
proc fit_indiv_1n_mc {sfile suff} { run_mcmc fit_indiv_1n $suff }


# --------------------------------------------------- #
# Fit the best in sub_2c to all spectra, individually #
proc fit_indiv_2c {sfile {suff ""} {sub 0}} {

    if {[string range $sfile 5 6] == "fi" } {
        # suzaku data #
        @fits/fit_indiv_2c_base
        da $sfile; ign 1:0.0-.5,10.-**
    } else {
        # load the relevant sub fit #
        @fits/fit_sub_2c
        # setup the model base #
        _fit_indiv_base $sfile $sub
    }
    
    # do other tweaking #
    freez 2 12 25-**; thaw 30 34 38
    foreach i {30 34 38} {
        set n [lindex [tcloutr par $i] 0]
        if {$n<1e-8} {set n 1e-8}
        new $i $n 0.01 1e-8 1e-8
    }
    fit 1000
    stepp best 8 1.0 2.5 35
    
    # switch to fluxes #
    add 6 cflux & 7 & 10 & -10
    new 19 1 -1; fit 1000
    add 8 cflux & 6.1 & 6.7 & -11.5
    new 30 1 -1; fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_2c$suff 1.0
}
proc fit_indiv_2c_mc {sfile suff} { run_mcmc fit_indiv_2c $suff }


# --------------------------------------------------- #
# Fit the best in sub_2d to all spectra, individually #
proc fit_indiv_2d {sfile {suff ""} {sub 0}} {

    # load the relevant sub fit #
    @fits/fit_sub_2d
    
    # setup the model base #
    _fit_indiv_base $sfile $sub
    
    # do other tweaking #
    freez 2 12 25-**; thaw 28 34 38 42
    foreach i {34 38 42} {
        set n [lindex [tcloutr par $i] 0]
        if {$n<1e-8} {set n 1e-8}
        new $i $n 0.01 1e-8 1e-8
    }
    set n [lindex [tcloutr par 28] 0]
    if {$n>-1e-8} {set n -1e-8}
    new 28 $n 0.01 -1 -1 -1e-8 -1e-8
    
    fit 1000
    stepp best 8 1.0 2.5 35
    
    # switch to fluxes #
    add 6 cflux & 7 & 10 & -10
    new 19 1 -1; fit 1000
    add 8 cflux & 6.1 & 6.7 & -11.5
    new 30 1 -1; fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_2d$suff 1.0
}
proc fit_indiv_2d_mc {sfile suff} { run_mcmc fit_indiv_2d $suff }



# --------------------------------------------------- #
# Fit the best in sub_1m to all spectra, individually #
# uses a variable and constant line at 6.4. keV       #
proc fit_indiv_2 {sfile {suff ""} {sub 0}} {

    # load the relevant sub fit #
    @fits/fit_sub_1m
    
    # setup the model base #
    _fit_indiv_base $sfile $sub

    # do other tweaking #
    freez 2 12 28-**; thaw 74 78
    #untie 24; new 24 ,,0.01 1e-8 1e-8
    fit 1000
    stepp best 8 1 3 70
    
    add 6 cflux & 7 & 10 & -10
    new 19 1 -1; fit 1000
    
    add 8 cflux & 2 & 10 & -11.5
    new 26 1 -1; fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_2$suff 1.0
}
proc fit_indiv_2_mc {sfile suff} { run_mcmc fit_indiv_2 $suff }


# --------------------------------------------------- #
# Fit the best in sub_1n to all spectra, individually #
# sigma of 6.4 line is 0; add another line at ~6.2    #
proc fit_indiv_3 {sfile {suff ""}} {

    # load the relevant sub fit #
    @fits/fit_sub_1n
    
    # setup the model base #
    _fit_indiv_base $sfile

    # do other tweaking #
    freez 2 12 21 28-**; thaw 74 78
    fit 1000
    stepp best 8 -2 3 130
    
    add 6 cflux & 7 & 10 & -10
    new 19 1 -1; fit
    
    editm tbab(ztbab*tbpcf*zxipc*zxipc*cflux*powerlaw + cflux*(zga+zga) + zga+po+zga+zga + \
        zga+zga+zga+zga+zga+zga+zga+zga+zga+zga) & 2 & 10 & -11.5
    new 26 1 -1; new 30 0.1
    
    # log-transform the norm of the 6.15 line #
    add 25 po & -1 0.1 -4 -4 0 0 & 0 -1
    new 30=10.**p85; fit
    
    az_calc_errors [az_free_params] fits/fit_indiv_3$suff 1.0
}
proc fit_indiv_3_mc {sfile suff} { run_mcmc fit_indiv_3 $suff }



# ----------------------------------------------------- #
# Fit the best in sub_6c to all spectra, individually #
proc fit_indiv_6c {sfile {suff ""} {sub 0}} {

    # load the relevant sub fit #
    @fits/fit_sub_6c
    
    # setup the model base #
    _fit_indiv_base $sfile $sub
    
    
    # do other tweaking #
    ign 1:0.0-2.5,10.-**
    freez 15 19 20 23 24
    fit 1000    
    az_calc_errors [az_free_params] fits/fit_indiv_6c$suff 1.0
}



