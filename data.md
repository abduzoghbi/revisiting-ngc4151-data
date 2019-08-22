
Using:
- `heasoft 6.25`
- `xspec 12.10.1`
- `xmmsas_20180620_1732-17.0.0`
- `ccf` update `20-12-2018`

# XMM Data

## Download the data
- Use xamin to create a list of all obsids and save them to data/xmm/obsids.txt, with obsids being the first column


```python
import numpy as np
import os
import subprocess as subp
import glob
import time
import re
import aztools as az
from ftplib import FTP
from astropy.io import fits as pyfits
import astropy.time as atime
```


```python
base_dir = '/home/abzoghbi/data/ngc4151/spec_analysis'
data_dir = 'data/xmm'
os.system('mkdir -p %s'%data_dir)
obsids = ['0112310101', '0112310501', '0112830201', '0112830501', '0112830601',
          '0143500101', '0143500201', '0143500301', '0402660101', '0402660201',
          '0402660301', '0657840101', '0657840201', '0657840301', '0657840401',
          '0657840501', '0679780101', '0679780201', '0679780301', '0679780401',
          '0679780501', '0761670101', '0761670201', '0761670301', '0761670401',
          '0761670501', '0761670601', '0761670701', '0761670801', '0761670901']
obsids = np.array(obsids)
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 30 observations
    0112310101, 0112310501, 0112830201, 0112830501, 0112830601, 0143500101, 0143500201, 0143500301, 0402660101, 0402660201, 0402660301, 0657840101, 0657840201, 0657840301, 0657840401, 0657840501, 0679780101, 0679780201, 0679780301, 0679780401, 0679780501, 0761670101, 0761670201, 0761670301, 0761670401, 0761670501, 0761670601, 0761670701, 0761670801, 0761670901


- We use `ftplib` to get the data from heasarc (may take some time)


```python
os.chdir(base_dir)
ftp = FTP('legacy.gsfc.nasa.gov', 'anonymous', 'anonymous@gmail.com')
ftp.cwd('xmm/data/rev0')
failed = []
for o in obsids:
    tar_file = '%s/%s.tar'%(data_dir, o)
    # download file only if not already downloaded
    if not os.path.exists(tar_file):
        try:
            ftp.retrbinary('RETR %s.tar'%o ,open(tar_file, 'wb').write)
        except:
            print('failed downloading %s'%o)
            os.system('rm %s >/dev/null 2>&1'%tar_file)
            failed.append(o)

```

    failed downloading 0657840501



```python
for f in failed:
    obsids = np.delete(obsids, np.argwhere(obsids==f)[0,0])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 29 observations
    0112310101, 0112310501, 0112830201, 0112830501, 0112830601, 0143500101, 0143500201, 0143500301, 0402660101, 0402660201, 0402660301, 0657840101, 0657840201, 0657840301, 0657840401, 0679780101, 0679780201, 0679780301, 0679780401, 0679780501, 0761670101, 0761670201, 0761670301, 0761670401, 0761670501, 0761670601, 0761670701, 0761670801, 0761670901


## Process the PN data
We use our shell script `xmm_process`. Split it into two parts so we can run things in parallel across observations. The first creates `ccf.cif`, and the second creates the event files


```python

os.chdir('%s/%s'%(base_dir, data_dir))
os.system('mkdir -p log')
procs = []
for o in obsids:
    if os.path.exists(o): continue
    os.system('tar -xf %s.tar'%o)
    os.chdir(o)
    os.system('rm -r 3XMM om_mosaic PPS >/dev/null 2>&1')
    os.system('mv ODF odf')
    os.chdir('odf')
    if not os.path.exists('ccf.cif'):
        os.system('gzip -d *gz')
        log_file = '../../log/%s_process.log'%o
        proc = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process > %s 2>&1'%log_file])
        procs.append(proc)
    os.chdir('../..')

# wait for the tasks to end
for p in procs: p.wait()
```


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in obsids:
    os.chdir(o)
    os.system('mkdir -p pn')
    os.chdir('pn')
    if len(glob.glob('*EVL*')) == 0 and len(glob.glob('pn.fits')) == 0:
        log_file = '../../log/%s_process_pn.log'%o
        p = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process pn > %s 2>&1'%log_file])
        procs.append(p)
    os.chdir('../..')

# wait for the tasks to end
for p in procs: p.wait()
```

## Spectral Extraction
### Standard Filtering & Region
- `xmm_filter.py` does standard background filtering and opens `ds9` and requrest a region file called `ds9.reg`, which contrains the source and background regions in **Physical** coordinate units, so `xmm_spec.py` can understand it.
- Here, we use an annular region for the source, so later we can check for pileup; at this stage we set the inner radius to 0 and outer radius to 50 arcsec. The background is circular with radius of 50 arcsec.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
for o in obsids:
    print('-- obs %s --'%o)
    os.chdir('%s/pn'%o)
    os.system('ln -s %s pn.fits'%glob.glob('*EVL*')[0])
    os.system('mkdir -p spec')
    os.chdir('spec')
    if not exists('pn_filtered.fits') or not exists('ds9.reg'):
        # check if we have a saved region file, or a temporary region file
        # for faster loading
        saved_reg = '../../../log/%s_ds9.reg'%o
        if exists(saved_reg):
            os.system('cp %s ds9.reg'%saved_reg)
            region = ''
        else:
            region = '--region'
            tmp_reg = '../../../log/ds9.reg'
            if exists(tmp_reg):
                os.system('cp %s tmp.reg'%tmp_reg)

        subp.call(['/bin/bash', '-i', '-c', 
                'sasinit; xmm_filter.py ../pn.fits pn --std %s'%region])
        if not exists(saved_reg):
            os.system('cp ds9.reg %s'%tmp_reg) 
    os.chdir('../../..')
    
```

    -- obs 0112310101 --
    -- obs 0112310501 --
    -- obs 0112830201 --
    -- obs 0112830501 --
    -- obs 0112830601 --
    -- obs 0143500101 --
    -- obs 0143500201 --
    -- obs 0143500301 --
    -- obs 0402660101 --
    -- obs 0402660201 --
    -- obs 0402660301 --
    -- obs 0657840101 --
    -- obs 0657840201 --
    -- obs 0657840301 --
    -- obs 0657840401 --
    -- obs 0679780101 --
    -- obs 0679780201 --
    -- obs 0679780301 --
    -- obs 0679780401 --
    -- obs 0679780501 --
    -- obs 0761670101 --
    -- obs 0761670201 --
    -- obs 0761670301 --
    -- obs 0761670401 --
    -- obs 0761670501 --
    -- obs 0761670601 --
    -- obs 0761670701 --
    -- obs 0761670801 --
    -- obs 0761670901 --


- Some of the observations had no data, so we don't consider them anymore:
```
0112310501, 0112830601, 0402660301, 0679780501: No PN science exposure
0657840101: High background
```
    



```python
no_data = ['0112310501', '0112830601', '0402660301', '0679780501', '0657840101']
obsids = np.array([o for o in obsids if not o in no_data])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 24 observations
    0112310101, 0112830201, 0112830501, 0143500101, 0143500201, 0143500301, 0402660101, 0402660201, 0657840201, 0657840301, 0657840401, 0679780101, 0679780201, 0679780301, 0679780401, 0761670101, 0761670201, 0761670301, 0761670401, 0761670501, 0761670601, 0761670701, 0761670801, 0761670901


### Check for pileup
We use `epatplot` to check for pileup starting with an annulus source region with inner radius of 0. A pileup is present if the fraction of expected to predicted single or doubles events deviates from 1 by more than 3 sigma. If there is pileup, we increase the inner radius of the annulus to 3 arcsec and repeat the increase in steps of 0.5 arcsec, until the fractions are consistent with 1.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
print('obs_num | obsid | fractions | radius | pileup?')
for iobs, o in enumerate(obsids):
    if os.path.exists('%s/pn/pileup'%o): continue
    os.chdir('%s/pn'%o)
    os.system('mkdir -p pileup')
    os.chdir('pileup')
    
    pileup = True
    radius = 0.0
    reg_lines = open('../spec/ds9.reg').readlines()
    if 'background' in reg_lines[-1]:
        reg_lines = reg_lines[:-2] + [reg_lines[-1], reg_lines[-2]]
    g = re.match('\(.*\)\n', reg_lines[-1])
    dum_rad = reg_lines[-1].split(',')[2]
    reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',%g,')
    reg_text = ''.join(reg_lines)
    

    while pileup:
        with open('ds9.reg', 'w') as fp: fp.write(reg_text%radius)
        subp.call(['/bin/bash', '-i', '-c', 
           'sasinit; xmm_spec.py ../pn.fits ds9.reg --check_pileup > pileup.log 2>&1'])
        line = [l for l in open('pileup.log').readlines() if '+/-' in l][0]
        frac = np.array(np.array(line.split())[[2,4,6,8]], np.double)
        pileup = (frac[0] > 1+3*frac[1]) or (frac[2] < 1-3*frac[3])
        # 1 arcsec = 20 pixels
        text = '{:4d} | {} | {} | {} | {:d}'.format(iobs+1, o, line[10:-1], radius/20., pileup)
        print(text)
        if pileup: radius = np.max([3.0*20, radius+10])
    os.chdir('../spec')
    os.system('cp ds9.reg ds9.reg.0')
    os.system('cp ../pileup/ds9.reg ds9.reg')
    os.chdir('../../../')
            
```

    obs_num | obsid | fractions | radius | pileup?
       1 | 0112310101 |  s: 1.005 +/- 0.008   d: 0.990 +/- 0.011 | 0.0 | 0
       2 | 0112830201 |  s: 0.993 +/- 0.005   d: 1.025 +/- 0.008 | 0.0 | 0
       3 | 0112830501 |  s: 0.983 +/- 0.008   d: 1.053 +/- 0.013 | 0.0 | 0
       4 | 0143500101 |  s: 1.006 +/- 0.007   d: 0.980 +/- 0.009 | 0.0 | 0
       5 | 0143500201 |  s: 1.013 +/- 0.007   d: 0.963 +/- 0.009 | 0.0 | 1
       5 | 0143500201 |  s: 1.001 +/- 0.007   d: 0.993 +/- 0.010 | 3.0 | 0
       6 | 0143500301 |  s: 1.009 +/- 0.006   d: 0.977 +/- 0.008 | 0.0 | 0
       7 | 0402660101 |  s: 1.005 +/- 0.007   d: 0.990 +/- 0.010 | 0.0 | 0
       8 | 0402660201 |  s: 1.023 +/- 0.006   d: 0.942 +/- 0.008 | 0.0 | 1
       8 | 0402660201 |  s: 1.012 +/- 0.006   d: 0.971 +/- 0.009 | 3.0 | 1
       8 | 0402660201 |  s: 1.011 +/- 0.007   d: 0.972 +/- 0.009 | 3.5 | 1
       8 | 0402660201 |  s: 1.010 +/- 0.007   d: 0.975 +/- 0.009 | 4.0 | 0
       9 | 0657840201 |  s: 1.003 +/- 0.012   d: 0.998 +/- 0.017 | 0.0 | 0
      10 | 0657840301 |  s: 1.019 +/- 0.009   d: 0.957 +/- 0.012 | 0.0 | 1
      10 | 0657840301 |  s: 1.009 +/- 0.010   d: 0.979 +/- 0.014 | 3.0 | 0
      11 | 0657840401 |  s: 1.004 +/- 0.010   d: 0.982 +/- 0.014 | 0.0 | 0
      12 | 0679780101 |  s: 1.015 +/- 0.011   d: 0.959 +/- 0.016 | 0.0 | 0
      13 | 0679780201 |  s: 1.022 +/- 0.011   d: 0.942 +/- 0.015 | 0.0 | 1
      13 | 0679780201 |  s: 1.006 +/- 0.012   d: 0.982 +/- 0.017 | 3.0 | 0
      14 | 0679780301 |  s: 1.012 +/- 0.012   d: 0.972 +/- 0.016 | 0.0 | 0
      15 | 0679780401 |  s: 1.010 +/- 0.012   d: 0.970 +/- 0.016 | 0.0 | 0
      16 | 0761670101 |  s: 1.018 +/- 0.007   d: 0.953 +/- 0.010 | 0.0 | 1
      16 | 0761670101 |  s: 1.008 +/- 0.008   d: 0.977 +/- 0.011 | 3.0 | 0
      17 | 0761670201 |  s: 1.018 +/- 0.006   d: 0.954 +/- 0.009 | 0.0 | 1
      17 | 0761670201 |  s: 1.005 +/- 0.007   d: 0.989 +/- 0.010 | 3.0 | 0
      18 | 0761670301 |  s: 1.017 +/- 0.006   d: 0.959 +/- 0.009 | 0.0 | 1
      18 | 0761670301 |  s: 1.003 +/- 0.007   d: 0.998 +/- 0.010 | 3.0 | 0
      19 | 0761670401 |  s: 1.016 +/- 0.006   d: 0.961 +/- 0.008 | 0.0 | 1
      19 | 0761670401 |  s: 1.005 +/- 0.006   d: 0.988 +/- 0.009 | 3.0 | 0
      20 | 0761670501 |  s: 1.012 +/- 0.007   d: 0.972 +/- 0.009 | 0.0 | 1
      20 | 0761670501 |  s: 1.004 +/- 0.007   d: 0.992 +/- 0.010 | 3.0 | 0
      21 | 0761670601 |  s: 1.017 +/- 0.006   d: 0.961 +/- 0.009 | 0.0 | 1
      21 | 0761670601 |  s: 1.009 +/- 0.007   d: 0.981 +/- 0.010 | 3.0 | 0
      22 | 0761670701 |  s: 1.010 +/- 0.006   d: 0.977 +/- 0.008 | 0.0 | 0
      23 | 0761670801 |  s: 1.019 +/- 0.006   d: 0.951 +/- 0.008 | 0.0 | 1
      23 | 0761670801 |  s: 1.004 +/- 0.006   d: 0.991 +/- 0.009 | 3.0 | 0
      24 | 0761670901 |  s: 1.006 +/- 0.006   d: 0.985 +/- 0.008 | 0.0 | 0


### Extract the Spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for iobs, o in enumerate(obsids):
    os.chdir('%s/pn/spec'%o)
    if len(glob.glob('spec*grp')) == 0:
        time.sleep(1)
        p = subp.Popen(['/bin/bash', '-i', '-c', 
           'sasinit; xmm_spec.py pn_filtered.fits ds9.reg -o spec_%d > spec.log 2>&1'%(iobs+1)])
        procs.append(p)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()
```

### Extract spectra from full region


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for iobs, o in enumerate(obsids):
    os.chdir('%s/pn'%o)
    os.system('mkdir -p spec0')
    os.chdir('spec0')
    
    if len(glob.glob('spec0*grp')) == 0:
        # use region from spec, but with no central region extracion #
        reg_lines = open('../spec/ds9.reg').readlines()
        if 'background' in reg_lines[-1]:
            reg_lines = reg_lines[:-2] + [reg_lines[-1], reg_lines[-2]]
        g = re.match('\\\\(.*\\\\)\\\\n', reg_lines[-1])
        dum_rad = reg_lines[-1].split(',')[2]
        reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',0,')
        with open('ds9.reg', 'w') as fp: fp.write(''.join(reg_lines))
        
    
        time.sleep(1)
        p = subp.Popen(['/bin/bash', '-i', '-c', 
           'sasinit; xmm_spec.py ../spec/pn_filtered.fits ds9.reg -o spec0_%d > spec.log 2>&1'%(iobs+1)])
        procs.append(p)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()
```

---
---
## Process MOS Data


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in obsids:
    os.chdir(o)
    os.system('mkdir -p mos')
    os.chdir('mos')
    if len(glob.glob('*EVL*')) == 0:
        log_file = '../../log/%s_process_mos.log'%o
        time.sleep(1)
        p = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process mos > %s 2>&1'%log_file])
        procs.append(p)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait()
```

## Spectral Extraction
### Standard filtering
Similar to pn case


```python
# wait until above is done! #
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
for o in obsids:
    print('-- obs %s --'%o)
    os.chdir('%s/mos'%o)
    os.system('ln -s %s m1.fits'%glob.glob('*M1*EVL*')[0])
    os.system('ln -s %s m2.fits'%glob.glob('*M2*EVL*')[0])
    os.system('mkdir -p spec')
    os.chdir('spec')
    for m in ['m1', 'm2']:
        # this line comes after first run; these 2 obs have now m1 data
        if m == 'm1' and o in ['0657840201', '0761670401']: continue
        if not exists('%s_filtered.fits'%m) or not exists('ds9_%s.reg'%m):
            # check if we have a saved region file, or a temporary region file
            # for faster loading
            saved_reg = '../../../log/%s_ds9_%s.reg'%(o, m)
            if exists(saved_reg):
                os.system('cp %s ds9_%s.reg'%(saved_reg, m))
                region = ''
            else:
                region = '--region'
                tmp_reg = '../../pn/spec/ds9.reg'
                if exists(tmp_reg):
                    os.system('cp %s tmp.reg'%tmp_reg)

            subp.call(['/bin/bash', '-i', '-c', 
                    'sasinit; xmm_filter.py ../%s.fits %s --std %s'%(m, m, region)]) 
    os.chdir('../../..')
```

    -- obs 0112310101 --
    -- obs 0112830201 --
    -- obs 0112830501 --
    -- obs 0143500101 --
    -- obs 0143500201 --
    -- obs 0143500301 --
    -- obs 0402660101 --
    -- obs 0402660201 --
    -- obs 0657840201 --
    -- obs 0657840301 --
    -- obs 0657840401 --
    -- obs 0679780101 --
    -- obs 0679780201 --
    -- obs 0679780301 --
    -- obs 0679780401 --
    -- obs 0761670101 --
    -- obs 0761670201 --
    -- obs 0761670301 --
    -- obs 0761670401 --
    -- obs 0761670501 --
    -- obs 0761670601 --
    -- obs 0761670701 --
    -- obs 0761670801 --
    -- obs 0761670901 --


```
No mos-1 data in 0657840201, 0761670401
```

### Check for pileup


```python
os.chdir('%s/%s'%(base_dir, data_dir))
print('obs_num | obsid | fractions | radius | pileup?')
for iobs, o in enumerate(obsids):
    os.chdir('%s/mos'%o)
    os.system('mkdir -p pileup')
    os.chdir('pileup')
    
    for m in ['m1', 'm2']:
        pileup = True
        radius = 0.0
        reg_file = '../spec/ds9_%s.reg'%m
        if not os.path.exists(reg_file): 
            print('no %s; skipping !'%reg_file)
            continue
        reg_lines = open(reg_file).readlines()
        if 'background' in reg_lines[-1]:
            reg_lines = reg_lines[:-2] + [reg_lines[-1], reg_lines[-2]]
        dum_rad = reg_lines[-1].split(',')[2]
        reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',%g,')
        reg_text = ''.join(reg_lines)

        while pileup:
            with open('ds9_%s.reg'%m, 'w') as fp: fp.write(reg_text%radius)
            subp.call(['/bin/bash', '-i', '-c', 
               ('sasinit; xmm_spec.py ../%s.fits ds9_%s.reg --check_pileup'
                ' > pileup_%s.log 2>&1')%(m,m,m)])
            line = [l for l in open('pileup_%s.log'%m).readlines() if '+/-' in l][0]
            frac = np.array(np.array(line.split())[[2,4,6,8]], np.double)
            pileup = (frac[0] > 1+3*frac[1]) or (frac[2] < 1-3*frac[3])
            # 1 arcsec = 20 pixels
            text = '{:4d}-{} | {} | {} | {} | {:d}'.format(
                iobs+1, m, o, line[10:-1], radius/20., pileup)
            print(text)
            if pileup: radius = np.max([3.0*20, radius+10])
            if radius/20 > 6: break
    os.chdir('../spec')
    os.system('cp ds9_%s.reg ds9_%s.reg.0'%(m,m))
    os.system('cp ../pileup/ds9_%s.reg ds9_%s.reg'%(m,m))
    os.chdir('../../../')
            
```

    obs_num | obsid | fractions | radius | pileup?
       1-m1 | 0112310101 |  s: 1.021 +/- 0.012   d: 0.944 +/- 0.020 | 0.0 | 0
       1-m2 | 0112310101 |  s: 1.032 +/- 0.012   d: 0.893 +/- 0.020 | 0.0 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.014   d: 0.891 +/- 0.022 | 3.0 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.014   d: 0.891 +/- 0.023 | 3.5 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.014   d: 0.893 +/- 0.023 | 4.0 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.015   d: 0.893 +/- 0.024 | 4.5 | 1
       1-m2 | 0112310101 |  s: 1.032 +/- 0.015   d: 0.898 +/- 0.025 | 5.0 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.016   d: 0.895 +/- 0.025 | 5.5 | 1
       1-m2 | 0112310101 |  s: 1.033 +/- 0.016   d: 0.896 +/- 0.026 | 6.0 | 1
       2-m1 | 0112830201 |  s: 1.021 +/- 0.009   d: 0.947 +/- 0.016 | 0.0 | 1
       2-m1 | 0112830201 |  s: 1.022 +/- 0.010   d: 0.942 +/- 0.017 | 3.0 | 1
       2-m1 | 0112830201 |  s: 1.022 +/- 0.010   d: 0.940 +/- 0.017 | 3.5 | 1
       2-m1 | 0112830201 |  s: 1.021 +/- 0.010   d: 0.942 +/- 0.018 | 4.0 | 1
       2-m1 | 0112830201 |  s: 1.021 +/- 0.011   d: 0.942 +/- 0.018 | 4.5 | 1
       2-m1 | 0112830201 |  s: 1.020 +/- 0.011   d: 0.946 +/- 0.018 | 5.0 | 0
       2-m2 | 0112830201 |  s: 1.040 +/- 0.010   d: 0.870 +/- 0.015 | 0.0 | 1
       2-m2 | 0112830201 |  s: 1.040 +/- 0.010   d: 0.864 +/- 0.016 | 3.0 | 1
       2-m2 | 0112830201 |  s: 1.041 +/- 0.010   d: 0.862 +/- 0.016 | 3.5 | 1
       2-m2 | 0112830201 |  s: 1.041 +/- 0.011   d: 0.864 +/- 0.017 | 4.0 | 1
       2-m2 | 0112830201 |  s: 1.041 +/- 0.011   d: 0.865 +/- 0.017 | 4.5 | 1
       2-m2 | 0112830201 |  s: 1.041 +/- 0.011   d: 0.865 +/- 0.018 | 5.0 | 1
       2-m2 | 0112830201 |  s: 1.041 +/- 0.011   d: 0.865 +/- 0.018 | 5.5 | 1
       2-m2 | 0112830201 |  s: 1.042 +/- 0.012   d: 0.862 +/- 0.018 | 6.0 | 1
       3-m1 | 0112830501 |  s: 1.024 +/- 0.016   d: 0.943 +/- 0.026 | 0.0 | 0
       3-m2 | 0112830501 |  s: 1.045 +/- 0.016   d: 0.848 +/- 0.025 | 0.0 | 1
       3-m2 | 0112830501 |  s: 1.046 +/- 0.017   d: 0.843 +/- 0.027 | 3.0 | 1
       3-m2 | 0112830501 |  s: 1.046 +/- 0.018   d: 0.845 +/- 0.027 | 3.5 | 1
       3-m2 | 0112830501 |  s: 1.044 +/- 0.018   d: 0.853 +/- 0.028 | 4.0 | 1
       3-m2 | 0112830501 |  s: 1.045 +/- 0.018   d: 0.849 +/- 0.029 | 4.5 | 1
       3-m2 | 0112830501 |  s: 1.043 +/- 0.019   d: 0.856 +/- 0.029 | 5.0 | 1
       3-m2 | 0112830501 |  s: 1.044 +/- 0.019   d: 0.852 +/- 0.030 | 5.5 | 1
       3-m2 | 0112830501 |  s: 1.044 +/- 0.020   d: 0.853 +/- 0.031 | 6.0 | 1
       4-m1 | 0143500101 |  s: 0.989 +/- 0.010   d: 1.063 +/- 0.019 | 0.0 | 0
       4-m2 | 0143500101 |  s: 1.010 +/- 0.011   d: 0.966 +/- 0.018 | 0.0 | 0
       5-m1 | 0143500201 |  s: 0.996 +/- 0.010   d: 1.018 +/- 0.018 | 0.0 | 0
       5-m2 | 0143500201 |  s: 1.008 +/- 0.010   d: 0.981 +/- 0.018 | 0.0 | 0
       6-m1 | 0143500301 |  s: 0.988 +/- 0.009   d: 1.054 +/- 0.016 | 0.0 | 0
       6-m2 | 0143500301 |  s: 1.007 +/- 0.009   d: 0.975 +/- 0.015 | 0.0 | 0
       7-m1 | 0402660101 |  s: 1.011 +/- 0.012   d: 0.994 +/- 0.020 | 0.0 | 0
       7-m2 | 0402660101 |  s: 1.035 +/- 0.011   d: 0.871 +/- 0.018 | 0.0 | 1
       7-m2 | 0402660101 |  s: 1.032 +/- 0.013   d: 0.884 +/- 0.021 | 3.0 | 1
       7-m2 | 0402660101 |  s: 1.032 +/- 0.013   d: 0.883 +/- 0.021 | 3.5 | 1
       7-m2 | 0402660101 |  s: 1.033 +/- 0.013   d: 0.880 +/- 0.022 | 4.0 | 1
       7-m2 | 0402660101 |  s: 1.032 +/- 0.014   d: 0.880 +/- 0.022 | 4.5 | 1
       7-m2 | 0402660101 |  s: 1.032 +/- 0.014   d: 0.878 +/- 0.023 | 5.0 | 1
       7-m2 | 0402660101 |  s: 1.031 +/- 0.014   d: 0.885 +/- 0.024 | 5.5 | 1
       7-m2 | 0402660101 |  s: 1.029 +/- 0.015   d: 0.895 +/- 0.024 | 6.0 | 1
       8-m1 | 0402660201 |  s: 1.011 +/- 0.011   d: 0.984 +/- 0.020 | 0.0 | 0
       8-m2 | 0402660201 |  s: 1.025 +/- 0.011   d: 0.913 +/- 0.019 | 0.0 | 1
       8-m2 | 0402660201 |  s: 1.026 +/- 0.012   d: 0.913 +/- 0.021 | 3.0 | 1
       8-m2 | 0402660201 |  s: 1.027 +/- 0.013   d: 0.911 +/- 0.021 | 3.5 | 1
       8-m2 | 0402660201 |  s: 1.025 +/- 0.013   d: 0.923 +/- 0.022 | 4.0 | 1
       8-m2 | 0402660201 |  s: 1.025 +/- 0.013   d: 0.924 +/- 0.022 | 4.5 | 1
       8-m2 | 0402660201 |  s: 1.025 +/- 0.014   d: 0.922 +/- 0.023 | 5.0 | 1
       8-m2 | 0402660201 |  s: 1.027 +/- 0.014   d: 0.918 +/- 0.024 | 5.5 | 1
       8-m2 | 0402660201 |  s: 1.027 +/- 0.015   d: 0.918 +/- 0.024 | 6.0 | 1
    no ../spec/ds9_m1.reg; skipping !
       9-m2 | 0657840201 |  s: 1.042 +/- 0.026   d: 0.848 +/- 0.040 | 0.0 | 1
       9-m2 | 0657840201 |  s: 1.042 +/- 0.027   d: 0.845 +/- 0.041 | 3.0 | 1
       9-m2 | 0657840201 |  s: 1.043 +/- 0.027   d: 0.845 +/- 0.042 | 3.5 | 1
       9-m2 | 0657840201 |  s: 1.043 +/- 0.027   d: 0.845 +/- 0.042 | 4.0 | 1
       9-m2 | 0657840201 |  s: 1.045 +/- 0.028   d: 0.840 +/- 0.043 | 4.5 | 1
       9-m2 | 0657840201 |  s: 1.042 +/- 0.028   d: 0.853 +/- 0.044 | 5.0 | 1
       9-m2 | 0657840201 |  s: 1.044 +/- 0.028   d: 0.843 +/- 0.044 | 5.5 | 1
       9-m2 | 0657840201 |  s: 1.045 +/- 0.029   d: 0.839 +/- 0.044 | 6.0 | 1
      10-m1 | 0657840301 |  s: 1.014 +/- 0.019   d: 0.956 +/- 0.033 | 0.0 | 0
      10-m2 | 0657840301 |  s: 1.031 +/- 0.019   d: 0.861 +/- 0.031 | 0.0 | 1
      10-m2 | 0657840301 |  s: 1.029 +/- 0.020   d: 0.873 +/- 0.032 | 3.0 | 1
      10-m2 | 0657840301 |  s: 1.027 +/- 0.020   d: 0.879 +/- 0.032 | 3.5 | 1
      10-m2 | 0657840301 |  s: 1.026 +/- 0.020   d: 0.886 +/- 0.033 | 4.0 | 1
      10-m2 | 0657840301 |  s: 1.024 +/- 0.020   d: 0.892 +/- 0.033 | 4.5 | 1
      10-m2 | 0657840301 |  s: 1.023 +/- 0.020   d: 0.895 +/- 0.034 | 5.0 | 1
      10-m2 | 0657840301 |  s: 1.023 +/- 0.021   d: 0.898 +/- 0.034 | 5.5 | 0
      11-m1 | 0657840401 |  s: 1.010 +/- 0.019   d: 0.973 +/- 0.033 | 0.0 | 0
      11-m2 | 0657840401 |  s: 1.022 +/- 0.019   d: 0.912 +/- 0.032 | 0.0 | 0
      12-m1 | 0679780101 |  s: 1.021 +/- 0.023   d: 0.937 +/- 0.038 | 0.0 | 0
      12-m2 | 0679780101 |  s: 1.036 +/- 0.023   d: 0.879 +/- 0.037 | 0.0 | 1
      12-m2 | 0679780101 |  s: 1.036 +/- 0.023   d: 0.877 +/- 0.037 | 3.0 | 1
      12-m2 | 0679780101 |  s: 1.034 +/- 0.023   d: 0.888 +/- 0.038 | 3.5 | 0
      13-m1 | 0679780201 |  s: 1.013 +/- 0.020   d: 0.993 +/- 0.035 | 0.0 | 0
      13-m2 | 0679780201 |  s: 1.029 +/- 0.020   d: 0.912 +/- 0.033 | 0.0 | 0
      14-m1 | 0679780301 |  s: 1.006 +/- 0.023   d: 1.018 +/- 0.041 | 0.0 | 0
      14-m2 | 0679780301 |  s: 1.021 +/- 0.023   d: 0.948 +/- 0.039 | 0.0 | 0
      15-m1 | 0679780401 |  s: 1.016 +/- 0.022   d: 0.954 +/- 0.038 | 0.0 | 0
      15-m2 | 0679780401 |  s: 1.035 +/- 0.023   d: 0.878 +/- 0.037 | 0.0 | 1
      15-m2 | 0679780401 |  s: 1.033 +/- 0.024   d: 0.892 +/- 0.039 | 3.0 | 0
      16-m1 | 0761670101 |  s: 1.015 +/- 0.012   d: 0.974 +/- 0.020 | 0.0 | 0
      16-m2 | 0761670101 |  s: 1.031 +/- 0.012   d: 0.900 +/- 0.020 | 0.0 | 1
      16-m2 | 0761670101 |  s: 1.031 +/- 0.013   d: 0.901 +/- 0.022 | 3.0 | 1
      16-m2 | 0761670101 |  s: 1.031 +/- 0.014   d: 0.897 +/- 0.022 | 3.5 | 1
      16-m2 | 0761670101 |  s: 1.031 +/- 0.014   d: 0.898 +/- 0.023 | 4.0 | 1
      16-m2 | 0761670101 |  s: 1.030 +/- 0.014   d: 0.904 +/- 0.024 | 4.5 | 1
      16-m2 | 0761670101 |  s: 1.029 +/- 0.015   d: 0.906 +/- 0.024 | 5.0 | 1
      16-m2 | 0761670101 |  s: 1.031 +/- 0.015   d: 0.900 +/- 0.025 | 5.5 | 1
      16-m2 | 0761670101 |  s: 1.032 +/- 0.016   d: 0.898 +/- 0.026 | 6.0 | 1
      17-m1 | 0761670201 |  s: 1.019 +/- 0.011   d: 0.952 +/- 0.019 | 0.0 | 0
      17-m2 | 0761670201 |  s: 1.028 +/- 0.011   d: 0.902 +/- 0.018 | 0.0 | 1
      17-m2 | 0761670201 |  s: 1.026 +/- 0.012   d: 0.915 +/- 0.020 | 3.0 | 1
      17-m2 | 0761670201 |  s: 1.026 +/- 0.012   d: 0.914 +/- 0.021 | 3.5 | 1
      17-m2 | 0761670201 |  s: 1.025 +/- 0.013   d: 0.920 +/- 0.021 | 4.0 | 1
      17-m2 | 0761670201 |  s: 1.026 +/- 0.013   d: 0.915 +/- 0.022 | 4.5 | 1
      17-m2 | 0761670201 |  s: 1.025 +/- 0.013   d: 0.917 +/- 0.022 | 5.0 | 1
      17-m2 | 0761670201 |  s: 1.025 +/- 0.014   d: 0.920 +/- 0.023 | 5.5 | 1
      17-m2 | 0761670201 |  s: 1.025 +/- 0.014   d: 0.922 +/- 0.023 | 6.0 | 1
      18-m1 | 0761670301 |  s: 1.017 +/- 0.010   d: 0.963 +/- 0.017 | 0.0 | 0
      18-m2 | 0761670301 |  s: 1.030 +/- 0.010   d: 0.892 +/- 0.017 | 0.0 | 1
      18-m2 | 0761670301 |  s: 1.031 +/- 0.011   d: 0.889 +/- 0.019 | 3.0 | 1
      18-m2 | 0761670301 |  s: 1.032 +/- 0.012   d: 0.886 +/- 0.019 | 3.5 | 1
      18-m2 | 0761670301 |  s: 1.031 +/- 0.012   d: 0.888 +/- 0.020 | 4.0 | 1
      18-m2 | 0761670301 |  s: 1.031 +/- 0.012   d: 0.890 +/- 0.020 | 4.5 | 1
      18-m2 | 0761670301 |  s: 1.030 +/- 0.013   d: 0.895 +/- 0.021 | 5.0 | 1
      18-m2 | 0761670301 |  s: 1.029 +/- 0.013   d: 0.897 +/- 0.021 | 5.5 | 1
      18-m2 | 0761670301 |  s: 1.028 +/- 0.013   d: 0.902 +/- 0.022 | 6.0 | 1
    no ../spec/ds9_m1.reg; skipping !
      19-m2 | 0761670401 |  s: 1.021 +/- 0.015   d: 0.936 +/- 0.026 | 0.0 | 0
      20-m1 | 0761670501 |  s: 1.006 +/- 0.011   d: 1.020 +/- 0.020 | 0.0 | 0
      20-m2 | 0761670501 |  s: 1.029 +/- 0.011   d: 0.902 +/- 0.018 | 0.0 | 1
      20-m2 | 0761670501 |  s: 1.030 +/- 0.012   d: 0.901 +/- 0.020 | 3.0 | 1
      20-m2 | 0761670501 |  s: 1.028 +/- 0.012   d: 0.907 +/- 0.021 | 3.5 | 1
      20-m2 | 0761670501 |  s: 1.027 +/- 0.013   d: 0.910 +/- 0.021 | 4.0 | 1
      20-m2 | 0761670501 |  s: 1.028 +/- 0.013   d: 0.907 +/- 0.022 | 4.5 | 1
      20-m2 | 0761670501 |  s: 1.030 +/- 0.013   d: 0.899 +/- 0.022 | 5.0 | 1
      20-m2 | 0761670501 |  s: 1.030 +/- 0.014   d: 0.899 +/- 0.023 | 5.5 | 1
      20-m2 | 0761670501 |  s: 1.031 +/- 0.014   d: 0.897 +/- 0.023 | 6.0 | 1
      21-m1 | 0761670601 |  s: 1.015 +/- 0.010   d: 0.977 +/- 0.018 | 0.0 | 0
      21-m2 | 0761670601 |  s: 1.031 +/- 0.011   d: 0.890 +/- 0.017 | 0.0 | 1
      21-m2 | 0761670601 |  s: 1.031 +/- 0.012   d: 0.889 +/- 0.019 | 3.0 | 1
      21-m2 | 0761670601 |  s: 1.030 +/- 0.012   d: 0.892 +/- 0.019 | 3.5 | 1
      21-m2 | 0761670601 |  s: 1.031 +/- 0.012   d: 0.887 +/- 0.020 | 4.0 | 1
      21-m2 | 0761670601 |  s: 1.032 +/- 0.012   d: 0.884 +/- 0.020 | 4.5 | 1
      21-m2 | 0761670601 |  s: 1.033 +/- 0.013   d: 0.881 +/- 0.021 | 5.0 | 1
      21-m2 | 0761670601 |  s: 1.033 +/- 0.013   d: 0.882 +/- 0.021 | 5.5 | 1
      21-m2 | 0761670601 |  s: 1.032 +/- 0.013   d: 0.883 +/- 0.022 | 6.0 | 1
      22-m1 | 0761670701 |  s: 1.007 +/- 0.010   d: 1.012 +/- 0.017 | 0.0 | 0
      22-m2 | 0761670701 |  s: 1.029 +/- 0.010   d: 0.900 +/- 0.016 | 0.0 | 1
      22-m2 | 0761670701 |  s: 1.028 +/- 0.010   d: 0.903 +/- 0.017 | 3.0 | 1
      22-m2 | 0761670701 |  s: 1.026 +/- 0.011   d: 0.911 +/- 0.018 | 3.5 | 1
      22-m2 | 0761670701 |  s: 1.027 +/- 0.011   d: 0.909 +/- 0.018 | 4.0 | 1
      22-m2 | 0761670701 |  s: 1.025 +/- 0.011   d: 0.916 +/- 0.019 | 4.5 | 1
      22-m2 | 0761670701 |  s: 1.025 +/- 0.011   d: 0.914 +/- 0.019 | 5.0 | 1
      22-m2 | 0761670701 |  s: 1.025 +/- 0.012   d: 0.913 +/- 0.019 | 5.5 | 1
      22-m2 | 0761670701 |  s: 1.026 +/- 0.012   d: 0.912 +/- 0.020 | 6.0 | 1
      23-m1 | 0761670801 |  s: 1.011 +/- 0.009   d: 0.982 +/- 0.016 | 0.0 | 0
      23-m2 | 0761670801 |  s: 1.026 +/- 0.010   d: 0.916 +/- 0.016 | 0.0 | 1
      23-m2 | 0761670801 |  s: 1.026 +/- 0.011   d: 0.918 +/- 0.018 | 3.0 | 1
      23-m2 | 0761670801 |  s: 1.026 +/- 0.011   d: 0.918 +/- 0.018 | 3.5 | 1
      23-m2 | 0761670801 |  s: 1.025 +/- 0.011   d: 0.922 +/- 0.019 | 4.0 | 1
      23-m2 | 0761670801 |  s: 1.024 +/- 0.011   d: 0.926 +/- 0.019 | 4.5 | 1
      23-m2 | 0761670801 |  s: 1.024 +/- 0.012   d: 0.924 +/- 0.020 | 5.0 | 1
      23-m2 | 0761670801 |  s: 1.025 +/- 0.012   d: 0.923 +/- 0.020 | 5.5 | 1
      23-m2 | 0761670801 |  s: 1.024 +/- 0.012   d: 0.924 +/- 0.021 | 6.0 | 1
      24-m1 | 0761670901 |  s: 1.009 +/- 0.010   d: 1.002 +/- 0.018 | 0.0 | 0
      24-m2 | 0761670901 |  s: 1.029 +/- 0.010   d: 0.900 +/- 0.016 | 0.0 | 1
      24-m2 | 0761670901 |  s: 1.030 +/- 0.010   d: 0.895 +/- 0.017 | 3.0 | 1
      24-m2 | 0761670901 |  s: 1.029 +/- 0.011   d: 0.898 +/- 0.017 | 3.5 | 1
      24-m2 | 0761670901 |  s: 1.028 +/- 0.011   d: 0.902 +/- 0.018 | 4.0 | 1
      24-m2 | 0761670901 |  s: 1.028 +/- 0.011   d: 0.904 +/- 0.018 | 4.5 | 1
      24-m2 | 0761670901 |  s: 1.026 +/- 0.011   d: 0.911 +/- 0.019 | 5.0 | 1
      24-m2 | 0761670901 |  s: 1.026 +/- 0.012   d: 0.912 +/- 0.020 | 5.5 | 1
      24-m2 | 0761670901 |  s: 1.026 +/- 0.012   d: 0.915 +/- 0.020 | 6.0 | 1


### Extract Spectra
this may need to be run a few times in case python cannot handle openning so many processes at once. Running the code again will skip files already created.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
print('mos 1')
for iobs, o in enumerate(obsids):
    os.chdir('%s/mos/spec'%o)
    if os.path.exists('ds9_m1.reg') and len(glob.glob('spec_m1*grp'))==0:
        p = subp.Popen(['/bin/bash', '-i', '-c', 
           ('sasinit; xmm_spec.py m1_filtered.fits ds9_m1.reg -o spec_m1_%d '
            '> spec_m1.log 2>&1')%(iobs+1)])
        procs.append(p)
        time.sleep(1)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()

procs = []
print('mos 2')
for iobs, o in enumerate(obsids):
    os.chdir('%s/mos/spec'%o)
    if os.path.exists('ds9_m2.reg') and len(glob.glob('spec_m2*grp'))==0:
        p = subp.Popen(['/bin/bash', '-i', '-c', 
           ('sasinit; xmm_spec.py m2_filtered.fits ds9_m2.reg -o spec_m2_%d '
            '> spec_m2.log 2>&1')%(iobs+1)])
        procs.append(p)
        time.sleep(1)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()
```

    mos 1
    mos 2


---
---
## Process RGS Data


```python
os.chdir('%s/%s'%(base_dir, data_dir))
procs = []
for o in obsids:
    os.chdir(o)
    os.system('mkdir -p rgs')
    os.chdir('rgs')
    if len(glob.glob('spec_rgs*')) == 0:
        log_file = '../../log/%s_process_rgs.log'%o
        p = subp.Popen(['/bin/bash', '-i', '-c', 'sasinit; xmm_process rgs > %s 2>&1'%log_file])
        time.sleep(1)
        procs.append(p)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait()
```

#### Rename the RGS spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s/rgs'%o)
    if len(glob.glob('spec_rgs_%d*'%iobs)) == 0:
        os.system('rename _rgs. _rgs_%d. spec_rgs*'%(iobs+1))
        root = 'spec_rgs_%d.'%(iobs+1) + '%s'
        with pyfits.open(root%'grp') as fp:
            fp[1].header['backfile'] = root%'bgd'
            fp[1].header['respfile'] = root%'rsp'
            os.system('rm tmp.grp > /dev/null 2>&1')
            fp.writeto('tmp.grp')
        os.system('mv %s _%s'%(root%'grp', root%'grp'))
        os.system('mv %s tmp.grp'%(root%'grp'))
        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
              'ftgrouppha tmp.grp %s opt respfile=%s')%(root%'grp', root%'rsp')
        subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')

```

## Summary of PN Spectral Data


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
print('{:5} | {:12} | {:10.8} | {:10.3} | {:10.3}'.format(
        'num', 'obsid', 'mjd', 'rate', 'exposure'))
spec_data = []
for iobs,o in enumerate(obsids):
    with pyfits.open('%s/pn/spec/spec_%d.grp'%(o, iobs+1)) as fp:
        exposure = fp[1].header['exposure']
        counts = fp[1].data.field('counts').sum()
        mjd = (atime.Time(fp[0].header['date_end']).mjd + 
               atime.Time(fp[0].header['date_obs']).mjd ) / 2
        spec_data.append([mjd, counts/exposure, exposure/1e3])
        text = '{:5} | {:12} | {:10.8} | {:10.3} | {:10.3}'.format(
                iobs+1, o, mjd, counts/exposure, exposure/1e3)
        print(text)
spec_data = np.array(spec_data)
```

    num   | obsid        | mjd        | rat        | exp       
        1 | 0112310101   |  51899.889 |       7.07 |       21.0
        2 | 0112830201   |  51900.806 |       6.82 |       50.9
        3 | 0112830501   |  51900.255 |       6.83 |       17.6
        4 | 0143500101   |  52784.177 |       24.3 |       11.1
        5 | 0143500201   |  52785.967 |       20.8 |       12.7
        6 | 0143500301   |  52786.746 |       31.5 |       12.7
        7 | 0402660101   |  53871.499 |        6.4 |       28.0
        8 | 0402660201   |  54069.028 |       7.53 |       22.9
        9 | 0657840201   |  55724.662 |       16.1 |       2.21
       10 | 0657840301   |  55890.235 |       18.3 |       5.61
       11 | 0657840401   |    55904.2 |       21.2 |       6.59
       12 | 0679780101   |   56060.25 |       19.4 |       6.29
       13 | 0679780201   |  56088.184 |       8.37 |       8.75
       14 | 0679780301   |  56245.789 |       16.1 |       3.83
       15 | 0679780401   |  56271.736 |       18.9 |       6.64
       16 | 0761670101   |  57338.846 |       6.38 |       22.8
       17 | 0761670201   |  57340.826 |       7.11 |       24.0
       18 | 0761670301   |  57342.825 |       6.25 |       30.6
       19 | 0761670401   |    57356.8 |       6.67 |       21.9
       20 | 0761670501   |  57346.783 |       8.01 |       26.9
       21 | 0761670601   |  57348.801 |       8.32 |       29.9
       22 | 0761670701   |  57372.711 |       12.3 |       29.4
       23 | 0761670801   |  57374.679 |       10.9 |       29.7
       24 | 0761670901   |  57378.669 |       11.7 |       30.4



```python
## keep only exposures > 5ks
igood = np.argwhere(spec_data[:,2] >= 5)[:,0]
spec_obsids = obsids[igood]
spec_data = spec_data[igood]
print('There are %d spec observations'%len(spec_obsids))
print(', '.join(spec_obsids))
```

    There are 22 spec observations
    0112310101, 0112830201, 0112830501, 0143500101, 0143500201, 0143500301, 0402660101, 0402660201, 0657840301, 0657840401, 0679780101, 0679780201, 0679780401, 0761670101, 0761670201, 0761670301, 0761670401, 0761670501, 0761670601, 0761670701, 0761670801, 0761670901


## Extract PN At Small Time steps
### Identify the time selection


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
# time_step in ks
time_step = 5
tselect, tselect_expr = [], []
print('{:5} | {:12} | {:10} | {:5}'.format('nu', 'obsid', 'exposure', 'nsub spec'))
for iobs,o in enumerate(spec_obsids):
    nsub_spec = np.int(spec_data[iobs, 2] // time_step)
    with pyfits.open('%s/pn/spec/pn_filtered.fits'%o) as fp:
        events = np.sort(fp[1].data.field('time'))
        evt_per_sub = len(events)//nsub_spec
        tcut = np.concatenate([[events[0]], events[np.arange(1, nsub_spec)*evt_per_sub], [events[-1]]])
        tselect_expr.append(['(TIME > %.10g) && (TIME <= %.10g)'%(x,y) 
                        for x,y in zip(tcut[:-1], tcut[1:])])
        tselect.append([[x,y] for x,y in zip(tcut[:-1], tcut[1:])])
    print('{:5} | {:12} | {:10.3} | {:5}'.format(iobs+1, o, spec_data[iobs, 2], nsub_spec))
```

    nu    | obsid        | exposure   | nsub spec
        1 | 0112310101   |       21.0 |     4
        2 | 0112830201   |       50.9 |    10
        3 | 0112830501   |       17.6 |     3
        4 | 0143500101   |       11.1 |     2
        5 | 0143500201   |       12.7 |     2
        6 | 0143500301   |       12.7 |     2
        7 | 0402660101   |       28.0 |     5
        8 | 0402660201   |       22.9 |     4
        9 | 0657840301   |       5.61 |     1
       10 | 0657840401   |       6.59 |     1
       11 | 0679780101   |       6.29 |     1
       12 | 0679780201   |       8.75 |     1
       13 | 0679780401   |       6.64 |     1
       14 | 0761670101   |       22.8 |     4
       15 | 0761670201   |       24.0 |     4
       16 | 0761670301   |       30.6 |     6
       17 | 0761670401   |       21.9 |     4
       18 | 0761670501   |       26.9 |     5
       19 | 0761670601   |       29.9 |     5
       20 | 0761670701   |       29.4 |     5
       21 | 0761670801   |       29.7 |     5
       22 | 0761670901   |       30.4 |     6


### Now extract the spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
# use tselect array from above #
procs, ispec, ispecs = [], 1, []
for iobs in range(len(spec_obsids)):
    dum = [ispec+i for i in range(len(tselect[iobs]))]
    ispecs.append(dum)
    ispec = dum[-1]+1
        
for iobs,o in enumerate(spec_obsids):
    os.chdir('%s/pn'%o)
    os.system('mkdir -p subspec')
    os.chdir('subspec')

    if len(glob.glob('spec*grp')) != len(tselect_expr[iobs]):
        cmd = 'sasinit'
        for isel,tsel in enumerate(tselect_expr[iobs]):
            ispec = ispecs[iobs][isel]
            cmd += (';xmm_spec.py ../spec/pn_filtered.fits ../spec/ds9.reg'
                ' --e_expr " && %s" -o spec_%d > spec_%d.log 2>&1')%(tsel, ispec, ispec)
        time.sleep(0.5)
        p = subp.Popen(['/bin/bash', '-i', '-c', cmd])    
        procs.append(p)
    os.chdir('../../..')
# wait for the tasks to end
for p in procs: p.wait()
```

### Save useful data for other notebooks


```python
os.chdir('%s/%s'%(base_dir, data_dir))
# save some useful data for other notebooks
np.savez('log/data.npz', obsids=obsids, spec_obsids=spec_obsids, spec_data=spec_data, 
         tselect=tselect, tselect_ispec=ispecs)
```

<br /> <br /> <br />

---
---
## SUZAKU data


```python
base_dir = '/u/home/abzoghbi/data/ngc4151/spec_analysis'
data_dir = 'data/suzaku'
os.system('mkdir -p %s'%data_dir)
obsids = ['906006010', '906006020', '707024010', '701034010']
obsids = np.sort(obsids)
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 4 observations
    701034010, 707024010, 906006010, 906006020


#### We use `ftplib` to get the data from heasarc (may take some time)


```python
os.chdir(base_dir)
ftp = FTP('legacy.gsfc.nasa.gov', 'anonymous', 'anonymous@gmail.com')
ftp.cwd('suzaku/data/obs')
failed = []
for o in obsids:
    tar_file = '%s/%s.tar'%(data_dir, o)
    ftp.cwd(o[0])
    # download file only if not already downloaded
    if not os.path.exists(tar_file):
        try:
            ftp.retrbinary('RETR %s.tar'%o ,open(tar_file, 'wb').write)
        except:
            print('failed downloading %s'%o)
            os.system('rm %s >/dev/null 2>&1'%tar_file)
            failed.append(o)
    ftp.cwd('..')
```


```python
for f in failed:
    obsids = np.delete(obsids, np.argwhere(obsids==f)[0,0])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 4 observations
    906006010, 906006020, 707024010, 707034010


## Process the XIS data
We use our shell script `suzaku_process`.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
os.system('mkdir -p log')
procs = []
for o in obsids:
    if os.path.exists(o): continue
    os.system('tar -xf %s.tar'%o)
    if not os.path.exists('%s_p'%o):
        log_file = 'log/%s_process.log'%o
        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
               'suzaku_process %s xis > %s 2>&1'%(o, log_file))
        proc = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(proc)
        time.sleep(0.2)

# wait for the tasks to end
for p in procs: p.wait()
```

## Spectral Extraction
- Use `suzaku_xis_spec.py` 
- Region size: 160''


```python
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
procs = []
for iobs,o in enumerate(obsids):
    print('-- obs %s --'%o)
    os.chdir('%s_p/'%o)
    os.system('mkdir -p spec')
    os.chdir('spec')
    if len(glob.glob('spec_xi*grp')) != 3:
        # check if we have a saved region file, or a temporary region file
        # for faster loading
        saved_reg = '../../log/%s_src.reg'%o
        if exists(saved_reg):
            os.system('cp %s src.reg'%saved_reg)
            os.system('cp %s bgd.reg'%(saved_reg.replace('_src.', '_bgd.')))
            region = ''
        else:
            region = '--create_region'

        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
              'suzaku_xis_spec.py -o spec_%d %s'%(iobs+1, region))
        p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(p)
        if not exists(saved_reg):
            os.system('cp src.reg %s'%saved_reg) 
            os.system('cp bgd.reg %s'%(saved_reg.replace('_src.', '_bgd.')))
        time.sleep(0.3)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait() 
```

    -- obs 701034010 --
    -- obs 707024010 --
    -- obs 906006010 --
    -- obs 906006020 --



```python
## add xi0, xi3 -> fi #
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s_p/spec'%o)
    txt = '\n'.join(['spec_xi0_{0}.{1} spec_xi3_{0}.{1}'.format(iobs+1, s) 
                        for s in ['src', 'bgd', 'rsp']])
    with open('tmp.dat', 'w') as fp: fp.write(txt+'\n')
    cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null; rm spec_fi_{0}.*;'
            'addascaspec tmp.dat spec_fi_{0}.src spec_fi_{0}.rsp spec_fi_{0}.bgd').format(iobs+1)
    p = subp.call(['/bin/bash', '-i', '-c', cmd])
    # group #
    cmd = 'ogrppha.py spec_fi_{0}.src spec_fi_{0}.grp -f 3 -s 6'.format(iobs+1)
    subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')
```

## Summary of Spectral Data


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
print('{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.3}'.format(
        'num', 'obsid', 'mjd_s', 'mjd_e', 'rate', 'exposure'))
spec_data = []
for iobs,o in enumerate(obsids):
    with pyfits.open('%s_p/spec/spec_xi0_%d.grp'%(o, iobs+1)) as fp:
        exposure = fp[1].header['exposure']
        counts = fp[1].data.field('counts').sum()
        tmid = np.array([fp[0].header['tstart'],  fp[0].header['tstop']])
        mref = fp[0].header['mjdrefi'] + fp[0].header['mjdreff']
        mjd = tmid / (24*3600) + mref
        spec_data.append([mjd, counts/exposure, exposure/1e3])
        text = '{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.5}'.format(
                iobs+1, o, mjd[0], mjd[1], counts/exposure, exposure/1e3)
        print(text)
spec_data = np.array(spec_data)
```

    num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
        1 | 701034010    |   54087.84 |  54090.385 |       1.11 |     124.98
        2 | 707024010    |  56242.862 |  56245.961 |       3.63 |     150.28
        3 | 906006010    |  55882.674 |  55883.975 |        4.8 |     61.665
        4 | 906006020    |   55913.67 |  55915.041 |       5.73 |     60.596


## Extract XIS At Small Time steps
### Identify the time selection


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
cmd_clean = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
# time_step in ks
time_step = 5
tselect, tselect_expr = [], []
print('{:5} | {:12} | {:10} | {:5}'.format('num', 'obsid', 'exposure', 'nsub spec'))
for iobs,o in enumerate(obsids):
    
    os.chdir('%s_p'%o)
    os.system('mkdir -p subspec')
    os.chdir('subspec')
    os.system('cp ../spec/*reg .')
    
    # extract light curve #
    lcfile = 'lc_008_xi0__1.src'
    if not os.path.exists(lcfile):
        cmd = cmd_clean + 'suzaku_xis_lc.py -t 8 --rootdir ../xis/event_cl'
        p = subp.call(['/bin/bash', '-i', '-c', cmd])

    ldata,dt = az.LCurve.read_fits_file(lcfile, min_exp=0.7)
    ldata = ldata[:, ldata[2]>0]
    nsub_spec = ldata.shape[1]//np.int(time_step * 1e3/dt)
    icut = np.array(np.linspace(0, ldata.shape[1]-1, nsub_spec+1), np.int)
    tcut = ldata[0, icut]
    with pyfits.open(lcfile) as fp:
        mjdref = fp[0].header['mjdrefi'] + fp[0].header['mjdreff']
    tcut = tcut/(24*3600) + mjdref
    
    tselect_expr.append(['mjd \\"%10.10g,%10.10g\\"'%(x,y) 
                        for x,y in zip(tcut[:-1], tcut[1:])])
    tselect.append([[x,y] for x,y in zip(tcut[:-1], tcut[1:])])
    print('{:5} | {:12} | {:10.5} | {:5}'.format(iobs+1, o, spec_data[iobs, 2], nsub_spec))
    os.chdir('../..')
```

    num   | obsid        | exposure   | nsub spec
        1 | 701034010    |     124.98 |    24
        2 | 707024010    |     150.28 |    29
        3 | 906006010    |     61.665 |    12
        4 | 906006020    |     60.596 |    12


### Now extract the spectra


```python
os.chdir('%s/%s'%(base_dir, data_dir))
# use tselect array from above #
procs, ispec, ispecs = [], 1, []
for iobs in range(len(obsids)):
    dum = [ispec+i for i in range(len(tselect[iobs]))]
    ispecs.append(dum)
    ispec = dum[-1]+1
        
for iobs,o in enumerate(obsids):
    os.chdir('%s_p/subspec'%o)

    
    cmd0 = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
    for isel,tsel in enumerate(tselect_expr[iobs]):
        ispec = ispecs[iobs][isel]
        logfile = '../../log/sub_%d.log'%ispec
        cmd = cmd0 + 'suzaku_xis_spec.py -o spec_%d --t_expr "%s" --noclean > %s 2>&1'%(
                    ispec, tsel, logfile)

        if len(glob.glob('spec_xi*_%d.grp'%ispec)) != 3:
            time.sleep(1)
            p = subp.Popen(['/bin/bash', '-i', '-c', cmd])    
            procs.append(p)

        if len(procs) == 30:
            for p in procs: p.wait()
            procs = []
    
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait()
```


```python
## add xi0, xi3 -> fi #
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s_p/subspec'%o)
    cmd0 = 'export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
    for isel,tsel in enumerate(tselect_expr[iobs]):
        ispec = ispecs[iobs][isel]
        if not os.path.exists('spec_fi_{0}.grp'.format(ispec)):
            txt = '\n'.join(['spec_xi0_{0}.{1} spec_xi3_{0}.{1}'.format(ispec, s) 
                            for s in ['src', 'bgd', 'rsp']])
            with open('tmp.dat', 'w') as fp: fp.write(txt)
            cmd = 'addascaspec tmp.dat spec_fi_{0}.src spec_fi_{0}.rsp spec_fi_{0}.bgd'.format(ispec)
            p = subp.call(['/bin/bash', '-i', '-c', cmd0+cmd])
            # group #
            cmd = 'ogrppha.py spec_fi_{0}.src spec_fi_{0}.grp -f 3 -s 6'.format(ispec)
            subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')
```


```python
os.chdir('%s/%s'%(base_dir, data_dir))
# save some useful data for other notebooks
np.savez('log/data.npz', obsids=obsids, spec_data=spec_data, 
         tselect=tselect, tselect_ispec=ispecs)
```

<br /> <br /> <br />

---
---
## NuSTAR data


```python
base_dir = '/u/home/abzoghbi/data/ngc4151/spec_analysis'
data_dir = 'data/nustar'
os.system('mkdir -p %s'%data_dir)
obsids = ['60001111002', '60001111003', '60001111005']
obsids = np.array(obsids)
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
```

    There are 3 observations
    60001111002, 60001111003, 60001111005


#### We use `ftplib` to get the data from heasarc (may take some time)


```python
os.chdir('%s/%s'%(base_dir, data_dir))
ftp = FTP('legacy.gsfc.nasa.gov', 'anonymous', 'anonymous@gmail.com')
ftp.cwd('nustar/data/obs')
failed = []
for o in obsids:
    ftp.cwd('%s/%s'%(o[1:3], o[0]))
    tar_file = '%s.tar'%o
    # download file only if not already downloaded
    if not os.path.exists(tar_file):
        try:
            ftp.retrbinary('RETR %s'%tar_file ,open(tar_file, 'wb').write)
        except:
            print('failed downloading %s'%o)
            os.system('rm %s >/dev/null 2>&1'%tar_file)
            failed.append(o)
    ftp.cwd('../..')

```


```python
for f in failed:
    if f in obsids:
        obsids = np.delete(obsids, np.argwhere(obsids==f)[0,0])
print('There are %d observations'%len(obsids))
print(', '.join(obsids))
nobs = len(obsids)
```

    There are 3 observations
    60001111002, 60001111003, 60001111005


### Process the NuSTAR data
We use our shell script `nustar_process`.


```python
os.chdir('%s/%s'%(base_dir, data_dir))
os.system('mkdir -p log')
procs = []
for o in obsids:
    if not os.path.exists(o):
        os.system('tar -xf %s.tar'%o)
    if not os.path.exists('%s_p'%o):
        # download large files by http
        
        log_file = 'log/%s_process.log'%o
        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
               'nustar_process %s > %s 2>&1'%(o, log_file))
        proc = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(proc)
        time.sleep(0.2)

# wait for the tasks to end
for p in procs: p.wait()
```

### Spectral Extraction
- Use `nustar_spec.py` 
- Region size: 150''


```python
os.chdir('%s/%s'%(base_dir, data_dir))
exists = os.path.exists
obsids = np.sort(obsids)
procs = []
for iobs,o in enumerate(obsids):
    print('-- obs %s --'%o)
    os.chdir('%s_p/'%o)
    os.system('mkdir -p spec')
    os.chdir('spec')
    if len(glob.glob('spec*grp')) != 2:
        # check if we have a saved region file, or a temporary region file
        # for faster loading
        saved_reg = '../../log/%s_src.reg'%o
        if exists(saved_reg):
            os.system('cp %s src.reg'%saved_reg)
            os.system('cp %s bgd.reg'%(saved_reg.replace('_src.', '_bgd.')))
            region = ''
        else:
            region = '--create_region'

        cmd = ('export HEADASNOQUERY=;export HEADASPROMPT=/dev/null;'
              'nustar_spec.py -o spec_%d %s'%(iobs+1, region))
        p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
        procs.append(p)
        if not exists(saved_reg):
            os.system('cp src.reg %s'%saved_reg) 
            os.system('cp bgd.reg %s'%(saved_reg.replace('_src.', '_bgd.')))
        time.sleep(0.3)
    os.chdir('../..')
# wait for the tasks to end
for p in procs: p.wait() 
```

    -- obs 60001111002 --
    -- obs 60001111003 --
    -- obs 60001111005 --



```python
## group the spectra #
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs,o in enumerate(obsids):
    os.chdir('%s_p/spec'%o)
    cmd = ('rm *grp; ogrppha.py spec_{0}_a_sr.pha spec_{0}_a.grp -f 3 -s 6;'
           'ogrppha.py spec_{0}_b_sr.pha spec_{0}_b.grp -f 3 -s 6').format(iobs+1)
    subp.call(['/bin/bash', '-i', '-c', cmd])
    os.chdir('../..')
```

### Summary of spectral data


```python
os.chdir('%s/%s/'%(base_dir, data_dir))
print('{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.3}'.format(
        'num', 'obsid', 'mjd_s', 'mjd_e', 'rate', 'exposure'))
spec_data = []
for iobs,o in enumerate(obsids):
    with pyfits.open('%s_p/spec/spec_%d_a.grp'%(o, iobs+1)) as fp:
        exposure = fp[1].header['exposure']
        counts = fp[1].data.field('counts').sum()
        tmid = np.array([fp[0].header['tstart'], fp[0].header['tstop']])
        mref = fp[0].header['mjdrefi'] + fp[0].header['mjdreff']
        mjd = tmid / (24*3600) + mref
        spec_data.append([mjd, counts/exposure, exposure/1e3])
        text = '{:5} | {:12} | {:10.8} | {:10.8} | {:10.3} | {:10.5}'.format(
                iobs+1, o, mjd[0], mjd[1], counts/exposure, exposure/1e3)
        print(text)
spec_data = np.array(spec_data)
```

    num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
        1 | 60001111002  |  56243.264 |  56243.765 |       7.64 |     21.864
        2 | 60001111003  |  56243.792 |  56245.045 |       7.17 |     57.036
        3 | 60001111005  |  56245.345 |   56246.73 |       8.26 |     61.531



```python
# summary #
os.chdir('%s/%s'%(base_dir, data_dir))
# save some useful data for other notebooks
np.savez('log/data.npz', obsids=obsids, spec_data=spec_data)
```

### Overlap between NuSTAR & Suzaku
Use `suzaku_2` with `nustar_2`


```python
# ---- SUZAKU ---- #
# num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
#     1 | 701034010    |   54087.84 |  54090.385 |       1.11 |     124.98
#     2 | 707024010    |  56242.862 |  56245.961 |       3.63 |     150.28
#     3 | 906006010    |  55882.674 |  55883.975 |        4.8 |     61.665
#     4 | 906006020    |   55913.67 |  55915.041 |       5.73 |     60.596

# ---- NUSTAR ---- #
# num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
#     1 | 60001111002  |  56243.264 |  56243.765 |       7.64 |     21.864
#     2 | 60001111003  |  56243.792 |  56245.045 |       7.17 |     57.036
#     3 | 60001111005  |  56245.345 |   56246.73 |       8.26 |     61.531

```
