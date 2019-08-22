
### Description
- Extract light curves for timing analysis


```python
import sys,os
base_dir = '/home/abzoghbi/data/ngc4151/spec_analysis'
sys.path.append(base_dir)
from spec_helpers import *
%load_ext autoreload
%autoreload 2
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload



```python
### Read useful data from data notebook
data_dir = 'data/xmm'
spec_dir = 'data/xmm_spec'
os.chdir('%s/%s'%(base_dir, data_dir))
data = np.load('log/data.npz')
spec_obsids = data['spec_obsids']
obsids = data['obsids']
spec_data = data['spec_data']
spec_ids = [i+1 for i,o in enumerate(obsids) if o in spec_obsids]
```

---
## Light curve Extraction



```python
# prepare for light curve extraction #
os.chdir('%s/%s'%(base_dir, data_dir))
for iobs, o in enumerate(obsids):
    os.chdir('%s/pn'%o)
    if os.path.exists('lc'):
        os.chdir('../..') 
        continue
    os.system('mkdir -p lc')
    os.chdir('lc')
    os.system('ln -s ../../odf/ccf.cif ccf.cif >/dev/null 2>&1')
    if not os.path.exists('pn_filtered.fits'):
        # use region from spec, but with no central region extracion #
        reg_lines = open('../spec/ds9.reg').readlines()
        g = re.match('\\(.*\\)\\n', reg_lines[-1])
        dum_rad = reg_lines[-1].split(',')[2]
        reg_lines[-1] = reg_lines[-1].replace(',' + dum_rad + ',', ',0,')
        reg_text = ''.join(reg_lines)
        with open('ds9.reg', 'w') as fp: fp.write(reg_text)
        
        p = subp.call(['/bin/bash', '-i', '-c', 
                       'sasinit; xmm_filter.py ../pn.fits pn --std --stdR 0.5'])
    os.chdir('../../..')
```

### Extract function


```python
def _extract_lc(lcdir, ebins, dt):
    os.chdir('%s/%s'%(base_dir, data_dir))
    procs = []
    for o in obsids:
        os.system('mkdir -p %s/pn/lc/%s'%(o, lcdir))
        os.chdir('%s/pn/lc/%s'%(o, lcdir))
        if len(glob.glob('lc_{:03g}*.fits'.format(dt))) != 3*(len(ebins.split())-1): 
            cmd = ('sasinit;xmm_lc.py ../pn_filtered.fits ../ds9.reg'
                   ' -e "%s" -t %g >lc.log 2>&1')%(ebins, dt)
            time.sleep(0.5)
            p = subp.Popen(['/bin/bash', '-i', '-c', cmd])
            procs.append(p)
        os.chdir('../../../..')
    # wait for the tasks to end
    for p in procs: p.wait()
```

### 1 bin in 2-10 keV: `1b`


```python
lcdir, ebins, dt = '1b', '2 10', 128
_extract_lc(lcdir, ebins, dt)
```

### 8 bins: `8l -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 9)])`


```python
lcdir, ebins, dt = '8l', '2 2.4 3 3.7 4.5 5.5 6.7 8.2 10', 128
_extract_lc(lcdir, ebins, dt)
```

### 16 bins: `16l -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 17)])`


```python
lcdir, ebins, dt = '16l', ('2 2.2 2.4 2.7 3 3.3 3.7 4 4.5 4.9 5.5 6 6.7 7.4 8.2 9 10'), 128
_extract_lc(lcdir, ebins, dt)
```

### 24 bins: `24l -> ' '.join(['{:2.2g}'.format(x) for x in logspace(log10(2), log10(10), 25)])`


```python
lcdir, ebins, dt = '24l', ('2 2.33 2.67 3 3.33 3.67 4 4.33 4.67 5 5.33 5.67 6 6.33 6.67 '
                           '7 7.33 7.67 8 8.33 8.67  9 9.33 9.67 10'), 128
_extract_lc(lcdir, ebins, dt)
```

### 32 bins: `32l -> ' '.join(['{:2.3g}'.format(x) for x in logspace(log10(2), log10(10), 33)])`


```python
lcdir, ebins, dt = '32l', ('2 2.1 2.21 2.33 2.45 2.57 2.7 2.84 2.99 3.14 3.31 3.48 '
        '3.66 3.85 4.04 4.25 4.47 4.7 4.95 5.2 5.47 5.75 6.05 6.36 6.69 7.03 7.4 7.78 '
        '8.18 8.6 9.04 9.51 10'), 256
_extract_lc(lcdir, ebins, dt)
```


```python

```
