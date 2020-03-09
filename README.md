## Revisiting The Spectral and Timing Properties of NGC 4151

These papes contain the data and codes associated with the article entitled: **[Revisiting The Spectral and Timing Properties of NGC 4151](https://arxiv.org/abs/1908.09862)**, published in the Astrophysical Journal.

### Abstract
> NGC 4151 is the brightest Seyfert 1 nucleus in X-rays. It was the first object to show short time delays in the Fe K band, which were attributed to relativistic reverberation, providing a new tool for probing regions at the black hole scale. Here, we report the results of a large XMM-Newton campaign in 2015 to study these short delays further. Analyzing high quality data that span time scales between hours and decades, we find that neutral and ionized absorption contribute significantly to the spectral shape. Accounting for their effects, we find no evidence for a relativistic reflection component, contrary to early work. Energy-dependent lags are significantly measured in the new data, but with an energy profile that does not resemble a broad iron line, in contrast to the old data. The complex lag-energy spectra, along with the lack of strong evidence for a relativistic spectral component, suggest that the energy-dependent lags are produced by absorption effects. The long term spectral variations provide new details on the variability of the narrow Fe Kα line . We find that its variations are correlated with, and delayed with respect to, the primary X-ray continuum. We measure a delay of tau = 3.3+1.8 −0.7 days, implying an origin in the inner broad line region (BLR). The delay is half the Hβ line delay, suggesting a geometry that differs slightly from the optical BLR.


### Description
The analysis is organized into several python notebooks, which sometimes call outside functions either from the my [toolset package `aztools`](https://zoghbi-a.github.io/aztools/), [the psd/lag calculation package](https://zoghbi-a.github.io/plag/) or the specefic helper scripts: 
- [`spec_helpers.py`](https://github.com/zoghbi-a/Revisiting-NGC-4151-Data/blob/master/spec_helpers.py): This contrain a collection of functions used in the spectral modeling. These are documented individually.
- [`timing_helpers.py`](https://github.com/zoghbi-a/Revisiting-NGC-4151-Data/blob/master/timing_helpers.pyy): Contains the functions used in the timing analysis, and these are mostly used in the [Timing](timing) (see bellow) notebook, and they are also documented individually
- [`fit.tcl`](https://github.com/zoghbi-a/Revisiting-NGC-4151-Data/blob/master/fit.tcl): This contrains a collection of tcl functions called from xspec to do the modeling. They are mostly called from the spectral modeling notebooks (see below).

A quick description of each notebooks is as follows:

- [Data](data.md): This is the data extraction tools, and contains the code for downloading, reducing the data, and extracting the spectra from XMM, Suzaku and NuSTAR used in the data.
- [Data Timing](data_timing.md): Contrains the code for extracting the light curves used in the timing analysis.
- [Prepare Spec](prepare_spec.md): This mostly does the gain correction discussed in the paper.
- [Explore Spec](explore_spec.md): A general exploratioin of the spectra.
- [Spec](spec): Starting from the the results in [Explore Spec](explore_spec), do the a consistent spectral modeling. This is where the spectral results in the paper come from.
-[Spec Narrow line](spec_narrow_line.md): Use the results from [Spec](spec) and look in more details at the variability of the narrow line. 
- [Subspec Narrow Line](subspec_narrow_line.md): Study the variability of the narrow line on the short (5 ks) time scales, using both the XMM and Suzaku data.
- [Spec Other Parameters](spec_other_params.md): Study the variability of parameters other than the narrow line resulting from the spectral modeling in [Spec](spec.md).
- [Timing](timing.md): Contains the code for extracting and modeling the power spectra and lags.
- [Lag Models](lag_models.md): Contains the code for producing the lag models presented in the discussion section of the paper.

### Data Products
All the data products are available through the Open Science Framework at the [following link](https://osf.io/x4jde/files/). There are three files:
- [`xmm_spec.tgz`](https://osf.io/a68e9/): contains the spectral products and modeling from individual observations.
- [`xmm_subspec.tgz`](https://osf.io/8sqcp/): contains the spectral products from 5 ks segments, used for the narrow line analysis.
- [`suzaku_subspec.tgz`](https://osf.io/f8buc/): contains the spectral products from 5 ks segments, used for the narrow line analysis.
- [`xmm_timing.tgz`](https://osf.io/86qbd/): contains all the timing products.
