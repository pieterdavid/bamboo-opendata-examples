# [bamboo](https://gitlab.cern.ch/cp3-cms/bamboo) examples based on [RDataFrame tutorials](https://root.cern/doc/master/group__tutorial__dataframe.html) with open data

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pieterdavid/bamboo-docker/master?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252Fpieterdavid%252Fbamboo-opendata-examples%26urlpath%3Dlab%252Ftree%252Fbamboo-opendata-examples%252Fhiggs4l_tutorial_CMSOpenData.py%26branch%3Dmaster)

These can be run directly on [Binder](https://mybinder.readthedocs.io/en/latest/) through the badge above, or installed locally (see the [documentation](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/install.html) for more details).

In case the link above does not work, [this one](https://mybinder.org/v2/gh/pieterdavid/bamboo-docker/master) can be used to launch the session, where this repository can be cloned with ``git clone https://github.com/pieterdavid/bamboo-opendata-examples.git``.

Currently, the following examples are available:

- CMS H->4l search, based on [the df103 tutorial](https://root.cern/doc/master/df103__NanoAODHiggsAnalysis_8C.html).
  The [default](higgs4l_tutorial_CMSOpenData.yml) configuration file uses the skimmed (51MB) input data,
  for performance testing the full (12GB) version
  [`higgs4l_tutorial_CMSOpenData_full.yml`](higgs4l_tutorial_CMSOpenData_full.yml)
  (if file access is fast enough multithreading, e.g. ``--threads 4``, works very nicely)
  ```sh
  bambooRun -m higgs4l_tutorial_CMSOpenData.py:Higgs4L higgs4l_tutorial_CMSOpenData.yml -o test_out/df103
  ```
- ATLAS W boson analysis, based on [the df105 tutorial](https://root.cern/doc/master/df105__WBosonAnalysis_8py.html)
  ```sh
  bambooRun -m wanalysis_tutorial_ATLASOpenData.py:WAnalysis wanalysis_tutorial_ATLASOpenData.yml -o test_out/df105
  ```
- the [IRIS-HEP analysis description language benchmarks](https://github.com/iris-hep/adl-benchmarks-index), also available in a [separate repository](https://github.com/pieterdavid/bamboo-adl-benchmarks)
  ```sh
  bambooRun -m adl_benchmarks.py:IRISHEP_ADLBenchmarks adl_benchmarks.yml -o test_out/adl_benchmarks
  ```

Some of these, especially of the ADL benchmarks, go through enough combinatorics to benefit from [implicit multithreading](https://doi.org/10.1088/1742-6596/898/7/072022) in [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) ([DOI 10.5281/zenodo.260230](https://doi.org/10.5281/zenodo.260230)); this can be enabled by passing `--threads N` (with `N` three or four).
