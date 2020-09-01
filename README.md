# [bamboo](https://gitlab.cern.ch/cp3-cms/bamboo) examples based on [RDataFrame tutorials](https://root.cern/doc/master/group__tutorial__dataframe.html) with open data

[![Binder](https://mybinder.org/badge_logo.svg)](https://github.com/pieterdavid/bamboo-dockerhub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Fpieterdavid%2Fbamboo-opendata-examples.git&urlpath=lab%2Ftree%2Fbamboo-opendata-examples%2Fhiggs4l_tutorial_CMSOpenData.py&branch=master)

These can be run directly on [Binder](https://mybinder.readthedocs.io/en/latest/) through the badge above, or installed locally (see the [documentation](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/install.html) for more details).

Currently, the following examples are available:

- CMS H->4l search, based on [the df103 tutorial](https://root.cern/doc/master/df103__NanoAODHiggsAnalysis_8C.html)
- ATLAS W boson analysis, based on [the df105 tutorial](https://root.cern/doc/master/df105__WBosonAnalysis_8py.html)
- the [IRIS-HEP analysis description language benchmarks](https://github.com/iris-hep/adl-benchmarks-index), also available in a [separate repository](https://github.com/pieterdavid/bamboo-adl-benchmarks)

Some of these, especially of the ADL benchmarks, go through enough combinatorics to benefit from [implicit multithreading](https://doi.org/10.1088/1742-6596/898/7/072022) in [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) ([DOI 10.5281/zenodo.260230](https://doi.org/10.5281/zenodo.260230)); this can be enabled by passing `--threads N` (with `N` three or four).
