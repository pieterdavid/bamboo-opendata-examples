"""
Implementation of the IRIS-HEP ADL benchmarks https://github.com/iris-hep/adl-benchmarks-index
using the bamboo framework https://gitlab.cern.ch/cp3-cms/bamboo

After installing bamboo (see the documentation, plotIt is not required) this can be run with
$ bambooRun -m bamboo/examples/adl_benchmarks.py:IRISHEP_ADLBenchmarks bamboo/examples/adl_benchmarks.yml -o adl_bamboo
the resulting histograms will be stored in `adl_bamboo/results/SingleMuon_test.root`.
The `--examples` argument can be used to run just some of the benchmarks (default is all);
since some of the benchmarks go through a lot of combinatorics, implicit multithreading
is useful to speed them up, and can be enabled by passing `--threads N` (with N=5 or so).
"""
__author__ = "Pieter David <pieter.david@cern.ch>"

from bamboo.analysismodules import NanoAODHistoModule
import bamboo.logging
logger = bamboo.logging.getLogger(__name__)

from bamboo.treedecorators import NanoAODDescription
description_CMSRun1OpenData_ROOTbenchmark = NanoAODDescription(groups=["HLT_", "PV_", "MET_"], collections=["nMuon", "nElectron", "nTau", "nPhoton", "nJet"])

class IRISHEP_ADLBenchmarks(NanoAODHistoModule):
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        ## disabled, no need for plotIt then
        logger.info("The resulting histograms are stored in the {0} directory".format(resultsdir))
    def mergeCounters(self, outF, infileNames, sample=None):
        pass ## disabled because the test file is not a full NanoAOD
    def addArgs(self, parser):
        super(IRISHEP_ADLBenchmarks, self).addArgs(parser)
        self.nExamples = 8
        parser.add_argument("--examples", type=str, default="all", help='Set of examples to run (comma-separated, any set of 1-{0:d}), or "all" (default)'.format(self.nExamples))
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        return super(IRISHEP_ADLBenchmarks, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, description=description_CMSRun1OpenData_ROOTbenchmark)
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, SummedPlot
        from bamboo.plots import EquidistantBinning as EqBin
        from bamboo import treefunctions as op

        if self.args.examples == "all":
            examples = list(i+1 for i in range(self.nExamples)) # 1-4 are fine, so is 7
        else:
            examples = list(set(int(tk) for tk in self.args.examples.split(",")))
        logger.info("Running the following examples: {0}".format(",".join(str(i) for i in examples)))
                
        plots = []
        if 1 in examples:
            ## Example 1: Plot the missing ET of all events.
            plots.append(Plot.make1D("Ex1_MET",
                tree.MET.pt, noSel,
                EqBin(100, 0., 2000.), title="MET (GeV)"))

        if 2 in examples:
            ## Example 2: Plot pT of all jets in all events.
            plots.append(Plot.make1D("Ex2_jetPt",
                op.map(tree.Jet, lambda j : j.pt), noSel,
                EqBin(100, 15., 60.), title="Jet p_{T} (GeV/c)"))

        if 3 in examples:
            ## Example 3: Plot pT of jets with |η| < 1.
            centralJets1 = op.select(tree.Jet, lambda j : op.abs(j.eta) < 1.)
            plots.append(Plot.make1D("Ex3_central1_jetPt",
                op.map(centralJets1, lambda j : j.pt), noSel,
                EqBin(100, 15., 60.), title="Jet p_{T} (GeV/c)"))

        if 4 in examples:
            ## Example 4: Plot the missing ET of events that have at least two jets with pT > 40 GeV.
            jets40 = op.select(tree.Jet, lambda j : j.pt > 40)
            hasTwoJets40 = noSel.refine("twoJets40", cut=(op.rng_len(jets40) >= 2))
            plots.append(Plot.make1D("Ex4_twoJets40_MET",
                tree.MET.pt, hasTwoJets40,
                EqBin(100, 0., 2000.), title="MET (GeV)"))

        if 5 in examples:
            ## Example 5: Plot the missing ET of events that have an opposite-sign muon pair with an invariant mass between 60 and 120 GeV.
            dimu_Z = op.combine(tree.Muon, N=2, pred=(lambda mu1, mu2 : op.AND(
                mu1.charge != mu2.charge,
                op.in_range(60., op.invariant_mass(mu1.p4, mu2.p4), 120.)
                )))
            hasDiMuZ = noSel.refine("hasDiMuZ", cut=(op.rng_len(dimu_Z) > 0))
            plots.append(Plot.make1D("Ex5_dimuZ_MET",
                tree.MET.pt, hasDiMuZ,
                EqBin(100, 0., 2000.), title="MET (GeV)"))

        if 6 in examples:
            ## Example 6: Plot pT of the trijet system with the mass closest to 172.5 GeV in each event and plot the maximum b-tagging discriminant value among the jets in the triplet.
            trijets = op.combine(tree.Jet, N=3)
            hadTop = op.rng_min_element_by(trijets,
                fun=lambda comb: op.abs((comb[0].p4+comb[1].p4+comb[2].p4).M()-172.5))
            hadTop_p4 = (hadTop[0].p4 + hadTop[1].p4 + hadTop[2].p4)
            hasTriJet = noSel.refine("hasTriJet", cut=(op.rng_len(trijets) > 0))
            plots.append(Plot.make1D("Ex6_trijet_topPt",
                hadTop_p4.Pt(), hasTriJet,
                EqBin(100, 15., 40.), title="Trijet p_{T} (GeV/c)"))
            plots.append(Plot.make1D("Ex6_trijet_maxbtag",
                op.max(op.max(hadTop[0].btag, hadTop[1].btag), hadTop[2].btag), hasTriJet,
                EqBin(100, 0., 1.), title="Trijet maximum b-tag"))
            if verbose:
                plots.append(Plot.make1D("Ex6_njets",
                    op.rng_len(tree.Jet), noSel,
                    EqBin(20, 0., 20.), title="Number of jets"))
                plots.append(Plot.make1D("Ex6_ntrijets",
                    op.rng_len(trijets), noSel,
                    EqBin(100, 0., 1000.), title="Number of 3-jet combinations"))
                plots.append(Plot.make1D("Ex6_trijet_mass",
                    hadTop_p4.M(), hasTriJet,
                    EqBin(100, 0., 250.), title="Trijet mass (GeV/c^{2})"))

        if 7 in examples:
            ## Example 7: Plot the sum of pT of jets with pT > 30 GeV that are not within 0.4 in ΔR of any lepton with pT > 10 GeV.
            el10  = op.select(tree.Electron, lambda el : el.pt > 10.)
            mu10  = op.select(tree.Muon    , lambda mu : mu.pt > 10.)
            cleanedJets30 = op.select(tree.Jet, lambda j : op.AND(
                j.pt > 30.,
                op.NOT(op.rng_any(el10, lambda el : op.deltaR(j.p4, el.p4) < 0.4 )),
                op.NOT(op.rng_any(mu10, lambda mu : op.deltaR(j.p4, mu.p4) < 0.4 ))
                ))
            plots.append(Plot.make1D("Ex7_sumCleanedJetPt",
                op.rng_sum(cleanedJets30, lambda j : j.pt), noSel,
                EqBin(100, 15., 200.), title="Sum p_{T} (GeV/c)"))

        if 8 in examples:
            ## Example 8: For events with at least three leptons and a same-flavor opposite-sign lepton pair, find the same-flavor opposite-sign lepton pair with the mass closest to 91.2 GeV and plot the transverse mass of the missing energy and the leading other lepton.
            # The plot is made for each of the different flavour categories (l+/- l-/+ l') and then summed,
            # because concatenation of containers is not (yet) supported.
            lepColl = { "El" : tree.Electron, "Mu" : tree.Muon }
            mt3lPlots = []
            for dlNm,dlCol in lepColl.items():
                dilep = op.combine(dlCol, N=2, pred=(lambda l1,l2 : op.AND(l1.charge != l2.charge)))
                hasDiLep = noSel.refine("hasDilep{0}{0}".format(dlNm), cut=(op.rng_len(dilep) > 0))
                dilepZ = op.rng_min_element_by(dilep, fun=lambda ll : op.abs(op.invariant_mass(ll[0].p4, ll[1].p4)-91.2))
                for tlNm,tlCol in lepColl.items():
                    if tlCol == dlCol:
                        hasTriLep = hasDiLep.refine("hasTrilep{0}{0}{1}".format(dlNm,tlNm),
                            cut=(op.rng_len(tlCol) > 2))
                        residLep = op.select(tlCol, lambda l : op.AND(l.idx != dilepZ[0].idx, l.idx != dilepZ[1].idx))
                        l3 = op.rng_max_element_by(residLep, lambda l : l.pt)
                    else:
                        hasTriLep = hasDiLep.refine("hasTriLep{0}{0}{1}".format(dlNm,tlNm),
                            cut=(op.rng_len(tlCol) > 0))
                        l3 = op.rng_max_element_by(tlCol, lambda l : l.pt)
                    mtPlot = Plot.make1D("Ex8_3lMT_{0}{0}{1}".format(dlNm,tlNm),
                        op.sqrt(2*l3.pt*tree.MET.pt*(1-op.cos(l3.phi-tree.MET.phi))), hasTriLep,
                        EqBin(100, 15., 250.), title="M_{T} (GeV/c^2)")
                    mt3lPlots.append(mtPlot)
                    plots.append(mtPlot)
            plots.append(SummedPlot("Ex8_3lMT", mt3lPlots))

        return plots
