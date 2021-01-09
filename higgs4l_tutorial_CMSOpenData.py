"""
Implementation of the H->ZZ*->4l analysis on CMS 2012 open data,
see http://opendata.cern.ch/record/5500 and the RDataFrame implementation in
https://root.cern/doc/master/df103__NanoAODHiggsAnalysis_8C.html
"""
from bamboo.analysismodules import NanoAODHistoModule
from functools import partial

from bamboo.treedecorators import NanoAODDescription
description_CMSRun1OpenData_ROOT_H4ltutorial = NanoAODDescription(groups=["PV_"], collections=["nMuon", "nElectron"])

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

class Higgs4L(NanoAODHistoModule):
    def mergeCounters(self, outF, infileNames, sample=None):
        pass ## disabled because the test file is not a full NanoAOD
    def isMC(self, sampleName):
        return sampleName.split("_")[0] not in ("DoubleMuParked", "DoubleElectron")
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        return super(Higgs4L, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, description=description_CMSRun1OpenData_ROOT_H4ltutorial)
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        muons = op.sort(op.select(tree.Muon, lambda mu : op.AND(
            mu.pt > 5,
            op.abs(mu.eta) < 2.4,
            op.abs(mu.pfRelIso04_all) < 0.40,
            op.abs(mu.dxy) < 0.5,
            op.abs(mu.dz ) < 1.,
            op.sqrt(mu.dxy**2 + mu.dz**2)/op.sqrt(mu.dxyErr**2+mu.dzErr**2) < 4, ## SIP3D
            )), lambda mu : -mu.pt)
        electrons = op.sort(op.select(tree.Electron, lambda el : op.AND(
            el.pt > 7.,
            op.abs(el.eta) < 2.5,
            op.abs(el.pfRelIso03_all) < 0.40,
            op.abs(el.dxy) < 0.5,
            op.abs(el.dz ) < 1.,
            op.sqrt(el.dxy**2 + el.dz**2)/op.sqrt(el.dxyErr**2+el.dzErr**2) < 4, ## SIP3D
            )), lambda el : -el.pt)

        plots += self.controlPlots_2l(noSel, muons, electrons)

        mZ = 91.1876

        def reco_4l(leptons, lName, baseSel):
            ## select events with four leptons, and find the best Z candidate
            ## shared between 4el and 4mu
            has4l = baseSel.refine(f"has4{lName}", cut=[
                op.rng_len(leptons) == 4,
                op.rng_sum(leptons, lambda l : l.charge) == 0,
                ])
            allZcand = op.combine(leptons, N=2, pred=lambda l1,l2 : l1.charge != l2.charge)
            bestZ = op.rng_min_element_by(allZcand, lambda ll : op.abs(op.invariant_mass(ll[0].p4, ll[1].p4)-mZ))
            otherLeptons = op.select(leptons, partial(lambda l,oz=None : op.AND(l.idx != oz[0].idx, l.idx != oz[1].idx), oz=bestZ))
            return has4l, bestZ, otherLeptons

        ## Mixed category: take leading two for each (as in the other implementations
        def cuts2lMixed(leptons):
            return [ op.rng_len(leptons) == 2,
                     leptons[0].charge != leptons[1].charge,
                     leptons[0].pt > 20.,
                     leptons[1].pt > 10.
                   ]
        has4lMixed = noSel.refine(f"has4lMixed", cut=cuts2lMixed(muons)+cuts2lMixed(electrons))
        mixed_mElEl = op.invariant_mass(electrons[0].p4, electrons[1].p4)
        mixed_mMuMu = op.invariant_mass(muons[0].p4, muons[1].p4)
        ## two categories: elel closest to Z or mumu closest to Z
        has2El2Mu = has4lMixed.refine(f"has2El2Mu", cut=(op.abs(mixed_mElEl-mZ) < op.abs(mixed_mMuMu-mZ)))
        has2Mu2El = has4lMixed.refine(f"has2Mu2El", cut=(op.abs(mixed_mElEl-mZ) > op.abs(mixed_mMuMu-mZ)))

        mH_cats = []
        for catNm, (has4l, bestZ, otherZ) in {
                "4Mu" : reco_4l(muons    , "Mu", noSel),
                "4El" : reco_4l(electrons, "El", noSel),
                "2El2Mu" : (has2El2Mu, electrons, muons),
                "2Mu2El" : (has2Mu2El, muons, electrons)
                }.items():
            bestZp4  = bestZ[0].p4  + bestZ[1].p4
            otherZp4 = otherZ[0].p4 + otherZ[1].p4
            hasZZ = has4l.refine(f"{has4l.name}ZZ", cut=[
                op.deltaR(bestZ[0].p4 , bestZ[1].p4 ) > 0.02,
                op.deltaR(otherZ[0].p4, otherZ[1].p4) > 0.02,
                op.in_range(40., bestZp4.M(), 120.),
                op.in_range(12., otherZp4.M(), 120.)
                ])
            plots += self.controlPlots_4l(hasZZ, bestZ, otherZ)
            m4l = (bestZ[0].p4+bestZ[1].p4+otherZ[0].p4+otherZ[1].p4).M()
            hasZZ_m4l70 = hasZZ.refine(f"{hasZZ.name}m4l70", cut=(m4l > 70.))
            p_mH = Plot.make1D(f"H_mass_{catNm}", m4l, hasZZ_m4l70, EqBin(36, 70., 180.),
                    plotopts={"show-overflow": False, "log-y": False, "y-axis-range": [0., 18.]})
            mH_cats.append(p_mH)
            plots.append(p_mH)
        plots.append(SummedPlot("H_mass", mH_cats))

        return plots

    def controlPlots_2l(self, noSel, muons, electrons):
        plots = [
            Plot.make1D("nEl", op.rng_len(electrons), noSel, EqBin(10, 0., 10.), xTitle="Number of tight electrons"),
            Plot.make1D("nMu", op.rng_len(muons), noSel, EqBin(10, 0., 10.), xTitle="Number of tight muons"),
            ]
        hasOSElEl = noSel.refine("hasOSElEl", cut=[ op.rng_len(electrons) >= 2,
            electrons[0].charge != electrons[1].charge, electrons[0].pt > 20., electrons[1].pt > 10. ])
        plots.append(Plot.make1D("massZto2e", op.invariant_mass(electrons[0].p4, electrons[1].p4),
            hasOSElEl, EqBin(120, 40., 120.), title="mass of Z to 2e",
            xTitle="Invariant Mass of Nelectrons=2 (in GeV/c^2)"))
        plots += [ 
            Plot.make1D("OSElEl_PTl1", electrons[0].pt, hasOSElEl, EqBin(50, 0., 100.)),
            Plot.make1D("OSElEl_PTl2", electrons[1].pt, hasOSElEl, EqBin(50, 0., 100.)),
            ]
        hasOSMuMu = noSel.refine("hasOSMuMu", cut=[ op.rng_len(muons) >= 2,
            muons[0].charge != muons[1].charge, muons[0].pt > 20., muons[1].pt > 10. ])
        plots.append(Plot.make1D("massZto2mu", op.invariant_mass(muons[0].p4, muons[1].p4),
            hasOSMuMu, EqBin(120, 40., 120.), title="mass of Z to 2mu",
            xTitle="Invariant Mass of Nmuons=2 (in GeV/c^2)"))
        plots += [ 
            Plot.make1D("OSMuMu_PTl1", muons[0].pt, hasOSMuMu, EqBin(50, 0., 100.)),
            Plot.make1D("OSMuMu_PTl2", muons[1].pt, hasOSMuMu, EqBin(50, 0., 100.)),
            ]
        return plots
    def controlPlots_4l(self, sel4l, bestZ, otherZ):
        prefix = sel4l.name[3:]
        plots = [
            Plot.make1D(f"{prefix}_Z1_l1PT", bestZ[0].pt, sel4l, EqBin(50, 0., 100.)),
            Plot.make1D(f"{prefix}_Z1_l2PT", bestZ[1].pt, sel4l, EqBin(50, 0., 100.)),
            Plot.make1D(f"{prefix}_Z2_l1PT", otherZ[0].pt, sel4l, EqBin(50, 0., 100.)),
            Plot.make1D(f"{prefix}_Z2_l2PT", otherZ[1].pt, sel4l, EqBin(50, 0., 100.)),
            ]
        return plots
