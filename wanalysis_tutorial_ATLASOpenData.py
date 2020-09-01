import bamboo.analysismodules

class ATLASOutreachOpenDataModule(bamboo.analysismodules.AnalysisModule):
    def isMC(self, sampleName):
        return not sampleName.startswith("data")
    def prepareTree(self, tree, sample=None, sampleCfg=None, lazyBackend=False):
        from bamboo.dataframebackend import DataframeBackend, LazyDataframeBackend
        backendCls = (LazyDataframeBackend if lazyBackend else DataframeBackend)
        t = decorateATLASOutreachOpenData(tree)
        be, noSel = backendCls.create(t)
        return t, noSel, be, None
class ATLASOutreachOpenDataHistoModule(ATLASOutreachOpenDataModule, bamboo.analysismodules.HistogramsModule):
    def __init__(self, args):
        super(ATLASOutreachOpenDataHistoModule, self).__init__(args)

def makePtEtaPhiEP4(pt, eta, phi, e):
    from bamboo.treeoperations import Construct
    return Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >", (pt, eta, phi, e)).result

class WAnalysis(ATLASOutreachOpenDataHistoModule):
    def prepareTree(self, tree, sample=None, sampleCfg=None, lazyBackend=False):
        t,noSel,be,lumiArgs = super(WAnalysis, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, lazyBackend=lazyBackend)
        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=[ t.scaleFactor.ELE, t.scaleFactor.MUON, t.scaleFactor.LepTRIGGER, t.scaleFactor.PILEUP, t.mcWeight ])
        return t,noSel,be,lumiArgs
    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, SummedPlot
        from bamboo.plots import EquidistantBinning as EqBin
        from bamboo import treefunctions as op
        plots = []
        metSel = noSel.refine("MET", cut=(tree.met.et > 30000))
        trigSel = metSel.refine("trig", cut=op.OR(tree.trigE, tree.trigM))
        goodLeptons = op.select(tree.lep, lambda l : op.AND(l.isTightID, l.pt > 35000., l.ptcone30/l.pt < 0.1, l.etcone20/l.pt < 0.1))
        oneLepSel = trigSel.refine("1goodlep", cut=(op.rng_len(goodLeptons) == 1))
        lep = goodLeptons[0]
        signalSel = oneLepSel.refine("signalRegion", cut=op.AND(
            op.abs(lep.z0*op.sin(lep.p4.theta())) < 0.5,
            op.multiSwitch(
                (lep.type == 11, op.AND(
                    op.abs(lep.eta) < 2.46, op.NOT(op.in_range(1.37, op.abs(lep.eta), 1.52)),
                    op.abs(lep.trackd0pvunbiased/lep.tracksigd0pvunbiased) < 5
                    )),
                (lep.type == 13, op.AND(
                    op.abs(lep.eta) < 2.5,
                    op.abs(lep.trackd0pvunbiased/lep.tracksigd0pvunbiased) < 3
                    )),
                op.c_bool(False)
            )))
        metp4 = makePtEtaPhiEP4(tree.met.et, op.c_float(0.), tree.met.phi, tree.met.et)
        plots.append(Plot.make1D("mt_w", (lep.p4+metp4).Mt()/1000., signalSel, EqBin(40, 60., 180.), title="m_{T,W}"))

        return plots


def decorateATLASOutreachOpenData(aTree):
    """ Simple decorator for ATLAS OpenData outreach trees - a simplified copy from bamboo.treedecorators.decorateNanoAOD (e.g. without systematics) """
    groups = ["scaleFactor_", "met_"]
    collections = ["lep_n", "jet_n", "photon_n", "tau_n", "largeRjet_n"]
    from functools import partial
    from bamboo.treedecorators import allLeafs, proxy, funProxy, itemProxy
    from bamboo.treeoperations import GetColumn, GetArrayColumn, Construct
    from bamboo.treeproxies import ContainerGroupItemProxy, ContainerGroupProxy, LeafGroupProxy, TreeBaseProxy

    def addP4ToObj(prefix, lvNms):
        return funProxy(partial( (lambda getEta,getE,inst:
                   Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >",
                             (inst.pt, getEta(inst), inst.phi, getE(inst))).result),
                   ((lambda inst: inst.eta) if f"{prefix}eta" in lvNms else (lambda inst: 0.)),
                   ((lambda inst: inst.E) if f"{prefix}E" in lvNms else (lambda inst: 0.)) ))

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm))) for lvNm,lv in allTreeLeafs.items()))
    tree_children = list()
    def setTreeAtt(name, proxy, setParent=True):
        tree_dict[name] = proxy
        if setParent:
            tree_children.append(proxy)

    ## non-collection branches to group
    grp_found = []
    for prefix in groups:
        if not any(lvNm.startswith(prefix) for lvNm in allTreeLeafs):
            logger.warning("No branch name starting with {0} in the tree - skipping group".format(prefix))
        else:
            grp_found.append(prefix)
    for prefix in grp_found:
        grpNm = prefix.rstrip("_")
        grp_dict = {
            "__doc__" : "{0} leaf group proxy class".format(grpNm)
            }
        grp_lvNms = set(lvNm for lvNm in allTreeLeafs.keys() if lvNm.startswith(prefix))
        grp_dict.update(dict((lvNm[len(prefix):], proxy(GetColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm))) for lvNm in grp_lvNms))
        if f"{prefix}pt" in grp_lvNms and f"{prefix}phi" in grp_lvNms:
            grp_dict["p4"] = addP4ToObj(prefix, grp_lvNms)
        grpcls = type("{0}LeafGroupProxy".format(grpNm), (LeafGroupProxy,), grp_dict)
        for lvNm in grp_lvNms:
            del tree_dict[lvNm]
        ## default group proxy, replaced below if needed
        grp_proxy = grpcls(grpNm, None)
        setTreeAtt(grpNm, grp_proxy)

    ## SOA, nanoAOD style (LeafCount, shared)
    cnt_found = []
    for sizeNm in collections:
        if sizeNm not in allTreeLeafs:
            logger.warning("{0} is not a branch in the tree - skipping collection".format(sizeNm))
        else:
            cnt_found.append(sizeNm)

    from bamboo.treeproxies import vecPat, makeProxy
    for sizeNm in cnt_found:
        grpNm = "_".join(sizeNm.split("_")[:-1]) ## strip "_n"
        prefix = "{0}_".format(grpNm)
        itm_dict = {
            "__doc__" : "{0} proxy class".format(grpNm)
            }
        itm_lvs = set(lvNm for lvNm,lv in allTreeLeafs.items() if lvNm.startswith(prefix) and lvNm != sizeNm)
        sizeOp = GetColumn(allTreeLeafs[sizeNm].GetTypeName(), sizeNm)
        for lvNm in itm_lvs:
            lvNm_short = lvNm[len(prefix):]
            m = vecPat.match(allTreeLeafs[lvNm].GetTypeName())
            if not m:
                raise RuntimeError("problem with interpreting as vector")
            col = GetArrayColumn(m.group("item"), lvNm, sizeOp).result
            itm_dict[lvNm_short] = itemProxy(col)
        ## create p4 branches (naive, but will be reused for variation case)
        if f"{prefix}pt" in itm_lvs and f"{prefix}phi" in itm_lvs:
            itm_dict["p4"] = addP4ToObj(prefix, itm_lvs)
        itm_bases = [ContainerGroupItemProxy]
        itmcls = type("{0}GroupItemProxy".format(grpNm), tuple(itm_bases), itm_dict)
        ## default collection proxy, replaced below if needed
        coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
        setTreeAtt(grpNm, coll_orig)
        for lvNm in itm_lvs:
            del tree_dict[lvNm]
        del tree_dict[sizeNm] ## go through op.rng_len

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)
    for pc in tree_children:
        pc._parent = treeProxy

    return treeProxy
