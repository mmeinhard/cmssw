import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.jets_cff import *
from PhysicsTools.NanoAOD.muons_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.electrons_cff import *
from PhysicsTools.NanoAOD.photons_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.extraflags_cff import *
from PhysicsTools.NanoAOD.ttbarCategorization_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.particlelevel_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.NanoAOD.triggerObjects_cff import *
from PhysicsTools.NanoAOD.isotracks_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.boosted_cff import *

from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_nanoAOD_92X_cff import run2_nanoAOD_92X

nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

linkedObjects = cms.EDProducer("PATObjectCrossLinker",
   jets=cms.InputTag("finalJets"),
   muons=cms.InputTag("finalMuons"),
   electrons=cms.InputTag("finalElectrons"),
   taus=cms.InputTag("finalTaus"),
   photons=cms.InputTag("finalPhotons"),
)

simpleCleanerTable = cms.EDProducer("NanoAODSimpleCrossCleaner",
   name=cms.string("cleanmask"),
   doc=cms.string("simple cleaning mask with priority to leptons"),
   jets=cms.InputTag("linkedObjects","jets"),
   muons=cms.InputTag("linkedObjects","muons"),
   electrons=cms.InputTag("linkedObjects","electrons"),
   taus=cms.InputTag("linkedObjects","taus"),
   photons=cms.InputTag("linkedObjects","photons"),
   jetSel=cms.string("pt>15"),
   muonSel=cms.string("isPFMuon && innerTrack.validFraction >= 0.49 && ( isGlobalMuon && globalTrack.normalizedChi2 < 3 && combinedQuality.chi2LocalPosition < 12 && combinedQuality.trkKink < 20 && segmentCompatibility >= 0.303 || segmentCompatibility >= 0.451 )"),
   electronSel=cms.string(""),
   tauSel=cms.string(""),
   photonSel=cms.string(""),
   jetName=cms.string("Jet"),muonName=cms.string("Muon"),electronName=cms.string("Electron"),
   tauName=cms.string("Tau"),photonName=cms.string("Photon")
)

btagSFdir="PhysicsTools/NanoAOD/data/btagSF/"

btagWeightTable = cms.EDProducer("BTagSFProducer",
    src = cms.InputTag("linkedObjects","jets"),
    cut = cms.string("pt > 25. && abs(eta) < 2.5"),
    discNames = cms.vstring(
        "pfCombinedInclusiveSecondaryVertexV2BJetTags",
        "pfDeepCSVJetTags:probb+pfDeepCSVJetTags:probbb",       #if multiple MiniAOD branches need to be summed up (e.g., DeepCSV b+bb), separate them using '+' delimiter
        "pfCombinedMVAV2BJetTags"        
    ),
    discShortNames = cms.vstring(
        "CSVV2",
        "DeepCSVB",
        "CMVA"
    ),
    weightFiles = cms.vstring(                                  #default settings are for 2017 94X. toModify function is called later for other eras.
        btagSFdir+"CSVv2_94XSF_V2_B_F.csv",
        btagSFdir+"DeepCSV_94XSF_V2_B_F.csv",                    
        "unavailable"                                           #if SFs for an algorithm in an era is unavailable, the corresponding branch will not be stored
    ),
    operatingPoints = cms.vstring("3","3","3"),                 #loose = 0, medium = 1, tight = 2, reshaping = 3
    measurementTypesB = cms.vstring("iterativefit","iterativefit","iterativefit"),     #e.g. "comb", "incl", "ttbar", "iterativefit"
    measurementTypesC = cms.vstring("iterativefit","iterativefit","iterativefit"),
    measurementTypesUDSG = cms.vstring("iterativefit","iterativefit","iterativefit"),
    sysTypes = cms.vstring("central","central","central")
)

run2_miniAOD_80XLegacy.toModify(btagWeightTable,                
    cut = cms.string("pt > 25. && abs(eta) < 2.4"),             #80X corresponds to 2016, |eta| < 2.4
    weightFiles = cms.vstring(                                  #80X corresponds to 2016 SFs
        btagSFdir+"CSVv2_Moriond17_B_H.csv",            
        "unavailable",                    
        btagSFdir+"cMVAv2_Moriond17_B_H.csv"                                            
    )
)

run2_nanoAOD_92X.toModify(btagWeightTable,                      #92X corresponds to MCv1, for which SFs are unavailable
    weightFiles = cms.vstring(
        "unavailable",
        "unavailable",                    
        "unavailable"                                            
    )
)                    

genWeightsTable = cms.EDProducer("GenWeightsTableProducer",
    genEvent = cms.InputTag("generator"),
    lheInfo = cms.InputTag("externalLHEProducer"),
    preferredPDFs = cms.VPSet( # see https://lhapdf.hepforge.org/pdfsets.html
        cms.PSet( name = cms.string("NNPDF31_nnlo_hessian_pdfas"), lhaid = cms.uint32(306000) ),
        cms.PSet( name = cms.string("PDF4LHC15_nnlo_30_pdfas"), lhaid = cms.uint32(91400) ),
        cms.PSet( name = cms.string("NNPDF30_nlo_as_0118"), lhaid = cms.uint32(260000) ), # for some 92X samples. Note that the nominal weight, 260000, is not included in the LHE ...
        cms.PSet( name = cms.string("NNPDF30_lo_as_0130"), lhaid = cms.uint32(262000) ), # some MLM 80X samples have only this (e.g. /store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root )
        cms.PSet( name = cms.string("NNPDF30_nlo_nf_4_pdfas"), lhaid = cms.uint32(292000) ), # some FXFX 80X samples have only this (e.g. WWTo1L1Nu2Q, WWTo4Q)
        cms.PSet( name = cms.string("NNPDF30_nlo_nf_5_pdfas"), lhaid = cms.uint32(292200) ), # some FXFX 80X samples have only this (e.g. DYJetsToLL_Pt, WJetsToLNu_Pt, DYJetsToNuNu_Pt)
    ),
    namedWeightIDs = cms.vstring(),
    namedWeightLabels = cms.vstring(),
    lheWeightPrecision = cms.int32(14),
    maxPdfWeights = cms.uint32(150), 
    debug = cms.untracked.bool(False),
)
lheInfoTable = cms.EDProducer("LHETablesProducer",
    lheInfo = cms.InputTag("externalLHEProducer"),
    precision = cms.int32(14),
    storeLHEParticles = cms.bool(True) 
)

l1bits=cms.EDProducer("L1TriggerResultsConverter", src=cms.InputTag("gtStage2Digis"), legacyL1=cms.bool(False))

content = cms.EDAnalyzer("EventContentAnalyzer")

nanoSequence = cms.Sequence(
        nanoMetadata + jetSequence + muonSequence + tauSequence + electronSequence+photonSequence+vertexSequence+metSequence+prefiringweight+boostedSequence+
        isoTrackSequence + # must be after all the leptons 
        linkedObjects  +
        jetTables + muonTables + tauTables + electronTables + photonTables +  globalTables +vertexTables+ metTables+simpleCleanerTable + triggerObjectTables + isoTrackTables +
	l1bits + boostedTables)

nanoSequenceMC = cms.Sequence(genParticleSequence + particleLevelSequence + nanoSequence + jetMC + muonMC + electronMC + photonMC + tauMC + metMC + ttbarCatMCProducers +  globalTablesMC + btagWeightTable + genWeightsTable + genParticleTables + particleLevelTables + lheInfoTable  + ttbarCategoryTable + L1PrefiringTable + boostedTables + boostedMC)


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
def nanoAOD_addDeepBTagFor80X(process):
    print "Updating process to run DeepCSV btag on legacy 80X datasets"
    updateJetCollection(
               process,
               jetSource = cms.InputTag('slimmedJets'),
               jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None'),
               btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc'], ## to add discriminators
               btagPrefix = ''
           )
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.looseJetId.src="selectedUpdatedPatJets"
    process.tightJetId.src="selectedUpdatedPatJets"
    process.tightJetIdLepVeto.src="selectedUpdatedPatJets"
    process.bJetVars.src="selectedUpdatedPatJets"
    process.slimmedJetsWithUserData.src="selectedUpdatedPatJets"
    process.qgtagger80x.srcJets="selectedUpdatedPatJets"
    patAlgosToolsTask = getPatAlgosToolsTask(process)
    patAlgosToolsTask .add(process.updatedPatJets)
    patAlgosToolsTask .add(process.patJetCorrFactors)
    process.additionalendpath = cms.EndPath(patAlgosToolsTask)
    return process



from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
def nanoAOD_recalibrateMETs(process,isData):
    runMetCorAndUncFromMiniAOD(process,isData=isData)
    process.nanoSequence.insert(process.nanoSequence.index(process.jetSequence),cms.Sequence(process.fullPatMetSequence))
    process.basicJetsForMetForT1METNano = process.basicJetsForMet.clone(
        src = process.updatedJetsWithUserData.src,
        skipEM = False,
        type1JetPtThreshold = 0.0,
        calcMuonSubtrRawPtAsValueMap = cms.bool(True),
    )
    process.jetSequence.insert(process.jetSequence.index(process.updatedJetsWithUserData),cms.Sequence(process.basicJetsForMetForT1METNano))
    process.updatedJetsWithUserData.userFloats.muonSubtrRawPt = cms.InputTag("basicJetsForMetForT1METNano:MuonSubtrRawPt")
    process.corrT1METJetTable.src = process.finalJets.src
    process.corrT1METJetTable.cut = "pt<15 && abs(eta)<9.9"
    for table in process.jetTable, process.corrT1METJetTable:
        table.variables.muonSubtrFactor = Var("1-userFloat('muonSubtrRawPt')/(pt()*jecFactor('Uncorrected'))",float,doc="1-(muon-subtracted raw pt)/(raw pt)",precision=6)
    process.metTables += process.corrT1METJetTable
#    makePuppiesFromMiniAOD(process,True) # call this before in the global customizer otherwise it would reset photon IDs in VID
#    runMetCorAndUncFromMiniAOD(process,isData=isData,metType="Puppi",postfix="Puppi",jetFlavor="AK4PFPuppi")
#    process.puppiNoLep.useExistingWeights = False
#    process.puppi.useExistingWeights = False
#    process.nanoSequenceCommon.insert(process.nanoSequenceCommon.index(jetSequence),cms.Sequence(process.puppiMETSequence+process.fullPatMetSequencePuppi))
    return process

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
def nanoAOD_runMETfixEE2017(process,isData):
    runMetCorAndUncFromMiniAOD(process,
        isData=isData,
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139},
        postfix = "ModifiedMET")
    process.nanoSequence.insert(process.nanoSequence.index(jetSequence),process.fullPatMetSequenceFixEE2017)

def nanoAOD_customizeCommon(process):
    run2_miniAOD_80XLegacy.toModify(process, nanoAOD_addDeepBTagFor80X)
    return process


def nanoAOD_customizeData(process):
    process = nanoAOD_customizeCommon(process)
    process = nanoAOD_nanoAOD_runMETfixEE2017(process,isData=True)
    #for modifier in run2_nanoAOD_94XMiniAODv1, run2_nanoAOD_94XMiniAODv2:
    #    modifier.toModify(process, lambda p: nanoAOD_runMETfixEE2017(p,isData=True))
    if hasattr(process,'calibratedPatElectrons80X'):
        process.calibratedPatElectrons80X.isMC = cms.bool(False)
        process.calibratedPatPhotons80X.isMC = cms.bool(False)
    return process

def nanoAOD_customizeMC(process):
    process = nanoAOD_customizeCommon(process)
    process = nanoAOD_runMETfixEE2017(process,isData=False)
    #for modifier in run2_nanoAOD_94XMiniAODv1, run2_nanoAOD_94XMiniAODv2:
    #    modifier.toModify(process, lambda p: nanoAOD_runMETfixEE2017(p,isData=False))
    if hasattr(process,'calibratedPatElectrons80X'):
        process.calibratedPatElectrons80X.isMC = cms.bool(True)
        process.calibratedPatPhotons80X.isMC = cms.bool(True)
    return process

#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#def correctMET(process):
#  print ("HEHRERE")
#  runMetCorAndUncFromMiniAOD (
#        process,
#        isData = True, # false for MC
#        fixEE2017 = True,
#        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
#        postfix = "ModifiedMET"
#  )
#  return process

### Era dependent customization
_80x_sequence = nanoSequence.copy()
#remove stuff 
_80x_sequence.remove(isoTrackTable)
_80x_sequence.remove(isoTrackSequence)
#add stuff
_80x_sequence.insert(_80x_sequence.index(jetSequence), extraFlagsProducers)
_80x_sequence.insert(_80x_sequence.index(l1bits)+1, extraFlagsTable)

run2_miniAOD_80XLegacy.toReplaceWith( nanoSequence, _80x_sequence)

	

