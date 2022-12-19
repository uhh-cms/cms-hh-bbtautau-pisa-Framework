import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll
from PhysicsTools.NanoAOD.custom_jme_cff import AddParticleNetAK4Scores

def nanoAOD_addDeepInfoAK4CHS(process, addDeepBTag, addDeepFlavour, addParticleNet):
  _btagDiscriminators=[]
  if addDeepBTag:
    print("Updating process to run DeepCSV btag")
    _btagDiscriminators += ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc']
  if addDeepFlavour:
    print("Updating process to run DeepFlavour btag")
    _btagDiscriminators += ['pfDeepFlavourJetTags:probb','pfDeepFlavourJetTags:probbb','pfDeepFlavourJetTags:problepb','pfDeepFlavourJetTags:probc']
  if addParticleNet:
    print("Updating process to run ParticleNet btag")
    _btagDiscriminators += _pfParticleNetAK4JetTagsAll
  if len(_btagDiscriminators)==0: return process
  print("Will recalculate the following discriminators: "+", ".join(_btagDiscriminators))
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None'),
    btagDiscriminators = _btagDiscriminators,
    postfix = 'WithDeepInfo',
  )
  process.load("Configuration.StandardSequences.MagneticField_cff")
  process.jetCorrFactorsNano.src="selectedUpdatedPatJetsWithDeepInfo"
  process.updatedJets.jetSource="selectedUpdatedPatJetsWithDeepInfo"
  return process

def customise(process):
  process.MessageLogger.cerr.FwkReport.reportEvery = 100
  process.finalGenParticles.select = cms.vstring(
    # PISA settings:
    # "drop *",
    # "keep++ abs(pdgId) == 15", # keep full decay chain for taus
    # "+keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15", #keep leptons, with at most one mother back in the history
    # "+keep+ abs(pdgId) == 6 || abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 35 || abs(pdgId) == 39  || abs(pdgId) == 9990012 || abs(pdgId) == 9900012",   # keep VIP particles
    # "drop abs(pdgId)= 2212 && abs(pz) > 1000", #drop LHC protons accidentally added by previous keeps

    # taken from https://github.com/cms-sw/cmssw/blob/9318eb0ff5982d62e93f385ffbe72e04887c93c8/PhysicsTools/NanoAOD/python/genparticles_cff.py
    "drop *",
    "keep++ abs(pdgId) == 15 & (pt > 15 ||  isPromptDecayed() )",#  keep full tau decay chain for some taus
    #"drop status==1 & pt < 1", #drop soft stable particle in tau decay
    "keep+ abs(pdgId) == 15 ",  #  keep first gen decay product for all tau
    "+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())", # keep gamma above 10 GeV (or all prompt) and its first parent
    "+keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15", #keep leptons, with at most one mother back in the history
    "drop abs(pdgId)= 2212 && abs(pz) > 1000", #drop LHC protons accidentally added by previous keeps
    "keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)", #keep all B and C hadrons
    "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",   # keep neutrinos
    "keep status == 3 || (status > 20 && status < 30)", #keep matrix element summary
    "keep isHardProcess() ||  fromHardProcessDecayed()  || fromHardProcessFinalState() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())",  #keep event summary based on status flags
    "keep  (status > 70 && status < 80 && pt > 15) ", # keep high pt partons right before hadronization
    "keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 37 ",   # keep VIP(articles)s
    #"keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ",                                                     # keep K0
    "keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)", #keep SUSY fiction particles
  )

  process = nanoAOD_addDeepInfoAK4CHS(process, False, False, True)
  process = AddParticleNetAK4Scores(process, 'jetTable')

  process.boostedTauTable.variables.dxy = Var("leadChargedHadrCand().dxy()", float,
    doc="d_{xy} of lead track with respect to PV, in cm (with sign)", precision=10)
  process.boostedTauTable.variables.dz = Var("leadChargedHadrCand().dz()", float,
    doc="d_{z} of lead track with respect to PV, in cm (with sign)", precision=14)
  return process
