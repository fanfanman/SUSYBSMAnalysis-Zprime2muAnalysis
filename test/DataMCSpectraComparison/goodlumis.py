import os
from FWCore.PythonUtilities.LumiList import LumiList

def combine(Sept17_ll, Prompt_ll):
    Prompt_ll.removeRuns(xrange(132440,144114+1))
    return Sept17_ll | Prompt_ll

Sept17_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_132440-144114_7TeV_Sep17ReReco_Collisions10_JSON.txt')
#Prompt_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v2.txt')  # gives 35/pb
Prompt_ll = LumiList('/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.txt')   # gives 33/pb
Run2010AB_ll = combine(Sept17_ll, Prompt_ll)

cmssw_base = os.environ['CMSSW_BASE']
Sept17MuonsOnly_ll = LumiList(os.path.join(cmssw_base, 'src/HeavyFlavorAnalysis/Onia2MuMu/certification/7TeV/Collisions10/Reprocessing/Cert_132440-144114_7TeV_Sep17ReReco_Collisions10_JSON_BPAG.txt'))
PromptMuonsOnly_ll = LumiList(os.path.join(cmssw_base, 'src/HeavyFlavorAnalysis/Onia2MuMu/certification/7TeV/Collisions10/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_BPAG_v2.txt'))
Run2010ABMuonsOnly_ll = combine(Sept17MuonsOnly_ll, PromptMuonsOnly_ll)

Run2010AB = Run2010AB_ll.getCMSSWString().split(',')
Run2010ABMuonsOnly = Run2010ABMuonsOnly_ll.getCMSSWString().split(',')


