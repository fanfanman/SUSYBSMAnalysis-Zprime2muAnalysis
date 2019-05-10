
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_SingleMuonRun2017F-31Mar2018-v1_data2017ebe'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_SingleMuonRun2017F-31Mar2018-v1_data2017ebe'
config.Data.outLFNDirBase = '/store/user/zhangfa/Data2017eventbyevent'
config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.whitelist = ["T2_US_*"]
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 400
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt'

