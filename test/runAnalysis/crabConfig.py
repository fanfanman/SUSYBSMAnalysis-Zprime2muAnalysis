
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_2018_dyPt_4Jets_2018_2018ebeVsPt'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_2018_dyPt_4Jets_2018_2018ebeVsPt'
config.Data.outLFNDirBase = '/store/user/zhangfa/DYMC2018ebeVsPt'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

