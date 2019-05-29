
from CRABClient.UserUtilities import config
config = config()
config.General.requestName = 'dileptonAna_resolution_2018_dy800to1400_2018ebeRelMass'
config.General.workArea = 'crab2'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmssw_cfg.py'   
config.Data.inputDataset =  '/ZToMuMu_NNPDF31_13TeV-powheg_M_800_1400/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.publication = False
config.Data.outputDatasetTag = 'dileptonAna_resolution_2018_dy800to1400_2018ebeRelMass'
config.Data.outLFNDirBase = '/store/user/zhangfa/DYMC2018ebeRelMass'
#config.Data.ignoreLocality = True
#config.General.instance = 'preprod' 
config.Site.storageSite = 'T3_US_FNALLPC'
config.JobType.maxMemoryMB  = 8000

config.Data.splitting = 'EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 500000

