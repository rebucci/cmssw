import datetime

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, StorageConfiguration, Workflow

project = 'Emulator_truncOFF_20180426'
#version = datetime.datetime.now().strftime('%Y%m%d_%H%M')

storage = StorageConfiguration(
    output=[
        "hdfs://eddie.crc.nd.edu:19000/store/user/$USER/" + project,
        "file:///hadoop/store/user/$USER/" + project,
        "root://deepthought.crc.nd.edu//store/user/$USER/" + project,
        "gsiftp://T3_US_NotreDame/store/user/$USER/" + project,
        "srm://T3_US_NotreDame/store/user/$USER/" + project,
    ]
)

workflows = []

TTBar_PU0 = Workflow(
    label='TTBar_PU0',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValTTbar_14TeV/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='ttbar0',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(TTBar_PU0)

TTBar_PU140 = Workflow(
    label='TTBar_PU140',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValTTbar_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='ttbar140',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(TTBar_PU140)

TTBar_PU200 = Workflow(
    label='TTBar_PU200',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValTTbar_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='ttbar200',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(TTBar_PU200)

QCD_PU0 = Workflow(
    label='QCD_PU0',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='qcd0',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(QCD_PU0)

QCD_PU140 = Workflow(
    label='QCD_PU140',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='qcd140',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(QCD_PU140)

QCD_PU200 = Workflow(
    label='QCD_PU200',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValQCD_Pt-15To7000_Flat_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='qcd200',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(QCD_PU200)

SingleMu_PU0 = Workflow(
    label='SingleMu_PU0',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValSingleMuPt10Extended/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='mu0',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(SingleMu_PU0)

SingleMu_PU140 = Workflow(
    label='SingleMu_PU140',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValSingleMuPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='mu140',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(SingleMu_PU140)

SingleMu_PU200 = Workflow(
    label='SingleMu_PU200',
    sandbox=cmssw.Sandbox(release='~/Tracklet_Emulator/CMSSW_9_3_2', include=['L1Trigger/TrackFindingTracklet/test/']),
    dataset=cmssw.Dataset(
        dataset='/RelValSingleMuPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/GEN-SIM-DIGI-RAW',
        events_per_task=50000
    ),
    category=Category(
        name='mu200',
        cores=1,
        memory=8000,
        disk=32000,
    ),
    command='cmsRun L1TrackNtupleMaker_cfg.py',
    publish_label='test',
    merge_size='2.5G',
    merge_command="hadd @outputfiles @inputfiles",
    outputs=['output.root']
)
workflows.append(SingleMu_PU200)

config = Config(
    workdir='/tmpscratch/users/$USER/' + project,
    plotdir='/afs/crc.nd.edu/user/r/rbucci/www/lobster/' + project,
    storage=storage,
    workflows=workflows,
    advanced=AdvancedOptions(
        bad_exit_codes=[127, 160],
        log_level=1,
        threshold_for_failure=50,
        threshold_for_skipping=50,
    )
)