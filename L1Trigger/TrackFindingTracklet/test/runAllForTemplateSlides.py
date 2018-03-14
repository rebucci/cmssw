import ROOT as r
import os.path

from optparse import OptionParser
parser = OptionParser()

# Use this for user specific label at the end of the filename
parser.add_option('--treeName', metavar='F', type='string', action='store',
                  default='',
                  dest='name',
                  help='')
parser.add_option('--runPGun', metavar='M', action='store_true',
                  default=False,
                  dest='runPGun',
                  help='')
parser.add_option('--run2GeV', metavar='M', action='store_true',
                  default=False,
                  dest='run2GeV',
                  help='')
(options, args) = parser.parse_args()
argv = []

r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L L1TrackNtuplePlot.C++");

treeName = ""
if options.name != "":
    treeName = "_"+options.name

truncationOptions = [
'WithTruncation',
'WithoutTruncation'
]

PUs = [
'PU0',
'PU140',
'PU200'
]

pGunSamples = {
	'MuonPt10' : 13,
	'MuonPt100' : 13,
	'ElectronPt10' : 11,
	'ElectronPt35' : 11,
}

run = r.gROOT.ProcessLine


print 'Running ttbar samples'
for PU in PUs:
	for truncation in truncationOptions:
		inputFile = 'TTbar_{PU}_{truncation}'.format( PU=PU, truncation = truncation)
		if not os.path.isfile( inputFile+'.root' ): 
			continue
		print '---> ',inputFile
                run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )
 		print '---------> Muons'
 		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,13,0,3)'.format(inputFile=inputFile, treeName=treeName) )
 		print '---------> Electrons'
 		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,11,0,3)'.format(inputFile=inputFile, treeName=treeName) )

		print 'In jet'
		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",1,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )
		print 'In high pt jet'
		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",2,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )
        
                if options.run2GeV:
                    print "----> 2 GeV threshold instead"
                    run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0)'.format(inputFile=inputFile, treeName=treeName) )
                    run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,13)'.format(inputFile=inputFile, treeName=treeName) )
                    run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,11)'.format(inputFile=inputFile, treeName=treeName) )
                    run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",1)'.format(inputFile=inputFile, treeName=treeName) )
                    run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",2)'.format(inputFile=inputFile, treeName=treeName) )

        print "Done"

if options.runPGun :
    print 'Running PGun samples'
    for sample, pdg in pGunSamples.iteritems():
        print '---> ',sample
        for PU in PUs:
            for truncation in truncationOptions:
                inputFile = '{sample}_{PU}_{truncation}'.format( sample=sample, PU=PU, truncation = truncation)
                if not os.path.isfile( inputFile+'.root' ): 
                    continue
                print '--------->',inputFile
                run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,{pdg},0,3)'.format(inputFile=inputFile,pdg=pdg, treeName=treeName) )

    print "Done"

