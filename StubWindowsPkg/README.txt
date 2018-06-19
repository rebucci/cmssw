~~ Stub Windows Analysis Package ~~

EXTRACTOR 
- Produce usable root file using the extractor
		cd DataProduction/test
	cmsRun SLHC_Extr.py
- Run macros in RecoExtractor/test over extracted files

ANALYSIS
1. In StubAnalysis directory, first set up the tools.
	- Unzip tar file
		tar -zxf tclap-1.2.1.tar.gz
	- Make AMA_ana
		cd AMA_ana
		make
	- Make FE_ana
		cd ../FE_ana
		make
2. Get the stub rates file from the AMA_ana
		./AM_ana -c rates -i INPUTFILE -o OUTPUTFILE
3. Get the efficiency file from AM_ana
		./AM_ana -c stub_eff -i INPUTFILE -o OUTPUTFILE -t PDGID
4. Use PlotStubRates.C and PlotStubEffs.C macros to make plots from the output files produced by AM_ana. Instructions are at the top of each macro.
