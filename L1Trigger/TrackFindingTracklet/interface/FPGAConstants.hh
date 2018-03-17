#ifndef FPGACONSTANTS_H
#define FPGACONSTANTS_H

//Uncomment if you want root output
//#define USEROOT

//Gemetry extensions
//static string geomext="D3";  //Use only detector region 3
//static string geomext="D4";  //Use only detector region 4
//static string geomext="D3D4";  //Use only detector region 3+4
//static string geomext="D5D6";  //Use forward disks
//static string geomext="full";  //Use full
static string geomext="new";  //Use full

static int TMUX = 6;

static string fitpatternfile="fitpattern.txt";

//If this string is non-empty we will write ascii file with
//processed events
static string skimfile="";
//static string skimfile="evlist_skim.txt";

//Debug options (should be false for 'normal' operation)
static bool dumppars=false;
static bool dumpproj=false;
static bool dumpmatch=false;


static bool writeInvTable=false; //Write out tables of drinv in tracklet calculator

static bool writeFitDerTable=false; //Write out track derivative tables


static bool writeDTCLinks=false;
static bool writeStubsLayer=false;
static bool writeStubsLayerperSector=false;
static bool writeLayerRouter=false;  //obsolete
static bool writeDiskRouter=false;   //obsolete
static bool writeAllStubs=false;
static bool writeVMOccupancyME=false;
static bool writeVMOccupancyTE=false;
static bool writeTE=false;
static bool writeTrackletCalculator=false;
static bool writeTrackletPars=false;
static bool writeNeighborProj=false;
static bool writeAllProjections=false;
static bool writeVMProjections=false;
static bool writeTrackProj=false;
static bool writeTrackProjOcc=false;
static bool writeME=false;
static bool writeMatchCalculator=false;
static bool writeResiduals=false;
static bool writeProjectionTransceiver=false;
static bool writeMatchTransceiver=false;
static bool writeFitTrack=false;
static bool writeChiSq=false;
static bool writeNMatches=false;
static bool writez0andrinv=false;
static bool writeDiskMatch1=false;
static bool writeHitEff=false;

static bool writeTETables=false;
static bool writeVMTables=false;


static bool writeHitPattern=false;
static bool writeTrackletParsOverlap=false;
static bool writeTrackletParsDisk=false;

static bool writeAllCT=false; //write out .dat file containing all output tracks in bitwise format

static bool writeResEff=false; //write files for making resolution & efficiency plots for standable code version
static bool writePars=false; //write files for making plots of track parameters



static bool writestubs=false;
static bool writestubs_in2=false;
static bool writeifit=false;
static bool padding=true;

static bool exactderivatives=false;  //for both the integer and float
static bool exactderivativesforfloating=true; //only for the floating point
static bool useapprox=true; //use approximate postion based on integer representation for floating point
static double alphamax=5.0/(65.0*65.0);
static int nbitsalpha=6;
static double kalpha=alphamax/(1<<(nbitsalpha-1));
static int alphaBitsTable=2; //For number of bits in track derivative table
static int nrinvBitsTable=3; //number of bits for tabulating rinv dependence

static bool writetrace=false; //Print out details about startup
static bool debug1=false; //Print information about tracking
static bool writeoutReal = false; 
static bool writemem=false; //Note that for 'full' detector this will open
                            //a LOT of files, and the program will run excruciatingly slow
static unsigned int writememsect=3;  //writemem only for this sector

static bool warnNoMem=false;  //If true will print out warnings about missing projection memories

//Program flow (should be true for normal operation)
//enables the stub finding in these layer/disk combinations
static bool doL1L2=true;
static bool doL3L4=true;
static bool doL5L6=true;

static bool doD1D2=true; 
static bool doD3D4=true;

static bool doL1D1=true;
static bool doL2D1=true;

static bool doMEMatch=true;

static bool allSector=false; //if true duplicate stubs in all sectors

static bool doProjections=true;

static const int MEBinsBits=3;
static const int MEBins=(1<<MEBinsBits);

//Geometry
static double zlength=120.0;

// can automatically determine the above values using script plotstub.cc:
// root -b -q 'plotstub.cc("evlist_MuPlus_1to10_D11_PU0")'

// these assume D11 geometry!
static double rmeanL1=25.1493;
static double rmeanL2=37.468;
static double rmeanL3=52.5977;
static double rmeanL4=68.7737;
static double rmeanL5=86.0591;
static double rmeanL6=110.844;

static double zmeanD1=131.18;
static double zmeanD2=155.0;
static double zmeanD3=185.34;
static double zmeanD4=221.619;
static double zmeanD5=265.0;


//now the half-way between rmindisk and rmaxdisk corresponds to PS/SS boundary 
static double rmindisk=0.0;
static double rmaxdisk=120.0;

static double rmindiskvm=27.0;
static double rmaxdiskvm=67.0;

static double rmaxdiskl1overlapvm=45.0;
static double rmindiskl2overlapvm=40.0;


// need separate lookup values for inner two vs outer three disks for 2S modules
// these assume D11 geometry!
static double rDSS[20] = {66.7728, 71.7967, 77.5409, 82.5584, 84.8736, 89.8953, 95.7791, 100.798, 102.495, 107.52,   // <=== these 10 are for inner 2 disks
			  65.1694, 70.1936, 75.6641, 80.6908, 83.9581, 88.9827, 94.6539, 99.6772, 102.494, 107.519}; // <=== these 10 are for outer 3 disks


static double drmax=(rmaxdisk-rmindisk)/32.0;

static double dzmax=zlength/32.0;

static double drdisk=rmaxdisk-rmindisk;

static double rmean[6]={rmeanL1,rmeanL2,rmeanL3,rmeanL4,rmeanL5,rmeanL6};

static double zmean[5]={zmeanD1,zmeanD2,zmeanD3,zmeanD4,zmeanD5};



static double rminL1=rmeanL1-drmax; 
static double rmaxL1=rmeanL1+drmax; 
static double rminL2=rmeanL2-drmax; 
static double rmaxL2=rmeanL2+drmax; 
static double rminL3=rmeanL3-drmax; 
static double rmaxL3=rmeanL3+drmax; 
static double rminL4=rmeanL4-drmax; 
static double rmaxL4=rmeanL4+drmax; 
static double rminL5=rmeanL5-drmax; 
static double rmaxL5=rmeanL5+drmax; 
static double rminL6=rmeanL6-drmax; 
static double rmaxL6=rmeanL6+drmax; 

static double zminD1=zmeanD1-dzmax; 
static double zmaxD1=zmeanD1+dzmax; 
static double zminD2=zmeanD2-dzmax; 
static double zmaxD2=zmeanD2+dzmax; 
static double zminD3=zmeanD3-dzmax; 
static double zmaxD3=zmeanD3+dzmax; 
static double zminD4=zmeanD4-dzmax; 
static double zmaxD4=zmeanD4+dzmax; 
static double zminD5=zmeanD5-dzmax; 
static double zmaxD5=zmeanD5+dzmax; 


static bool   enstubbend = false; 
static double two_pi=8.0*atan(1.0);

static double ptcut=1.91; //Minimum pt
static double rinvcut=0.01*0.3*3.8/ptcut; //0.01 to convert to cm-1
static double ptcutte=1.6; //Minimum pt in TE
static double rinvcutte=0.01*0.3*3.8/ptcutte; //0.01 to convert to cm-1 in TE
static double bendcut=1.5;
static double bendcutdisk=3.0;
static double z0cut=15.0;


static unsigned int NSector=27; 
static int Nphibits=2;         //Number of bits required to label the phi VM
static int L1Nphi=(1<<Nphibits)-1; //Number of odd layer VMs
static int Nzbits=3;         //Number of bits required to label the z VM
static int L1Nz=(1<<Nzbits); //Number of z VMs in odd layers
static int VMzbits=4;        //Number of bits for the z position in VM
static int L2Nphi=(1<<Nphibits); //Number of even layer VMs
static int L2Nz=(1<<Nzbits); //Number of z VMs in even layers
static int VMrbits=2;        //Number of bits for r position 'in VM'
static int VMphibits=3;      //Number of bits for phi position in VM


static const unsigned int NLONGVMBITS=3; 
static const unsigned int NLONGVMRBITS=3;   //4 bins on each side (+-z)
static const unsigned int NLONGVMBINS=(1<<NLONGVMBITS);
static const unsigned int NLONGVMRBINS=(1<<NLONGVMRBITS);
static const unsigned int NLONGVMODDLAYERBITS=6;
static const unsigned int NLONGVMODDDISKBITS=6;
static const double rinnerdisk=22.0;
static const double routerPSdisk=65.0;



//limits per FED region
//static int NMAXstub  = 250;
//static int NMAXroute = 250;

static unsigned int MAXOFFSET=10000; //set to 0 for regular truncation

static unsigned int MAXSTUBSLINK = 36 + MAXOFFSET; //Max stubs per link
static unsigned int MAXLAYERROUTER = 36 + MAXOFFSET; //Max stubs handled by layer router
static unsigned int MAXDISKROUTER = 36 + MAXOFFSET; //Max stubs handled by disk router
static unsigned int MAXVMROUTER = 36 + MAXOFFSET; //Max stubs handled by VM router
static unsigned int MAXTE = 36 + MAXOFFSET; //Maximum number of stub pairs to try in TE 
static unsigned int MAXTC = 36 + MAXOFFSET; //Maximum number of tracklet parameter calculations
//static unsigned int MAXPROJECTIONTRANSCEIVER = 10000; //Maximum number of projections to neighbor
static unsigned int MAXPROJROUTER = 36 + MAXOFFSET; //Maximum number of projections to route
static unsigned int MAXME = 36 + MAXOFFSET; //Maximum number of stub-projection matches to try
static unsigned int MAXMC = 36 + MAXOFFSET; //Maximum number of match calculations
static unsigned int MAXFIT = 36 + MAXOFFSET; //Maximum number of track fits


static double dphisector=two_pi/NSector;

//Constants for defining stub representations
static int nbitsrL123=7;
static int nbitsrL456=8;

static int nbitszL123=12;
static int nbitszL456=8;

static int nbitsphistubL123=14;
static int nbitsphistubL456=17;

static int nrbitsdisk=12;
static int nzbitsdisk=7;

static int nrbitsprojdisk=12;
static int nrbitsprojderdisk=9;


static int nrbitsdiskvm=5;
static int nzbitsdiskvm=2;

static int Nrbitsdisk=2;         //Number of bits required to label the r VM


static int nbitsphiprojL123=nbitsphistubL123;
static int nbitsphiprojL456=nbitsphistubL456;

static int nbitszprojL123=12;
static int nbitszprojL456=8;

static int nbitsphiprojderL123=7;
static int nbitsphiprojderL456=8;

static int nbitszprojderL123=8+1;
static int nbitszprojderL456=7+2;






//Bits used to store track parameter in tracklet
static int nbitsrinv=14;
static int nbitsphi0=18;
static int nbitst=14;
static int nbitsz0=10;

//Minimal ranges for track parameters
static double maxrinv=0.006;
static double maxphi0=0.59;
static double maxt=9.0;
static double maxz0=28.0;

static double rmin[6]={rminL1,rminL2,rminL3,rminL4,rminL5,rminL6};

//These are constants for track paramter calculation. Some will be 
//used in the actual calculation. Others are just needed to compare 
//the results to the floating point calculations

static int idrinvbits=19;
static int it1shift=23;
static int it2shift=15;
static int it3shift=13;
static int it4bits=9;
static int it5bits=9;
static int irinvshift=10;
static int it7shift=10;
static int it7tmpbits=39;
static int it7tmpshift=20;
static int it7shift2=5;
static int it9bits=12;
static int itshift=5;
static int it12shift=8;

static int it1shiftdisk=23;
static int it2shiftdisk=15;
static int it3shiftdisk=13;
static int it4bitsdisk=9;
static int it5bitsdisk=9;
static int irinvshiftdisk=10;
static int it7shiftdisk=10;
static int it7tmpbitsdisk=39;
static int it7tmpshiftdisk=20;
static int it7shift2disk=5;
static int it9bitsdisk=12;
static int itshiftdisk=5;
static int it12shiftdisk=8;


static double kphi=two_pi/((0.75*NSector)*(1<<nbitsphistubL123));
static double kphi1=two_pi/((0.75*NSector)*(1<<nbitsphistubL456));
static double kz=2*zlength/(1<<nbitszL123);
static double kr=2*drmax/(1<<nbitsrL456);			
static double kdrinv=1.0/(kr*(1<<idrinvbits));
static double kdelta=(1<<(it1shift+it3shift-2*(idrinvbits-it2shift)))*kphi1*kphi1;	
static double kt5=1.0/(1<<it5bits);	
static double krinv=kphi1/(kr*(1<<idrinvbits));
static double kt=kz*kdrinv;
static double kt7=kphi1*(1<<(irinvshift+it7shift-idrinvbits));
static int it7tmpfactor=(1<<(it7tmpbits-idrinvbits+irinvshift))*kphi1/sqrt(6.0);

//static double kt9=1.0/(1<<(2*it7tmpbits-2*it7tmpshift-2*it7shift));
static double kt9=1.0/(1<<it9bits);
static double kt12=kz/(1<<(idrinvbits-itshift-it12shift));


static double kzdisk=2.0*dzmax/(1<<nzbitsdisk);
static double krdisk=drdisk/(1<<nrbitsdisk);
static double kdeltadisk=(1<<(it1shiftdisk+it3shiftdisk-2*(idrinvbits-it2shiftdisk)))*kphi1*kphi1;	
static double kdrinvdisk=1.0/(krdisk*(1<<idrinvbits));
static double ktdisk=kzdisk*kdrinvdisk;
static double kt12disk=kz/(1<<(idrinvbits-itshiftdisk-it12shiftdisk));
static double krinvdisk=kphi1/(krdisk*(1<<idrinvbits));
static double kt7disk=kphi1/(1<<(irinvshiftdisk+it7shiftdisk-idrinvbits));
static int it7tmpfactordisk=(1<<(it7tmpbitsdisk-idrinvbits+irinvshiftdisk))*kphi1/sqrt(6.0);
static double kt9disk=1.0/(1<<(2*it7tmpbitsdisk-2*it7tmpshiftdisk-2*it7shiftdisk)); 


static int rinvbitshift=13; //(int)(1.0+log((maxrinv/(1<<(nbitsrinv-1)))/krinv)/log(2.0));
static int phi0bitshift=1; //(int)(1.0+log((maxphi0/(1<<(nbitsphi0-1)))/kphi1)/log(2.0));
static int tbitshift=9; //(int)(1.0+log((maxt/(1<<(nbitst-1)))/kt)/log(2.0));
static int z0bitshift=0; //(int)(1.0+log((maxz0/(1<<(nbitsz0-1)))/kz)/log(2.0));


static int rinvbitshiftdisk=13; //(int)(1.0+log((maxrinv/(1<<(nbitsrinv-1)))/krinvdisk)/log(2.0));
static int phi0bitshiftdisk=1; //(int)(1.0+log((maxphi0/(1<<(nbitsphi0-1)))/kphi1)/log(2.0));


static double krinvpars=krinv*(1<<rinvbitshift);
static double kphi0pars=kphi1*(1<<phi0bitshift);
static double kphi0parsdisk=kphi1*(1<<phi0bitshiftdisk);
static double ktpars=kt*(1<<tbitshift);
static double kzpars=kz;

static int t2bits=23;
static int t3shift=8;
static int t4shift=8;
static int t4shift2=8;
static int t6bits=12;

static double krinvparsdisk=krinvdisk*(1<<rinvbitshiftdisk);
static double ktparsdisk=ktdisk*(1<<tbitshift);
static double kt2disk=(1.0/ktparsdisk)/(1<<t2bits);
static double kt3disk=kt2disk*kzdisk*(1<<t3shift);
static double kt4disk=kt3disk*krinvparsdisk*(1<<t4shift);
static double kt5disk=1.0/(1<<it5bitsdisk);	
static double kt6disk=1.0/(1<<t6bits);
static double krprojdisk=kt3disk*kt6disk*(1<<t6bits);
static double krprojderdisk=kt2disk;
static double kphiprojderdisk=kt2disk*krinvparsdisk;

static double kst5disk=kt4disk*kt4disk*(1<<(2*t4shift2));

static int rprojdiskbitshift=6;
static int phiderdiskbitshift=20;
static int rderdiskbitshift=7;
static double krprojshiftdisk=krprojdisk*(1<<rprojdiskbitshift);

static double kphiprojdisk=kphi0parsdisk*4.0;
static double kphiprojderdiskshift=kphiprojderdisk*(1<<phiderdiskbitshift);
static double krprojderdiskshift=krprojderdisk*(1<<rderdiskbitshift);


//Parameters for projections

static int is1shift=9;
static int is2shift=23;
static int is3bits=10;
static int is5shift=9;

static double ks1=kr*krinvpars*(1<<is1shift);
static double ks2=ks1*ks1*(1<<is2shift);
static double ks3=1.0/(1<<is3bits);
static double ks4=ks1*ks3*(1<<(idrinvbits+phi0bitshift-is1shift+is3bits-rinvbitshift));
static double ks5=kt*kr*(1<<tbitshift)*(1<<is5shift);
static double ks6=ks5*ks3;

static int phiderbitshift=7;
//static int zderbitshift=6-3;
static int zderbitshift=6; //changed to tracklet 1.0 version to run on ttbar PU=200


static double kphiproj123=kphi0pars*4;
static double kphiproj456=kphi0pars/2;
static double kzproj=kz;
static double kphider=krinvpars*(1<<phiderbitshift);
static double kzder=ktpars*(1<<zderbitshift);

static int phiresidbits=12; 
static int zresidbits=9;
static int rresidbits=7;
//static int rresidbits=8; //tracklet 1.0 version


//Trackfit
static int fitrinvbitshift=9;  //6 OK?
static int fitphi0bitshift=6;  //4 OK?
static int fittbitshift=10;     //4 OK? //lower number gives rounding problems
static int fitz0bitshift=8;    //6 OK?

static int chisqphifactbits=14;
static int chisqzfactbits=14;

//Duplicate Removal
static int minIndStubs=3;
static bool AdjacentRemoval=true;
static string RemovalType="ichi";
//"ichi" (pairwise, keep track with best ichisq), "nstub" (pairwise, keep track with more stubs), "grid" (TMTT-like removal)

#endif




