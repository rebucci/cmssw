#include "jobparams.h"

//tclap
#include <tclap/CmdLine.h>
using namespace TCLAP;

jobparams::jobparams(int argc, char** argv){
  


   try {
     // command line parser
     CmdLine cmd("Command option", ' ', "0.9");

     ValueArg<std::string> option("c","case","type of analysis (rates/sectors/rate_n_sec/sec_n_test/stub_eff/PR_eff)",
				false, "rates", "string");
     cmd.add(option);

     ValueArg<std::string> pattfile("d","data","name of the pattern data file",
				   false, "", "string");
     cmd.add(pattfile);

     ValueArg<std::string> testfile("f","fortest","name of the input file used to test",
				    false, "test.root", "string");
     cmd.add(testfile);

     ValueArg<bool> dbg("g","gdb","debug mode (standalone PR used) or not",
			false, 0, "bool");
     cmd.add(dbg);

     ValueArg<std::string> inputfile("i","input","path and name of the input file",
				false, "/scratch/viret/data.root", "string");
     cmd.add(inputfile);
       
     ValueArg<std::string> inputfileQCD("j","jinput","path and name of the input QCD file",
                                        false, "/scratch/viret/data.root", "string");
     cmd.add(inputfileQCD);
       
     ValueArg<int> nevt("n","nevt","number of events for the eff test?",
			false, 0, "int");
     cmd.add(nevt);

     ValueArg<std::string> outfile("o","output","name of the output file",
				false, "output.root", "string");
     cmd.add(outfile);
         
     ValueArg<float> ptmin("p","ptmin","minimal pt value for primary selection (in GeV/c)",
			   false, 2., "float");
     cmd.add(ptmin);
       
     ValueArg<float> sptminmu("y","sptminmu","minimal pt value for muon for SW tuning (in GeV/c)",
                            false, 2., "float");
     cmd.add(sptminmu);
  
     ValueArg<float> sptmaxmu("z","sptmaxmu","maximal pt value for muon for SW tuning (in GeV/c)",
                            false, 6., "float");
     cmd.add(sptmaxmu);
       
     ValueArg<float> sptminele("k","sptminele","minimal pt value for electron for SW tuning (in GeV/c)",
                            false, 5., "float");
     cmd.add(sptminele);
       
     ValueArg<float> sptmaxele("l","sptmaxele","maximal pt value for electron for SW tuning (in GeV/c)",
                            false, 20., "float");
     cmd.add(sptmaxele);
       
     ValueArg<float> maxlosses("m","maxlosses","max proportion of CIC segments with good stub losses accepted (in %)",
			   false, 1., "float");
     cmd.add(maxlosses);

 
     ValueArg<float> qmax("q","qmax","maximal d0 value for primary selection (in cm)",
			  false, 0.1, "float");
     cmd.add(qmax);

     ValueArg<int> type("t","type","PDG id of the particle to analyze",
			false, 13, "int");
     cmd.add(type);

     ValueArg<int> rate("r","rate","L1 rate (for pattern gen, in kHz)",
			false, 100, "int");
     cmd.add(rate);

     // parse
     cmd.parse(argc, argv);
     
     m_inputfile    = inputfile.getValue();
     m_inputfileQCD = inputfileQCD.getValue();
     m_outfile      = outfile.getValue();
     m_testfile     = testfile.getValue();
     m_pattfile     = pattfile.getValue();
     m_opt          = option.getValue();
     m_nevt         = nevt.getValue();
     m_dbg          = dbg.getValue();
     m_type         = type.getValue();
     m_rmax         = qmax.getValue();
     m_ptmin        = ptmin.getValue();
     m_ptminmu      = sptminmu.getValue();
     m_ptmaxmu      = sptmaxmu.getValue();
     m_ptminele     = sptminele.getValue();
     m_ptmaxele     = sptmaxele.getValue();
     m_maxlosses    = maxlosses.getValue();
     m_rate         = rate.getValue();
   }
   catch (ArgException &e){ // catch exception from parse
     std::cerr << "ERROR: " << e.error() << " for arg " << e.argId()  << std::endl;
     abort();
   }
}
