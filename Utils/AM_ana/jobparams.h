#ifndef jobparams_H
#define jobparams_H

#include <string>
#include <cstdio>
#include <tclap/CmdLine.h>
using namespace TCLAP;

class jobparams{

 public:

  /** constructer reading from standard arg */
  jobparams(int argc, char** argv);

  /** default constructer */
  jobparams(){}

  /** copy constructer */
  jobparams(const jobparams& tparams); 

  /** destructer */
  ~jobparams(){}

  /** return value */

  bool        dbg() const;
  std::string option() const;
  std::string inputfile() const;
  std::string inputfileQCD() const;
  std::string testfile() const;
  std::string outfile() const;
  std::string pattfile() const;
  float       qmax() const;
  float       ptmin() const;
  float       sptminmu() const;
  float       sptmaxmu() const;
  float       sptminele() const;
  float       sptmaxele() const;
  float       maxlosses() const;
  int         nevt() const;
  int         type() const;
  int         rate() const;

 private:

  bool         m_dbg;   
  std::string  m_opt;  
  std::string  m_inputfile;
  std::string  m_inputfileQCD;
  std::string  m_outfile;
  std::string  m_testfile; 
  std::string  m_pattfile;  
  float        m_rmax;
  float        m_ptmin;
  float        m_ptminmu;
  float        m_ptmaxmu;
  float        m_ptminele;
  float        m_ptmaxele;
  float        m_maxlosses;
  int          m_nevt;
  int          m_type;
  int          m_rate;
};

inline std::string jobparams::option() const{
  return m_opt;
}

inline std::string jobparams::inputfile() const{
  return m_inputfile;
}

inline std::string jobparams::inputfileQCD() const{
    return m_inputfileQCD;
}

inline std::string jobparams::outfile() const{
  return m_outfile;
}

inline std::string jobparams::testfile() const{
  return m_testfile;
}

inline std::string jobparams::pattfile() const{
  return m_pattfile;
}

inline bool jobparams::dbg() const{
  return m_dbg;
}

inline float jobparams::qmax() const{
  return m_rmax;
}

inline float jobparams::ptmin() const{
  return m_ptmin;
}

inline float jobparams::sptminmu() const{
    return m_ptminmu;
}

inline float jobparams::sptmaxmu() const{
    return m_ptmaxmu;
}

inline float jobparams::sptminele() const{
    return m_ptminele;
}

inline float jobparams::sptmaxele() const{
    return m_ptmaxele;
}

inline float jobparams::maxlosses() const{
  return m_maxlosses;
}

inline int jobparams::nevt() const{
  return m_nevt;
}

inline int jobparams::type() const{
  return m_type;
}

inline int jobparams::rate() const{
  return m_rate;
}

#endif
