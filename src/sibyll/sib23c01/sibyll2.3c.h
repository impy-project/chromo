#ifndef _include_sib23c_interface_h_
#define _include_sib23c_interface_h_
//----------------------------------------------
//  C++ interface for the SIBYLL event generator
//----------------------------------------------
//wrapper
extern"C"{

  typedef char s_name[6];

  // SIBYLL particle stack (FORTRAN COMMON)
  // variables are: np : numer of particles on stack
  //                 p : 4momentum + mass of particles on stack
  //             llist : id of particles on stack
  extern struct { 
    double p[5][8000]; 
    int llist[8000]; 
    int np; 
  }s_plist_;


  // additional particle stack for the mother particles of unstable particles
  // stable particles have entry zero
  extern struct { 
    int llist1[8000]; 
  }s_plist1_;


  // tables with particle properties
  // charge, strangeness and baryon number
  extern struct {
    int ichp[99]; 
    int istr[99]; 
    int ibar[99]; 
  }s_chp_;

  // tables with particle properties
  // mass and mass squared
  extern struct {
    double am[99]; 
    double am2[99]; 
  }s_mass1_;

  // table with particle names
  extern struct {char namp[6][99];}s_cnam_;

  // debug info
  extern struct {int ncall; int ndebug; int lun; }s_debug_;

  // lund random generator setup
  //extern struct {int mrlu[6]; float rrlu[100]; }ludatr_;


  // sibyll main subroutine
  void sibyll_(int&,int&,double&);
    
  // subroutine to initiate sibyll
  void sibyll_ini_();

  // subroutine to SET DECAYS
  void dec_ini_();

  // subroutine to initiate random number generator
  void rnd_ini_();
  
  // print event
  void sib_list_(int&);

  // decay routine
  void decsib_();

  double s_rndm_(int&);

  // phojet random generator setup
  void pho_rndin_(int&,int&,int&,int&);
}
#endif

