#include<iostream>
#include<cmath>

#include"sibyll2.3c.h"

int main(int argc, char *argv[]){

  // redirect table output to file
  s_debug_.lun = 7;

  rnd_ini_();
  
  std::cout << "sequence of random numbers."  << std::endl;
  int a = 0;
  for(int i=0; i<5; ++i)
    std::cout << i << " " << s_rndm_(a) << std::endl;

  // reset generator
  // set seed for random number generator to different sequence
  int na,nb,nc,nd;
  na = 15; // default would be 12
  nb = 34;
  nc = 56;
  nd = 78;
  pho_rndin_(na,nb,nc,nd);

  std::cout << "another sequence of random numbers."  << std::endl;
  for(int i=0; i<5; ++i)
    std::cout << i << " " << s_rndm_(a) << std::endl;

  // reset to initial sequence
  rnd_ini_();
  std::cout << "again first sequence of random numbers. seed:" << std::endl;

  for(int i=0; i<5; ++i)
    std::cout << i << " " << s_rndm_(a) << std::endl;

  sibyll_ini_();

  dec_ini_();
 
  int ibeam = 13;

  int itarget = 1;

  double ecm = 7000;

  int Nevt = 1;

  // links to particle stack
  int *id = s_plist_.llist;
  double *px = s_plist_.p[0];
  double *py = s_plist_.p[1];
  double *pz = s_plist_.p[2];
  double *en = s_plist_.p[3];
  double *mass = s_plist_.p[4];

  int pid;
  int iout = 6;
  int particleCharge;

  for( int ievt=0; ievt < Nevt; ++ievt ){

    sibyll_( ibeam, itarget , ecm );
    decsib_();

    sib_list_(iout);

    std::cout << "number of particles: " << s_plist_.np << std::endl;
    // loop over final state particles
    for( int ipart= 0; ipart < s_plist_.np; ++ipart ){

      // select particle properties (charged)
      // fortran indexing: 1 .. upper_bound !!!
      // c++ indexing: 0 .. upper_bound - 1 
      // need to shift id when passing to array of particle properties

      // get rid of unstable flag
      pid = id[ ipart ] % 10000;

      //std::cout << ipart << " " << pid << std::endl;
      particleCharge = s_chp_.ichp[ int(std::abs( pid )) - 1 ];
      if( 0 != particleCharge ){
    	std::cout << "(id,px,chrg): "  << pid << " "
		  << px[ ipart ] << " "
		  << particleCharge << std::endl;
      }
      
    }
    
  }

}
