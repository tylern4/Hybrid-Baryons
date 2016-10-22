#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"


using namespace std;

//----------------------------------------------------
//WE ASSUMED THAT IN REACTION A + B = 1 + 2 + 3
//A - VIRTUAL PHOTON
//B - INITIAL PROTON
//1 - PARTICLE WITH MASS m1
//2 - PARTICLE WITH MASS m2
//3 - PARTICLE WITH MASS m3

//We assume as global variables: MP - initial proton mass (is taken from pdg-table), E_beam - beam energy (is taken from input file).

//Here all derivations are performed via special root functions. If you need further level of explanation please see anti_rot_explanation, where all calculations are shown in more details and rotations are performed via matrices.

 void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2,Float_t &phi_rot2,Float_t E_beam,TLorentzVector  P4_E_prime, TLorentzVector  &P4_E_prime_boosted) {
 
 Float_t theta_fermi, phi_fermi;
 TLorentzVector P4_EL, P4_in_Prot, P4_gamma;
 
 P4_EL.SetXYZT(0.,0.,E_beam,E_beam);
 P4_in_Prot.SetXYZT(px_fermi,py_fermi,pz_fermi,sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 
 P4_gamma = P4_EL -  P4_E_prime;
 
  Float_t Q2 = -P4_gamma.Mag2();
//  cout <<   P4_gamma[0] << " "<< P4_gamma[1] << " "<< P4_gamma[2] << " "<< P4_gamma[3] << " q2\n";
//cout << P4_EL[0] << " "<<  P4_EL[1] << " "<<  P4_EL[2] << " "<< P4_EL[3] << "\n";


// cout << P4_E_prime[0] << " "<<  P4_E_prime[1] << " "<<  P4_E_prime[2] << " "<< P4_E_prime[3] << "   qqq6\n";
 
// cout << P4_gamma_test[3]<< " qwqw\n";

 theta_fermi = acos(pz_fermi/sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi));
 phi_fermi = acos(fabs(px_fermi)/sqrt(px_fermi*px_fermi+py_fermi*py_fermi));
 
 if ((px_fermi < 0.)&&(py_fermi > 0.)) phi_fermi = M_PI-phi_fermi;
 if ((px_fermi < 0.)&&(py_fermi < 0.)) phi_fermi = phi_fermi + M_PI;
 if ((px_fermi > 0.)&&(py_fermi < 0.)) phi_fermi = 2.*M_PI - phi_fermi;
 
Float_t eps_t1 = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
Float_t eps_l1 = Q2*eps_t1/P4_gamma[3]/P4_gamma[3];

  P4_EL.RotateZ(-phi_fermi);
  P4_EL.RotateY(-theta_fermi);
  
  P4_gamma.RotateZ(-phi_fermi);
  P4_gamma.RotateY(-theta_fermi);
   
  P4_E_prime.RotateZ(-phi_fermi);
  P4_E_prime.RotateY(-theta_fermi);
   
  P4_in_Prot.RotateZ(-phi_fermi);
  P4_in_Prot.RotateY(-theta_fermi);
  

//cout << P4_EL[0] << "   "<<  P4_EL[1] << "   "<<  P4_EL[2] << "   "<< P4_EL[3] << " ioi\n";
//cout << R[0] << "   "<<  R[1] << "   "<< R[2] <<  "\n";

  
// cout << "rrr   "<< P4_in_Prot[0] << "   "<<  P4_in_Prot[1] << "   "<<  P4_in_Prot[2] << "   "<< P4_in_Prot[3]<<"   " <<sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi) << "\n";
 

//------------------------------------
 
 Float_t beta;
 
 beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/sqrt(MP*MP+px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi);
// beta = sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi)/MP;
 P4_in_Prot.Boost(0,0,-beta);
 P4_EL.Boost(0,0,-beta);
 P4_E_prime.Boost(0,0,-beta);
 P4_gamma.Boost(0,0,-beta);
 
//cout << "rrr   "<< P4_in_Prot[0] << "   "<<  P4_in_Prot[1] << "   "<<  P4_in_Prot[2] << "   "<< P4_in_Prot[3]<<"   " <<sqrt(px_fermi*px_fermi+py_fermi*py_fermi+pz_fermi*pz_fermi) << "\n";

//-----------------------

//ROT2


  theta_rot2 = acos(P4_EL[2]/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));


  phi_rot2 = acos(fabs(P4_EL[0])/sqrt(P4_EL[0]*P4_EL[0]+P4_EL[1]*P4_EL[1]+P4_EL[2]*P4_EL[2]));
 
 if ((P4_EL[0] < 0.)&&(P4_EL[1] > 0.)) phi_rot2 = M_PI-phi_rot2;
 if ((P4_EL[0] < 0.)&&(P4_EL[1] < 0.)) phi_rot2 = phi_rot2 + M_PI;
 if ((P4_EL[0] > 0.)&&(P4_EL[1] < 0.)) phi_rot2 = 2.*M_PI - phi_rot2;


 P4_EL.RotateY(theta_rot2);
 P4_EL.RotateZ(phi_rot2);
 

 P4_gamma.RotateY(theta_rot2);
 P4_gamma.RotateZ(phi_rot2);
 
 E_beam_fermi = P4_EL[3];

//cout <<"rrrtrrrr                  "<< P4_EL[0] << "   "<<  P4_EL[1] << "   "<<  P4_EL[2] << "   "<< P4_EL[3] << "\n"; 
 
 P4_E_prime.RotateY(theta_rot2);
 P4_E_prime.RotateZ(phi_rot2);
 
  P4_in_Prot.RotateY(theta_rot2);
  P4_in_Prot.RotateZ(phi_rot2);
 //----------------------------
 P4_E_prime_boosted=P4_E_prime;
 
Float_t eps_t2 = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_EL.Vect()).Cross(P4_E_prime.Vect())).Mag2()));
Float_t eps_l2 = Q2*eps_t2/P4_gamma[3]/P4_gamma[3];

//cout << eps_t1<< " "<< eps_t2<<" "<< eps_t1/eps_t2<<"  \n";
//cout << eps_l2<<"  f\n";

 /*
  Float_t W_f = (P4_in_Prot + P4_EL - P4_E_prime).Mag();
  Float_t Q2 = -P4_gamma.Mag2();
 
// cout << (W_f*W_f+Q2-MP*MP)/2./MP <<" "<< P4_gamma[3]<<"  \n";
TRotation rot;
TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_EL.Vect().Cross(P4_E_prime.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();

 P4_gamma.Transform(rot);


 
P4_gamma.Boost(0,0,-sqrt(P4_gamma[3]*P4_gamma[3]+Q2)/(P4_gamma[3]+MP));
*/
//cout << (W_f*W_f-Q2-MP*MP)/2./W_f << " "<< P4_gamma[3]<< "  cms\n";
 };
