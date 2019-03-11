#ifndef ROT_H
#define ROT_H
#include <iomanip>
#include <string>
#include <stdio.h>
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

void rot(Float_t Q2, Float_t E_beam, TLorentzVector P4_E_fin, TLorentzVector P4_1, TLorentzVector P4_2, TLorentzVector P4_3, Float_t &M12, Float_t &M23, Float_t &theta_hadr,Float_t &alpha_hadr, Float_t &phi_hadr);

  
#endif
