#ifndef FERMI_ROT_H
#define FERMI_ROT_H

#include <iomanip>
#include <string>
#include <stdio.h>
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include "global.h"

void fermi_rot(Float_t &E_beam_fermi, Float_t &theta_rot2, Float_t E_beam,TLorentzVector P4_E_prime,TLorentzVector  &P4_E_prime_boosted);
  
#endif
