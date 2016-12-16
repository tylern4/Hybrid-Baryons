#include <iomanip>
#include <string>
#include <stdio.h> 
#include <stdlib.h>  
#include <math.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include "global.h"
#include <sstream>
#include "hist_def.h"
#include "interpol_int.h"




void hist_write(){
 
   //Dividing by Delta cos to obtain d(-cos theta) distributions for 1.225 < W < 2.15 GeV.
  for (Int_t ik=0;ik<=36;ik++){
  for (Int_t ij=1;ij<=100;ij++){
h_odn_theta_2[ik] ->SetBinContent (ij,      h_odn_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };
  
  //Dividing by Delta cos to obtain d(-cos theta) distributions for 2.1625 < W < 3.0375 GeV.  
  for (Int_t ik=0;ik<18;ik++){
  for (Int_t ij=1;ij<=100;ij++){
h_odn_wwide_theta_2[ik] ->SetBinContent (ij, h_odn_wwide_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };
  
  //Dividing by Delta cos to obtain d(-cos theta) distributions for 3.0875 < W < 4.5375 GeV.  
  for (Int_t ik=0;ik<15;ik++){
  for (Int_t ij=1;ij<=100;ij++){
h_odn_wgt3_theta_2[ik] ->SetBinContent (ij, h_odn_wgt3_theta[ik]->GetBinContent(ij)/(cos((h_odn_theta[ik]->GetBinLowEdge(ij)))-cos((h_odn_theta[ik]->GetBinLowEdge(ij)+h_odn_theta[ik]->GetBinWidth(ij)))));
     };
  };    
     
     
     
TFile *file1 = TFile::Open("out_hist_test.root","RECREATE");
file1->cd();
 file1->mkdir("1diff_w_1225_215");
 file1->mkdir("1diff_w_21625_30375");
 file1->mkdir("1diff_w_30875_45375");
 file1->mkdir("int_w_for_q2_01_13");
 file1->mkdir("int_q2_for_w_125_2075");
 file1->mkdir("int_cr_sect_norm");
  
h_zel->SetMinimum(0.);
h_phi_e->SetMinimum(0.);
h_Q2vsW->SetMinimum(0.); 
h_nu->SetMinimum(0.);
h_dalitz->SetMinimum(0.);


h_W->Write("", TObject::kOverwrite);
h_Q2->Write("", TObject::kOverwrite);
h_phi_e->Write("", TObject::kOverwrite);
h_zel->Write("", TObject::kOverwrite);


h_Q2vsW->Write("", TObject::kOverwrite);
h_Q2vsW2->Write("", TObject::kOverwrite);

h_nu->Write("", TObject::kOverwrite);
h_dalitz->Write("", TObject::kOverwrite);

h_0_miss->Write("", TObject::kOverwrite);
h_0_miss_2->Write("", TObject::kOverwrite);
h_pim_miss->Write("", TObject::kOverwrite);
h_pim_miss_2->Write("", TObject::kOverwrite);
h_0_miss_en->Write("", TObject::kOverwrite);
h_0_miss_en_2->Write("", TObject::kOverwrite);
 

h_0_miss_fermi->Write("", TObject::kOverwrite);
h_0_miss_fermi_2->Write("", TObject::kOverwrite);
h_pim_miss_fermi->Write("", TObject::kOverwrite);
h_pim_miss_fermi_2->Write("", TObject::kOverwrite);
h_0_miss_en_fermi->Write("", TObject::kOverwrite);
h_0_miss_en_fermi_2->Write("", TObject::kOverwrite);

h_eradgam->Write("", TObject::kOverwrite);
h_fermi_bonn->Write("", TObject::kOverwrite);




file1->cd("1diff_w_1225_215");
     
for (Short_t ii=0; ii<=36; ii++) {
h_odn_inv_m12[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<=36; ii++) {
h_odn_inv_m23[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<=36; ii++) {
h_odn_alpha[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<=36; ii++) {
h_odn_theta[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<=36; ii++) {
h_odn_theta_2[ii]->Write("", TObject::kOverwrite); 
};


file1->cd("1diff_w_21625_30375");
for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_inv_m12[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_inv_m23[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_alpha[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_theta[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<18; ii++) {
h_odn_wwide_theta_2[ii]->Write("", TObject::kOverwrite); 
};


file1->cd("1diff_w_30875_45375");
for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_inv_m12[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_inv_m23[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_alpha[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_theta[ii]->Write("", TObject::kOverwrite); 
};
for (Short_t ii=0; ii<15; ii++) {
h_odn_wgt3_theta_2[ii]->Write("", TObject::kOverwrite); 
};


file1->cd("int_w_for_q2_01_13");
for (Short_t ii=0; ii<12; ii++) {
h_odn_w_dep_tot[ii]->SetMinimum(0.);
h_odn_w_dep_t[ii]->SetMinimum(0.);
h_odn_w_dep_l[ii]->SetMinimum(0.);
h_odn_w_dep_l2[ii]->SetMinimum(0.);

h_odn_w_dep_t[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_l[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_l2[ii]->Write("", TObject::kOverwrite); 
h_odn_w_dep_tot[ii]->Write("", TObject::kOverwrite);

};

file1->cd("int_q2_for_w_125_2075");
for (Short_t ii=0; ii<33; ii++) {

h_odn_q2_dep_t[ii]->SetMinimum(0.);
h_odn_q2_dep_l[ii]->SetMinimum(0.);
h_odn_q2_dep_l2[ii]->SetMinimum(0.);
h_odn_q2_dep_tot[ii]->SetMinimum(0.);


h_odn_q2_dep_t[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_l[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_l2[ii]->Write("", TObject::kOverwrite); 
h_odn_q2_dep_tot[ii]->Write("", TObject::kOverwrite); 

};

/*
file1->cd("int_cr_sect_norm");
for (Short_t ii=0; ii<27; ii++) {
h_int_crsect_t[ii]->Scale(1./h_int_crsect_t[ii]->GetEntries()*71.);
h_int_crsect_l[ii]->Scale(1./h_int_crsect_l[ii]->GetEntries()*71.);
 };  
   
   
for (Short_t ii=0; ii<27; ii++) {
h_int_crsect_t[ii]->Write("", TObject::kOverwrite); 
h_int_crsect_l[ii]->Write("", TObject::kOverwrite); 
};
*/


file1->Write(); 


};
