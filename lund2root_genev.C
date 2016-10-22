#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void lund2root_genev(string input_filename){
  
gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);

const double DEG=180./3.1415926;

char the_filename[200];
sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());

char output_filename[200];
sprintf(output_filename, "%s.root",the_filename);
TFile *outputfile=new TFile(output_filename, "recreate");

TTree *T = new TTree("h2","h2");

  Float_t xsec;
  Float_t e_e,p_e,th_e,phi_e;
  Float_t e_prot,p_prot,th_prot,phi_prot;
  Float_t e_pp1,p_pp1,th_pp1,phi_pp1;
  Float_t e_pm1,p_pm1,th_pm1,phi_pm1;
TBranch *B1 = T->Branch("xsec",&xsec,"data/F");
TBranch *B2 = T->Branch("e_e",&e_e,"data/F");
TBranch *B3 = T->Branch("p_e",&p_e,"data/F");
TBranch *B4 = T->Branch("th_e",&th_e,"data/F");
TBranch *B5 = T->Branch("phi_e",&phi_e,"data/F");
TBranch *B6 = T->Branch("e_prot",&e_prot,"e_prot/F");
TBranch *B7 = T->Branch("p_prot",&p_prot,"data/F");
TBranch *B8 = T->Branch("th_prot",&th_prot,"data/F");
TBranch *B9 = T->Branch("phi_prot",&phi_prot,"data/F");
TBranch *B10 = T->Branch("e_pp1",&e_pp1,"data/F");
TBranch *B11 = T->Branch("p_pp1",&p_pp1,"data/F");
TBranch *B12 = T->Branch("th_pp1",&th_pp1,"data/F");
TBranch *B13 = T->Branch("phi_pp1",&phi_pp1,"data/F");
TBranch *B14 = T->Branch("e_pm1",&e_pm1,"data/F");
TBranch *B15 = T->Branch("p_pm1",&p_pm1,"data/F");
TBranch *B16 = T->Branch("th_pm1",&th_pm1,"data/F");
TBranch *B17 = T->Branch("phi_pm1",&phi_pm1,"data/F");

  ifstream input(input_filename.c_str());
    
  if (!input.good()) {cout << "file doesn't exist" << endl;  return;}
    
  char textline[200];  
  
  double tmp;
  double px_e,py_e,pz_e;
  double px_prot,py_prot,pz_prot;  
  double px_pp1,py_pp1,pz_pp1;  
  double px_pm1,py_pm1,pz_pm1;  
  
  int i=0;
   
//   while(i<2){  
  while(!input.eof()){
//     input >> tmp;
//     if(input.eof()) break; 
//   input  >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  tmp  >>  tmp  >>  tmp  >>  xsec;    
  input >> tmp >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  tmp  >>  tmp  >>  tmp  >>  xsec;
//   input.getline(textline,200);  
  input >> tmp >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  px_e  >>  py_e  >>  pz_e  >>  e_e >> tmp >> tmp  >>  tmp   >> tmp;  
//   input.getline(textline,200);  
  input >> tmp >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  px_prot  >>  py_prot  >>  pz_prot  >>  e_prot >> tmp >> tmp  >>  tmp   >> tmp;  
//   input.getline(textline,200);  
  input >> tmp >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  px_pp1  >>  py_pp1  >>  pz_pp1  >>  e_pp1 >> tmp >> tmp  >>  tmp   >> tmp;  
//   input.getline(textline,200);  
  input >> tmp >> tmp  >>  tmp  >>  tmp  >>  tmp   >> tmp >>  px_pm1  >>  py_pm1  >>  pz_pm1  >>  e_pm1 >> tmp >> tmp  >>  tmp   >> tmp;  
//   input.getline(textline,200);  

  if(input.eof()) break; 
 
//   cout << xsec << endl;
//   cout << px_e << " " << py_e << " " << pz_e << endl;
//   cout << px_prot << " " << py_prot << " " << pz_prot << endl;
//   cout << px_pp1 << " " << py_pp1 << " " << pz_pp1 << endl;
//   cout << px_pm1 << " " << py_pm1 << " " << pz_pm1 << endl;  
  
  p_e=sqrt(px_e*px_e+py_e*py_e+pz_e*pz_e);
  th_e=acos(pz_e/p_e)*DEG;
  phi_e=atan2(py_e,px_e)*DEG;
  if (phi_e<0) phi_e=360+phi_e;
  p_prot=sqrt(px_prot*px_prot+py_prot*py_prot+pz_prot*pz_prot);
  th_prot=acos(pz_prot/p_prot)*DEG;
  phi_prot=atan2(py_prot,px_prot)*DEG;
  if (phi_prot<0) phi_prot=360+phi_prot;
  p_pp1=sqrt(px_pp1*px_pp1+py_pp1*py_pp1+pz_pp1*pz_pp1);
  th_pp1=acos(pz_pp1/p_pp1)*DEG;
  phi_pp1=atan2(py_pp1,px_pp1)*DEG;
  if (phi_pp1<0) phi_pp1=360+phi_pp1;
  p_pm1=sqrt(px_pm1*px_pm1+py_pm1*py_pm1+pz_pm1*pz_pm1);
  th_pm1=acos(pz_pm1/p_pm1)*DEG;
  phi_pm1=atan2(py_pm1,px_pm1)*DEG;
  if (phi_pm1<0) phi_pm1=360+phi_pm1;

  //   B1->Fill();   B2->Fill();  B3->Fill();  B4->Fill();  B5->Fill();  B6->Fill();  B7->Fill();  B8->Fill();  B9->Fill();  B10->Fill();  B11->Fill();  B12->Fill();  B13->Fill();  B14->Fill();  B15->Fill();  B16->Fill();  B17->Fill();
  T->Fill();
  
  i++;
   cout << i << "\r";
//    cout << i << "\n";
   
//   if (i==1000000) break;   
  }

  cout << "\n total events " << i << endl;
  
T->Write();
outputfile->Write();
outputfile->Close();  

}
