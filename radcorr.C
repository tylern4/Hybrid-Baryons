#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>



void radcorr() {
gStyle->SetOptStat(0);
Int_t i;
Float_t R;



//Opening file with histogram 

//TFile *file = new TFile("out_norad_2_039.root");
//TFile *file = new TFile("out_norad_2_039_wupto19.root");
//TFile *file = new TFile("out_norad_1_515.root");
//TFile *file = new TFile("out_norad_5_754_q2_20_24.root");
TFile *file = new TFile("out_norad_5_754_q2_w_upto_3.root");

TH2F* q2vsw_norad = (TH2F*)file->Get("Q2vsW"); 
TH2F* q2vsw_norad2 = (TH2F*)file->Get("Q2vsW2");
//for (Int_t i=1; i<=q2vsw_norad2->GetNbinsX(); i++){
//for (Int_t j=1; j<=q2vsw_norad2->GetNbinsY(); j++){

//q2vsw_norad->SetBinError(i,j,0.);
//cout << q2vsw_norad->GetBinContent(i,j) << "\n";
//};
//};

//q2vsw_norad->Scale(1./(5250000));


//TFile *file1 = new TFile("out_rad_2_039.root");
//TFile *file1 = new TFile("out_rad_2_039_wupto19.root");
//TFile *file1 = new TFile("out_rad_1_515.root");
//TFile *file1 = new TFile("out_rad_5_754_q2_20_24.root");
TFile *file1 = new TFile("out_rad_5_754_q2_w_upto_3.root");

//TH2F* q2vsw_rad = (TH2F*)file1->Get("h105");
 TH2F* q2vsw_rad = (TH2F*)file1->Get("Q2vsW");
  TH2F* q2vsw_rad2 = (TH2F*)file1->Get("Q2vsW2");
////for (Int_t i=1; i<=q2vsw_rad2->GetNbinsX(); i++){
//for (Int_t j=1; j<=q2vsw_rad2->GetNbinsY(); j++){

//q2vsw_rad->SetBinError(i,j,0.);

//};
//};
//


//q2vsw_rad->Scale(1./(5250000));
//q2vsw_rad->Scale(1./(6237.*250000.));

q2vsw_norad->Divide(q2vsw_rad);

for (Int_t i=1; i<=q2vsw_norad->GetNbinsX(); i++){
for (Int_t j=1; j<=q2vsw_norad->GetNbinsY(); j++){


//q2vsw_norad->SetBinError(i,j,q2vsw_norad->GetBinContent(i,j)*sqrt(1./q2vsw_rad2->GetBinContent(i,j)+1./q2vsw_norad2->GetBinContent(i,j)));
q2vsw_norad->SetBinError(i,j,0.);
//cout << q2vsw_norad->GetBinError(i,j) << "\n";
};
};

//q2vsw_norad->Draw("LEGO2");

TCanvas *c = new TCanvas("c","c",500,500);
c->cd();
//c->cd()->SetBottomMargin(0.2);
c->cd()->SetLeftMargin(0.15);
TH1 *odn = q2vsw_norad->ProjectionX("odn",1,4);

odn->SetMarkerStyle(20);
//odn->Scale(0.1);
odn->SetMinimum(0.95);
odn->SetMaximum(1.4);
odn->SetLineColor(kBlack);
odn->SetMarkerColor(kBlack);
odn->SetTitle("E_{beam} = 5.754 GeV");

odn->GetYaxis()->SetLabelSize(0.04);
odn->GetXaxis()->SetLabelSize(0.04);
odn->GetXaxis()->SetTitle("W, GeV");
odn->GetYaxis()->SetTitle("factor");
odn->GetYaxis()->SetTitleOffset(1.5);
odn->GetXaxis()->SetNdivisions(10);
odn->GetYaxis()->SetNdivisions(15);
odn->GetXaxis()->SetTitleSize(0.04);
odn->GetYaxis()->SetTitleSize(0.04);

odn->Scale(1./4.);
odn->Draw("E1P");



//TH1 *odnq = q2vsw_rad->ProjectionX("odnq",2,2);
//odnq->SetMarkerStyle(20);
//odnq->SetMarkerColor(kRed);
//odnq->Draw("E1PX0 same");




TH1 *odn2 = q2vsw_norad->ProjectionX("odn2",2,2);
odn2->SetMarkerStyle(20);
odn2->SetLineColor(kBlue);
odn2->SetMarkerColor(kBlue);
//odn2->Draw("E1PX0 same");


TH1 *odn3 = q2vsw_norad->ProjectionX("odn3",3,3);
odn3->SetMarkerStyle(20);
odn3->SetLineColor(kGreen);
odn3->SetMarkerColor(kGreen);
//odn3->Draw("E1PX0 same");

TH1 *odn4 = q2vsw_norad->ProjectionX("odn4",4,4);
odn4->SetMarkerStyle(20);
odn4->SetLineColor(kRed);
odn4->SetMarkerColor(kRed);
//odn4->Draw("E1PX0 same");


TH1 *odn5 = q2vsw_norad->ProjectionX("odn5",5,5);
odn5->SetMarkerStyle(20);
odn5->SetLineColor(6);
odn5->SetMarkerColor(6);
//odn5->Draw("E1PX0 same");


TH1 *odn6 = q2vsw_norad->ProjectionX("odn6",10,10);
odn6->SetMarkerStyle(20);
odn6->SetLineColor(9);
odn6->SetMarkerColor(9);
//odn6->Draw("E1PX0 same");


leg = new TLegend(0.1,0.3,0.4,0.4);
leg->AddEntry(odn,"0.4 GeV^{2} < Q^{2} < 0.45 GeV^{2} ","p");
leg->AddEntry(odn2,"0.45 GeV^{2} < Q^{2} < 0.5 GeV^{2} ","p");
leg->AddEntry(odn3,"0.5 GeV^{2} < Q^{2} < 0.55 GeV^{2} ","p");
leg->AddEntry(odn4,"0.55 GeV^{2} < Q^{2} <  0.6 GeV^{2} ","p");
//leg->AddEntry(odn,"2. GeV^{2} < Q^{2} < 2.1 GeV^{2} ","p");
//leg->AddEntry(odn2,"2.1 GeV^{2} < Q^{2} < 2.2 GeV^{2} ","p");
//leg->AddEntry(odn3,"2.3 GeV^{2} < Q^{2} < 2.3 GeV^{2} ","p");
//leg->AddEntry(odn4,"2.3 GeV^{2} < Q^{2} <  2.4 GeV^{2} ","p");
//leg->AddEntry(odn5,"0.75 GeV^{2} < Q^{2} < 0.8 GeV^{2} ","p");
//leg->Draw();


};
