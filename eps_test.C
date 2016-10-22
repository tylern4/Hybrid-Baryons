void eps_test() {
gStyle->SetOptStat(0);

Float_t W = 1.6;
Float_t Q2 = 0.45;
Float_t E_beam = 2.039;
Float_t MP = 0.939272;
Float_t phi_e = 0.;
Float_t M_PI = 3.141592653589793;

TLorentzVector P4_Eini, P4_Efin, P4_gamma;

Float_t En_gamma = (W*W+Q2-MP*MP)/2./MP;

Float_t En_Efin = E_beam - En_gamma;

Float_t Theta_Efin = acos(1.-Q2/E_beam/En_Efin/2.);


P4_Eini.SetXYZT(0.,0.,E_beam, E_beam);

P4_Efin.SetXYZT(En_Efin*cos(phi_e)*sin(Theta_Efin),En_Efin*sin(phi_e)*sin(Theta_Efin),En_Efin*cos(Theta_Efin),En_Efin);

P4_gamma = P4_Eini - P4_Efin;

Float_t eps_t = 1./(1.+ 2.*(1. + En_gamma*En_gamma/Q2)*tan(Theta_Efin/2.)*tan(Theta_Efin/2.));



Float_t eps_t2 = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_Eini.Vect()).Cross(P4_Efin.Vect())).Mag2()));



 TRotation rot;
 TVector3 uz = P4_gamma.Vect().Unit();
 TVector3 ux = (P4_Eini.Vect().Cross(P4_Efin.Vect())).Unit();
 ux.Rotate(3.*M_PI/2,uz);
 rot.SetZAxis(uz,ux).Invert();
P4_Eini.Transform(rot);
P4_Efin.Transform(rot);
 P4_gamma.Transform(rot);

//Float_t eps_t_after = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_Eini.Vect()).Cross(P4_Efin.Vect())).Mag2()));

P4_Eini.Boost(0.,0.,-sqrt(En_gamma*En_gamma+Q2)/(En_gamma+MP));
P4_Efin.Boost(0.,0.,-sqrt(En_gamma*En_gamma+Q2)/(En_gamma+MP));
P4_gamma.Boost(0.,0.,-sqrt(En_gamma*En_gamma+Q2)/(En_gamma+MP));

//cout << P4_Eini[3]<<" 2\n";
//cout << P4_gamma.Mag2()<< " \n";
Float_t En_Efin_after = P4_Efin[3];
Float_t Theta_Efin_after = P4_Efin.Theta();

//Float_t eps_t_after = 1./(1.+ 2.*(1. + En_Efin_after*En_Efin_after/Q2)*tan(Theta_Efin_after/2.)*tan(Theta_Efin_after/2.));

Float_t eps_t_after = 1./(1.+Q2*((P4_gamma.Vect()).Mag2())/2./(((P4_Eini.Vect()).Cross(P4_Efin.Vect())).Mag2()));
cout<< eps_t<<" "<< eps_t2<<" "<<eps_t_after<<" \n";
};
