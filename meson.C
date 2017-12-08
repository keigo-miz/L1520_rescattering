/* based on Titov's paper (2003) */

#include "spinor_mat.h"

/* iPart 1: pion,  2: eta */

/* form factor */
double FF(int FFid, double t, int iPart) {
  /* FFid 1: F_{\varphi NN},  2: F_{V\gamma\varphi} */
  /* Ryu (2014)'s parameters */
  double M[2] = {Mpi0, Meta};
  double L;
  if (FFid == 1) {
         if (iPart == 1) L = 0.7;  // [GeV]
    else if (iPart == 2) L = 1.0;  // [GeV]
  } else if (FFid == 2) {
         if (iPart == 1) L = 0.77;  // [GeV]
    else if (iPart == 2) L = 0.9;  // [GeV]
  }
  return ((sq(L)-sq(M[iPart-1])) / (sq(L)-t));
}

TComplex Meson_amp(double Eg, double costh, int hel_gamma, int hel_phi, int mi, int mf, int iPart) {
  TLorentzVector k(0, 0, Eg, Eg), p(0,0,0,Mp);  // lab system
  TLorentzVector W = k + p;
  /* lab --> CM */
  k.Boost(-W.BoostVector());
  p.Boost(-W.BoostVector());

  double pdk = TMath::Sqrt((W.M2()-sq(Mp+Mphi))*(W.M2()-sq(Mp-Mphi)))/(2*W.M());
  TLorentzVector q;
  q.SetXYZM(pdk*TMath::Sqrt(1-sq(costh)),0,pdk*costh,Mphi);
  /* CM --> V rest frame */
  k.Boost(-q.BoostVector());
  p.Boost(-q.BoostVector());
  q.Boost(-q.BoostVector());

  TLorentzVector pp = k + p - q;
  double s = W.M2();
  double t = (k-q).M2();

  double M, gg, gN;
  if (iPart == 1) {
    M = Mpi0;
    gg = -0.141;
    gN = 13.26;
  } else if (iPart == 2) {
    M = Meta;
    gg = -0.707;
    gN = 3.527;
  }
  double e = TMath::Sqrt(4*pi/137.04);
  double A = FF(1, t, iPart)*FF(2, t, iPart)*e*gg*gN/(t-sq(M))/Mphi;

  /* polarization vector */
  TComplex eps_gamma[4], eps_V[4], eps_Vc[4];
  FillPolVector(k, hel_gamma, eps_gamma);
  FillPolVector(k, hel_phi, eps_V);  // GJ frame
  for (int i=0; i<4; i++) eps_Vc[i] = TComplex::Conjugate(eps_V[i]);

  /* deal with the difference of metric definition */
  TComplex kk[4] = {TComplex(k(3),0), TComplex(k(0),0), TComplex(k(1),0), TComplex(k(2),0)};
  TComplex qq[4] = {TComplex(q(3),0), TComplex(q(0),0), TComplex(q(1),0), TComplex(q(2),0)};

  TComplex sum = TComplex(0,0);
  sum += qq[0]*eps_Vc[1]*kk[2]*eps_gamma[3];
  sum -= qq[0]*eps_Vc[1]*kk[3]*eps_gamma[2];
  sum -= qq[0]*eps_Vc[2]*kk[1]*eps_gamma[3];
  sum += qq[0]*eps_Vc[2]*kk[3]*eps_gamma[1];
  sum += qq[0]*eps_Vc[3]*kk[1]*eps_gamma[2];
  sum -= qq[0]*eps_Vc[3]*kk[2]*eps_gamma[1];
  sum -= qq[1]*eps_Vc[0]*kk[2]*eps_gamma[3];
  sum += qq[1]*eps_Vc[0]*kk[3]*eps_gamma[2];
  sum += qq[1]*eps_Vc[2]*kk[0]*eps_gamma[3];
  sum -= qq[1]*eps_Vc[2]*kk[3]*eps_gamma[0];
  sum -= qq[1]*eps_Vc[3]*kk[0]*eps_gamma[2];
  sum += qq[1]*eps_Vc[3]*kk[2]*eps_gamma[0];
  sum += qq[2]*eps_Vc[0]*kk[1]*eps_gamma[3];
  sum -= qq[2]*eps_Vc[0]*kk[3]*eps_gamma[1];
  sum -= qq[2]*eps_Vc[1]*kk[0]*eps_gamma[3];
  sum += qq[2]*eps_Vc[1]*kk[3]*eps_gamma[0];
  sum += qq[2]*eps_Vc[3]*kk[0]*eps_gamma[1];
  sum -= qq[2]*eps_Vc[3]*kk[1]*eps_gamma[0];
  sum -= qq[3]*eps_Vc[0]*kk[1]*eps_gamma[2];
  sum += qq[3]*eps_Vc[0]*kk[2]*eps_gamma[1];
  sum += qq[3]*eps_Vc[1]*kk[0]*eps_gamma[2];
  sum -= qq[3]*eps_Vc[1]*kk[2]*eps_gamma[0];
  sum -= qq[3]*eps_Vc[2]*kk[0]*eps_gamma[1];
  sum += qq[3]*eps_Vc[2]*kk[1]*eps_gamma[0];

  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);
  TComplex B = barDot(u(pp,mf),Gamma5*u(p,mi));

  return (TComplex(0,A)*B*sum);
}
