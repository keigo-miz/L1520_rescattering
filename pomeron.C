/* based on Titov's paper (2003) */

#include "spinor_mat.h"

/* isoscalar EM FF of the nucleon */
double F1(double t) {
  double t0 = 0.7;  // [GeV^2]
  double a_N = 2.;  // Titov = TMath::Sqrt(2.8);
  double num = 4*sq(Mp) - sq(a_N)*t;
  double den = (4*sq(Mp) - t)*sq(1-t/t0);
  return (num/den);
}

/* FF of V-photon-Pomeron */
double Fv(double t) {
  double sq_mu0 = 1.1;  // [GeV^2]
  double num = 2*sq_mu0;
  double den = (1-t/sq(Mphi))*(2*sq_mu0+sq(Mphi)-t);
  return (num/den);
}

/* Pomeron trajectory */
double alpha_P(double t) {
  return (1.08 + 0.25*t);
}

/* scalar function */
TComplex M_P(double s, double t) {
  double sp = 4.;  // [GeV^2]
  double Cp = 0.618052;  // [GeV^{-2}]

  double real = Cp*F1(t)*Fv(t)*TMath::Power(s/sp,alpha_P(t)-1.);
  return (real * TComplex::Exp(-TComplex(0,(pi/2)*alpha_P(t))));
}

/* vertex function h */
Mat h_P(int i, int j, TLorentzVector k, TLorentzVector q, TLorentzVector pbar) {
  /* metric */
  double g[4][4] = {};
  g[0][0] = 1.; g[1][1] = -1.; g[2][2] = -1.; g[3][3] = -1.;
  /* deal with the difference of metric definition */
  double kk[4] = {k(3), k(0), k(1), k(2)};
  double qq[4] = {q(3), q(0), q(1), q(2)};
  double ppbar[4] = {pbar(3), pbar(0), pbar(1), pbar(2)};

  Mat kslash = G[0]*kk[0] - G[1]*kk[1] - G[2]*kk[2] - G[3]*kk[3];
  Mat qslash = G[0]*qq[0] - G[1]*qq[1] - G[2]*qq[2] - G[3]*qq[3];
  Mat term1 = kslash*(g[i][j]-qq[i]*qq[j]/q.M2());
  Mat term2 = G[j]*(kk[i]-qq[i]*(k*q)/q.M2());

  double qbar_j = qq[j] - ppbar[j]*(k*q)/(pbar*k);
  Mat term3 = (G[i] - qslash*(qq[i]/q.M2()))*qbar_j;
  return (term1 - term2 - term3);
}

TComplex Pomeron_amp(double Eg, double costh, int hel_gamma, int hel_phi, int mi, int mf) {
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
  TLorentzVector pbar = 0.5*(p+pp);
  double s = W.M2();
  double t = (k-q).M2();

  /* polarization vector */
  TComplex eps_gamma[4], eps_V[4];
  FillPolVector(k, hel_gamma, eps_gamma);
  FillPolVector(k, hel_phi, eps_V);  // GJ frame
  /* sup --> sub */
  for (int i=1; i<4; i++) {
    eps_gamma[i] *= TComplex(-1,0);
    eps_V[i] *= TComplex(-1,0);
  }

  TComplex Gamma = TComplex(0,0);
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      Gamma += TComplex::Conjugate(eps_V[i])*barDot(u(pp,mf),h_P(i,j,k,q,pbar)*u(p,mi))*eps_gamma[j];
    }
  }
  TComplex I = -M_P(s,t)*Gamma;
  return I;
}
