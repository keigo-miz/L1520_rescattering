/* based on Ryu's paper (2014) */

#include "../spinor_mat.h"
#define MLs (1.5195)  // L(1520) mass
#define Mk (0.493677)  // K+ mass

/* spin-3/2 spinor */
TComplex ee(int mu, int lambda, TLorentzVector p) {  // Ryu's thesis (A.13)
  TComplex ehat[3][3] = {  // ehat[0] = ehat_-, ehat[1] = ehat_0, ehat[2] = ehat_+ (A.14)
    {TComplex(1/TMath::Sqrt(2),0), TComplex(0,-1/TMath::Sqrt(2)), TComplex(0,0)}, // -
    {TComplex(0,0), TComplex(0,0), TComplex(1,0)},                                // 0
    {TComplex(-1/TMath::Sqrt(2),0), TComplex(0,-1/TMath::Sqrt(2)), TComplex(0,0)} // +
  };
  TComplex dot = TComplex(0,0);
  for (int i=0; i<3; i++) dot += ehat[lambda+1][i]*TComplex(p(i),0);
  dot *= TComplex(1/MLs,0);

  TComplex ret;
  if (mu==0) {
    ret = dot;
  } else {
    ret = ehat[lambda+1][mu-1] + TComplex(p(mu-1)/(p.E()+MLs),0)*dot;
  }
  return ret;
}

Spi u32(int mu, TLorentzVector p, int S) {
  Spi ret;
  if (S==3) {
    ret = u(1,p)*ee(mu,1,p);
  } else if (S==1) {
    ret = u(1,p)*(ee(mu,0,p)*TComplex(TMath::Sqrt(2./3.),0)) + u(-1,p)*(ee(mu,1,p)*TComplex(TMath::Sqrt(1./3.),0));
  } else if (S==-1) {
    ret = u(1,p)*(ee(mu,-1,p)*TComplex(TMath::Sqrt(1./3.),0)) + u(-1,p)*(ee(mu,0,p)*TComplex(TMath::Sqrt(2./3.),0));
  } else if (S==-3) {
    ret = u(-1,p)*ee(mu,-1,p);
  } else {
    std::cerr << "[error] invalid S in u32" << std::endl;
    exit(0);
  }
  return ret;
}

TComplex M_L_s(double Eg, int hel_gamma, int mi, TLorentzVector k2, int sign1, int sign2) {
  TLorentzVector k1(0, 0, Eg, Eg), p1(0,0,0,Mp);  // lab system
  TVector3 b = (k1 + p1).BoostVector();
  /* lab --> CM */
  k1.Boost(-b);
  p1.Boost(-b);

  TLorentzVector qs = k1 + p1;
  TLorentzVector p2 = k1 + p1 - k2;

  /* constants */
  double e = TMath::Sqrt(4*pi/137.04);
  double g = 11.;
  double kappa = 1.79;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(hel_gamma, k1, eps);

  /* deal with the difference of metric definitions */
  double kk1[4] = {k1(3), k1(0), k1(1), k1(2)};
  double kk2[4] = {k2(3), k2(0), k2(1), k2(2)};
  double qqs[4] = {qs(3), qs(0), qs(1), qs(2)};

  /* slashed matrices */
  Mat k1slash = G[0]*kk1[0] - G[1]*kk1[1] - G[2]*kk1[2] - G[3]*kk1[3];
  Mat qsslash = G[0]*qqs[0] - G[1]*qqs[1] - G[2]*qqs[2] - G[3]*qqs[3];
  Mat eslash  = G[0]*eps[0] - G[1]*eps[1] - G[2]*eps[2] - G[3]*eps[3];

  Spi tmp0 = (Gamma5*qsslash+Gamma5*Mp)*(1/(qs.M2()-sq(Mp)))*eslash*u(mi,p1);
  int S = 2*hel_gamma + mi;
  TComplex tmp1 = barDot(u32(0,p2,S),tmp0*kk2[0])
                - barDot(u32(1,p2,S),tmp0*kk2[1])
                - barDot(u32(2,p2,S),tmp0*kk2[2])
                - barDot(u32(3,p2,S),tmp0*kk2[3]);
  TComplex term1 = TComplex(0,sign1*e*g/Mk)*tmp1;

  Spi tmp2 = (Gamma5*qsslash+Gamma5*Mp)*(1/(qs.M2()-sq(Mp)))*eslash*k1slash*u(mi,p1);
  TComplex tmp3 = barDot(u32(0,p2,S),tmp2*kk2[0])
                - barDot(u32(1,p2,S),tmp2*kk2[1])
                - barDot(u32(2,p2,S),tmp2*kk2[2])
                - barDot(u32(3,p2,S),tmp2*kk2[3]);
  TComplex term2 = TComplex(0,sign2*e*kappa*g/(2*Mk*Mp))*tmp3;

  return (term1 + term2);
}

TComplex M_L_t(double Eg, int hel_gamma, int mi, TLorentzVector k2, int sign) {
  TLorentzVector k1(0, 0, Eg, Eg), p1(0,0,0,Mp);  // lab system
  TVector3 b = (k1 + p1).BoostVector();
  /* lab --> CM */
  k1.Boost(-b);
  p1.Boost(-b);

  TLorentzVector p2 = k1 + p1 - k2;
  TLorentzVector qt = k1 - k2;

  /* constants */
  double e = TMath::Sqrt(4*pi/137.04);
  double g = 11.;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(hel_gamma, k1, eps);

  /* deal with the difference of metric definition */
  TComplex kk2[4] = {TComplex(k2(3),0), TComplex(k2(0),0), TComplex(k2(1),0), TComplex(k2(2),0)};
  double qqt[4] = {qt(3), qt(0), qt(1), qt(2)};

  TComplex tmp0 = kk2[0]*eps[0] - kk2[1]*eps[1] - kk2[2]*eps[2] - kk2[3]*eps[3];
//  tmp0 = TComplex(1,0);

  int S = 2*hel_gamma + mi;
  TComplex tmp1 = barDot(u32(0,p2,S),Gamma5*qqt[0]*u(mi,p1))
                - barDot(u32(1,p2,S),Gamma5*qqt[1]*u(mi,p1))
                - barDot(u32(2,p2,S),Gamma5*qqt[2]*u(mi,p1))
                - barDot(u32(3,p2,S),Gamma5*qqt[3]*u(mi,p1));

  return (TComplex(0,sign*2*e*g/(Mk*(qt.M2()-sq(Mk))))*tmp1*tmp0);
}

TComplex M_L_c(double Eg, int hel_gamma, int mi, TLorentzVector k2, int sign) {
  TLorentzVector k1(0, 0, Eg, Eg), p1(0,0,0,Mp);  // lab system
  TVector3 b = (k1 + p1).BoostVector();
  /* lab --> CM */
  k1.Boost(-b);
  p1.Boost(-b);

  TLorentzVector p2 = k1 + p1 - k2;

  /* constants */
  double e = TMath::Sqrt(4*pi/137.04);
  double g = 11.;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(hel_gamma, k1, eps);

  int S = 2*hel_gamma + mi;
  TComplex tmp0 = barDot(u32(0,p2,S),Gamma5*eps[0]*u(mi,p1))
                - barDot(u32(1,p2,S),Gamma5*eps[1]*u(mi,p1))
                - barDot(u32(2,p2,S),Gamma5*eps[2]*u(mi,p1))
                - barDot(u32(3,p2,S),Gamma5*eps[3]*u(mi,p1));

  return (TComplex(0,sign*e*g/Mk)*tmp0);
}

TComplex FRL(double s, double t, double L1, double L2) {
  int n1, n2;
  n1 = 1;
  n2 = 1;
//    L1 = 1.5;  // [GeV]
//    L2 = 1.5;  // [GeV]
  double nL4_1 = n1*L1*L1*L1*L1;
  double nL4_2 = n2*L2*L2*L2*L2;
  double left = TMath::Power(nL4_1 / (nL4_1 + sq(s-sq(Mp))), (Double_t) n1);
  double right = TMath::Power(nL4_2 / (nL4_2 + sq(t)), (Double_t) n2);
  return TComplex(left*right,0);
}

TComplex M_L(double Eg, int hel_gamma, int mi, TLorentzVector k2, int sign1, int sign2, int sign3, int sign4, double L1, double L2) {
  double s = 2*Mp*Eg + sq(Mp);
  double Eg_cm = (s - sq(Mp)) / (2*TMath::Sqrt(s));
  TLorentzVector k1(0,0,Eg_cm,Eg_cm);
  double t = (k1 - k2).M2();
  return ((M_L_s(Eg,hel_gamma,mi,k2,sign1,sign2) + M_L_t(Eg,hel_gamma,mi,k2,sign3) + M_L_c(Eg,hel_gamma,mi,k2,sign4))*FRL(s,t,L1,L2));
}
