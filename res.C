/* based on Ryu's paper (2014) */

#include "../amp_calc/spinor_mat.h"
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
    ret = u(p,1)*ee(mu,1,p);
  } else if (S==1) {
    ret = u(p,1)*(ee(mu,0,p)*TComplex(TMath::Sqrt(2./3.),0)) + u(p,-1)*(ee(mu,1,p)*TComplex(TMath::Sqrt(1./3.),0));
  } else if (S==-1) {
    ret = u(p,1)*(ee(mu,-1,p)*TComplex(TMath::Sqrt(1./3.),0)) + u(p,-1)*(ee(mu,0,p)*TComplex(TMath::Sqrt(2./3.),0));
  } else if (S==-3) {
    ret = u(p,-1)*ee(mu,-1,p);
  } else {
    std::cerr << "[error] invalid S in u32" << std::endl;
    exit(0);
  }
  return ret;
}

TComplex I_Ls(double Eg, int hel_gamma, int mi, TLorentzVector k2) {
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
  FillPolVector(k1, hel_gamma, eps);

  /* deal with the difference of metric definitions */
  double kk1[4] = {k1(3), k1(0), k1(1), k1(2)};
  double kk2[4] = {k2(3), k2(0), k2(1), k2(2)};
  double qqs[4] = {qs(3), qs(0), qs(1), qs(2)};

  /* slashed matrices */
  Mat k1slash = G[0]*kk1[0] - G[1]*kk1[1] - G[2]*kk1[2] - G[3]*kk1[3];
  Mat qsslash = G[0]*qqs[0] - G[1]*qqs[1] - G[2]*qqs[2] - G[3]*qqs[3];
  Mat eslash  = G[0]*eps[0] - G[1]*eps[1] - G[2]*eps[2] - G[3]*eps[3];

  Spi tmp0 = (Gamma5*qsslash+Gamma5*Mp)*(1/(qs.M2()-sq(Mp)))*eslash*u(p1,mi);
  int S = 2*hel_gamma + mi;
  TComplex tmp1 = barDot(u32(0,p2,S),tmp0*kk2[0])
                - barDot(u32(1,p2,S),tmp0*kk2[1])
                - barDot(u32(2,p2,S),tmp0*kk2[2])
                - barDot(u32(3,p2,S),tmp0*kk2[3]);
  TComplex term1 = TComplex(0,-e*g/Mk)*tmp1;

  Spi tmp2 = (Gamma5*qsslash+Gamma5*Mp)*(1/(qs.M2()-sq(Mp)))*eslash*k1slash*u(p1,mi);
  TComplex tmp3 = barDot(u32(0,p2,S),tmp2*kk2[0])
                - barDot(u32(1,p2,S),tmp2*kk2[1])
                - barDot(u32(2,p2,S),tmp2*kk2[2])
                - barDot(u32(3,p2,S),tmp2*kk2[3]);
  TComplex term2 = TComplex(0,e*kappa*g/(2*Mk*Mp))*tmp3;

  return (term1 + term2);
}

TComplex I_Lt(double Eg, int hel_gamma, int mi, TLorentzVector k2) {
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
  FillPolVector(k1, hel_gamma, eps);

  /* deal with the difference of metric definition */
  TComplex kk2[4] = {TComplex(k2(3),0), TComplex(k2(0),0), TComplex(k2(1),0), TComplex(k2(2),0)};
  double qqt[4] = {qt(3), qt(0), qt(1), qt(2)};

  TComplex tmp0 = kk2[0]*eps[0] - kk2[1]*eps[1] - kk2[2]*eps[2] - kk2[3]*eps[3];

  int S = 2*hel_gamma + mi;
  TComplex tmp1 = barDot(u32(0,p2,S),Gamma5*qqt[0]*u(p1,mi))
                - barDot(u32(1,p2,S),Gamma5*qqt[1]*u(p1,mi))
                - barDot(u32(2,p2,S),Gamma5*qqt[2]*u(p1,mi))
                - barDot(u32(3,p2,S),Gamma5*qqt[3]*u(p1,mi));

  return (TComplex(0,2*e*g/(Mk*(qt.M2()-sq(Mk))))*tmp1*tmp0);
}

TComplex I_Lc(double Eg, int hel_gamma, int mi, TLorentzVector k2) {
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
  FillPolVector(k1, hel_gamma, eps);

  int S = 2*hel_gamma + mi;
  TComplex tmp0 = barDot(u32(0,p2,S),Gamma5*eps[0]*u(p1,mi))
                - barDot(u32(1,p2,S),Gamma5*eps[1]*u(p1,mi))
                - barDot(u32(2,p2,S),Gamma5*eps[2]*u(p1,mi))
                - barDot(u32(3,p2,S),Gamma5*eps[3]*u(p1,mi));

  return (TComplex(0,e*g/Mk)*tmp0);
}

TComplex I_Rs(double Eg, double costh, int hel_phi, int mf, TLorentzVector k1) {
  double W = TMath::Sqrt(2*Mp*Eg+sq(Mp));
  double Eg_cm = (sq(W) - sq(Mp)) / (2*W);
  TLorentzVector pPhoton(0,0,Eg_cm,Eg_cm);
  TLorentzVector qs(0,0,0,W);
  /* 4-mom of outgoing phi-meson */
  double E_phi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));
  TLorentzVector k2(p_phi*TMath::Sqrt(1-sq(costh)),0,p_phi*costh,E_phi);
  /* CM --> phi-rest frame */
  TVector3 b = k2.BoostVector();
  pPhoton.Boost(-b);
  qs.Boost(-b);
  k1.Boost(-b);
  k2.Boost(-b);
  TLorentzVector p1 = qs - k1;
  TLorentzVector p2 = qs - k2;

  /* constants */
  double g2 = 11. * 0.25;
  double kappa = 0.2;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(pPhoton, hel_phi, eps);

  /* deal with the difference of metric definition */
  double kk1[4] = {k1(3), k1(0), k1(1), k1(2)};
  double kk2[4] = {k2(3), k2(0), k2(1), k2(2)};
  double qqs[4] = {qs(3), qs(0), qs(1), qs(2)};

  /* slashed matrices */
  Mat qsslash = G[0]*qqs[0] - G[1]*qqs[1] - G[2]*qqs[2] - G[3]*qqs[3];
  Mat k2slash = G[0]*kk2[0] - G[1]*kk2[1] - G[2]*kk2[2] - G[3]*kk2[3];
  Mat eslash = G[0]*TComplex::Conjugate(eps[0])
             - G[1]*TComplex::Conjugate(eps[1])
             - G[2]*TComplex::Conjugate(eps[2])
             - G[3]*TComplex::Conjugate(eps[3]);

  int S = 2*hel_phi + mf;
  Spi k1u = u32(0,p1,S)*kk1[0]
          - u32(1,p1,S)*kk1[1]
          - u32(2,p1,S)*kk1[2]
          - u32(3,p1,S)*kk1[3];

  Spi tmp0 = (eslash*qsslash + eslash*Mp)*(1/(qs.M2()-sq(Mp)))*Gamma5*k1u;
  TComplex term1 = TComplex(0,g2/Mk)*barDot(u(p2,mf),tmp0);
  TComplex term2 = TComplex(0,-g2*kappa/(2*Mk*Mp))*barDot(u(p2,mf),k2slash*tmp0);

  return (term1 + term2);
}

TComplex I_Rt(double Eg, double costh, int hel_phi, int mf, TLorentzVector k1) {
  double W = TMath::Sqrt(2*Mp*Eg+sq(Mp));
  double Eg_cm = (sq(W) - sq(Mp)) / (2*W);
  TLorentzVector pPhoton(0,0,Eg_cm,Eg_cm);
  TLorentzVector qs(0,0,0,W);
  /* 4-mom of outgoing phi-meson */
  double E_phi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));
  TLorentzVector k2(p_phi*TMath::Sqrt(1-sq(costh)),0,p_phi*costh,E_phi);
  /* CM --> phi-rest frame */
  TVector3 b = k2.BoostVector();
  pPhoton.Boost(-b);
  qs.Boost(-b);
  k1.Boost(-b);
  k2.Boost(-b);
  TLorentzVector p1 = qs - k1;
  TLorentzVector p2 = qs - k2;
  TLorentzVector qt = k1 - k2;

  /* constants */
  double g2 = 11. * 4.7;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(pPhoton, hel_phi, eps);

  /* deal with the difference of metric definition */
  double qqt[4] = {qt(3), qt(0), qt(1), qt(2)};
  TComplex kk1[4] = {TComplex(k1(3),0), TComplex(k1(0),0), TComplex(k1(1),0), TComplex(k1(2),0)};

  int S = 2*hel_phi + mf;
  Spi tmp0 = Gamma5*qqt[0]*u32(0,p1,S)
           - Gamma5*qqt[1]*u32(1,p1,S)
           - Gamma5*qqt[2]*u32(2,p1,S)
           - Gamma5*qqt[3]*u32(3,p1,S);

  TComplex tmp1 = kk1[0]*TComplex::Conjugate(eps[0])
                - kk1[1]*TComplex::Conjugate(eps[1])
                - kk1[2]*TComplex::Conjugate(eps[2])
                - kk1[3]*TComplex::Conjugate(eps[3]);

  return (TComplex(0,2*g2/(Mk*(qt.M2()-sq(Mk))))*tmp1*barDot(u(p2,mf),tmp0));
}

TComplex I_Rc(double Eg, double costh, int hel_phi, int mf, TLorentzVector k1) {
  double W = TMath::Sqrt(2*Mp*Eg+sq(Mp));
  double Eg_cm = (sq(W) - sq(Mp)) / (2*W);
  TLorentzVector pPhoton(0,0,Eg_cm,Eg_cm);
  TLorentzVector qs(0,0,0,W);
  /* 4-mom of outgoing phi-meson */
  double E_phi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));
  TLorentzVector k2(p_phi*TMath::Sqrt(1-sq(costh)),0,p_phi*costh,E_phi);
  /* CM --> phi-rest frame */
  TVector3 b = k2.BoostVector();
  pPhoton.Boost(-b);
  qs.Boost(-b);
  k1.Boost(-b);
  k2.Boost(-b);
  TLorentzVector p1 = qs - k1;
  TLorentzVector p2 = qs - k2;

  /* constants */
  double g2 = 11. * 4.7;
  Mat Gamma5 = G[0]*G[1]*G[2]*G[3]*TComplex(0,1);

  /* polarization vector */
  TComplex eps[4];
  FillPolVector(pPhoton, hel_phi, eps);

  int S = 2*hel_phi + mf;
  Spi tmp0 = u32(0,p1,S)*TComplex::Conjugate(eps[0])
           - u32(1,p1,S)*TComplex::Conjugate(eps[1])
           - u32(2,p1,S)*TComplex::Conjugate(eps[2])
           - u32(3,p1,S)*TComplex::Conjugate(eps[3]);

  return (TComplex(0,g2/Mk)*barDot(u(p2,mf),Gamma5*tmp0));
}

TComplex FRL(int IsR, double s, double t) {
  int n1, n2;
  double L1, L2;
  if (IsR == 1) {
    n1 = 1;
    n2 = 1;
    L1 = 1.59;  // [GeV]
    L2 = 1.5;  // [GeV]
  } else {
    n1 = 1;
    n2 = 1;
    L1 = 1.478;  // [GeV]
    L2 = 1.73;  // [GeV]
  }
  double nL4_1 = n1*L1*L1*L1*L1;
  double nL4_2 = n2*L2*L2*L2*L2;
  double left = TMath::Power(nL4_1 / (nL4_1 + sq(s-sq(Mp))), (Double_t) n1);
  double right = TMath::Power(nL4_2 / (nL4_2 + sq(t)), (Double_t) n2);
  return TComplex(left*right,0);
}

TComplex I_L(double Eg, int hel_gamma, int mi, TLorentzVector k2) {
  double s = 2*Mp*Eg + sq(Mp);
  double Eg_cm = (s - sq(Mp)) / (2*TMath::Sqrt(s));
  TLorentzVector k1(0,0,Eg_cm,Eg_cm);
  double t = (k1 - k2).M2();
  return ((I_Ls(Eg,hel_gamma,mi,k2) + I_Lt(Eg,hel_gamma,mi,k2) + I_Lc(Eg,hel_gamma,mi,k2))*FRL(0,s,t));
}

TComplex I_R(double Eg, double costh, int hel_phi, int mf, TLorentzVector k1) {
  double s = 2*Mp*Eg + sq(Mp);
  double E_phi = (s + sq(Mphi) - sq(Mp)) / (2*TMath::Sqrt(s));
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));
  TLorentzVector k2(p_phi*TMath::Sqrt(1-sq(costh)),0,p_phi*costh,E_phi);
  double t = (k1 - k2).M2();
  return ((I_Rs(Eg,costh,hel_phi,mf,k1) + I_Rt(Eg,costh,hel_phi,mf,k1) + I_Rc(Eg,costh,hel_phi,mf,k1))*FRL(1,s,t));
}

TComplex ImM(double Eg, double costh, int hel_gamma, int hel_phi, int mi, int mf) {
  if (2*hel_gamma+mi != 2*hel_phi+mf) return TComplex(0,0);
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);

  double dx = 0.05;
  TComplex sum = TComplex(0,0);
  TLorentzVector pK;
  /* integration */
//  for (double cosK=-1.; cosK<1.; cosK+=dx) {
//    for (double phi=0.; phi<2*pi; phi+=dx) {
//      pK = TLorentzVector(r*TMath::Sqrt(1-sq(cosK))*TMath::Cos(phi), r*TMath::Sqrt(1-sq(cosK))*TMath::Sin(phi), r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));
//      sum += I_L(Eg, hel_gamma, mi, pK)*TComplex::Conjugate(I_R(Eg, costh, hel_phi, mf, pK));
//    }
//  }
//  sum *= TComplex(sq(dx),0);
  /* fast integration (available only when t=t_min) */
  /* get costh at which, amp. is infinity */
  double E_k = TMath::Sqrt(sq(r) + sq(Mk));
  double E_phi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));
  double magic_cos = (2*E_k*E_phi - sq(Mphi)) / (2*p_phi*r);
  if (magic_cos < 1.-dx/2) {
    for (double cosK=magic_cos+dx/2; cosK<1.; cosK+=dx) {
      pK = TLorentzVector(r*TMath::Sqrt(1-sq(cosK)), 0., r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));
      sum += I_L(Eg, hel_gamma, mi, pK)*TComplex::Conjugate(I_R(Eg, costh, hel_phi, mf, pK));
    }
    for (double cosK=magic_cos-dx/2; cosK>-1.; cosK-=dx) {
      pK = TLorentzVector(r*TMath::Sqrt(1-sq(cosK)), 0., r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));
      sum += I_L(Eg, hel_gamma, mi, pK)*TComplex::Conjugate(I_R(Eg, costh, hel_phi, mf, pK));
    }
  } else {
    for (double cosK=-1.; cosK<1.; cosK+=dx) {
      pK = TLorentzVector(r*TMath::Sqrt(1-sq(cosK)), 0., r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));
      sum += I_L(Eg, hel_gamma, mi, pK)*TComplex::Conjugate(I_R(Eg, costh, hel_phi, mf, pK));
    }
  }
  sum *= TComplex(2*pi*dx,0);

  return (sum * TComplex(0, r / (W*32*sq(pi))));
}

double dsigma_dt_R(double Eg, double costh) {
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          sum += sq(TComplex::Abs(ImM(Eg,costh,hel_g,hel_V,mi,mf)));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}
