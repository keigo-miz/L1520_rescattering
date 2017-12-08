#ifndef spinor_mat_h_
#define spinor_mat_h_

#include <cmath>
#define sq(x) ((x)*(x))

#define pi (M_PI)
#define Mp (0.938272)  // proton mass
#define Mphi (1.01946)  // phi mass
#define Mpi0 (0.1349770)  // pi0 mass
#define Meta (0.547862)  // eta mass
#define hbarc2 (389.379304)  // GeV^2 ubarn

/* spinor */
class Spi {
  public:
    TComplex e[4];  // should be private ..
    Spi();
    Spi(TComplex a0, TComplex a1, TComplex a2, TComplex a3);
    void disp();
    Spi operator+ (Spi x);
    Spi operator- (Spi x);
    Spi operator* (TComplex k);
};

Spi::Spi() {
  e[0] = TComplex(0,0); e[1] = TComplex(0,0); e[2] = TComplex(0,0); e[3] = TComplex(0,0);
}

Spi::Spi(TComplex a0, TComplex a1, TComplex a2, TComplex a3) {
  e[0] = a0; e[1] = a1; e[2] = a2; e[3] = a3;
}

void Spi::disp() {
  std::cout << "(" << e[0].Re() << "," << e[0].Im() << ") "
            << "(" << e[1].Re() << "," << e[1].Im() << ") "
            << "(" << e[2].Re() << "," << e[2].Im() << ") "
            << "(" << e[3].Re() << "," << e[3].Im() << ") "
            << std::endl;
}

Spi Spi::operator+ (Spi x) {
  Spi y = Spi();
  for (int i=0; i<4; i++) { y.e[i] = e[i] + x.e[i]; }
  return y;
}

Spi Spi::operator- (Spi x) {
  Spi y = Spi();
  for (int i=0; i<4; i++) { y.e[i] = e[i] - x.e[i]; }
  return y;
}

Spi Spi::operator* (TComplex k) {
  Spi y = Spi();
  for (int i=0; i<4; i++) { y.e[i] = k*e[i]; }
  return y;
}

/* 4x4 matrix (element = complex number) */
class Mat {
  private:
    TComplex e[4][4];
  public:
    Mat();
    Mat(TComplex a00, TComplex a01, TComplex a02, TComplex a03,
        TComplex a10, TComplex a11, TComplex a12, TComplex a13,
        TComplex a20, TComplex a21, TComplex a22, TComplex a23,
        TComplex a30, TComplex a31, TComplex a32, TComplex a33);
    void disp();
    Mat operator+ (Mat x);
    Mat operator- (Mat x);
    Mat operator* (double k);
    Mat operator* (TComplex k);
    Mat operator* (Mat x);
    Spi operator* (Spi x);
};

Mat::Mat() {
  e[0][0] = TComplex(0,0); e[0][1] = TComplex(0,0); e[0][2] = TComplex(0,0); e[0][3] = TComplex(0,0);
  e[1][0] = TComplex(0,0); e[1][1] = TComplex(0,0); e[1][2] = TComplex(0,0); e[1][3] = TComplex(0,0);
  e[2][0] = TComplex(0,0); e[2][1] = TComplex(0,0); e[2][2] = TComplex(0,0); e[2][3] = TComplex(0,0);
  e[3][0] = TComplex(0,0); e[3][1] = TComplex(0,0); e[3][2] = TComplex(0,0); e[3][3] = TComplex(0,0);
}

Mat::Mat(TComplex a00, TComplex a01, TComplex a02, TComplex a03,
         TComplex a10, TComplex a11, TComplex a12, TComplex a13,
         TComplex a20, TComplex a21, TComplex a22, TComplex a23,
         TComplex a30, TComplex a31, TComplex a32, TComplex a33) {
  e[0][0] = a00; e[0][1] = a01; e[0][2] = a02; e[0][3] = a03;
  e[1][0] = a10; e[1][1] = a11; e[1][2] = a12; e[1][3] = a13;
  e[2][0] = a20; e[2][1] = a21; e[2][2] = a22; e[2][3] = a23;
  e[3][0] = a30; e[3][1] = a31; e[3][2] = a32; e[3][3] = a33;
}

void Mat::disp() {
  for (int i=0; i<4; i++) {
    std::cout << "(" << e[i][0].Re() << "," << e[i][0].Im() << ") "
              << "(" << e[i][1].Re() << "," << e[i][1].Im() << ") "
              << "(" << e[i][2].Re() << "," << e[i][2].Im() << ") "
              << "(" << e[i][3].Re() << "," << e[i][3].Im() << ") "
              << std::endl;
  }
}

Mat Mat::operator+ (Mat x) {
  Mat y = Mat();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) { y.e[i][j] = e[i][j] + x.e[i][j]; }
  }
  return y;
}

Mat Mat::operator- (Mat x) {
  Mat y = Mat();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) { y.e[i][j] = e[i][j] - x.e[i][j]; }
  }
  return y;
}

Mat Mat::operator* (double k) {
  Mat y = Mat();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) { y.e[i][j] = k*e[i][j]; }
  }
  return y;
}

Mat Mat::operator* (TComplex k) {
  Mat y = Mat();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) { y.e[i][j] = k*e[i][j]; }
  }
  return y;
}

Mat Mat::operator* (Mat x) {
  Mat y = Mat();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      for (int k=0; k<4; k++) { y.e[i][j] += e[i][k]*x.e[k][j]; }
    }
  }
  return y;
}

Spi Mat::operator* (Spi x) {
  Spi y = Spi();
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) { y.e[i] += e[i][j]*x.e[j]; }
  }
  return y;
}

/* gamma matrices */
Mat G[4] = {
  Mat(TComplex(1,0), TComplex(0,0), TComplex(0,0), TComplex(0,0),
      TComplex(0,0), TComplex(1,0), TComplex(0,0), TComplex(0,0),
      TComplex(0,0), TComplex(0,0), TComplex(-1,0), TComplex(0,0),
      TComplex(0,0), TComplex(0,0), TComplex(0,0), TComplex(-1,0)),
  Mat(TComplex(0,0), TComplex(0,0), TComplex(0,0), TComplex(1,0),
      TComplex(0,0), TComplex(0,0), TComplex(1,0), TComplex(0,0),
      TComplex(0,0), TComplex(-1,0), TComplex(0,0), TComplex(0,0),
      TComplex(-1,0), TComplex(0,0), TComplex(0,0), TComplex(0,0)),
  Mat(TComplex(0,0), TComplex(0,0), TComplex(0,0), TComplex(0,-1),
      TComplex(0,0), TComplex(0,0), TComplex(0,1), TComplex(0,0),
      TComplex(0,0), TComplex(0,1), TComplex(0,0), TComplex(0,0),
      TComplex(0,-1), TComplex(0,0), TComplex(0,0), TComplex(0,0)),
  Mat(TComplex(0,0), TComplex(0,0), TComplex(1,0), TComplex(0,0),
      TComplex(0,0), TComplex(0,0), TComplex(0,0), TComplex(-1,0),
      TComplex(-1,0), TComplex(0,0), TComplex(0,0), TComplex(0,0),
      TComplex(0,0), TComplex(1,0), TComplex(0,0), TComplex(0,0))
};

/* DiracAdjoint(x)*y */
TComplex barDot(Spi x, Spi y) {
  Spi z = G[0]*y;
  TComplex sum;
  for (int i=0; i<4; i++) sum += TComplex::Conjugate(x.e[i])*z.e[i];
  return sum;
}

/* Dirac spinor of the nucleon */
Spi u(TLorentzVector p, int m) {
  TComplex chi_r[2];
  if (m==1) {
   chi_r[0] = TComplex(1,0);
   chi_r[1] = TComplex(0,0);
  } else if (m==-1) {
   chi_r[0] = TComplex(0,0);
   chi_r[1] = TComplex(1,0);
  }
  TComplex u2 = (TComplex(p(0),0)*chi_r[1] + TComplex(0,-p(1))*chi_r[1] + TComplex(p(2),0)*chi_r[0]) / (TComplex(p(3)+p.M(),0));
  TComplex u3 = (TComplex(p(0),0)*chi_r[0] + TComplex(0,p(1))*chi_r[0] + TComplex(-p(2),0)*chi_r[1]) / (TComplex(p(3)+p.M(),0));

  TComplex N = TComplex(TMath::Sqrt(p.E()+p.M()),0);  // normalization factor
  return (Spi(chi_r[0],chi_r[1],u2,u3) * N);
}

/* polarization vector */
void FillPolVector(TLorentzVector p, int hel, TComplex *eps) {
  if (hel == 0) {
    /* helicity frame */
    // eps[0] = TComplex((p.Vect()).Mag()/p.M(),0);
    // eps[1] = TComplex(p.Px() / (p.M() * p.Beta()),0);
    // eps[2] = TComplex(0,0);
    // eps[3] = TComplex(p.Gamma()*p.CosTheta(),0);
    /* GJ frame */
    eps[0] = TComplex(0,0);
    eps[1] = TComplex(p.Px() / (p.Vect().Mag()),0);
    eps[2] = TComplex(0,0);
    eps[3] = TComplex(p.CosTheta(),0);
  } else {
    eps[0] = TComplex(0,0);
    eps[1] = TComplex((-hel)*p.CosTheta()/TMath::Sqrt(2),0);
    eps[2] = TComplex(0,-1/TMath::Sqrt(2));
    eps[3] = TComplex((hel/TMath::Sqrt(2))*p.Px()/(p.Vect()).Mag(),0);
  }
}

#endif
