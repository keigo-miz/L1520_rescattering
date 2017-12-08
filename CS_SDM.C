#include "pomeron.C"
#include "meson.C"
#include "res.C"

double dsigma_dt(double Eg, double costh) {
  TComplex I;
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          I = TComplex(0,0);
          I += Pomeron_amp(Eg,costh,hel_g,hel_V,mi,mf);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,1);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,2);
          I += iImI(Eg,costh,hel_g,hel_V,mi,mf);
          sum += sq(TComplex::Abs(I));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}

double dsigma_dt_T(double Eg, double costh) {
  TComplex I;
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          I = TComplex(0,0);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,1);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,2);
          sum += sq(TComplex::Abs(I));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}

double dsigma_dt_P(double Eg, double costh) {
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          sum += sq(TComplex::Abs(Pomeron_amp(Eg,costh,hel_g,hel_V,mi,mf)));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}

TComplex SDM_N(double Eg, double costh) {
  TComplex I;
  TComplex N = TComplex(0,0);
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          I = TComplex(0,0);
          I += Pomeron_amp(Eg,costh,hel_g,hel_V,mi,mf);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,1);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,2);
          I += iImI(Eg,costh,hel_g,hel_V,mi,mf);
          N += I*TComplex::Conjugate(I);
        }
      }
    }
  }
  return N;
}

TComplex rho0(double Eg, double costh, int lambda, int lambdap) {
  TComplex I1, I2;
  TComplex sum = TComplex(0,0);
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        I1 = TComplex(0,0);
        I1 += Pomeron_amp(Eg,costh,hel_g,lambda,mi,mf);
        I1 += Meson_amp(Eg,costh,hel_g,lambda,mi,mf,1);
        I1 += Meson_amp(Eg,costh,hel_g,lambda,mi,mf,2);
        I1 += iImI(Eg,costh,hel_g,lambda,mi,mf);
        I2 = TComplex(0,0);
        I2 += Pomeron_amp(Eg,costh,hel_g,lambdap,mi,mf);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,1);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,2);
        I2 += iImI(Eg,costh,hel_g,lambdap,mi,mf);
        sum += I1*TComplex::Conjugate(I2);
      }
    }
  }
  return (sum / SDM_N(Eg, costh));
}

TComplex rho1(double Eg, double costh, int lambda, int lambdap) {
  TComplex I1, I2;
  TComplex sum = TComplex(0,0);
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        I1 = TComplex(0,0);
        I1 += Pomeron_amp(Eg,costh,-hel_g,lambda,mi,mf);
        I1 += Meson_amp(Eg,costh,-hel_g,lambda,mi,mf,1);
        I1 += Meson_amp(Eg,costh,-hel_g,lambda,mi,mf,2);
        I1 += iImI(Eg,costh,hel_g,lambda,mi,mf);
        I2 = TComplex(0,0);
        I2 += Pomeron_amp(Eg,costh,hel_g,lambdap,mi,mf);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,1);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,2);
        I2 += iImI(Eg,costh,hel_g,lambdap,mi,mf);
        sum += I1*TComplex::Conjugate(I2);
      }
    }
  }
  return (sum / SDM_N(Eg, costh));
}

TComplex rho2(double Eg, double costh, int lambda, int lambdap) {
  TComplex I1, I2;
  TComplex sum = TComplex(0,0);
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        I1 = TComplex(0,0);
        I1 += Pomeron_amp(Eg,costh,-hel_g,lambda,mi,mf);
        I1 += Meson_amp(Eg,costh,-hel_g,lambda,mi,mf,1);
        I1 += Meson_amp(Eg,costh,-hel_g,lambda,mi,mf,2);
        I1 += iImI(Eg,costh,hel_g,lambda,mi,mf);
        I2 = TComplex(0,0);
        I2 += Pomeron_amp(Eg,costh,hel_g,lambdap,mi,mf);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,1);
        I2 += Meson_amp(Eg,costh,hel_g,lambdap,mi,mf,2);
        I2 += iImI(Eg,costh,hel_g,lambdap,mi,mf);
        sum += TComplex(hel_g,0)*I1*TComplex::Conjugate(I2);
      }
    }
  }
  return (TComplex(0,1) * sum / SDM_N(Eg, costh));
}
