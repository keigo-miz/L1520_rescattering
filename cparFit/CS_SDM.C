#include "../meson.C"
#include "pomeron.C"
#include "res.C"

double dsigma_dt(double Eg, double costh, double Cp, double n1, double n2, double L1, double L2) {
  TComplex I;
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          I = TComplex(0,0);
          I += Pomeron_amp(Eg,costh,hel_g,hel_V,mi,mf,Cp);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,1);
          I += Meson_amp(Eg,costh,hel_g,hel_V,mi,mf,2);
          I += iImI(Eg,costh,hel_g,hel_V,mi,mf,n1,n2,L1,L2);
          sum += sq(TComplex::Abs(I));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}

double dsigma_dt_P(double Eg, double costh, double Cp) {
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          sum += sq(TComplex::Abs(Pomeron_amp(Eg,costh,hel_g,hel_V,mi,mf,Cp)));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}

double dsigma_dt_R(double Eg, double costh, double n1, double n2, double L1, double L2) {
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int mf=-1; mf<=+1; mf+=2) {
      for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
        for (int hel_V=-1; hel_V<=+1; hel_V++) {
          sum += sq(TComplex::Abs(iImI(Eg,costh,hel_g,hel_V,mi,mf,n1,n2,L1,L2)));
        }
      }
    }
  }
  double den = 64*pi*sq(2*Eg*Mp);
  return (hbarc2*sum/den);
}
