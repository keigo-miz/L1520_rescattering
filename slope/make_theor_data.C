#include "CS_SDM.C"

double dsdt(double Eg, double tp) {
  double W = TMath::Sqrt(2*Eg*Mp + sq(Mp));
  double EgCM = (sq(W) - sq(Mp)) / (2*W);
  double Ephi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double Pphi = TMath::Sqrt((sq(W) - sq(Mphi+Mp))*(sq(W) - sq(Mphi-Mp))) / (2*W);
  double tmin = 2*EgCM*Pphi - 2*EgCM*Ephi + sq(Mphi);
  double t = tp + tmin;
  double cost = (t + 2*EgCM*Ephi - sq(Mphi)) / (2*EgCM*Pphi);

  double ret = dsigma_dt(Eg,cost);
  return ret;
}

void make_theor_data(int bin = 1) {
  double bin2Eg[9] = {1.77,1.97,2.17,2.32,2.42,2.52,2.62,2.72,2.82};
  double Eg = bin2Eg[bin-1];
  std::cout << "(bin,Eg): (" << bin << ", " << Eg << ")" << std::endl;
  double t = -0.000001;
  while (t>-0.1) {
    std::cout << t << "  " << dsdt(Eg,t) << std::endl;
    t += -0.01;
  }
  while (t>-0.61) {
    std::cout << t << "  " << dsdt(Eg,t) << std::endl;
    t += -0.05;
  }
}
