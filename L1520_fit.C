#include "res2.C"

double dsigma_dcos(double Eg, double costh, int sign1, int sign2, int sign3, int sign4,double L1, double L2) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector k2(r*TMath::Sqrt(1-sq(costh)), 0., r*costh, TMath::Sqrt(sq(r)+sq(Mk)));
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
      sum += sq(TComplex::Abs(M_L(Eg, hel_g, mi, k2,sign1,sign2,sign3,sign4,L1,L2)));
    }
  }
  double p_i = (sq(W)-sq(Mp))/(2*W);
  sum *= (1./(64*sq(pi)*sq(W)))*(r/p_i)*(1./4.);
  return (hbarc2*sum*1000.);  // [nbarn]
}

double sigma(double Eg, int sign1, int sign2, int sign3, int sign4,double L1, double L2) {
  double dx = 0.05;
  double sum = 0.0;
  for (double cos=-1.; cos<1.; cos+=dx) {
    sum += dsigma_dcos(Eg, cos, sign1, sign2, sign3, sign4,L1,L2);
  }
  return (2*pi*sum*dx);
}

void L1520_fit() {
  TGraph *tg0 = new TGraph("fig3.txt");

  TCanvas *c1 = new TCanvas();
  c1->SetGrid();
  TH1 *frame1 = c1->DrawFrame(1.6,0.,5,2000.);
  TF1 *cs0 = new TF1("cs0", "sigma(x, 1, -1, 1, -1, [0], [1])", 1.7, 5.);
  tg0->Draw("l");
  cs0->SetParameter(0,1.0);
  cs0->SetParameter(1,1.0);
  cs0->SetParLimits(0,0.5,2.0);
  cs0->SetParLimits(1,0.5,2.0);
  tg0->Fit("cs0","N");
  cs0->Draw("same");
}
