#include "res.C"

double dsigma_dcos(double Eg, double costh) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector k2(r*TMath::Sqrt(1-sq(costh)), 0., r*costh, TMath::Sqrt(sq(r)+sq(Mk)));
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
      sum += sq(TComplex::Abs(M_L(Eg, hel_g, mi, k2)));
    }
  }
  double p_i = (sq(W)-sq(Mp))/(2*W);
  sum *= (1./(64*sq(pi)*sq(W)))*(r/p_i)*(1./4.);
  return (hbarc2*sum*1000.);  // [nbarn]
}

double sigma(double Eg) {
  double dx = 0.05;
  double sum = 0.0;
  for (double cos=-1.; cos<1.; cos+=dx) {
    sum += dsigma_dcos(Eg, cos);
  }
  return (2*pi*sum*dx);
}

void L1520() {
  TGraph *tg0 = new TGraph("fig3.txt");

  TCanvas *c1 = new TCanvas();
  c1->SetGrid();
  TH1 *frame1 = c1->DrawFrame(1.6,0.,5,1500.);
  TF1 *cs0 = new TF1("cs0", "sigma(x)", 1.7, 5.);
  cs0->Draw("same");
  tg0->Draw("l");
}
