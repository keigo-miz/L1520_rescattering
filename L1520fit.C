#include "res_for_L1520fit.C"

double dsigma_dcos(double Eg, double costh, double n1, double n2, double L1, double L2) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector k2(r*TMath::Sqrt(1-sq(costh)), 0., r*costh, TMath::Sqrt(sq(r)+sq(Mk)));
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
      sum += sq(TComplex::Abs(I_L(Eg, hel_g, mi, k2, n1, n2, L1, L2)));
    }
  }
  double p_i = (sq(W)-sq(Mp))/(2*W);
  sum *= (1./(64*sq(pi)*sq(W)))*(r/p_i)*(1./4.);
  return (hbarc2*sum*1000.);  // [nbarn]
}

double sigma(double Eg, double n1, double n2, double L1, double L2) {
  double dx = 0.05;
  double sum = 0.0;
  for (double cos=-1.; cos<1.; cos+=dx) {
    sum += dsigma_dcos(Eg, cos, n1, n2, L1, L2);
  }
  return (2*pi*sum*dx);
}

void L1520fit() {
  TGraphErrors *tg0 = new TGraphErrors("dat/clas.dat");
  tg0->SetMarkerStyle(20);

  TCanvas *c1 = new TCanvas();
  c1->SetFillStyle(4000);
  c1->SetFrameFillStyle(4000);
  TH1 *frame0 = c1->DrawFrame(1.0,0.,5.0,1200.);
  frame0->SetFillStyle(4000);
  tg0->Draw("p");

  /* draw theor curve */
  double Eth = (sq(Mk + MLs) - sq(Mp))/(2*Mp);
  TF1 *cs0 = new TF1("cs0", "sigma(x,[0],[1],[2],[3])", Eth, 5.);
  cs0->SetParameters(1,1,1.,1.);
  cs0->FixParameter(0,0.6);
  cs0->FixParameter(1,8.0);
  cs0->FixParameter(2,1.157);
//  cs0->SetParLimits(0,0.6,8.);
//  cs0->SetParLimits(1,0.8,8.);
  cs0->SetParLimits(2,0.1,5.);
  cs0->SetParLimits(3,0.1,5.);
  tg0->Fit("cs0");
  TLatex *latex = new TLatex();
  latex->DrawLatexNDC(0.6,0.75,Form("n_{1} = %4.3lf", cs0->GetParameter(0)));
  latex->DrawLatexNDC(0.6,0.7,Form("n_{2} = %4.3lf", cs0->GetParameter(1)));
  latex->DrawLatexNDC(0.6,0.65,Form("#Lambda_{1} = %4.3lf", cs0->GetParameter(2)));
  latex->DrawLatexNDC(0.6,0.6,Form("#Lambda_{2} = %4.3lf", cs0->GetParameter(3)));
}
