#include "../res.C"

double dsigma_dcos(double Eg, double costh) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector k2(r*TMath::Sqrt(1-sq(costh)), 0., r*costh, TMath::Sqrt(sq(r)+sq(Mk)));
  double sum = 0.;
  for (int mi=-1; mi<=+1; mi+=2) {
    for (int hel_g=-1; hel_g<=+1; hel_g+=2) {
      sum += sq(TComplex::Abs(I_L(Eg, hel_g, mi, k2)));
    }
  }
  double p_i = (sq(W)-sq(Mp))/(2*W);
  sum *= (1./(64*sq(pi)*sq(W)))*(r/p_i)*(1./4.);
  return (hbarc2*sum);  // [ubarn]
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
  TGraphErrors *tg0 = new TGraphErrors("dat/clas.dat");
  TGraphErrors *tg1 = new TGraphErrors("dat/saphir.dat");
  TGraphErrors *tg2 = new TGraphErrors("dat/lamp2.dat");
  TGraph *tg3 = new TGraph("dat/fig3.txt");
  tg0->SetMarkerStyle(20);
  tg1->SetMarkerStyle(24);
  tg2->SetMarkerStyle(25);
  tg3->SetLineColor(3);

  TCanvas *c1 = new TCanvas();
  c1->SetFillStyle(4000);
  c1->SetFrameFillStyle(4000);
  TH1 *frame0 = c1->DrawFrame(1.50001,0.,5.0,1.5);
  frame0->SetFillStyle(4000);
  frame0->GetYaxis()->SetNdivisions(503);
  frame0->GetYaxis()->CenterTitle(kTRUE);
  frame0->GetYaxis()->SetTitleSize(0.06);
  frame0->GetYaxis()->SetTitleOffset(0.6);
  frame0->GetYaxis()->SetTitle("#sigma (#mub)");
  frame0->GetXaxis()->CenterTitle(kTRUE);
  frame0->GetXaxis()->SetTitleSize(0.05);
  frame0->GetXaxis()->SetTitleOffset(0.9);
  frame0->GetXaxis()->SetTitle("E_{#gamma} (GeV)");

  /* draw theor curve */
  double Eth = (sq(Mk + MLs) - sq(Mp))/(2*Mp);
  TF1 *cs0 = new TF1("cs0", "sigma(x)", Eth-0.01662, 5.);
  tg3->Draw("l");
  cs0->Draw("same");
  tg0->Draw("p");
  tg1->Draw("p");
  tg2->Draw("p");

  /* legend */
  TLegend *leg = new TLegend(0.6,0.5,0.9,0.9);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(tg0,"CLAS 2013","lp");
  leg->AddEntry(tg1,"SAPHIR 2011","lp");
  leg->AddEntry(tg2,"LAMP2 1980","lp");
  leg->AddEntry(tg3,"Ryu #font[52]{et al.} (2014)","l");
  leg->AddEntry(cs0,"Fit result","lp");
  leg->Draw();

  c1->Print("../../dt/pic/5disc/L1520fit.pdf");
}
