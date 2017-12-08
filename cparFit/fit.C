#include "CS_SDM.C"

void fit() {
  TGraphErrors *tg0 = new TGraphErrors("dat/leps.dat");
  TGraphErrors *tg1 = new TGraphErrors("dat/cpar1.dat");
  TGraphErrors *tg2 = new TGraphErrors("dat/cpar_mibe.dat");
  tg0->SetMarkerStyle(20);
  tg1->SetMarkerStyle(24);
  tg2->SetMarkerStyle(25);

  TF1 *f1 = new TF1("f1","dsigma_dt(x,1.,[0],[1],[2],[3],[4])",  1.573,3.0);
  TF1 *f2 = new TF1("f2","dsigma_dt_R(x,1.,[0],[1],[2],[3])",1.6,3.0);
  f2->SetLineStyle(2);
  f2->SetLineColor(3);
  double pp[5] = {0.62, 2., 8., 4., 0.51};
  f1->SetParameters(pp[0],pp[1],pp[2],pp[3],pp[4]);
//  f2->SetParameters(1,1,1,1);

//  f1->FixParameter(1,2.);
//  f1->FixParameter(2,8.);
//  f1->FixParameter(3,4.);
//  f1->FixParameter(4,0.6);

  f1->SetParLimits(0,0.5,1.);
  /* FF par limits */
  f1->SetParLimits(1,1.5,5);
  f1->SetParLimits(2,6.,15);
  f1->SetParLimits(3,2.,10.);
  f1->SetParLimits(4,0.5,0.8);

  TCanvas *c1 = new TCanvas();
  TH1 *frame0 = c1->DrawFrame(1.5,0,3.,1.6);
  tg0->Draw("p");
  tg0->Fit("f1");
}
