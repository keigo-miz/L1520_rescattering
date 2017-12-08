#include "CS_SDM.C"

void draw() {
  TGraphErrors *tg1 = new TGraphErrors("dat/cpar1.dat");
  TGraphErrors *tg2 = new TGraphErrors("dat/cpar_mibe.dat");
  tg1->SetMarkerStyle(24);
  tg2->SetMarkerStyle(25);

  TF1 *f1 = new TF1("f1","dsigma_dt(x,1.,[0],[1],[2],[3],[4])",  1.573,3.0);
  TF1 *f2 = new TF1("f2","dsigma_dt_R(x,1.,[0],[1],[2],[3])",1.6,3.0);
  f2->SetLineStyle(2);
  f2->SetLineColor(3);

//  double pp[5] = {0.621686, 20., 20., 2.79486, 0.507631};
//  double pp[5] = {0.583984, 2., 8., 4., 0.6};
//  double pp[5] = {0.62138, 5., 15., 2.8139, 0.503733};
//  double pp[5] = {0.616674, 3., 10., 3, 0.516387};
  double pp[5] = {0.617702, 5., 15., 2.85635, 0.529565};
  f1->SetParameters(pp[0],pp[1],pp[2],pp[3],pp[4]);
  f2->SetParameters(pp[1],pp[2],pp[3],pp[4]);

  TCanvas *c1 = new TCanvas();
  TH1 *frame0 = c1->DrawFrame(1.5,0,3.,1.6);
  tg1->Draw("p");
  tg2->Draw("p");
  f1->Draw("same");
  f2->Draw("same");
}
