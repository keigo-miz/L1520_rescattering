#include "CS_SDM.C"

void draw0() {
  TGraphErrors *tg1 = new TGraphErrors("dat/cpar1.dat");
  TGraphErrors *tg2 = new TGraphErrors("dat/cpar_mibe.dat");
  tg1->SetMarkerStyle(24);
  tg2->SetMarkerStyle(25);

  const Int_t n = 5;
  TF1 *ff[n];
  for (int i=0; i<n; i++) {
    ff[i] = new TF1(Form("ff%d",i),"dsigma_dt_R(x,1.,[0],[1],[2],[3])+0.00001",1.8,3.0);
  }

  ff[0]->SetParameters(2.,  8., 4., 0.6);
  ff[1]->SetParameters(2., 10., 4., 0.6);
  ff[2]->SetParameters(1.,  8., 4., 0.6);
  ff[3]->SetParameters(2.,  8., 2., 0.6);
  ff[4]->SetParameters(2.,  8., 1., 0.6);

  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  TH1 *frame0 = c1->DrawFrame(1.5,0.00001,3.,0.3);
  tg1->Draw("p");
  tg2->Draw("p");
  for (int i=0; i<n; i++) {
    ff[i]->SetLineColor(i+1);
    ff[i]->Draw("same");
  }
}
