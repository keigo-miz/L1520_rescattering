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

void show() {
  const Int_t ne = 9;
  TGraphErrors *cs[ne];
  for (int i=0; i<ne; i++) {
    cs[i] = new TGraphErrors(Form("dat/cs%d.dat",i+1));
    cs[i]->SetMarkerStyle(20);
  }

  /* mibe-san */
  TGraphErrors *mibe[4];
  for (int i=0; i<4; i++) {
    mibe[i] = new TGraphErrors(Form("dat/mibe%d.dat",i+1));
    mibe[i]->SetMarkerStyle(24);
  }

  TLegend *leg = new TLegend(0.25, 0.3, 0.8, 0.75);
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->AddEntry(cs[0], "This work", "lp");
  leg->AddEntry(mibe[0], "LEPS (2005)", "lp");

  TF1 *ff[ne];
  ff[0] = new TF1("ff0","dsigma_dt(1.77,x)",-0.6,0.0);
  ff[1] = new TF1("ff1","dsigma_dt(1.97,x)",-0.6,0.0);
  ff[2] = new TF1("ff2","dsigma_dt(2.17,x)",-0.6,0.0);
  ff[3] = new TF1("ff3","dsigma_dt(2.32,x)",-0.6,0.0);
  ff[4] = new TF1("ff4","dsigma_dt(2.42,x)",-0.6,0.0);
  ff[5] = new TF1("ff5","dsigma_dt(2.52,x)",-0.6,0.0);
  ff[6] = new TF1("ff6","dsigma_dt(2.62,x)",-0.6,0.0);
  ff[7] = new TF1("ff7","dsigma_dt(2.72,x)",-0.6,0.0);
  ff[8] = new TF1("ff8","dsigma_dt(2.82,x)",-0.6,0.0);

  TCanvas *c1 = new TCanvas("c","c",800,800);
  c1->SetBottomMargin(0.18);
  c1->SetTopMargin(0.02);
  c1->SetLeftMargin(0.18);
  c1->SetRightMargin(0.02);
  c1->SetFillStyle(4000);
  c1->Divide(3,3,0,0);
  TH1 *frame1[ne];
  for (int i=0; i<ne; i++) {
    c1->cd(i+1)->SetFillStyle(4000);
    c1->cd(i+1)->SetFrameFillStyle(4000);
    if (i<6) {
      frame1[i] = c1->cd(i+1)->DrawFrame(-0.59999,0.0001,-0.0001,1.19999);
    } else {
      frame1[i] = c1->cd(i+1)->DrawFrame(-0.59999,0.,-0.0001,1.19999);
    }
    frame1[i]->GetXaxis()->SetNdivisions(503);
    frame1[i]->GetXaxis()->SetLabelSize(0.);
    frame1[i]->GetYaxis()->SetLabelSize(0.);
    ff[i]->Draw("same");
    if (i<4) mibe[i]->Draw("p");
    cs[i]->Draw("p");
    if (i==0) leg->Draw();
  }
  frame1[6]->GetXaxis()->SetLabelSize(0.08);
  frame1[6]->GetXaxis()->SetLabelOffset(0.01);
  frame1[6]->GetYaxis()->SetLabelSize(0.08);
  frame1[6]->GetYaxis()->SetLabelOffset(0.01);
  frame1[6]->GetXaxis()->SetTitle("t #minus t_{min} (GeV^{2})");
  frame1[6]->GetXaxis()->CenterTitle(kTRUE);
  frame1[6]->GetXaxis()->SetTitleSize(0.1);
  frame1[6]->GetXaxis()->SetTitleOffset(0.9);
  frame1[6]->GetYaxis()->SetTitle("d#sigma/dt (#mub/GeV^{2})");
  frame1[6]->GetYaxis()->CenterTitle(kTRUE);
  frame1[6]->GetYaxis()->SetTitleSize(0.1);
  frame1[6]->GetYaxis()->SetTitleOffset(1.);

  frame1[0]->GetYaxis()->SetLabelSize(0.1);
  frame1[0]->GetYaxis()->SetLabelOffset(0.01);
  frame1[3]->GetYaxis()->SetLabelSize(0.1);
  frame1[3]->GetYaxis()->SetLabelOffset(0.01);
  frame1[7]->GetXaxis()->SetLabelSize(0.085);
//  frame1[7]->GetXaxis()->SetLabelOffset(0.08);
  frame1[8]->GetXaxis()->SetLabelSize(0.085);
//  frame1[8]->GetXaxis()->SetLabelOffset(0.08);

  /* pad0 */
  TLatex *latex = new TLatex();
  c1->cd();
  TPad *pad0 = new TPad("pad0","",0.,0.,1.,1.);
  pad0->SetFillStyle(4000);
  pad0->Draw();
  pad0->SetNumber(10);
  c1->cd(10);
  latex->SetTextSize(0.04);
  for (int i=0; i<ne; i++) {
    double xx = 0.10 + (i%3)*0.31;
    double yy = 0.95 - (i/3)*0.31;
    double elow, ehigh;
    if (i<3) {
      elow  = 1.67 + 0.2*i;
      ehigh = 1.87 + 0.2*i;
    } else {
      elow  = 1.97 + 0.1*i;
      ehigh = 2.07 + 0.1*i;
    }
    latex->DrawLatexNDC(xx,yy,Form("%4.2lf < E_{#gamma} < %4.2lf",elow,ehigh));
  }
  c1->Print("dsigdt.pdf");
}
