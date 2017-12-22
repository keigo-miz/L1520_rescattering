#include "CS_SDM.C"

double rr1(double Eg) {
  return rho0(Eg, 1., 0, 0).Re();
}

double rr2(double Eg) {
  return 0;
//  return rho0(Eg, 1., 1, -1).Re();
}

double rr3(double Eg) {
  return ((rho1(Eg,1.,1,-1).Re() - rho2(Eg,1.,1,-1).Im())/2.);
}

double rr4(double Eg) {
  return 0;
//  return ((rho1(Eg,1.,1,-1).Re() + rho2(Eg,1.,1,-1).Im())/2.);
}

double rr5(double Eg) {
  return 0;
//  return (2*rho1(Eg,1.,1,1).Re() + rho1(Eg,1.,0,0).Re());
}

void show_SDME() {
  gStyle->SetLabelFont(132,"XY");
  gStyle->SetLabelSize(0.042,"XY");
  gStyle->SetTitleFont(132,"XY");
  double Elow = 1.7;
  double Ehigh = 2.9;
  /* SDME */
  TF1 *rho[5];
  for (int i=0; i<5; i++) {
    rho[i] = new TF1(Form("rho%d",i+1), Form("rr%d(x)",i+1),Elow,Ehigh);
    rho[i]->SetLineColor(2);
  }

  TGraphErrors *tg0[5], *chang[5];
  for (int i=0; i<5; i++) {
    tg0[i] = new TGraphErrors(Form("dat/rho%dVe.dat",i+1));
    tg0[i]->SetMarkerStyle(20);
    chang[i] = new TGraphErrors(Form("dat/chang%d.dat",i+1));
    chang[i]->SetMarkerStyle(24);
  }

  /* title */
  TString title[5] = {"#rho^{0}_{00}",
                      "Re#rho^{0}_{1#minus1}",
                      "#bar{#rho}^{1}_{1#minus1}",
                      "#Delta_{1#minus1}",
                      "2#rho^{1}_{11}#plus#rho^{1}_{00}"};
  double pos[5] = {0.8,0.65,0.8,0.77,0.55};

  /* Mibe-san's */
  double xx[2],yy[2],xe[2],ye[2];
  xx[0] = 2.072; yy[0] = 0.197; xe[0] = 0.; ye[0] = 0.0372;
  xx[1] = 2.272; yy[1] = 0.189; xe[1] = 0.; ye[1] = 0.0247;
  TGraphErrors *mibe = new TGraphErrors(2,xx,yy,xe,ye);
  mibe->SetMarkerStyle(25);

  /* Legend setting */
  TLegend *leg = new TLegend(0.2, 0.5, 0.8, 0.9);
  leg->SetTextFont(132);
  leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->AddEntry(tg0[0], "This work", "lp");
  leg->AddEntry(chang[0], "LEPS (2010)", "lp");

  TCanvas *c1 = new TCanvas();
  c1->SetBottomMargin(0.18);
  c1->SetTopMargin(0.02);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.08);
  c1->SetFillStyle(4000);
  c1->Divide(3,2,0,0);
  TH1 *frame1[5];
  TLatex *latex = new TLatex();
  latex->SetTextFont(132);
  latex->SetTextSize(0.12);
  TLine *line0 = new TLine(Elow,0.,Ehigh,0.);
  line0->SetLineStyle(2);
  for (int i=0; i<5; i++) {
    c1->cd(i+1)->SetFillStyle(4000);
    c1->cd(i+1)->SetFrameFillStyle(4000);
    frame1[i] = c1->cd(i+1)->DrawFrame(Elow,-0.09999,Ehigh,0.49999);
    frame1[i]->GetXaxis()->SetLabelSize(0.);
    frame1[i]->GetYaxis()->SetLabelSize(0.);
//    if (i==2) mibe->Draw("p");
    line0->Draw();
    rho[i]->Draw("same");
    chang[i]->Draw("p");
    tg0[i]->Draw("p");
    if (i<3) latex->DrawLatexNDC(pos[i],0.85,title[i]);
    else     latex->DrawLatexNDC(pos[i],0.87,title[i]);
  }
  c1->cd(1); leg->Draw();
  c1->cd(6)->SetFillStyle(4000);
  c1->cd(6)->SetFrameFillStyle(4000);
  frame1[3]->GetXaxis()->SetTitle("#font[12]{E_{#gamma}} (GeV)");
  frame1[3]->GetXaxis()->CenterTitle(kTRUE);
  frame1[3]->GetXaxis()->SetTitleSize(0.1);
  frame1[3]->GetXaxis()->SetTitleOffset(0.9);

  /* label settings */
  frame1[0]->GetYaxis()->SetLabelSize(0.09);
  frame1[3]->GetXaxis()->SetLabelSize(0.08);
  frame1[3]->GetYaxis()->SetLabelSize(0.08);
  frame1[3]->GetYaxis()->SetLabelOffset(0.01);
  frame1[0]->GetYaxis()->SetLabelOffset(0.01);
  frame1[4]->GetXaxis()->SetLabelSize(0.09);

  /* pad */
  c1->cd();
  TPad *pad0 = new TPad("pad0","",0.,0.,1.,1.);
  pad0->SetFillStyle(4000);
  pad0->Draw();
  pad0->SetNumber(7);
  c1->cd(7);
  latex->SetTextSize(0.042);
  double start = 0.672, diff = 0.2620;
  latex->DrawLatexNDC(start+diff*0/5,0.505,"1.8");
  latex->DrawLatexNDC(0.734,0.505,"2");
  latex->DrawLatexNDC(start+diff*2/5,0.505,"2.2");
  latex->DrawLatexNDC(start+diff*3/5,0.505,"2.4");
  latex->DrawLatexNDC(start+diff*4/5,0.505,"2.6");
  latex->DrawLatexNDC(start+diff*5/5,0.505,"2.8");

  c1->Print("../dt/pic/5disc/rhoVe_rescatt.pdf");
}
