#include "CS_SDM.C"

void show_cpar() {
  gStyle->SetLabelFont(132,"XY");
  gStyle->SetLabelSize(0.042,"XY");
  gStyle->SetTitleFont(132,"XY");
  double Eth = Mphi + sq(Mphi)/(2*Mp);

  /* theoretical curves */
  TF1 *theor_curve1 = new TF1("theor_curve1", "dsigma_dt(x,1.)",   Eth, 3.0);
  TF1 *theor_curve2 = new TF1("theor_curve2", "dsigma_dt_R(x,1.)", Eth, 3.0);
  TF1 *theor_curve3 = new TF1("theor_curve3", "dsigma_dt_P(x,1.)", Eth, 3.0);
  TF1 *theor_curve4 = new TF1("theor_curve4", "dsigma_dt_T(x,1.)", Eth, 3.0);
  theor_curve1->SetLineColor(2);
  theor_curve2->SetLineColor(3);
  theor_curve3->SetLineColor(6);
  theor_curve4->SetLineColor(7);
  theor_curve2->SetLineStyle(2);
  theor_curve3->SetLineStyle(2);
  theor_curve4->SetLineStyle(2);

  /* data */
  TGraphErrors *tgC1 = new TGraphErrors("cparFit/dat/cpar1.dat");
  TGraphErrors *tgC2 = new TGraphErrors("cparFit/dat/cpar_mibe.dat");
  tgC1->SetMarkerStyle(20);
  tgC2->SetMarkerStyle(24);

  /* Legend setting */
  TLegend *leg1 = new TLegend(0.5, 0.2, 0.95, 0.65);
  leg1->SetTextFont(132);
  leg1->SetFillStyle(0); leg1->SetBorderSize(0);
  leg1->AddEntry(tgC1, "This work", "lp");
  leg1->AddEntry(tgC2, "LEPS (2005)", "lp");
  leg1->AddEntry(theor_curve3, "Pomeron", "l");
  leg1->AddEntry(theor_curve4, "#pi^{0}#plus#eta", "l");
  leg1->AddEntry(theor_curve2, "#font[12]{K^{#plus}}#Lambda(1520) rescatt.", "l");
  leg1->AddEntry(theor_curve1, "Total", "l");

  TCanvas *c1 = new TCanvas("c","c",600,400); // default 700x500
  c1->SetFillStyle(4000);
  c1->SetFrameFillStyle(4000);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.08);
  TH1 *frame1 = c1->DrawFrame(1.5,0.,3.,1.3);
  frame1->GetXaxis()->SetTitle("#font[12]{E_{#gamma}} (GeV)");
  frame1->GetXaxis()->CenterTitle(kTRUE);
  frame1->GetXaxis()->SetTitleSize(0.05);
  frame1->GetYaxis()->SetTitle("(#font[12]{d#sigma}/#font[12]{dt})_{#font[12]{t}=#font[12]{t}_{min}} (#mub/GeV^{2})");
  frame1->GetYaxis()->CenterTitle(kTRUE);
  frame1->GetYaxis()->SetTitleSize(0.05);
  /* label settings */
  frame1->GetXaxis()->SetLabelSize(0.05);
  frame1->GetYaxis()->SetLabelSize(0.05);

  theor_curve1->Draw("same");
  theor_curve2->Draw("same");
  theor_curve3->Draw("same");
  theor_curve4->Draw("same");
  tgC2->Draw("p");
  tgC1->Draw("p");
  leg1->Draw();

  c1->Print("cpar_rescatt.pdf");
}
