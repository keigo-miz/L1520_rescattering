#include "res.C"

TComplex preL(int iC, double Eg, int hel_gamma, int mi, double cosK, double phi) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector pK(r*TMath::Sqrt(1-sq(cosK))*TMath::Cos(phi), r*TMath::Sqrt(1-sq(cosK))*TMath::Sin(phi), r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));

  TComplex ret;
  if (iC == 1) {  // s-channel
    ret = M_L_s(Eg, hel_gamma, mi, pK);
  } else if (iC == 2) {  // t-channel
    ret = M_L_t(Eg, hel_gamma, mi, pK);
  } else if (iC == 3) {  // contact term
    ret = M_L_c(Eg, hel_gamma, mi, pK);
  }
  return ret;
}

TComplex preR(int iC, double Eg, double costh, int hel_phi, int mf, double cosK, double phi) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  TLorentzVector pK(r*TMath::Sqrt(1-sq(cosK))*TMath::Cos(phi), r*TMath::Sqrt(1-sq(cosK))*TMath::Sin(phi), r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));

  TComplex ret;
  if (iC == 1) {  // s-channel
    ret = M_R_s(Eg, costh, hel_phi, mf, pK);
  } else if (iC == 2) {  // t-channel
    ret = M_R_t(Eg, costh, hel_phi, mf, pK);
  } else if (iC == 3) {  // contact term
    ret = M_R_c(Eg, costh, hel_phi, mf, pK);
  }
  return ret;
}

double calc_costh(double Eg) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);
  double E_k = TMath::Sqrt(sq(r) + sq(Mk));
  double E_phi = (sq(W) + sq(Mphi) - sq(Mp)) / (2*W);
  double p_phi = TMath::Sqrt(sq(E_phi) - sq(Mphi));

  double magic_cos = (2*E_k*E_phi - sq(Mphi)) / (2*p_phi*r);
  return magic_cos;
}

void show_amp() {
  TF1 *f0 = new TF1("f0","preL(1, 2.2, +1,-1, x,0.).Im()", -1,1);
  TF1 *f1 = new TF1("f1","preL(2, 2.2, +1,+1, x,0.).Im()", -1,1);
  TF1 *f2 = new TF1("f2","preL(2, 2.2, +1,-1, x,0.).Im()", -1,1);
  TF1 *f3 = new TF1("f3","preL(3, 2.2, +1,+1, x,0.).Im()", -1,1);
  TF1 *f4 = new TF1("f4","preL(3, 2.2, +1,-1, x,0.).Im()", -1,1);
  f1->SetLineColor(3);
  f2->SetLineColor(4);
  f3->SetLineColor(5);
  f4->SetLineColor(6);

  TLegend *leg1 = new TLegend(0.6,0.6,0.9,0.9);
  leg1->SetFillStyle(0); leg1->SetBorderSize(0);
  leg1->AddEntry(f0, "s (+1,-1)", "l");
  leg1->AddEntry(f1, "t (+1,+1)", "l");
  leg1->AddEntry(f2, "t (+1,-1)", "l");
  leg1->AddEntry(f3, "c (+1,+1)", "l");
  leg1->AddEntry(f4, "c (+1,-1)", "l");

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);
  TH1 *frame1 = c1->cd(1)->DrawFrame(-1, -9, 1, 6);
  f0->Draw("same");
  f1->Draw("same");
  f2->Draw("same");
  f3->Draw("same");
  f4->Draw("same");
  leg1->Draw();

  TF1 *g0 = new TF1("g0","preR(1, 2.2,1.0, +1,-1, x,0.).Im()", -1,1);
  TF1 *g1 = new TF1("g1","preR(1, 2.2,1.0, +0,+1, x,0.).Im()", -1,1);
  TF1 *g2 = new TF1("g2","preR(2, 2.2,1.0, +1,+1, x,0.).Im()", -1,1);
  TF1 *g3 = new TF1("g3","preR(2, 2.2,1.0, +1,-1, x,0.).Im()", -1,1);
  TF1 *g4 = new TF1("g4","preR(2, 2.2,1.0, +0,+1, x,0.).Im()", -1,1);
  TF1 *g5 = new TF1("g5","preR(3, 2.2,1.0, +1,+1, x,0.).Im()", -1,1);
  TF1 *g6 = new TF1("g6","preR(3, 2.2,1.0, +1,-1, x,0.).Im()", -1,1);
  TF1 *g7 = new TF1("g7","preR(3, 2.2,1.0, +0,+1, x,0.).Im()", -1,1);
  g1->SetLineColor(3);
  g2->SetLineColor(4);
  g3->SetLineColor(5);
  g4->SetLineColor(6);
  g5->SetLineColor(7);
  g6->SetLineColor(8);
  g7->SetLineColor(9);
  g4->SetLineStyle(2);
  g5->SetLineStyle(2);
  g6->SetLineStyle(2);
  g7->SetLineStyle(2);

  TLegend *leg2 = new TLegend(0.6,0.6,0.9,0.9);
  leg2->SetFillStyle(0); leg2->SetBorderSize(0);
  leg2->AddEntry(g0, "s (+1,-1)", "l");
  leg2->AddEntry(g1, "s (+0,+1)", "l");
  leg2->AddEntry(g2, "t (+1,+1)", "l");
  leg2->AddEntry(g3, "t (+1,-1)", "l");
  leg2->AddEntry(g4, "t (+0,+1)", "l");
  leg2->AddEntry(g5, "c (+1,+1)", "l");
  leg2->AddEntry(g6, "c (+1,-1)", "l");
  leg2->AddEntry(g7, "c (+0,+1)", "l");

  TH1 *frame2 = c1->cd(2)->DrawFrame(-1, -80, 1, 80);
  g0->Draw("same");
  g1->Draw("same");
  g2->Draw("same");
  g3->Draw("same");
  g4->Draw("same");
  g5->Draw("same");
  g6->Draw("same");
  g7->Draw("same");
  leg2->Draw();

  TCanvas *c2 = new TCanvas();
  TH1 *frame3 = c2->DrawFrame(0.8, -1000, 1, 1000);
  TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9);
  leg3->SetFillStyle(0); leg3->SetBorderSize(0);
  leg3->AddEntry(g2, "t (+1,+1)", "l");
  leg3->AddEntry(g3, "t (+1,-1)", "l");
  leg3->AddEntry(g4, "t (+0,+1)", "l");
  g2->Draw("same");
  g3->Draw("same");
  g4->Draw("same");
  leg3->Draw();

  TF1 *fcos = new TF1("fcos", "calc_costh(x)", 1.7, 3.);
  TCanvas *c3 = new TCanvas();
  c3->SetGrid();
  fcos->Draw();
  std::cout << "magic cos (Eg=2.0) = " << calc_costh(2.0) << std::endl;
  std::cout << "magic cos (Eg=2.5) = " << calc_costh(2.5) << std::endl;
  std::cout << "magic cos (Eg=3.0) = " << calc_costh(3.0) << std::endl;
}
