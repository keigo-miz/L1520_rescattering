#include "res.C"

TComplex ImM0(double Eg, double costh, int hel_gamma, int hel_phi, int mi, int mf, double cosK) {
  double W = TMath::Sqrt(2*Mp*Eg + sq(Mp));
  double r = TMath::Sqrt((sq(W)-sq(MLs+Mk))*(sq(W)-sq(MLs-Mk)))/(2*W);

  TLorentzVector pK(r*TMath::Sqrt(1-sq(cosK)), 0., r*cosK, TMath::Sqrt(sq(r)+sq(Mk)));
  return M_L(Eg, hel_gamma, mi, pK)*TComplex::Conjugate(M_R(Eg, costh, hel_phi, mf, pK));
}

void show_ImM() {
  TF1 *ImMf0 = new TF1("ImMf0","ImM0(2.2,1., 1,1,  1, 1, x).Re()", -1,1);
  TF1 *ImMf1 = new TF1("ImMf1","ImM0(2.2,1., 1,1,  1,-1, x).Re()", -1,1);
  TF1 *ImMf2 = new TF1("ImMf2","ImM0(2.2,1., 1,1, -1, 1, x).Re()", -1,1);
  TF1 *ImMf3 = new TF1("ImMf3","ImM0(2.2,1., 1,1, -1,-1, x).Re()", -1,1);
  ImMf1->SetLineColor(3);
  ImMf2->SetLineColor(4);
  ImMf3->SetLineColor(5);

  TLegend *leg3 = new TLegend(0.6,0.6,0.9,0.9);
  leg3->SetFillStyle(0); leg3->SetBorderSize(0);
  leg3->AddEntry(ImMf0, "(1,1,  1, 1,)", "l");
  leg3->AddEntry(ImMf1, "(1,1,  1,-1,)", "l");
  leg3->AddEntry(ImMf2, "(1,1, -1, 1,)", "l");
  leg3->AddEntry(ImMf3, "(1,1, -1,-1,)", "l");

  TCanvas *c2 = new TCanvas();
  c2->Divide(1,2);
  TH1 *frame3 = c2->cd(1)->DrawFrame(-1, -40, 1, 40);
  ImMf0->Draw("same");
  ImMf1->Draw("same");
  ImMf2->Draw("same");
  ImMf3->Draw("same");
  leg3->Draw();

  TF1 *ImMg0 = new TF1("ImMg0","ImM0(2.2,1., 1,0,  1, 1, x).Re()", -1,1);
  TF1 *ImMg1 = new TF1("ImMg1","ImM0(2.2,1., 1,0,  1,-1, x).Re()", -1,1);
  TF1 *ImMg2 = new TF1("ImMg2","ImM0(2.2,1., 1,0, -1, 1, x).Re()", -1,1);
  TF1 *ImMg3 = new TF1("ImMg3","ImM0(2.2,1., 1,0, -1,-1, x).Re()", -1,1);
  ImMg1->SetLineColor(3);
  ImMg2->SetLineColor(4);
  ImMg3->SetLineColor(5);

  TLegend *leg4 = new TLegend(0.6,0.6,0.9,0.9);
  leg4->SetFillStyle(0); leg4->SetBorderSize(0);
  leg4->AddEntry(ImMg0, "(1,0,  1, 1)", "l");
  leg4->AddEntry(ImMg1, "(1,0,  1,-1)", "l");
  leg4->AddEntry(ImMg2, "(1,0, -1, 1)", "l");
  leg4->AddEntry(ImMg3, "(1,0, -1,-1)", "l");

  TH1 *frame4 = c2->cd(2)->DrawFrame(-1, -40, 1, 40);
  ImMg0->Draw("same");
  ImMg1->Draw("same");
  ImMg2->Draw("same");
  ImMg3->Draw("same");
  leg4->Draw();
}
