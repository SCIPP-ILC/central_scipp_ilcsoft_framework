void eazy(){
  TCanvas c;
  c.SetLogy();
  az->Draw();
  ez->SetLineColor(2);
  ez->Draw("SAME");
}
