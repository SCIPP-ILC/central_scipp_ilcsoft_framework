void alphaBeta(){
  TCanvas c;
  c.SetLogy();
  alpha->Draw();
  beta->SetLineColor(2);
  beta->Draw("SAME");
  auto l=new TLegend(0,.2,0,.2);
  l->AddEntry(alpha);
  l->AddEntry(beta);
}
