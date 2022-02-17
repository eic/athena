void drawHits() {
  TFile *f = new TFile("test.root","READ");
  TTree *t = (TTree*) f->Get("events");
  TCanvas *c = new TCanvas();
  t->Draw("DRICHHits.position.x");
  c->SaveAs("test.png");
  f->Close();
};
