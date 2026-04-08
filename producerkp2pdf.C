{
  TH2F * h2 = new TH2F("h2","h2",100,-150,50,100,-50,100);
  h2->SetStats(0);
  Tout->Draw("RayY:RayX>>h2","DetAmp>100 ","colz", 149425, 0);
  
  c1->SetTitle("All Pads");
  h2->SetTitle("All Pads");
  c1->SetGridx();
  c1->SetGridy();
  
  c1->Print("RKP2.pdf(","Title:All Pads");

  for(int ich=0;ich<128;ich++)
    {
      Tout->Draw("RayY:RayX>>h2",TString::Format("DetAmp>100 && DetX==%d",ich),"colz", 149425, 0);  
      c1->SetTitle(TString::Format("Pad %d",ich));
      h2->SetTitle(TString::Format("Pad %d",ich));
      c1->Print("RKP2.pdf",TString::Format("Title:Pad %d",ich));
    }



  Tout->Draw("RayY:RayX>>h2","DetAmp>100 ","colz", 149425, 0);
  c1->SetTitle("All Pads");
  h2->SetTitle("All Pads");
  c1->Print("RKP2.pdf)","Title:All Pads");
}
