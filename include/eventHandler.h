#include "RQ_OBJECT.h"
#include "TFile.h"
#include "TGraph.h" 
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TQObject.h"
#include <iostream> 

class eventHandler : public TQObject {


  public: 
 eventHandler(TCanvas &eventCanvas) : tc(eventCanvas) { }; 
  ~eventHandler(); 
  void print();
  void eventSlot(Int_t event, Int_t x, Int_t y, TObject *selected);
  Double_t threshold() { return _threshold; }; 
  private:
  Double_t _threshold; 
  TCanvas &tc; 

  RQ_OBJECT("eventHandler"); 
}; 
