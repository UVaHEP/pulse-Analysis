#ifndef DARKPEAKER_H
#define DARKPEAKER_H

#include "TH1F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TString.h"


class DarkPeaker {
public:
  DarkPeaker();
  ~DarkPeaker(){;}
  void SetBuffer(TH1F *newbuf, double sampleTime);
  int AnalyzePeaks();
  void FindNoise();
  void FindBackground();
  TH1F* GetBackground();
  int GetNPeaks();
  double CalcDarkRate();
  Double_t* GetPositionX();
  Double_t* GetPositionY();
  Double_t* GetBkgdCorrectedY();
  Double_t* GetDeltaX();
private:
  TSpectrum *tspectrum;
  TH1F *buf;
  TH1F *hbkg;
  TH1F *hdist;
  TF1  *tfNoise;
  int npeaks;
  double snglPeak;
  Double_t *bkgCorrectedY;
  Double_t *deltaX;
  double dT;
};


TString getoutput(TString cmd);


#endif
