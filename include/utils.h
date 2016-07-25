#ifndef DARKPEAKER_H
#define DARKPEAKER_H

#include "TH1F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TCanvas.h"

#define TSPECTFLOAT Float_t  // ROOT 5
//#define TSPECTFLOAT Double_t // ROOT 6


class DarkPeaker {
public:
  DarkPeaker(double peThreshold=-1e20);
  ~DarkPeaker(){;}
  void SetBuffer(TH1F *newbuf, double sampleTime);
  int AnalyzePeaks();
  void FindNoise();
  void FindBackground();
  TH1F* GetBackground();
  int GetNPeaks();
  double CalcDarkRate();
  void SetSearchThreshold(double peThreshold) {_peThreshold = peThreshold;}
  TSPECTFLOAT* GetPositionX();
  TSPECTFLOAT* GetPositionY();
  TSPECTFLOAT* GetBkgdCorrectedY();
  TSPECTFLOAT* GetDeltaX();
private:
  TSpectrum *tspectrum;
  TH1F *buf;
  TH1F *hbkg;
  TH1F *hdist;  // frequeny distribution of ADC samples
  TH1F *hscan;  // count of samples above threshold
  TF1  *tfNoise;
  int npeaks;
  double _peThreshold;
  double snglPeak;
  TSPECTFLOAT *bkgCorrectedY;
  TSPECTFLOAT *deltaX;
  double dT;
  TCanvas *tcD;
};


TString getoutput(TString cmd);


#endif
