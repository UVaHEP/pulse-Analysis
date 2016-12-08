#ifndef DARKPEAKER_H
#define DARKPEAKER_H

#include "TH1F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TCanvas.h"

//#define TSPECTFLOAT Float_t  // ROOT 5
#define TSPECTFLOAT Double_t // ROOT 6


//class DarkPeaker : public TSpectrum { 
class DarkPeaker{ 
public:
  DarkPeaker(double peThreshold=-1e20);
  ~DarkPeaker(){;}
  void Reset();
  void SetBuffer(TH1F *newbuf, double sampleTime);
  int AnalyzePeaks();
  void FindNoise();
  void FindBackground();
  TH1F* GetBackground();
  int GetNPeaks();
  double CalcDarkRate(); /// simple peak rate in MHz
  void SetSearchThreshold(double peThreshold) {_peThreshold = peThreshold;}
  TH1F* GetHdist() {return (TH1F*)(hdist->Clone());}
  bool haveAnalysis;
  TSPECTFLOAT* GetTSpectrumX();  // get arrays from tspectrum
  TSPECTFLOAT* GetTSpectrumY();
  void GetPoint(int i, double &x, double &y) const;  // return time ordered peaks data
  
 private:
  TSpectrum *tspectrum;
  TH1F *buf;
  TH1F *hbkg;
  TH1F *hdist;  // frequency distribution of ADC samples
  TH1F *hscan;  // count of samples above threshold
  TF1  *tfNoise;
  int npeaks;
  double _peThreshold;
  double snglPeak;
  double *peaksX;
  double *peaksY;
  TSPECTFLOAT *bkgCorrectedY;
  TSPECTFLOAT *deltaT;
  double dT;
  TCanvas *tcD;
};


TString getoutput(TString cmd);
void BinLogX(TH1*h);

#endif
