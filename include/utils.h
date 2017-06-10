#ifndef DARKPEAKER_H
#define DARKPEAKER_H

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TCanvas.h"
#include <vector>
using std::vector;

// hack to handle changes between ROOT versions
#if __GNUC__ < 6
#define TSPECTFLOAT Float_t  // ROOT 5
#else
#define TSPECTFLOAT Double_t // ROOT 6
#endif

//class DarkPeaker : public TSpectrum { // alternative design
class DarkPeaker{ 
public:
  DarkPeaker();
  ~DarkPeaker(){;}

  int AnalyzePeaks(double peThreshold);
  double CalcDarkRate(); /// simple peak rate in MHz
  void DumpPeaks();
  void FindNoise();
  void FindBackground();
  /// calc integrals and FWHM for isolated peaks
  /// nLow, nHigh are bin ranges to left,right of peak
  void Integrate(int nLow, int nHigh, bool selectIsolated=true);
  
  TH1F* GetBackground() const {return hbkg;}
  TH1F* GetBuffer() const {return buf;}
  double GetFWHM(int i) const;
  TH1F* GetFWHMhist() const;
  TH1F* GetHdist() {return (TH1F*)(hdist->Clone());}
  TH1F* GetHscan() {return (TH1F*)(hscan->Clone());}
  void GetIntegral(int i, double &x, double &y) const; 
  int GetNIntegrals(){ return pIntegrals.size(); }
  int GetNPeaks() const {return npeaks;}
  void GetPoint(int i, double &x, double &y) const;  // return time ordered peaks data
  double GetSearchThreshold() {return _peThreshold;}
  TSPECTFLOAT* GetTSpectrumX();  // get arrays from tspectrum
  TSPECTFLOAT* GetTSpectrumY();
  
  void Reset();
  void SetBuffer(TH1F *newbuf, double sampleTime);
  void SetSearchThreshold(double peThreshold) {_peThreshold = peThreshold;}

private:
  double CalcFWHM(int i) const; /// calc FWHM for peak i
  
 private:
  bool haveAnalysis;
  TSpectrum *tspectrum;
  TH1F *buf;    // the picoscope buffer
  TH1F *hbkg;   // background fit from TSpectrum
  TH1F *hdist;  // frequency distribution of ADC samples
  TH1F *hscan;  // count of samples above threshold
  TH1F *hFWHM;  // distribution of pulse widths (for no overlapping pulses)
  TF1  *tfNoise;
  int npeaks;
  double _peThreshold;
  double _noiseCut;
  double snglPeak;
  vector<double> peaksX;
  vector<double> peaksY;
  vector<double> bkgCorrectedY;
  vector<double> pIntegrals;
  vector<double> pFWHM;
  TSPECTFLOAT *deltaT;
  double dT;
  TCanvas *tcD;
};

// helper functions
TString getoutput(TString cmd);
void BinLogX(TH1*h);

// fit a double exponential to separate DCR and afterpulsing in Delta_time distribution
void fit2Exp(TH1F *h, TString opt="L,Q");

// improved fit model
void fcnDcrAp(TH1F *h, TString opt="L,Q");
class FitDcrAp {
 public:
  FitDcrAp();
  void Fit(TH1F *h, TString opt="L,Q");
  TF1 *GetExpFit(){return expoFit;}
  TF1 *GetApFit(){return afterPulseFit;}
  TF1 *GetDcrFcn();
  // DCR rate based on sample time in ns
  // HACK - Only works if x-axis is in units of ns
  double GetDCR() {return 1e9/(afterPulseFit->GetParameter(3));}  // DCR in Hz
  double GetAPrate() {return afterPulseFit->GetParameter(1);}
  /*
  double GetDCRerr();

  double GetAPerr();
  */
 private:
  TF1 *expoFit;
  TF1 *afterPulseFit;
  TF1 *dcrFcn;
};




// calculate DCR based on double exponential fit
void calcDcrAp(TH1F *hdTime, double &fitDCR, double &errDCR, double &rateAP, double errAP);


double apCalc2D(TH2F *h, double dcr, double pe, double sigma);

#endif
