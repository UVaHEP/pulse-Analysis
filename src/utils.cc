#include "utils.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TAxis.h"
#include <stdlib.h>
#include <pthread.h>
#include <iostream>

using std::cout;
using std::endl;



//DarkPeaker::DarkPeaker(double peThreshold) : TSpectrum(5000) {
DarkPeaker::DarkPeaker(double peThreshold) {
  tspectrum=new TSpectrum(5000,2);
  haveAnalysis=false;
  Reset();
  if (peThreshold) _peThreshold=peThreshold;
  tcD = new TCanvas("tcD","Dark Peak Analysis");
  // x-axis limits for these will be defined in FindNoise()
  hdist=new TH1F("hdist","Background subtracted ADC distribution",100,0,100);
  hscan=new TH1F("hscan","Threshold scan",100,0,100);  // samples above ADC count threshold
}

void DarkPeaker::Reset(){
  // check if we have run the analyze method and reset storage
  if (haveAnalysis){
    delete hbkg;
    peaksX.clear();
    peaksY.clear();
    bkgCorrectedY.clear();
    pIntegrals.clear();
    pFWHM.clear();
    delete deltaT;
    haveAnalysis=false;
  }
  buf=0;
  hbkg=0;
  deltaT=0;
  npeaks=0;
  haveAnalysis=false;
}


void DarkPeaker::SetBuffer(TH1F *newbuf, double sampleTime){
  Reset();
  buf=newbuf;
  dT=sampleTime;
}


void DarkPeaker::DumpPeaks(){
  for (int i=0; i<=npeaks; i++){
    cout << peaksX[i] << " " << peaksY[i] << endl;
  }
}

TSPECTFLOAT* DarkPeaker::GetTSpectrumX(){ return tspectrum->GetPositionX();}
TSPECTFLOAT* DarkPeaker::GetTSpectrumY(){ return tspectrum->GetPositionY();}


int DarkPeaker::AnalyzePeaks(){
  if (!buf) {
    cout << "DarkPeaker::AnalyzePeaks : No buffer set!" << endl;
    return 1;
  }
  // estimate background
  hbkg=(TH1F*)buf->Clone("background");
  hbkg->SetTitle("background");
  FindBackground();

  // estimate noise and single peak height
  FindNoise();

  // find peaks
  double threshold;
  if (_peThreshold>-1e19) threshold = _peThreshold / buf->GetMaximum();
  else threshold = (snglPeak/2) / buf->GetMaximum();
  cout << "snglPeak/2 : " << snglPeak/2 << endl;
  //
  double sigma=2; // this can/should be optimzed
  npeaks=tspectrum->Search(buf,sigma,"nobackground,nomarkov,nodraw",threshold);
  //npeaks=tspectrum->Search(buf,sigma,"nomarkov,nodraw",threshold);
  //int npeaks=tspectrum->Search(buf,sigma,"nomarkov",threshold);
  cout << "Found " << npeaks << " peaks" << endl;
  cout << "Approximate DCR: " << CalcDarkRate() << " MHz" << endl;

  // retrieve the peaks and sort by time order
  Int_t *index=new Int_t[npeaks];
  TSPECTFLOAT *xpeaks = tspectrum->GetPositionX();
  TSPECTFLOAT *ypeaks = tspectrum->GetPositionY();
  TMath::Sort(npeaks, xpeaks, index, kFALSE);  // index sort by timestamp
  for (int i=0;i<npeaks;i++) {
    peaksX.push_back( xpeaks[index[i]] );
    peaksY.push_back( ypeaks[index[i]] );
  }  
  haveAnalysis=true;
  return 0;
}

void DarkPeaker::GetPoint(int i, double &x, double &y) const{
  if (i>=npeaks) {
    cout << "Warning: npoints= " << npeaks
	 << " index: " << i << " requested" << endl;
  }
  x=peaksX[i];
  y=peaksY[i];
}

void DarkPeaker::GetIntegral(int i, double &x, double &y) const{
  if (i>=pIntegrals.size()) {
    cout << "Warning: nIntegrals= " << pIntegrals.size()
	 << " index: " << i << " requested" << endl;
    x=0;
    y=0;
    return;
  }
  x=peaksX[i];
  y=pIntegrals[i];
}

double DarkPeaker::GetFWHM(int i) const{
  if (i>=pFWHM.size()) {
    cout << "Warning: nFWHM= " << pFWHM.size()
	 << " index: " << i << " requested" << endl;
    return 0;
  }
  return pFWHM[i];
}



/// Calculate simple integrals of peaks by adding bkg-subtracted bins in range
/// peak-iLow to peak+iHigh
/// for now we select only isolated peaks (eg. no other peaks in the integration window)
/// we can update this with fits, or something that extends to include overlapping peaks
void DarkPeaker::Integrate(int iLow, int iHigh, bool selectIsolated){
  for (int i=0; i<peaksX.size(); i++){
    double x=peaksX[i];
    int binX=buf->FindBin(x);
    int bLow=binX-iLow;
    int bHigh=binX+iHigh;
    if ( (binX<iLow+1)||(binX+iHigh > buf->GetNbinsX()) ) continue;  // we are too close to an edge to integrate
    // enforce isolation here by searching for peaks to right/left side
    bool leftPeak = i>0 && binX-buf->FindBin(peaksX[i-1]) <= iLow;
    bool rightPeak = i<peaksX.size()-1 && buf->FindBin(peaksX[i+1])-binX <= iHigh;
    if (rightPeak || leftPeak){
      //cout << "****** not isolated ******" << endl;
      continue;  // this peak is not isolated
    }
    double pulseInteg=buf->Integral(bLow,bHigh) - hbkg->Integral(bLow,bHigh);
    pIntegrals.push_back( pulseInteg );
    //cout << "integral " << x << " " << pulseInteg << endl;
    double fwhm=CalcFWHM(i);
    pFWHM.push_back( fwhm );  // for each pulse integral calculate corresponding FWHM
  }
}

double DarkPeaker::CalcFWHM(int i) const{
  if (i>=npeaks) {
    cout << "(FWHM) Warning: npoints= " << npeaks
	 << " index: " << i << " requested" << endl;
  }
  int ipeak=buf->FindBin(peaksX[i]);
  double ypeak=peaksY[i]-hbkg->GetBinContent(ipeak);
  int ilow=1,ihigh=buf->GetNbinsX();
  double xlow,xhigh;
  // walk left, find first bin < ypeak/2
  for (int i=ipeak-1; i>0 ;i--){
    if ( buf->GetBinContent(i)-hbkg->GetBinContent(i) < ypeak/2 ) {
      ilow=i;
      break;
    }
  }
  double x1=buf->GetBinCenter(ilow);
  double x2=buf->GetBinCenter(ilow+1);
  double y1=buf->GetBinContent(ilow)-hbkg->GetBinContent(ilow);
  double y2=buf->GetBinContent(ilow+1)-hbkg->GetBinContent(ilow+1);
  if (y2<y1) xlow=x1;
  else { //interpolate
    xlow = (ypeak/2-y1)*(x2-x1)/(y2-y1) + x1;
  }
  // walk right, find first bin < ypeak/2
  for (int i=ipeak+1; i<=buf->GetNbinsX(); i++){
    if ( buf->GetBinContent(i)-hbkg->GetBinContent(i) < ypeak/2 ) {
      ihigh=i;
      break;
    }
  }
  x1=buf->GetBinCenter(ihigh-1);
  x2=buf->GetBinCenter(ihigh);
  y1=buf->GetBinContent(ihigh-1)-hbkg->GetBinContent(ihigh-1);
  y2=buf->GetBinContent(ihigh)-hbkg->GetBinContent(ihigh);
  if (y2>y1) xhigh=x2;
  else {
    xhigh = (ypeak/2-y2)*(x2-x1)/(y2-y1) + x1;
  }
  //cout << "xlow/x/xhigh " << xlow <<"/"<< peaksX[i]<<"/"<<xhigh<<endl;
  return xhigh-xlow;
}

void DarkPeaker::FindBackground(){
  tcD->cd();
  int nbins=buf->GetNbinsX();
  TSPECTFLOAT *bksource=new TSPECTFLOAT[nbins];
  for (int i = 0; i < nbins; i++) bksource[i]=buf->GetBinContent(i + 1);
  /*
  // https://stackoverflow.com/questions/4364823/how-do-i-obtain-the-frequencies-of-each-value-in-an-fft
  TH1 *hfft;
  hfft=buf->FFT(0,"RE");
  double sum=0, sumw=0;
  for (int i = 1; i < nbins/2; i++){
    double mag=hfft->GetBinContent(i);
    sum += hfft->GetBinCenter(i) * mag;
    sumw+= mag;
    //    cout << i << " " << hfft->GetBinContent(i) << endl;
    // cout << i << " " << sum << " " << sumw << endl;
  }
  cout << "FFT estimate of baseline range: " << sum/sumw << endl;
  cout << "FFT estimate of baseline range: " << hfft->GetMaximumBin() << endl;
  delete hfft;
  */
  //tspectrum->Background(bksource,nbins,180
  int niter=40;  // Magic #! Width of clipping window is set to ~twice #bins spanned by signal pulse
  tspectrum->Background(bksource,nbins,niter,TSpectrum::kBackDecreasingWindow,
   			TSpectrum::kBackOrder4,kTRUE,
  			TSpectrum::kBackSmoothing3,kFALSE);
  // replace 180 w/ max fourrier component or 1/2 of it, etc
  for (int i = 0; i < nbins; i++) hbkg->SetBinContent(i + 1,bksource[i]);
  hbkg->SetLineColor(kGreen+2);
  hbkg->SetLineWidth(3);
  delete bksource;
}


// find noise from distribution pulse height spectrum
void DarkPeaker::FindNoise(){
  // background subtracted sample data
  hdist->SetBins(200,0,buf->GetMaximum());  // distribution of ADC counts
  hscan->SetBins(200,0,buf->GetMaximum()/3);  // ADC samples above threshold
  
  for (int i = 1; i <= buf->GetNbinsX(); i++){
    double val=buf->GetBinContent(i) - hbkg->GetBinContent(i);
    hdist->Fill(val);
    for (int b=1; b<=hscan->GetNbinsX(); b++){
      if (val > hscan->GetBinLowEdge(b))
	hscan->SetBinContent(b,hscan->GetBinContent(b)+1);
    }
  }		 
  
  //Using hscan to locate and eliminate noise
  //written by Grace E. Cummings, 22 July 2016
  double binwid=hdist->GetBinWidth(1);
 
  //limits and relevant values
  double maxheightthresh = hscan->GetMaximum();
  double xmaxThres = hscan->GetXaxis()->GetXmax();
  double xminThres = hscan->GetXaxis()->GetXmin();
  double xendfit= hscan->GetBinCenter(hscan->GetNbinsX());
  
  //Finds bin with contents <= XX% of maximum. This serves as end of fit
  for (int i = hscan->GetMaximumBin(); i<=hscan->GetNbinsX(); i++){
    if (hscan->GetBinContent(i)<=0.33*maxheightthresh){ // || i>=maxbinthresh+10){
      xendfit = hscan->GetBinCenter(i);
      break;
    }
  }

  //Creates Exponential fit of hscan to find lambda value. Fits the background contributions
  TF1 *thresFit = new TF1("thresFit", "expo", xminThres, xmaxThres);
  gStyle->SetOptFit(0001);
  hscan->Fit("thresFit","","",hscan->GetBinCenter(2),xendfit);
  tcD->cd();
  hscan->DrawCopy();
  //gStyle->SetOptStat(0);
  gPad->Update();
  double threshSlope = TMath::Abs(thresFit->GetParameter(1));

  //Threshold now will be 4 times 1/slope value. This should select out most background. 4 times was determined visually
  _peThreshold=4.5/threshSlope;  // BH tweak to 4.5 times
  
  //hscan option
  std::cout <<
    "1PE peak not found, estimate noise to be above "<< _peThreshold << std::endl;
}
// dark count rate in MHz
double DarkPeaker::CalcDarkRate(){
  if (buf) return npeaks / (dT*1e6*buf->GetNbinsX());
  return 0;
}



TString getoutput(TString cmd)
{
  // setup
  int MAX_BUFFER=5000;
  TString data;
  FILE *stream;
  char buffer[MAX_BUFFER];

  // do it
  stream = popen(cmd.Data(), "r");  // or c_str() for C++ string
  while ( fgets(buffer, MAX_BUFFER, stream) != NULL )
    data.Append(buffer);
  pclose(stream);

  // exit
  // return trim(data);
  return data;
}


// Create the histogram to range from log10(min) to log10(max) 

void BinLogX(TH1*h) {
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];
   
   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
     // cout << new_bins[i] << endl;
   }
   axis->Set(bins, new_bins);
   delete new_bins;
}

/*
// set up logarithmic bins
void BinLogX(TH1*h) {
  
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  from = TMath::Log10(from);
  to = TMath::Log10(to);
  
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
    cout << new_bins[i] << endl;
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}
*/
