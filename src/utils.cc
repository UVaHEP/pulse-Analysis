#include "utils.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TList.h"
#include "TBox.h"
#include <stdlib.h>
#include <pthread.h>
#include <iostream>

using std::cout;
using std::endl;



//DarkPeaker::DarkPeaker(double peThreshold) : TSpectrum(5000) {
DarkPeaker::DarkPeaker() {
  tspectrum=new TSpectrum(5000,2);
  haveAnalysis=false;
  Reset();
  //tcD = new TCanvas("tcD","Dark Peak Analysis");
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
  _peThreshold=-1;
  _noiseCut=0;
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


int DarkPeaker::AnalyzePeaks(double peThreshold){
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
  if (peThreshold>0) _peThreshold=peThreshold/ buf->GetMaximum();
  else threshold = _noiseCut / buf->GetMaximum();

  double sigma= 30;// Add adjustment to options //2; // this can/should be optimzed

  if (threshold>1) {
    npeaks=0;
    return 0;
  }
  npeaks=tspectrum->Search(buf,sigma,"nobackground,nomarkov,nodraw",threshold);
  //npeaks=tspectrum->Search(buf,sigma,"nomarkov,nodraw",threshold);
  //int npeaks=tspectrum->Search(buf,sigma,"nomarkov",threshold);
  //cout << "Found " << npeaks << " peaks" << endl;
  //cout << "Approximate DCR: " << CalcDarkRate() << " MHz" << endl;

  // retrieve the peaks and sort by time order
  Int_t *index=new Int_t[npeaks];
  TSPECTFLOAT *xpeaks = tspectrum->GetPositionX();
  TSPECTFLOAT *ypeaks = tspectrum->GetPositionY();
  TMath::Sort(npeaks, xpeaks, index, kFALSE);  // index sort by timestamp
  for (int i=0;i<npeaks;i++) {
    peaksX.push_back( xpeaks[index[i]] );
    peaksY.push_back( ypeaks[index[i]] );
  }
  delete index, xpeaks, ypeaks;
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
void DarkPeaker::Integrate(int nLow, int nHigh, bool selectIsolated){
  for (int i=0; i<peaksX.size(); i++){
  
    // find bins for integral
    double x=peaksX[i];
    int binX=buf->FindBin(x);
    int bLow=binX-nLow;
    int bHigh=binX+nHigh;

    if ( (binX<nLow+1)||(binX+nHigh > buf->GetNbinsX()) ) continue;  // pulse is too close to edge of histogram

    double fwhm=CalcFWHM(i);
    pFWHM.push_back( fwhm );  // for each pulse integral calculate corresponding FWHM
    
    // enforce ~5sigma isolation here by searching for peaks to right/left side 
    double isoX = fwhm*5.0/2;
       
    bool hasEarlyPeak = i>0 && binX-buf->FindBin(peaksX[i-1]) <= nLow;
    bool hasLatePeak = i<peaksX.size()-1 && buf->FindBin(peaksX[i+1])-binX <= nHigh;
    if (hasEarlyPeak || hasLatePeak) continue;  // this peak is not isolated

    double pulseInteg=buf->Integral(bLow,bHigh) - hbkg->Integral(bLow,bHigh);
    pIntegrals.push_back( pulseInteg );

  }
}

// FWHM in same units as x-axis
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
  return xhigh-xlow;
}

void DarkPeaker::FindBackground(){
  //tcD->cd();
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


// Find noise from distribution pulse height spectrum
// Set pe search threshold based on noise distribution, this can be changed by the user when
// calling the AnalyzePeaks method
void DarkPeaker::FindNoise(){
  // background subtracted sample data
  hdist->SetBins(256,0,1<<14);  // distribution of ADC counts
  hscan->SetBins(256,0,1<<14);  // ADC samples above threshold
  
  for (int i = 1; i <= buf->GetNbinsX(); i++){
    double val=buf->GetBinContent(i) - hbkg->GetBinContent(i);
    hdist->Fill(val);
    for (int b=1; b<=hscan->GetNbinsX(); b++){
      if (val > hscan->GetBinLowEdge(b))
	hscan->SetBinContent(b,hscan->GetBinContent(b)+1);
    }
  }		 
  

  //  double binwid=hdist->GetBinWidth(1);
 
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

  // Exponential fit for hscan to find lambda value.
  TF1 *thresFit = new TF1("thresFit", "expo", xminThres, xmaxThres);
  gStyle->SetOptFit(0001);
  hscan->Fit("thresFit","Q0","",hscan->GetBinCenter(2),xendfit);
  //tcD->cd();
  //hscan->DrawCopy();
  //gStyle->SetOptStat(0);
  gPad->Update();
  double threshSlope = TMath::Abs(thresFit->GetParameter(1));

  // This should select out most background. Ad hoc determination of magic number=4.5
  _noiseCut=4.5/threshSlope;
  
  //hscan option
  //std::cout <<
  //  "1PE peak not found, estimate noise to be above "<< _peThreshold << std::endl;
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


// fit a double exponential to separate DCR and afterpulsing in Delta_time distribution
void fit2Exp(TH1F *h, TString opt){
  // 1st fit for DCR using long delta time tail
  int startBin = 8;  // HACK!!! Start this fit after AP peak
  h->Fit("expo","LQ","",h->GetBinCenter(startBin),h->GetXaxis()->GetXmax());
  TF1 *fun1=new TF1(*(h->GetFunction("expo")));
  fun1->SetName("expoDCR");
  double a1,b1;  //fit parameters for DCR
  a1=fun1->GetParameter(0);
  b1=fun1->GetParameter(1); 
  
  // fit parameters for afterpulse component
  // we start the combined fit at the maximum bin, a minimum bin value is set to
  // prevent crazy fits at low statistics (~flat timing distribution)
  double xAtMax=h->GetBinCenter( TMath::Min(3,h->GetMaximumBin()) );
  double startx=h->GetBinCenter(xAtMax);
  double apNorm=h->GetMaximum()-fun1->Eval(xAtMax);
  double a2,b2; // fit parameters for APulse
  a2=TMath::Log(apNorm/10);  // Assume ~10% AP at lowest dT to start
  b2=20*b1;  // Assume much faster fall off of AP spectrum
  
  TF1 *afterPulseFit = new TF1("afterPulseFit","exp([0]+[1]*(x-[4]))+exp([2]+[3]*(x-[4]))",xAtMax, h->GetXaxis()->GetXmax());
  afterPulseFit->SetParNames("A_{1}","Lamba_{DCR}","A_{2}","Lamba_{2}");
  afterPulseFit->SetParameters(a1,b1,a2,b2,xAtMax);
  afterPulseFit->FixParameter(4,xAtMax);
  h->Fit("afterPulseFit",opt,"",xAtMax,h->GetXaxis()->GetXmax());

  // update DCR fit function parameters and range for drawing
  fun1->SetParameters( afterPulseFit->GetParameter(0), afterPulseFit->GetParameter(1) );
  fun1->SetParError(0, afterPulseFit->GetParError(0));
  fun1->SetParError(1, afterPulseFit->GetParError(1));
  fun1->SetRange(xAtMax,h->GetXaxis()->GetXmax());
  fun1->SetLineStyle(2);
  fun1->DrawCopy("same");
  
  h->GetListOfFunctions()->Add(fun1);
  h->GetListOfFunctions()->Add(afterPulseFit);
}


Double_t fcnDcrAp(Double_t *x, Double_t *par){
  double norm=par[0];
  double pa=par[1];  // AP fraction
  double ta=par[2];  // AP time constant
  double td=par[3];  // Dark count time constant
  return norm * ( pa/ta*exp(-x[0]/ta) + (1.0-pa)/td*exp(-x[0]/td) );
}


Double_t fcnDcr(Double_t *x, Double_t *par){
  double norm=par[0];
  double pa=par[1];  // AP fraction
  double ta=par[2];  // AP time constant
  double td=par[3];  // Dark count time constant
  return norm * ( (1.0-pa)/td*exp(-x[0]/td) );
}

FitDcrAp::FitDcrAp(){
  afterPulseFit=new TF1("fcnDcrAp",fcnDcrAp,0,1,4);
  afterPulseFit->SetLineColor(kRed);
  afterPulseFit->SetParNames("norm","fAP","tauAP","tauDCR");
  dcrFcn=new TF1("fcnDcr",fcnDcr,0,1,4);
  dcrFcn->SetLineColor(kGreen);
}



void FitDcrAp::Fit(TH1F *h, TString opt){
  // start with exponential fit to whole distribution
  double startx=h->GetBinCenter( h->GetMaximumBin() );
  double endx=h->GetXaxis()->GetXmax();
  h->Fit("expo","LQ","",startx,endx);
  expoFit=(TF1*)(h->GetFunction("expo")->Clone());
  expoFit->SetLineStyle(2);
  double a1=expoFit->GetParameter(0);
  double b1=expoFit->GetParameter(1);
  double norm=h->Integral(h->GetMaximumBin(), h->GetNbinsX());
  norm*=h->GetBinWidth(1);
  
  afterPulseFit->SetRange(startx,endx);
  afterPulseFit->SetParNames("norm","fAP","tauAP","tauDCR");
  dcrFcn->SetRange(startx,endx);
  afterPulseFit->SetParameter(0,norm);
  afterPulseFit->SetParameter(1,0.05);  // starting guess at 5% for afterpulse
  afterPulseFit->SetParLimits(1,0.0,0.5);
  afterPulseFit->SetParameter(2,-0.1/b1);
  afterPulseFit->SetParLimits(2,0,-0.9/b1);
  afterPulseFit->SetParameter(3,-1/b1);
  h->Fit("fcnDcrAp",opt,"",startx,endx);
  afterPulseFit->DrawCopy("same");
  expoFit->DrawCopy("same");
}


void calcDcrAp(TH1F *hdTime, double &fitDCR, double &errDCR, double &rateAP, double errAP){
  fit2Exp(hdTime,"L");
  TF1 *tf_expoDCR=hdTime->GetFunction("expoDCR");
  double par[2];
  tf_expoDCR->GetParameters(par);
  double meanDt = -2/par[1];  // Delta_t in [ns]
  fitDCR = 1 / meanDt * 1000;  // in MHz
  errDCR = fitDCR*tf_expoDCR->GetParError(1)/par[1];
  rateAP=errAP=0;
  if (fitDCR<0 || TMath::Abs(errDCR/fitDCR)>0.25){
    std::cout << "Warning Fit yields negative DCR or large error" << std::endl;
    std::cout << "Afterpulse calculation skipped" << std::endl;
  }
  else {  // calculate afterpulse rate
    //Written in part by Grace E. Cummings, 30 July 2016
    TF1 *tf_apFcn = hdTime->GetFunction("afterPulseFit");
    double xmaxBinTime = hdTime->GetBinCenter( hdTime->GetMaximumBin() );
    double xmaxTime = hdTime->GetXaxis()->GetXmax();

    //ap Rate calculated as ratio of fit of afterpulses to fit of dark counts
    //Fit of afterpulses is the afterPulseFit-expoDCR
    double aPFint = tf_apFcn->Integral(xmaxBinTime,xmaxTime);
    double expDCRint = tf_expoDCR->Integral(xmaxBinTime,xmaxTime);
    double binWidthhdTime = hdTime->GetBinWidth(1);
    
    rateAP = (aPFint-expDCRint)/expDCRint;
    
    //Error in afterpulse rate. Sqrt(excess)~sigma of excess
    double excess = (aPFint-expDCRint)/binWidthhdTime;
    errAP=sqrt(excess)/(expDCRint/binWidthhdTime);    
  }
}


TF1* FitDcrAp::GetDcrFcn(){
  double par[4];
  afterPulseFit->GetParameters(par);
  dcrFcn->SetParameters(par);
  return dcrFcn;
}

// estimate afterpulses from 2D plot
int apCalc2D(TH2F *h, double dcr, double pe, double sigma, bool show){
  int ymaxAP=h->GetYaxis()->FindBin(pe-2.5*sigma);
  int ymax=h->GetNbinsY();
  int xmaxAP=h->GetXaxis()->FindBin(1/dcr);  // DCR in Hz
  int xmax=h->GetNbinsX();
  double countAP=0;
  for (int bx=1; bx<xmax;bx++){
    for (int by=1; by<ymax;by++){
      int gbin=h->GetBin(bx,by);
      double val=h->GetBinContent(gbin);
      if (bx<xmaxAP && by<ymaxAP) countAP+=val;
    }
  }
  if (show){
    TCanvas *tc=new TCanvas("tctemp","2D AP measurement");
    tc->cd(1)->SetLogx();
    h->DrawCopy("col");
    TBox *box=new TBox(h->GetXaxis()->GetBinLowEdge(1),
		       h->GetYaxis()->GetBinLowEdge(1),
		       h->GetXaxis()->GetBinLowEdge(xmaxAP),
		       h->GetYaxis()->GetBinLowEdge(ymaxAP));
    box->SetLineColor(2);
    box->SetFillStyle(0);
    box->Draw("l");
  }
  return countAP;
}



