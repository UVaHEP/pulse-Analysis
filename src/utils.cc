#include "utils.h"
#include "TMath.h"
#include <stdlib.h>
#include <pthread.h>
#include <iostream>

using std::cout;
using std::endl;



DarkPeaker::DarkPeaker(double peThreshold){
  tspectrum=new TSpectrum(5000,2);
  buf=0;
  hbkg=0;
  bkgCorrectedY=0;
  hdist=0;
  if (peThreshold) _peThreshold=peThreshold;
  tcD = new TCanvas("tcD","Dark Peak Analysis"); 
}

void DarkPeaker::SetBuffer(TH1F *newbuf, double sampleTime){
  buf=newbuf;
  dT=sampleTime;
  if (hbkg) delete hbkg;
  if (bkgCorrectedY) {
    delete bkgCorrectedY;
    bkgCorrectedY=0;
  }
  if (hdist) {
    delete hdist;
    hdist=0;
  }
  npeaks=0;
}
int DarkPeaker::GetNPeaks(){return npeaks;}
TH1F* DarkPeaker::GetBackground(){return hbkg;}
TSPECTFLOAT* DarkPeaker::GetPositionX(){ return tspectrum->GetPositionX();}
TSPECTFLOAT* DarkPeaker::GetPositionY(){ return tspectrum->GetPositionY();}
TSPECTFLOAT* DarkPeaker::GetBkgdCorrectedY(){
  if (bkgCorrectedY) return bkgCorrectedY;
  bkgCorrectedY = new TSPECTFLOAT[npeaks];
  TSPECTFLOAT *xpeaks = tspectrum->GetPositionX();
  TSPECTFLOAT *ypeaks = tspectrum->GetPositionY();
  for (int i=0; i<npeaks; i++)
    bkgCorrectedY[i]=ypeaks[i]-hbkg->GetBinContent((int)xpeaks[i]+1);
  return bkgCorrectedY;
}
TSPECTFLOAT* DarkPeaker::GetDeltaX(){
  Int_t *index=new Int_t[npeaks];
  TSPECTFLOAT *xpeaks = tspectrum->GetPositionX();
  TMath::Sort(npeaks, xpeaks, index, kFALSE);  // index sort by timestamp
  deltaX=new TSPECTFLOAT[npeaks-1];
  for (int i=0;i<npeaks-1;i++) {
    int idx1=index[i];
    int idx2=index[i+1];
    deltaX[i]=xpeaks[idx2]-xpeaks[idx1];
  }
  delete index;
  return deltaX;
}


int DarkPeaker::AnalyzePeaks(){
  if (!buf) {
    cout << "DarkPeaker::AnalyzePeaks : No buffer set!" << endl;
    return 1;
  }
  // estimate background
  hbkg=(TH1F*)buf->Clone("background");
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
  cout << "Dark pulse rate: " << CalcDarkRate() << " MHz" << endl;
  return 0;
}

void DarkPeaker::FindBackground(){
  tcD->cd();
  int nbins=buf->GetNbinsX();
  TSPECTFLOAT *bksource=new TSPECTFLOAT[nbins];
  for (int i = 0; i < nbins; i++) bksource[i]=buf->GetBinContent(i + 1);
  /*
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
  tspectrum->Background(bksource,nbins,180,TSpectrum::kBackDecreasingWindow,
   			TSpectrum::kBackOrder4,kTRUE,
  			TSpectrum::kBackSmoothing3,kFALSE);
  // replace 180 w/ max fourrier component or 1/2 of it, etc
  for (int i = 0; i < nbins; i++) hbkg->SetBinContent(i + 1,bksource[i]);
  delete bksource;
}


// find noise from distribution pulse height spectrum
void DarkPeaker::FindNoise(){
  // background subtracted sample data
  hdist=new TH1F("hdist","ADC distribution",
		       100,0,buf->GetMaximum());
  hscan=new TH1F ("hscan","Threshold scan",50,0,buf->GetBinContent(buf->GetMaximumBin()));
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
  int maxbinthresh = hscan->GetMaximumBin();
  double maxheightthresh = hscan->GetBinContent(maxbinthresh);
  double xmaxThres = hscan->GetXaxis()->GetXmax();
  double xminThres = hscan->GetXaxis()->GetXmin();
  double xendfit=hscan->GetBinCenter(hscan->GetNbinsX());
  
  //Finds bin with contents <= 10% of maximum. This serves as end of fit
  for (int i = maxbinthresh; i<=hscan->GetNbinsX(); i++){
    if (hscan->GetBinContent(i)<=0.1*maxheightthresh || i>=maxbinthresh+10){
      xendfit = hscan->GetBinCenter(i);
      break;
    }
  }

  //Creates Exponential fit of hscan to find lambda value. Fits the background contributions
  TF1 *thresFit = new TF1("thresFit", "expo", xminThres, xmaxThres);
  tcD->cd();
  hscan->Fit("thresFit","","",0,xendfit);
  hscan->DrawCopy();
  tcD->Update();
  double threshSlope = TMath::Abs(thresFit->GetParameter(1));

  //Threshold now will be 4 times 1/slope value. This should select out most background. 4 times was determined visually
  _peThreshold=4/threshSlope;
  

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
