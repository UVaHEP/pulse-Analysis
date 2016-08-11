#include "lightutils.h"
#include "TMath.h"
#include <stdlib.h>
#include <pthread.h>
#include <iostream>

using std::cout;
using std::endl;


LightPeaker::LightPeaker(double peThreshold){
  tspectrum=new TSpectrum(5000,2);
  buf=0;
  hbkg=0;
  bkgCorrectedY=0;
  hdist=0;
  if (peThreshold) _peThreshold=peThreshold;
  //  tcD = new TCanvas("tcD","Light Peak Analysis"); 
}

void LightPeaker::SetBuffer(TH1F *newbuf, double sampleTime){
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
int LightPeaker::GetNPeaks(){return npeaks;}
TH1F* LightPeaker::GetBackground(){return hbkg;}
TSPECTFLOAT* LightPeaker::GetPositionX(){ return tspectrum->GetPositionX();}
TSPECTFLOAT* LightPeaker::GetPositionY(){ return tspectrum->GetPositionY();}

TSPECTFLOAT* LightPeaker::GetBkgdCorrectedY(){
  if (bkgCorrectedY) return bkgCorrectedY;
  bkgCorrectedY = new TSPECTFLOAT[npeaks];
  TSPECTFLOAT *xpeaks = tspectrum->GetPositionX();
  TSPECTFLOAT *ypeaks = tspectrum->GetPositionY();
  for (int i=0; i<npeaks; i++)
    bkgCorrectedY[i]=ypeaks[i]-hbkg->GetBinContent((int)xpeaks[i]+1);
  return bkgCorrectedY;
}

int LightPeaker::AnalyzePeaks(){
  if (!buf) {
    cout << "LightPeaker::AnalyzePeaks : No buffer set!" << endl;
    return 1;
  }
  // estimate background
  hbkg=(TH1F*)buf->Clone("background");
  FindBackground();

  // estimate Threshold for the greatest pe peak
  FindThreshold();

  // find peaks
  double threshold;
  if (_peThreshold>-1e19) threshold = _peThreshold / buf->GetMaximum();
  else threshold = (snglPeak/2) / buf->GetMaximum();


  double sigma=2; // this can/should be optimzed
  npeaks=tspectrum->Search(buf,sigma,"nobackground,noMarkov,nodraw",threshold);
  //npeaks=tspectrum->Search(buf,sigma,"nomarkov,nodraw",threshold);
  //int npeaks=tspectrum->Search(buf,sigma,"nomarkov",threshold);
  //cout << "Found " << npeaks << " peaks" << endl;
  return 0;
}

void LightPeaker::FindBackground(){
  /* tcD = new TCanvas("tcD","Light Peak Analysis"); 
  std::cout << "tcD's Name:" << tcD->GetName() << std::endl;
  tcD->cd();*/
  int nbins=buf->GetNbinsX();
  TSPECTFLOAT *bksource=new TSPECTFLOAT[nbins];
  for (int i = 0; i < nbins; i++) bksource[i]=buf->GetBinContent(i + 1);
    tspectrum->Background(bksource,nbins,180,TSpectrum::kBackDecreasingWindow,
   			TSpectrum::kBackOrder4,kTRUE,
  			TSpectrum::kBackSmoothing3,kFALSE);
  // replace 180 w/ max fourier component or 1/2 of it, etc
  for (int i = 0; i < nbins; i++) hbkg->SetBinContent(i + 1,bksource[i]);
  delete bksource;
}


// find noise from distribution pulse height spectrum
void LightPeaker::FindThreshold(){
  double maxHeightBkgCorrected = buf->GetMaximum();

  //right now, these histograms are not being used.
  //staying defined because I will probably think of a better way to select the threshold
  //background subtracted sample data
  /*  hdist=new TH1F("hdist","Pulse Height Distribution",
		       100,0,buf->GetMaximum());
  hscan=new TH1F ("hscan","Threshold scan",50,0,buf->GetBinContent(buf->GetMaximumBin()));
  for (int i = 1; i <= buf->GetNbinsX(); i++){
    double val=buf->GetBinContent(i) - hbkg->GetBinContent(i);
    hdist->Fill(val);
    for (int b=1; b<=hscan->GetNbinsX(); b++){
      if (val > hscan->GetBinLowEdge(b))
	hscan->SetBinContent(b,hscan->GetBinContent(b)+1);
    }
    }*/		 
  
  _peThreshold = maxHeightBkgCorrected*.5;//only want the highest peaks
  

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
