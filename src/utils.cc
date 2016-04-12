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
    hdist->Fill( buf->GetBinContent(i) - hbkg->GetBinContent(i)  );
    for (int b=1; b<=hscan->GetNbinsX(); b++){
      if (val > hscan->GetBinLowEdge(b))
	hscan->SetBinContent(b,hscan->GetBinContent(b)+1);
    }
  }		 

  double xmin=hdist->GetXaxis()->GetXmin();
  double xmax=hdist->GetXaxis()->GetXmax();

  TF1 *tf01 = new TF1("f01","gaus",xmin,xmax);

  // fit the largest peak to estimate the noise width 
  int noiseXbin=hdist->GetMaximumBin();
  double noiseX=hdist->GetBinCenter(noiseXbin);
  double noiseY=hdist->GetBinContent(noiseXbin);

  // starting parameter for noise width
  double binwid=hdist->GetBinWidth(noiseXbin);
  double fwhm=binwid;
  for (int i=noiseXbin-1; i>=0; i--){
    if (hdist->GetBinContent(i)<=noiseY/2) {
      fwhm=noiseX-hdist->GetBinCenter(i);
      break;
    }
  }
  tf01->SetParameters(noiseY,noiseX,fwhm);

  // fit the noise
  cout << "Fit for noise"<< endl;
  hdist->Fit("f01","0","",noiseX-fwhm,noiseX+fwhm);
  double par[3];  // norm, mean, sigma for noise
  tf01->GetParameters(par);
  //for (int b=1; b<=hscan->GetNbinsX(); b++){
  //  if (hscan->GetBinCenter(b)<par[1]+2*par[2]) hscan->SetBinContent(b,0);
  //}

  // use TSprectrum to find the 0,1 PE peaks
  TSpectrum *ts2=new TSpectrum();
  int npeaks=ts2->Search(hdist,par[2]/binwid,"nodraw");
  if (npeaks<2) { // can't find 1PE peak
    // place 1PE 6S.D. over noise, so we count 3S.D. fluctuations
    //snglPeak=par[1]+6*par[2];
    snglPeak=6*par[2];
    std::cout <<
      "1PE peak not found, estimate at noise+6 S.D. " <<
      snglPeak << std::endl;
  }
  else { // use 1PE from TSpectrum
    Int_t *index=new Int_t[npeaks];
    TSPECTFLOAT *xpeaks= ts2->GetPositionX();
    TSPECTFLOAT *ypeaks= ts2->GetPositionY();
    TMath::Sort(npeaks, xpeaks, index, kFALSE);  // index sort by x
    // not sure if absolute height is best or [height - mean noise]
    // use height for now
    snglPeak=xpeaks[index[1]]; 
    //snglPeak=xpeaks[index[1]]-xpeaks[index[0]];
    std::cout <<
      "1PE estimated at " <<
      snglPeak << std::endl;
  }

  //hdist->Write();
  //hscan->Write();

  delete ts2;
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
