// Pulse Average

// By TA

// This will capture a user specified number of waveforms from a picoscope and then generate and average of those waveforms using ROOT


#include "picoscopeInterface.h"
#include "ps6000.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TString.h"
#include "TTimer.h"
#include "TSystem.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
using std::vector;


using namespace picoscope;

struct FitResults {
  double baseline;
  double Qtotal;
  int ncells;
  double M;
  double tauRise;
  double tauFall;
  double A;  // pulse amplitude
};


struct FitResults pulseFitter(TProfile *h, double dt, int ncells=1);

void setupPicoscope(ps6000 &dev, chRange range, float trigger, int samples, int nwaveforms, bool invert = false) {

  dev.open();
  dev.setChCoupling(picoscope::A, picoscope::DC_50R);
  dev.setChRange(picoscope::A, range);//allows us to set range on Picoscope
  dev.enableChannel(picoscope::A);

  dev.setTimebase(1); //400 ns
  /*  
  if (invert) 
    dev.setSimpleTrigger(picoscope::A, dev.mVToADC(trigger, range), trgFalling, 0, 0);
  else
    dev.setSimpleTrigger(picoscope::A, dev.mVToADC(trigger,range), trgRising, 0, 0); 
  */
  //dev.setSimpleTrigger(picoscope::A, dev.mVToADC(trigger,range), trgRising, 0, 0); 
  dev.setSimpleTrigger(AUX, -100, trgFalling, 0, 0);//When triggering off anything else Jul 10 changed to rising from falling

  dev.setSamples(samples); 
  //dev.setPreTriggerSamples(samples/10);
  //dev.setPostTriggerSamples(samples/90);
  dev.setCaptureCount(nwaveforms);
  dev.prepareBuffers();  

  
}

/*
chRange autoRange(ps6000 &dev, chName name, chRange startRange, int samples) {

  chRange ranges[] = { PS_50MV, PS_100MV, PS_200MV, PS_500MV, PS_1V, PS_2V, PS_5V};
  int pos = 0; 
  while (pos < 7) {
    if (ranges[pos] == startRange)
      break;
    pos++; 
  }
  dev.setChRange(name, startRange); 
  dev.setCaptureCount(1);
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  dev.prepareBuffers();
  bool searching = true; 
  while(searching) {
    dev.setChRange(name, ranges[pos]);
    dev.captureBlock();
    searching = false; 
    vector<vector<short> > & data = dev.getWaveforms();
    for (auto &w : data) {
      for (int i = 0; i < w.size(); i++) {
	if (w[i] > 0.9*dev.getMaxADC() or w[i] < 0.9 * dev.getMinADC()) {
	  // If we've got values > 90% of our max adc value, time to increase
	  pos++; 
	  searching = true;
	  if (pos >= 7) { 
	    std::cout << "Signal greater than our possible ranges, returning 5V range" << std::endl;
	    searching = false;
	    return PS_5V; 
	  }
	  break;
	}
      }
    }
  }

  std::cout << "Range set to:" << ranges[pos] << std::endl; 
  return ranges[pos];
    
}
*/

chRange autoRange(ps6000 &dev, int nbuf){
  TFile autorangeF("autorange.root", "UPDATE");
  
  vector <vector<short> > data;
  int mvRange[]={10,20,50,100,200,500,1000,2000,5000};
  chRange ranges[] = {PS_10MV,PS_20MV, PS_50MV, PS_100MV, PS_200MV, PS_500MV, PS_1V, PS_2V, PS_5V};
  //dev.setChRange(picoscope::A, PS_50MV);
  //dev.setChCoupling(picoscope::A, picoscope::DC_50R);
  dev.setCaptureCount(nbuf);
  gSystem->Sleep(1000);
  chRange autoRange=PS_50MV;
  int mvScale=0;
  chRange range;
  for ( int psRange=PS_50MV; psRange <= PS_5V; psRange++ ){
    std::cout<<"Autoranging pass: " << mvRange[psRange] << "mV range" << std::endl;
    mvScale=mvRange[psRange];
    range=ranges[psRange];
    autoRange=(chRange)psRange;
    // set picoscope to test voltage range
    dev.setChRange(picoscope::A, autoRange);
    dev.prepareBuffers();
    dev.captureBlock(); 
    data = dev.getWaveforms();   // can we call ps5000aMaximumValue etc on this object?
    dev.prepareBuffers();
    dev.captureBlock(); 
    data = dev.getWaveforms();
    bool overThresh=false;
    int wcount = 0;
    TH1F *autoWave = NULL;
    for (auto &waveform : data) {
      //print first 5
      if (wcount < 2) {
	autoWave = new TH1F("hdata_autorange", "Autorange Waveforms",
			    waveform.size(), 0, waveform.size());
	for (int i = 0; i < waveform.size(); i++) {
	  autoWave->SetBinContent(i+1, waveform[i]);
	}
	autoWave->Write();
	
	std::sort(waveform.begin(), waveform.end(), std::greater<short>());
	for (int i = 0; i < 5; i++) {
	  std::cout << waveform[i] << "\t";
	}
	
	std::cout << std::endl;
	std::cout << std::endl;
	

      wcount++;
      }


      
      for (int i = 0; i < waveform.size(); i++) {

	if (abs(waveform[i])>26000){      // ~ 80% of ADC range
	  //std::cout << "Overthresh" << std::endl;
	  overThresh=true;
	  break;
	}
      }
    } // end of buffer loop
    if (!overThresh) return range;  // CLEAN ME UP!
  }
  return range;
}


int main(int argc, char **argv) {

  TString outfn        = "lightPulseMeasurement_defaultFileName.root";
  int samples          = 1000;       // number of samples per waveform
  int nbuffers         = 10000;     // number of waveforms per capture cycle
  int ncells           = 1;         // number of cells if applying fit
  float trigger        = 50;
  int opt;
  int freq             = 0;//frequency of laser
  float vexcess        = 99;
  bool invert = false;
  bool doFit = false;
  while ((opt = getopt(argc, argv, "Fis:b:o:n:t:f:r:")) != -1) {
    switch (opt) {
    case 's':
      samples = atoi(optarg);
      break;
    case 'b':
      nbuffers=atoi(optarg);
      break;
    case 'o':
      outfn = optarg;
      std::cout<<"set outfn " << outfn<<std::endl;
      break;
    case 'i':
      invert = true;
      break;
    case 'n':
      std::cout << "optarg:" << optarg << std::endl;
      ncells=atoi(optarg);
      break;
    case 't':
      trigger = atof(optarg);
      break;
    case 'F':
      doFit = true;
      break;
    case 'f':
      freq = atoi(optarg);
      break;
    case 'r':
      vexcess = atof(optarg);
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-s samples to capture, -b nbuffers[10000],  -o output[lightPulseMeasurement_defaultFileName.root\n");
      exit(EXIT_FAILURE);
    }
  }
      std::cout <<"over voltage " << vexcess << std::endl;
      std::cout <<"frequency  " << freq<< std::endl;
      std::cout <<"ncells     " << ncells <<std::endl;

  
  TApplication theApp("App", &argc, argv, 0, -1);
  ps6000 dev;
  chRange range = PS_50MV;     //range on picoscope, will capture amplitudes over

  setupPicoscope(dev, range, trigger, samples, nbuffers,invert);
  range=autoRange(dev,1000);
  //setupPicoscope(dev, range, trigger, samples, nbuffers,invert);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();
  dev.captureBlock();
  int nSaves=20;
  vector<TH1F *> wavesSave;
  vector <vector<short> > &data = dev.getWaveforms();
  Int_t waveSize  = data.front().size();
  float adcmV=dev.adcToMv(1, range);
  dev.close();
  float timebase = dev.timebaseNS();
  std::cout << "Captured:" << data.size() << " waveforms." << std::endl;
  std::cout << "Buf Size:" << waveSize << " samples." << std::endl;
  TFile *f = new TFile(outfn, "RECREATE");  
  
  TH1F *hdata     = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  TH1F *hdataMv   = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);
  TProfile *hprof = new TProfile("hprof","Profile of all pulses, mV measurements", waveSize, 0, waveSize, "S");

  for (auto &w : data) {
    for (int i = 0; i < w.size(); i++) {
      hdata->SetBinContent(i, -1*w[i]);
      hdataMv->SetBinContent(i,dev.adcToMv(-1*w[i],range));
    }
    if (wavesSave.size()<nSaves){
      wavesSave.push_back((TH1F*)hdata->Clone(TString::Format("wave%02d",wavesSave.size()+1)));
    }
    for (int i=1; i<=hdataMv->GetNbinsX();i++){
      hprof->Fill(hdataMv->GetBinCenter(i),
		  hdataMv->GetBinContent(i));
    }
    hdata->Clear();
  }
  TCanvas *tc=new TCanvas("Pulse Profile","Pulse Profile",1600,400);
  tc->cd();
  TTimer *timer = new TTimer();
  timer->Connect("Timeout()","TCanvas",tc,"Close()");
  timer->Start(2000,kTRUE);
  
  hprof->DrawCopy();
  
  struct FitResults result;
  if (doFit) {
    std::cout << "Fitting!" << std::endl;
    result=pulseFitter(hprof,timebase,ncells);
  }
  for (int i=0; i<wavesSave.size(); i++) wavesSave[i]->Write();

  
  hprof->Write();
  TH1F *dT  = new  TH1F("dT", "Time Steps [s]", 1,-0.5,0.5);
  dT->SetBinContent(1, timebase);
  dT->Write();
  TH1F *hRange  = new  TH1F("hRange", "Picoscope Range Setting", 2,0,2);
  hRange->SetBinContent(1, range);
  hRange->SetBinContent(2, adcmV);
  hRange->Write();
  TH1F *hGain  = new  TH1F("hGain", "Calculated Gain", 1,-0.5,0.5);
  hGain->SetBinContent(1, result.M);
  hGain->Write();
  TH1F *hAmpl  = new  TH1F("hAmpl", "Mean pulse amplitude", 1,-0.5,0.5);
  hAmpl->SetBinContent(1, result.A);
  hAmpl->Write();
  TH1F *hFreq  = new TH1F("hFreq", "Frequency of laser", 1,-0.5,0.5);
  hFreq->SetBinContent(1,freq);
  hFreq->Write();
  TH1F *hVexcess = new TH1F("hVexcess", "Overvolage (V)",1,-0.5,0.5);
  std::cout<<"excess voltage "<< vexcess<< std::endl;
  hVexcess->SetBinContent(1,vexcess);
  hVexcess->Write(); 
  
  f->Close(); 

  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  
  theApp.Run(true);



  
  return 0;


}

// Y-range of TProfile is expected to be in millivolts
struct FitResults pulseFitter(TProfile *hprof, double dt, int ncells){
  double rLoad=50;
  double q_e=1.602e-19;
  int nsamples=0;

  // force positive polarity
  double xmin=hprof->GetMinimum();
  double xmax=hprof->GetMaximum();
  if (TMath::Abs(xmin)>TMath::Abs(xmax)) hprof->Scale(-1);

  // use first ~10% of range as baseline estimate
  double BLSrange=0.1;
  double baseline=0;
  int nbins=hprof->GetNbinsX();
  for (int i=1; i<=(int)(nbins*BLSrange); i++){
    nsamples++;
    baseline+=hprof->GetBinContent(i);
  }
  baseline/=nsamples;
  
  // Q = I*t = (V-baseline)/rLoad * t
  // V in volts
  double Q=(hprof->Integral()-hprof->GetNbinsX()*baseline)/rLoad*dt/1000;
  Q=TMath::Abs(Q);
  // calclate gain
  double M = Q/ncells/q_e;

  // do fit to falling edge y=exp(A + B*x), B=1/tau
  double max=hprof->GetMaximum();
  int blow=hprof->FindLastBinAbove(0.8*(max-baseline)+baseline);
  int bhigh=hprof->FindLastBinAbove(0.2*(max-baseline)+baseline);
  xmin=hprof->GetBinCenter(blow);
  xmax=hprof->GetBinCenter(bhigh);

  hprof->Fit("expo","Lq","",xmin,xmax);
  TF1 *fnFall= (TF1*)(hprof->GetFunction("expo")->Clone("fnFall"));
  double tauFall=TMath::Abs(1/hprof->GetFunction("expo")->GetParameter(1))*dt;


  // fit rising edge
  // note: the rising edge may be dominated by scope rise time
  blow=hprof->FindFirstBinAbove(0.2*(max-baseline)+baseline)-1;
  bhigh=hprof->FindFirstBinAbove(0.8*(max-baseline)+baseline)-1;
  xmin=hprof->GetBinCenter(blow-1);
  xmax=hprof->GetBinCenter(bhigh+1);
  hprof->Fit("expo","Lq","",xmin,xmax);
  double tauRise=TMath::Abs(1/hprof->GetFunction("expo")->GetParameter(1))*dt;


  hprof->DrawCopy();
  fnFall->DrawCopy("lsame");
  std::cout << "baseline estimate: " << baseline << " mV" << std::endl;
  std::cout << "total charge in pulse: " << Q << std::endl;
  std::cout << "tau(rising) " << tauRise << std::endl;
  std::cout << "tau(falling) " << tauFall << std::endl;
  std::cout << "Gain: " << M << std::endl;
  std::cout << "Mean pulse height " << hprof->GetMaximum()-baseline << " mV" << std::endl;
  if (ncells==1) cout << "ncells set to 1, divide gain by ncells" << std::endl;
  // check that pulse is well centered and not truncated
  int iRise=(1+tauRise/dt)*5;
  int iFall=(1+tauFall/dt)*5;
  if ( hprof->GetBinLowEdge(hprof->GetMaximumBin())-iRise < (int)(nbins*BLSrange) ) {
    std::cout << "Baseline calculation may be biased.  Is the pulse too close to left edge of histogram? " << std::endl;
  }
  if ( hprof->GetBinLowEdge(hprof->GetMaximumBin())+iFall > nbins ){
    std::cout << "Pulse may be truncated at right side of histogram" << std::endl;
  }
  struct FitResults result;
  result.baseline=baseline;
  result.Qtotal=Q;
  result.ncells=ncells;
  result.M=M;
  result.tauRise=tauRise;
  result.tauFall=tauFall;
  result.A=hprof->GetMaximum()-baseline;
  return result;
}

