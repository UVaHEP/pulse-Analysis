#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "lightutils.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>

using namespace picoscope;


//Talks to the picoscope
void setupPicoscope(ps5000a &dev, chRange range, int samples, int nbuffers) {

  dev.open(picoscope::PS_12BIT);//allows for ADC conversion and resolution
  dev.setChCoupling(picoscope::A, picoscope::DC);
  dev.setChRange(picoscope::A, range);//allows us to set range on Picoscope
  dev.enableChannel(picoscope::A);

  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1);
  dev.setSimpleTrigger(EXT, 18000, trgFalling, 0, 0);//When triggering off anything else
  //dev.setSimpleTrigger(EXT, -10000, trgFalling, 0, 0);//When triggering off laser 
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();  
}

//Captures height of first buffer to be used to set limits on the ranges of subsequent histograms
int estimateDistLimits(vector <vector<short> > &d, bool volts, chRange range, ps5000a &dev) {
  int topOfFirstPulse;
  Int_t waveSize = d.front().size();
  TH1F* h1 = new TH1F("h1","First Pulse", waveSize, 0, waveSize);
  TH1F* h2 = new TH1F("h2","First Pulse", waveSize, 0, waveSize);
  for (auto &waveform : d.front()) {
    for (int i = 0; i < d.front().size(); i++){
      h1->SetBinContent(i, -1*d.front()[i]);
      h2->SetBinContent(i,dev.adcToMv(-1*d.front()[i],range));
    }
    if (!volts) {
      topOfFirstPulse = h1->GetMaximum();
    }
    else {
      topOfFirstPulse = h2->GetMaximum();
    }
  }

   return topOfFirstPulse;
}

int main(int argc, char **argv) {

  ps5000a dev;
  chRange range = PS_1V;   //range on picoscope, will capture amplitudes over 100pe
  int samples = 500;       // number of samples per waveform
  int nbuffers = 10000;    // number of waveforms per capture cycle
  int xlow = samples/2;    // starting point for pulse integral
  int xwid = 15;           // bins to integrate
  double theshLowLimit = 0;//just so it works in LightPeaker
  int opt;
  bool millivolts = true;  //uses mV not ADC for analysis and output
  bool quit = false;
  bool quiet = false;
  TString outfn = "lightPulseMeasurement_defaultFileName.root";
  
  while ((opt = getopt(argc, argv, "b:o:hq0x:w:z:u")) != -1) {
    switch (opt) {
    case 'b':
      nbuffers=atoi(optarg);
      break;
    case 'o':
      outfn = optarg;
      std::cout<<"set outfn " << outfn<<std::endl;
      break;
    case 'x':
      xlow=atoi(optarg);
      break;
    case 'w':
      xwid=atoi(optarg);
      break;
    case 'q':   // exit when finished
      quit=true;
      break;
    case '0':   // turn off graphics
      quiet=true;  
      break;	
    case 'u':
      millivolts=false;
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-b nbuffers[10000] -w bins to integrate -o output[lightPulseMeasurement_defaultFileName.root\n");
      fprintf(stderr, "-x starting bin of pulse integral[50]\n");
      fprintf(stderr, "-q exit when finished\n");
      fprintf(stderr, "-u use ADC instead of default mV\n");
      exit(EXIT_FAILURE);
    }
  }
  TApplication theApp("App", &argc, argv, 0, -1);

  //Take the data
  setupPicoscope(dev, range, samples, nbuffers); 
  dev.captureBlock();
  vector <vector<short> > &data = dev.getWaveforms();
  dev.close();
  
  float timebase = dev.timebaseNS();

  LightPeaker *lPk = new LightPeaker(theshLowLimit);

  Int_t waveSize = data.front().size();
  TH1F *hdata = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  TH1F *hdataMv = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);

  
  //Estimate limits, prepare peak distribution histogram
  int repPeakHeight = estimateDistLimits(data, millivolts, range, dev);
  TH1F *hPeaksDist = new TH1F("hPeakDist","Distribution of Peak heights;mV",repPeakHeight,repPeakHeight/2,repPeakHeight+repPeakHeight/2);
  TH1F *hSigmaOMean = new TH1F("hSigmaOMean","Sigma/Mean Value",1,-1,1);
  TH1F *hMeanPeakHeight = new TH1F("hMeanPeakHeight","Mean Value;; mV",1,-1,1);
  TH1F *hSigmaPeakHeight = new TH1F("hSigmaPeakHeight","Sigma for Peak Distribution",1,-1,1);
  int buffNum = 0;

  //Fills histogram with data: makes new histogram for each waveform
  for (auto &waveform : data) {
    buffNum++;
    std::cout << "Processing Buffer: " << buffNum << std::endl;
    for (int i = 0; i < waveform.size(); i++) {
      hdata->SetBinContent(i, -1*waveform[i]);
      hdataMv->SetBinContent(i,dev.adcToMv(-1*waveform[i],range));
    }

    //Prepares analyzer and picks units
    if (!millivolts) {
      lPk->SetBuffer(hdata, timebase);
    }
    else {
      lPk->SetBuffer(hdataMv, timebase);
    }

    //Using TSpectrum to find the peak
    lPk->AnalyzePeaks();
    
    //This object type is due to TSpectrum expecting more peaks/buff
    //We only have 1peak/buff, and the TSpectrum is remade each time
    //Height is always first value, as there is always only one value
    Float_t *peakHeight = lPk->GetBkgdCorrectedY();
    hPeaksDist->Fill(peakHeight[0]);
    
    //Resets histogram to increase speed, fix memory leak
    hdata->Reset();
    hdataMv->Reset();

  }

  TF1 *hPeaksDistFit = new TF1("hPeaksDistFit","gaus",repPeakHeight/2,repPeakHeight+repPeakHeight/2);
  hPeaksDist->Fit("hPeaksDistFit","","",repPeakHeight/2,repPeakHeight+repPeakHeight/2);
  
  double meanPeakHeight = hPeaksDistFit->GetParameter(1);
  double sigmaPeakHeight = hPeaksDistFit->GetParameter(2);
  double sigmaOmean = sigmaPeakHeight/meanPeakHeight;
  hSigmaOMean->SetBinContent(1,sigmaOmean);
  hMeanPeakHeight->SetBinContent(1,meanPeakHeight);
  hSigmaPeakHeight->SetBinContent(1,sigmaPeakHeight);
  
  TFile *f = new TFile(outfn, "RECREATE");  

  TCanvas *tc=new TCanvas("tc","Information about Pulses",1200,400);
  tc->Divide(2,1);
  gStyle->SetOptStat(0);
  tc->cd(1);
  hPeaksDist->DrawCopy();
  tc->cd(2);
  hSigmaOMean->DrawCopy();
  
  hPeaksDist->Write();
  hSigmaOMean->Write();
  hMeanPeakHeight->Write();
  hSigmaPeakHeight->Write();
  
  f->Close();

  std::cout<< "Close TCanvas: Information about Pulses to exit" << std::endl;
  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  theApp.Run(true);
  
  return 0;
}



