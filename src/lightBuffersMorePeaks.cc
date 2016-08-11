#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TProfile.h"
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

//Uses a rough estimate of mean as a basis for histogram and fit limits
int baseDistLimits(Float_t *heights, int nbuffers) {
  int maximum = 0;
  int minimum = maximum;
  int roughMiddle = 0;
  for (int i =0; i < nbuffers; i++) {
    if (heights[i] > maximum){
      maximum = heights[i];
	}
  }
  return maximum;
}

//Trying to make a function to make it easier to interact with the heights

LightPeaker *lPk = new LightPeaker(0);

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

  //Get's peak heights, and make array
  Float_t *peakHeight;
  Int_t waveSize = data.front().size();
  int buffNum = 0;
  TSPECTFLOAT *heightsOfPeaks = new TSPECTFLOAT[nbuffers];
  TH1F *hdata = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  TH1F *hdataMv = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);
  TH1F *hsum = new TH1F("hsum","Sum of all pulses, mV measurements", waveSize, 0, waveSize);
  TProfile *hprof = new TProfile("hprof","Profile of all pulses, mV measurements", waveSize, 0, waveSize,"S");

  for (auto &waveform : data) {
    buffNum++;
    std::cout << "Processing Buffer: " << buffNum << std::endl;
    for (int i = 0; i < waveform.size(); i++) {
      hdata->SetBinContent(i, -1*waveform[i]);
      hdataMv->SetBinContent(i,dev.adcToMv(-1*waveform[i],range));
    }

    //Sum of all data
    hsum->Add(hdataMv);
    for (int i=1; i<=hdataMv->GetNbinsX();i++){
      hprof->Fill(hdataMv->GetBinCenter(i),
		 hdataMv->GetBinContent(i));
    }
    
    //Prepares analyzer and picks units
    if (!millivolts) {
      lPk->SetBuffer(hdata, timebase);
    }
    else {
      lPk->SetBuffer(hdataMv, timebase);
    }

    lPk->AnalyzePeaks();
    peakHeight = lPk->GetBkgdCorrectedY();

    heightsOfPeaks[buffNum-1] = peakHeight[0];

    hdata->Reset();
    hdataMv->Reset();
  }

    
  int repPeakHeight = baseDistLimits(heightsOfPeaks,nbuffers);

  //Estimate limits, prepare peak distribution histogram
  TH1F *hPeaksDist = new TH1F("hPeakDist","Distribution of Peak heights;mV",repPeakHeight,repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10) ;
  TH1F *hSigmaOMean = new TH1F("hSigmaOMean","Sigma/Mean Value",1,-1,1);
  TH1F *hMeanPeakHeight = new TH1F("hMeanPeakHeight","Mean Value;; mV",1,-1,1);
  TH1F *hSigmaPeakHeight = new TH1F("hSigmaPeakHeight","Sigma for Peak Distribution",1,-1,1);

  for (int i = 0; i < nbuffers; i++) {
    hPeaksDist->Fill(heightsOfPeaks[i]);
      }

  
  TF1 *hPeaksDistFit = new TF1("hPeaksDistFit","gaus",repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10);
  hPeaksDist->Fit("hPeaksDistFit","","",repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10);
  
  double meanPeakHeight = hPeaksDistFit->GetParameter(1);
  double sigmaPeakHeight = hPeaksDistFit->GetParameter(2);
  double sigmaOmean = sigmaPeakHeight/meanPeakHeight;
  hSigmaOMean->SetBinContent(1,sigmaOmean);
  hMeanPeakHeight->SetBinContent(1,meanPeakHeight);
  hSigmaPeakHeight->SetBinContent(1,sigmaPeakHeight);
  
  TFile *f = new TFile(outfn, "RECREATE");  

  TCanvas *tc=new TCanvas("tc","Information about Pulses",1200,800);
  tc->Divide(2,2);
  gStyle->SetOptStat(0);
  tc->cd(1);
  hPeaksDist->DrawCopy();
  tc->cd(2);
  hSigmaOMean->DrawCopy();
  tc->cd(3);
  hsum->DrawCopy();
  tc->cd(4);
  hprof->DrawCopy();
  
  hPeaksDist->Write();
  hSigmaOMean->Write();
  hMeanPeakHeight->Write();
  hSigmaPeakHeight->Write();
  hsum->Write();
  hprof->Write();
  
  
  f->Close();

  std::cout<< "Close TCanvas: Information about Pulses to exit" << std::endl;
  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  theApp.Run(true);
  
  return 0;
}



