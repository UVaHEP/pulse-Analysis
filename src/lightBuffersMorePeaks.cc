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


int main(int argc, char **argv) {

  ps5000a dev;
  chRange range = PS_1V;//range on picoscope, will capture amplitudes over 100pe
  int samples = 500;   // number of samples per waveform
  int nbuffers=10000;   // number of waveforms per capture cycle
  int nrepeat=5;        // number of capture cycles
  int xlow=samples/2;   // starting point for pulse integral
  int xwid=15;          // bins to integrate
  double theshLowLimit = 0;//just so it works in LightPeaker
  int opt;
  bool millivolts = true;//uses mV not ADC for analysis and output
  bool quit=false;
  bool quiet=false;
  TString outfn="lightPulseMeasurement_defaultFileName.root";
  
  while ((opt = getopt(argc, argv, "b:c:S:n:o:hq0x:w:z:u")) != -1) {
    switch (opt) {
    case 'b':
      nbuffers=atoi(optarg);
      break;
    case 'c':
      nrepeat=atoi(optarg);
      break;
      //    case 'S':
      //nsave=atoi(optarg);
      //break;
    case 'o':
      outfn = optarg;
      std::cout<<"set outfn " << outfn<<std::endl;
      break;
      //case 'z':
      //z0=atoi(optarg);
      //break;
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
      fprintf(stderr, "-b nbuffers[10000] -c nrepeat[5] -w nsave[20] -o output[darkBuffers.root]\n");
      fprintf(stderr, "-x starting bin of pulse integral[50]\n");
      fprintf(stderr, "-q exit when finished\n");
      fprintf(stderr, "-u use ADC instead of default mV\n");
      exit(EXIT_FAILURE);
    }
  }
  TApplication theApp("App", &argc, argv, 0, -1);
  
  setupPicoscope(dev, range, samples, nbuffers); 

  
  //for (int i = 0; i < nrepeat; i++) {
  // std::cout << "Capturing Block:" << i << std::endl;     
  dev.captureBlock();
     
  // }

  vector <vector<short> > &data = dev.getWaveforms();
  float timebase = dev.timebaseNS();
  
  TH1F *hdata = NULL;
  TH1F *hdataMv = NULL;

  LightPeaker *lPk = new LightPeaker(theshLowLimit);

  Int_t waveSize = data.front().size();
  hdata = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  hdataMv = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);

    //Fills histogram with data: makes new histogram for each waveform
  for (auto &waveform : data) {
    
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

    lPk->Analyze();
    
    hdata->Reset();
    hdataMv->Reset();

  }


  TFile *f = new TFile(outfn, "RECREATE");  
  dev.close();

  
  TCanvas *tc=new TCanvas("tc","Information about Pulses",1200,400);
  tc->Divide(3,1);
  gStyle->SetOptStat(0);
  tc->cd(1);
  hdataMv->DrawCopy();
  tc->cd(2);
  hdata->DrawCopy();
  
     
  f->Close();

  std::cout<< "Close TCanvas: Information about Pulses to exit" << std::endl;
  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  theApp.Run(true);
  
  return 0;
}



