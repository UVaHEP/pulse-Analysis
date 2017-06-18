// Pulse Average

// By TA

// This will capture a user specified number of waveforms from a picoscope and then generate and average of those waveforms using ROOT


#include "picoscopeInterface.h"
#include "ps6000.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>


using namespace picoscope;

void setupPicoscope(ps6000 &dev, chRange range, int samples, int nwaveforms, bool invert = false) {

  dev.open();
  dev.setChCoupling(picoscope::A, picoscope::DC_50R);
  dev.setChRange(picoscope::A, range);//allows us to set range on Picoscope
  dev.enableChannel(picoscope::A);
  dev.setTimebase(1); //400 ns
  if (invert) 
    dev.setSimpleTrigger(picoscope::A, -8000, trgFalling, 0, 0);
  else
    dev.setSimpleTrigger(picoscope::A, 8000, trgRising, 0, 0); 
  
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  dev.setCaptureCount(nwaveforms);
  dev.prepareBuffers();  

  
}

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



int main(int argc, char **argv) {

  TString outfn        = "lightPulseMeasurement_defaultFileName.root";
  int samples          = 3000;       // number of samples per waveform
  int nbuffers         = 10000;     // number of waveforms per capture cycle
  int opt;
  bool invert = false;
  while ((opt = getopt(argc, argv, "s:b:o:i")) != -1) {
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
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-s samples to capture, -b nbuffers[10000],  -o output[lightPulseMeasurement_defaultFileName.root");
      exit(EXIT_FAILURE);
    }
  }

  
  TApplication theApp("App", &argc, argv, 0, -1);
  ps6000 dev;
  chRange range        = PS_1V;     //range on picoscope, will capture amplitudes over
  setupPicoscope(dev, range, samples, nbuffers,invert);
  range = autoRange(dev, picoscope::A, range, samples);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();
  dev.captureBlock();
  vector <vector<short> > &data = dev.getWaveforms();
  Int_t waveSize  = data.front().size();
  dev.close();
  float timebase = dev.timebaseNS();
  std::cout << "Captured:" << data.size() << " waveforms." << std::endl; 
  TFile *f = new TFile(outfn, "RECREATE");  
  
  TH1F *hdata     = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  TH1F *hdataMv   = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);
  TProfile *hprof = new TProfile("hprof","Profile of all pulses, mV measurements", waveSize, 0, waveSize,"S");
  
  for (auto &w : data) {
    for (int i = 0; i < w.size(); i++) {
      hdata->SetBinContent(i, -1*w[i]);
      hdataMv->SetBinContent(i,dev.adcToMv(-1*w[i],range));
    }
    for (int i=1; i<=hdataMv->GetNbinsX();i++){
      hprof->Fill(hdataMv->GetBinCenter(i),
		  hdataMv->GetBinContent(i));
    }

  }
  TCanvas *tc=new TCanvas("Pulse Profile","Pulse Profile",1600,400);
  tc->cd();
  
  hprof->DrawCopy();

  hprof->Write();
  f->Close(); 

  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  
  theApp.Run(true);



  
  return 0;


}
