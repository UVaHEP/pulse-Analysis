#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>

using namespace picoscope;
using std::cout;
using std::endl;


void PHD(TH1F* h, TH1F* p, int firstbin, int width){
  p->Fill(h->Integral(firstbin,firstbin+width-1));
}


void setupPicoscope(ps5000a &dev, chRange range, int samples, int nbuffers) {

  dev.open(picoscope::PS_12BIT);
  dev.setChCoupling(picoscope::A, picoscope::DC);
  dev.setChRange(picoscope::A, range);
  dev.enableChannel(picoscope::A);

  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1);
  dev.setSimpleTrigger(EXT, 18000, trgFalling, 0, 0);
  //dev.setSimpleTrigger(EXT, -10000, trgFalling, 0, 0); 
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();  
}

int main(int argc, char **argv) {

  ps5000a dev;
  chRange range = PS_50MV;
  int samples = 100;    // number of ADC samples per waveform
  int nbuffers=10000;   // number of waveforms per capture cycle
  int nrepeat=5;        // number of capture cycles
  int nsave=20;         // number of waveforms to save for samples
  int z0=0;              // starting point for background integral
  int xlow=samples/2;   // starting point for pulse integral
  int xwid=15;          // bins to integrate
  const int nBaseline=20; // number of bins to use to estimate baseline
  int opt;
  bool findWindow=false;   // find the integration window automatically
  bool quit=false;
  bool quiet=false;
  TString outfn="pulsed.root";
  
  while ((opt = getopt(argc, argv, "b:c:S:n:o:hq0xa:w:z:")) != -1) {
    switch (opt) {
    case 'b':
      nbuffers=atoi(optarg);
      break;
    case 'c':
      nrepeat=atoi(optarg);
      break;
    case 'S':
      nsave=atoi(optarg);
      break;
    case 'o':
      outfn = optarg;
      cout<<"set outfn " << outfn<<endl;
      break;
    case 'z':
      z0=atoi(optarg);
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
    case 'a':
      findWindow="true";
      break;
    case '0':   // turn off graphics
      quiet=true;  
      break;	
    case 'h':
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-b nbuffers[10000] -c nrepeat[5] -w nsave[20] -o output[darkBuffers.root]\n");
      fprintf(stderr, "-x starting bin of pulse integral[50]\n");
      fprintf(stderr, "-q exit when finished\n");	    
      exit(EXIT_FAILURE);
    }
  }
  TApplication theApp("App", &argc, argv, 0, -1);
  
  setupPicoscope(dev, range, samples, nbuffers); 


  TH2F *hpersist=new TH2F("hpersist","Persistence Display",samples,
			0,samples,250,-0.5,15000-0.5);
  hpersist->GetXaxis()->SetTitle("Sample time [2ns/div]");
  
  TH2F *hpersistCopy = (TH2F*)hpersist->Clone();
  hpersistCopy->SetTitle("Threshold Scan");
  hpersistCopy->GetYaxis()->SetTitle("Threshold [ADC]");
   
  TH1F *hsum=new TH1F("hsum","Sum of wave data",samples,0,samples);
  TH1F* hpulses1=new TH1F("hpulses1","Pulse area distribution",2000,-20000,200000);
  hpulses1->GetXaxis()->SetTitle("Pulse Area");
  TH1F* hpulses0=new TH1F("hpulses0","Pulse area distribution",2000,-20000,200000);
  hpulses0->SetLineColor(kRed);

  TH1F* waveForms[nsave];
  int wavSaved=0;
  
  for (int i = 0; i < nrepeat; i++) {
    cout << "Capturing Block:" << i << endl;     
    dev.captureBlock();

    vector <vector<short> > &data = dev.getWaveforms();
    
    for (auto &waveform : data) {
      TH1F *hsamp = new TH1F("hsamp","Samples", waveform.size(), 0, waveform.size());
      for (int i = 0; i < waveform.size(); i++) {
	hpersist->Fill(i, -1*waveform[i]);
	hpersistCopy->Fill(i, -1*waveform[i]);
	hsamp->SetBinContent(i, -1*waveform[i]);
      }
      hsum->Add(hsamp);


      if (i==0 && findWindow){  //use first buffer to figure out the pulse integration window
	double baseline=hsum->Integral(1,nBaseline)/nBaseline;
	int maxbin=hsum->GetMaximumBin();
	double max=hsum->GetBinContent(maxbin)-baseline; 
	// search find locations of ~10% peak height
	int left=hsum->FindFirstBinAbove(max/10+baseline);
	int right=hsum->FindLastBinAbove(max/10+baseline);
	// set window to 2x this "width"
	xwid=(right-left)*2;
	xlow=maxbin-(maxbin-left)*2;
	if (xlow <= nBaseline) cout << "WARNING: pulse is too close to background area, move to the right" << endl;
      }
      if (i==0)
	cout << "Setting integration window ( " << xlow << " , " << xwid << " )" << endl;
      
      //Note! Currently this is in ADC  counts, not anything else
      PHD(hsamp,hpulses1,xlow,xwid);  // integrals are inclusive over bins limits
      PHD(hsamp,hpulses0,z0,xwid);
      if (wavSaved<nsave) {
	waveForms[wavSaved]=
	  (TH1F*)hsamp->Clone(TString::Format("wave%02d",wavSaved));
	wavSaved++;
      }
      delete hsamp; 
    }
  }

  dev.close(); 
  TFile *f = new TFile(outfn, "RECREATE");
  
  TCanvas *tc=new TCanvas("tc","Pulse heights",1200,400);
  tc->Divide(3,1);
  gStyle->SetOptStat(0);
  
  // projections of persistance histogram to show counts vs threshold
  tc->cd(1)->SetLogy();
  hpersistCopy->ProjectionY("_py1",xlow,xlow+xwid)->DrawCopy();
  hpersistCopy->ProjectionY("_py0",z0,z0+xwid)->DrawCopy("same");
  hpersistCopy->GetXaxis()->SetTitle("Threshold [ADC]");

  // plot the persistance histogram
  tc->cd(2);
  hpersist->DrawCopy("col");
  hsum->Scale( 0.9 * hpersist->GetYaxis()->GetXmax() / hsum->GetMaximum() );
  hsum->SetLineColor(kGray);
  hsum->SetLineWidth(3);
  hsum->DrawCopy("same");

  // pulse height distributions
  tc->cd(3);
  hpulses1->DrawCopy(); 
  hpulses0->Scale( hpulses1->GetMaximum()/hpulses0->GetMaximum() );
  hpulses0->DrawCopy("same");
  tc->Update();

  hpersist->Write();
  hpulses0->Write();
  hpulses1->Write(); 
  hsum->Write();
  // save integration parameters
  TVectorD vInt(2);
  vInt[0]=xlow;
  vInt[1]=xwid;
  vInt.Write("vInt");
  
  f->Close();

  cout << "Hit any ^c to exit" << endl;
  theApp.Run(true);
  
  return 0;
}



