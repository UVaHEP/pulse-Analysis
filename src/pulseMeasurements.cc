#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "picoutils.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include "TTimer.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <vector>

using namespace picoscope;
using std::cout;
using std::endl;
using std::vector;

double PHD(TH1F* h, TH1F* p, int firstbin, int width, double baseline){
  double val=h->Integral(firstbin,firstbin+width-1)-baseline*width;
  p->Fill(val);
  return val;
}


void setupPicoscope(ps5000a &dev, chRange range, int samples, int nbuffers) {

  dev.open(picoscope::PS_12BIT);
  dev.setChCoupling(picoscope::A, picoscope::DC);
  dev.setChRange(picoscope::A, range);
  dev.enableChannel(picoscope::A);

  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1);
  //dev.setSimpleTrigger(EXT, 18000, trgFalling, 0, 0);
  dev.setSimpleTrigger(EXT, -5000, trgFalling, 0, 0);//triggering off laser
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();  
}

bool timeThat=false;   //can set timer to close TCanvas

int main(int argc, char **argv) {

  ps5000a dev;
  chRange range = PS_50MV;
  int samples = 150;    // number of ADC samples per waveform (300 ns)
  int nbuffers=10000;   // number of waveforms per capture cycle
  int nrepeat=5;        // number of capture cycles
  int nsave=20;         // number of waveforms to save for samples
  int x0=1;             // starting bin for background integral
  int xlow=samples/2;   // starting bin for pulse integral
  int xwid=15;          // number of bins to integrate
  const int nBaseline=20; // number of bins to use to estimate baseline
  int opt;
  bool findWindow=false;  // find the integration window automatically\

  bool quit=false;
  bool quiet=false;
  TString outfn="pulsed.root";
  
  while ((opt = getopt(argc, argv, "b:c:S:n:o:hq0xaTw:z:")) != -1) {
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
      x0=atoi(optarg);
      break;
    case 'x':
      xlow=atoi(optarg);
      break;
    case 'T':
      timeThat=true;
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
      fprintf(stderr, "Usage: %s ",argv[0]);
      fprintf(stderr, "-b nbuffers[10000] -c nrepeat[5] -w nsave[20] -z0 starting bin for background integration[0]");
      fprintf(stderr, "-o output[pulsed.root]\n");
      fprintf(stderr, "-a automatically find integration window\n");
      fprintf(stderr, "-T times out the TCanvas to close it");
      fprintf(stderr, "-x starting bin of pulse integral[50]\n");
      fprintf(stderr, "-q exit when finished\n");	    
      exit(EXIT_FAILURE);
    }
  }
  TApplication theApp("App", &argc, argv, 0, -1);
  
  setupPicoscope(dev, range, samples, nbuffers);
  int mvScale = autoRange(dev,50000);
  TH1F *hRange    = new TH1F("hRange","Picoscope range setting",2,0,2);
  hRange->SetBinContent(1,mvScale);
  hRange->SetBinContent(2,1.0*mvScale/32767);
  dev.setCaptureCount(nbuffers);
   
  vector<TH1F *> wavesSave;
  
  TH2F *hpersist=new TH2F("hpersist","Persistence Display",samples,
			  0,samples,512,-2048,32768);  // ADC counts are reported in steps of 20 units
  hpersist->GetXaxis()->SetTitle("Sample time [2ns/div]");
  
  
  TH1F *hsamp = new TH1F("hsamp","Samples", samples, 0, samples);  // use for reading buffer
  TProfile *hprof=new TProfile("hprof","Wave data profile",samples,0,samples);
  TProfile *hBuf0;
  TH1F* hpulses1=new TH1F("hpulses1","Pulse area distribution",2500,-40000,460000);
  hpulses1->GetXaxis()->SetTitle("Pulse Area");
  TH1F* hpulses0=new TH1F("hpulses0","Pulse area distribution",2500,-40000,460000);
  hpulses0->SetLineColor(kRed);
  double baseline=0;
  
  for (int iblock = 0; iblock < nrepeat; iblock++) {
    cout << "Capturing Block:" << iblock << endl;     
    dev.captureBlock();

    vector <vector<short> > &data = dev.getWaveforms();

    if (iblock==0) {  // use first block to figure out the pulse integration window
      for (auto &waveform : data) {
	hsamp->Reset();
	for (int i = 0; i < samples; i++) {  // samples = waveform.size()
	  hprof->Fill(i,-1*waveform[i]);
	  hsamp->SetBinContent(i, -1*waveform[i]);
	}
	if (wavesSave.size()<nsave) {
	  int nsav=wavesSave.size()+1;
	  wavesSave.push_back((TH1F*)hsamp->Clone(TString::Format("wave%02d",nsav)));
	}
      }
      if (findWindow){
	baseline=hprof->Integral(x0,x0+nBaseline-1)/nBaseline;
	int maxbin=hprof->GetMaximumBin();
	double dYmax=hprof->GetBinContent(maxbin)-baseline;
	// search find locations of ~10% peak height
	int left=hprof->FindFirstBinAbove(dYmax/10+baseline);
	int right=hprof->FindLastBinAbove(dYmax/10+baseline);
	// set window to this "width"
	xwid=(right-left);  // consider adding a couple bins to account
	xlow=left;          // for pulse shape
	if (xlow <= x0+nBaseline-1)
	  cout << "WARNING: pulse is too close to background area, move to the right" << endl;
      }
      hBuf0 = (TProfile*) hprof->Clone("hbuf0");
    }
    // 2nd loop over buffer 0 is inefficient, can be improved?
    for (auto &waveform : data) {
      hsamp->Reset();
      for (int i = 0; i < samples; i++) {  // samples = waveform.size()
	if (iblock>0) hprof->Fill(i,-1*waveform[i]);  // don't duplicate entries
	hpersist->Fill(i, -1*waveform[i]);  
	hsamp->SetBinContent(i, -1*waveform[i]);
      }

      //Note! Currently this is in ADC  counts, not anything else
      PHD(hsamp,hpulses1,xlow,xwid,baseline);  // integrals are inclusive over bins limits
      double val=PHD(hsamp,hpulses0,x0,xwid,baseline);
      //cout << val << endl;
    }
  }

  dev.close();
  
  TFile *f = new TFile(outfn, "RECREATE");

  hBuf0->Write();
  std::cout << "Writing output to: " << outfn << std::endl; 
  TCanvas *tc=new TCanvas("tc","Pulse heights",1200,400);
  tc->Divide(3,1);
  gStyle->SetOptStat(0);

  // save sample waveforms
  for (int i=0; i<wavesSave.size(); i++) wavesSave[i]->Write();
  
  //Timer to turn off so can run multiple scans !!!not done!!!
  if (timeThat==1){
    TTimer *timer = new TTimer();
    timer->Connect("Timeout()","TCanvas",tc,"Close()");
    timer->Start(2000,kTRUE);
  }
  
  // projections of persistance histogram to show counts vs threshold
  tc->cd(1)->SetLogy();
  TH1D *hIntime = hpersist->ProjectionY("_py1",xlow,xlow+xwid);
  hIntime->SetName("hIntInT");
  hIntime->SetTitle("Threshold Scan;ADC");
  TH1D *hOotime = hpersist->ProjectionY("_py0",x0,x0+xwid);
  hOotime->SetName("hIntOoT");
  hOotime->SetTitle("Threshold Scan;ADC");
  if (hOotime->GetMaximum() > hIntime->GetMaximum()) hIntime->SetMaximum(hOotime->GetMaximum()*1.10);
  hIntime->Write();
  hOotime->Write();
  hIntime->Smooth();
  hOotime->Smooth();
  hIntime->DrawCopy("");
  hOotime->SetLineColor(kRed);
  hOotime->DrawCopy("same");
  TH1F *hlimits=new TH1F("hLimits","xlow,xwid,x0",3,0,3);
  hlimits->SetBinContent(1,xlow);
  hlimits->SetBinContent(2,xwid);
  hlimits->SetBinContent(3,x0);
  hlimits->Write();
  TH1F *hbaseline=new TH1F("hbaseline","baseline",1,0,1);
  hbaseline->SetBinContent(1,baseline);
  hbaseline->Write();
  
  // plot the persistance histogram
  tc->cd(2);
  hpersist->DrawCopy("col");
  hprof->Scale( 0.9 * hpersist->GetYaxis()->GetXmax() / hprof->GetMaximum() );
  hprof->SetLineColor(kGray);
  hprof->SetLineWidth(3);
  hprof->DrawCopy("same");

  // pulse height distributions
  tc->cd(3);
  hpulses1->DrawCopy(); 
  hpulses0->Scale( hpulses1->GetMaximum()/hpulses0->GetMaximum() );
  hpulses0->DrawCopy("same");
  tc->Update();

  hpersist->Write();
  hpulses0->Write();
  hpulses1->Write(); 
  hprof->Write();
  hRange->Write();
  // save integration parameters
  TVectorD vInt(2);
  vInt[0]=xlow;
  vInt[1]=xwid;
  vInt.Write("vInt");
  
  f->Close();

  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  
  //cout << "Hit any ^c to exit" << endl;
  theApp.Run(true);
  
  return 0;
}



