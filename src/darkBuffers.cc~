#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TGraph.h" 
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "utils.h"
#include <TApplication.h>
#include <getopt.h>
#include <stdio.h>

using namespace picoscope;

void fitExp(TH1F *h, TString opt="L"){
  h->Fit("expo",opt,"",h->GetBinCenter(8),h->GetXaxis()->GetXmax());
  TF1 *fun2=new TF1(*(h->GetFunction("expo")));
  fun2->SetName("expo2");
  fun2->SetRange(h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  fun2->SetLineStyle(2);
  fun2->Draw("same");
  h->GetListOfFunctions()->Add(fun2);
}


int main(int argc, char **argv) {

  TString outfn="darkBuffers.root";
  int samples = 40000;
  int nbuffers = 50;
  int opt;
  bool quit=false;
  bool quiet=false;
  while ((opt = getopt(argc, argv, "b:n:o:hq0")) != -1) {
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
    case 'q':   // exit when finished
      quit=true;
      break;
    case '0':   // turn off graphics
      quiet=true;  
      break;	
    case 'h':
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-s nsamples[40000] -b nbuffers[50] -o output[darkBuffers.root]\n");
      exit(EXIT_FAILURE);
    }
  }

  TApplication theApp("App", &argc, argv);


  ps5000a dev;
  chRange range = PS_20MV;
  dev.open(picoscope::PS_12BIT);
  dev.setChCoupling(picoscope::A, picoscope::DC);
  dev.setChRange(picoscope::A, range);
  dev.enableChannel(picoscope::A);

  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1);
  //dev.setSimpleTrigger(EXT, 18000, trgRising, 0, 0); 
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
  //  dev.setCaptureCount(10000);
  dev.setCaptureCount(nbuffers);
  dev.prepareBuffers();
  dev.captureBlock(); 
  dev.close();

  float timebase = dev.timebaseNS();
  std::cout << "Timebase: " << timebase << std::endl; 
  vector <vector<short> > data = dev.getWaveforms();

  vector<float> graphWaveform(data[0].size());
  vector<float> graphtime(graphWaveform.size());
  float timebaseStart = timebase*samples/2*-1;

  for (int i = 0; i < data[0].size(); i++) {
    graphtime[i] = timebaseStart+i*timebase;
  }


  TFile f(outfn, "RECREATE");
  std::cout << "writing waveforms..." << std::endl;
  TH1F *hist = NULL;

  TH1F *dT  = new  TH1F("dT", "Time Steps [ns]", 1,0,1);
  dT->Fill(0.0, timebase);
  TH1F *dV = new TH1F("dV", "Voltage Steps[mV]", 1,0,1);
  dV->Fill(0.0, dev.adcToMv(1, range));
  std::cout << "dV:" << dev.adcToMv(1, range) << std::endl; 
  dT->Write();
  dV->Write();

  delete dT;
  delete dV;


  // delta time distribution
  TH1F *hdTime=new TH1F("hdTime","Delta times;x [2 ns]",101,-2.5,502.5);
  // pulse height distribution
  TH1F* hpeaks=new TH1F("hpeaks","Peaks",200,0,15000);
  TCanvas *tc=new TCanvas("tc","Samples",50,20,1200,400);
  TCanvas *tc1=new TCanvas("tc1","Peaks and Time distribution",0,450,1200,400);
  tc1->Divide(2,1);
  gStyle->SetOptStat(0);
  
  DarkPeaker *dPk = new DarkPeaker();

  bool first=false;
  int nbuf=0;
  for (auto &waveform : data) {
    nbuf++;
    std::cout << "Processing buffer: " << nbuf << std::endl;
    hist = new TH1F("pulses", "pulses", waveform.size(), 0, waveform.size());
    for (int i = 0; i < waveform.size(); i++) {
      hist->SetBinContent(i, -1*waveform[i]);
    }
    if (first) {
      int iymax=(int)hist->GetMaximum();
      iymax=iymax*1.1;
      iymax-=iymax%1000;
      hpeaks->SetBins(200,0,iymax);
      first=false;
    }
    
    dPk->SetBuffer(hist,timebase);
    dPk->AnalyzePeaks();
    
    // retrieve peak heights
    Float_t *yvals = dPk->GetBkgdCorrectedY();
    int npeaks=dPk->GetNPeaks();

    // Fill pulse height histogram
    for (int i=0;i<npeaks;i++) hpeaks->Fill(yvals[i]);
    
    // fill delta time distro
    Float_t *deltaX = dPk->GetDeltaX();
    for (int i=0;i<npeaks-1;i++) hdTime->Fill(deltaX[i]);

    // draw samples buffer with peaks and background
    tc->cd();
    hist->DrawCopy();
    TH1F* bkg=dPk->GetBackground();
    bkg->SetLineColor(kGreen+2);
    bkg->SetLineWidth(3);
    bkg->Draw("same");
    tc->Update();

    // draw pulse analysis plots
    tc1->cd(1);
    hpeaks->DrawCopy();
    tc1->cd(2);
    fitExp(hdTime,"LQ");
    hdTime->DrawCopy();
    tc1->Update();

    hist->Write(); 
    delete hist;
  }
  hdTime->Write();
  hpeaks->Write();
  fitExp(hdTime,"L");
  hdTime->DrawCopy();
  tc1->Update();
  
  // Get the average dark pulse rate
  TH1F *hRate=new TH1F("hRate","Dark Pulse Rate;;MHz",1,-1,1);
  double meanDt = -2/hdTime->GetFunction("expo")->GetParameter(1);  // Dt in [ns]
  double rate = 1 / meanDt * 1000;  // in MHz
  std::cout << "Average dark pulse rate: " << rate << std::endl;
  double par[2];
  TF1 *fcn=hdTime->GetFunction("expo");
  fcn->GetParameters(par);
  hRate->SetBinContent(1,rate);
  hRate->SetBinError(1,rate*hdTime->GetFunction("expo")->GetParError(1)/
		     par[1]);
  hRate->Write();

  // Find the afterpulsing rate
  int maxbin=hdTime->GetMaximumBin();

  // rough estimate of excess over exponential fit, can improve
  double excess=0;
  for (int i=maxbin; i<=hdTime->GetNbinsX(); i++){
    double x=hdTime->GetBinCenter(i);
    double y=hdTime->GetBinContent(i);
    double ey=hdTime->GetBinError(i);
    if ( (y-fcn->Eval(x)) < ey/2 ) break;
    excess += y-fcn->Eval(x);
  }
  double arate=excess/hdTime->Integral();
  /*
  int xmax=hdTime->GetXaxis()->GetXmax();
  int xmin=hdTime->GetBinLowEdge(maxbin);
  double histInt=hdTime->Integral(maxbin,hdTime->GetNbinsX()) *
    hdTime->GetBinWidth(maxbin);
  double funcInt=exp(par[0])/par[1] * ( exp(par[1]*xmax) - 1 );
  double arate = (histInt-funcInt)/funcInt;
  */
  
  TH1F *hAp=new TH1F("hAp","After Pulse Rate",1,-1,1);
  hAp->SetBinContent(1,arate);
  hAp->SetBinError(1,arate/sqrt(excess)); // rough estimate, can improve
  hAp->Write();
  
  f.Close();

  std::cout << "===============================" << std::endl;
  std::cout << "Fit dark pulse rate: " << rate << std::endl;
  std::cout << "After pulse rate: " << arate << std::endl;
  std::cout << "===============================" << std::endl;
  
  if (quit) return 0;
  std::cout << "Hit any ^c to exit" << std::endl;
  theApp.Run(true);
  return 0;


}

