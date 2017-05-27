
#include "picoutils.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "eventHandler.h"




void setupScope(ps5000a &dev, chRange &range, int samples) { 
  dev.setChCoupling(picoscope::A, picoscope::DC);
  dev.setChRange(picoscope::A, range);
  dev.enableChannel(picoscope::A);
  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1);
  //dev.setSimpleTrigger(EXT, 18000, trgRising, 0, 0); 
  dev.setSamples(samples); 
  dev.setPreTriggerSamples(samples/2);
  dev.setPostTriggerSamples(samples/2);
}


void autoRange(ps5000a &dev){
  vector <vector<short> > data;
  int mvRange[]={10,20,50,100,200,500,1000,2000,500};
  dev.setChRange(picoscope::A, PS_10MV);
  dev.setCaptureCount(1);
  chRange autoRange=PS_10MV;
  for ( int psRange=PS_10MV; psRange < PS_1V; psRange++ ){
    std::cout<<"Autoranging pass: " << mvRange[psRange] << "mV range" << std::endl;
    autoRange=(chRange)psRange;
    dev.prepareBuffers();
    dev.captureBlock(); 
    data = dev.getWaveforms();
    bool overThresh=false;
    for (auto &waveform : data) {
      for (int i = 0; i < waveform.size(); i++) {
	if (abs(waveform[i])>26000){      // ~ 80% of ADC range
	  overThresh=true;
	  break;
	}
      }
    } // end of buffer loop
    if (!overThresh) return;
    // set picoscope to next highest voltage
    dev.setChRange(picoscope::A, autoRange);
  }
}

// fix me to work for a specified number of buffers and remove some code below
void acquireBuffers(ps5000a &dev, vector <vector<short> > &data){
  dev.prepareBuffers();
  dev.captureBlock(); 
  data = dev.getWaveforms();
}

// read buffers from a file instead of acquiring from a scope
int readBuffers(TString filename, vector <vector<short> > &data, double &timebase, double &adc2mV){
  int nbuffers = 0;
  TFile infile(filename);
  TIter nextkey(gDirectory->GetListOfKeys());
  TKey *key;
  vector<short> buf;
  while ( (key = (TKey*)nextkey()) ) {
    TObject *obj = key->ReadObj();
    if ( TString(obj->GetName())=="pulses" ){
      TH1F *hbuf=(TH1F*)obj;
      buf.clear();
      nbuffers+=1;
      for (int i=1; i<=hbuf->GetNbinsX();i++)
	buf.push_back((short)(-1.0*hbuf->GetBinContent(i)));
      data.push_back(buf);
    }
  }
  timebase=((TH1F*)(infile.Get("dT")))->GetBinContent(1);
  adc2mV=((TH1F*)(infile.Get("dV")))->GetBinContent(1);
}


Double_t userThresholdFn(ps5000a &dev, int samples, TApplication &app) {
  dev.setCaptureCount(1);
  dev.prepareBuffers();
  dev.captureBlock();
  TCanvas *tc=new TCanvas("tc","Samples",50,20,1200,400);
  eventHandler eH(*tc); 
  TH1F* hpeaks=new TH1F("hpeaks","Peaks",200,0,27000);
  TH1F *hist = NULL;
  TH1F *hpeakthresh = new TH1F("hpeakthresh","Signal Peak Threshold",50,0,27000);

  float timebase = dev.timebaseNS();
  std::cout << "Timebase: " << timebase << std::endl; 
  vector <vector<short> > data = dev.getWaveforms();

  vector<float> graphWaveform(data[0].size());
  vector<float> graphtime(graphWaveform.size());
  float timebaseStart = timebase*samples/2*-1;

  auto waveform = data[0]; 
  
  for (int i = 0; i < waveform.size(); i++) {
    graphtime[i] = timebaseStart+i*timebase;
  }
  
  hist = new TH1F("pulses", "pulses;x[2 ns]", waveform.size(), 0, waveform.size());
  for (int i = 0; i < waveform.size(); i++) {
    hist->SetBinContent(i, -1*waveform[i]);
  }
  int iymax=(int)hist->GetMaximum();
  iymax=iymax*1.1;
  iymax-=iymax%1000;
  hpeaks->SetBins(200,0,iymax);
  hist->Draw();
  tc->Update(); 
  std::cout << "Please select the single peak threshold" << std::endl;

  tc->Connect("ProcessedEvent(Int_t,Int_t, Int_t, TObject *)", "eventHandler", &eH, "eventSlot(Int_t, Int_t, Int_t, TObject *)"); 
  tc->Connect("Closed()", "TApplication", &app, "Terminate()"); 
  app.SetReturnFromRun(true); 
  app.Run(true); 

  return eH.threshold(); 
}
