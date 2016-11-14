#include "picoscopeInterface.h"
#include "ps6000.h"
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
void setupPicoscope(ps6000 &dev, chRange range, int samples, int nbuffers) {

  //  dev.open(picoscope::PS_12BIT);//allows for ADC conversion and resolution
  dev.open();//allows for ADC conversion and resolution
  dev.setChCoupling(picoscope::A, picoscope::DC_50R);
  dev.setChRange(picoscope::A, range);//allows us to set range on Picoscope
  dev.enableChannel(picoscope::A);

  //  dev.enableBandwidthLimit(picoscope::A); 
  dev.setTimebase(1); //400 ns
  dev.setSimpleTrigger(AUX, 18000, trgFalling, 0, 0);//When triggering off anything else
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
  for (int i =0; i < nbuffers; i++) {
    if (heights[i] > maximum){
      maximum = heights[i];
	}
  }
  return maximum;
}


LightPeaker *lPk = new LightPeaker(0);

int main(int argc, char **argv) {

  ps6000 dev;
  chRange range        = PS_100MV;     //range on picoscope, will capture amplitudes over 100pe
  int samples          = 3000;       // number of samples per waveform
  int nbuffers         = 10000;     // number of waveforms per capture cycle
  double theshLowLimit = 0;         // just so it works in LightPeaker
  double onePePeakNorm = 1.0;       // In case there isn't normalization given
  int opt;
  bool millivolts      = true;       //uses mV not ADC for analysis and output
  bool pePeakNormalization = false;  //can also do analysis in pe peaks (requires millivolts option if true)
  bool quit            = false;
  bool quiet           = false;
  TString outfn        = "lightPulseMeasurement_defaultFileName.root";
  
  while ((opt = getopt(argc, argv, "s:b:o:hq0up::n:")) != -1) {
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
    case 'u':
      millivolts=false;
      break;
    case 'p':
      pePeakNormalization = true;
      break;
    case 'n':
      onePePeakNorm=atof(optarg);
      std::cout<< "Our normalization is: " << onePePeakNorm << std::endl;
      break;
    default: /* '?' */
      fprintf(stderr, "Usage: %s",argv[0]);
      fprintf(stderr, "-b nbuffers[10000] -w bins to integrate -o output[lightPulseMeasurement_defaultFileName.root -n 1 pe normalization\n");
      fprintf(stderr, "-x starting bin of pulse integral[50]\n");
      fprintf(stderr, "-q exit when finished\n");
      fprintf(stderr, "-u use ADC instead of default mV\n");
      exit(EXIT_FAILURE);
    }
  }
  TApplication theApp("App", &argc, argv, 0, -1);


  
  // Auto range in case of overflow before taking data 
  /*
  setupPicoscope(dev, range, samples, 1000);
  dev.captureBlock();
  vector <vector<short> > &data = dev.getWaveforms();
  dev.close();
  Int_t waveSize  = data.front().size();
  int waveMax=0;
  for (auto &waveform : data) {
    for (int i = 0; i < waveform.size(); i++) {
      if (-1*waveform[i] > waveMax) waveMax=-1*waveform[i];
    }
  }
  if (waveMax>32000) range=(chRange)(range+1);
  */
  
  // Take the data
  setupPicoscope(dev, range, samples, nbuffers); 
  dev.captureBlock();
  vector <vector<short> > &data = dev.getWaveforms();
  //data = dev.getWaveforms();
  Int_t waveSize  = data.front().size();
  dev.close();
  
  float timebase = dev.timebaseNS();

  //Get's peak heights and makes array of said heights
  Float_t *peakHeight;
  
  int buffNum     = 0;
  TSPECTFLOAT *heightsOfPeaks     = new TSPECTFLOAT[nbuffers];
  TSPECTFLOAT *normHeightsOfPeaks = new TSPECTFLOAT[nbuffers];
  TH1F *hdata     = new TH1F("hdata","The Pulses", waveSize, 0, waveSize);
  TH1F *hdataMv   = new TH1F("hdataMv","Pulses in mV", waveSize, 0, waveSize);
  TH1F *hsum      = new TH1F("hsum","Sum of all pulses, mV measurements", waveSize, 0, waveSize);
  TProfile *hprof = new TProfile("hprof","Profile of all pulses, mV measurements", waveSize, 0, waveSize,"S");
  TH1F *h1PePeak  = new TH1F("h1PePeak","Height in mV of 1Pe peak",1,-1,1);
 
  for (auto &waveform : data) {
    buffNum++;
    
    if (buffNum % 1000 == 0){
      std::cout << "Processing Buffer: " << buffNum << std::endl;
    }
    
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

    heightsOfPeaks[buffNum-1]     = peakHeight[0];
    normHeightsOfPeaks[buffNum-1] = peakHeight[0]/onePePeakNorm;

    hdata->Reset();
    hdataMv->Reset();
  }

  TFile *f = new TFile(outfn, "RECREATE");  

  //Finds something to base the histogram limits off of 
  int repPeakHeight = baseDistLimits(heightsOfPeaks,nbuffers);

  //Prepare peak distribution histogram
  TH1F *hPeaksDist = new TH1F("hPeakDist","Distribution of Peak heights; mV",repPeakHeight,repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10) ;

  for (int i = 0; i < nbuffers; i++) {
    hPeaksDist->Fill(heightsOfPeaks[i]);
      }
  
  TH1F *hSigmaOMeanPe = new TH1F("hSigmaOMeanPe","Sigma/Mean Value with 1 pe normalization",1,-1,1);
  TH1F *hMeanPePeaks  = new TH1F("hMeanPePeaks","Mean Value in 1 Pe peaks",1,-1,1);//taken from a fit of peaks normalized to 1 pe peaks
  TH1F *hSigmaPeakHeightPe = new TH1F("hSigmaPeakHeightPe","Sigma for Peak Distribution",1,-1,1);
    
  //Creates a complete pe peak normalized regime
  if (pePeakNormalization == true && millivolts == true){
    int repNormHeight  = baseDistLimits(normHeightsOfPeaks,nbuffers);
    TH1F *hPeaksDistPe = new TH1F("hPeakDistPe","Distribution of Peak heights, normalized to 1 Pe peak",repNormHeight,repNormHeight-repNormHeight*.75,repNormHeight+repNormHeight*0.10) ;

    for (int i = 0; i < nbuffers; i++) {
    hPeaksDistPe->Fill(normHeightsOfPeaks[i]);
      }

    TF1 *hPeaksDistPeFit = new TF1("hPeaksDistPeFit","gaus",repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10);
    hPeaksDistPe->Fit("hPeaksDistPeFit","","",repNormHeight-repNormHeight*.75,repNormHeight+repNormHeight*0.10);

    double meanPeakHeightPe  = hPeaksDistPeFit->GetParameter(1);
    std::cout<< "The mean of the pe peak distribution is :" << meanPeakHeightPe<<std::endl;
    double sigmaPeakHeightPe = hPeaksDistPeFit->GetParameter(2);
    double sigmaOmeanPe      = sigmaPeakHeightPe/meanPeakHeightPe;

    hSigmaOMeanPe->SetBinContent(1,sigmaOmeanPe);
    hMeanPePeaks->SetBinContent(1,meanPeakHeightPe);
    hSigmaPeakHeightPe->SetBinContent(1,sigmaPeakHeightPe);

    hSigmaOMeanPe->Write();
    hMeanPePeaks->Write();
    hSigmaPeakHeightPe->Write();
    hPeaksDistPe->Write();
    
  }

  //Everything is in millivolts
  TH1F *hSigmaOMean      = new TH1F("hSigmaOMean","Sigma/Mean Value",1,-1,1);
  TH1F *hMeanPeakHeight  = new TH1F("hMeanPeakHeight","Mean Value in mV;; mV",1,-1,1);
  TH1F *hSigmaPeakHeight = new TH1F("hSigmaPeakHeight","Sigma for Peak Distribution",1,-1,1);

  TF1 *hPeaksDistFit = new TF1("hPeaksDistFit","gaus",repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10);
  hPeaksDist->Fit("hPeaksDistFit","","",repPeakHeight-repPeakHeight*.75,repPeakHeight+repPeakHeight*0.10);
   
  double meanPeakHeight  = hPeaksDistFit->GetParameter(1);
  double sigmaPeakHeight = hPeaksDistFit->GetParameter(2);
  double sigmaOmean      = sigmaPeakHeight/meanPeakHeight;
  
  hSigmaOMean->SetBinContent(1,sigmaOmean);
  hMeanPeakHeight->SetBinContent(1,meanPeakHeight);
  hSigmaPeakHeight->SetBinContent(1,sigmaPeakHeight);
  h1PePeak->SetBinContent(1,onePePeakNorm);

  
  TCanvas *tc=new TCanvas("tc","Information about Pulses",1600,400);
  tc->Divide(4,1);
  gStyle->SetOptStat(0);
  tc->cd(1);
  hPeaksDist->DrawCopy();
  tc->cd(2);
  hMeanPeakHeight->DrawCopy();
  tc->cd(3);
  hMeanPePeaks->DrawCopy();
  tc->cd(4);
  hprof->DrawCopy();
  
  hPeaksDist->Write();
  hSigmaOMean->Write();
  hMeanPeakHeight->Write();
  hSigmaPeakHeight->Write();
  hsum->Write();
  hprof->Write();
  h1PePeak->Write();
  
  
  f->Close();

  std::cout<< "Close TCanvas: Information about Pulses to exit" << std::endl;
  tc->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  theApp.Run(true);
  
  return 0;
}



