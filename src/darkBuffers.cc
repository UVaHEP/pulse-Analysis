#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TGraph.h" 
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "utils.h"
#include "eventHandler.h"
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

void fit2Exp(TH1F *h, TString opt="L"){
  // fit for DCR using long delta time tail
  h->Fit("expo",opt,"",h->GetBinCenter(8),h->GetXaxis()->GetXmax());
  TF1 *fun1=new TF1(*(h->GetFunction("expo")));
  fun1->SetName("expoDCR");
  double a1,b1;  //fit parameters for DCR
  a1=fun1->GetParameter(0);
  b1=fun1->GetParameter(1); 
  
  // estimate parameter for afterpulse component
  int maxbinhdTime;
  double apNorm;
  double xmaxStarthdTime;
  maxbinhdTime = TMath::Min(3,h->GetMaximumBin());
  xmaxStarthdTime = h->GetBinCenter(maxbinhdTime);
  apNorm=h->GetMaximum()-fun1->Eval(xmaxStarthdTime);
  double a2,b2; // fit parameters for APulse
  a2=TMath::Log(apNorm);
  b2=-1/6.; //HACK!
  TF1 *afterPulseFit = new TF1("afterPulseFit","exp([0]+[1]*x)+exp([2]+1*[3]*(x-[4]))",xmaxStarthdTime, h->GetXaxis()->GetXmax());
  afterPulseFit->SetParameters(a1,b1,a2,b2,xmaxStarthdTime);
  afterPulseFit->FixParameter(4,xmaxStarthdTime);
  h->Fit("afterPulseFit",opt,"",xmaxStarthdTime,h->GetXaxis()->GetXmax());
  a1=fun1->GetParameter(0);
  b1=fun1->GetParameter(1);
  //a2=fun1->GetParameter(2);
  fun1->SetParameters(a1,b1);

  fun1->SetRange(h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
  fun1->SetLineStyle(2);
  fun1->Draw("same");
  h->GetListOfFunctions()->Add(fun1);
  h->GetListOfFunctions()->Add(afterPulseFit);
}

void setupScope(ps5000a &dev,   chRange &range, int samples) { 
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

void acquireBuffers(ps5000a &dev, vector <vector<short> > &data){
  dev.prepareBuffers();
  dev.captureBlock(); 
  data = dev.getWaveforms();
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

int main(int argc, char **argv) {

  TString outfn="darkBuffers.root";
  int samples = 40000;
  int nbuffers = 50;
  double peThreshold=0;
  int opt;
  bool quit=false;
  bool quiet=false;
  bool userThreshold = false; 
  chRange range = PS_20MV;
  TString fileToOpen;
  while ((opt = getopt(argc, argv, "s:b:o:P:R:f:uhq0")) != -1) {
    switch (opt) {
    case 's':
      samples = atoi(optarg);
      break;
    case 'b':
      nbuffers=atoi(optarg);
      break;
    case 'u':
      // Capture first waveform, then allow the user to select a threshold via GUI
      userThreshold = true;
      break;
    case 'o':
      outfn = optarg;
      std::cout<<"set output file " << outfn<<std::endl;
      break;
    case 'P':
      peThreshold=atof(optarg);
      std::cout<<"1PE value " << peThreshold<<std::endl;
    case 'R':
      if (TString(optarg)=="PS_50MV") {
	range = PS_50MV;
	std::cout<<"setting range to PS_50MV"<<std::endl;
      }
      else if (TString(optarg)=="PS_20MV") {
	std::cout<<"setting range to PS_20MV"<<std::endl;
      }
      else std::cout<<"Unknown range, defaulting to PS_20MV"<<std::endl;
      break;
    case 'q':   // exit when finished
      quit=true;  // not implemented
      break;
    case '0':   // turn off graphics
      quiet=true;    // not implemented
      break;
    case 'f':   // loads in a file, ignores picoscope stuff
      fileToOpen = optarg;
      std::cout<<"reading data from: " <<fileToOpen<<std::endl;      
      break;
    case 'h':
    default: /* '?' */
      fprintf(stderr, "\nUsage: %s [options]\n",argv[0]);
      fprintf(stderr, " -s nsamples[40000] : number of samples per buffer\n");
      fprintf(stderr, " -b nbuffers[50] : number of buffers\n");
      fprintf(stderr, " -u : use GUI to select 1PE threshold, default is auto threshold\n");
      fprintf(stderr, " -o output[darkBuffers.root] : Output filename\n");
      fprintf(stderr, " -R Range[PS_20MV] : Voltage range selection [PS_20MV,PS_50MV]\n");
      fprintf(stderr, " -P [ADC] : User setting to 1PE threshold in ADC counts\n");
      fprintf(stderr, " -f filename : do not acquire data, process data from file\n");
      exit(EXIT_FAILURE);
    }
  }

  TApplication theApp("App", &argc, argv);

  ps5000a dev;
  vector <vector<short> > data;
  float timebase;
  
  if (fileToOpen.Length()==0){  
    dev.open(picoscope::PS_12BIT);
    setupScope(dev, range, samples); 
    timebase = dev.timebaseNS();
    std::cout << "Timebase: " << timebase << std::endl;
    
    // run GUI to pick 1PE threshold
    if (userThreshold) {
      peThreshold = userThresholdFn(dev, samples, theApp);
    }
    std::cout << "pe Threshold: " << peThreshold << std::endl; 
  
    dev.setCaptureCount(nbuffers);
    acquireBuffers(dev,data);
    //dev.prepareBuffers();
    //dev.captureBlock(); 
    dev.close();
    //data = dev.getWaveforms();
  }
  else { // reaf buffers from file
    nbuffers = 0;
    TFile infile(fileToOpen);
    TIter nextkey(gDirectory->GetListOfKeys());
    TKey *key;
    vector<short> buf;
    timebase = ((TH1F*)(infile.Get("dT")))->GetBinContent(1);
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

  // We use these histograms to write results to the TFile
  // delta time distribution
  TH1F *hdTime=new TH1F("hdTime","Delta times;x [2 ns]",101,-2.5,502.5);
  // pulse height distribution
  float maxPeakRange=35000; // HACK!!!!
  TH1F* hpeaks=new TH1F("hpeaks","Peaks",200,0,maxPeakRange);
  //Threshold for counted peaks, maximum needs to be dynamic
  TH1F *hpeakthresh=new TH1F("hpeakthresh","Scan of peaks;Threshold",500,0,maxPeakRange);
  TCanvas *tc=new TCanvas("tc","Samples",50,20,1200,400);
  TCanvas *tc1=new TCanvas("tc1","Peaks and Time distributions",0,450,1200,400);
  tc1->Divide(3,1);
  gStyle->SetOptStat(0);
  int totalPeaks=0;
  
  DarkPeaker *dPk = new DarkPeaker(peThreshold/2);

  bool first=true;
  int nbuf=0;
  
  for (auto &waveform : data) {
    nbuf++;
    std::cout << "Processing buffer: " << nbuf << std::endl;
    //std::cout << "Bins in hist: " << waveform.size() << std::endl;
    hist = new TH1F("pulses", "pulses", waveform.size(), 0, waveform.size());
    for (int i = 0; i < waveform.size(); i++) {
      hist->SetBinContent(i, -1*waveform[i]);
    }

    dPk->SetBuffer(hist,timebase);
    dPk->AnalyzePeaks();
    
    if (first) {
      int iymax=(int)hist->GetMaximum();
      iymax=iymax*1.1;
      iymax-=iymax%1000;
      hpeaks->SetBins(200,0,iymax);
      dPk->GetHdist()->Write();  // save a copy of data used for noise estimation
      first=false;
    }
    
    // retrieve peak heights
    TSPECTFLOAT *yvals = dPk->GetBkgdCorrectedY();
    int npeaks=dPk->GetNPeaks();
    totalPeaks+=npeaks;
    
    // Fill pulse height histogram, pulse threshold histogram
    for (int i=0;i<npeaks;i++){
      hpeaks->Fill(yvals[i]);
      for (int b = 1; b<=hpeakthresh->GetNbinsX();b++){
	if (yvals[i]>hpeakthresh->GetBinLowEdge(b))
	  hpeakthresh->SetBinContent(b,hpeakthresh->GetBinContent(b)+1);	
      }
    }
    
    // fill delta time distro
    TSPECTFLOAT *deltaT = dPk->GetDeltaX();
    for (int i=0;i<npeaks-1;i++) hdTime->Fill(deltaT[i]);

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
    //peak threshold histogram
    tc1->cd(2);
    hpeakthresh->DrawCopy();
    tc1->cd(3);
    fit2Exp(hdTime,"LQ");  
    hdTime->DrawCopy();
    tc1->Update();

    hist->Write(); 
    delete hist;
  }
  
  hdTime->Write();
  hpeaks->Write();
  hpeakthresh->Write();
  fit2Exp(hdTime,"L");
  hdTime->DrawCopy();
  tc1->Update();

  //Total time for all buffers
  double timeTotal = nbuffers*samples*timebase;
  std::cout <<"Number of buffers: "<<nbuffers<<std::endl;
  
  // Dark pulse rate information
  TH1F *hRate=new TH1F("hRate","Dark Pulse Rate;;MHz",2,-1,2); // bin1 DCR fit, bin2 DCR count
  TH1F *hCount=new TH1F("hCount","Dark Pulse Count;;",1,-1,1);
  TH1F *hTtot=new TH1F("hTtot","Total time of samples;;",1,-1,1);
  TH1F *hAp=new TH1F("hAp","After Pulse Rate",1,-1,1);
  // simple counts
  hCount->SetBinContent(1,totalPeaks);
  hTtot->SetBinContent(1,timeTotal);
  double dcr=totalPeaks/timeTotal/1000000; // in MHz, not corrected for AP
  hRate->SetBinContent(2,totalPeaks/timeTotal/1000000);
  hRate->SetBinError(2,TMath::Sqrt(totalPeaks)/timeTotal);
  // extracted from exponential fit
  TF1 *tf_expoDCR=hdTime->GetFunction("expoDCR");
  double par[2];
  tf_expoDCR->GetParameters(par);
  double meanDt = -2/par[1];  // Delta_t in [ns]
  double dcrFit = 1 / meanDt * 1000;  // in MHz
  double dcrErr = dcrFit*tf_expoDCR->GetParError(1)/par[1];
  hRate->SetBinContent(1,dcrFit);
  hRate->SetBinError(1,dcrErr);

  std::cout << "Average dark pulse rate (uncorrected): " << dcr << std::endl;

  double aPrate=0;
  if (dcrFit<0 || TMath::Abs(dcrErr/dcrFit)>0.25){
    std::cout << "Warning Fit yields negative DCR or large error" << std::endl;
    std::cout << "Afterpulse calculation skipped" << std::endl;
  }
  else {  // calculate afterpulse rate
    //Written in part by Grace E. Cummings, 30 July 2016
    TF1 *tf_apFcn = hdTime->GetFunction("afterPulseFit");
    double maxBinTime = hdTime->GetMaximumBin();
    double xmaxBinTime = hdTime->GetBinCenter(maxBinTime);
    double xmaxTime = hdTime->GetXaxis()->GetXmax();

    //ap Rate calculated as ratio of fit of afterpulses to fit of dark counts
    //Fit of afterpulses is the afterPulseFit-expoDCR
    double aPFint = tf_apFcn->Integral(xmaxBinTime,xmaxTime);
    double expDCRint = tf_expoDCR->Integral(xmaxBinTime,xmaxTime);
    double binWidthhdTime = hdTime->GetBinWidth(1);
    
    aPrate = (aPFint-expDCRint)/expDCRint;
    
    hAp->SetBinContent(1,aPrate);
    //Error in afterpulse rate. Sqrt(excess)~sigma of excess
    double excess = (aPFint-expDCRint)/binWidthhdTime;
    hAp->SetBinError(1,sqrt(excess)/(expDCRint/binWidthhdTime));    
  }
    
  hRate->Write();
  hCount->Write();
  hTtot->Write();
  hAp->Write();
    
		     
  //Find the crosstalk fraction
  //written by Grace E. Cummings, 26 July 2016
  int maxPeaksBin = hpeaks->GetMaximumBin(); 
  double maxPeaksHeight = hpeaks->GetBinContent(maxPeaksBin); 
  double xminPeaks = hpeaks->GetXaxis()->GetXmin();
  double xmaxPeaks = hpeaks->GetXaxis()->GetXmax();
  double xendPeaksFit = hpeaks->GetBinCenter(hpeaks->GetNbinsX());
  double sigmaGuessBin = 5+maxPeaksBin;//parameter for fit backup
  
  //Find point where histogram value drops bellow 7% of max. Used to end fit
  for (int i = maxPeaksBin; i<=hpeaks->GetNbinsX();i++){
    if (hpeaks->GetBinContent(i)<=0.07*maxPeaksHeight){
      xendPeaksFit = hpeaks->GetBinCenter(i);
      break;
    }
  }

  //Creates Gaussian Fit and plots
  double sigmaGuess = (hpeaks->GetBinCenter(sigmaGuessBin))-hpeaks->GetBinCenter(maxPeaksBin);
  TF1 *peaksFit = new TF1("peaksFit","[0]*exp(-0.5*((x-[1])/[2])**2)",xminPeaks,xmaxPeaks);
  peaksFit->SetParameter(0,hpeaks->GetMaximum());
  peaksFit->SetParameter(1,hpeaks->GetBinCenter(maxPeaksBin));
  peaksFit->SetParameter(2,sigmaGuess);
  TCanvas *tc2 =new TCanvas("tc2","Samples",50,20,400,400);
  tc2->cd();
  hpeaks->Fit("peaksFit","","",0,xendPeaksFit);
  hpeaks->DrawCopy();
  double meanhPeaks = peaksFit->GetParameter(1);
  double meanhPeakmV = dev.adcToMv(meanhPeaks,range);
  double uncMeanhPeaks = peaksFit->GetParError(1);
  double sigmahPeaks = peaksFit->GetParameter(2);
  double uncSigmahPeaks = peaksFit->GetParError(2);
  double sigmaOmeanhPeaks = sigmahPeaks/meanhPeaks;

  //calculate Crosstalk fraction. Crosstalk peaks are what wasn't fitted
  double crosstalkFraction = 0;
  double crosstalkPeaks = 0;
  double ourcrosstalkPeaks = 0;
  double ourcrosstalkFraction = 0;
  
  ourcrosstalkPeaks = (totalPeaks-peaksFit->Integral(0,xmaxPeaks)/hpeaks->GetBinWidth(1));
  ourcrosstalkFraction = ourcrosstalkPeaks/totalPeaks;

  for (int i = 1; i<=hpeaks->GetNbinsX();i++){
    if (hpeaks->GetBinCenter(i)>=meanhPeaks*1.5){
    crosstalkPeaks+=hpeaks->GetBinContent(i);
    }
  }
  crosstalkFraction = crosstalkPeaks/totalPeaks;
  
  //Save the Crosstalk Fraction (lightspin way)
  TH1F *hCrossTalk = new TH1F("hCrossTalk","Crosstalk Fraction",1,-1,1);
  hCrossTalk->SetBinContent(1,crosstalkFraction);
  //Error in Crosstalk Fraction. Sqrt(multi-photon peaks)~sigma of the crosstalk/totalPeaks
  hCrossTalk->SetBinError(1,sqrt(crosstalkPeaks)/totalPeaks);
  hCrossTalk->Write();

  //Save crosstalk fraction our way
  TH1F *hOurMethodCrossTalk = new TH1F("hOurMethodCrossTalk","CrossTalk Fraction by subtracting 1Pe Peak",1,-1,1);
  hOurMethodCrossTalk->SetBinContent(1,ourcrosstalkFraction);
  hOurMethodCrossTalk->Write();
  
  //Save mean of 1PE peak
  TH1F *h1PePeak = new TH1F("h1PePeak","1 PE Peak Mean Value",1,-1,1);
  h1PePeak->SetBinContent(1,meanhPeaks);
  h1PePeak->SetBinError(1,uncMeanhPeaks);
  h1PePeak->Write();
  TH1F *h1PePeakmV = new TH1F("h1PePeakmV","1 PE Peak Mean Value, mV",1,-1,1);
  h1PePeakmV->SetBinContent(1,meanhPeakmV);
  h1PePeakmV->Write();
  
  //Save sigma of 1Pe peak
  TH1F *h1PePeakSigma = new TH1F("h1PePeakSigma","1 PE Peak Sigma Value",1,-1,1);
  h1PePeakSigma->SetBinContent(1,sigmahPeaks);
  h1PePeakSigma->SetBinError(1,uncSigmahPeaks);
  h1PePeakSigma->Write();

  //Save sigma/1Pe peak
  TH1F *hSigmaOMean = new TH1F("hSigmaOMean","Sigma/Mean",1,-1,1);
  hSigmaOMean->SetBinContent(1,sigmaOmeanhPeaks);
  hSigmaOMean->Write();
  
  f.Close();
 
  std::cout << "===============================" << std::endl;
  std::cout << "Dark pulse rate (uncorrected): " << dcr << " MHz" << std::endl;
  std::cout << "Dark pulse fit: " << dcrFit << " MHz" << std::endl;
  std::cout << "Afterpulse probability:  " << aPrate << std::endl;
  std::cout << "Crosstalk Fraction (lightspin method): " << crosstalkFraction << std::endl;
  //std::cout << "Our crosstalk fraction method: : "<< ourcrosstalkFraction <<std::endl;
  std::cout << "1Pe Peak value (mV): " << meanhPeakmV << std::endl;
  std::cout << "===============================" << std::endl;
  
  std::cout<< "Close TCanvas: Peaks and Time distributions to exit" << std::endl;
  tc1->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  if (quit)
    return 0;
  else
    theApp.Run(true);
  
  return 0;

}

