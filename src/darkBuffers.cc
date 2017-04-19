#include "picoscopeInterface.h"
#include "ps5000a.h"
#include "TFile.h"
#include "TGraph.h" 
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "utils.h"
#include "picoutils.h"
#include "TVector.h"
#include <TApplication.h>
#include <getopt.h>
#include <stdio.h>

using namespace picoscope;

class peakData {
public:
  peakData(int b, int i, double x, double h) :
    buffer(b),index(i),xpeak(x),height(h){;}
  int buffer;
  int index;
  double xpeak;
  double height;
  void Print() {
    std::cout << "(ibuf, idx, x, ycor) " << buffer << " " << index
	      << " " << xpeak << " " << height << std::endl;
  }
};

void usage(char **argv){
  fprintf(stderr, "\nUsage: %s [options]\n",argv[0]);
  fprintf(stderr, " -s nsamples[40000] : number of samples per buffer\n");
  fprintf(stderr, " -b nbuffers[50] : number of buffers\n");
  fprintf(stderr, " -a write out all buffers. 10 are written as default\n");
  fprintf(stderr, " -u : use GUI to select 1PE threshold, default is auto threshold\n");
  fprintf(stderr, " -o output[darkBuffers.root] : Output filename\n");
  fprintf(stderr, " -R Range[PS_20MV] : Voltage range selection [PS_20MV,PS_50MV,PS_100MV]\n");
  fprintf(stderr, " -P [ADC] : User setting to 1PE threshold in ADC counts\n");
  fprintf(stderr, " -f filename : do not acquire data, process data from file\n");
}

int main(int argc, char **argv) {
  TString outfn="darkBuffers.root";
  int samples = 40000;  // default samples and buffer numbers
  int nbuffers = 50;    // default number of buffers to take
  int nbuffersWrite=10; // default number of buffers to write
  int iLimitL=5;  // pulse integration limits in bin counts [5]
  int iLimitH=18; // eg peak-ilimitL to peak+iLimitH [18]
  int nbufUser=0;
  double peThreshold=0;
  int opt;
  bool quit=false;
  bool quiet=false;
  bool userThreshold = false; 
  chRange range = PS_20MV;
  TString fileToOpen;
  while ((opt = getopt(argc, argv, "s:b:o:P:R:f:uhq0a")) != -1) {
    switch (opt) {
    case 's':
      samples = atoi(optarg);
      break;
    case 'b':
      nbufUser=atoi(optarg);
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
      if (TString(optarg)=="PS_100MV") {
	range = PS_100MV;
	std::cout<<"setting range to PS_100MV"<<std::endl;
      }
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
    case 'a':
      nbuffersWrite=1000000; // write "all" buffers
    case 'f':   // loads in a file, ignores picoscope stuff
      fileToOpen = optarg;
      std::cout<<"reading data from: " <<fileToOpen<<std::endl;      
      break;
    case 'h':
    default: /* '?' */
      usage(argv);
      exit(EXIT_FAILURE);
    }
  }

  TApplication theApp("App", &argc, argv);

  // acquire or read data
  ps5000a dev;
  vector <vector<short> > data;
  double timebase, adc2mV;
  
  if (nbufUser>0) nbuffers=nbufUser;
  if (fileToOpen.Length()==0){
    dev.open(picoscope::PS_12BIT);
    setupScope(dev, range, samples); 
    timebase = dev.timebaseNS();  // despite the name this returns units of seconds
    
    // run GUI to pick 1PE threshold
    if (userThreshold) {
      peThreshold = userThresholdFn(dev, samples, theApp);
    }
    std::cout << "pe Threshold: " << peThreshold << std::endl; 
  
    dev.setCaptureCount(nbuffers);
    acquireBuffers(dev,data);
    dev.close();
    adc2mV = dev.adcToMv(1, range);
  }
  else { // read buffers from file
    nbuffers=readBuffers(fileToOpen,data,timebase,adc2mV);
    if (nbuffers==0) {
      std::cout << "No buffers found in input file.  Exiting..." << std::endl;
      return 1;
    }
    samples=data[0].size();
  }
  
  std::cout << "Timebase: " << timebase << std::endl;
  std::cout << "Samples/buffer: " << samples << std::endl;
  double timebaseNS=timebase*1e9;  
  
  // open output file and setup storage elements
  TFile f(outfn, "RECREATE");
  TH1F *hist = 0;

  TH1F *dT  = new  TH1F("dT", "Time Steps [s]", 1,0,1);
  dT->Fill(0.0, timebase);
  TH1F *dV = new TH1F("dV", "Voltage Steps [mV]", 1,0,1);
  dV->Fill(0.0, adc2mV);
  std::cout << "dV:" << adc2mV << std::endl; 
  dT->Write();
  dV->Write();

  // delta time distribution
  // 4ns bins 701,-2,2802
  // 8ns bins 351,-4,2804
  TH1F *hdTime=new TH1F("hdTime","Delta times;#Delta time [ns]",351,-4,2804);
  // pulse height distribution
  const float maxPeakRange=32768;  // roughly the max ADC value
  TH1F* hpeaks=new TH1F("hpeaks","Peaks;ADC",256,0,maxPeakRange);  // y axis is reset below
  TH1F* hintegrals=new TH1F("hintegrals","Area of isolated peaks",256,0,maxPeakRange*20); // also reset below
  TH1F* hFWHM=new TH1F("hFWHM","FWHM of peaks in bins",100,0,20);
  // 2D plot of pulse heights vs delta time
  TH2F *hdPT=new TH2F("hdPeakvTime","Peak vs Delta times;#Delta time [s];ADC counts",
		      101,-2.5,502.5,512,0,maxPeakRange); //x-axis is reset below
  double xmin=6e-9;
  double xmax=5e-6;
  int bins=(int)(TMath::Log10(xmax/xmin)*10);
  hdPT->SetBins(bins,TMath::Log10(xmin),TMath::Log10(xmax),400,0,maxPeakRange);
  BinLogX(hdPT);
  // Threshold for counted peaks
  TH1F *hpeakScan=new TH1F("hpeakScan","Threshold Scan;Threshold [ADC]",256,0,maxPeakRange);
  // Diagnostic histogram to keep track of TSpectrum search threshold based on noise
  TH1F *hsearchThresh=new TH1F("hsearchThresh","1PE search threshold;Threshold (ADC)",64,0,maxPeakRange/16);

  gStyle->SetOptStat(0);
  TCanvas *tc=new TCanvas("tc","Samples",50,20,1200,400);
  TCanvas *tc1=new TCanvas("tc1","Peaks and Time distributions",0,450,1200,400);
  TCanvas *tcPT=new TCanvas("tcPT","Peaks v. time",0,600,600,400);
  tc1->Divide(3,1);
  int totPeaks=0;
  
  DarkPeaker *dPk = new DarkPeaker(peThreshold/2);

  bool first=true;
  int nbuf=0;  // buffers processed 
  TString buftitle;
  vector<peakData> *vPeaks = new vector<peakData>; // keep track of all peak info
  FitDcrAp *dcrFitter = new FitDcrAp();
  
  // loop through buffers 
  for (auto &waveform : data) {   
    nbuf++;
    if (nbufUser>0 && nbuf>=nbufUser) break; 
    buftitle.Form("Buffer[%d];sample time [x%.1e]; ADC",nbuf,timebase);
    hist = new TH1F("pulses", buftitle, waveform.size(), 0, waveform.size());
    for (int i = 0; i < waveform.size(); i++) { // translate buffer to histogram
      hist->SetBinContent(i, -1*waveform[i]);
    }
    
    dPk->SetBuffer(hist,timebase);
    dPk->AnalyzePeaks();
    
    int npeaks=dPk->GetNPeaks();
    totPeaks+=npeaks;
    
    if (first) {
      dPk->GetHdist()->Write();  
      dPk->GetHscan()->Write();  // save a copy of data used for noise estimation in 1st buffer
      first=false;
    }
    hsearchThresh->Fill(dPk->GetSearchThreshold());

    // retrieve time ordered peak data and fill histograms
    double x,y,prev;
    TH1F *bkg=dPk->GetBackground();
    for (int i=0;i<npeaks;i++){
      dPk->GetPoint(i,x,y);
      double ycor = y - bkg->Interpolate(x);  // correct bin count for baseline shift
      vPeaks->push_back(peakData(nbuf,i,x,ycor));
      hpeaks->Fill(ycor);
      for (int b = 1; b<=hpeakScan->GetNbinsX();b++){
	if (ycor>hpeakScan->GetBinLowEdge(b))
	  hpeakScan->SetBinContent(b,hpeakScan->GetBinContent(b)+1);	
      }
      if (i>0) {
	hdTime->Fill((x-prev)*timebaseNS);
	//hdPT->Fill(x-prev,y);  // if using linear x-binning
	hdPT->Fill((x-prev)*timebase,y); // for log x-binning 
      }
      prev = x;
    }

    // now fill the histogram of integrated peaks
    // to do: add code to automatically calculate the optimal integration window
    dPk->Integrate(iLimitL,iLimitH);
    for (int i=0; i<dPk->GetNIntegrals(); i++){
      dPk->GetIntegral(i,x,y);
      hintegrals->Fill(y);
      hFWHM->Fill(dPk->GetFWHM(i));
    }
    
    // draw samples buffer with peaks and background
    tc->cd();
    hist->DrawCopy();
    bkg->DrawCopy("same");
    tc->Update();

    // draw pulse analysis plots
    tc1->cd(1);
    hpeaks->DrawCopy();

    // peak threshold scan histogram
    tc1->cd(2);
    hpeakScan->DrawCopy();

    // delta_t histogram and fits for DCR
    tc1->cd(3);
    //fit2Exp(hdTime,"LQ");
    dcrFitter->Fit(hdTime,"LQ");
    hdTime->DrawCopy();
    dcrFitter->GetApFit()->DrawCopy("same");
    dcrFitter->GetDcrFcn()->DrawCopy("same");
    dcrFitter->GetExpFit()->DrawCopy("same");
    tc1->Update();

    // 2D plot of peaks vs delta time
    tcPT->cd()->SetLogx();
    int bx,by,bz;
    int bmax= hdPT->GetMaximumBin();
    hdPT->GetBinXYZ(bmax,bx,by,bz);
    int maxYbin=3.2*by;
    //    float maxYbin=hdPT->ProjectionY()->FindLastBinAbove(0);
    hdPT->GetYaxis()->SetRange(1,maxYbin);
    hdPT->DrawCopy("col");
    tcPT->Update();
        
    if (nbuf<=nbuffersWrite) hist->Write(); 
    if (nbuffersWrite>0 && nbuf==1) dPk->GetBackground()->Write();

    delete hist;
  }
  ///////////////////////////////////////////////////////////////////////////
  // end of loop over buffers
  ///////////////////////////////////////////////////////////////////////////
  hdPT->Write();
  hdTime->Write();
  hintegrals->Write();
  hpeakScan->Write();
  hsearchThresh->Write();
  hFWHM->Write();

  //Total time for all buffers
  double timeTotal = nbuf*samples*timebase;
  TH1F *hTtot=new TH1F("hTtot","Total time of samples;;",1,-1,1);
  hTtot->SetBinContent(1,timeTotal);
  hTtot->Write();

  
  TCanvas *tc2=new TCanvas("tc2","Peaks and Time distributions (0.1PE cut)",0,450,800,400);
  tc2->Divide(2,1);


  // Fit the 1PE peak, then refine it
  hpeaks->Fit("gaus","0Q");
  TF1 *peFcn=hpeaks->GetFunction("gaus");
  double mu=peFcn->GetParameter(1);
  double sig=peFcn->GetParameter(2);
  std::cout << "** Fit to peak height distribution" << std::endl;
  hpeaks->Fit("gaus","0","",mu-2*sig,mu+2*sig);
  peFcn=hpeaks->GetFunction("gaus");
  double onePE=peFcn->GetParameter(1);  // 1PE peak in ADC, bkg corrected
  tc1->cd(1);
  hpeaks->DrawCopy();
  peFcn->DrawCopy("same");
  hpeaks->Write();

  // loop through peaks and remove everything < 0.1(5)PE cleans up grass
  vector<peakData> *vPeaks01 = new vector<peakData>;
  vector<peakData> *vPeaks05 = new vector<peakData>;
  
  for (auto &pkData : *vPeaks){
    if (pkData.height<onePE/10) continue;
    vPeaks01->push_back(pkData);
    if (pkData.height<onePE/2) continue;
    vPeaks05->push_back(pkData);
  }
  int totPeaks01= vPeaks01->size();
  int totPeaks05= vPeaks05->size();
  
  ///////////////////////////////////////////////////////////////////////////
  // now clean up the calculations/plots to using 1PE/2 threshold
  ///////////////////////////////////////////////////////////////////////////
  TH1F *hdTime01=new TH1F(*hdTime);
  hdTime01->Reset();
  hdTime01->SetName("hdTimeCut");
  hdTime01->SetTitle("Delta times (0.1PE cut)");
  TH1F *hdTime05=new TH1F(*hdTime);
  hdTime05->Reset();
  hdTime05->SetName("hdTimeCut");
  hdTime05->SetTitle("Delta times (0.5PE cut)");

  TH2F *hdPT01=new TH2F(*hdPT);
  hdPT01->Reset();
  hdPT01->SetName("hdPeakvTimeCut");
  hdPT01->SetTitle("Peaks vs Delta times (0.1PE cut)");
  hdPT01->SetMaximum(hdPT->GetMaximum());  // fix the color maps to same values

  // loop over 0.1 PE peaks
  for (int i=0; i<vPeaks01->size(); i++){
    peakData &pk=(*vPeaks01)[i];
    if (i==0) continue;
    peakData &last=(*vPeaks01)[i-1];
    double dT=0;
    if (last.buffer==pk.buffer) {
      hdTime01->Fill((pk.xpeak-last.xpeak)*timebaseNS);
      hdPT01->Fill((pk.xpeak-last.xpeak)*timebase,pk.height);
    }
  }
  
  // loop over 0.5 PE peaks
  int num1_5=0;  // for cross talk
  for (int i=0; i<vPeaks05->size(); i++){
    peakData &pk=(*vPeaks05)[i];
    if ( pk.height>=onePE*1.5 ) num1_5++;
    if (i==0) continue;
    peakData &last=(*vPeaks05)[i-1];
    double dT=0;
    if (last.buffer==pk.buffer) {
      hdTime05->Fill((pk.xpeak-last.xpeak)*timebaseNS);
    }
  }

  
  TCanvas *tcPT2=new TCanvas("tcPT2","Peaks v. time (0.1PE cut)",600,600,600,400);
  tcPT2->cd()->SetLogx();
  hdPT01->DrawCopy("col");

  
  // redo the fit for DCR
  tc2->cd(1);
  std::cout << "** Extraction of DCR and Afterpulsing 0.1 PE cut" << std::endl;
  dcrFitter->Fit(hdTime01,"L");
  hdTime01->DrawCopy();
  dcrFitter->GetApFit()->DrawCopy("same");
  dcrFitter->GetDcrFcn()->DrawCopy("same");
  dcrFitter->GetExpFit()->DrawCopy("same");
  tc2->Update();
  hdTime01->Write();
  TH1F *hRate=new TH1F("hRate","Dark Pulse Rate;;MHz",2,-1,2); // bin1 DCR fit, bin2 DCR count
  // DCR rate from counting
  double dcrCount=1.0*totPeaks05/timeTotal; // in Hz
  hRate->SetBinContent(2,dcrCount);
  hRate->SetBinError(2,TMath::Sqrt(vPeaks05->size())/timeTotal);
  // DCR extracted from exponential fit
  double dcrFit01=dcrFitter->GetDCR();
  double dcrErr01=0;
  double rateAP01=dcrFitter->GetAPrate();
  double errAP01=0;
  hRate->SetBinContent(1,dcrFit01);
  hRate->SetBinError(1,dcrErr01);
  hRate->Write();
  TH1F *hAp=new TH1F("hAp","After Pulse Rate",1,-1,1);
  hAp->SetBinContent(1,rateAP01);
  hAp->SetBinError(1,errAP01);
  hAp->Write();

  tc2->cd(2);
  std::cout << "** Extraction of DCR and Afterpulsing 0.5PE cut" << std::endl;
  dcrFitter->Fit(hdTime05,"L");
  hdTime05->DrawCopy();
  dcrFitter->GetApFit()->DrawCopy("same");
  dcrFitter->GetDcrFcn()->DrawCopy("same");
  dcrFitter->GetExpFit()->DrawCopy("same");
  double dcrFit05=dcrFitter->GetDCR();
  double dcrErr05=0;
  double rateAP05=dcrFitter->GetAPrate();
  double errAP05=0;
      
  // cross talk  
  double xTalkFrac = 1.0*num1_5/vPeaks05->size();
  TH1F *hCrossTalk = new TH1F("hCrossTalk","Crosstalk Fraction",1,-1,1);
  hCrossTalk->SetBinContent(1,xTalkFrac);
  //Error in Crosstalk Fraction. Sqrt(multi-photon peaks)~sigma of the crosstalk/totPeaks05
  hCrossTalk->SetBinError(1,sqrt(num1_5)/totPeaks05);
  hCrossTalk->Write();
      
  // Save various derived information
  TH1F *hCount=new TH1F("hCount","Dark Pulse Count 0.5PE;;",1,-1,1);
  hCount->SetBinContent(1,totPeaks05);
  hCount->Write();
  
  double onePEadc = peFcn->GetParameter(1);
  double onePEmV = onePEadc*adc2mV;
  double onePEadcErr = peFcn->GetParError(1);
  double onePEadcSigma = peFcn->GetParameter(2);
  double onePEadcSigmaErr = peFcn->GetParError(2);
  double onePEsigmaOmean = onePEadcSigma/onePEadc;

   
  //Save mean of 1PE peak
  TH1F *h1PePeak = new TH1F("h1PePeak","1 PE Peak Mean Value",1,-1,1);
  h1PePeak->SetBinContent(1,onePEadc);
  h1PePeak->SetBinError(1,onePEadcErr);
  h1PePeak->Write();
  TH1F *h1PePeakmV = new TH1F("h1PePeakmV","1 PE Peak Mean Value, mV",1,-1,1);
  h1PePeakmV->SetBinContent(1,onePEmV);
  h1PePeakmV->Write();
  
  //Save sigma of 1PE peak
  TH1F *h1PePeakSigma = new TH1F("h1PePeakSigma","1 PE Peak Sigma Value",1,-1,1);
  h1PePeakSigma->SetBinContent(1,onePEadcSigma);
  h1PePeakSigma->SetBinError(1,onePEadcSigmaErr);
  h1PePeakSigma->Write();

  //Save sigma/1PE peak
  TH1F *hSigmaOMean = new TH1F("hSigmaOMean","Sigma/Mean",1,-1,1);
  hSigmaOMean->SetBinContent(1,onePEsigmaOmean);
  hSigmaOMean->Write();
  
 
  std::cout << "===============================" << std::endl;
  std::cout << "Number of buffers: "<< nbuf << " total time: " << timeTotal << std::endl;
  std::cout << "Number of peaks (minimal cut): " << totPeaks << std::endl;
  std::cout << "Number of peaks >0.1PE " << totPeaks01 << std::endl;
  std::cout << "Number of peaks >0.5PE " << totPeaks05 << std::endl;
  std::cout << "Dark pulse rate (counting): " << dcrCount/1e6 << " MHz" << std::endl;
  std::cout << "Dark pulse rate (fit01): " << dcrFit01/1e6 << " MHz" << std::endl;
  std::cout << "Afterpulse probability (fit01):  " << rateAP01 << std::endl;
  std::cout << "Dark pulse rate (fit05): " << dcrFit05/1e6 << " MHz" << std::endl;
  std::cout << "Afterpulse probability (fit05):  " << rateAP05 << std::endl;
  std::cout << "Afterpulse probability (2D):  " << apCalc2D(hdPT,1e-6/dcrCount,mu,sig) << std::endl;
  std::cout << "Crosstalk Fraction (1.5 PE threshold): " << xTalkFrac << std::endl;
  std::cout << "1Pe Peak value: " << onePEadc << " (ADC)  " << onePEmV << " (mV)" << std::endl;
  std::cout << "===============================" << std::endl;

  f.Close();


  
  std::cout<< "Close TCanvas: Peaks and Time distributions to exit" << std::endl;
  tc1->Connect("TCanvas","Closed()","TApplication",gApplication,"Terminate()");
  if (quit)
    return 0;
  else
    theApp.Run(true);
  
  return 0;

}

