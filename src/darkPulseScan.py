#!/usr/bin/python

import argparse,sys
import subprocess,time

# keep ROOT TApplication from grabbing -h flag
from ROOT import PyConfig
PyConfig.IgnoreCommandLineOptions = True
from ROOT import *


class bcolors:
    HEADER = '\033[95m' # purple
    INFO = '\033[94m'  #blue
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#######################
# main
#######################

parser = argparse.ArgumentParser(description='Take dark pulse data')
parser.add_argument('-v','--voltage', type=float, default = 0,
                    help="starting voltage")
parser.add_argument('-s','--nsteps', type=int, default=0,
                    help="number of voltage steps")
parser.add_argument('-S','--stepsize', type=float, default = 0.0,
                    help="size of voltage steps")
parser.add_argument('-x', '--vmax', type=float, default = 0,
                    help="max voltage")
parser.add_argument('-R','--range', type=str, default="PS_20MV",
                    help="set voltage scale on picoscope [PS_20MV]")
parser.add_argument('-z', '--zero', action='store_true',
                    help="zero voltage measure")
parser.add_argument('-o', '--output', type=str, default = "darkBuffers",
                    help="outputfile template")
parser.add_argument('-b', '--nbuf', type=int, default = 50,
                    help="number of buffers to take [5]")
parser.add_argument('-f', '--usefile', dest="usefile", default=None, action="store_true",
                    help="Use this flag if reading input from root files")
parser.add_argument('-0', '--quiet', dest="quiet", default=None, action="store_true",
                    help="Use this flag to disable graphics")
args = parser.parse_args()


if args.voltage: voltage=args.voltage
if not args.voltage and not args.zero:
    print "No voltage given:",args.voltage
    sys.exit()
if args.zero: voltage=0
#if  args.voltage<=0:
#    print "Invalid voltage given:",args.voltage
#    sys.exit()
vmax=args.vmax
outname=args.output

nsteps = args.nsteps
stepsize = args.stepsize


if nsteps>0:
    if vmax>0: stepsize=float(vmax-voltage)/nsteps
    elif stepsize==0:
        print "Can't determine voltage steps"
        sys.exit()
elif stepsize>0 and vmax>0:
    nsteps=int(round(float(vmax-voltage)/stepsize))
    

#loop over voltages
vend=vmax

if args.usefile: tfScan=TFile(outname+"_.root","recreate")
else: tfScan=TFile(outname+".root","recreate")
tgDCR=TGraphErrors(); tgDCR.SetTitle("Dark count rate vs Voltage;V;MHz"); tgDCR.SetName("gDCR")
tgDCR.SetMaximum(10)
tgAP=TGraphErrors(); tgAP.SetTitle("Afterpulse probability vs Voltage"); tgAP.SetName("gAPRate")
tgXT=TGraphErrors(); tgXT.SetTitle("Crosstalk fraction vs. Voltage;V"); tgXT.SetName("gXtalkFrac")
tg1P=TGraphErrors(); tg1P.SetTitle("1PE Peak vs. Voltage;V"); tg1P.SetName("gOnePE")
tgSMM=TGraphErrors(); tgSMM.SetTitle("Sigma/Mean vs. Mean"); tgSMM.SetName("gSoMvM")
tgSMV=TGraphErrors(); tgSMV.SetTitle("Sigma/mean vs. Voltage;V"); tgSMV.SetName("gSoMvV")
tIV=TGraph(); tIV.SetTitle("I vs V"); tIV.SetName("gIV")
tc=TCanvas("cgr","Pulse Data",1500,1000)
tc.Divide(3,2)

#print "nsteps",nsteps
nbuf=args.nbuf

vsteps=[]
for i in range(nsteps+1):
    vsteps.append(round(voltage+i*stepsize,3))

for i in range(nsteps+1):
    v=vsteps[i]
    print bcolors.HEADER+"\nStarting data at V= "+str(v)+bcolors.ENDC
    
    # take pulses
    filename=outname+"_"+str(v)+".root"
    flags=""
    if args.quiet: flags="-0"
    if args.usefile:
        filename2=filename.replace(".root","_.root")
        subprocess.call(["darkBuffers",flags,"-qf"+filename, "-o"+filename2])
        resultsFile=filename2
        iVal=0  # no readback voltage
    else:
        #set voltage
        subprocess.call(["setVoltage.py","-pqv"+str(v)])
        iReading=subprocess.check_output(["setVoltage.py","-pqv"+str(v)])
        iVal=iReading.split("measure:")[1]
        iVal=float(iVal.split()[0])
        print "Readback current",iVal    
        subprocess.call(["darkBuffers",flags,"-aqb"+str(nbuf), "-o"+filename])
        resultsFile=filename

    print "Reading",resultsFile
    tf=TFile(resultsFile)
    version=tf.Get("hVersion").GetBinContent(1);
    darkRate=tf.Get("hRate").GetBinContent(1);
    error=tf.Get("hRate").GetBinError(1);
    tgDCR.SetPoint(tgDCR.GetN(),v,darkRate/1e6);
    tgDCR.SetPointError(tgDCR.GetN()-1,0,error/1e6);
    if version==1:
        apRate=tf.Get("hAp").GetBinContent(1);
        apErr=tf.Get("hAp").GetBinError(1);
    else:
        apRate=tf.Get("hAp").GetBinContent(3);
        apErr=tf.Get("hAp").GetBinError(3);
    error=tf.Get("hAp").GetBinError(1);
    tgAP.SetPoint(tgAP.GetN(),v,apRate);
    tgAP.SetPointError(tgAP.GetN()-1,0,apErr);
    crossTalk=tf.Get("hCrossTalk").GetBinContent(1);
    error=tf.Get("hCrossTalk").GetBinError(1);
    tgXT.SetPoint(tgXT.GetN(),v,crossTalk);
    tgXT.SetPointError(tgXT.GetN()-1,0,error);
    pePeak=tf.Get("h1PePeak").GetBinContent(1);
    error=tf.Get("h1PePeak").GetBinError(1);
    tg1P.SetPoint(tg1P.GetN(),v,pePeak);
    tg1P.SetPointError(tg1P.GetN()-1,0,error);
    sigmaOmean=tf.Get("hSigmaOMean").GetBinContent(1);
    tgSMM.SetPoint(tgSMM.GetN(),pePeak,sigmaOmean);
    tgSMV.SetPoint(tgSMV.GetN(),v,sigmaOmean);
    tIV.SetPoint(tIV.GetN(),v,iVal);
    tc.cd(1);
    tgDCR.Draw("ALP*")
    tc.cd(2)
    tgAP.Draw("ALP*")
    tc.cd(3);
    tgXT.Draw("ALP*")
    tc.cd(4);
    tg1P.Draw("ALP*")
    tc.cd(5);
    tgSMM.Draw("ALP*")
    tc.cd(6);
    tgSMV.Draw("ALP*")
    tc.Update()
    tf.Close()
    
if args.usefile==None: subprocess.call(["setVoltage.py"])

#time.sleep(2)
tc.SaveAs(outname+".pdf")
tfScan.cd()
tgDCR.Write()
tgAP.Write()
tgXT.Write()
tg1P.Write()
tgSMM.Write()
tgSMV.Write()
tIV.Write()
tfScan.Write()
tfScan.Close()

