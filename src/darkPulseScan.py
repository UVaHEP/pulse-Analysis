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
tg=TGraphErrors(); tg.SetTitle("Dark count rate vs Voltage;V;MHz"); tg.SetName("gDCR")
tg.SetMaximum(10)
tga=TGraphErrors(); tga.SetTitle("Afterpulse rate vs Voltage"); tga.SetName("gAPRate")
tgb=TGraphErrors(); tgb.SetTitle("Crosstalk fraction vs. Voltage;V"); tgb.SetName("gXtalkFrac")
tgc=TGraphErrors(); tgc.SetTitle("1PE Peak vs. Voltage;V"); tgc.SetName("gOnePE")
tgd=TGraphErrors(); tgd.SetTitle("Sigma/Mean vs. Mean"); tgd.SetName("gSoM1")
tge=TGraphErrors(); tge.SetTitle("Sigma/mean vs. Voltage;V"); tge.SetName("gSom2")
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
    darkRate=tf.Get("hRate").GetBinContent(1);
    error=tf.Get("hRate").GetBinError(1);
    tg.SetPoint(tg.GetN(),v,darkRate/1e6);
    tg.SetPointError(tg.GetN()-1,0,error/1e6);
    afterRate=tf.Get("hAp").GetBinContent(1);
    error=tf.Get("hAp").GetBinError(1);
    tga.SetPoint(tga.GetN(),v,afterRate);
    tga.SetPointError(tga.GetN()-1,0,error);
    crossTalk=tf.Get("hCrossTalk").GetBinContent(1);
    error=tf.Get("hCrossTalk").GetBinError(1);
    tgb.SetPoint(tgb.GetN(),v,crossTalk);
    tgb.SetPointError(tgb.GetN()-1,0,error);
    pePeak=tf.Get("h1PePeak").GetBinContent(1);
    error=tf.Get("h1PePeak").GetBinError(1);
    tgc.SetPoint(tgc.GetN(),v,pePeak);
    tgc.SetPointError(tgc.GetN()-1,0,error);
    sigmaOmean=tf.Get("hSigmaOMean").GetBinContent(1);
    tgd.SetPoint(tgd.GetN(),pePeak,sigmaOmean);
    tge.SetPoint(tge.GetN(),v,sigmaOmean);
    tIV.SetPoint(tIV.GetN(),v,iVal);
    tc.cd(1);
    tg.Draw("ALP*")
    tc.cd(2)
    tga.Draw("ALP*")
    tc.cd(3);
    tgb.Draw("ALP*")
    tc.cd(4);
    tgc.Draw("ALP*")
    tc.cd(5);
    tgd.Draw("ALP*")
    tc.cd(6);
    tge.Draw("ALP*")
    tc.Update()
    tf.Close()
    
if args.usefile==None: subprocess.call(["setVoltage.py"])

#time.sleep(2)
tc.SaveAs(outname+".pdf")
tfScan.cd()
tg.Write()
tga.Write()
tgb.Write()
tgc.Write()
tgd.Write()
tge.Write()
tIV.Write()
tfScan.Write()
tfScan.Close()

