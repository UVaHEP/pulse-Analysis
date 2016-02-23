#!/usr/bin/python

import argparse,sys
import subprocess,time


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
parser.add_argument('-0', '--zero', action='store_true',
                    help="zero voltage measure")
parser.add_argument('-o', '--output', type=str, default = "darkBuffers",
                    help="outputfile template")
parser.add_argument('-b', '--nbuf', type=int, default = 5,
                    help="number of buffers to take [5]")
args = parser.parse_args()

from ROOT import *

if args.voltage: voltage=args.voltage
if not args.voltage and not args.zero:
    print "No voltage given:",args.voltage
    sys.exit()
if args.zero: voltage=0
if  args.voltage<=0:
    print "Invalid voltage given:",args.voltage
    sys.exit()
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

tg=TGraphErrors()
tg.SetTitle("Dark pulse rate vs Voltage;V;MHz")
tga=TGraphErrors()
tga.SetTitle("Afterpulse rate vs Voltage")
tc=TCanvas("cgr","Pulse Data",1000,500)
tc.Divide(2,1)

#print "nsteps",nsteps
nbuf=args.nbuf

for i in range(nsteps+1):
    v=round(voltage+i*stepsize,3)
    print bcolors.HEADER+"\nStarting data at V= "+str(v)+bcolors.ENDC
    # set voltage
    subprocess.call(["setVoltage.py","-pqv"+str(v)])
    # takepulses
    filename=outname+"_"+str(v)+".root"
    print "Saving data to",filename
    subprocess.call(["darkBuffers","-qb"+str(nbuf), "-o"+filename])
    tf=TFile(filename)
    darkRate=tf.Get("hRate").GetBinContent(1);
    error=tf.Get("hRate").GetBinError(1);
    tg.SetPoint(tg.GetN(),v,darkRate);
    tg.SetPointError(tg.GetN()-1,0,error);
    afterRate=tf.Get("hAp").GetBinContent(1);
    error=tf.Get("hAp").GetBinError(1);
    tga.SetPoint(tga.GetN(),v,afterRate);
    tga.SetPointError(tga.GetN()-1,0,error);
    tc.cd(1);
    tg.Draw("ALP*")
    tc.cd(2)
    tga.Draw("ALP*")
    tc.Update()
    
subprocess.call(["setVoltage.py"])

time.sleep(2)
tc.SaveAs(outname+".pdf")


