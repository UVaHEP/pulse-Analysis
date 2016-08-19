#!/usr/bin/python

import argparse,sys
import subprocess,time


class bcolors:
    HEADER = '\033[95m' # purple
    INFO = '\033[94m'  #blue
    OKGREEN = '\033[92m'
    WARNING = '\033[93m' #Yellow
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#######################
# main
#######################

parser = argparse.ArgumentParser(description='Take light pulse data')

parser.add_argument('-v','--voltage', type=float, default = 0,
                    help="voltage applied for pulses")
parser.add_argument('-0', '--zero', action='store_true',
                    help="zero voltage measure")
parser.add_argument('-o', '--output', type=str, default = "lightPulseScans_defaultFileName",
                    help="outputfile template")
parser.add_argument('-b', '--nbuf', type=int, default = 10000,
                    help="number of buffers to take [10000]")
parser.add_argument('-I','--IntensitySteps', type=int, default = 0,
                    help="one less than the number of different intensities you plan to take")
parser.add_argument('-p','--pePeak', type=float, default = 0,
                    help="if you manually want to input a 1 pe peak")
#the way to get output in ADC currently cannot be called from here.
#Everything is labeled in mV anyway and I didn't feel like changing it
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
outname=args.output

tg=TGraphErrors()
tg.SetTitle("Sigma/Mean vs. Mean;mean (n*1PePeak)")
tg2=TGraphErrors()
tg2.SetTitle("Mean (1pe normalization) vs. Intensity; Laser Intensity")
tc=TCanvas("tc","Peak height Distribution data",900,450)
tc.Divide(2,1)

nbuf=args.nbuf
inten_steps = args.IntensitySteps
pe_peak = args.pePeak

#set voltage
subprocess.call(["setVoltage.py","-pqv"+str(voltage)])

#Find 1 pe peak
if pe_peak != 0:
    print "You opted not to take dark buffers. Your 1 pe peak value is"+str(pe_peak)
    onePePeak = pe_peak
else:
    print bcolors.HEADER+"Taking dark pulses to find 1 Pe peak"+bcolors.ENDC
    darkBuffname=outname+"_1peNormalization_darkCalculation.root"
    subprocess.call(["./darkBuffers","-b50","-o"+darkBuffname])
    dB=TFile(darkBuffname)
    onePePeak=dB.Get("h1PePeakmV").GetBinContent(1);

for i in range(inten_steps+1):
    print "The normalization to 1 pe peak is "+str(onePePeak)
    print bcolors.HEADER+"\nAbout to take pulses. Please go to the laser and set an intensity"+bcolors.ENDC
    print bcolors.HEADER+"\nPlease record the intensity (user input): "+bcolors.ENDC
    inten_level = raw_input()
    print bcolors.HEADER+"\nYour intensity is "+str(inten_level)+bcolors.ENDC
    filename=outname+"_inten"+str(inten_level)+"_"+str(voltage)+"V.root"
    #take pulses
    subprocess.call(["./lightBuffersMorePeaks","-b"+str(nbuf), "-o"+filename, "-p", "-n"+str(onePePeak)])
    # Create the graph of sigma/mean
    tf=TFile(filename)
    sigma_over_mean = tf.Get("hSigmaOMeanPe").GetBinContent(1);
    mean_penorm_height = tf.Get("hMeanPePeaks").GetBinContent(1);
    tg.SetPoint(tg.GetN(),mean_penorm_height,sigma_over_mean);
    tg2.SetPoint(tg2.GetN(),float(inten_level),mean_penorm_height);
    tc.cd(1);
    tg.Draw("ALP*")
    tc.cd(2);
    tg2.Draw("ALP*")
    tc.Update()
    print bcolors.WARNING+"When done, turn off the laser... peasant!"+bcolors.ENDC
    print "Your last used intensity was "+str(inten_level)

subprocess.call(["setVoltage.py"])

time.sleep(2)
tc.SaveAs(outname+".pdf")


