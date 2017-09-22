#ifndef PICOUTILS_H
#define PICOUTILS_H

#include "ps5000a.h"
#include "TApplication.h"
#include <vector>

using std::vector;
using namespace picoscope;

void setupScope(ps5000a &dev, chRange &range, int samples, unsigned int timebase=1);
int autoRange(ps5000a &dev, int nbuf=1);
void acquireBuffers(ps5000a &dev, vector <vector<short> > &data);
Double_t userThresholdFn(ps5000a &dev, int samples, TApplication &app, bool invertFlag=false);
int readBuffers(TString filename, vector <vector<short> > &data, double &timebase, double &adc2mV);

#endif


