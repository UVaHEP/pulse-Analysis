import os 


env = Environment()
include = Dir('#/include')
picoscopeDriver = Dir('/opt/uvahep/picoscopeInterface/include')
env = Environment(LIBS=['ps6000'], CCFLAGS='-std=c++11 -Wl,rpath=./', CPPPATH=[include], LIBPATH=['/opt/picoscope/lib', './lib'])

Export('env')

SConscript(['src/SConscript'])

