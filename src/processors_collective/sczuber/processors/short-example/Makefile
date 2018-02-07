CPPFLAGS = -Wall -g --std=c++11 \
		   -I/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/FastJet/3.1.2/include \
		   -I/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/CLHEP/2.1.4.1/include \
		   -I/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/lcio/v02-07/include \
		   -Wl,-rpath,/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/FastJet/3.1.2/lib \
		   -Wl,-rpath,/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/CLHEP/2.1.4.1/lib \
		   $(shell root-config --cflags)

LDFLAGS = -lm \
			-L/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/CLHEP/2.1.4.1/include \
			-L/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/FastJet/3.1.2/lib \
			-lfastjettools \
			-lfastjet \
			-lfastjetplugins \
			-L/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/lcio/v02-07/lib \
			-llcio \
			-lsio \
			-lsiscone_spherical \
			-lsiscone \
			-L/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-09/root/5.34.30/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -rdynamic


EXECUTABLES = validateMCParticles fastJetMCParticles fastJetReconstructedParticles MCParticleROOT testEvent


short-example-2.x :
	g++ $(CPPFLAGS) $(LDFLAGS) short-example-2.cc -o short-example-2.x

short-example.x :
	g++ $(CPPFLAGS) $(LDFLAGS) short-example.cc -o short-example.x


#short-example.o: short-example.cc short-example.hh
#	g++ $(CPPFLAGS) -c short-example.cc

.PHONY: clean

clean : 
	rm *.o *.x

#MyDict.cxx: Event.hh EventLinkdef.h
#	rootcint -f $@ -c -p $^


#libEvent.so: MyDict.cxx Event.cc
#	g++ -g -shared -o$@ `root-config --ldflags` -fPIC -O2 -I$(ROOTSYS)/include $^


#$(EXECUTABLES): %: short-example.o %.cc libEvent.so
#	gcc $*.cc $(CPPFLAGS) $(LDFLAGS) short-example.o libEvent.so -o $@
