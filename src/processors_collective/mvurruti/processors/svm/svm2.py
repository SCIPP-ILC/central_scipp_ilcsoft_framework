import time
import argparse
import ROOT

#list of root files as input
#loop through input files and plot based of the types (S,V,M)
#Start w/ vector true alone.
#plot all true vectors of diff root files in same plot

roots = [
    "bbsignal.root","bwsignal.root","signal_all.root"
    ]

c = ROOT.TCanvas()

for file_name in roots:
    print file_name
    file = ROOT.TFile(file_name)
    name = "V_Tru"
    graph = file.Get(name)
    
    try:
        graph.Draw("SAME")
    except ReferenceError:
        print("Tried to find a graph (%s) in the root file that was not there."%name)
        print("This is what is in the root file:")
        file.ls()
        exit()

    
time.sleep(4)
