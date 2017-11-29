import time
import argparse
import ROOT

#list of root files as input
#loop through input files and plot based of the types (S,V,M)
#Start w/ vector true alone.
#plot all true vectors of diff root files in same plot

roots = ["bbsignal.root","bwsignal.root","wbsignal.root","wwsignal.root"]

files=[]
for file_name in roots:
    files.append(ROOT.TFile(file_name))

names = ["V_Tru","V_Dtd","V_Dbl","S_Dtd","S_Dbl","S_Tru","M_Dtd","M_Dbl","M_Tru"]

for name in names:
    c = ROOT.TCanvas()
    legend = ROOT.TLegend(.75,.75,1,1)
    
    c.SetLogy()
    i = 2
    for f in files:
        f.cd()
        graph = f.Get(name)
    
        try:
            graph.SetLineColor(i)
            i+=1
            graph.Draw("same")
            legend.AddEntry(name,"test","l")
        except ReferenceError:
            print("Tried to find a graph (%s) in the root file that was not there."%name)
            print("This is what is in the root file:")
            file.ls()
            exit()
        except AttributeError:
            print("Misspelled root file")
            exit()
    legend.Draw()
    c.SaveAs("./%s.png"%(name))   

