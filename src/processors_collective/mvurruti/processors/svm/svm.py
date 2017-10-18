import time
import argparse
import ROOT

parser = argparse.ArgumentParser(description = "plots SVM root files.")

parser.add_argument("filename", metavar = "filename", nargs = "?", help = "input a root file name to plot.")

args = parser.parse_args()
prefix = args.filename
if prefix is None:
    raise Exception("Dude, you gotta supply a filename.")

if ".root" in prefix:
    prefix = prefix[-5:]
print("%s.root" % prefix)
f = ROOT.TFile("%s.root"%prefix)
c = ROOT.TCanvas()
c.SetLogy()

def plot(name): 
    graph = f.Get(name)
    try:
        graph.Draw()
    except ReferenceError:
        print("Tried to find a graph (%s) in the root file that was not there."%name)
        print("This is what is in the root file:")
        f.ls()
        exit()
    graph.GetYaxis().SetTitle("# of Events")
    if name == "V_Dtd" or name == "V_Dbl" or name == "V_Tru" or name == "S_Dtd" or name =="S_Dbl" or name == "S_Tru":
        graph.GetXaxis().SetTitle("Momentum (GeV)")
    if name == "M_Dtd" or name == "M_Dbl" or name == "M_Tru":
        graph.GetXaxis().SetTitle("Mass (GeV)")    
    time.sleep(4)
    c.SaveAs("./%s_%s.png"%(prefix,name))

graphs = ["V_Dtd","V_Dbl","V_Tru","S_Dtd","S_Dbl","S_Tru","M_Dtd","M_Dbl","M_Tru"]
for g in graphs:
    plot(g)
