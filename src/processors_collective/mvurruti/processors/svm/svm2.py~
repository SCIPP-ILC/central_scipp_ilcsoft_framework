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
f = ROOT.TFile("%s.root"%prefix)
c = ROOT.TCanvas()
c.SetLogy()

def plot(name): 
    graph = f.Get(name)
    graph.GetYaxis().SetTitle("# of Events")
    if name == "S" or name == "V":
        graph.GetXaxis().SetTitle("Momentum (GeV)")
    if name == "M":
        graph.GetXaxis().SetTitle("Mass (GeV)")    
    graph.Draw()
    time.sleep(4)
    c.SaveAs("./%s_%s.png"%(prefix,name))

graphs = ["S","V","M"]
for g in graphs:
    plot(g)
