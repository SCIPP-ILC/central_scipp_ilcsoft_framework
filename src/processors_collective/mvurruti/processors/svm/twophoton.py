from ROOT import *
f = TFile("twophoton.root")
c = TCanvas()
c.SetLogy()

def plot(name): 
    graph = f.Get(name)
    graph.Draw()
    c.SaveAs("./%s.png"%name)

graphs = ["S","V","M"]
for g in graphs:
    plot(g)
