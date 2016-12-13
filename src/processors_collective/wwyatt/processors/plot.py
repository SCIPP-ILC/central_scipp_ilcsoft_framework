import matplotlib.pyplot as plt
from ROOT import gROOT, TCanvas, TFile
import numpy as n
f = TFile("Bhabha_aa.root")
names = [
    "phi_s",
    "phi_m",
    "phi_l",
]
graphs = {}
for name in names:
    graphs[name] = f.Get(name)
    graphs[name].Draw()
