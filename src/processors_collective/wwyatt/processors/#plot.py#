import matplotlib.pyplot as plt
from ROOT import gROOT, TCanvas, TFile
import numpy as n
import time

hitmiss_files = [
"aa_positron_hitmiss.root",
"aa_electron_hitmiss.root",
"a_positron_hitmiss.root",
"a_electron_hitmiss.root",
"z_electron_hitmiss.root",
"z_positron_hitmiss.root",
]
def plotHitMiss(g, color):
    g.SetLineColor(color)
    g.SetMarkerColor(color)

def initGraphs():
    f = TFile("aa_positron_hitmiss.root")
    c = TCanvas()
    names = [
        "hh",
        "mm",
        "hm",
        ]
    graphs = {}
    color = 2

    for name in names:
        graph = f.Get(name)
        graph.Draw("SAME")
        plotHitMiss(graph, color)
        color+=1

    c.SaveAs("/export/home/wwyatt/Downloads/%s.png"%"marsietn")    

