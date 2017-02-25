import matplotlib.pyplot as plt
from ROOT import gROOT, TCanvas, TFile
import numpy as n
import time

print("plot initialized")

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

def initGraphs(_file):
    f = TFile(_file)
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

    c.SaveAs("/export/home/wwyatt/Downloads/%s.png"%_file)    

def init(files):
    for f in files:
        initGraphs(f)

init(hitmiss_files)
