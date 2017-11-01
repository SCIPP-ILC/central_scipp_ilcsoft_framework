import os
import ROOT
from ROOT import TFile, TCanvas, TH1D, gStyle, TPaveText

Outpath = os.environ["ILC"]+"/output/"

#Masses_eLpR = [ '39113', '39117', '39121', '39125', '39129', '39133' ]
Polarizations = [ 'eB.pB.I39215', 'eB.pW.I39214', 'eW.pB.I39213', 'eW.pW.I39212' ]


Plots = [ 
        #  'T']  #with neutrinos, all angles
        #  'DAB']  #without neutrinos, all angles
         'DED' ] #without neutrinos, forward angles
Pol = {}
Pol['eB.pB.I39215'] = 'BB'
Pol['eB.pW.I39214'] = 'BW'
Pol['eW.pB.I39213'] = 'WB'
Pol['eW.pW.I39212'] = 'WW'

Plotnames = {}
Plotnames['T'] = 'True'
Plotnames['DAB'] = 'Detectable'
Plotnames['DED'] = 'Detected'

ptype = 'TV'

def plot(Polarization):
    Rootfiles = []
    for pol in Polarization:
        #print pol
        for plotstring in Plots:
            print plotstring
            Rootfiles.append( TFile(Outpath+"TwoPhotonThrust_"+pol+"._"+plotstring+".root") )

    for plotstring in Plots:
        print plotstring
        canvas1 = TCanvas("c1","c1")
        firsthist = True
        for i,rootfile in enumerate(Rootfiles):
            #print i
             
            exec('plot = rootfile.'+ptype+'_'+plotstring) 
            if i+1>= 5: plot.SetLineColor(i+2)
            else : plot.SetLineColor(i+1)
            plot.SetLineWidth(2)
            plot.SetLineStyle(1)
            #if ptype == 'S': plot.SetLineStyle(1)
            #else: plot.SetLineStyle(2)
            plot.SetTitle(Pol[Polarization[i]])
            plot.Draw("same")
            #plot.SetLogy()

            if firsthist:
                plot.SetStats(False)
                plot.SetXTitle("Thrust")
                plot.GetXaxis().CenterTitle()
                plot.SetYTitle("Events")
                plot.GetYaxis().CenterTitle()
                plot.GetXaxis().SetRangeUser(-3,4)
                plot.GetYaxis().SetRangeUser(0,10000)
                
                firsthist = False
        #STOP rootfiles loop
        gStyle.SetOptTitle(0); #turn off default title
        title = TPaveText(0.25, 0.9, 0.75, 1.0, "brNDC"); 
        title.AddText(Plotnames[plotstring] + ' Thrust Value');
        title.Draw("same");

        canvas1.SetGrid(1,0)
        canvas1.BuildLegend() 
        canvas1.SaveAs(Outpath+ptype+'_'+plotstring+'_'+'aalowpt'+".png")
        canvas1.SaveAs(Outpath+ptype+'_'+plotstring+'_'+'aalowpt'+".ps")
#END PLOT



plot(Polarizations)



