from ROOT import TFile, TH1, TCanvas
from rutil import TCachingFile, Frame, canvas, nostat, keep, parseOpt


collection = "offlinePrimaryVertices4D"

tf = []
for rootfile, marker in (
        ("pvt_RUN4NOPUTTBAR_new_A_000.root", "opencircle"),
        ("pvt_RUN4PU200TTBAR_new_A_000.root","fullcircle"),
        ("pvt_RUN4PU200ZMM_new_A_000.root", "fulltriangle")):
    tcf = TCachingFile(rootfile) 
    tcf.cd(collection)
    tcf.marker = marker
    tf.append(tcf)

    
c1 = canvas("landscape")
F = F1 = Frame("time residuals, Signal vertex", "vertex trec - tsim [ns]",-0.1, 0.1)
for f in tf:
    F.add(f.Get("matchedvtxsel/trecsimSignal"), f.marker)
F.drawLegend()
F.Print("Signal-vertex-time-residual.png")


c2 = canvas("landscape")
F = F2 = Frame("time residuals, PU", "vertex trec - tsim [ns]",-0.1, 0.1)
for f in tf:
    F.add(f.Get("matchedvtxsel/trecsimPU"), f.marker)
F.drawLegend()
F.Print("PU-vertex-time-residual.png")

    
c3 = canvas("landscape")
F = F3 = Frame("normalized time residual, Signal vertex", "vertex (trec - tsim) / error",-10., 10.)
for f in tf:
    F.add(f.Get("matchedvtxsel/trecsimpullSignal"), f.marker)
F.drawLegend()
F.Print("Signal-vertex-time-pull.png")


c4 = canvas("landscape")
F = F4 = Frame("normalized time residuals, PU", "vertex (trec - tsim) / error",-10., 10.)
for f in tf:
    F.add(f.Get("matchedvtxsel/trecsimpullPU"), f.marker)
F.drawLegend()
F.Print("PU-vertex-time-pull.png")


c5 = canvas("landscape")
F = F4 = Frame("vertex time error vs number of tracks", "tracks with timing", 0. ,200, "error (ns)", 0.,  0.04)
for f in tf:
    F.add(f.Get("matchedvtxsel/trecerrvsntrk"), f.marker)
F.drawLegend()
F.Print("vertex-time-errror-vs-ntrk.png")
raw_input()
    

