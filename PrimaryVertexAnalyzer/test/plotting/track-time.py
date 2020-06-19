from ROOT import TFile, TH1, TCanvas
from rutil import TCachingFile, Frame, canvas, nostat, keep



tf = []
for rootfile, marker in (
        ("pvt_RUN4NOPUTTBAR_new_A_000.root", "opencircle"),
        ("pvt_RUN4PU200TTBAR_new_A_000.root","fullcircle"),
        ("pvt_RUN4PU200ZMM_new_A_000.root", "fulltriangle")):
    tcf = TCachingFile(rootfile) 
    tcf.cd("tracks")
    tcf.marker = marker
    tf.append(tcf)

    
c1 = canvas("landscape")
F = F1 = Frame("matched track time residuals", "track trec - tsim [ns]",-0.5, 0.5)
for f in tf:
    h = f.Get("trestrk_selmatched")
    h.Scale(1./h.Integral())
    #h.Fit("gaus", "same")
    F.add(h, f.marker)
    keep(h)
F.drawLegend()
F.Print("track-time-residuals.png")


c2 = canvas("landscape")
F = F2 = Frame("normalized time residuals", "track (trec - tsim) / error",-10., 10.)
for f in tf:
    h = f.Get("tpulltrk_selmatched")
    h.Scale(1./h.Integral())
    #h.Fit("gaus", "", "same")
    F.add(h, f.marker)
    keep(h)
F.drawLegend()
F.Print("track-time-pulls.png")


c3 = canvas("landscape")
F = F3 = Frame("normalized time residuals for pions", "track (trec - tsim) / error",-10., 10.)
for f in tf:
    h = f.Get("tpulltrk_selpion")
    h.Scale(1./h.Integral())
    #h.Fit("gaus", "", "same")
    F.add(h, f.marker)
    keep(h)
F.drawLegend()
F.Print("pion-track-time-pulls.png")

c4 = canvas("landscape")
F = F4 = Frame("normalized time residuals for kaons", "track (trec - tsim) / error",-10., 10.)
for f in tf:
    h = f.Get("tpulltrk_selkaon")
    h.Scale(1./h.Integral())
    #h.Fit("gaus", "", "same")
    F.add(h, f.marker)
    keep(h)
F.drawLegend()
F.Print("kaon-track-time-pulls.png")

c5 = canvas("landscape")
F = F4 = Frame("track time errors", "track error (ns)",0., 0.1)
for f in tf:
    h = f.Get("terrtrk_rec_selmatched")
    h.Scale(1./h.Integral())
    F.add(h, f.marker)
    keep(h)
F.drawLegend()
F.Print("track-time-errors.png")


raw_input()
