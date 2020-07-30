from ROOT import TFile, TH1, TCanvas
from rutil import TCachingFile, Frame, canvas, nostat, keep, parseOpt

args = parseOpt(dir="v21")

def show(key):
    return ("all" in args.args) or (str(key) in args.args)

def addplots(F, tf, hid):
    for f in tf:
        h = f.Get(hid)
        h.Scale(1./h.Integral())
        F.add(h, f.marker)
        keep(h)
    F.drawLegend()
 
tf = []
for rootfile, marker in (
        ("pvt_RUN4NOPUTTBAR_new_A.root", "opencircle"),
        ("pvt_RUN4PU200TTBAR_new_A.root","fullcircle"),
        #("pvt_RUN4PU200ZMM_new_A.root", "fulltriangle")
     ):
    tcf = TCachingFile(args.dir + "/" + rootfile) 
    tcf.cd("tracks")
    tcf.marker = marker
    tf.append(tcf)


if show(1):
    c1 = canvas("landscape")
    F = F1 = Frame("matched track time residuals", "track trec - tsim [ns]",-0.5, 0.5)
    for f in tf:
        h = f.Get("trestrk_selmatched")
        h.Scale(1./h.Integral())
        h.Fit("gaus", "N", "", -0.05, 0.05)
        F.add(h, f.marker)
        keep(h)
    F.drawLegend()
    F.Print("track-time-residuals.png")


if show(2):
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


if show("2a"):
    c2a = canvas("landscape")
    F = F2a = Frame("normalized time residuals (#sigma_{t}<0.1)", "track (trec - tsim) / error",-10., 10.)
    for f in tf:
        h = f.Get("tpulltrk_sigmat01_selmatched")
        h.Scale(1./h.Integral())
        #h.Fit("gaus", "", "same")
        F.add(h, f.marker)
        keep(h)
        F.drawLegend()
        F.Print("track-time-pulls-small-sigmat.png")

if show(3):
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

if show("3a"):
    c3a = canvas("landscape")
    F = F3a = Frame("normalized time residuals for pions (#sigma_{t}<0.1)", "track (trec - tsim) / error",-10., 10.)
    for f in tf:
        h = f.Get("tpulltrk_sigmat01_selpion")
        h.Scale(1./h.Integral())
        #h.Fit("gaus", "", "same")
        F.add(h, f.marker)
        keep(h)
    F.drawLegend()
    F.Print("pion-track-time-pulls-small-sigmat.png")

if show("4"):
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

if show("4a"):
    c4a = canvas("landscape")
    F = F4a = Frame("normalized time residuals for kaons (#sigma_{t}<0.1)", "track (trec - tsim) / error",-10., 10.)
    for f in tf:
        h = f.Get("tpulltrk_sigmat01_selkaon")
        h.Scale(1./h.Integral())
        #h.Fit("gaus", "", "same")
        F.add(h, f.marker)
        keep(h)
        F.drawLegend()
        F.Print("kaon-track-time-pulls-small-sigmat.png")

if show(5):
    c5 = canvas("landscape")
    F = F4 = Frame("track time errors", "track error (ns)",0., 0.1, logy=True, ymax=0.4)
    for f in tf:
        h = f.Get("terrtrk_rec_selmatched")
        h.Scale(1./h.Integral())
        F.add(h, f.marker)
        keep(h)
        F.drawLegend()
        F.Print("track-time-errors.png")


if show(6):
    c = c6 = canvas("landscape", "", 2, 2)
    c.cd(1)
    F61 = Frame("","",-1.1, 1.1)
    addplots(F61, tf, "tqualtrk_rec_selmatched_barrel_hipt")

    c.cd(2)
    F62 = Frame("","",-1.1, 1.1)
    addplots(F62, tf, "tqualtrk_rec_selmatched_endcap_hipt")

    c.cd(3)
    F61 = Frame("","",-1.1, 1.1)
    addplots(F63, tf, "tqualtrk_rec_selmatched_barrel_lopt")

    c.cd(4)
    F62 = Frame("","",-1.1, 1.1)
    addplots(F64, tf, "tqualtrk_rec_selmatched_endcap_lopt")
    
    c.Print("yo.png")

raw_input()
