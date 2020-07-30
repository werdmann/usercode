from rutil import Frame, canvas, TCachingFile, keep, parseOpt
from ROOT import TMath, TF1
import sys



args = parseOpt(pufile="v22/pvt_RUN4PU200TTBAR__A.root", nopufile="v22/pvt_RUN4NOTTBAR__A.root", name="")

#tcf = TCachingFile("v18/pvt_RUN4PU200TTBAR_new_A.root")
tcfpu = TCachingFile(args.pufile)
tcfnopu = TCachingFile(args.nopufile)
scale = 1.0*tcfpu.Get("offlinePrimaryVertices/sigmaZ").GetEntries()/tcfnopu.Get("offlinePrimaryVertices/sigmaZ").GetEntries()
tcfnopu.norm = scale


def get(tcf, hid, style=1,rebin=2, scale=None):
    h = tcf.Get(hid)
    h.SetLineWidth(2)
    h.SetLineStyle(style)
    if rebin > 1:
        h.Rebin(rexbin)
    if scale is not None:
        h.Scale(scale)
    keep(h)
    print hid,h.Integral(), h.GetEntries()
    return h

c1 = canvas("large")
F = F1 = Frame("4D vertices time residuals", "t_{rec}-t_{sim} [ns]", -0.1, 0.1, "# reconstructed vertices")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/trecsimSignal",1,1), "hist", color=2,  label="t#bar{t} 200PU")
F.add(get(tcfnopu, "offlinePrimaryVertices4D/matchedvtxsel/trecsimSignal",2,1), "hist", color=2,  label="t#bar{t} 0 PU")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/trecsimPU",1,1, 0.01), "hist", color=4,  label="PU 4d 200PU")
F.add(get(tcfpu,   "offlinePrimaryVertices/matchedvtxsel/trecsim_fromtracksSignal",1,1), "hist", color=1,  label="3d 200 PU")
F.add(get(tcfnopu, "offlinePrimaryVertices/matchedvtxsel/trecsim_fromtracksSignal",2,1), "hist", color=1,  label="3d 0 PU")
F.drawLegend(0.65, 0.65, 0.95, 0.95)
F1.Print("figures/4d-time-residuals.png")


#c2 = canvas("large")
#F = F2 = Frame("PU", "t_{rec}-t_{sim} [ns]", -0.1, 0.1, "# reconstructed vertices")
#F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/trecsimPU",1,1), "hist", color=2,  label="4d 200PU")
#F.add(get(tcfnopu, "offlinePrimaryVertices4D/matchedvtxsel/trecsimPU",2,1), "hist", color=2,  label="4d 0 PU")
#F.add(get(tcfpu,   "offlinePrimaryVertices/matchedvtxsel/trecsim_fromtracksPU",1,2), "hist", color=1,  label="3d 200 PU")
#F.add(get(tcfnopu, "offlinePrimaryVertices/matchedvtxsel/trecsim_fromtracksPU",2,2), "hist", color=1,  label="3d 0 PU")
#F.drawLegend(0.65, 0.65, 0.95, 0.95)


c3 = canvas("large")
F = F3 = Frame("4D vertices normalized time residuals", "(t_{rec}-t_{sim})/#sigma_{t} [ns]", -10., 10., "# reconstructed vertices")
F.add(get(tcfnopu, "offlinePrimaryVertices4D/matchedvtxsel/trecsimpullSignal",2,1), "hist", color=2,  label="t#bar{t} 0 PU")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/trecsimpullSignal",1,1), "hist", color=2,  label="t#bar{t} 200PU")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/trecsimpullPU",1,1, 0.01), "hist", color=4,  label="PU 200PU")
#F1.add(tcfnopu.Get("offlinePrimaryVertices4D/trecsimSignal"), "hist", color=1,  label="4d no PU")
#F.add(tcfpu.Get("offlinePrimaryVertices4D/matchedvtxsel/trecsimpull_fromtracksSignal"), "hist", color=1,  label="3d 200 PU")
F.drawLegend(0.65, 0.65, 0.95, 0.95)
F.Print("figures/4d-time-pulls.png")

raw_input()
sys.exit(0)

c5 = canvas("large")
F = F5 = Frame("t#bar{t}", "z_{rec}-z_{sim} [cm]", -0.01, 0.01, "# reconstructed vertices")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/zrecsimHRSignal",1,1), "hist", color=2,  label="4d 200PU")
F.add(get(tcfnopu, "offlinePrimaryVertices4D/matchedvtxsel/zrecsimHRSignal",2,1), "hist", color=2,  label="4d 0 PU")
F.add(get(tcfpu,   "offlinePrimaryVertices/matchedvtxsel/zrecsimHRSignal",1,1), "hist", color=1,  label="3d 200 PU")
F.add(get(tcfnopu, "offlinePrimaryVertices/matchedvtxsel/zrecsimHRSignal",2,1), "hist", color=1,  label="3d 0 PU")
F.drawLegend(0.65, 0.65, 0.95, 0.95)

c6 = canvas("large")
F = F6 = Frame("PU", "z_{rec}-z_{sim} [cm]", -0.01, 0.01, "# reconstructed vertices")
F.add(get(tcfpu,   "offlinePrimaryVertices4D/matchedvtxsel/zrecsimHRPU",1,1), "hist", color=2,  label="4d 200PU")
F.add(get(tcfpu,   "offlinePrimaryVertices/matchedvtxsel/zrecsimHRPU",1,1), "hist", color=1,  label="3d 200 PU")
F.drawLegend(0.65, 0.65, 0.95, 0.95)



raw_input()
