from rutil import Frame, canvas, TCachingFile, parseOpt
from ROOT import TMath, TF1

args = parseOpt(file="v22/pvt_RUN4PU200TTBAR__A.root", name="")

tcf = TCachingFile(args.file)
hsel3d = tcf.Get("offlinePrimaryVertices/ntpselvssimPU")
hsel4d = tcf.Get("offlinePrimaryVertices4D/ntpselvssimPU")
hreal3d = tcf.Get("offlinePrimaryVertices/ntprealselvssimPU")
hreal4d = tcf.Get("offlinePrimaryVertices4D/ntprealselvssimPU")
hsplit3d = tcf.Get("offlinePrimaryVertices/ntpsplitselvssimPU")
hsplit4d = tcf.Get("offlinePrimaryVertices4D/ntpsplitselvssimPU")
hother3d = tcf.Get("offlinePrimaryVertices/ntpotherfakeselvssimPU")
hother4d = tcf.Get("offlinePrimaryVertices4D/ntpotherfakeselvssimPU")
hfake3d = tcf.Get("offlinePrimaryVertices/ntpfakeselvssimPU")
hfake4d = tcf.Get("offlinePrimaryVertices4D/ntpfakeselvssimPU")

c1 = canvas("large")
F = F1 = Frame("200 PU t#bar{t}", "PU", 0., 320., "# reconstructed vertices", 0., 160)
F1.add(hreal3d, "fullcircle", color=1,  label="3d matched")
F1.add(hreal4d, "fullcircle", color=2,  label="4d matched")
#F1.add(hsplit3d, "fullcircle",  label="3d", add_to_legend = False)
#F1.add(hsplit4d, "fullcircle",  label="4d")
#F1.add(hother3d, "fullcircle",  label="3d")
#F1.add(hother4d, "fullcircle",  label="4d")
F1.add(hfake4d, "opencircle", color=2, label="4d fake")
F1.add(hfake3d, "opencircle", color=1, label="3d fake")

if True:
    sigmaz=4.26 # beamspot
    dzeff = 0.018
    alpha=TMath.Erf(dzeff/sigmaz/2.)
    epsilon = 0.70
    print alpha
        
    #lf3d = TF1("f","x*%f-0.5*%f*x*x"%(epsilon, epsilon**2*alpha),0., 300.)
    lf3d = TF1("f","x*[0]-0.5*[1]*x*x",0., 300.)
    lf3d.SetParameter(0, epsilon)
    lf3d.SetParameter(1,  epsilon**2*alpha)
    lf3d.SetLineColor(1)
    lf3d.SetLineWidth(1)
    lf3d.SetLineStyle(2)
    hreal3d.Fit(lf3d, "")
    #lf3d.Draw("same")
    eff3d_200 =  lf3d.GetParameter(0) -0.5 * lf3d.GetParameter(1)*200
    
    dzeff = 0.0185
    alpha=TMath.Erf(dzeff/sigmaz/2.)
    epsilon = 0.67
    print alpha
    
    #lf4d = TF1("f","x*%f-0.5*%f*x*x"%(epsilon, epsilon**2*alpha),0., 300.)
    lf4d = TF1("f","x*[0]-0.5*[1]*x*x",0., 300.)
    lf4d.SetParameter(0, epsilon)
    lf4d.SetParameter(1,  epsilon**2*alpha)
    lf4d.SetLineColor(2)
    lf4d.SetLineWidth(1)
    lf4d.SetLineStyle(1)
    #lf4d.Draw("same")
    hreal4d.Fit(lf4d, "N")
    eff4d_200 =  lf4d.GetParameter(0) -0.5 * lf4d.GetParameter(1)*200
    
F1.drawLegend(0.65, 0.3, 0.95, 0.65)
if args.name == "":
    F1.Print("numvtvsnsim.png")
else:
    F1.Print("numvtvsnsim_%s.png"%args.name)


print "3d efficiency(@200) ",eff3d_200
print "4d efficiency(@200) ",eff4d_200
raw_input()
