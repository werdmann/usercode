from ROOT import TFile, TH1F, TH2F, TLegend, gPad, TCanvas, TGaxis, gStyle, TGraph, TGraphErrors, TObject, TLatex, TLine, TArrow

import array
from types import *

import random,math,re,sys,os
import hashlib

def keep(*items):
    global repository_
    try:
        n=len(repository_)
    except:
        repository_={}
    for i in items:
        randomkey=random.random()
        repository_[randomkey]=i

def nostat():
    gStyle.SetOptStat(0)

def randid():
    return "".join([random.sample("ABCDEFGHIJKLMNIOQRSTUVWXYZ",1)[0] for n in range(20)])


class TCachingFile:
    files = []
    def __init__(self, path, color=1, label=None, norm=None, caching=True, defaultpath=""):
        if not os.path.exists(path):
            if not os.path.exists( os.path.join(defaultpath, path) ):
                print "TCachingFile file not found ",path
            path = os.path.join(defaultpath, path)

        self.path=path
        self.tfile=TFile(path)
        self.color=color
        self.label=label
        self.cache={}
        self.cwd=""
        self.norm=None
        self.caching = caching
        if norm is not None:
            if isinstance(norm, float):
                self.norm=norm
            elif isinstance(norm, FunctionType):
                self.norm=norm(self.tfile)
                
        info=self.tfile.Get("Info") #TObjString
        if info:
            self.info=str(info.String())
        else:
            self.info=None
          
        build = self.tfile.Get("Build")
        if build:
            self.build=str(build.String())
        else:
            self.build=None
        
        self.id = len(self.files)
        self.files.append(self)

    def Get(self, hid):
        if self.cwd == "":
            p = hid
        else:
            p = self.cwd+"/"+hid
        m = hashlib.md5()
        m.update(self.cwd)
        cp=m.hexdigest()+":"+p
        if not cp in self.cache:
            h0 = self.tfile.Get(p)
            if h0:
                h = h0.Clone(cp)
            else:
                print "histogram not found : ",p
                return None
            if self.norm:
                h.Scale(self.norm)
            if self.caching: self.cache[cp]=h
            h.SetLineColor(self.color)
            h.SetMarkerColor(self.color)
            h.file=self
            return h
        else:
            return self.cache[cp]


    def cd(self, dir=""):
        self.cwd=dir

    def close(self):
        self.tfile.Close()

    def __getitem__(self, hid):
        return self.Get(hid)



def normalize(h):
    if h.Integral()>0:
        h.Scale(1./h.Integral())
        



class Frame:

    markers={ "opencircle":(24,1), "smallopencircle":(24,0.5),
              "fullcircle":(20,1), "smallfullcircle":(20,0.5),
              "opensquare":(25,1), "smallopensquare":(25,0.5),
              "fullsquare":(21,1), "smallfullsquare":(21,0.5),
              "opentriangle":(26,1), "smallopentriangle":(26,0.5),
              "fulltriangle":(22,1), "smallfulltriangle":(22,0.5),
              "opencross":(28,1), "smallopencross":(28,0.5),
              "fullcross":(34,1), "smallfullcross":(34,0.5),
              "opendiamond":(27,1), "smallopendiamond":(27,0.5),
              "fulldiamond":(33,1), "smallfulldiamond":(33,0.5)
              }
    #http://root.cern.ch/root/html/TAttMarker.html
    
    def __init__(self,title, xtitle="",xmin=0,xmax=1.,ytitle="", ymin=0, ymax=0., option="", nxbin=2, nybin=2, logy=False, logx=False):
        randomid="%d"%(random.randint(10000000000,99999999999))
        self.ymin=ymin
        self.ymax=ymax
        self.xmin=xmin
        self.xmax=xmax
        self.palette = (1,2,4,8,9,42, 46, 38, 31)
        self.logy = logy
        if logy and self.ymin <= 0:
            self.ymin= 1.e-5

        self.logx = logx
        if logx and self.xmin <= 0:
            self.xmin= 1.e-5
            
        if ymax > ymin:
            self.hF=TH2F(randomid,title,nxbin,xmin,xmax,nybin,ymin,ymax)
        else:
            self.hF=TH2F(randomid,title,nxbin,xmin,xmax,nybin,ymin,ymin+1)
            ymax = ymin + 1

        self.yaxis2 = None
        self.hF.GetXaxis().SetTitle(xtitle)
        self.hF.GetYaxis().SetTitle(ytitle)

        if logy:
            gPad.SetLogy()
        if logx:
            gPad.SetLogx()
        self.hF.SetStats(0)
        self.hF.Draw()
        self.pad = gPad.cd() # yes, the cd() is needed, gPad is special

        self.items=[]
        self.options=[]
        self.labels=[]
        # list of options (e.g. logy), comma or white-space separated
        self.frame_options = option.replace(","," ").split()

        self.tls = []
        self.tlines = []

    def process_options(self, h, opts):
        """ FIXME, work in progress, to be used in add and drawLegend consistently """
        if opts == "":
            optstring = h.GetDrawOption()
            if optstring == "same":
                optsring = "L"
        else:
            optstring = ""
            for o in opts.split("+"):
                if o=="hist":
                    optstring += "L"
                elif o.lower() in self.markers:
                    optstring +="P"
                elif o == "L":
                    optstring += "L"
                elif o == "C":
                    optstring += "C"
                else:
                    print "unknown option : ",o
        print "optstring ",optstring, "opts" ,opts
            
    def auto_color(self):
        return self.palette[ (len(self.items) - 1) % len(self.palette) ]

    
    def auto_marker_color(self):
        markers = ("fullcircle", "opencircle", "fullsquare", "opensquare", "fulltriangle", "opentriangle")
        n = len(self.items)-1
        m = n % len(markers)
        c = (m + n // len(markers)) % len(self.palette)
        return "small"+markers[m], self.palette[c]
        
        
    def add(self, item, drawOption="", label=None, color=None, size=None, f=None, axis=1, add_to_legend=True):
        self.items.append(item)
        self.options.append(drawOption)

        if add_to_legend == False:
            self.labels.append(None)
        elif label:
            self.labels.append(label)
        else:
            try:
                if item.file.label:
                    self.labels.append(item.file.label)
                elif item.file.info:
                    self.labels.append(item.file.info)            
                else:
                    self.labels.append("")
            except AttributeError:
                # histogram may not have a file attribute
                self.labels.append("")
                
        if color is None and drawOption == "auto":
            drawOption, color = self.auto_marker_color()
            self.options[-1] = drawOption

        elif color is None:
            color = self.auto_color()
            
        item.SetLineColor(color)
        item.SetMarkerColor(color)

        if size is not None:
            item.SetMarkerSize(size)


        if drawOption=="P" and item.GetMarkerStyle()==1:
            item.SetMarkerStyle(20)
            
        if drawOption.lower() in self.markers:
            item.SetMarkerStyle(self.markers[drawOption.lower()][0])
            item.SetMarkerSize(self.markers[drawOption.lower()][1])
            drawOption="P"

        if drawOption.startswith("L+") and drawOption[2:].lower() in self.markers:
            item.SetMarkerStyle(self.markers[drawOption[2:].lower()][0])
            item.SetMarkerSize(self.markers[drawOption[2:].lower()][1])
            drawOption="PL"

            
        if drawOption == "fill":
            drawOption = "hist"
            if color:
                item.SetFillColor(color)
        
        if drawOption == "box":
            item.SetFillColor(0)
            item.SetLineColor(color)
            

        if f is not None:
            f(item)

        if axis == 2 and self.yaxis2 is not None:
            yscale = (self.ymax - self.ymin) / (self.ymax2 - self.ymin2)
            # rescale for the second axis
            if type(item) is TGraph:
                xdata = [x for x in item.GetX()]
                ydata = [self.ymin +(y - self.ymin2) * yscale for y in item.GetY()]
                for n in range(len(ydata)):
                    item.SetPoint(n, xdata[n], ydata[n])
            
        self.pad.cd()
        item.Draw("same0 "+drawOption)
        self.pad.Update()

    def addLatex(self, x, y, text, textsize=0.3, textcolor=1):
        tl = TLatex(x, y, text)
        tl.SetTextSize(textsize)
        tl.SetTextColor(textcolor)
        tl.Draw("same")
        self.tls.append(tl)
        return tl

    def addLine(self, x1, y1, x2, y2, arrow="", arrowsize=0.1):
        if arrow == "":
            tl = TLine(x1, y1, x2, y2)
        elif arrow == "<->":
            tl = TArrow(x1, y1, x2, y2, arrowsize, "<|>")
        elif arrow == "->":
            print x2, y2, arrowsize
            tl = TArrow(x1, y1, x2, y2, arrowsize, "|>")
            tl.SetAngle(40)
            tl.SetLineWidth(1)
            tl.SetFillColor(1)
        tl.Draw("same")
        self.tlines.append(tl)
        return tl
            
    def drawLegend(self,xl=0.8,yl=0.8,xh=0.95,yh=0.95,title="",columns=1):
        self.pad.cd()
        self.legend=TLegend(xl, yl, xh, yh, title)
        self.legend.SetNColumns(columns)
        self.legend.SetFillColor(0)
        for h,l,o in zip(self.items,self.labels,self.options):
            if l is None or l=="":
                continue
            if o=="":
                o=h.GetDrawOption()
                if o=="same":
                    o="L"
            if o=="hist" or o=="H":
                o="L"
            elif o.lower() in self.markers:
                o="P"
            self.legend.AddEntry(h,l,o)
        self.legend.Draw()
        self.Draw()
        
    def DrawLegend(self,xl=0.8,yl=0.8,xh=0.95,yh=0.95,title="",columns=1):
        self.drawLegend(xl, yl, xh, yh, title, columns)

    def Draw(self):
        ymax = self.ymax
        if self.ymax <= self.ymin:
            for h in self.items:
                try:
                    y = h.GetBinContent(h.GetMaximumBin())*1.10
                except AttributeError:
                    try: 
                        y = max( h.GetY() )
                    except AttributeError:
                        y = ymax

                if y > ymax: ymax = y
        
        if ymax <= self.ymin:
            ymax = self.ymin+1
        self.hF.GetYaxis().SetLimits(self.ymin,ymax)

        if "logy" in self.frame_options:
            self.pad.SetLogy()
        if "logx" in self.frame_options:
            self.pad.SetLogx()
        if "grid" in self.frame_options or "gridxy" in self.frame_options:
            self.pad.SetGrid()
        if "gridx" in self.frame_options:
            self.pad.SetGridx()
        if "gridy" in self.frame_options:
            self.pad.SetGridy()

        self.pad.Update()

    def Print(self,filename=None):
        if filename:
            self.pad.Print(filename)
       
    def addVerticalAxis(self, title, ymin, ymax):
        self.yaxis2 = TGaxis(self.xmax,self.ymin, self.xmax, self.ymax,
                             ymin, ymax,510,"+L")
        self.yaxis2.SetName("axis")
        self.yaxis2.SetLabelOffset(0.01)
        self.yaxis2.SetTitle(title)
        self.yaxis2.Draw()
        self.ymin2 = ymin
        self.ymax2 = ymax
        return self.yaxis2
        
    def getXaxis(self):
        """ return the x-axis (for setting attributes) """
        return self.hF.GetXaxis()
    
    def getYaxis(self):
        """ return the y-axis (for setting attributes) """
        return self.hF.GetYaxis()

    def getYaxis2(self):
        """ return the second y-axis (note: this is a TGaxis, not a TAxis)"""
        return self.axis2
                                 
def frame(title, xtitle="",xmin=0,xmax=1.,ytitle="", ymin=0, ymax=0.):
    randomid="%d"%(random.randint(10000000000,99999999999))
    hF=TH2F(randomid,title,2,xmin,xmax,2,ymin,ymax)
    hF.GetXaxis().SetTitle(xtitle)
    hF.GetYaxis().SetTitle(ytitle)
    hF.Draw()
    keep(hF)
    return hF


def tgraph_from_tuples( xylist, sort=False ):
    if sort:
        xylist = sorted(xylist)
    x,y=( array.array('d',l) for l in zip(*xylist) )
    tg = TGraph(len(xylist), x, y)
    keep(tg)
    tg.SetFillColor(0)
    return tg

def tgraph_from_list( a ):
    n=len(a)
    tg = TGraph(n, array.array('d', range(n)), array.array('d', a) )
    keep(tg)
    tg.SetFillColor(0)
    return tg

def tgraph_from_lists( x, y, sort = False):
    n=min(len(x), len(y))
    if sort:
        xylist = sorted( zip(x,y) )
        x,y=( array.array('d',l) for l in zip(*xylist) )
    tg = TGraph(n, array.array('d', x[:n]), array.array('d', y[:n]) )
    keep(tg)
    tg.SetFillColor(0)
    return tg

def tgrapherrors_from_lists( x, y, ex, ey ):
    if ex is None and ey is None:
        n=min(len(x), len(y))
        tg = TGraphErrors(n, array.array('d', x[:n]),
                          array.array('d', y[:n]),
                          array.array('d', [0]*n),
                          array.array('d', [0]*n))
    elif ex is None:
        n=min(len(x), len(y), len(ey))
        tg = TGraphErrors(n, array.array('d', x[:n]),
                          array.array('d', y[:n]),
                          array.array('d', [0]*n),
                          array.array('d', ey[:n]))
    elif ey is None:
        n=min(len(x), len(ex), len(y))
        tg = TGraphErrors(n, array.array('d', x[:n]),
                          array.array('d', y[:n]),
                          array.array('d', ex[:n]),
                          array.array('d', [0]*n))
    else:
        n=min(len(x), len(ex), len(y), len(ey))
        tg = TGraphErrors(n, array.array('d', x[:n]),
                          array.array('d', y[:n]),
                          array.array('d', ex[:n]),
                          array.array('d', ey[:n]))
    keep(tg)
    tg.SetFillColor(0)
    return tg


def canvas(option, title="", nx=1, ny=1, n=0):
    """ create a canvas and optionally divide it into subcanvasses
    option determines the size/shape
    tile   title
    nx,ny  divide the canvas into nx by ny subpads
    n      determine a suitable subdivision for n pads if n > 0
    """
    if n > 0:
        nx = int( math.ceil(math.sqrt( n )))
        ny = int( math.ceil(float(n) / nx))

    size=1200
    global canvascounter
    try:
        canvascounter +=1
    except:
        canvascounter =1 
    cname = "c%d"%canvascounter
    if option.startswith("x"):
        option = option[1:]
        size = size * 2
    if option=="landscape":
        c=TCanvas(cname, title,int(size),int(size/math.sqrt(2)))
    elif option=="wide":
        c=TCanvas(cname,title,int(size),int(size/2))
    elif option=="portrait":
        c=TCanvas(cname,title,int(size/math.sqrt(2)),int(size))
    elif option=="small":
        c=TCanvas(cname,title,int(size/2),int(size/2))
    elif option=="postcard":
        c=TCanvas(cname,title,int(size/2),int(size/4))
    elif option=="smallwide":
        c=TCanvas(cname,title,int(size/2),int(size/4))
        gStyle.SetLabelSize(.08, "XY");
    elif option=="large":
        c=TCanvas(cname,title,int(size),int(size))
    elif option=="xlarge":
        c=TCanvas(cname,title,int(size*1.5),int(size*1.5))
    elif option=="xlandscape":
        c=TCanvas(cname,title,int(size*math.sqrt(2)),size)
    else:
        print "canvas format option not understood : ",option
        return None
        
    keep(c)
    if nx>1 or ny>1:
        c.Divide(nx,ny)
        c.cd(1)
    else:
        c.cd()
    return c

class ArgContainer:
    def __init__(self):
        self.keys=[]
        self.args=[] # non-keyword arguments
    def set(self, key, value):
        self.keys.append(key)
        setattr(self, key, value)
    def __getitem__(self, key):
        return getattr(self, key)
    
def parseOpt(**opt):
    """ parse command line arguments and returns a dummy object
    that has attributes corresponing to the arguments
    
    example
        args=parseOpt(filename="pv.root")


    creates an attribute named filename. If a commandline argument -filename=foo has been given,
    the value of filename will be foo, if not, it will be pv.rootn

    normally argument types are string, but "True" or "False" are 
    converted to boolean


    """
    arg=ArgContainer()
    
    opt_t={}
    for key in opt:
        opt_t[key]=type(opt[key])
    
    for a in sys.argv:
        if a in ("help","-h","--h","--help"):
            print sys.argv[0],
            for key in opt: print "(%s=%s)"%(key,opt[key]),
            print
            sys.exit(0)
            
    for a in sys.argv[1:]:
        m=re.match("-*(\S+)=(\S+)",a)
        if m:
            keyword,value=m.groups()
            if keyword in opt:
                if keyword in opt_t:
                    opt[keyword]=opt_t[keyword](value)
                else:
                    opt[keyword]=value
            else:
                print "unknown option",a 
        else:
            arg.args.append(a)

    for o in opt:
        #globals()[o]=opt[o]
        if opt[o]=="True":
            #setattr(arg,o,True)
            arg.set(o,True)
        elif opt[o]=="False":
            #setattr(arg,o,False)
            arg.set(o,False)
        else:
            #setattr(arg,o,opt[o])
            arg.set(o,opt[o])
        #print o,opt[o]

    return arg



def parseCfg(filename="default.cfg"):
    """ parse a config file like an opt list, one line per option
    
    example
        args=parseCfg(filename="hello.cfg")

     where hello.cfg looks like
     rootfile hello.root
     bla 5
     bla 6
     

     The latter will create an array bla=[5,6]
    """

    arg=ArgContainer()
    try:
        f=open(filename)
    except IOERROR:
        print "rutil.parseCfg: file not found ",filname
        return

    opt={}
    for a in f.readlines():

        if a.startswith("#"):
            continue
        
        m=re.match("(\S+) (.*)",a)
        if m:
            keyword,value=m.groups()
            if value=="True":
                value=True
            elif value=="False":
                value=False

                
            if keyword in opt:
                opt[keyword].append(value)
            else:
                opt[keyword]=[value]
                
        else:
            print "option not understood: ",a
                
    for o in opt:
        if len(o)==1:
            setattr(arg,o,opt[o][0])
        else:
            setattr(arg,o,opt[o])

    return arg



class Plot(object):
    """ container holding information about a plot (a histogram)
        mostly useful when used in a PlotGroup

    """
    def __init__(self, file_path, hid, label, color, option, scale=None):
        self.path = file_path
        self.hid = hid  # may contain a path
        self.label = label
        self.color = color
        self.option = option
        self.scale = scale  # None or a string descriptor
        self.tfile = None  # defined when openend (PlotGroup has the full path)
        self.hist = None
        
    def get_norm(self,root_folder=""):
        if self.scale == None:
            return 1.
        try:
            return float(self.scale)
        except ValueError:
            pass

        try:
            hid, quantity = self.scale.split(".")
            h = self.tfile.Get(root_folder + hid)
            if quantity == "entries":
                return h.GetEntries()
        except ValueError:
            print "require path.quantity format, don't understand ",self.scale
            return 1.

        
    def get(self, dir_path=".", root_folder=""):
        if not os.path.exists(os.path.join(dir_path, self.path)):
            print "*"*80
            print "file not found: ",os.path.join(dir_path, self.path)
            print "*"*80
            return TH1F("","",1,0.,1.)
        
        self.tfile = TCachingFile(os.path.join(dir_path, self.path))
        self.hist = self.tfile.Get(root_folder + self.hid)
        if self.scale is not None:
            norm = self.get_norm(root_folder)
            if norm > 0:
                self.hist.Scale(1/norm)
        return self.hist

            
class PlotGroup(object):
    def __init__(self, description="", dir_path=".", root_folder="", group_label = "", plots=[]):
        self.description = description

        if root_folder == "" or root_folder.endswith("/"):
            self.root_folder = root_folder
        else:
            self.root_folder =  root_folder + "/"

        if group_label == "" or group_label.startswith(" "):
            self.group_label = group_label
        else:
            self.group_label = " " + group_label
        
        self.path = dir_path
        self.plots = []
        self.tf = []
        for p in plots:
            self.add_plot(p)
        
    def add(self, p):
        self.plots.append(p)

    def plot_on_frame(self, f, add_to_legend=True):
        for p in self.plots:
            hist = p.get(self.path, self.root_folder)
            f.add(hist, p.option, p.label + self.group_label, p.color,
                  add_to_legend = add_to_legend)
                         
