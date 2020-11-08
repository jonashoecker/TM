from ROOT import *

atlasStyle = TStyle("ATLAS", "Atlas style")

#  use plain black on white colors
icol = 0
atlasStyle.SetFrameBorderMode(icol)
atlasStyle.SetCanvasBorderMode(icol)
atlasStyle.SetPadBorderMode(icol)
atlasStyle.SetPadColor(icol)
atlasStyle.SetCanvasColor(icol)
atlasStyle.SetStatColor(icol)
# atlasStyle.SetFillColobr(icol)

#  set the paper & margin sizes
atlasStyle.SetPaperSize(20, 26)
atlasStyle.SetPadTopMargin(0.05)
atlasStyle.SetPadRightMargin(0.05)
atlasStyle.SetPadBottomMargin(0.16)
atlasStyle.SetPadLeftMargin(0.10)

#  use large fonts
font = 42
tsize = 0.04
atlasStyle.SetTextFont(font)


atlasStyle.SetTextSize(tsize)
atlasStyle.SetLabelFont(font, "x")
atlasStyle.SetTitleFont(font, "x")
atlasStyle.SetLabelFont(font, "y")
atlasStyle.SetTitleFont(font, "y")
atlasStyle.SetLabelFont(font, "z")
atlasStyle.SetTitleFont(font, "z")

atlasStyle.SetLabelSize(tsize, "x")
atlasStyle.SetTitleSize(tsize, "x")
atlasStyle.SetLabelSize(tsize, "y")
atlasStyle.SetTitleSize(tsize, "y")
atlasStyle.SetLabelSize(tsize, "z")
atlasStyle.SetTitleSize(tsize, "z")


# use bold lines and markers
atlasStyle.SetMarkerStyle(0)
atlasStyle.SetMarkerSize(1.2)
atlasStyle.SetHistLineWidth(2)
atlasStyle.SetLineStyleString(2, "[12 12]")  # postscript dashes

# get rid of X error bars and y error bar caps
# atlasStyle.SetErrorX(0.001)

# do not display any of the standard histogram decorations
atlasStyle.SetOptTitle(0)
# atlasStyle.SetOptStat(1111)
atlasStyle.SetOptStat(0)
# atlasStyle.SetOptFit(1111)
atlasStyle.SetOptFit(0)

#  put tick marks on top and RHS of plots
atlasStyle.SetPadTickX(1)
atlasStyle.SetPadTickY(1)

gROOT.SetStyle("Plain")

# gStyle.SetPadTickX(1)
# gStyle.SetPadTickY(1)
gROOT.SetStyle("ATLAS")
gROOT.ForceStyle()
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
#  overwrite atlas styles

atlasStyle.SetMarkerSize(1.0)
atlasStyle.SetPadLeftMargin(0.14)
atlasStyle.SetPadRightMargin(0.03)
atlasStyle.SetPadBottomMargin(0.12)
atlasStyle.SetPadTopMargin(0.02)
atlasStyle.SetFrameFillColor(0)
