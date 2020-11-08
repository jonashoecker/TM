import math as m
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from array import array
import rootStyle

mmu = 105.7 / 1000. # mass of muon
c = 299792458 # m / s
tau = 2.1969811e-6  # seconds



def FillArray( p ):

	betagamma = p/mmu
	beta = betagamma/m.sqrt(1+betagamma**2)
	gamma = 1./m.sqrt(1-beta**2)

	x = 15000. # m
	dx = 100.0 # m
	ddx = 0. # m
	xg, yg = array( 'd' ), array( 'd' )
	cont = True
	while (cont):
		xg.append(x/1000.)
		proba = m.exp( - ddx / ( c * gamma * tau ) )
		yg.append(proba)
		x -= dx
		ddx += dx
		if x < 0: cont = False

	return xg, yg


def PlotHisto():

	x10, y10 = FillArray(10)
	x6, y6 = FillArray(6)
	x4, y4 = FillArray(4)
	x2, y2 = FillArray(2)
	x1, y1 = FillArray(1)

	c = TCanvas("c1","c1",0,0,600,500)
	frame = TH2F("frame","frame",100,0,18,10,0,1)
	frame.Draw()
	frame.GetXaxis().SetTitle("Altitude (km)")
	frame.GetYaxis().SetTitle("Probabilit#acute{e} de survie (%)")

	gra = [ TGraph( len(x10), x10, y10 ), TGraph( len(x6), x6, y6 ), TGraph( len(x4), x4, y4 ), TGraph( len(x2), x2, y2 ), TGraph( len(x1), x1, y1) ]
	coa = [ TColor.GetColor( '#4C1FC4' ), TColor.GetColor( '#1F97C4' ), TColor.GetColor( '#1EC22A' ), TColor.GetColor( '#F7A126' ), TColor.GetColor( '#F74226' ) ]
	naa = [ '10 GeV/c', '6 GeV/c', '4 GeV/c', '2 GeV/c', '1 GeV/c' ]

	legend = TLegend(0.58,0.17,1.03,0.47)
	legend.SetHeader('Impulsion de d#acute{e}part (GeV/c)')
	legend.SetTextFont()
	legend.SetBorderSize(0)
	legend.SetFillStyle(0)

	ic = 0
	for gr in gra:
		gr.SetLineWidth(2)
		gr.SetLineColor(coa[ic])
		gr.Draw("c")
		legend.AddEntry( gr, naa[ic], 'l' )
		ic += 1

	legend.Draw()

	# indicate creation of muons
	line = TLine(15, 0.88, 15, 1)
	line.SetLineWidth(1)
	line.SetLineColor(1)
	line.SetLineStyle(3)
	line.Draw()

	text = TLatex()
	text.SetTextSize(0.03)
	text.DrawLatex(14.5,0.85,"Point de ")
	text.DrawLatex(14.5,0.82,"cr#acute{e}ation ")
	text.DrawLatex(14.5,0.79,"des muons")

	c.Update()
	c.Print("boa.pdf")


def main():
	PlotHisto()

if __name__ == "__main__":
	main()