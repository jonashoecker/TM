import math as m
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from array import array
import rootStyle
import muonsurvival as ms
import b


def Plot():
	c = TCanvas("c","c",0,0,600,500)

	x, y = array( 'd' ), array ( 'd' )
	ex, ey = array( 'd' ), array ( 'd' )

	x.append(413.) # altitude gnv
	x.append(1030) # altitude chamonix
	x.append(1719) # altitude recule
	x.append(2317.) # altitude intermediate station
	x.append(3790.) # altitude sommet midi
	for i in range(5): ex.append(0.)
	counts = [ 5768., 4855., 6342., 5008., 12961. ]
	time   = [ 2400., 1800., 1800., 1200., 1800. ]
	for i in range(5):
		y.append(counts[i] / time[i])
		ey.append(m.sqrt(counts[i]) / time[i])
		print i, x[i],y[i], ey[i]

	h = TH1F("frame","frame",100,0,5000.)
	h.GetXaxis().SetTitle("Altitude [m]")
	h.GetYaxis().SetTitle("Taux [Hz]")
	h.SetMaximum(8)
	h.SetMinimum(0.0)
	h.Draw()

	gr = TGraphErrors( 5, x, y, ex, ey )
	gr.SetMarkerStyle(20)
	gr.SetMarkerSize(0.45)
	gr.SetLineWidth(1)

	# fit
	fit = TF1("it", "[0]*exp([1]*x)+[2]"); 
	fit.SetLineColor(2)
	fit.SetLineWidth(1)
	gr.Fit(fit); 
	gr.Draw('p')
	c.SetGrid()

	# show lifetime expectation
	lifetime = TF1("lt","7.2*exp(-(3790.-x)/(299792458.*(4./0.105)*2.2e-6))",0.,3790.)
	lifetime.SetLineWidth(2)
	lifetime.SetLineColor(4)
	lifetime.SetLineStyle(2)
	#lifetime.Draw("same")

	c.Update()
	c.Print('CAENair.pdf')

	c.Print("rvsdair.pdf")

def main():
	Plot()



if __name__ == "__main__":
	main()