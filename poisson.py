import math as m
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from array import array
import rootStyle


#def probability(ev):
	#res = landa**ev * m.e**(-landa) / m.factorial(ev)
	#return res

#print 'the probability is :', 	probability(ev)

landa = 10.54

def PlotGraph():
	ev = 0

	x, y = array( 'd' ), array( 'd' )
	cont = True
	while(cont):
		px = landa**ev * m.e**(-landa) / m.factorial(ev)
		ev += 1
		x.append(ev)
		y.append(px*100)
		if ev > 20:
			cont = False

	c = TCanvas("c1","c1",0,0,600,500)
	c.SetGrid()

	gr = TGraph( len(x), x, y )
	gr.SetMaximum(20)
	gr.SetMinimum(0)
	gr.Draw("al")
	gr.SetLineWidth(2)
	c_Bethe = TColor.GetColor( '#1A52A9' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )

	gr.GetXaxis().SetTitle("event per 10s intervals")
	gr.GetYaxis().SetTitle("probability in %")
	c.Update()
	c.Print("alfa.pdf")


def main():

	PlotGraph()


if __name__ == "__main__":
	main()
