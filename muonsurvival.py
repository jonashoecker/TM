import math as m
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from array import array
import rootStyle

# parametrisation of momentum spectrum 
# Atmospheric muon flux at sea level, underground, and underwater
# E. V. Bugaev, A. Misaki, V. A. Naumov, T. S. Sinegovskaya, S. I. Sinegovsky, and N. Takahashi
# Phys. Rev. D 58, 054001 Published 15 July 1998

gam = [ [ 0.3061, 1.2743, -0.2630, 0.0252 ], [ 1.7910, 0.3040, 0, 0 ], [ 3.6720, 0, 0, 0 ], [ 4, 0, 0, 0 ] ]
prange = [ 1., 927., 1588., 416250. ]
const = [ 2.95e-3, 1.781e-2, 14.35, 1.e3 ]
def getIndex(p):
	ic = 0
	if p > prange[3]: 
		ic = 3
	elif p > prange[2]: 
		ic =2
	elif p > prange[1]: 
		ic =1
	return ic	

def pSpectrumPrimitive( p, c, g0, g1, g2, g3 ): # in GeV
	print 'coefficients:',p, c, g0, g1, g2, g3
	sgn = 1.
	if g2 < 0: sgn = -1.0
	g2 = m.fabs(g2)
	t0 = ( g1 * m.log(p) - 1. )**2 / (4. * sgn*g2 * m.log(p))
	t1 = sgn*m.sqrt(m.pi) * c * p**(-g0 - g3 * m.log(p)**3) * m.e**(t0) 
	t4 = (m.log(p) * (g1 + sgn*2. * g2 * m.log(p)) - 1.) / (2. * m.sqrt(g2) * m.sqrt(m.log(p)))
	t5 =  (2. * m.sqrt(g2) * m.sqrt(m.log(p)))
	print 't0,t1,t4,t5',t0,t1,t4,t5
	print 'p, log(p)',p,m.log(p)
	return t1 * m.erf(t4) / t5


def pSpectrumIntegral( p1, p2 ): # in GeV
	ic1 = getIndex(p1)
	ic2 = getIndex(p2)
	pP1 = pSpectrumPrimitive( p1, const[ic1], gam[ic1][0], gam[ic1][1], gam[ic1][2], gam[ic1][3] )
	pP2 = pSpectrumPrimitive( p2, const[ic2], gam[ic2][0], gam[ic2][1], gam[ic2][2], gam[ic2][3] )
	return pP2 - pP1

def pSpectrumNumIntegral( p1, p2 ): # in GeV
	npoints = 100000
	delta_p = (p2 - p1)/float(npoints)
	i_tot = 0.0
	for i in range(npoints):
		p = float(i)*delta_p + p1
		i_tot += pSpectrum(p)
	return i_tot * delta_p

def poly(p, c, g0, g1, g2, g3): # p in GeV/c
	lnp = m.log( p, 10 )
	retval = c * p**( - ( g0 + g1 * lnp + g2 * lnp**2 + g3 * lnp**3 ) )
	return retval

def pSpectrum(p): # p in GeV/c
	ic = getIndex(p)
	return poly(p, const[ic], gam[ic][0], gam[ic][1], gam[ic][2], gam[ic][3])
	
def PlotHistRateVsDepth():

	c = TCanvas("c","c",0,0,600,500)
	x, y = array( 'd' ), array ( 'd' )
	npoints = 1000
	dmin = 0
	dmax = 100
	delta_d = (dmax - dmin)/float(npoints)
	rate_0 = pSpectrumNumIntegral( 1.0, 3.e2 )

	for i in range(npoints):
		d = float(i)*delta_d + dmin
		p = (d - 0.54)/1.58
		if p < 1: p = 1
		rate_d = pSpectrumNumIntegral( p, 3.e2 )
		print 'in RateVsDepthHist: ',d,p,rate_d,rate_d/rate_0
		x.append(d)
		y.append(100*rate_d/rate_0)

	c.SetLogy()

	gr = TGraph( len(x), x, y)
	gr.SetMaximum(dmax*1.1)
	gr.SetMinimum(dmin*0.9)
	gr.Draw("ac")
	gr.SetLineWidth(4)
	c_Bethe = TColor.GetColor( '#BD265D' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )


	gr.GetXaxis().SetTitle("profondeur [m]")
	gr.GetYaxis().SetTitle("taux relatif [%]")

	c.Update()
	c.Print("rvsd.pdf")




def PlotHisto():

	c = TCanvas("c","c",0,0,600,500)

	x, y = array( 'd' ), array ( 'd' )
	nmax = 0.0
	nmin = 1e10
	pmin = 2.
	pmax = 1.e2
	npoints = 10000
	delta_p = (pmax - pmin)/float(npoints)
	i_tot = 0.0
	p_part = 30 # in GeV
	i_part = 0.0
	p_ave = 0.
	for i in range(npoints):
		p = float(i)*delta_p + pmin
		nmuon = pSpectrum(p)
		# print p, nmuon, nmuon*p**3
		x.append(p)
		y.append(nmuon)
		p_ave += p*nmuon
		i_tot += nmuon*delta_p
		if p > p_part: i_part += nmuon*delta_p
		if y[-1] > nmax: nmax = y[-1]
		if y[-1] < nmin: nmin = y[-1]

		
	print 'total numerical integral  :',i_tot,pSpectrumNumIntegral( pmin, pmax )
	### i_a = pSpectrumIntegral( pmin, pmax )
	### print 'total analytical integral :',i_a
	print 'average momentum:',p_ave/float(i_tot)
	print 'integral for p > ',p_part,' GeV:', i_part,' ==> fraction = ',i_part/i_tot
	#c.SetLogx()
	c.SetLogy()
	#c.SetLogx()

	gr = TGraph( len(x), x, y)
	gr.SetMaximum(nmax*1.1)
	gr.SetMinimum(nmin*0.9)
	gr.Draw("al")
	gr.SetLineWidth(4)
	c_Bethe = TColor.GetColor( '#1A52A9' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )

	gr.GetXaxis().SetTitle("impulsion (GeV/c)")
	gr.GetYaxis().SetTitle("nombre de muons (echelle arbitraire)")

	c.Update()
	c.Print("ms.pdf")



def main():
	#PlotHisto()
	PlotHistRateVsDepth()


if __name__ == "__main__":
	main()
