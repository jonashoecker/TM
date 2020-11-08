import math as m
from random import seed
from random import random
from ROOT import *
from array import array
import rootStyle
import muonsurvival as muon

a = 15. # size of scintillator in cm
rx = a/2.0
dy = 8.0 # distance between scintillators in cm
y0 = dy

seed(1) # set random seed

def getWeight( phi ):
	theta = (random()-0.5)*m.pi
	return theta, m.cos(theta + phi*m.pi/180.0)**2

def hitBoth(x, theta):
	tanth = m.tan(theta)
	# upper
	y = y0 - dy
	dx = y*tanth
	if m.fabs(x + dx) > a/2.: return False
	# lower
	dx = y0*tanth
	if m.fabs(x + dx) > a/2.: return False
	return True

def simMuon( phi ):
	ntracks = 1000000
#	ntracks = 10000
	sum_w = 0.
	for i in range(ntracks):
		# generate x coordinate
		x = 2.0*(random()-0.5)*rx
		# generate angle and weight
		theta, w = getWeight(phi)
		# compute coordinate of lower plate
		if hitBoth(x, theta): sum_w += w
	return sum_w/float(ntracks)

def phiLoop():

	c = TCanvas("c","c",0,0,600,500)
	x, y = array( 'f' ), array ( 'f' )
	phi = array( 'f' )

	nphi = int(90.0/2.)
	dphi = 1./nphi
	w0 = simMuon( 0 )

	dist0 = 20.0
	rate_0 = muon.pSpectrumNumIntegral( 1.0, 3.e2 )
	ave_surv = 0.0
	sum_weight = 0.0
	sum_weight1 = 0.0

	for i in range(nphi+1):
		phi.append( 90*(i)/float(nphi) )
		w = simMuon( phi[i] )
		x.append(phi[i])
		y.append(100.0*w/w0)

		if i > 0:
			alpha = (phi[i-1] + phi[i])/2.
			weight = (y[i-1] + y[i])/2./100.
			dist = dist0/m.cos(alpha*m.pi/180.0)
			pmin = (dist - 0.54)/1.58
			rate_d = muon.pSpectrumNumIntegral( pmin, 3.e2 )
			p_surv = rate_d/rate_0
			ave_surv += p_surv*weight
			sum_weight += weight
			if (dist < 19): sum_weight1 += weight

			print alpha, weight, dist, pmin, p_surv

	ave_surv /= sum_weight
	print 'sum of weights:',sum_weight,sum_weight1,sum_weight1/sum_weight
	print '--------------------------------------'
	print 'Average survival probability: ',ave_surv
	print '--------------------------------------'

	# plot frame
	h = TH1F("frame","frame",100,0,90)
	h.GetXaxis().SetTitle("Angle des scintillateurs [deg]")
	h.GetYaxis().SetTitle("Taux relatif [%]")
	h.SetMaximum(110)
	h.SetMinimum(0.0)
	h.Draw()

	f = TF1("fun","100*cos(x*3.14159/180.)",0,90)
	f.SetLineColor(1)
	f.SetLineStyle(2)
	f.Draw("same")

	# add graph
	gr = TGraph( len(x), x, y)
	gr.Draw("c")
	gr.SetLineWidth(4)
	c_Bethe = TColor.GetColor( '#BD265D' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )

	# add data points
	xp = array( 'f', [0, 20, 45, 70, 90] ) # angle in deg
	ncount = array ( 'f', [1950., 3698., 7428., 17453., 2433.] )
	itime  = array ( 'f', [1200., 2400., 7200., 31200., 6000.] )
	yp = array ( 'f', [0, 0, 0, 0, 0] )  # number of counts / time interval (seconds)
	exp = array( 'f', [0, 0, 0, 0, 0] )
	eyp = array( 'f', [0, 0, 0, 0, 0] )
	# normalise to 0 deg
	n = len(xp)
	xnorm = ncount[0]/itime[0]
	for i in range(0,n):
		yp[i] = ncount[i]/itime[i]/xnorm
		eyp[i] = m.sqrt( yp[i]*(1.0 - yp[i])/ncount[i] + yp[i]*(1.0 - yp[i])/ncount[0] )
		yp[i] *= 100.0
		eyp[i] *= 100.0
		print 'Datapoint ',i,' (',xp[i],' deg) :',yp[i],' +- ',eyp[i],' %'
	grp = TGraphErrors( n, xp, yp, exp, eyp )
	grp.SetMarkerStyle(20)
	grp.SetMarkerSize(0.8)
	c_p = TColor.GetColor( '#2222FF' )
	grp.SetLineColor(c_p)
	grp.SetLineWidth(2)
	grp.SetMarkerColor(c_p)
	grp.Draw("p")

	c.Update()
	c.Print("musim.pdf")



def main():
	phiLoop()

if __name__ == "__main__":
	main()
