import math as m
from math import pi
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
from array import array
import rootStyle
 
z = -1.0 # muon charge
Z = 7. # atomic number (nitrogen)
A = 14. # atomic mass (nitrogen) g / m
#Z = 14.0 # atomic number (silicon)
#A = 28.0 # atomic mass (silicon) g / m
K = 0.307075 # MeV cm**2 / mol
me = 0.511 # mass of electron in MeV
mass = 105.7 # mass of muon
c = 299792458 # m / s
I = 12.0 * Z # eV
I *= 1e-6 # in MeV
delta = 0.0 # density effect at large distance
dx = 1,5e+6
rho_air = 1.293 # kg/m^3 at sea level
rho_air *= 1.e3/100**3
rho_limestone = 2720 # kg/m^3 at sea level
rho_limestone *= 1.e3/100**3


def rho_alt(altitude):
	return rho_air*m.exp(-altitude / 10400.0) # 10.4 in m


def dE(dx,beta,gamma):
	Wmax = 2.0 * me * beta**2 * gamma**2 / ( 1.0 + 2.0*gamma*me/mass + (me/mass)**2 )
	p1 = K * z**2 * Z / ( A * beta**2 )
	p2 = 0.5 * m.log( 2 * me * beta**2 * gamma**2 * Wmax / (I**2) )
	p3 = - beta**2
	p4 = - delta / 2.0
	res = - dx * p1 * ( p2 + p3 + p4 ) 
	return res

def plotGraph():
	dx = 0.5

	x, y = array( 'd' ), array( 'd' )
	p = 0.02
	while p < 100:
		betagamma = p*1.0e3/mass
		beta = betagamma/m.sqrt(1+betagamma**2)
		gamma = 1./m.sqrt(1-beta**2)
		print beta,gamma,beta*gamma,p
		dEdx = dE(dx, beta, gamma)/dx
		x.append(betagamma)
		y.append(-dEdx)
		p += 0.01

	c = TCanvas("c1","c1",0,0,600,500)
	c.SetLogx(); 
	c.SetGrid()

	gr = TGraph( len(x), x, y )
	gr.SetMaximum(10)
	gr.SetMinimum(0)
	gr.Draw("al")
	gr.SetLineWidth(2)
	c_Bethe = TColor.GetColor( '#1A52A9' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )

	gr.GetXaxis().SetTitle("Impulsion (GeV/c)")
	gr.GetYaxis().SetTitle("#frac{-dE}{dx#upoint#rho}  (MeV g^{-1} cm^{2})")
	c.Update()
	c.Print("c.pdf")

def plotGraph2():

	x, y = array( 'd' ), array ( 'd' )
	alt = 15000.0 # m
	dx = 10 # m
	p = 6 # GeV
	while alt >= 0:
		rho = rho_alt(alt)
		eloss = calcEnergyLoss(p, dx, rho)
		print alt,rho/(1.e3/100**3),eloss,p
		p += eloss
		alt -= dx
		x.append(p)
		y.append(eloss)

	c = TCanvas("c","c",0,0,600,500)
	c.SetLogx();
	c.SetGrid()

	gr = TGraph( len(x), x, y)
	gr.SetMaximum(7)
	gr.SetMinimum(0)
	gr.Draw("al")
	gr.SetLineWidth(4)
	c_Bethe = TColor.GetColor( '#1A52A9' )
	gr.SetLineColor(c_Bethe)
	gr.GetYaxis().SetLabelOffset( 0.01 )

	gr.GetXaxis().SetTitle("dx")
	gr.GetYaxis().SetTitle("dE")
	c.Update()
	c.Print("e.png")



def calcEnergyLoss(p, dx, rho):

	betagamma = p/mass
	beta = betagamma/m.sqrt(1+betagamma**2)
	gamma = 1./m.sqrt(1-beta**2)
	energyLoss = dE(dx, beta, gamma)*rho/1000
	return energyLoss

def energylossAir():
	alt = 15000.0 # m
	dx = 10 # m
	p = 6 # GeV
	while alt >= 0:
		rho = rho_alt(alt)
		eloss = calcEnergyLoss(p, dx, rho)
		print alt,rho/(1.e3/100**3),eloss,p
		p += eloss
		alt -= dx

def energylossRock():
	thickness = 10 # m
	p = 4 # GeV
	x = thickness
	dx = 0.1
	while x >= 0:
		eloss = calcEnergyLoss(p, dx, rho_limestone)
		print dx,rho_limestone/(1.e3/100**3),p,eloss
		p += eloss
		x -= dx



def main():

	plotGraph()
	#plotGraph2()
	#energylossAir()
	#energylossRock()
	#print calcEnergyLoss(1.3585, 1, rho_limestone)


if __name__ == "__main__":
	main()




