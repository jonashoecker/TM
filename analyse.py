import math
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *
import atlasstyle

class DataPoint:
	def __init__ ( self, date, time, event, atime, adccount, siPMsignal, deadtime, temp, name ):
		self.date            = date 
		self.time            = time
		self.event           = float(event)
		self.atime           = float(atime)       # atime given in seconds
		self.adccount        = float(adccount)
		self.siPMsignal      = float(siPMsignal)
		self.deadtime        = float(deadtime)
		self.temp            = float(temp)
		self.name            = name

	def printContent( self ):
		print 'This object has time: ',self.time

def readdata ( fname1 ): 
	dataList = []
	f = open( fname1,'r' )
	if not f: 
		print "could not open file", f
		exit(1)
	flines = f.readlines()
	for line in flines:
		if line[0] == '2':
			l = line.split()
			dp =  DataPoint( l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8] )
			dataList.append( dp )

	return dataList

#def readdata ( fname2 ): 
	dataList = []
	f = open( fname2, 'r' )
	if not f: 
		print "could not open file", f
		exit(1)
	flines = f.readlines()
	for line in flines:
		if line[0] == '2':
			l = line.split()
			dp =  DataPoint( l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8] )
			dataList.append( dp )

	return dataList

def computeFrequency( dl ):
	freq = len(dl)/(dl[-1].atime - dl[0].atime)*1000
	return freq
	
def adcCountDistribution( dl1, dl2, nbins ): 
	
	c = TCanvas( 'c', 'c', 0, 0, 800, 600 )

	# define binning
	xmin = 0.0
	xmax = 700.0
	hist1 = TH1F( 'ADC counts', 'ADC counts', nbins, xmin, xmax )
	hist2 = TH1F( 'ADC counts', 'ADC counts', nbins, xmin, xmax )
	hist2.SetLineColor(2)
	hist2.SetLineStyle(2)
	hist1.SetLineWidth(2)
	hist2.SetLineWidth(2)


	# count number of events in each ADC count interval
	# and compute normalisation
	norm1 = 1.
	for dp in dl1: 
		hist1.Fill(dp.adccount)
		if dp.adccount < 120: norm1 += 1
	norm2 = 1.
	for dp in dl2: 
		hist2.Fill(dp.adccount)	
		if dp.adccount < 120: norm2 += 1


	# normalise histograms
	hist1.Scale( 1./norm1 )
	hist2.Scale( 1./norm2 )

	hist1.Draw('hist')
	hist2.Draw('histsame')
	c.Update()
	c.Print('ADCcounts.pdf')

def siPMsignalDistribution( dl, nbins ): 
	
	c = TCanvas( 'c', 'c', 0, 0, 800, 600 )

	# define binning
	xmin = 0
	xmax = 400
	hist = TH1F( 'siPMsignal', 'siPMsignal', nbins, xmin, xmax )

	
	# count number of events in each ADC count interval
	for dp in dl: 
		hist.Fill(dp.siPMsignal)

	hist.Draw('')
	c.Update()
	c.Print('SiPMSignal.pdf')




def eventDistribution( dl, interval, mergeBins ): 	# interval given in seconds 
	interval *= 1000                   				# convert to seconds
	t0 = dl[0].atime
	t1 = dl[-1].atime
	ninterval = int((t1 - t0) / interval + 1)
	intCount = [0.0] * ninterval
	

	# count number of events in each time interval
	for dp in dl: 
		iInterval = int((dp.atime - t0) / interval)
		intCount[iInterval] += 1.0

	# find maximum number of events in interval
	maxEvents = 0
	numEvents = 0
	for iC in intCount:
		numEvents += iC
		if iC > maxEvents : maxEvents = iC
	print 'Total number of Events:',int(numEvents)
	print 'Maximum number of events in interval:',int(maxEvents)

	# define binning
	bins = []
	for i in xrange(0,int((maxEvents)+2)/mergeBins):
		bins.append(i*mergeBins)
	print bins

	# fill histogram array (list)
	histArray = [0.0] * (int(maxEvents) + 2)
	for iC in intCount:	
		histArray[int(iC)] += 1.0


	print histArray

	# histogram
	n, bins, patches = plt.hist( x=intCount, bins=bins, color='#1589ff', alpha=0.8, rwidth=1)
	plt.grid( axis='y', alpha=0.75 )
	plt.xlabel('Number of counts in interval')
	plt.ylabel('Number of intervals')
	plt.title('')
	maxfreq = n.max()
	# Set a clean upper y-axis limit.
	plt.ylim( ymax=np.ceil( maxfreq*1.1 ) )
	plt.show()

def main():
	#fname2 = "GREEN-Vertical"
	fname1 = "CW_data-horizontal-24dec2019.txt"
	fname2 = "CW_data-vertical-24dec2019.txt"
	dataList1 = readdata( fname1 )
	dataList2 = readdata( fname2 )
	print 'Number of events per second in file ',fname1,':',computeFrequency( dataList1 ),' Hz (total number of events:', len(dataList1),')'
	#print 'Number of events per second in file ',fname2,':',computeFrequency( dataList2 ),' Hz (total number of events:', len(dataList2),')'
	timeInterval = 10 # in seconds
	mergeBins = 1
	#eventDistribution( dataList1, timeInterval, mergeBins )
	nbins = 50
	adcCountDistribution( dataList1, dataList2, nbins )
	#siPMsignalDistribution( dataList1, nbins )

if __name__ == "__main__":
	main()


