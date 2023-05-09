#! /usr/bin/env python

import sys
import os
from array import *
import numpy as np
import ROOT
import yaml
import math

ROOT.gROOT.SetBatch(True)

class AliHFInvMassFitter:
	ETypeOfBkg = {"kExpo": 0,
				  "kLin": 1,
				  "kPol2": 2,
				  "kNoBk": 3,
				  "kPow": 4,
				  "kPowEx": 5}
	ETypeOfSgn = {"kGaus": 0,
				  "k2Gaus": 1,
				  "k2GausSigmaRatioPar":2}

	def __init__(self, histoToFit, minvalue, maxvalue, fittypeb=ETypeOfBkg["kExpo"], fittypes=ETypeOfSgn["kGaus"]):
		self.fMinMass = 1.72
		self.fMaxMass = 2.05
		self.fTypeOfFit4Bkg = fittypeb
		self.fTypeOfFit4Sgn = fittypes

		self.fPolDegreeBkg = 4
		self.fCurPolDegreeBkg = -1
		self.fMassParticle = 1.864
		self.fMass = 1.865
		self.fMassErr = 0
		self.fSigmaSgn = 0.012
		self.fSigmaSgnErr = 0
		self.fSigmaSgn2Gaus = 0.012
		self.fFixedMean = False
		self.fBoundMean = False
		self.fMassLowerLim = 0
		self.fMassUpperLim = 0
		self.fFixedSigma = False
		self.fBoundSigma = False
		self.fSigmaVar = 0.012
		self.fParSig = 0.1
		self.fFixedSigma2Gaus = False
		self.fFixedRawYield = -1
		self.fFrac2Gaus = 0.2
		self.fFixedFrac2Gaus = False
		self.fRatio2GausSigma = 0
		self.fFixedRatio2GausSigma = False
		self.fNParsSig = 3
		self.fNParsBkg = 2
		self.fOnlySideBands = False
		self.fNSigma4SideBands = 4
		self.fCheckSignalCountsAfterFirstFit = True
		self.fFitOption = "L,E"
		self.fRawYield = 0
		self.fRawYieldErr = 0
		self.fSigFunc = 0x0
		self.fBkgFuncSb = 0x0
		self.fBkgFunc = 0x0
		self.fBkgFuncRefit = 0x0
		self.fReflections = False
		self.fNParsRfl = 0
		self.fRflOverSig = 0
		self.fFixRflOverSig = False
		self.fHistoTemplRfl = 0x0
		self.fSmoothRfl = False
		self.fRawYieldHelp = 0
		self.fRflFunc = 0x0
		self.fBkRFunc = 0x0
		self.fSecondPeak = False
		self.fNParsSec = 0
		self.fSecMass = -999
		self.fSecWidth = 9999
		self.fFixSecMass = False
		self.fFixSecWidth = False
		self.fSecFunc = 0x0
		self.fTotFunc = 0x0
		
		self.fHistoInvMass = histoToFit.Clone("fHistoInvMass")
		self.fHistoInvMass.SetDirectory(0)
		self.SetNumberOfParams()


	#------------------------------------------#
	def SetNumberOfParams(self):
		switchNParsBkg = {0: 2,
				  1: 2,
				  2: 3,
				  3: 1,
				  4: 2,
				  5: 2,
				  6: self.fPolDegreeBkg+1
				 }
		switchNParsSig = {0: 3,
					 1: 5,
					 2: 5
				 }

		self.fNParsBkg = switchNParsBkg.get(self.fTypeOfFit4Bkg,"Error in computing fNParsBkg: check fTypeOfFit4Bkg")
		self.fNParsSig = switchNParsSig.get(self.fTypeOfFit4Sgn,"Error in computing fNParsSig: check fTypeOfFit4Sgn")


		if(self.fReflections):
			self.fNParsRfl = 1
		else:
			self.fNParsRfl = 0
			
		if(self.fSecondPeak):
			self.fNParsSec = 3
		else:
			self.fNParsSec = 0

	#------------------------------------------#
	def MassFitter(self, draw):
		# Main function to fit the invariant mass distribution
		# returns 0 if fit fails 
		# returns 1 if the fit succeeds
		# returns 2 if there is no signal and the fit is performed only w/bkg

		ROOT.TVirtualFitter.SetDefaultFitter("Minuit")
		#ROOT.TFitter.SetPrecision(1e-2)

		integralHisto = self.fHistoInvMass.Integral(self.fHistoInvMass.FindBin(self.fMinMass),self.fHistoInvMass.FindBin(self.fMaxMass),"width")
		self.fOnlySideBands = True
		self.fBkgFuncSb = self.CreateBackgroundFitFunction("funcbkgsb",integralHisto)
		status = -1
		print("-- First fit with only background on the sides bands - exclusion region = {} sigma --".format(self.fNSigma4SideBands))
		if self.fTypeOfFit4Bkg == 6:
			if self.PrepareHighPolFit(self.fBkgFuncSb):
				status = 0
		else:
			fitoption = "R,S,{},+,0".format(self.fFitOption)
			status = int(self.fHistoInvMass.Fit("funcbkgsb",fitoption ))
		self.fBkgFuncSb.SetLineColor(ROOT.kGray+1) # kGray + 1
		if status != 0:
			print("   --=> Failed first fit with only background, minuit status = {}".format(status))
			return 0
	
		fOnlySideBands = False
		if not self.fBkgFunc:
			self.fBkgFunc = self.CreateBackgroundFitFunction("funcbkg",integralHisto)
			for ipar in range(0, self.fNParsBkg):
				self.fBkgFunc.SetParameter(ipar,self.fBkgFuncSb.GetParameter(ipar))
			self.fBkgFunc.SetLineColor(ROOT.kGray+1) # kGray + 1
		
		print("--- Estimate signal counts in the peak region ---")
		estimSignal=self.CheckForSignal(self.fMass,self.fSigmaSgn)
		doFinalFit=True
		if(self.fCheckSignalCountsAfterFirstFit and estimSignal<0):
			if(draw): 
				self.DrawFit()
			estimSignal=0.
			doFinalFit=False
			print("Abandon fit: no signal counts after first fit")
			
		fRawYieldHelp=estimSignal # needed for reflection normalization
		if not self.fBkgFuncRefit:
			self.fBkgFuncRefit = self.CreateBackgroundFitFunction("funcbkgrefit",integralHisto)
			for ipar in range(0, self.fNParsBkg):
				self.fBkgFuncRefit.SetParameter(ipar,self.fBkgFunc.GetParameter(ipar))
		self.fBkgFuncRefit.SetLineColor(2)
		self.fSigFunc = self.CreateSignalFitFunction("fsigfit",estimSignal)
		if self.fSecondPeak:
			print("   --=> Final fit includes a second inv. mass peak")
			estimSec=self.CheckForSignal(self.fSecMass,self.fSecWidth)
			self.fSecFunc = self.CreateSecondPeakFunction("fsecpeak",estimSec)
		if self.fReflections:
			print("   --=> Final fit includes reflections")
			self.fRflFunc = self.CreateReflectionFunction("freflect")
			self.fBkRFunc = self.CreateBackgroundPlusReflectionFunction("fbkgrfl")
		self.fTotFunc = self.CreateTotalFitFunction("funcmass")  

		if(doFinalFit):
			print("--- Final fit with signal+background on the full range ---")
			fitoption = "R,{},,0".format(self.fFitOption)
			status=int(self.fHistoInvMass.Fit("funcmass",fitoption))
			if (status != 0):
				print("   --. Failed fit with signal+background, minuit status = {}".format(status))
				#return 0
	
		for ipar in range (0, self.fNParsBkg):
			self.fBkgFuncRefit.SetParameter(ipar,self.fTotFunc.GetParameter(ipar))
			self.fBkgFuncRefit.SetParError(ipar,self.fTotFunc.GetParError(ipar))
		for ipar in range(0, self.fNParsSig):
			self.fSigFunc.SetParameter(ipar,self.fTotFunc.GetParameter(ipar+self.fNParsBkg))
			self.fSigFunc.SetParError(ipar,self.fTotFunc.GetParError(ipar+self.fNParsBkg))
		if self.fSecondPeak:
			for ipar in range(0, self.fNParsSec):
				self.fSecFunc.SetParameter(ipar,self.fTotFunc.GetParameter(ipar+self.fNParsBkg+self.fNParsSig))
				self.fSecFunc.SetParError(ipar,self.fTotFunc.GetParError(ipar+self.fNParsBkg+self.fNParsSig))				
			self.fSecFunc.SetLineColor(ROOT.kMagenta+1)
			self.fSecFunc.SetLineStyle(3)
		if self.fReflections:
			for ipar in range(0, self.fNParsRfl):
				self.fRflFunc.SetParameter(ipar,self.fTotFunc.GetParameter(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec))
				self.fRflFunc.SetParError(ipar,self.fTotFunc.GetParError(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec))
				self.fBkRFunc.SetParameter(ipar+self.fNParsBkg,self.fTotFunc.GetParameter(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec))
				self.fBkRFunc.SetParError(ipar+self.fNParsBkg,self.fTotFunc.GetParError(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec))
			for ipar in range(0, self.fNParsBkg):
				self.fBkRFunc.SetParameter(ipar,self.fTotFunc.GetParameter(ipar))
				self.fBkRFunc.SetParError(ipar,self.fTotFunc.GetParError(ipar))
			self.fRflFunc.SetLineColor(ROOT.kGreen+1)
			self.fBkRFunc.SetLineColor(ROOT.kRed+1)
			self.fBkRFunc.SetLineStyle(7)

		self.fMass=self.fSigFunc.GetParameter(1)
		print("-------------")
		print(self.fMass)
		self.fMassErr=self.fSigFunc.GetParError(1)
		self.fSigmaSgn=self.fSigFunc.GetParameter(2)
		self.fSigmaSgnErr=self.fSigFunc.GetParError(2)
		self.fTotFunc.SetLineColor(4)
		self.fRawYield=self.fTotFunc.GetParameter(self.fNParsBkg)/self.fHistoInvMass.GetBinWidth(1)
		self.fRawYieldErr=self.fTotFunc.GetParError(self.fNParsBkg)/self.fHistoInvMass.GetBinWidth(1)
		self.fRawYieldHelp=self.fRawYield
		if(draw): 
			self.DrawFit()
		if(doFinalFit): 
			return 1
		else: 
			return 2

	# ------------------------------------------#
	def DrawFit(self):
		c0 = ROOT.TCanvas("c0")
		self.DrawHere(c0)

	def DrawHere(self, c, nsigma=3, writeFitInfo=0):
		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetCanvasColor(0)
		ROOT.gStyle.SetFrameFillColor(0)
		c.cd()
		self.fHistoInvMass.GetXaxis().SetRangeUser(self.fMinMass,self.fMaxMass)
		self.fHistoInvMass.SetMarkerStyle(20)
		self.fHistoInvMass.SetMinimum(0.)
		self.fHistoInvMass.SetMaximum(400)
		self.fHistoInvMass.Draw("PE")
		if self.fBkgFunc: 
			self.fBkgFunc.Draw("same")
		if self.fBkgFuncRefit: 
			self.fBkgFuncRefit.Draw("same")
		if self.fBkRFunc: 
			self.fBkRFunc.Draw("same")
		if self.fRflFunc: 
			self.fRflFunc.Draw("same")
		if self.fSecFunc: 
			self.fSecFunc.Draw("same")
		if self.fTotFunc: 
			self.fTotFunc.Draw("same")

		if writeFitInfo > 0:
			pinfos = ROOT.TPaveText(0.12,0.65,0.47,0.89,"NDC")
			pinfom = ROOT.TPaveText(0.6,0.7,1.,.87,"NDC")
			pinfos.SetBorderSize(0)
			pinfos.SetFillStyle(0)
			pinfom.SetBorderSize(0)
			pinfom.SetFillStyle(0)
			if(self.fTotFunc):
				pinfom.SetTextColor(ROOT.kBlue)
				for ipar in range(1, self.fNParsSig):
					str = "{} = {:.3f} +/- {:.3f}".format(self.fTotFunc.GetParName(ipar+self.fNParsBkg),self.fTotFunc.GetParameter(ipar+self.fNParsBkg),self.fTotFunc.GetParError(ipar+self.fNParsBkg))  
					if (self.fTotFunc.GetParameter(self.fNParsBkg+self.fNParsSig-2)<0.2):
						str = "{} = {:.3f} +/- {:.3f}".format(self.fTotFunc.GetParName(ipar+self.fNParsBkg),self.fTotFunc.GetParameter(ipar+self.fNParsBkg)*1000,self.fTotFunc.GetParError(ipar+self.fNParsBkg)*1000)    
					pinfom.AddText(str)

				if writeFitInfo >= 1:
					pinfom.Draw()

				bkg = ROOT.Double()
				errbkg = ROOT.Double()
				self.Background(nsigma,bkg,errbkg)
				signif = ROOT.Double()
				errsignif = ROOT.Double()
				self.Significance(nsigma,signif,errsignif)

				pinfos.AddText("S = {:.0f} \u00B1 {:.0f}".format(self.fRawYield,self.fRawYieldErr))
				pinfos.AddText("B ({:.0f}\u03C3) = {:.0f} \u00B1 {:.0f}".format(nsigma,bkg,errbkg))
				if bkg==0:
					pinfos.AddText("S/B ({:.0f}\u03C3) = 0 ".format(nsigma))
				else:
					pinfos.AddText("S/B ({:.0f}\u03C3) = {:.4f} ".format(nsigma,self.fRawYield/bkg))					
				if(self.fRflFunc):  pinfos.AddText("Refl/Sig =  {:.3f} \u00B1 {:.3f} ".format(self.fRflFunc.GetParameter(0),self.fRflFunc.GetParError(0)))
				pinfos.AddText("Signif ({:.0f}\u03C3) = {:.1f} \u00B1 {:.1f} ".format(nsigma,signif,errsignif))
				if writeFitInfo>=2: pinfos.Draw()

		c.Update()
		c.SaveAs("testing.png")
		return


	#------------------------------------------#
	def CheckForSignal(self, mean, sigma):
		# Checks if there are signal counts above the background
		# in the invariant mass region of the peak

		minForSig=mean-4.*sigma
		maxForSig=mean+4.*sigma
		binForMinSig=self.fHistoInvMass.FindBin(minForSig)
		binForMaxSig=self.fHistoInvMass.FindBin(maxForSig)
		sumpeak=0.
		sumback=0.
		self.fBkgFunc.Print()
		for ibin in range(binForMinSig, binForMaxSig+1):
			sumpeak = sumpeak + self.fHistoInvMass.GetBinContent(ibin)
			sumback = sumback + self.fBkgFunc.Eval(self.fHistoInvMass.GetBinCenter(ibin))
		diffUnderPeak=(sumpeak-sumback)
		print("   --=> IntegralUnderHisto={}  IntegralUnderBkgFunc={}   EstimatedSignal={}".format(sumpeak,sumback,diffUnderPeak))
		if (diffUnderPeak/math.sqrt(sumpeak)) < 1:
			print("   --=> (Tot-Bkg)/sqrt(Tot)={} --=> Likely no signal".format(diffUnderPeak/math.sqrt(sumpeak)))
			return -1
		
		return diffUnderPeak*self.fHistoInvMass.GetBinWidth(1)		


	#------------------------------------------#
	def CreateBackgroundFitFunction(self, fname, integral):
		# Creates the background fit function
		self.SetNumberOfParams()
		funcbkg = ROOT.TF1(fname,self.FitFunction4Bkg, self.fMinMass, self.fMaxMass, self.fNParsBkg)
		if self.fTypeOfFit4Bkg == 0:
			funcbkg.SetParNames("BkgInt","Slope")
			funcbkg.SetParameters(integral,-2.)
		elif self.fTypeOfFit4Bkg == 1:
			funcbkg.SetParNames("BkgInt","Slope")
			funcbkg.SetParameters(integral,-100.)
		elif self.fTypeOfFit4Bkg == 2:
			funcbkg.SetParNames("BkgInt","Coef1","Coef2")
			funcbkg.SetParameters(integral,-10.,5)
		elif self.fTypeOfFit4Bkg == 3:
			funcbkg.SetParNames("Const")
			funcbkg.SetParameter(0,0.)
			funcbkg.FixParameter(0,0.)
		elif self.fTypeOfFit4Bkg == 4:
			funcbkg.SetParNames("BkgInt","Coef1")
			funcbkg.SetParameters(integral,0.5)
		elif self.fTypeOfFit4Bkg == 5:
			funcbkg.SetParNames("Coef1","Coef2")
			funcbkg.SetParameters(-10.,5.)
		elif self.fTypeOfFit4Bkg == 6:		
			for j in range(0,self.fNParsBkg):
				funcbkg.SetParName(j, "Coef{}".format(j))
				funcbkg.SetParameter(j,0)
				funcbkg.SetParameter(0,integral)
		else:
			print("Wrong choice of fTypeOfFit4Bkg {}",self.fTypeOfFit4Bkg)
			funcbkg.Delete()
			return 0x0
		funcbkg.SetLineColor(ROOT.kBlue+3) # kBlue+3

		return funcbkg
		
		
	#------------------------------------------#
	def CreateSecondPeakFunction(self, fname, integsig):
		# Creates a function for a gaussian peak in the background fit function
		# Can be used e.g. to include the D+=>KKpi peak in the D_s inv. mass fit

		funcsec =  ROOT.TF1(fname,self.FitFunction4SecPeak,self.fMinMass,self.fMaxMass,3)
		funcsec.SetParameter(0,integsig)
		funcsec.SetParameter(1,self.fSecMass)
		if self.fFixSecMass: funcsec.FixParameter(1,self.fSecMass)
		funcsec.SetParameter(2,self.fSecWidth)
		if self.fFixSecWidth: funcsec.FixParameter(2,self.fSecWidth)
		funcsec.SetParNames("SecPeakInt","SecMean","SecSigma")
		return funcsec


	#------------------------------------------#
	def CreateReflectionFunction(self, fname):
		# Creates a function for reflections contribution
		# in the D0.Kpi inv. mass distribution
		funcrfl =  ROOT.TF1(fname,self.FitFunction4Refl,self.fMinMass,self.fMaxMass,1)
		funcrfl.SetParameter(0,self.fRflOverSig)
		funcrfl.SetParLimits(0,0.,1.)
		if self.fFixRflOverSig: 
			funcrfl.FixParameter(0,self.fRflOverSig)
		funcrfl.SetParNames("ReflOverS")
		return funcrfl

		
	#------------------------------------------#
	def CreateBackgroundPlusReflectionFunction(self, fname):
		# Creates the function with sum of background and reflections
		#
		self.SetNumberOfParams()
		totParams = self.fNParsBkg+self.fNParsRfl
		fbr = ROOT.TF1(fname,self.FitFunction4BkgAndRefl,self.fMinMass,self.fMaxMass,totParams)
		for ipar in range(0,self.fNParsBkg):
				fbr.SetParameter(ipar,self.fBkgFunc.GetParameter(ipar))
				fbr.SetParName(ipar,self.fBkgFunc.GetParName(ipar))
		for ipar in range(0,self.fNParsRfl):
				fbr.SetParameter(ipar+self.fNParsBkg,self.fRflFunc.GetParameter(ipar))
				fbr.SetParName(ipar+self.fNParsBkg,self.fRflFunc.GetParName(ipar))
				# par limits not set because this function is not used for fitting but only for drawing
		return fbr


	#------------------------------------------#
	def CreateSignalFitFunction(self, fname, integsig):
		# Creates the fit function for the signal peak

		self.SetNumberOfParams()
		funcsig = ROOT.TF1(fname,self.FitFunction4Sgn,self.fMinMass,self.fMaxMass,self.fNParsSig)
		if self.fTypeOfFit4Sgn==self.ETypeOfSgn["kGaus"]:
			funcsig.SetParameter(0,integsig)
			if self.fFixedRawYield>-0.1: funcsig.FixParameter(0,self.fFixedRawYield)
			funcsig.SetParameter(1,self.fMass)
			if self.fFixedMean: funcsig.FixParameter(1,self.fMass)
			if self.fBoundMean: funcsig.SetParLimits(1,self.fMassLowerLim, self.fMassUpperLim)
			funcsig.SetParameter(2,self.fSigmaSgn)
			if self.fFixedSigma: funcsig.FixParameter(2,self.fSigmaSgn)
			if self.fBoundSigma: funcsig.SetParLimits(2,self.fSigmaVar*(1-self.fParSig), self.fSigmaVar*(1+self.fParSig))
			funcsig.SetParNames("SgnInt","Mean","Sigma")
		elif self.fTypeOfFit4Sgn==self.ETypeOfSgn["k2Gaus"]:
			funcsig.SetParameter(0,integsig)
			if self.fFixedRawYield>-0.1: funcsig.FixParameter(0,self.fFixedRawYield)
			funcsig.SetParameter(1,self.fMass)
			if self.fFixedMean: funcsig.FixParameter(1,self.fMass)
			funcsig.SetParameter(2,self.fSigmaSgn)
			funcsig.SetParLimits(2,0.004,0.05)
			if self.fFixedSigma: funcsig.FixParameter(2,self.fSigmaSgn)
			funcsig.SetParameter(3,self.fFrac2Gaus)
			if self.fFixedFrac2Gaus: 
				funcsig.FixParameter(3,self.fFrac2Gaus)
			else: 
				funcsig.SetParLimits(3,0.,1.)
			funcsig.SetParameter(4,self.fSigmaSgn2Gaus)
			if self.fFixedSigma2Gaus: 
				funcsig.FixParameter(4,self.fSigmaSgn2Gaus)
			else: 
				funcsig.SetParLimits(4,0.004,0.05)
			funcsig.SetParNames("SgnInt","Mean","Sigma1","Frac","Sigma2")
		elif self.fTypeOfFit4Sgn==self.ETypeOfSgn["k2GausSigmaRatioPar"]:
			funcsig.SetParameter(0,integsig)
			if self.fFixedRawYield>-0.1: funcsig.FixParameter(0,self.fFixedRawYield)
			funcsig.SetParameter(1,self.fMass)
			if self.fFixedMean: funcsig.FixParameter(1,self.fMass)
			funcsig.SetParameter(2,self.fSigmaSgn)
			funcsig.SetParLimits(2,0.004,0.05)
			if self.fFixedSigma: funcsig.FixParameter(2,self.fSigmaSgn)
			funcsig.SetParameter(3,self.fFrac2Gaus)
			if self.fFixedFrac2Gaus: 
				funcsig.FixParameter(3,self.fFrac2Gaus)
			else:
				funcsig.SetParLimits(3,0.,1.)
			funcsig.SetParameter(4,self.fSigmaSgn2Gaus)
			if self.fFixedRatio2GausSigma:
				funcsig.FixParameter(4,self.fRatio2GausSigma)
			else:
				funcsig.SetParLimits(4,0.,20.)
			funcsig.SetParNames("SgnInt","Mean","Sigma1","Frac","RatioSigma12")
	
		return funcsig

		
	#------------------------------------------#
	def CreateTotalFitFunction(self, fname):
		# Creates the total fit fucntion (signal+background+possible second peak)
		#

		self.SetNumberOfParams()
		totParams=self.fNParsBkg+self.fNParsRfl+self.fNParsSec+self.fNParsSig
		ftot=ROOT.TF1(fname,self.FitFunction4Mass,self.fMinMass,self.fMaxMass,totParams)
		parmin = ROOT.Double()
		parmax = ROOT.Double()
		for ipar in range(0, self.fNParsBkg):
			ftot.SetParameter(ipar,round(self.fBkgFunc.GetParameter(ipar),3))
			ftot.SetParName(ipar,self.fBkgFunc.GetParName(ipar))
			print("Background parameter {}".format(round(self.fBkgFunc.GetParameter(ipar),3)))
		#ftot.SetParLimits(0,0,round(self.fBkgFunc.GetParameter(0),3)*2)
		for ipar in range(0, self.fNParsSig):
			ftot.SetParameter(ipar+self.fNParsBkg,round(self.fSigFunc.GetParameter(ipar),3))
			ftot.SetParName(ipar+self.fNParsBkg,self.fSigFunc.GetParName(ipar))
			print("Signal parameter {}".format(round(self.fSigFunc.GetParameter(ipar),3)))
			#self.fSigFunc.GetParLimits(ipar,parmin,parmax)
		ftot.SetParLimits(2+self.fNParsBkg,0.004,0.05)
		ftot.SetParLimits(1+self.fNParsBkg,1.84,1.9)
			#print("par limits {} to {} ".format(parmin,parmax))
		if(self.fSecondPeak and self.fSecFunc):
			for ipar in range(0, self.fNParsSec):
				ftot.SetParameter(ipar+self.fNParsBkg+self.fNParsSig,self.fSecFunc.GetParameter(ipar))
				ftot.SetParName(ipar+self.fNParsBkg+self.fNParsSig,self.fSecFunc.GetParName(ipar))
				self.fSecFunc.GetParLimits(ipar,parmin,parmax)
				ftot.SetParLimits(ipar+self.fNParsBkg+self.fNParsSig,parmin,parmax)
		if(self.fReflections and self.fRflFunc):
			for ipar in range(0, self.fNParsRfl):
				ftot.SetParameter(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec,self.fRflFunc.GetParameter(ipar))
				ftot.SetParName(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec,self.fRflFunc.GetParName(ipar))
				self.fRflFunc.GetParLimits(ipar,parmin,parmax)
				ftot.SetParLimits(ipar+self.fNParsBkg+self.fNParsSig+self.fNParsSec,parmin,parmax)
		return ftot
		

	#------------------------------------------#	
	def FitFunction4Bkg(self, x, par):
		# Fit function for the background

		maxDeltaM = self.fNSigma4SideBands * self.fSigmaSgn
		if(self.fOnlySideBands and (abs(x[0] - self.fMass) < maxDeltaM)):
			ROOT.TF1.RejectPoint()
			return 0
		if(self.fOnlySideBands and self.fSecondPeak and (abs(x[0] - self.fSecMass) < (self.fNSigma4SideBands * self.fSecWidth))):
			ROOT.TF1.RejectPoint()
			return 0
		total = 0
		mpi = 139.57039
		if self.fTypeOfFit4Bkg == 0:
			#exponential
			#exponential = A*exp(B*x) => integral(exponential)=A/B*exp(B*x)](min,max)
			#=> A = B*integral/(exp(B*max)-exp(B*min)) where integral can be written
			#as integralTot- integralGaus (=par [2])
			#Par:
			# * [0] = integralBkg
			# * [1] = B
			#exponential = [1]*[0]/(exp([1]*max)-exp([1]*min))*exp([1]*x)
			total = par[0]*par[1]/(math.exp(par[1]*self.fMaxMass)-math.exp(par[1]*self.fMinMass))*math.exp(par[1]*x[0])
		elif self.fTypeOfFit4Bkg == 1:
			#linear
			#y=a+b*x => integral = a(max-min)+1/2*b*(max^2-min^2) => a = (integral-1/2*b*(max^2-min^2))/(max-min)=integral/(max-min)-1/2*b*(max+min)
			# * [0] = integralBkg
			# * [1] = b
			total = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))
		elif self.fTypeOfFit4Bkg == 2:
			#parabola
			#y=a+b*x+c*x**2 => integral = a(max-min) + 1/2*b*(max^2-min^2) +
			#+ 1/3*c*(max^3-min^3) =>
			#a = (integral-1/2*b*(max^2-min^2)-1/3*c*(max^3-min^3))/(max-min)
			# * [0] = integralBkg
			# * [1] = b
			# * [2] = c
			total = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))+par[2]*(x[0]*x[0]-1/3.*(self.fMaxMass*self.fMaxMass*self.fMaxMass-self.fMinMass*self.fMinMass*self.fMinMass)/(self.fMaxMass-self.fMinMass))
		elif self.fTypeOfFit4Bkg == 3:
			total = par[0]
		elif self.fTypeOfFit4Bkg == 4:
			#power function
			#y=a(x-m_pi)^b => integral = a/(b+1)*((max-m_pi)^(b+1)-(min-m_pi)^(b+1))
			#
			#a = integral*(b+1)/((max-m_pi)^(b+1)-(min-m_pi)^(b+1))
			# * [0] = integralBkg
			# * [1] = b
			# a(power function) = [0]*([1]+1)/((max-m_pi)^([1]+1)-(min-m_pi)^([1]+1))*(x-m_pi)^[1]
			total = par[0]*(par[1]+1.)/(math.pow(self.fMaxMass-mpi,par[1]+1.)-math.pow(self.fMinMass-mpi,par[1]+1.))*math.pow(x[0]-mpi,par[1])
		elif self.fTypeOfFit4Bkg == 5:
			#power function wit exponential
			#y=a*Sqrt(x-m_pi)*exp(-b*(x-m_pi))
			total = par[0]*math.sqrt(x[0] - mpi)*math.exp(-1*par[1]*(x[0]-mpi))
		elif self.fTypeOfFit4Bkg == 6:
			#     # pol 3, following convention for pol 2
			#     #y=a+b*x+c*x**2+d*x**3 => integral = a(max-min) + 1/2*b*(max^2-min^2) +
			#     #+ 1/3*c*(max^3-min^3) + 1/4 d * (max^4-min^4) =>
			#     #a = (integral-1/2*b*(max^2-min^2)-1/3*c*(max^3-min^3) - 1/4 d * (max^4-min^4) )/(max-min)
			#     # * [0] = integralBkg
			#     # * [1] = b
			#     # * [2] = c
			#     # * [3] = d
			for i in range(1, self.fPolDegreeBkg + 1):
				total = total + par[i]*math.pow(x[0]-self.fMassParticle,i)/math.factorial(i)
		else:
			print("fTypeOfFit4Bkg out of range")
			return 0

		return total


	#------------------------------------------#
	def FitFunction4Sgn(self, x, par):
		# Fit function for the signal

		# AliInfo("Signal function set to: Gaussian")
		sigval = 0
		g1 = 0
		g2 = 0
		if self.fTypeOfFit4Sgn == 0:
			#gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
			#Par:
			# * [0] = integralSgn
			# * [1] = mean
			# * [2] = sigma
			sigval=par[0]/math.sqrt(2*math.pi)/par[2]*np.exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2])
		elif self.fTypeOfFit4Sgn == 1:
			#double gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
			#Par:
			# * [0] = integralSgn
			# * [1] = mean
			# * [2] = sigma1
			# * [3] = 2nd gaussian ratio
			# * [4] = deltaSigma
			g1=(1.-par[3])/math.sqrt(2.*math.pi)/par[2]*math.exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])
			g2=par[3]/math.sqrt(2.*math.pi)/par[4]*math.exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[4]/par[4])
			sigval=par[0]*(g1+g2)
		elif self.fTypeOfFit4Sgn == 2:
			#double gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
			#Par:
			# * [0] = integralSgn
			# * [1] = mean
			# * [2] = sigma1
			# * [3] = 2nd gaussian ratio
			# * [4] = ratio sigma12
			g1=(1.-par[3])/math.sqrt(2.*math.pi)/par[2]*math.exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])
			g2=par[3]/math.sqrt(2.*math.pi)/(par[4]*par[2])*math.exp(-(x[0]-par[1])*(x[0]-par[1])/2./(par[4]*par[2])/(par[4]*par[2]))
			sigval=par[0]*(g1+g2)

		self.fRawYieldHelp = par[0]/self.fHistoInvMass.GetBinWidth(1)
		return sigval


	#------------------------------------------#
	def FitFunction4Refl(self, x, par):
		# Fit function for reflections
		# D0=>Kpi decays with swapped mass assignment to pion and kaon decay tracks
		if not self.fHistoTemplRfl: return 0
		
		xbin =self.fHistoTemplRfl.FindBin(x[0])
		value=self.fHistoTemplRfl.GetBinContent(xbin)
		binmin=self.fHistoTemplRfl.FindBin(self.fMinMass*1.00001)
		binmax=self.fHistoTemplRfl.FindBin(self.fMaxMass*0.99999)
		norm=self.fHistoTemplRfl.Integral(binmin,binmax)*self.fHistoTemplRfl.GetBinWidth(xbin)
		if (abs(value)<1.e-14 and self.fSmoothRfl):
			# very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
			value = value + self.fHistoTemplRfl.GetBinContent(xbin-1)+self.fHistoTemplRfl.GetBinContent(xbin+1)
			value = value/3.

		return par[0]*value/norm*self.fRawYieldHelp*self.fHistoInvMass.GetBinWidth(1)


	#------------------------------------------#
	def FitFunction4BkgAndRefl(self, x, par):
		 # Fit function with the sum of background and reflections
		 #
		bkg = 0
		if self.fTypeOfFit4Bkg == 0:
			bkg = par[0]*par[1]/(np.exp(par[1]*self.fMaxMass)-np.exp(par[1]*self.fMinMass))*np.exp(par[1]*x[0])
		elif self.fTypeOfFit4Bkg == 1:
			bkg = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))
		elif self.fTypeOfFit4Bkg == 2:
			bkg = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))+par[2]*(x[0]*x[0]-1/3.*(self.fMaxMass*self.fMaxMass*self.fMaxMass-self.fMinMass*self.fMinMass*self.fMinMass)/(self.fMaxMass-self.fMinMass))
		elif self.fTypeOfFit4Bkg == 3:
			bkg = par[0]
		elif self.fTypeOfFit4Bkg == 4:
			bkg = par[0]*(par[1]+1.)/(math.pow(self.fMaxMass-mpi,par[1]+1.)-math.pow(self.fMinMass-mpi,par[1]+1.))*math.pow(x[0]-mpi,par[1])
		elif self.fTypeOfFit4Bkg == 5:
			bkg = par[0]*math.sqrt(x[0] - mpi)*np.exp(-1*par[1]*(x[0]-mpi))
		elif self.fTypeOfFit4Bkg == 6:
			for i in range(1, self.fPolDegreeBkg + 1):
				bkg = bkg + par[i]*math.pow(x[0]-self.fMassParticle,i)/math.factorial(i)
		else:
			print("fTypeOfFit4Bkg out of range")
			return 0
		 
		refl=0
		if self.fReflections: 
			if not self.fHistoTemplRfl:
				return 0

			binx=self.fHistoTemplRfl.FindBin(x[0])
			value=self.fHistoTemplRfl.GetBinContent(binx)
			binmin=self.fHistoTemplRfl.FindBin(self.fMinMass*1.00001)
			binmax=self.fHistoTemplRfl.FindBin(self.fMaxMass*0.99999)
			norm=self.fHistoTemplRfl.Integral(binmin,binmax)*self.fHistoTemplRfl.GetBinWidth(binx)
			if(abs(value)<1.e-14 and self.fSmoothRfl):
				#very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
				value = value + self.fHistoTemplRfl.GetBinContent(binx-1)+self.fHistoTemplRfl.GetBinContent(binx+1)
				value = value/3

			refl = par[0]*value/norm*self.fRawYieldHelp*self.fHistoInvMass.GetBinWidth(1)		
	
		return bkg+refl


	#------------------------------------------#
	def FitFunction4SecPeak(self, x, par):
		#/ Fit function for a second gaussian peak
		#/ To be used, e.g., for D+=>KKpi in the Ds mass spectrum

		#gaussian = A/(sigma*sqrt(2*pi))*exp(-(x-mean)^2/2/sigma^2)
		#Par:
			# * [0] = integralSgn
			# * [1] = mean
			# * [2] = sigma
		secgaval=par[0]/math.sqrt(2.*math.pi)/par[2]*math.exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])
		return secgaval


	#------------------------------------------#
	def FitFunction4Mass(self, x, par):
		# Total fit function (signal+background+possible second peak)
		#
		bkg = 0
		if self.fTypeOfFit4Bkg == 0:
			bkg = par[0]*par[1]/(np.exp(par[1]*self.fMaxMass)-np.exp(par[1]*self.fMinMass))*np.exp(par[1]*x[0])
		elif self.fTypeOfFit4Bkg == 1:
			bkg = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))
		elif self.fTypeOfFit4Bkg == 2:
			bkg = par[0]/(self.fMaxMass-self.fMinMass)+par[1]*(x[0]-0.5*(self.fMaxMass+self.fMinMass))+par[2]*(x[0]*x[0]-1/3.*(self.fMaxMass*self.fMaxMass*self.fMaxMass-self.fMinMass*self.fMinMass*self.fMinMass)/(self.fMaxMass-self.fMinMass))
		elif self.fTypeOfFit4Bkg == 3:
			bkg = par[0]
		elif self.fTypeOfFit4Bkg == 4:
			bkg = par[0]*(par[1]+1.)/(math.pow(self.fMaxMass-mpi,par[1]+1.)-math.pow(self.fMinMass-mpi,par[1]+1.))*math.pow(x[0]-mpi,par[1])
		elif self.fTypeOfFit4Bkg == 5:
			bkg = par[0]*math.sqrt(x[0] - mpi)*np.exp(-1*par[1]*(x[0]-mpi))
		elif self.fTypeOfFit4Bkg == 6:
			for i in range(1, self.fPolDegreeBkg + 1):
				bkg = bkg + par[i]*math.pow(x[0]-self.fMassParticle,i)/math.factorial(i)
		else:
			print("fTypeOfFit4Bkg out of range")
			return 0
		
		sig = 0
		g1 = 0
		g2 = 0
		if self.fTypeOfFit4Sgn == 0:
			sig = par[self.fNParsBkg]/math.sqrt(2*math.pi)/par[self.fNParsBkg+2]*np.exp(-(x[0]-par[self.fNParsBkg+1])*(x[0]-par[self.fNParsBkg+1])/2/par[self.fNParsBkg+2]/par[self.fNParsBkg+2])
		elif self.fTypeOfFit4Sgn == 1:
			g1=(1.-par[self.fNParsBkg+3])/math.sqrt(2.*math.pi)/par[self.fNParsBkg+2]*np.exp(-(x[0]-par[self.fNParsBkg+1])*(x[0]-par[self.fNParsBkg+1])/2./par[self.fNParsBkg+2]/par[self.fNParsBkg+2])
			g2=par[self.fNParsBkg+3]/math.sqrt(2.*math.pi)/par[self.fNParsBkg+4]*np.exp(-(x[0]-par[self.fNParsBkg+1])*(x[0]-par[self.fNParsBkg+1])/2./par[self.fNParsBkg+4]/par[self.fNParsBkg+4])
			sig = par[self.fNParsBkg+0]*(g1+g2)
		elif self.fTypeOfFit4Sgn == 2:
			g1=(1.-par[self.fNParsBkg+3])/math.sqrt(2.*math.pi)/par[self.fNParsBkg+2]*np.exp(-(x[0]-par[self.fNParsBkg+1])*(x[0]-par[self.fNParsBkg+1])/2./par[self.fNParsBkg+2]/par[self.fNParsBkg+2])
			g2=par[self.fNParsBkg+3]/math.sqrt(2.*math.pi)/(par[self.fNParsBkg+4]*par[self.fNParsBkg+2])*np.exp(-(x[0]-par[self.fNParsBkg+1])*(x[0]-par[self.fNParsBkg+1])/2./(par[self.fNParsBkg+4]*par[self.fNParsBkg+2])/(par[self.fNParsBkg+4]*par[self.fNParsBkg+2]))
			sig = par[self.fNParsBkg+0]*(g1+g2)

		sec = 0
		if(self.fSecondPeak):
			sec = par[self.fNParsBkg+self.fNParsSig+0]/math.sqrt(2.*math.pi)/par[self.fNParsBkg+self.fNParsSig+2]*np.exp(-(x[0]-par[self.fNParsBkg+self.fNParsSig+1])*(x[0]-par[self.fNParsBkg+self.fNParsSig+1])/2./par[self.fNParsBkg+self.fNParsSig+2]/par[self.fNParsBkg+self.fNParsSig+2])

		refl = 0
		if(self.fReflections):
			if not self.fHistoTemplRfl:
				return 0

			binx=self.fHistoTemplRfl.FindBin(x[0])
			value=self.fHistoTemplRfl.GetBinContent(binx)
			binmin=self.fHistoTemplRfl.FindBin(self.fMinMass*1.00001)
			binmax=self.fHistoTemplRfl.FindBin(self.fMaxMass*0.99999)
			norm=self.fHistoTemplRfl.Integral(binmin,binmax)*self.fHistoTemplRfl.GetBinWidth(binx)
			if(abs(value)<1.e-14 and self.fSmoothRfl):
				#very rough, assume a constant trend, much better would be a pol1 or pol2 over a broader range
				value = value + self.fHistoTemplRfl.GetBinContent(binx-1)+self.fHistoTemplRfl.GetBinContent(binx+1)
				value = value/3
  
			refl = par[0]*value/norm*self.fRawYieldHelp*self.fHistoInvMass.GetBinWidth(1)

		return bkg+sig+sec+refl


	#------------------------------------------#
	def Signal(self, nOfSigma, signal, errsignal):
		# Return signal integral in mean +- n sigma
		#

		minMass=self.fMass-nOfSigma*self.fSigmaSgn
		maxMass=self.fMass+nOfSigma*self.fSigmaSgn
		self.SignalinRange(minMass,maxMass,signal,errsignal)
		return

	
	#------------------------------------------#
	def SignalinRange(self, minM, maxM, signal, errsignal):
		# Return signal integral in a range
		#

		signal=self.fSigFunc.Integral(minM, maxM)/self.fHistoInvMass.GetBinWidth(1)
		errsignal=(self.fRawYieldErr/self.fRawYield)*signal #assume relative error is the same as for total integral
		return


	#------------------------------------------#
	def Background(self, nOfSigma, background, errbackground):
		# Return background integral in mean +- n sigma
		#
		
		minMass=self.fMass-nOfSigma*self.fSigmaSgn
		maxMass=self.fMass+nOfSigma*self.fSigmaSgn
		self.BackgroundinRange(minMass,maxMass,background,errbackground)
		return


	#------------------------------------------#
	def BackgroundinRange(self, minM, maxM, background, errbackground):
		# Return background integral in a range
		#
			
		#funcbkg 
		if(self.fBkgFuncRefit): 
			funcbkg=self.fBkgFuncRefit
		elif(self.fBkgFunc):
			funcbkg=self.fBkgFunc
		if not funcbkg:
			print("Bkg function not found!")
			return
			
		intB=funcbkg.GetParameter(0)
		intBerr=funcbkg.GetParError(0)

		#relative error evaluation: from histo

		leftBand=self.fHistoInvMass.FindBin(self.fMass-self.fNSigma4SideBands*self.fSigmaSgn)
		rightBand=self.fHistoInvMass.FindBin(self.fMass+self.fNSigma4SideBands*self.fSigmaSgn)
		intB=self.fHistoInvMass.Integral(1,leftBand)+self.fHistoInvMass.Integral(rightBand,self.fHistoInvMass.GetNbinsX())
		sum2=0
		for i in range(1,leftBand+1):
			sum2 = sum2 + self.fHistoInvMass.GetBinError(i)*self.fHistoInvMass.GetBinError(i)
  
		for i in range (rightBand, (self.fHistoInvMass.GetNbinsX())+1):
			sum2 = sum2 + self.fHistoInvMass.GetBinError(i)*self.fHistoInvMass.GetBinError(i)
		
		intBerr=math.sqrt(sum2)

		background=funcbkg.Integral(minM,maxM)/(self.fHistoInvMass.GetBinWidth(1))
		errbackground=intBerr/intB*background

		return


	#------------------------------------------#
	def Significance(self, nOfSigma, significance, errsignificance):
		# Return significance in mean +- n sigma
		#
		
		minMass=self.fMass-nOfSigma*self.fSigmaSgn
		maxMass=self.fMass+nOfSigma*self.fSigmaSgn
		self.SignificanceinRange(minMass, maxMass, significance, errsignificance)

		return

	#------------------------------------------#
	def SignificanceinRange(self, minM, maxM, significance, errsignificance):
		# Return significance integral in a range
		#

		background = ROOT.Double()
		errbackground = ROOT.Double()
		self.BackgroundinRange(minM,maxM,background,errbackground)
		
		if (self.fRawYield+background <= 0.):
				significance=-1
				errsignificance=0
				return

		self.ComputeSignificance(self.fRawYield,self.fRawYieldErr,background,errbackground,significance,errsignificance)

		return


	#------------------------------------------#
	def ComputeSignificance(self, signal, errsignal, background ,errbackground ,significance, errsignificance):
		errSigSq=errsignal*errsignal
		errBkgSq=errbackground*errbackground
		sigPlusBkg=signal+background
		if(sigPlusBkg>0. and signal>0.):
			significance =  signal/math.sqrt(signal+background)
			errsignificance = significance*math.sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal)
		else:
			significance=0.
			errsignificance=0.
		return


	#------------------------------------------#
	def PrepareHighPolFit(self, fback):
		#/ Perform intermediate fit steps up to fPolDegreeBkg-1
		#/ in case of fit with a polynomial with degree > 2 (fTypeOfFit4Bkg=6)

		estimatecent=0.5*(self.fHistoInvMass.GetBinContent(self.fHistoInvMass.FindBin(self.fMass-3.5*self.fSigmaSgn))+self.fHistoInvMass.GetBinContent(self.fHistoInvMass.FindBin(self.fMass+3.5*self.fSigmaSgn)))# just a first rough estimate
		estimateslope=(self.fHistoInvMass.GetBinContent(self.fHistoInvMass.FindBin(self.fMass+3.5*self.fSigmaSgn))-self.fHistoInvMass.GetBinContent(self.fHistoInvMass.FindBin(self.fMass-3.5*self.fSigmaSgn)))/(7*self.fSigmaSgn)# first rough estimate

		self.fCurPolDegreeBkg=2
		funcbkg = 0x0
		funcPrev=0x0
		while(self.fCurPolDegreeBkg <= self.fPolDegreeBkg):
			fname = "temp{}".format(self.fCurPolDegreeBkg)
			funcbkg = ROOT.TF1(fname,self.BackFitFuncPolHelper,self.fMinMass,self.fMaxMass,self.fCurPolDegreeBkg+1)
			if(funcPrev):
				for j in range(0,self.fCurPolDegreeBkg): # now is +1 degree w.r.t. previous fit funct
					funcbkg.SetParameter(j,funcPrev.GetParameter(j))
				funcPrev.delete
			else:
				funcbkg.SetParameter(0,estimatecent)
				funcbkg.SetParameter(1,estimateslope)
			print("   --=> Pre-fit of background with pol degree {} ---".format(self.fCurPolDegreeBkg))
			self.fHistoInvMass.Fit(funcbkg,"REMN","")
			funcPrev=funcbkg.Clone("ftemp")
			funcbkg.delete
			self.fCurPolDegreeBkg = self.fCurPolDegreeBkg + 1
  
		for j in range(0,self.fPolDegreeBkg+1):
			fback.SetParameter(j,funcPrev.GetParameter(j))
			fback.SetParError(j,funcPrev.GetParError(j))
		print("   --=> Final background fit with pol degree {} ---",self.fPolDegreeBkg)
		fitf = "R,{},+,0".format(self.fFitOption)
		self.fHistoInvMass.Fit(fback,fitf)# THIS IS JUST TO SET NOT ONLY THE PARAMETERS BUT ALSO chi2, etc...

		# The following lines might be useful for debugging
		#   TCanvas *cDebug=new TCanvas()
		#   cDebug.cd()
		#   fHistoInvMass.Draw()
		#   TString strout=Form("Test%d.root",(Int_t)fhistoInvMass.GetBinContent(fhistoInvMass.FindBin(fMass)))
		#   cDebug.Print(strout.Data())
		# delete cDebug

		funcPrev.delete
		return True


	#------------------------------------------#
	def BackFitFuncPolHelper(self, x, par):
		# Helper function for polynomials with degree>2
		#

		maxDeltaM = self.fNSigma4SideBands*self.fSigmaSgn
		if(self.fOnlySideBands and abs(x[0]-self.fMass) < self.maxDeltaM):
			ROOT.TF1.RejectPoint()
			return 0
  
		back=par[0]
		for it in range(1,self.fCurPolDegreeBkgit+1):
			back = back + par[it]*math.power(x[0]-self.fMassParticle,it)/math.factorial(it)
  
		return back


	#------------------------------------------#
	def SetTemplateReflections(self, h, fopt, minRange, maxRange):
		#/ Method to create the reflection invariant mass distributions from MC templates
		#/ option could be:
		#/    "template"                use MC histograms
		#/    "1gaus" ot "singlegaus"   single gaussian function fit to MC templates
		#/    "2gaus" ot "doublegaus"   double gaussian function fit to MC templates
		#/    "pol3"                    3rd order polynomial fit to MC templates
		#/    "pol6"                    6th order polynomial fit to MC templates
		
		if not h:
			self.fReflections=False
			return 0x0
  
		self.fHistoTemplRfl=h.Clone("hTemplRfl")
		fopt = fopt.lower()
		self.fReflections=True
		print("--- Reflection templates from simulation ---")
		if "templ" in fopt:
			print("   --=> Reflection contribution using directly the histogram from simulation")
			return self.fHistoTemplRfl

		f=0x0
		isPoissErr=True
		xMinForFit=h.GetBinLowEdge(1)
		xMaxForFit=h.GetXaxis().GetBinUpEdge(h.GetNbinsX())
		if(minRange>=0 and maxRange>=0):
			xMinForFit=max(minRange,h.GetBinLowEdge(1))
			xMaxForFit=min(maxRange,h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))

		if(fopt == "1gaus" or fopt == "singlegaus"):
			print("   --=> Reflection contribution from single-Gaussian fit to histogram from simulation")
			f=ROOT.TF1("mygaus","gaus",xMinForFit,xMaxForFit)
			f.SetParameter(0,h.GetMaximum())
			#    f.SetParLimits(0,0,100.*h.Integral())
			f.SetParameter(1,1.865)
			f.SetParameter(2,0.050)
			self.fHistoTemplRfl.Fit(f,"REM0","")#,h.GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
		elif(fopt == "2gaus" or fopt == "doublegaus"):
			print("   --=> Reflection contribution from double-Gaussian fit to histogram from simulation")
			f=ROOT.TF1("my2gaus","[0]*([3]/( sqrt(2.*pi)*[2])*exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+(1.-[3])/( sqrt(2.*pi)*[5])*exp(-(x-[4])*(x-[4])/(2.*[5]*[5])))",xMinForFit,xMaxForFit)
			f.SetParameter(0,h.GetMaximum())
			#    f.SetParLimits(0,0,100.*h.Integral())
			f.SetParLimits(3,0.,1.)
			f.SetParameter(3,0.5)
			f.SetParameter(1,1.84)
			f.SetParameter(2,0.050)
			f.SetParameter(4,1.88)
			f.SetParameter(5,0.050)
			self.fHistoTemplRfl.Fit(f,"REM0","")#,h.GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))		
		elif(fopt == "pol3"):
			print("   --=> Reflection contribution from pol3 fit to histogram from simulation")
			f=ROOT.TF1("mypol3","pol3",xMinForFit,xMaxForFit)
			f.SetParameter(0,h.GetMaximum())
			#    f.SetParLimits(0,0,100.*h.Integral())
			# Hard to initialize the other parameters...
			self.fHistoTemplRfl.Fit(f,"REM0","")
			#    Printf("We USED %d POINTS in the Fit",f.GetNumberFitPoints())
		elif(fopt == "pol6"):
			print("   --=> Reflection contribution from pol6 fit to histogram from simulation")
			f=ROOT.TF1("mypol6","pol6",xMinForFit,xMaxForFit)
			f.SetParameter(0,h.GetMaximum())
			#    f.SetParLimits(0,0,100.*h.Integral())
			# Hard to initialize the other parameters...
			self.fHistoTemplRfl.Fit(f,"RLEMI0","")#,h.GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
		else:
			# no good option passed
			print("   --=> Bad option for reflection configuration . reflections will not be included in the fit")
			self.fReflections=False
			self.fHistoTemplRfl.delete
			self.fHistoTemplRfl=0x0
			return 0x0
  
		# Fill fHistoTemplRfl with values of fit function
		if f:
			for j in range(1,self.fHistoTemplRfl.GetNbinsX()+1):
				self.fHistoTemplRfl.SetBinContent(j,f.Integral(self.fHistoTemplRfl.GetBinLowEdge(j),self.fHistoTemplRfl.GetXaxis().GetBinUpEdge(j))/self.fHistoTemplRfl.GetBinWidth(j))
				if(self.fHistoTemplRfl.GetBinContent(j)>=0 and abs(h.GetBinError(j)*h.GetBinError(j)-h.GetBinContent(j))>0.1*h.GetBinContent(j)):
					isPoissErr=False
			for j in range(1,self.fHistoTemplRfl.GetNbinsX()+1):
				if(isPoissErr):
					if(self.fHistoTemplRfl.GetBinContent(j)>0):
						self.fHistoTemplRfl.SetBinError(j,math.sqrt(self.fHistoTemplRfl.GetBinContent(j)))
					else:
						self.fHistoTemplRfl.SetBinError(j,0)
				else:
					self.fHistoTemplRfl.SetBinError(j,0.001*self.fHistoTemplRfl.GetBinContent(j))
			self.fReflections=True
			return self.fHistoTemplRfl
		else:
			print("   --=> Fit to MC template for reflection failed => reflections will not be included in the fit")
			self.fReflections=False
			self.fHistoTemplRfl.delete
			self.fHistoTemplRfl=0x0
			return 0x0
		
		return 0x0


	#-----------------------------------------#
	def SetInitialReflOverS(self, rovers):
		self.fRflOverSig = rovers

	#-----------------------------------------#
	def SetFixReflOverS(self, ros):
		self.SetInitialReflOverS(ros)
		self.fFixRflOverSig = True

	#------------------------------------------#
	def GetRawYieldBinCountingByPDG(self, errRyBC, nOfSigma, option, pdgCode):
		#/ Method to compute the signal using inv. mass histo bin counting
		#/ => interface method to compute yield in nsigma range around peak
		#/ pdgCode: if==411,421,413,413 or 4122: range defined based on PDG mass
		#           else (default) mean of gaussian fit

		massVal=self.fMass
		switchpdgCode = {411: 1.86965,
						 421: 1.86483,
						 431: 1.96834,
						 4122: 2.28646,
						 413: 0.14543
					 }
		massVal = switchpdgCode.get(pdgCode, self.fMass)
  
		minMass=massVal-nOfSigma*self.fSigmaSgn
		maxMass=massVal+nOfSigma*self.fSigmaSgn
		return self.GetRawYieldBinCounting(errRyBC,minMass,maxMass,option)


	#------------------------------------------#
	def GetRawYieldBinCounting(self, errRyBC, minMass, maxMass, option):
		#/ Method to compute the signal using inv. mass histo bin counting
		#/ after background subtraction from background fit function
		#/   option=0: background fit function from 1st fit step (only side bands)
		#/   option=1: background fit function from 2nd fit step (S+B)
		
		minBinSum=self.fHistoInvMass.FindBin(minMass)
		maxBinSum=self.fHistoInvMass.FindBin(maxMass)
		if(minBinSum<1):
			print("Left range for bin counting smaller than allowed by histogram axis, setting it to the lower edge of the first histo bin")
			minBinSum=1
		
		if(maxBinSum>self.fHistoInvMass.GetNbinsX()):
			print("Right range for bin counting larger than allowed by histogram axis, setting it to the upper edge of the last histo bin")
			maxBinSum=self.fHistoInvMass.GetNbinsX()
		
		cntSig=0.
		cntErr=0.
		errRyBC=0
		fbackground=self.fBkgFunc
		if(option==1):
			fbackground=self.fBkgFuncRefit
		if not fbackground:
			return 0.

		for jb in range(minBinSum,maxBinSum+1):
			cntTot=self.fHistoInvMass.GetBinContent(jb)
			cntBkg=fbackground.Integral(self.fHistoInvMass.GetBinLowEdge(jb),self.fHistoInvMass.GetBinLowEdge(jb)+self.fHistoInvMass.GetBinWidth(jb))/self.fHistoInvMass.GetBinWidth(jb)
			cntRefl=0
			if(option==1 and self.fRflFunc):
				cntRefl=self.fRflFunc.Integral(self.fHistoInvMass.GetBinLowEdge(jb),self.fHistoInvMass.GetBinLowEdge(jb)+self.fHistoInvMass.GetBinWidth(jb))/self.fHistoInvMass.GetBinWidth(jb)
			#Double_t cntBkg=fbackground.Eval(fHistoInvMass.GetBinCenter(jb))
			cntSecPeak=0
			if(option==1 and self.fSecondPeak and self.fSecFunc):
				cntSecPeak=self.fSecFunc.Integral(self.fHistoInvMass.GetBinLowEdge(jb),self.fHistoInvMass.GetBinLowEdge(jb)+self.fHistoInvMass.GetBinWidth(jb))/self.fHistoInvMass.GetBinWidth(jb)
			cntSig = cntSig + (cntTot-cntBkg-cntRefl-cntSecPeak)
			cntErr = cntErr + (self.fHistoInvMass.GetBinError(jb)*self.fHistoInvMass.GetBinError(jb))
  
		errRyBC=math.sqrt(cntErr)
		return cntSig


	#------------------------------------------#
	def GetResidualsAndPulls(self, hPulls, hResidualTrend, hPullsTrend, minrange, maxrange, option):
		#/ fill and return the residual and pull histos

		binmi=self.fHistoInvMass.FindBin(self.fMinMass*1.001)
		binma=self.fHistoInvMass.FindBin(self.fMaxMass*0.9999)
		if(maxrange>minrange):
			binmi=self.fHistoInvMass.FindBin(minrange*1.001)
			binma=self.fHistoInvMass.FindBin(maxrange*0.9999)
		
		if(hResidualTrend):
			#fHistoInvMass.Copy(hResidualTrend)
			hResidualTrend.SetBins(self.fHistoInvMass.GetNbinsX(),self.fHistoInvMass.GetXaxis().GetXmin(),self.fHistoInvMass.GetXaxis().GetXmax())
			hResidualTrend.SetName("{}_residualTrend".format(self.fHistoInvMass.GetName()))
			hResidualTrend.SetTitle("{}  (Residuals)".format(self.fHistoInvMass.GetTitle()))
			hResidualTrend.SetMarkerStyle(20)
			hResidualTrend.SetMarkerSize(1.0)
			hResidualTrend.Reset()
  
		if(hPullsTrend):
			hPullsTrend.SetBins(self.fHistoInvMass.GetNbinsX(),self.fHistoInvMass.GetXaxis().GetXmin(),self.fHistoInvMass.GetXaxis().GetXmax())
			hPullsTrend.Reset()
			hPullsTrend.SetName("{}_pullTrend".format(self.fHistoInvMass.GetName()))
			hPullsTrend.SetTitle("{} (Pulls)".format(self.fHistoInvMass.GetTitle()))
			hPullsTrend.SetMarkerStyle(20)
			hPullsTrend.SetMarkerSize(1.0)
			
		if(hPulls):
			hPulls.SetName("{}_pulls".format(self.fHistoInvMass.GetName()))
			hPulls.SetTitle("{} ; Pulls".format(self.fHistoInvMass.GetTitle()))
			hPulls.SetBins(40,-10,10)
			hPulls.Reset()
  
		res=-1.e-6
		minr=1.e+12
		maxr=-1.e+12
		arval=ROOT.TArrayD(binma-binmi+1)
		for jst in range(1,self.fHistoInvMass.GetNbinsX()+1):
			integFit=0
			if(option==0):
				integFit=self.fTotFunc.Integral(self.fHistoInvMass.GetBinLowEdge(jst),self.fHistoInvMass.GetBinLowEdge(jst)+self.fHistoInvMass.GetBinWidth(jst))
			else:
				integFit=self.fBkgFuncRefit.Integral(self.fHistoInvMass.GetBinLowEdge(jst),self.fHistoInvMass.GetBinLowEdge(jst)+self.fHistoInvMass.GetBinWidth(jst))
				if(option==2):
					integFit = integFit + self.fRflFunc.Integral(self.fHistoInvMass.GetBinLowEdge(jst),self.fHistoInvMass.GetBinLowEdge(jst)+self.fHistoInvMass.GetBinWidth(jst))
			res=self.fHistoInvMass.GetBinContent(jst)-integFit/self.fHistoInvMass.GetBinWidth(jst)
			if(jst>=binmi and jst<=binma):
				arval.AddAt(res,jst-binmi)
				if(res<minr):
					minr=res
				if(res>maxr):
					maxr=res
			#Printf("Res = %f from %f - %f",res,fHistoInvMass.GetBinContent(jst),fTotFunc.Integral(fHistoInvMass.GetBinLowEdge(jst),fHistoInvMass.GetBinLowEdge(jst)+fHistoInvMass.GetBinWidth(jst))/fHistoInvMass.GetBinWidth(jst))
			if(hResidualTrend):
				hResidualTrend.SetBinContent(jst,res)
				hResidualTrend.SetBinError(jst,self.fHistoInvMass.GetBinError(jst))    
			if(hPulls):
				if(jst>=binmi and jst<=binma):
					hPulls.Fill(res/self.fHistoInvMass.GetBinError(jst))
			if(hPullsTrend):
				hPullsTrend.SetBinContent(jst,res/self.fHistoInvMass.GetBinError(jst))
				hPullsTrend.SetBinError(jst,0.0001)
    
		if(hResidualTrend):
			hResidualTrend.GetXaxis().SetRange(binmi,binma)
			if(option!=0):
				fgauss=ROOT.TF1("signalTermForRes","[0]/math.sqrt(2.*math.pi)/[2]*math.exp(-(x-[1])*(x-[1])/2./[2]/[2])",self.fHistoInvMass.GetBinLowEdge(1),self.fHistoInvMass.GetBinLowEdge(self.fHistoInvMass.GetNbinsX()+1))
				fgauss.SetParameter(0,self.fRawYield*self.fHistoInvMass.GetBinWidth(1))
				fgauss.SetParameter(1,self.fMass)
				fgauss.SetParameter(2,self.fSigmaSgn)
				fgauss.SetLineColor(ROOT.kBlue)
				hResidualTrend.GetListOfFunctions().Add(fgauss)
	
		if(hPullsTrend):
			hPullsTrend.GetXaxis().SetRange(binmi,binma)
			hPullsTrend.SetMinimum(-7)
			hPullsTrend.SetMaximum(+7)
  
		if(abs(minr)>abs(maxr)):
			maxr=minr

		hname = "{}_residuals".format(self.fHistoInvMass.GetName())
		htitle = "{} residuals".format(self.fHistoInvMass.GetTitle())
		hout=ROOT.TH1F(hname, htitle, 25,-abs(maxr)*1.5,abs(maxr)*1.5)
		for j in range(0,binma-binmi+1):
			hout.Fill(arval.At(j))
		hout.Sumw2()
		hout.Fit("gaus","LEM","",-abs(maxr)*1.2,abs(maxr)*1.2)

		if(hPulls):
			hPulls.Sumw2()
			hPulls.Fit("gaus","LEM","",-3,3)

		arval.delete
		return hout


	#------------------------------------------#
	def GetOverBackgroundResidualsAndPulls(self, hPulls, hResidualTrend, hPullsTrend, minrange, maxrange):
		return self.GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend,minrange,maxrange,1)


	#------------------------------------------#
	def GetOverBackgroundPlusReflResidualsAndPulls(self, hPulls, hResidualTrend, hPullsTrend, minrange, maxrange):
		return self.GetResidualsAndPulls(hPulls,hResidualTrend,hPullsTrend,minrange,maxrange,2)


	#------------------------------------------#
	def PrintFunctions(self):
		#/ dump the function parameters
		#/
		if(self.fBkgFunc):
			printf("--- Background function in 1st step fit to side bands ---\n")
			self.fBkgFunc.Print()
		if(self.fTotFunc):
			printf("--- Total fit function from 2nd step fit ---\n")
			self.fTotFunc.Print()
		if(self.fBkgFuncRefit):
			printf("--- Background function in 2nd step fit ---\n")
			self.fBkgFuncRefit.Print()
		if(self.fRflFunc):
			printf("--- Reflections ---\n")
			self.fRflFunc.Print()
		if(self.fBkRFunc):
			printf("--- Background + reflections ---\n")
			self.fBkRFunc.Print()
		if(self.fSecFunc):
			printf("--- Additional Gaussian peak ---\n")
			self.fSecFunc.Print()
			


#h = ROOT.TH1F("h", "h", 20, 1.7, 2)
#for i in range(0,1000):
#	h.Fill(ROOT.gRandom.Gaus(1.86,0.02))
#h.FillRandom("gaus",1000)
#a1 = AliHFInvMassFitter(h, 1.7, 2)
#a1.fReflections = True
#a1.MassFitter(True)
#sig = h.Integral()
#a1.fBkgFunc = a1.CreateBackgroundFitFunction("bkgfit",sig)
#a1.fSigFunc = a1.CreateSignalFitFunction("sigfit",sig)
#totfit = a1.CreateTotalFitFunction("totfit")
#a1.SetNumberOfParams()
#totparams = a1.fNParsBkg + a1.fNParsSig
#totfit = ROOT.TF1("totfit",a1.FitFunction4Bkg,1.7,2,a1.fNParsBkg) 
#h.Fit("totfit")
