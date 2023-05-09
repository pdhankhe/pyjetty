#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as mp
import ROOT
import yaml
import math
ROOT.gROOT.SetBatch(True)
from pyjetty.mputils import treereader
from bisect import bisect_right

global pt_low, pt_high, pt_mid, pt_edge, nbins
pt_low = [3, 4, 5, 6, 7, 8, 10, 12, 16, 24]
pt_high = [4, 5, 6, 7, 8, 10, 12, 16, 24, 36]
pt_mid = [3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 14, 20, 30]
pt_edge = array('d',[3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 36])
nbins = len(pt_low)

global cut_names, cut_type, cut_value, branch_names, gen_branch_names

branch_names = ['jet_reco', 'Dpt', 'Dm','D_cos_t_star', 'cos_p', 'D_norm_dlxy','D_imp_par_prod','Dmeson_cand_type']
gen_branch_names= ['jet_pt', 'jet_eta', 'D_pt', 'D_eta', 'Dmeson_cand_type']
cut_names = ['D_cos_t_star', 'cos_p', 'D_norm_dlxy','D_imp_par_prod']

cut_type = ['abs', 'greater', 'greater', 'less']

cut_value = [
[0.8, 0.95, 5, -0.0003],
[0.8, 0.95, 5, -0.00015],
[0.8, 0.95, 4, -0.0001],
[0.8, 0.95, 4, -0.00008],
[0.8, 0.95, 4, -0.00008],
[0.9, 0.95, 3, -0.00005],
[0.9, 0.95, 3, -0.00005],
[1, 0.95, 3, 0.001],
[1, 0.90, 3, 0.001],
[1, 0.90, 3, 0.001]]

    

global num_of_sig

num_of_sig = {"peaklow": -2,
					 "leftsidelow": -8,
					 "leftsidehigh": -4,
					 "peakhigh": 2,
					 "rightsidelow": 4,
					 "rightsidehigh": 8}

global prompt, fd, ref, ctype
prompt = 0b1000
fd = 0b10000
ref = 0b100000
ctype = {
        "selected":   0b000000001,
        "signal":     0b000000010,
        "background": 0b000000100,
        "prompt":     0b000001000,
        "fd":         0b000010000,
        "reflection": 0b000100000,
        "selTopo":    0b001000000,
        "selPID":     0b010000000,
        "selTrack":   0b100000000,
        "all":        0b111111111
}

def make_var_hist(name, bins, xarray):
        h = ROOT.TH1F(name, name, bins, xarray)
        return h


def make_cut(tree, cut, type, value):
	if type == 'greater':
		if getattr(tree, cut)[0] > value:
			return True
	if type == 'less':
		if getattr(tree, cut)[0] < value:
			return True
	if type == 'abs':
		if abs(getattr(tree, cut)[0]) < value:
			return True
	return False


def make_hist(name, nxbins, xlow, xhigh, nybins=-1, ylow=-1, yhigh=-1):
	if nybins < 0:
		hist = ROOT.TH1F(name, name, nxbins, xlow, xhigh)
	else:
		hist = ROOT.TH2F(name, name, nxbins, xlow, xhigh, nybins, ylow, yhigh)
	hist.Sumw2()
	return hist


def get_hist_from_file(rootfile, histname):
	hist = rootfile.Get(histname)
	hist.Sumw2()
	hist.SetDirectory(0)
	return hist

def get_bin(hist, val, x=True):
	if x:
		bin = hist.GetXaxis().FindBin(val)
	else:
		bin = hist.GetYaxis().FindBin(val)
	return bin


def get_bin_edge(hist, val, low=True, x=True):
	if x:
		if low:
			edge = hist.GetXaxis().GetBinLowEdge(get_bin(hist, val, x))
		else:
			edge = hist.GetXaxis().GetBinUpEdge(get_bin(hist, val, x))
	else:
		if low:
			edge = hist.GetYaxis().GetBinLowEdge(get_bin(hist, val, x))
		else:
			edge = hist.GetYaxis().GetBinUpEdge(get_bin(hist, val, x))
	return edge


def calc_alpha(hist, bg_fitfunc, mean, sigma):
	bg = [0, 0, 0]
	section = ["peak", "leftside", "rightside"]
	for i in range(3):
		bg[i] = bg_fitfunc.Integral(get_bin_edge(hist, mean+num_of_sig[section[i]+"low"], True), get_bin_edge(hist, mean+num_of_sig[section[i]+"high"], False))
	try:
		alpha = bg[0]/(bg[1]+bg[2])
	except ZeroDivisionError:
		alpha = 0
	return alpha

def get_raw_counts(hist, mean, sigma, ptlow, pthigh, ct, cterr):
	ct = hist.Integral(get_bin(hist, mean-2*sigma), get_bin(hist, mean+2*sigma), ptlow, pthigh)
	cterr = math.sqrt(sig)


#def calc_sig_and_bg(hist, bgfunc, mean, sigma, ptlow, pthigh, bgcount, bgcounterr, sig, sigerr):
def calc_sig_and_bg(hist, bgfunc, mean, sigma, ptlow, pthigh):
	bgcount = hist.Integral(get_bin(hist, mean-8*sigma), get_bin(hist, mean-4*sigma), ptlow, pthigh) + hist.Integral(get_bin(hist, mean+4*sigma), get_bin(hist, mean+8*sigma), ptlow, pthigh)
	bgcounterr = math.sqrt(bgcount)

	bglow = bgfunc.Integral(get_bin_edge(hist, mean-8*sigma), get_bin_edge(hist, mean-4*sigma, False))
	bghigh = bgfunc.Integral(get_bin_edge(hist, mean+4*sigma), get_bin_edge(hist, mean+8*sigma, False))
	bgpeak = bgfunc.Integral(get_bin_edge(hist, mean-2*sigma), get_bin_edge(hist, mean+2*sigma, False))
	if((bglow + bghigh)==0):
		alpha = 0
	else:
		alpha = bgpeak/(bglow + bghigh)

	bgcount = bgcount*alpha
	bgcounterr = bgcounterr*alpha

	get_raw_counts(hist, mean, sigma, ptlow, pthigh, sig, sigerr)	

	yields = {
		"rawcounts": sig,
		"yielderr": sigerr,
		"bg": bgcount,
		"bgerr": bgcounterr,
		"subyield": ((sig-bgcount)/0.9545)
	}

	return yields
        

def calc_sob_and_sig(bg, signal, sob, sig):
	if(bgcount==0):
		sob = 0
	else:
		sob = signal/bg

	sig = signal/(math.sqrt(signal+bg))


def setup_hist_to_draw(hist, title, xtitle, ytitle, xrange=[], yrange=[]):
	hist.SetTitle(title)
	hist.SetXTitle(xtitle)
	hist.SetYTitle(ytitle)
	if xrange:
		hist.GetXaxis().SetRangeUser(xrange[0], xrange[1])
	if yrange:
		hist.GetYaxis().SetRangeUser(yrange[0], yrange[1])


def fill_hist_from_tree(treename, branchnames, filename, hist, cutnames, cuttype, cutvalue, candtype=["all"]):
	tr = treereader.RTreeReader(tree_name=treename, branches=branchnames, file_name=filename)
	for i in range(tr.tree.GetEntries()):
		tr.tree.GetEntry(i)

		if(i%10000==0):
			print(i)

		p = bisect_right(pt_low, tr.pt_cand[0]) - 1
		
		if((p < 0) | (tr.pt_cand[0] > 36)):
			continue

		for c in range(len(cutvalue[0])):
			flag = make_cut(tr, cutnames[c], cuttype[c], cutvalue[p][c])
			if not flag:
				break

		if flag:
			for ct in range(len(candtype)):
				type = int(tr.cand_type[0]) & ctype.get(candtype[ct], 0b111111111)
				if(type):
					hist[ct].Fill(tr.inv_mass[0], tr.pt_cand[0])
 

def set_fit_function(func, func_name, xrange=[], param=[], parrange=[]):
	ff = {
		"gaus": "gaus",
		"pol0": "pol0",
		"pol1": "pol1",
		"pol2": "pol2",
		"expo": "expo",
		"dgaus": "gaus(0) + gaus(3)",
		"totexp": "gaus(0) + expo(3)",
		"totpol": "gaus(0) + pol2(3)"
	}

	f1 = ROOT.TF1(func_name, ff[func], xrange[0], xrange[1])
	if param:
		par_arr = array('d',param)
		f1.SetParameters(par_arr)
	if parrange:
		for parval in parrange:
			f1.SetParLimits(parval[0], parval[1], parval[2])
	return f1


def fit_and_get_parameters(hist, fitfunc, numpar):
	hist.Fit(fitfunc, "R")
	fit_par_and_err = [[],[]]
	for i in range(numpar):
		fit_par_and_err[0].insert(i, fitfunc.GetParameter(i))
		fit_par_and_err[1].insert(i, fitfunc.GetParError(i))
	return fit_par_and_err


def run_over_tree(hnamelist, root_filename, input_file, cutvalues):
	#root_filename = root_dir + ofile
	rootfile = ROOT.TFile(root_filename, 'RECREATE')
	rootfile.Close()

	h_invmass_2D = []
	cand_list = []	

	for h in range(len(hnamelist)):
		name = "InvMass_vs_pt_" + hnamelist[h]
		h_invmass_2D.insert(h, make_hist(name, 55, 1.6, 2.15, 40, 0, 40))
	
	fill_hist_from_tree('d0', branch_names, input_file, h_invmass_2D, cut_names, cut_type, cutvalues, hnamelist)

	rootfile = ROOT.TFile(root_filename, 'UPDATE')
	for h in range(len(hnamelist)):
		h_invmass_2D[h].Write()
	rootfile.Close()


def open_hist(histname, root_filename):
	#root_filename = root_dir + ofile
	rootfile = ROOT.TFile.Open(root_filename, "READ")

	name = "InvMass_vs_pt_" + histname
	hist = get_hist_from_file(rootfile, name)

	rootfile.Close()
	return hist


def make_1D_hist(hist2D, histname, axisrange, x=True):
	if x:
		hist = hist2D.ProjectionY(histname, get_bin(hist2D, axisrange[0]), get_bin(hist2D, axisrange[1]) - 1)
	else:
		hist = hist2D.ProjectionX(histname, get_bin(hist2D, axisrange[0], False), get_bin(hist2D, axisrange[1], False) - 1)
	return hist


def integrate_hist(hist, xvalues=[], yvalues=[]):
	if xvalues:
		if yvalues:
			integral_val = hist.Integral(get_bin(hist, xvalues[0]), get_bin(hist, xvalues[1]), get_bin(hist, yvalues[0], False), get_bin(hist, yvalues[1], False))
		else:
			integral_val = hist.Integral(get_bin(hist, xvalues[0]), get_bin(hist, xvalues[1]))
	else:
		integral_val = hist.Integral()

	return integral_val


def sethistcontent(hist, xpoint, val, valerr):
	bin = get_bin(hist, xpoint)
	hist.SetBinContent(bin, val)
	hist.SetBinError(bin, valerr)


def make_hists(hnames=[]):
	hists =	[]
	for i, name in enumerate(hnames):
		hists.insert(i,	make_var_hist(name, nbins, pt_edge))
		hists[i].Sumw2()
	return hists
