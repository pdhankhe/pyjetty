#!/usr/bin/env python3

import numpy as np
import argparse
import os
from array import array
import pyjetty.alihfjets.dev.hfjet.process_io_event_hf as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter, jet_analysis

import yaml
import fastjet as fj
import fjext
import fjtools
import fjcontrib

import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)

		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()

		self.hNevents = ROOT.TH1F("hNevents", "hNevents", 2, array('f', [-0.5, 0.5, 1.5]) )
		self.hNevents.Sumw2()


	
		print(self.config_file)
	
		with open(self.config_file, 'r') as stream:
			self.config = yaml.load(stream,Loader=yaml.FullLoader)


	def analysis(self, df):
		#print(df)
		self.hNevents.Fill(1,df["ev_id"].nunique())


	def finalize(self):

		self.hNevents.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-c', '--configFile', help='Path of config file for analysis', type=str,metavar='configFile',default='config/configcuts.yaml', required=True)
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()
	
	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)


	hfaio = hfdio.HFAnalysisIO()
	
	hfa = HFAnalysisInvMass(config_file=args.configFile, name = args.output)

	hfa.event_selection.add_selection_range_abs('z_vtx_reco', 10)
	hfa.event_selection.add_selection_equal('is_ev_rej', 0)	
	
	#topomatic cut suggested by D2H

	hfaio.add_analysis(hfa)	
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)
	
	hfa.finalize()
