#!/usr/bin/env python3

from pyjetty.mputils import MPBase, pwarning, pinfo, perror, treewriter
import random
import uproot
import pandas as pd
import numpy as np
import fastjet as fj
import fjext
import os
import tqdm
import argparse


class HFAIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0', 
								 track_tree_name='PWGHF_TreeCreator/tree_Particle',
								 event_tree_name='PWGHF_TreeCreator/tree_event_char',
								 enable_jet=True,
								 enable_d0=True,
								 offset_parts=0)
		super(HFAIO, self).__init__(**kwargs)
		self.analyses = []
		self.df_grouped = None
		self.df_events = None

		# temp output
		self.tw = treewriter.RTreeWriter(name = 'd0j', file_name = 'djet_tout.root')

	def __del__(self):
		self.tw.write_and_close()

	def pd_tree(self, path, tname, squery=None):
		try:
			tree = uproot.open(path)[tname]
		except:
			pwarning('error getting', tname, 'from file:', path)
			return None
		if not tree:
			perror('Tree {} not found in file {}'.format(tname, path))
			return None
		df = tree.pandas.df()
		if squery:
			df = df.query(squery)
			df.reset_index(drop=True)
		return df


	def process_file(self, fname):
		_ev_cuts = "is_ev_rej == 0 & abs(z_vtx_reco) < 10"
		self.event_df = self.pd_tree(path=fname, tname=self.event_tree_name, squery=_ev_cuts)
		if self.event_df is None:
			return False
		pinfo('events from', fname, len(self.event_df))

		_d0cuts_base = "(pt_cand > 2.0 & pt_prong0 > 0.5 & pt_prong1 > 0.5 & abs(eta_cand) < 0.8) & "
		_d0cuts_kpi = _d0cuts_base
		_d0cuts_kpi += "(abs(nsigTPC_Pi_0) < 3. & (abs(nsigTOF_Pi_0) < 3. | nsigTOF_Pi_0 < -900) & abs(nsigTPC_K_1) < 3. & (abs(nsigTOF_K_1) < 3. | nsigTOF_K_0 < -900)) | "
		_d0cuts_kpi += "(abs(nsigTPC_Pi_1) < 3. & (abs(nsigTOF_Pi_1) < 3. | nsigTOF_Pi_1 < -900) & abs(nsigTPC_K_0) < 3. & (abs(nsigTOF_K_0) < 3. | nsigTOF_K_0 < -900))"
		self.d0_df = self.pd_tree(path=fname, tname=self.d0_tree_name, squery=_d0cuts_kpi)
		if self.d0_df is None:
			return False
		pinfo('d0s from', fname, len(self.d0_df))

		self.track_df = self.pd_tree(path=fname, tname=self.track_tree_name)
		# self.track_df = _track_df.groupby(['run_number','ev_id'])
		if self.track_df is None:
			return False
		pinfo('tracks from', fname, len(self.track_df))

		self.event_df.apply(self.process_event, axis=1)

	def process_event(self, df):
		if 'ev_id_ext' in df.keys():
			_ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'], df['ev_id'], df['ev_id_ext'])
		else:
			_ev_query = "run_number == {} & ev_id == {}".format(df['run_number'], df['ev_id'])
		_df_d0 = self.d0_df.query(_ev_query)
		_df_d0.reset_index(drop=True)
		_n_d0s = len(_df_d0.index)
		if _n_d0s < 1:
			return
		_df_tracks = self.track_df.query(_ev_query)
		_df_tracks.reset_index(drop=True)
		_parts = fjext.vectorize_pt_eta_phi(_df_tracks['ParticlePt'].values, _df_tracks['ParticleEta'].values, _df_tracks['ParticlePhi'].values)
		_d0s = fjext.vectorize_pt_eta_phi(_df_d0['pt_cand'].values, _df_d0['eta_cand'].values, _df_d0['phi_cand'].values)
		_d0s_gh = [p * 1.e-6 for p in _d0s]
		# now the jet finding with ghost D0s
		# write jets with constituents ... and/or fill histograms

		# some debug info
		# pinfo(_df_d0.keys())
		# if 'ev_id_ext' in df.keys():
		# 	pinfo('run#:', df['run_number'], 'ev_id:', df['ev_id'], 'ev_id_ext:', df['ev_id_ext'], 'ntrk:', len(_parts), 'nd0s:', _n_d0s)
		# else:
		# 	pinfo('run#:', df['run_number'], 'ev_id:', df['ev_id'], 'ntrk:', len(_parts), 'nd0s:', _n_d0s)

		self.tw.fill_branches(dpsj = _d0s)
		self.tw.fill_branches(dpsjgh = _d0s_gh)
		self.tw.fill_branches(minv = _df_d0['inv_mass'].values.tolist())
		self.tw.fill_tree()


def main():
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='single root file or a file with a list of files to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()

	if '.root' in args.flist:
		hfa = HFAIO()
		hfa.process_file(args.flist)


if __name__ == '__main__':
	main()