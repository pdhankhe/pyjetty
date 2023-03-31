"""
  Analysis IO class for jet analysis with HF and track dataframe.
  Each instance of the class handles the IO of a *single* track tree.
  Authors: Mateusz Ploskon
           Preeti Dhankher
"""


from pyjetty.mputils.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
from copy import deepcopy
import pandas as pd
import os
import tqdm
import numpy as np
from pyjetty.alihfjets.dev.hfjet.process.base.bitwise import filter_bit_df, tag_bit_df
# from numba import jit
# import numexpr

class DFSelection(MPBase):
	def __init__(self, **kwargs):
		super(DFSelection, self).__init__(**kwargs)
		self.selection = []
		self.df_selection = None
		self.query_strings = []
		self.query_string = ''


	def add_selection_equal(self, what, val):
		self.selection.append([what, val, None, 0])
		self.query_strings.append('({} == {})'.format(what, val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range(self, what, minv, maxv):
		self.selection.append([what, minv, maxv, 1])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, minv, what, maxv))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range_abs(self, what, val):
		self.selection.append([what, val, None, 2])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, -val, what, +val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_cond(self, what, val1, val2):
		self.selection.append([what, val1, val2, 3])
		self.query_strings.append('((({} > {}) & ({} < {})) | ({} < {}))'.format(what, -val2, what, +val2, what, val1))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'


	def add_selection_nsig(self, tpcp0, tofp0, tpck1, tofk1, tpcp1, tofp1, tpck0, tofk0, val1, val2):
	#	self.selection.append([tpcp0, tofp0, None, 4])
		self.query_strings.append('((abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})) | (abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})))'.format(tpcp0, val1, tofp0, val1, tofp0, val2, tpck1, val1, tofk1, val1, tofk1, val2, tpcp1, val1, tofp1, val1, tofp1, val2, tpck0, val1, tofk0, val1, tofk0, val2))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def compile_selection(self, df):
		self.df_selection = True
		for c in self.selection:
			if c[3] == 0:
				self.df_selection = (self.df_selection) & (df[c[0]] == c[1])
			if c[3] == 1:
				self.df_selection = (self.df_selection) & (df[c[0]] > c[1]) & (df[c[0]] < c[2])
			if c[3] == 2:
				self.df_selection = (self.df_selection) & (df[c[0]] > -c[1]) & (df[c[0]] < c[1])


class HFAnalysis(MPBase):
	def __init__(self, **kwargs):
		super(HFAnalysis, self).__init__(**kwargs)
		self.callback = None
		self.event_selection = DFSelection()
		self.d0_selection = DFSelection()
		self.d0_gen_selection = DFSelection()
	
	def analyze(self, df):
		self.compile_selection(df)
		_df = df[self.df_selection]
		#self.analysis(_df)
		if self.callback is not None:
			self.callback(df['ev_id'].values[0])

	def exec_analysis_event_df(self, df):
		self.pbar.update(1)
		_n_d0s = len(df)
		if _n_d0s < 1:
			return
		self.process_events(df)


	def exec_analysis_d0_df(self, df,isMC):
		self.pbar.update(1)
		#print(df)
		_n_d0s = len(df)
		
		if _n_d0s < 1:
			return
		if 'ev_id_ext' in list(self.event_df):
			_ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'].values[0], df['ev_id'].values[0], df['ev_id_ext'].values[0])
		else:
			_ev_query = "run_number == {} & ev_id == {}".format(df['run_number'].values[0], df['ev_id'].values[0])
		if isMC:
			self.cand_identifier=self.unique_identifier+['cand_type','pt_cand','eta_cand','phi_cand']
			self.df_particles = self.track_gen_df
		else:
			self.cand_identifier=self.unique_identifier+['cand_type','pt_cand','eta_cand','phi_cand','inv_mass']
			self.df_particles = self.track_df

		self.jet_identifier=['jet_pt','jet_eta','jet_phi']	
		self.df_tracks = self.df_particles.query(_ev_query) 
		self.df_tracks.reset_index(drop=True)
		pd_rec= self.analysis(df,isMC)
		self.df_tracks = None
		return pd_rec
			
		
	def selectfidacc(self, arr_pt, arr_y):
		array_is_sel = []
		pt_cand=arr_pt.to_numpy()
		y_cand=arr_y.to_numpy()

		for icand, pt in enumerate(pt_cand):
			if pt > 5:
				if abs(y_cand[icand]) < 0.8:
					array_is_sel.append(True)
				else:
					array_is_sel.append(False)
			else:
				yfid = -0.2/15 * pt**2 + 1.9/15 * pt + 0.5
				if abs(y_cand[icand]) < yfid:
					array_is_sel.append(True)
				else:
					array_is_sel.append(False)
		return array_is_sel


	def special_cuts_np(self, arr_pt_cand, arr_eta_cand, arr_imp_par_prong0, arr_imp_par_prong1,arr_imp_par_err_prong0, arr_imp_par_err_prong1, arr_d_len, arr_norm_dl):
                #implementation of AliRDHFCutsD0toKpi::IsSelectedSpecialCuts
                pt_cand = arr_pt_cand.to_numpy()
                eta_cand = arr_eta_cand.to_numpy()
                imp_par_prong0 = arr_imp_par_prong0.to_numpy()
                imp_par_prong1 = arr_imp_par_prong1.to_numpy()
                imp_par_err_prong0 = arr_imp_par_err_prong0.to_numpy()
                imp_par_err_prong1 = arr_imp_par_err_prong1.to_numpy()
                d_len = arr_d_len.to_numpy()
                norm_dl = arr_norm_dl.to_numpy()

                # normalised prong impact parameter
                normd0Cut = 0.5
                d0norm_prong_0 = np.divide(imp_par_prong0, imp_par_err_prong0)
                d0norm_prong_1 = np.divide(imp_par_prong1, imp_par_err_prong1)
                sel_d0norm_prong_0 = np.greater_equal(np.abs(d0norm_prong_0), normd0Cut)
                sel_d0norm_prong_1 = np.greater_equal(np.abs(d0norm_prong_1), normd0Cut)

                # decay length
                p_cand = np.multiply(pt_cand, np.cosh(eta_cand))
                decLengthCut = np.minimum(p_cand * 0.0066 + 0.01, 0.06)
                sel_d_len = np.greater_equal(d_len, decLengthCut)
                # normalised decay length
                normDecLengthCut = 1.
                sel_norm_dlen = np.greater_equal(np.abs(norm_dl), normDecLengthCut)


                sel_tot = np.logical_and(sel_d0norm_prong_0, np.logical_and(sel_d0norm_prong_1,np.logical_and(sel_d_len, sel_norm_dlen)))

                return sel_tot

	

	def apply_cut_fiducial_acceptance(self, df_):
                df_["is_sel_special"] = self.selectfidacc(df_["pt_cand"], df_["y_cand"])
                df_sel = df_[df_["is_sel_special"] == 1]
                return df_sel

	def apply_cut_special_np(self, df_):
                df_["is_sel_special"] = self.special_cuts_np(df_["pt_cand"], df_["eta_cand"],\
                df_["imp_par_prong0"], df_["imp_par_prong1"], df_["imp_par_err_prong0"], \
                df_["imp_par_err_prong1"], df_["d_len"], df_["norm_dl"])
                df_sel = df_[df_["is_sel_special"] == 1]
                return df_sel

	def apply_cuts_ptbin(self, df_):
		self.ptcand_binmin=self.config["analysisD0"]["pt_cand_binning"].get("sel_an_binmin")
		self.ptcand_binmax=self.config["analysisD0"]["pt_cand_binning"].get("sel_an_binmin")
		self.npt_cand_bins = len(self.ptcand_binmin)
		pt_cand_cuts = self.config["analysisD0"].get("cuts")
		self.analysis_cuts = deepcopy(pt_cand_cuts)
		print(self.ptcand_binmin)
		d0ev_df_custom_cuts=[]
		for ipt in range(self.npt_cand_bins):
			print(ipt)
			d0ev_df_custom_cuts.append(df_.query(self.analysis_cuts[ipt]))
		d0ev_df_custom_cuts = pd.concat(d0ev_df_custom_cuts)
		return(d0ev_df_custom_cuts)


	def analyze_slower(self):
		isMC = False
		#process total number of events
		_tot_event_df = self.event_df.copy()
		_tot_event_df.query(self.event_selection.query_string, inplace=True)
		_tot_event_df_ev_grouped= _tot_event_df.groupby(['run_number','ev_id'])
		#with tqdm.tqdm(total=len(_tot_event_df_ev_grouped)) as self.pbar:
		#	_tmp = _tot_event_df_ev_grouped.apply(self.exec_analysis_event_df)
		#self.pbar.close()
		
		print("number of events = " , len(_tot_event_df_ev_grouped))
		self.hNevents.Fill(1,len(_tot_event_df_ev_grouped))

		self.b_mcsigprompt = self.config["analysisD0"]["bitmap_sel"]["ismcprompt"]
		self.b_mcsigfd =self.config["analysisD0"]["bitmap_sel"]["ismcfd"]
		self.b_mcbkg = self.config["analysisD0"]["bitmap_sel"]["ismcbkg"]
		self.b_mcrefl = self.config["analysisD0"]["bitmap_sel"]["ismcrefl"]

		self.v_bitvar = self.config["analysisD0"]["bitmap_sel"]["var_name"]
		self.v_isstd = self.config["analysisD0"]["bitmap_sel"]["var_isstd"]
		self.v_ismcsignal = self.config["analysisD0"]["bitmap_sel"]["var_ismcsignal"]
		self.v_ismcprompt = self.config["analysisD0"]["bitmap_sel"]["var_ismcprompt"]
		self.v_ismcfd = self.config["analysisD0"]["bitmap_sel"]["var_ismcfd"]
		self.v_ismcbkg = self.config["analysisD0"]["bitmap_sel"]["var_ismcbkg"]
		self.v_ismcrefl = self.config["analysisD0"]["bitmap_sel"]["var_ismcrefl"]


		# generated D0 candidate dataframe
		d0_gen_df_copy=self.d0_gen_df.query(self.d0_gen_selection.query_string,engine="python")
		#reconstructed D0 candidate dataframe after applying D0 selection cuts
		_d0_df = self.d0_df.query(self.d0_selection.query_string,engine="python")
		
		self.unique_identifier =  ['run_number', 'ev_id']

		if 'ev_id_ext' in list(self.event_df):
			self.unique_identifier += ['ev_id_ext']

		d0ev_df = pd.merge(_d0_df, self.event_df,on=self.unique_identifier)
		d0ev_gen_df = pd.merge(d0_gen_df_copy, self.event_df,on=self.unique_identifier)
		

		#after merging the dataframe with event, apply the event selection cut on z vertex, event rejected..
		d0ev_df.query(self.event_selection.query_string, inplace=True)
		d0ev_gen_df.query(self.event_selection.query_string, inplace=True)

		pinfo('N generated d0s with d0 selection cuts and event cuts', len(d0ev_gen_df.index))
	
		#need to remove events id's at reconstructed level not present at generated level
		#coming from fake D0 due to looser selection cuts 	

		d0ev_df.sort_values(by=self.unique_identifier, inplace=True)
		d0ev_gen_df.sort_values(by=self.unique_identifier, inplace=True)
		

		df_d0runs = d0ev_df[self.unique_identifier].copy()
		df_gend0runs = d0ev_gen_df[self.unique_identifier].copy()

		#find reconstructed event matching to generator
		df_runs = pd.merge(df_d0runs, df_gend0runs, on=self.unique_identifier)
		df_runs.drop_duplicates(keep='first', inplace=True)

		d0ev_df = pd.merge(d0ev_df, df_runs, on=self.unique_identifier) 
		
		#print(df_runs)
		#apply fiducial cut on D candidate
		d0ev_df =self.apply_cut_fiducial_acceptance(d0ev_df)
		#apply special cuts for low pt
		d0ev_df = self.apply_cut_special_np(d0ev_df)
		#apply pt dependent custom selection cuts
		d0ev_df=self.apply_cuts_ptbin(d0ev_df)

		#apply fiducial cut on gen D candidate
		d0ev_gen_df=self.apply_cut_fiducial_acceptance(d0ev_gen_df)
		

		d0ev_df_grouped = d0ev_df.groupby(self.unique_identifier)
		d0ev_gen_df_grouped = d0ev_gen_df.groupby(self.unique_identifier)
	
		#print(d0ev_gen_df.head(6))
		pinfo('d0s after event cuts ', len(d0ev_df.index))
		pinfo('N d0 groups after event cuts ', len(d0ev_df_grouped))

		pinfo('generated d0s after event cuts ', len(d0ev_gen_df.index))
		pinfo('N generated d0 groups after event cuts ', len(d0ev_gen_df_grouped))

		#print(*d0ev_df_grouped)	
		with tqdm.tqdm(total=len(d0ev_df_grouped)) as self.pbar:
			_pd_reco=d0ev_df_grouped.apply(self.exec_analysis_d0_df,False)
			_pd_reco=_pd_reco.dropna()
			_pd_reco=_pd_reco.reset_index(drop=True)
			_pd_reco=_pd_reco.set_index(self.unique_identifier)
			_pd_reco[self.v_ismcfd] = np.array(tag_bit_df(_pd_reco, self.v_bitvar,self.b_mcsigfd), dtype=int)
			_pd_reco[self.v_ismcprompt] = np.array(tag_bit_df(_pd_reco, self.v_bitvar,self.b_mcsigprompt), dtype=int)
			_pd_reco[self.v_ismcrefl] = np.array(tag_bit_df(_pd_reco, self.v_bitvar,self.b_mcrefl), dtype=int)	
			print(_pd_reco)
		self.pbar.close()
		
		#d0ev_gen_df_grouped.apply(self.exec_analysis_d0_gen_df)
		with tqdm.tqdm(total=len(d0ev_gen_df_grouped)) as self.pbar:
			_pd_gen = d0ev_gen_df_grouped.apply(self.exec_analysis_d0_df,True)
			_pd_gen=_pd_gen.dropna()
			_pd_gen=_pd_gen.reset_index(drop=True)
			_pd_gen=_pd_gen.set_index(self.unique_identifier)	
			_pd_gen[self.v_ismcprompt] = np.array(tag_bit_df(_pd_gen, self.v_bitvar,self.b_mcsigprompt), dtype=int)
			_pd_gen[self.v_ismcfd] = np.array(tag_bit_df(_pd_gen, self.v_bitvar,self.b_mcsigfd), dtype=int)
			_pd_gen[self.v_ismcrefl] = np.array(tag_bit_df(_pd_gen, self.v_bitvar,self.b_mcrefl), dtype=int)
			#print(_pd_gen)
		self.pbar.close()
	
		self.EffandResponse(_pd_reco,_pd_gen)	
	
	# analysis on the single data frame
	# this is something specific to user - overload this one

	def analysis(self, df):
		if len(df) > 0:
			print (df)


class HFAnalysisIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0',
					d0_gen_tree_name='PWGHF_TreeCreator/tree_D0_gen',
					track_tree_name='PWGHF_TreeCreator/tree_Particle',
					track_gen_tree_name='PWGHF_TreeCreator/tree_Particle_gen',
					event_tree_name='PWGHF_TreeCreator/tree_event_char')
		super(HFAnalysisIO, self).__init__(**kwargs)
		self.analyses = []
		self.d0_df_grouped = None
		self.d0_gen_df_grouped = None

	def reset_analyses_list(self):
		self.analyses = []

	def add_analysis(self, a):
		self.analyses.append(a)

	def pd_tree(self, path, tname, squery=None):
		try:
			tree = uproot.open(path)[tname]
		except:
			pwarning('error getting', tname, 'from file:', path)
			return None
		if not tree:
			perror('Tree {} not found in file {}'.format(tname, path))
			return None
		df = uproot.concatenate(tree, library="pd")

		if squery:
			#df.query(squery, inplace=True)
			df = df.query(squery)
			df.reset_index(drop=True)
		return df

	def load_file(self, path):
		if not os.path.exists(path):
			pwarning('[w] file', path, 'does not exists.')
			return

		self.event_df_orig = self.pd_tree(path, self.event_tree_name)
		#select only run number, event id, zvtx and event rejected info
		self.event_df =self.event_df_orig[['run_number', 'ev_id','ev_id_ext','z_vtx_reco','is_ev_rej']].copy() 
		if self.event_df_orig is None:
			return False
		self.event_df.reset_index(drop=True)


		self.d0_df = self.pd_tree(path, self.d0_tree_name)
		if self.d0_df is None:
			return False

		self.track_df = self.pd_tree(path, self.track_tree_name)
		if self.track_df is None:
			return False

		self.d0_gen_df = self.pd_tree(path, self.d0_gen_tree_name)
		if self.d0_gen_df is None:
			return False

		self.track_gen_df = self.pd_tree(path, self.track_gen_tree_name)
		if self.track_gen_df is None:
			return False

		return True

	def execute_analyses(self):
		for ana in self.analyses:
			# update_df_references
			ana.event_df 	= self.event_df
			ana.d0_df 		= self.d0_df
			ana.track_df 	= self.track_df

			ana.d0_gen_df   = self.d0_gen_df
			ana.track_gen_df     = self.track_gen_df
			ana.analyze_slower()
			ana.event_df 	= None
			ana.d0_df 		= None
			ana.track_df 	= None
			ana.d0_gen_df         = None
			ana.track_gen_df     = None
		# drop the dfs
		self.event_df = None
		self.d0_df = None
		self.track_df = None
		self.d0_gen_df = None
		self.track_gen_df = None

	def update_status(self, mark):
		if mark != self.pbar2_mark:
			self.pbar2_mark = mark
			self.pbar2.update(1)

	def execute_analyses_on_inputfile(self, inputFile):
		print()
		if os.path.exists(inputFile):
			pinfo('file:', inputFile)
			if self.load_file(inputFile):
				self.execute_analyses()
		else:
			perror('file does not exist', inputFile)
		pinfo('done.')

	def __def__(self):
		self.d0_df_grouped = None
		self.d0_gen_df_grouped = None
