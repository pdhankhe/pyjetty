from pyjetty.mputils.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
from copy import deepcopy
import pandas as pd
import os
import fjext
import fjtools
import fjcontrib
import tqdm
import math as ma
import numpy as np
import fastjet as fj
from pyjetty.alihfjets.dev.hfjet.process.base.bitwise  import filter_bit_df, tag_bit_df
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
        
    def exec_analysis_d0_df(self, df, m, random_mass, min_pt, IsGen):
        self.pbar.update(1)
        _n_d0s = len(df)
        if _n_d0s < 1:
            return
   
        if 'ev_id_ext' in list(self.event_df):
            _ev_query = "run_number == {} & ev_id == {} & ev_id_ext == {}".format(df['run_number'].values[0], df['ev_id'].values[0], df['ev_id_ext'].values[0])
        else:
            _ev_query = "run_number == {} & ev_id == {}".format(df['run_number'].values[0], df['ev_id'].values[0])
            
        if IsGen:
            self.cand_identifier=self.unique_identifier+['cand_type','pt_cand','eta_cand','phi_cand']
            self.df_particles = self.track_gen_df
        else:
            self.cand_identifier=self.unique_identifier+['cand_type','pt_cand','eta_cand','phi_cand','inv_mass']
            self.df_particles = self.track_df

        self.df_tracks = self.df_particles.query(_ev_query)
        self.df_tracks.reset_index(drop=True)
        #print(self.df_tracks)
    
        df_tracks_accepted = self.df_tracks[self.df_tracks.ParticlePt > min_pt]
        
        m_array = np.full((df_tracks_accepted['ParticlePt'].values.size), m)
        
        if random_mass:
            rand_val = np.random.random((len(m_array)))
            K_mass = 0.4937     # GeV/c^2
            p_mass = 0.938272   # GeV/c^2
            # (p + pbar) / (pi+ + pi-) ~ 5.5%
            # (K+ + K-) / (pi+ + pi-) ~ 13%
            # But these are numbers with respect to the final _unreplaced_ pions, so there is
            # an additional factor of 1/(1 + 5.5% + 13%) to get things right
            K_factor = 0.13; p_factor = 0.055
            K_prob = K_factor / (1 + K_factor + p_factor)
            p_prob = 1 - p_factor / (1 + K_factor + p_factor)   # 1- just to look at diff random vals
            m_array = np.where(rand_val < K_prob, K_mass, m_array)
            m_array = np.where(rand_val > p_prob, p_mass, m_array)

        
        fj_particles = fjext.vectorize_pt_eta_phi_m(df_tracks_accepted['ParticlePt'].values, df_tracks_accepted['ParticleEta'].values,df_tracks_accepted['ParticlePhi'].values, m_array)
        self.df_tracks = None
        self.df_tracks_accepted = None
        
        return fj_particles
    
    def replace_daughter_group_fj_particles(self, df, IsGen):
        djmm = fjtools.DJetMatchMaker()
        
        if IsGen:
            m_cand_gen_array = np.full((df['pt_cand'].values.size), 1.864)
            djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values,m_cand_gen_array)
          
        else:
            djmm.set_Ds_pt_eta_phi_m(df['pt_cand'].values, df['eta_cand'].values, df['phi_cand'].values, df['inv_mass'].values)
            djmm.set_daughters0_pt_eta_phi(df['pt_prong0'].values, df['eta_prong0'].values, df['phi_prong0'].values)
            djmm.set_daughters1_pt_eta_phi(df['pt_prong1'].values, df['eta_prong1'].values, df['phi_prong1'].values)
            
        
        _D0_fj_particles = df[['run_number', 'ev_id','ev_id_ext', 'ev_id_long', 'ismcfd', 'ismcprompt', 'ismcrefl']].copy()
        
        array_index_df= _D0_fj_particles.index.values
     
        _D0_fj_particles['fj_D0_particles'] = _D0_fj_particles.apply(lambda _: '', axis=1)
        #_D0_fj_particles['fj_new_particles'] = _D0_fj_particles['fj_new_particles'].astype('object')
        for id0, d0 in enumerate(djmm.Ds):
            #replacing daughter tracks with matched D0 candidate
            djmm.ch = df['fj_particles'].iloc[id0]
            if IsGen:
                _parts_and_ds=djmm.ch
                for i in range(0,len(_parts_and_ds)):
                    for j in range(i+1,len(_parts_and_ds)):
                        daughtersSum=_parts_and_ds[i]+_parts_and_ds[j]
                        diff=daughtersSum-d0
                        if(ma.sqrt((diff.px()*diff.px()))<0.001 and ma.sqrt((diff.py()*diff.py()))<0.001 and ma.sqrt((diff.pz()*diff.pz()))<0.001):
                            _parts_and_ds[i]=_parts_and_ds[i]* 1.e-6
                            _parts_and_ds[j]=_parts_and_ds[j]* 1.e-6
             
        
            else:
                _parts_and_ds = djmm.ch
                _parts_and_ds = djmm.match(0.005, id0)
    
            #including D0
            _parts_and_ds.push_back(d0)
            # converting list to std fj vector
            fj_particles = fj.vectorPJ()
            for index in range(len(_parts_and_ds)):
                fj_particles.append(_parts_and_ds[index])
            _D0_fj_particles.at[id0, 'fj_D0_particles'] = fj_particles
                
        return _D0_fj_particles
        

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
            #print(ipt)
            d0ev_df_custom_cuts.append(df_.query(self.analysis_cuts[ipt]))
        d0ev_df_custom_cuts = pd.concat(d0ev_df_custom_cuts)
        return(d0ev_df_custom_cuts)


    def analyze_d0_df(self, m, offset_indices,
                group_by_evid, random_mass, min_pt, IsGen):

        #process total number of events
        _tot_event_df = self.event_df.copy()
        _tot_event_df.query(self.event_selection.query_string, inplace=True)
        _tot_event_df_ev_grouped= _tot_event_df.groupby(['run_number','ev_id'])
        print("number of events = " , len(_tot_event_df_ev_grouped))

        ### configuration for trimming dataframe
        
        self.IsFiducialCut = self.config["IsFiducialCut"]
        self.IsCustomSelectionCuts = self.config["IsCustomSelectionCuts"]
        self.IsSpecialLowpTCuts = self.config["IsSpecialLowpTCuts"]
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
        ############################################################################
        
        
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
        if self.IsFiducialCut:
            print("aplying fiducial cuts")
            d0ev_df =self.apply_cut_fiducial_acceptance(d0ev_df)
            
        #apply special cuts for low pt
        if self.IsSpecialLowpTCuts:
            print("aplying special cuts")
            d0ev_df = self.apply_cut_special_np(d0ev_df)
        
        #apply pt dependent custom selection cuts
        if self.IsCustomSelectionCuts:
            print("aplying pt dependent cuts")
            d0ev_df=self.apply_cuts_ptbin(d0ev_df)
            
        #apply fiducial cut on gen D candidate
        if self.IsFiducialCut:
            self.d0ev_gen_df=self.apply_cut_fiducial_acceptance(d0ev_gen_df)
        
        # tag prompt, fd and reflection signal
        d0_ev_tagged_df = self.tagging_prompt_fd_reflection(d0ev_df, IsGen)
        d0_ev_gen_tagged_df = self.tagging_prompt_fd_reflection(self.d0ev_gen_df, True)
    
        if IsGen:
            df_fjparticles = self.group_df_fjparticles(d0_ev_gen_tagged_df, group_by_evid, m, random_mass, min_pt, IsGen)
        else:
            df_fjparticles = self.group_df_fjparticles(d0_ev_tagged_df, group_by_evid, m, random_mass, min_pt, IsGen)
        
        d0_ev_tagged_df = None
        d0_ev_gen_tagged_df = None
        
        return df_fjparticles

   
        
    def group_df_fjparticles(self, D0_df, group_by_evid, m, random_mass, min_pt, IsGen):
    
        d0ev_df_grouped = None
        d0ev_df_grouped = D0_df.groupby(self.unique_identifier)

        if len(D0_df) > 0: #checking dataframe to not be empty
            # reconstructed level
            with tqdm.tqdm(total=len(d0ev_df_grouped)) as self.pbar:
                fj_particles_lists = d0ev_df_grouped.apply(self.exec_analysis_d0_df, m, random_mass,  min_pt, IsGen)
                fj_particles_lists=fj_particles_lists.dropna()
            self.pbar.close() 
        
            df_fj_particle = fj_particles_lists.to_frame()
            df_fj_particle.columns = ['fj_particles']
        
            #merge D0 and fj particles into single dataframe
            df_fjparticles = pd.merge(D0_df , df_fj_particle , on=self.unique_identifier)

            df_D0_daughter_fjparticles = self.replace_daughter_group_fj_particles(df_fjparticles, IsGen)
       
            #rename generated tages, needed fot matching the candidate
            if IsGen:
                df_D0_daughter_fjparticles.rename(columns = {'ismcfd':'ismcfd_gen', 'ismcprompt':'ismcprompt_gen', 'ismcrefl':'ismcrefl_gen', 'fj_D0_particles' : 'fj_D0_particles_gen' }, inplace = True)

            #tagging with event label
        
            df_fjparticle_reindex = df_D0_daughter_fjparticles.set_index(self.unique_identifier)
       
            #print(df_fjparticle_reindex)
    
            return df_fjparticle_reindex
    
        else: #return none
            pinfo('Number of D0 in the event is = ', len(D0_df.index))
            return d0ev_df_grouped


    def tagging_prompt_fd_reflection(self, _df, isGen):
        # save only relevant column
        if isGen:
            _df =_df[['run_number', 'ev_id','ev_id_ext','ev_id_long','cand_type','pt_cand','eta_cand','phi_cand']].copy()
        else:
            _df =_df[['run_number', 'ev_id','ev_id_ext','ev_id_long','cand_type','pt_cand','eta_cand','phi_cand', 'inv_mass', 'pt_prong0', 'eta_prong0', 'phi_prong0', 'pt_prong1', 'eta_prong1', 'phi_prong1']].copy()
        
        # tag bits
        _df[self.v_ismcfd] = np.array(tag_bit_df(_df, self.v_bitvar, self.b_mcsigfd), dtype=int)
        _df[self.v_ismcprompt] = np.array(tag_bit_df(_df, self.v_bitvar, self.b_mcsigprompt), dtype=int)
        _df[self.v_ismcrefl] = np.array(tag_bit_df(_df, self.v_bitvar, self.b_mcrefl), dtype=int)
        
        return _df

    def analysis(self, df):
        if len(df) > 0:
            print (df)


class HFAnalysisIO(MPBase):
    def __init__(self, input_file, **kwargs):
        self.configure_from_args(d0_tree_name='PWGHF_TreeCreator/tree_D0',
                    d0_gen_tree_name='PWGHF_TreeCreator/tree_D0_gen',
                    track_tree_name='PWGHF_TreeCreator/tree_Particle',
                    track_gen_tree_name='PWGHF_TreeCreator/tree_Particle_gen',
                    event_tree_name='PWGHF_TreeCreator/tree_event_char')
        super(HFAnalysisIO, self).__init__(**kwargs)
        self.analyses = []
        self.d0_df_grouped = None
        self.d0_gen_df_grouped = None
        self.input_file = input_file
        
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

    def load_file(self, input_file, Is_treff_sys, reject_tracks_fraction):
        if not os.path.exists(input_file):
            pwarning('[w] file', input_file, 'does not exists.')
            return

        self.event_df_orig = self.pd_tree(input_file, self.event_tree_name)
        #select only run number, event id, zvtx and event rejected info
        self.event_df =self.event_df_orig[['run_number', 'ev_id','ev_id_ext','z_vtx_reco','is_ev_rej']].copy()
        if self.event_df_orig is None:
            return False
        self.event_df.reset_index(drop=True)


        self.d0_df = self.pd_tree(input_file, self.d0_tree_name)
        if self.d0_df is None:
            return False

        self.track_df = self.pd_tree(input_file, self.track_tree_name)
        if self.track_df is None:
            return False

        print("number of particles before throwing out = "+str(len(self.track_df)))
       
        if Is_treff_sys:
            reject_tracks_fraction = 0.03
            n_remove = int(reject_tracks_fraction * len(self.track_df.index))
            np.random.seed()
            indices_remove = np.random.choice(self.track_df.index, n_remove, replace=False)
            self.track_df.drop(indices_remove, inplace=True)

        print("number of particles after throwing out = "+str(len(self.track_df)))

        self.d0_gen_df = self.pd_tree(input_file, self.d0_gen_tree_name)
        if self.d0_gen_df is None:
            return False

        self.track_gen_df = self.pd_tree(input_file, self.track_gen_tree_name)
        if self.track_gen_df is None:
            return False

        return True
                
    def execute_analyses(self, m, reject_tracks_fraction, offset_indices,
                group_by_evid, random_mass, min_pt, IsGen):
        for ana in self.analyses:
            # update_df_references
            ana.event_df 	= self.event_df
            ana.d0_df 		= self.d0_df
            ana.track_df 	= self.track_df

            ana.d0_gen_df   = self.d0_gen_df
            ana.track_gen_df = self.track_gen_df
            df_fjparticles = ana.analyze_d0_df(m, offset_indices,
                group_by_evid, random_mass, min_pt, IsGen)
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
        return df_fjparticles
    

    def update_status(self, mark):
        if mark != self.pbar2_mark:
            self.pbar2_mark = mark
            self.pbar2.update(1)

  #---------------------------------------------------------------
  # Convert ROOT TTree to SeriesGroupBy object of fastjet particles per event.
  # Optionally, define the mass assumption used in the jet reconstruction;
  #             remove a certain random fraction of tracks;
  #             randomly assign proton and kaon mass to some tracks
  #---------------------------------------------------------------
  
    def load_data(self, m=0.1396, offset_indices=False,
                group_by_evid=True, random_mass=False, min_pt=0., Is_treff_sys  = False, reject_tracks_fraction=0., IsGen=False):

        if os.path.exists(self.input_file):
            if self.load_file(self.input_file, Is_treff_sys , reject_tracks_fraction):
                df_fjparticles = self.execute_analyses(m, reject_tracks_fraction, offset_indices,
                group_by_evid, random_mass, min_pt, IsGen)
                return df_fjparticles
        else:
            perror('file list does not exist', self.input_file)
        pinfo('done.')

    def __def__(self):
        self.d0_df_grouped = None
        self.d0_gen_df_grouped = None
