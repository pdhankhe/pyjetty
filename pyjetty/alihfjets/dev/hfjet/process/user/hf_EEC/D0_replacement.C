class D0_reco
{
  private:
    // evt info
    TVector3 evt_vtx_true;
    TVector3 evt_vtx_reco;

    double x_true;
    double Q2_true;
    double nu_true;

    TLorentzVector hadron_beam;
    TLorentzVector struck_quark;

    // negative charged track true info
    vector<TLorentzVector> negl_p_true;
    vector<TVector3> negl_vtx_true;
    vector<int> negl_id_true;
    vector<int> negl_parent_index_true;
    vector<int> negl_parent_id_true;
    vector<TLorentzVector> negl_quark_p;

    // positive charged track true info
    vector<TLorentzVector> posl_p_true;
    vector<TVector3> posl_vtx_true;
    vector<int> posl_id_true;
    vector<int> posl_parent_index_true;
    vector<int> posl_parent_id_true;
    vector<TLorentzVector> posl_quark_p;

    // negative charged track reco info
    vector<TLorentzVector> negl_p_reco;
    vector<TVector3> negl_vtx_reco;
    vector<float> negl_prob_e;
    vector<float> negl_prob_pi;
    vector<float> negl_prob_K;
    vector<float> negl_prob_p;

    // positive charged track reco info
    vector<TLorentzVector> posl_p_reco;
    vector<TVector3> posl_vtx_reco;
    vector<float> posl_prob_e;
    vector<float> posl_prob_pi;
    vector<float> posl_prob_K;
    vector<float> posl_prob_p;

    // single track cut info
    float TRK_P_LO;
    float TRK_DCA;
    int ID_OPTION; // no PID, PID w/o assumption, PID w/ assumption
    int SMEAR_OPTION;
    int BFIELD_TYPE;

    // pair cut info
    float PAIR_DCA;
    float DECAY_L;
    float D0_DCA;
    float D0_COSTHETA;

    // if D0 vertex issue needs to be corrected, 1 if yes, 0 if no
    int do_correct_vertex;

    // mass pair vs pt
    TH2D* fg2d_Kpimass_vs_p[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpimass_vs_p[chargebin][etabin]; // 0: K-pi+

    // mass pair vs z
    TH2D* fg2d_Kpimass_vs_z[chargebin][etabin]; // 0: K-pi+
    TH2D* bg2d_Kpimass_vs_z[chargebin][etabin]; // 0: K-pi+

    TH2D* h2d_K_D0_p_vs_eta[etabin][pptbin];
    TH2D* h2d_pi_D0_p_vs_eta[etabin][pptbin];

    TH2D* h2d_D0_pt_vs_eta[Q2bin][xbin];
    TH2D* h2d_D0_z_vs_eta[Q2bin][xbin];

    TH2D* h2d_D0_pt_vs_eta_gen[Q2bin][xbin];
    TH2D* h2d_D0_z_vs_eta_gen[Q2bin][xbin];

    TH2D* h2d_ztheo_vs_zjet;

  public:
    D0_reco()
    {
      x_true = 1E1; // unphysical
      Q2_true = 1E-5; // out of range
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0; // 0--no smearing, 1--DM smearing, 1--LBL smearing, 2--Hybrid smearing
      BFIELD_TYPE = 0; // 0--Barbar, 1--Beast

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      D0_DCA = -9999;
      D0_COSTHETA = -9999;

      do_correct_vertex = 0; // corrects D0 verticies

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta] = new TH2D(Form("fg2d_Kpimass_vs_p_%d_%d",icharge,ieta),"Kpi mass vs Kpi p",480,0.6,3,100,0,10);
          fg2d_Kpimass_vs_p[icharge][ieta]->Sumw2();

          bg2d_Kpimass_vs_p[icharge][ieta] = new TH2D(Form("bg2d_Kpimass_vs_p_%d_%d",icharge,ieta),"Kpi mass vs Kpi p",480,0.6,3,100,0,10);
          bg2d_Kpimass_vs_p[icharge][ieta]->Sumw2();

          fg2d_Kpimass_vs_z[icharge][ieta] = new TH2D(Form("fg2d_Kpimass_vs_z_%d_%d",icharge,ieta),"Kpi mass vs Kpi z",480,0.6,3,100,0,1);
          fg2d_Kpimass_vs_z[icharge][ieta]->Sumw2();

          bg2d_Kpimass_vs_z[icharge][ieta] = new TH2D(Form("bg2d_Kpimass_vs_z_%d_%d",icharge,ieta),"Kpi mass vs Kpi z",480,0.6,3,100,0,1);
          bg2d_Kpimass_vs_z[icharge][ieta]->Sumw2();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_K_D0_p_vs_eta_%d_%d",ieta,ipt),"K p vs eta",100,0,10,160,-4,4);
          h2d_K_D0_p_vs_eta[ieta][ipt]->Sumw2();
          h2d_pi_D0_p_vs_eta[ieta][ipt] = new TH2D(Form("h2d_pi_D0_p_vs_eta_%d_%d",ieta,ipt),"pi p vs eta",100,0,10,160,-4,4);
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Sumw2();
        }
      }


      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix] = new TH2D(Form("h2d_D0_pt_vs_eta_%d_%d",iQ2,ix),"D0 pt vs eta",100,0,10,160,-4,4);
          h2d_D0_pt_vs_eta[iQ2][ix]->Sumw2();

          h2d_D0_z_vs_eta[iQ2][ix] = new TH2D(Form("h2d_D0_z_vs_eta_%d_%d",iQ2,ix),"D0 z vs eta",100,0,1,160,-4,4);
          h2d_D0_z_vs_eta[iQ2][ix]->Sumw2();

          h2d_D0_pt_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_D0_pt_vs_eta_gen_%d_%d",iQ2,ix),"D0 pt vs eta",100,0,10,160,-4,4);
          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Sumw2();

          h2d_D0_z_vs_eta_gen[iQ2][ix] = new TH2D(Form("h2d_D0_z_vs_eta_gen_%d_%d",iQ2,ix),"D0 z vs eta",100,0,1,160,-4,4);
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Sumw2();
        }
      }

      h2d_ztheo_vs_zjet = new TH2D("h2d_z_frag","z_{theo} vs z_{jet}",100,0,1,100,0,1);
      h2d_ztheo_vs_zjet->Sumw2();
    }
    virtual ~D0_reco()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          delete fg2d_Kpimass_vs_p[icharge][ieta];
          delete bg2d_Kpimass_vs_p[icharge][ieta];

          delete fg2d_Kpimass_vs_z[icharge][ieta];
          delete bg2d_Kpimass_vs_z[icharge][ieta];
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          delete h2d_K_D0_p_vs_eta[ieta][ipt];
          delete h2d_pi_D0_p_vs_eta[ieta][ipt];
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          delete h2d_D0_pt_vs_eta[iQ2][ix];
          delete h2d_D0_z_vs_eta[iQ2][ix];

          delete h2d_D0_pt_vs_eta_gen[iQ2][ix];
          delete h2d_D0_z_vs_eta_gen[iQ2][ix];
        }
      }

      delete h2d_ztheo_vs_zjet;
    };

    void Reset()
    {
      x_true = 1E1;
      Q2_true = 1E-5;
      nu_true = -9999;

      TRK_P_LO = -9999;
      TRK_DCA = -9999;
      ID_OPTION = -1; // -1--no PID, 0--DM PID, 1--PID w/o assumption, 2--PID w/ assumption
      SMEAR_OPTION = 0;
      BFIELD_TYPE = 0;

      PAIR_DCA = -9999; // 120um in unit of mm
      DECAY_L = -9999;
      D0_DCA = -9999;
      D0_COSTHETA = -9999;

      do_correct_vertex = 0;

      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta]->Reset("ICESM");
          bg2d_Kpimass_vs_p[icharge][ieta]->Reset("ICESM");

          fg2d_Kpimass_vs_z[icharge][ieta]->Reset("ICESM");
          bg2d_Kpimass_vs_z[icharge][ieta]->Reset("ICESM");
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt]->Reset("ICESM");
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Reset("ICESM");
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix]->Reset("ICESM");
          h2d_D0_z_vs_eta[iQ2][ix]->Reset("ICESM");

          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Reset("ICESM");
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Reset("ICESM");
        }
      }

      h2d_ztheo_vs_zjet->Reset("ICESM");
    }

    void SetLowP(const float mom_thr) { TRK_P_LO = mom_thr; }

    void SetDCACuts()
    {
      // single
      TRK_DCA = -9999; //0.03; // in units of mm

      // pair
      PAIR_DCA = 0.13; //0.12; // 120um in unit of mm
      DECAY_L = 0.04;
      D0_DCA = -9999; //0.01;
      D0_COSTHETA = 0.98;
    }

    void SetIDCuts(const int id_opt) { ID_OPTION = id_opt; }

    void SetSmearType(const int smear_opt) { SMEAR_OPTION = smear_opt; }

    void SetBFieldType(const int b_opt) { BFIELD_TYPE = b_opt; }

    void SetVectTrue(TVector3& _evt_vtx_true) { evt_vtx_true = _evt_vtx_true; }
    void SetVectReco(TVector3& _evt_vtx_reco) { evt_vtx_reco = _evt_vtx_reco; }

    void SetQ2True(double _Q2_true) { Q2_true = _Q2_true; }
    void SetXTrue(double _x_true) { x_true = _x_true; }

    void SetNuTrue(double _nu_true) { nu_true = _nu_true; };

    void SetDoCorrectVertex(int _do_correct) { do_correct_vertex = _do_correct; };

    void ClearTracks()
    {
      negl_p_true.clear();
      negl_vtx_true.clear();
      negl_id_true.clear();
      negl_parent_index_true.clear();
      negl_parent_id_true.clear();
      negl_quark_p.clear();

      negl_p_reco.clear();
      negl_vtx_reco.clear();
      negl_prob_e.clear();
      negl_prob_pi.clear();
      negl_prob_K.clear();
      negl_prob_p.clear();

      posl_p_true.clear();
      posl_vtx_true.clear();
      posl_id_true.clear();
      posl_parent_index_true.clear();
      posl_parent_id_true.clear();
      posl_quark_p.clear();

      posl_p_reco.clear();
      posl_vtx_reco.clear();
      posl_prob_e.clear();
      posl_prob_pi.clear();
      posl_prob_K.clear();
      posl_prob_p.clear();
    }

    void FillGenKin(erhic::EventPythia* py_evt)
    {
      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (abs(part->Id())!=421) continue; // FIX ME: for now fill all the D0 (maybe fill the D0->Kpi later)
        // NB: however at generation level there shouldn't be kin difference for different decay channels (it can be different once including the detector acceptance)

        TLorentzVector D0_mom4_gen = part->Get4Vector();
        double frag_z = hadron_beam.Dot(D0_mom4_gen)/(nu_true*hadron_beam.M());

        if (verbosity>1) std::cout << "D0 with pt " << D0_mom4_gen.Pt() << " z " << frag_z << " eta " << D0_mom4_gen.PseudoRapidity() << std::endl;

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (verbosity>1) std::cout << "Q2, x " << Q2_true << ", " << x_true << " bin " << iQ2bin << ", " << ixbin << std::endl;

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_D0_pt_vs_eta_gen[iQ2bin][ixbin]->Fill(D0_mom4_gen.Pt(),D0_mom4_gen.PseudoRapidity());
          h2d_D0_z_vs_eta_gen[iQ2bin][ixbin]->Fill(frag_z,D0_mom4_gen.PseudoRapidity());

          h2d_D0_pt_vs_eta_gen[Q2bin-1][xbin-1]->Fill(D0_mom4_gen.Pt(),D0_mom4_gen.PseudoRapidity()); // the last bin holds all entries - fill this too
          h2d_D0_z_vs_eta_gen[Q2bin-1][xbin-1]->Fill(frag_z,D0_mom4_gen.PseudoRapidity());
        }
      }
    }

    void FillSingleTracks(erhic::EventPythia* py_evt)
    {
      ClearTracks();

      erhic::ParticleMC* proton = py_evt->GetTrack(1);
      if (proton!=NULL)
      {
        assert(abs(proton->Id())!=2212);
        hadron_beam = proton->Get4Vector();
      }
      else return; // if incoming proton not found, skip the whole event

      if (do_correct_vertex == 1) correct_D0_verticies(py_evt);

      for(int ipart = 0; ipart < py_evt->GetNTracks(); ipart++)
      {
        erhic::ParticleMC* part = py_evt->GetTrack(ipart);

        if (part->GetStatus()!=1) continue; // only loop through final stable particles

        if (abs(part->Id())!=2212 && abs(part->Id())!=211 && abs(part->Id())!=321 && abs(part->Id())!=11 && abs(part->Id())!=13) continue; // proton, pion, kaon, electron, muon

        if (TRK_P_LO>-99 && part->GetPt()<TRK_P_LO) continue;

        TLorentzVector track_mom4_true = part->Get4Vector();
        TLorentzVector track_mom4_reco = track_mom4_true;

        TVector3 track_vtx_true = part->GetVertex();
        TVector3 track_vtx_reco = track_vtx_true;

        if (SMEAR_OPTION==0);
        if (SMEAR_OPTION==1)
        {
          track_mom4_reco = smearMomDM(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosDM(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==2)
        {
          track_mom4_reco = smearMomLBL(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosLBL(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==3)
        {
          track_mom4_reco = smearMomHybrid(track_mom4_true, BFIELD_TYPE);
          track_vtx_reco = smearPosHybrid(track_mom4_true.Vect(), track_vtx_true);
        }
        if (SMEAR_OPTION==4)
        {
          track_mom4_reco = smearMomATHENA(track_mom4_true);
          track_vtx_reco = smearPosATHENA(track_mom4_true.Vect(), track_vtx_true);
        }
        if (track_mom4_reco.E()>1000 || track_vtx_reco.Mag()>1000) continue; // outside eta or momentum range

        // single track DCA cut
        double track_dca = dcaSigned(track_mom4_reco.Vect(),track_vtx_reco,evt_vtx_reco);
        if (TRK_DCA>-99 && fabs(track_dca)<TRK_DCA) continue;

        //==========================
        //      PID selection
        //==========================
        bitset<4> track_binary_id;
        if (ID_OPTION==-1);
        else if (ID_OPTION==0) passing_DM_PID(track_mom4_reco,track_binary_id);
        else if (ID_OPTION==1 || ID_OPTION==2)
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        else
        {
          passing_DIRC(track_mom4_reco,track_binary_id,BFIELD_TYPE,0.5);
          passing_hside_dRICH(track_mom4_reco,track_binary_id);
          passing_eside_dRICH(track_mom4_reco,track_binary_id);
        }
        if (verbosity>2) cout << "track_binary_id " <<track_binary_id.to_ulong() << endl;


        //==============================================================================
        //    Assumption: Cherekov detectors (mass ordering, no mu/pi separation)
        // assign probability according to the multiplicity (NB: qualitative for now)
        //==============================================================================
        float prob_pi = -999;
        float prob_K = -999;
        float prob_e = -999;
        float prob_p = -999;
        identify_charged_hadrons(part->Id(), track_binary_id, prob_e, prob_pi, prob_K, prob_p);

        int quark_index = find_quark_origin(py_evt,part);
        erhic::ParticleMC* quark_orig = py_evt->GetTrack(quark_index-1);

        //============================
        //  opposite charge selection
        //============================
        bool flag_pos = false, flag_neg = false;
        // no charge info in the input tree, use id to separate charge
        if (abs(part->Id())==11 && part->Id()<0) flag_pos = true;
        else if (abs(part->Id())==11 && part->Id()>0) flag_neg = true;
        else if (part->Id()>0) flag_pos = true;
        else flag_neg = true;

        if (flag_neg)
        {
          negl_p_true.push_back(track_mom4_true);
          negl_vtx_true.push_back(track_vtx_true);
          negl_id_true.push_back(part->Id());
          negl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL) negl_parent_id_true.push_back(parent->Id());
          else negl_parent_id_true.push_back(-9999);

          negl_p_reco.push_back(track_mom4_reco);
          negl_vtx_reco.push_back(track_vtx_reco);
          negl_prob_e.push_back(prob_e);
          negl_prob_pi.push_back(prob_pi);
          negl_prob_K.push_back(prob_K);
          negl_prob_p.push_back(prob_p);

          if (quark_orig!=NULL) negl_quark_p.push_back(quark_orig->Get4Vector());
          else negl_quark_p.push_back(TLorentzVector(0,0,0,0));
        }

        if (flag_pos)
        {
          posl_p_true.push_back(track_mom4_true);
          posl_vtx_true.push_back(track_vtx_true);
          posl_id_true.push_back(part->Id());
          posl_parent_index_true.push_back(part->GetParentIndex());

          erhic::ParticleMC* parent = py_evt->GetTrack(part->GetParentIndex()-1);
          if (parent!=NULL) posl_parent_id_true.push_back(parent->Id());
          else posl_parent_id_true.push_back(-9999);

          posl_p_reco.push_back(track_mom4_reco);
          posl_vtx_reco.push_back(track_vtx_reco);
          posl_prob_e.push_back(prob_e);
          posl_prob_pi.push_back(prob_pi);
          posl_prob_K.push_back(prob_K);
          posl_prob_p.push_back(prob_p);

          if (quark_orig!=NULL) posl_quark_p.push_back(quark_orig->Get4Vector());
          else posl_quark_p.push_back(TLorentzVector(0,0,0,0));
        }
      }
    }

    void fill_Kpi_mass(const int charge_type, const int ipart1, const int ipart2)
    {
      if (charge_type>1) return;
      if (ipart1<0 || ipart2<0) return;

      bool is_SG = true;
      TLorentzVector kaon_p, pion_p;
      TLorentzVector quark_p;
      if (charge_type==0)
      { // K-pi+
        if (fabs(negl_parent_id_true[ipart1])!=421) is_SG = false;
        if (fabs(posl_parent_id_true[ipart2])!=421) is_SG = false;
        if (negl_parent_index_true[ipart1]!=posl_parent_index_true[ipart2]) is_SG = false;

        kaon_p = negl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(negl_id_true[ipart1])!=321) is_SG = false;

        pion_p = posl_p_reco[ipart2];
        pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
        if (fabs(posl_id_true[ipart2])!=211) is_SG = false;

        quark_p = negl_quark_p[ipart1];
      }
      else
      { // K+pi-
        if (fabs(posl_parent_id_true[ipart1])!=421) is_SG = false;
        if (fabs(negl_parent_id_true[ipart2])!=421) is_SG = false;
        if (posl_parent_index_true[ipart1]!=negl_parent_index_true[ipart2]) is_SG = false;

        kaon_p = posl_p_reco[ipart1];
        kaon_p.SetXYZM(kaon_p.X(),kaon_p.Y(),kaon_p.Z(),KMASS);
        if (fabs(posl_id_true[ipart1])!=321) is_SG = false;

        pion_p = negl_p_reco[ipart2];
        pion_p.SetXYZM(pion_p.X(),pion_p.Y(),pion_p.Z(),PIMASS);
        if (fabs(negl_id_true[ipart2])!=211) is_SG = false;

        quark_p = posl_quark_p[ipart1];
      }

      TLorentzVector pair(kaon_p+pion_p);
      int ietabin = -9999;
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (pair.PseudoRapidity()>=eta_lo[ieta] && pair.PseudoRapidity()<eta_hi[ieta]) ietabin = ieta;
      }
      if (ietabin<0) return;

      double frag_z = hadron_beam.Dot(pair)/(nu_true*hadron_beam.M());  // z= Pp/Pq where Pq=nuM

      cout<<"D0 pt "<<pair.Pt()<<" eta "<<pair.PseudoRapidity()<<endl;

      fg2d_Kpimass_vs_p[charge_type][ietabin]->Fill(pair.M(),pair.Pt());
      fg2d_Kpimass_vs_p[2][ietabin]->Fill(pair.M(),pair.Pt());
      fg2d_Kpimass_vs_z[charge_type][ietabin]->Fill(pair.M(),frag_z);
      fg2d_Kpimass_vs_z[2][ietabin]->Fill(pair.M(),frag_z);
      if (!is_SG)
      {
        bg2d_Kpimass_vs_p[charge_type][ietabin]->Fill(pair.M(),pair.Pt());
        bg2d_Kpimass_vs_p[2][ietabin]->Fill(pair.M(),pair.Pt());
        bg2d_Kpimass_vs_z[charge_type][ietabin]->Fill(pair.M(),frag_z);
        bg2d_Kpimass_vs_z[2][ietabin]->Fill(pair.M(),frag_z);
      }
      else
      {
        // cout<<"D0 pt = "<<pair.Pt()<<" eta = "<<pair.PseudoRapidity()<<" z = "<<frag_z<<endl;
        // if (struck_quark.E()>0.001)
        // {
        //   cout<<"D0 energy = "<<pair.E()<<" struck quark energy = "<<struck_quark.E()<<" frac = "<<pair.E()/struck_quark.E()<<endl;
        //   // cout<<"D0 Pt = "<<pair.Pt()<<" struck quark Pt = "<<struck_quark.Pt()<<" frac = "<<
        //   cout<<"frac = "<<(pair.Vect()).Dot(struck_quark.Vect())/(struck_quark.Vect()).Dot(struck_quark.Vect())<<endl;
        //   cout<<"**********************"<<endl;
        // }

        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          if (pair.Pt()>=ppt_lo[ipt] && pair.Pt()<ppt_hi[ipt])
          {
            h2d_K_D0_p_vs_eta[ietabin][ipt]->Fill(kaon_p.P(),kaon_p.PseudoRapidity());
            h2d_pi_D0_p_vs_eta[ietabin][ipt]->Fill(pion_p.P(),pion_p.PseudoRapidity());
          }
        }

        int iQ2bin = -9999;
        for (int iQ2 = 0; iQ2 < Q2bin-1; ++iQ2)
        {
          if (Q2_true>=Q2_lo[iQ2] && Q2_true<Q2_hi[iQ2]) iQ2bin = iQ2;
        }
        int ixbin = -9999;
        for (int ix = 0; ix < xbin-1; ++ix)
        {
          if (x_true>=x_lo[ix] && x_true<x_hi[ix]) ixbin = ix;
        }

        if (iQ2bin>=0 && ixbin>=0)
        {
          h2d_D0_pt_vs_eta[iQ2bin][ixbin]->Fill(pair.Pt(),pair.PseudoRapidity());
          h2d_D0_z_vs_eta[iQ2bin][ixbin]->Fill(frag_z,pair.PseudoRapidity());

          h2d_D0_pt_vs_eta[Q2bin-1][xbin-1]->Fill(pair.Pt(),pair.PseudoRapidity());
          h2d_D0_z_vs_eta[Q2bin-1][xbin-1]->Fill(frag_z,pair.PseudoRapidity());
        }

        h2d_ztheo_vs_zjet->Fill((pair.Vect()).Dot(quark_p.Vect())/(quark_p.Vect()).Dot(quark_p.Vect()),frag_z);
      }
    }

    void FillD0Pairs()
    {
      for (int ineg = 0; ineg < negl_id_true.size(); ++ineg)
      {
        for (int ipos = 0; ipos < posl_id_true.size(); ++ipos)
        {
          // ==========================
          //    decay topology cut
          // ==========================
          // pair DCA < cut value
          TVector3 dca_pair = negl_vtx_reco[ineg]-posl_vtx_reco[ipos];
          dca_pair.SetZ(0);
          // cout << "dca.pair " << dca_pair.Mag() << endl;
          if (PAIR_DCA>-99 && dca_pair.Mag()>PAIR_DCA) continue;
          // cout << "double check dca.pair " << dca_pair.Mag() << endl;

          // Decay length > cut value
          TVector3 decay_vtx_reco = (negl_vtx_reco[ineg]+posl_vtx_reco[ipos])*0.5;
          TVector3 decay_l = decay_vtx_reco-evt_vtx_reco;
          decay_l.SetZ(0);
          if (DECAY_L>-99 && decay_l.Mag()<DECAY_L) continue;

          // D0 DCA > cut value
          TVector3 D0_vec = negl_p_reco[ineg].Vect()+posl_p_reco[ipos].Vect();
          double D0_dca = dcaXY(D0_vec, decay_vtx_reco, evt_vtx_reco);
          if (D0_DCA>-99 && fabs(D0_dca)<D0_DCA) continue;

          // D0 cos(theta) > cut value
          TVector3 D0_vecT = D0_vec;
          D0_vecT.SetZ(0);
          double D0_costheta = TMath::Cos(D0_vecT.Angle(decay_l));
          if (D0_COSTHETA>-99 && D0_costheta<D0_COSTHETA) continue;

          if (ID_OPTION==-1)
          { // no hID (but with eID)
            if (fabs(negl_id_true[ineg])!=11 && fabs(posl_id_true[ipos])!=11)
            {
              fill_Kpi_mass(0,ineg,ipos);
              fill_Kpi_mass(1,ipos,ineg);
            }
          }

          if (ID_OPTION==0)
          { // PID with no low momentum cutoff
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }

          if (ID_OPTION==1)
          { // PID with low momentum cutoff & some mis-identified pi, K
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]>0.5)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]>0.5)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }

          if (ID_OPTION==2)
          { // PID with low momentum cutoff & all identified pi, K
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }

          if (ID_OPTION==3)
          { // PID with low momentum cutoff & all identified pi, K (DIRC at 50cm)
            if (posl_prob_pi[ipos]==1 && negl_prob_K[ineg]==1)
            {
              fill_Kpi_mass(0,ineg,ipos);
            }
            if (negl_prob_pi[ineg]==1 && posl_prob_K[ipos]==1)
            {
              fill_Kpi_mass(1,ipos,ineg);
            }
          }
        }
      }
    }

    void Write()
    {
      for (int icharge = 0; icharge < chargebin; ++icharge)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          fg2d_Kpimass_vs_p[icharge][ieta]->Write();
          bg2d_Kpimass_vs_p[icharge][ieta]->Write();
          fg2d_Kpimass_vs_z[icharge][ieta]->Write();
          bg2d_Kpimass_vs_z[icharge][ieta]->Write();
        }
      }

      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        for (int ipt = 0; ipt < pptbin; ++ipt)
        {
          h2d_K_D0_p_vs_eta[ieta][ipt]->Write();
          h2d_pi_D0_p_vs_eta[ieta][ipt]->Write();
        }
      }

      for (int iQ2 = 0; iQ2 < Q2bin; ++iQ2)
      {
        for (int ix = 0; ix < xbin; ++ix)
        {
          h2d_D0_pt_vs_eta[iQ2][ix]->Write();
          h2d_D0_z_vs_eta[iQ2][ix]->Write();

          h2d_D0_pt_vs_eta_gen[iQ2][ix]->Write();
          h2d_D0_z_vs_eta_gen[iQ2][ix]->Write();
        }
      }

      h2d_ztheo_vs_zjet->Write();
    }
};



void D0_replacement(const char* inFile = "ep_allQ2.20x100.small.root", const char* outFile = "hist.root", int nevt = 0, const int smear_option = 0, const int Bfield_type = 0, const int PID_option = 0)
{ // smear_2nd_vtx & momentum: 0--no smearing, 1--DM smearing, 2--LBL smearing, 3--Hybrid smearing, 4--ATHENA smearing
  // Bfield_type: 0--Barbar, 1--Beast
  // 0--no hID (but with eID), 1--PID with no low momentum cutoff, 2--PID with low momentum cutoff & some mis-identified pi, K, 3--PID with low momentum cutoff & all identified pi, K
  // DCA_cut: 0--no cut, 1--cut on DCA
  // do_correct_vertex: 0--no correction applied, 1--D0 and Lc vertex decay length recalculated, modifications applied to input root files and data used in above analysis

  // PDG data table
  pdg = new TDatabasePDG();

  //Event Class
  erhic::EventPythia *event(NULL); //Note that I use Pointer

  //Particle Class
  erhic::ParticleMC *particle(NULL); //Also use Pointer

  // erhic::ParticleMC *child(NULL);

  //Load ROOT File
  TFile *f = new TFile(inFile);

  //Get EICTree Tree
  TTree *tree = (TTree*)f->Get("EICTree");

  Int_t nEntries = tree->GetEntries();
  cout<<"-------------------------------"<<endl;
  cout<<"Total Number of Events = "<<nEntries<<endl<<endl;

  //Access event Branch
  tree->SetBranchAddress("event",&event); //Note &event, even with event being a pointer

  // ATHENA smeaing
  TFile* f_athena_tracks = new TFile("ATHENA_Resolutions_r.root","READ");
  TFile* f_athena_vertex = new TFile("VertexRes_ATHENA.root","READ");
  if (smear_option==4)
  {
    cout << "Setup ATHENA smearing parameters" << endl;
    setup_ATHENA_single_track_smearing(f_athena_tracks);
    setup_ATHENA_PV_smearing(f_athena_vertex);
  }

  //Define Some Variables
  Float_t Q2(0);
  Int_t nParticles(0);

  D0_reco ana_D0;
  ana_D0.SetLowP(0.1);
  ana_D0.SetDCACuts();
  ana_D0.SetIDCuts(PID_option);
  ana_D0.SetSmearType(smear_option);
  ana_D0.SetBFieldType(Bfield_type);
  // ana_D0.SetDoCorrectVertex(do_correct_vertex);

  Lc_reco ana_Lc;
  ana_Lc.SetDCACuts();
  ana_Lc.SetIDCuts(PID_option);
  ana_Lc.SetSmearType(smear_option);
  ana_Lc.SetBFieldType(Bfield_type);
  // ana_Lc.SetDoCorrectVertex(do_correct_vertex);

  //Loop Over Events
  if (nevt == 0) nevt = nEntries;
  for(Int_t ievt = 0; ievt < nevt; ievt++)
  {
    if (ievt%10000==0) cout<<"Processing event = "<<ievt<<"/"<<nevt<<endl;

    tree->GetEntry(ievt);

    //Write Out Q2
    Q2 = (Float_t) event->GetQ2(); //Can also do event->QSquared
    // printf("For Event %d, Q^2 = %.3f GeV^2!\n",ievt,Q2);
    // if (Q2>10) continue; // process low Q2

    //Get Total Number of Particles
    nParticles = event->GetNTracks();
    // printf("For Event %d, we have %d particles!\n",ievt,nParticles);

    TVector3 evt_vtx_true(0,0,0); // not sure how to get event vertex
    TVector3 evt_vtx_reco = evt_vtx_true;
    if (smear_option==4)
    {
      float multi_charged_hadron_35 = 0;
      float multi_charged_hadron_10 = 0;
      for (int ipart = 0; ipart < nParticles; ++ipart)
      {
        particle = event->GetTrack(ipart);
        if (particle->GetStatus()!=1) continue; // only loop through final stable particles

        if (abs(particle->Id())!=2212 && abs(particle->Id())!=211 && abs(particle->Id())!=321) continue; // proton, pion, kaon
        if(fabs(particle->GetEta())<3.5 && particle->GetP()>0.5) multi_charged_hadron_35++;
        if(fabs(particle->GetEta())<1.0 && particle->GetP()>0.5) multi_charged_hadron_10++;
      }
      evt_vtx_reco = smearPVTATHENA(evt_vtx_true, multi_charged_hadron_35); // smear transverse direction
      // comment out the longgitunal smearing for now
      // evt_vtx_reco = smearPVZATHENA(evt_vtx_true, multi_charged_hadron_10); // smear longitudinal direction
    }
    if (evt_vtx_reco.X()<-999 || evt_vtx_reco.Y()<-999)
    {
      // cout << evt_vtx_reco.X() << ", " << evt_vtx_reco.Y() << ", " << evt_vtx_reco.Z() << " mm" << endl;
      continue; // No valid PV smearing in transverse direction, jump to next event
    }

    ana_D0.SetVectTrue(evt_vtx_true);
    ana_D0.SetVectReco(evt_vtx_reco); // NB: no smear on primary vertex yet

    ana_D0.SetQ2True(event->GetQ2());
    ana_D0.SetXTrue(event->GetX());

    ana_D0.SetNuTrue(event->GetNu());

    ana_D0.FillGenKin(event);

    ana_D0.FillSingleTracks(event);
    ana_D0.FillD0Pairs();

    ana_Lc.SetVectTrue(evt_vtx_true);
    ana_Lc.SetVectReco(evt_vtx_reco); // NB: no smear on primary vertex yet

    ana_Lc.SetQ2True(event->GetQ2());
    ana_Lc.SetXTrue(event->GetX());

    ana_Lc.SetNuTrue(event->GetNu());

    ana_Lc.FillGenKin(event);

    ana_Lc.FillSingleTracks(event);
    ana_Lc.FillLcTriplets();
  }

  TFile* fout = new TFile(outFile,"recreate");

  ana_D0.Write();
  ana_Lc.Write();

  fout->Write();



  
}