
from __future__ import print_function

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjext
import ecorrel

import os
import argparse
import ROOT
import pyhepmc # this is the new
# import pyhepmc_ng # this is the old

from pyjetty.mputils import *
from pyjetty.mputils.mputils import pinfo, pwarning

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import numpy as np
import array as arr


import hepmc2antuple_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)
# Automatically set Sumw2 when creating new histograms
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()


class HepMC_quark_gluon(hepmc2antuple_base.HepMC2antupleBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(HepMC_quark_gluon, self).__init__(**kwargs)
    self.init()
    print(self)


    # PDG ID values for quarks and gluons
    self.quark_pdg_ids = [1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, -4, -5, -6, -7, -8]
    self.down_pdg_ids = [1, -1]
    self.up_pdg_ids = [2, -2]
    self.strange_pdg_ids = [3, -3]
    self.charm_pdg_ids = [4, -4]
    self.beauty_pdg_ids = [5, -5]
    self.gluon_pdg_ids = [9, 21] 

    # hadron level - ALICE tracking restriction
    self.max_eta_hadron = 0.9

    # analysis variables
    self.use_ptRL = False #change to parameters later
    
    
  #---------------------------------------------------------------
  def main(self):

    self.hepmc = 2

    # Open the HepMC file
    filename = "/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/herwig_onthefly/LHC_13000_MPI_10events.hepmc"

    # if self.hepmc == 3:
    #   input_hepmc = pyhepmc_ng.ReaderAscii(self.input)
    if self.hepmc == 2:
      # input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(self.input)
      input_hepmc = pyhepmc.ReaderAsciiHepMC2(self.input)
    #   with pyhepmc.open("modified_event.hepmc2", "w", format="hepmc2", precision=3) as f:
    # f.write(event)

    if input_hepmc.failed():
      print ("[error] unable to read from {}".format(self.input))
      sys.exit(1)

    # initialize histograms and jet tools
    self.initialize_hists()
    self.init_jet_tools()

    # start to analyze file
    # event_hepmc = pyhepmc_ng.GenEvent()
    event_hepmc = pyhepmc.GenEvent()
    print("checkpoint 1")

    while not input_hepmc.failed():
      ev = input_hepmc.read_event(event_hepmc)
      # print("input_hepmc", self.ev_id, input_hepmc.failed())
      print("checkpoint 2")
      if input_hepmc.failed():
        ("failed! error!")
        break
      print("checkpoint 3")
      self.fill_event(event_hepmc)
      # if self.ev_id >= 0:
      #   break
      self.increment_event()
      if self.nev > 0 and self.ev_id > self.nev:
        break
      
    # self.finish()

  #---------------------------------------------------------------
  # Initiate histograms
  #---------------------------------------------------------------
  def initialize_hists(self):
    observable = "EEC" #change to access config list later
    partontypeslist = ["charm", "light", "gluon", "inclusive"]

    # Store a list of all the histograms just so that we can rescale them later
    jetR = 0.4 # fix this block to config list later!!
    hist_list_name = "hist_list_R%s" % str(jetR).replace('.', '')
    setattr(self, hist_list_name, [])

    if (self.use_ptRL):
      self.obs_bins_EEC = np.logspace(np.log10(1E-4), np.log10(100), 51)
    else:
      self.obs_bins_EEC = np.logspace(np.log10(1E-4), np.log10(1), 51)

    obs_bins = getattr(self, "obs_bins_" + observable)
    # Use more finely binned pT bins for TH2s than for the RMs
    pt_bins = arr.array('d', list(range(0, 201, 1)))
    rapi_bins = np.linspace(-5,5,201)
    z_bins = np.linspace(0, 1.01, 102)


    dim = 5
    nbins  = [len(pt_bins)-1, len(pt_bins)-1, len(rapi_bins)-1, len(z_bins)-1, 50]
    min_li = [pt_bins[0],     pt_bins[0],      rapi_bins[0],      obs_bins[0],      z_bins[0]]
    max_li = [pt_bins[-1],    pt_bins[-1],     rapi_bins[-1],     obs_bins[-1],     z_bins[-1]]

    nbins = (nbins)
    xmin = (min_li)
    xmax = (max_li)
    
    nbins_array = arr.array('i', nbins)
    xmin_array = arr.array('d', xmin)
    xmax_array = arr.array('d', xmax)


    for parton_type in partontypeslist:
      if (self.use_ptRL):
        title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'y', 'z', '#it{p}_{T}#it{R}_{L}' ]
      else:
        title = [ '#it{p}_{T}^{ch jet}', '#it{p}_{T}^{D^{0}}', 'y', 'z', '#it{R}_{L}' ]

      # make THnSparse for parton EECs
      name = ('h_%s_JetPt_%s_R%s' % (observable, parton_type, jetR))
      hsparse = ROOT.THnSparseD(name,"%s-init_hsparsejet; #it{p}_{T,%s}^{ch jet}; #it{p}_{T}^{D^{0}}; y;R_{L}^{%s}" %(parton_type[0], parton_type[0] + "-init", parton_type[0] + "-init"), dim,  nbins_array, xmin_array, xmax_array)
      hsparse.Sumw2()

      # make another of THnSparse for the jet level (above is pair level)
      name_jetpt = ('h_JetPt_%s_R%s_jetlevel' % (parton_type, jetR))
      hsparse_jetpt = ROOT.THnSparseD(name_jetpt,"%s-init_hsparsejet_jetlevel; #it{p}_{T,%s}^{ch jet}; #it{p}_{T}^{D^{0}}; y" %(parton_type[0], parton_type[0] + "-init"), dim-1,  nbins_array[:-1], xmin_array[:-1], xmax_array[:-1])
      hsparse_jetpt.Sumw2()

      # assign THnSparse axes titles
      for i in range(0,dim):
        hsparse.GetAxis(i).SetTitle(title[i])
        if i<dim-1:
          hsparse_jetpt.GetAxis(i).SetTitle(title[i])
        if i == 0 or i == 1:
          hsparse.SetBinEdges(i, pt_bins)
          hsparse_jetpt.SetBinEdges(i, pt_bins)
        if i == 2:
          hsparse.SetBinEdges(i, rapi_bins)
          hsparse_jetpt.SetBinEdges(i, rapi_bins)
        if i == 3:
          hsparse.SetBinEdges(i, z_bins)
          hsparse_jetpt.SetBinEdges(i, z_bins)
        if i == 4:
          hsparse.SetBinEdges(i, obs_bins)
      
      setattr(self, name, hsparse)
      getattr(self, hist_list_name).append(hsparse)
      setattr(self, name_jetpt, hsparse_jetpt)
      getattr(self, hist_list_name).append(hsparse_jetpt)


  #---------------------------------------------------------------
  # Initiate jet defs, selectors, and sd (if required)
  #---------------------------------------------------------------
  def init_jet_tools(self):

    self.jetR_list = [0.4] # fix this to config file later!!
    for jetR in self.jetR_list:
      jetR_str = str(jetR).replace('.', '')

      # set up our jet definition and a jet selector
      jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
      setattr(self, "jet_def_R%s" % jetR_str, jet_def)

    pwarning('max eta for particles after hadronization set to', self.max_eta_hadron)
    parts_selector_h = fj.SelectorAbsEtaMax(self.max_eta_hadron)
    track_selector_ch = fj.SelectorPtMin(0.15) & parts_selector_h #ALICE parameters
    setattr(self, "track_selector_ch", track_selector_ch)

    for jetR in self.jetR_list:
      jetR_str = str(jetR).replace('.', '')

      jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
      #jet_selector = fj.SelectorPtMin(0.) & fj.SelectorAbsEtaMax(self.max_eta_hadron - jetR)
      setattr(self, "jet_selector_R%s" % jetR_str, jet_selector)

      count1 = 0  # Number of partonic parents which match to >1 ch-jets
      setattr(self, "count1_R%s" % jetR_str, count1)
      count2 = 0  # Number of partonic parents which match to zero ch-jets
      setattr(self, "count2_R%s" % jetR_str, count2)




  #---------------------------------------------------------------
  # Fill histograms for an event
  #---------------------------------------------------------------
  def fill_event(self, event_hepmc):

    # Loop through each event
    print("Event", event_hepmc.event_number, "has", len(event_hepmc.particles), "particles.")
    [print(ind, x.id, x.status, x.pid, x) for ind, x in enumerate(event_hepmc.particles) if ind<10]

    print("TESTING", event_hepmc.pdf_info, event_hepmc.pdf_info.parton_id1, event_hepmc.pdf_info.parton_id2)
    print("VERTICES", len(event_hepmc.vertices))
    # print("BEGIN", event_hepmc.vertices_begin)
    # print("END", event_hepmc.vertices_end)
    for ind, vertex in enumerate(event_hepmc.vertices):
    # for vertex in event_hepmc.vertices_begin
	  #     v != evt->vertices_end(); ++v):
      # Get the position of the vertex (assuming it's in space-time coordinates)
      # position = vertex.position()
      # Print the position (x, y, z, t)
      print("VERTEX", vertex.attributes, "//")
      print(vertex.id, vertex.in_event, vertex.status, vertex) #, vertex.parent_event) #, f"Vertex position: {position}")
      # Loop over the particles produced at the vertex (children of the vertex)
      print("PARTICLES IN", len(vertex.particles_in)) #, vertex.particles_in)
      for part in vertex.particles_in:
        print("  ", part.id, part.pid, part.status, part.momentum)
        # print(f"  Produced particle PDG ID: {particle.pid}, Momentum: {particle.momentum()}")
      # Loop over the particles that caused the vertex (parents of the vertex)
      print("PARTIFLES OUT", len(vertex.particles_out)) #, vertex.particles_out)
      for part in vertex.particles_out:
        print("  ", part.id, part.pid, part.status, part.momentum)
        # print(f"  Parent particle PDG ID: {particle.pid}")
      # if ind == 100:
      #   break

    vert_status_not0 = 0
    for ind, vertex in enumerate(event_hepmc.vertices):
      if vertex.status != 0:
        vert_status_not0+=1
    print("VERTSTATUSNOT0 = ", vert_status_not0)

    # in each event, we want to first identify the outgoing partons 
    parent_parton_1 = event_hepmc.particles[5]
    parent_parton_2 = event_hepmc.particles[6]
    self.parents = [parent_parton_1, parent_parton_2]
    self.parent_ids = [parent_parton_1.pid, parent_parton_2.pid]

    # then we want to find the jets
    parts_herwig_hch = pythiafjext.vectorize_select(event_hepmc, [pythiafjext.kFinal, pythiafjext.kVisible, pythiafjext.kCharged], 0, True)

    # then we want to see if any of the jets can be matched to the parent partons
    # then we take all those jets and add to the EEC histogram



    # # Loop through particles
    # for particle in event_hepmc.particles:
    #   pid = particle.pid  # PDG ID
    #   # print("PID", pid)

      
      
    

    

if __name__ == '__main__':
  
  parser = argparse.ArgumentParser(description='hepmc to ALICE Ntuple format', prog=os.path.basename(__file__))
  parser.add_argument('-i', '--input', help='input file', default='', type=str, required=True)
  parser.add_argument('-o', '--output', help='output root file', default='', type=str, required=True)
  parser.add_argument('--as-data', help='write as data - tree naming convention', action='store_true', default=False)
  parser.add_argument('--hepmc', help='what format 2 or 3', default=2, type=int)
  parser.add_argument('--nev', help='number of events', default=-1, type=int)
  parser.add_argument('-g', '--gen', help='generator type: pythia, herwig, jewel, jetscape, martini, hybrid', default='pythia', type=str, required=True)
  parser.add_argument('--no-progress-bar', help='whether to print progress bar', action='store_true', default=False)
  parser.add_argument('-p', '--include-parton', help='include additional tree of final-state partons', action='store_true', default=False)
  parser.add_argument('-d', '--include-D0', help='include additional tree of D0 information and mother IDs', action='store_true', default=False)
  args = parser.parse_args()
  
  converter = HepMC_quark_gluon(input = args.input, output = args.output, as_data = args.as_data, hepmc = args.hepmc, nev = args.nev, gen = args.gen, no_progress_bar = args.no_progress_bar, include_parton = args.include_parton, include_D0 = args.include_D0)
  converter.main()