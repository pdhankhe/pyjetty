#!/usr/bin/env python

from __future__ import print_function

import os
import argparse

# import pyhepmc
import pyhepmc_ng

import hepmc2antuple_base

# jit improves execution time by 18% - tested with jetty pythia8 events
# BUT produces invalid root file
# from numba import jit
# @jit

################################################################
class HepMC2antuple(hepmc2antuple_base.HepMC2antupleBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(HepMC2antuple, self).__init__(**kwargs)
    self.init()
    print(self)
    
  #---------------------------------------------------------------
  def main(self):
  
    if self.hepmc == 3:
      # input_hepmc = pyhepmc.io.ReaderAscii(self.input)
      input_hepmc = pyhepmc_ng.ReaderAscii(self.input)
    if self.hepmc == 2:
      # input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(self.input)
      input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(self.input)

    if input_hepmc.failed():
      print ("[error] unable to read from {}".format(self.input))
      sys.exit(1)

    # event_hepmc = pyhepmc.GenEvent()
    event_hepmc = pyhepmc_ng.GenEvent()

    while not input_hepmc.failed():
      ev = input_hepmc.read_event(event_hepmc)
      # print("input_hepmc", self.ev_id, input_hepmc.failed())
      if input_hepmc.failed():
        break
      self.fill_event(event_hepmc)
      self.increment_event()
      if self.nev > 0 and self.ev_id > self.nev:
        break
      
    self.finish()
    # print("finish event")

  #---------------------------------------------------------------
  def fill_event(self, event_hepmc):
    # print("in fill event!!")

    self.t_e.Fill(self.run_number, self.ev_id, 0, 0)

    for part in event_hepmc.particles:

      # get mother here --> need PID
      # if not looking at D0s or there is != 1 mother, then set motherPID=0
      motherPID = 0
      if self.include_D0:
        # print("PART", part)
        # print("PARENTS", part.parents)
        mothers = part.parents
        if len(mothers) == 1:
          mother = mothers[0]
          motherPID = mothers[0].pid
          if motherPID == 421 or motherPID == -421:
            # print("motherPID", motherPID) #, self.pdg.GetParticle(motherPID).GetName())
            
            if self.isD0toKpidecay(mother) and abs(part.pid)==321: # TODO: check that D0 goes to Kpi, and only add the D0 if the daughter is a kaon
              mothers_of_D0 = mother.parents
              # print("D0 MOTHERS", [mod.pid for mod in mothers_of_D0])
              if (len(mothers_of_D0) == 1):
                mother_of_D0 = mothers_of_D0[0]
                mother_of_D0_PID = mother_of_D0.pid
              else:
                mother_of_D0_PID = 0
              self.t_D.Fill(self.run_number, self.ev_id, mother.momentum.pt(), mother.momentum.eta(), mother.momentum.phi(), mother.momentum.rap(), mother.pid, mother_of_D0_PID)
              # print("rapidity", part.momentum.rap())
            elif self.isD0toKpidecay(mother) == False: # closing if D0->kpi
              motherPID = 0
          else: # closing if motherpid = 421
            motherPID = 0
            
    
      '''
      #if status == 1, then end_vertex is None
      if not (part.status == 1 and part.end_vertex == None):
        print("checking; PART:", part, "STATUS:", part.status, "END VTX:", part.end_vertex, "PID:", part.pid, "PDG:") #, self.pdg, "GEN:", self.gen)
      '''

      if self.accept_particle(part, part.status, part.end_vertex, part.pid, self.pdg, self.gen):
        # print("in here!")
      
        self.particles_accepted.add(self.pdg.GetParticle(part.pid).GetName())          
        self.t_p.Fill(self.run_number, self.ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid, motherPID)          
      
      elif self.include_parton and self.accept_particle(part, part.status, part.end_vertex, part.pid, self.pdg, self.gen, parton=True):

        self.partons_accepted.add(self.pdg.GetParticle(part.pid).GetName())
        self.t_pp.Fill(self.run_number, self.ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid)
      
        
#---------------------------------------------------------------
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
  
  converter = HepMC2antuple(input = args.input, output = args.output, as_data = args.as_data, hepmc = args.hepmc, nev = args.nev, gen = args.gen, no_progress_bar = args.no_progress_bar, include_parton = args.include_parton, include_D0 = args.include_D0)
  converter.main()
