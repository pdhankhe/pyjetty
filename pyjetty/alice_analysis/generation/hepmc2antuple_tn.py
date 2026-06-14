#!/usr/bin/env python

from __future__ import print_function

import os
import argparse
import re

import pyhepmc # use this for perlmutter -- also need to change in /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/select_particles.py
# import pyhepmc_ng # use this for hiccup

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
      input_hepmc = pyhepmc.io.ReaderAscii(self.input)
      # input_hepmc = pyhepmc_ng.ReaderAscii(self.input)
    if self.hepmc == 2:
      input_hepmc = pyhepmc.io.ReaderAsciiHepMC2(self.input)
      # input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(self.input)

    if input_hepmc.failed():
      print ("[error] unable to read from {}".format(self.input))
      sys.exit(1)

    event_hepmc = pyhepmc.GenEvent()
    # event_hepmc = pyhepmc_ng.GenEvent()

    while not input_hepmc.failed():
      ev = input_hepmc.read_event(event_hepmc)
      # print("input_hepmc", self.ev_id, input_hepmc.failed())
      if input_hepmc.failed():
        break
      self.fill_event(event_hepmc)
      self.increment_event()
      if self.nev > 0 and self.ev_id > self.nev:
        break
    
    if self.include_herwig_parton:
        self.fill_herwig_parton_info()
      
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
        # if (self.ev_id > 750 and self.ev_id < 800): # it was event 783 
        #   print(self.ev_id, "particle pid", part.pid)
        #   print("particle pid name", self.pdg.GetParticle(part.pid))
        if (part.pid == -14122):
          self.particles_accepted.add("Lambda_c+")
        else:
          self.particles_accepted.add(self.pdg.GetParticle(part.pid).GetName())  

        if self.include_D0:
          self.t_p.Fill(self.run_number, self.ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid, motherPID)
        else:
          self.t_p.Fill(self.run_number, self.ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid)        

        if self.for_jse:
          self.t_j.Fill(self.run_number, self.ev_id, part.momentum.px, part.momentum.py, part.momentum.pz, part.momentum.e, part.pid)  
      
      elif self.include_parton and self.accept_particle(part, part.status, part.end_vertex, part.pid, self.pdg, self.gen, parton=True):

        self.partons_accepted.add(self.pdg.GetParticle(part.pid).GetName())
        self.t_pp.Fill(self.run_number, self.ev_id, part.momentum.pt(), part.momentum.eta(), part.momentum.phi(), part.pid)
      


  def parse_log(self, path):
    """
    Parse a Herwig event log. For each event, extract ONLY the primary
    sub-process incoming/outgoing hard partons.
    Returns: event_number -> {'incoming': [...], 'outgoing': [...]}
    where each entry is (pid, px, py, pz, e).
    """
    events = {}
    with open(path) as fh:
        lines = fh.readlines()

    evt_re = re.compile(r"Event number\s+(\d+)")
    hdr_re = re.compile(r"^\s*\d+\s+\S+\s+(-?\d+)\b")
    mom_re = re.compile(
        r"^\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    )

    cur_evt = None
    in_primary = False     # are we inside the Primary sub-process block?
    section = None         # 'incoming' | 'outgoing' | None (we skip intermediates)
    pending_pid = None

    for line in lines:
      s = line.strip()

      m = evt_re.search(line)
      if m:
        cur_evt = int(m.group(1))
        events[cur_evt] = {'incoming': [], 'intermediates': [], 'outgoing': []}
        in_primary = False
        section = None
        pending_pid = None
        continue

      if cur_evt is None:
        continue

      # Enter the primary sub-process block
      if "Primary sub-process performed by" in line:
        in_primary = True
        section = None
        pending_pid = None
        continue

      # ANY of these ends the primary block. After the primary outgoing
      # section, the next thing is a "------" divider, then "Step 1", then
      # "Secondary sub-process". Stop at the first of them.
      if in_primary:
        if (s.startswith("Step")
            or "Secondary sub-process" in line
            or "performed by EventHandler" in line):
          in_primary = False
          section = None
          pending_pid = None
          continue

      if not in_primary:
        continue

      # --- section markers within the primary block ---
      if "--- incoming:" in line:
        section = 'incoming'; pending_pid = None; continue
      if "--- outgoing:" in line:
        section = 'outgoing'; pending_pid = None; continue
      if "--- intermediates:" in line:
        section = 'intermediates'; pending_pid = None; continue

      # The "------" divider closes the primary block (it appears right
      # after the outgoing section).
      if s.startswith("---") and set(s) <= {"-"}:
        in_primary = False
        section = None
        pending_pid = None
        continue
      # the "------------------..." full divider:
      if set(s) <= {"-"} and len(s) > 10:
        in_primary = False
        section = None
        pending_pid = None
        continue

      # only collect incoming/outgoing
      if section not in ('incoming', 'intermediates', 'outgoing'):
        pending_pid = None
        continue

      if pending_pid is None:
        h = hdr_re.match(line)
        if h:
          pending_pid = int(h.group(1))
        continue
      else:
        mm = mom_re.match(line)
        if mm:
          px, py, pz, e = (float(mm.group(k)) for k in range(1, 5))
          events[cur_evt][section].append((pending_pid, px, py, pz, e))
        pending_pid = None
        continue

    return events
    
  def fill_herwig_parton_info(self):
    print("Filling parton info from log file for Herwig events...")

    STATUS = {'incoming': 0, 'intermediates': 1, 'outgoing': 2}

    # ---- fill from the parsed log ----
    log = self.parse_log(self.herwig_log_file)   # {event_number: {'incoming':[...], 'outgoing':[...]}}

    # evt_no counts from 1, so in the Fill use evt_no-1
    for evt_no in sorted(log.keys()):
      for sect in ('incoming', 'intermediates', 'outgoing'):
      # outgoing = log[evt_no]['outgoing']

        for (this_pid, this_px, this_py, this_pz, this_e) in log[evt_no][sect]:
          self.t_parentparton.Fill(self.run_number, evt_no-1, float(this_px), float(this_py), float(this_pz), float(this_e), float(this_pid), STATUS[sect])
          if evt_no < 10:
            print(f"Event {evt_no}: Parton {this_pid} with momentum ({this_px}, {this_py}, {this_pz}), energy {this_e}")

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
  parser.add_argument('--jse', help='include additional tree for JSE information', action='store_true', default=False)
  parser.add_argument('--add-herwig-parton', help='save the herwig outgoing parton information', action='store_true', default=False)
  parser.add_argument('-l', '--herwig-log', help='path to the Herwig log file', type=str, required=False)
  args = parser.parse_args()
  
  converter = HepMC2antuple(input = args.input, output = args.output, as_data = args.as_data, hepmc = args.hepmc, nev = args.nev, gen = args.gen, no_progress_bar = args.no_progress_bar, include_parton = args.include_parton, include_D0 = args.include_D0, for_jse = args.jse, include_herwig_parton = args.add_herwig_parton, herwig_log_file = args.herwig_log)
  converter.main()
