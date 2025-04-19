import os
import argparse
import ROOT
import pyhepmc_ng

import hepmc2antuple_base


class HepMC_testing_herwig(hepmc2antuple_base.HepMC2antupleBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, **kwargs):
    super(HepMC_testing_herwig, self).__init__(**kwargs)
    self.init()
    print(self)
    
  #---------------------------------------------------------------
  def main(self):

    self.hepmc = 2

    # Open the HepMC file
    filename = "/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alihfjets/dev/hfjet/process/user/hf_EEC/herwig_onthefly/LHC_13000_MPI_10events.hepmc"
    # with open(filename, "r") as f:
        # reader = pyhepmc.ReaderAscii(f)
    print("hello")
    if self.hepmc == 3:
      input_hepmc = pyhepmc_ng.ReaderAscii(self.input)
    if self.hepmc == 2:
      input_hepmc = pyhepmc_ng.ReaderAsciiHepMC2(self.input)

    if input_hepmc.failed():
      print ("[error] unable to read from {}".format(self.input))
      sys.exit(1)

    # event_hepmc = pyhepmc.GenEvent()
    event_hepmc = pyhepmc_ng.GenEvent()
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

  def fill_event(self, event_hepmc):
    print("here!")

    # Loop through each event
    print("Event", event_hepmc.event_number, "has", len(event_hepmc.particles), "particles.")

    # Loop through particles
    particle_counter = 0
    finalstate_particles_counter = 0
    statuscodes_hist = ROOT.TH1I("statuscodes_hist", "Status Codes;Status Code;Counts", 20, 0, 20)
    pidcodes_hist = ROOT.TH1I("pidcodes_hist", "PID < 30 Codes;PID;Counts", 30, 0, 30)
    for particle in event_hepmc.particles:
      pid = particle.pid  # PDG ID
      # print("PID", pid)
      if particle_counter == 1 and self.ev_id == 0:
        statuscodes_hist.Fill(particle.status)
        pidcodes_hist.Fill(pid)
      
      # print(f"Particle: PDG {pid}, Parent(s): {parents}")
      if particle.status == 4:
        print("INGOING PARTICLE")
        [print(child.pid, child.status) for child in particle.children]
      if particle.status == 1:
        self.get_parent(particle, 0) #recursive function to print all generations of parents
        # for parent in particle.parents:
        #   print("  parent", parent.pid, parent.status)
        finalstate_particles_counter += 1
        print("#######")
      if finalstate_particles_counter == 5:
        break

      # if (particle_counter == 10):
        # break
      particle_counter += 1
    
    if particle_counter == 1 and self.ev_id == 0:
      c1 = ROOT.TCanvas("c1", "c1", 800, 600)
      statuscodes_hist.Draw()
      c1.SaveAs("statuscodes.pdf")

      c2 = ROOT.TCanvas("c2", "c2", 800, 600)
      pidcodes_hist.Draw()
      c2.SaveAs("pidcodes.pdf")

  def get_parent(self, particle, reverse_generation):
    print("  parent", particle.pid, "status", particle.status, "index", particle.id, "rev gen", reverse_generation)

    # base case
    # if (parent.status == )
    if (len(particle.parents) == 0):
      return
    else:
      for p in particle.parents:
        # print("index? ", p)
        self.get_parent(p, reverse_generation+1)
    

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
  
  converter = HepMC_testing_herwig(input = args.input, output = args.output, as_data = args.as_data, hepmc = args.hepmc, nev = args.nev, gen = args.gen, no_progress_bar = args.no_progress_bar, include_parton = args.include_parton, include_D0 = args.include_D0)
  converter.main()