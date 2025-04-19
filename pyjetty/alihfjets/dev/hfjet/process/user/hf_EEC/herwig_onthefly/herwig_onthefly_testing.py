

# need to module load herwig_with_deps for this file to work!!

import Herwig # seems like this works but how?
import subprocess




class HerwigRunner:
    def __init__(self, config="1/LHC_13000_MPI.run", events=100, debug=1):
        self.config = config
        self.events = events
        self.debug = debug

    def run(self):
        self.full_config_file = "/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/herwig/run/" + self.config
        command = ["Herwig", "run", self.full_config_file, "-N", str(self.events), "-d", str(self.debug)]
        
        # output_settings = ["--set /Herwig/Analysis/HepMCFile:Format=Root", "--set /Herwig/Analysis/HepMCFile:Filename=output.root"]
        # --verbose 2
        # command.extend(output_settings)
        print("COMMAND!", command)

        result = subprocess.run(command, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("Error:", result.stderr)

# Usage:
runner = HerwigRunner(events=10)
runner.run()



# herwig = Herwig.

# pythia = pythia8.Pythia()
# pythia.readString(s)
# if pythia.init():
#     print ('[i] pythia initialized with', config_strings)
#     # return pythia