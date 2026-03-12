#!/bin/bash

# should be sitting in pyjetty/pyjetty/alice_analysis/generation/pythia/config/

for i in {1..10}; do
    echo "$i"
    cd $i

    cp settings_HF_onlypromptDtoKpi.cmnd settings_HF_onlyprompt.cmnd
    # sed -i 's|^PDF:pSet*$|# PhaseSpace:etaMin = -0.9\n# PhaseSpace:etaMax = 0.9|' settings_5020_limiteta.cmnd
    
    cd ..
done