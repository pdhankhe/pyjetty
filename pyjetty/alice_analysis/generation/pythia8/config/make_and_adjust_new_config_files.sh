#!/bin/bash

# should be sitting in pyjetty/pyjetty/alice_analysis/generation/pythia/config/

for i in {1..20}; do
    echo "$i"
    cd $i

    cp settings_5020_CTEQ5L.cmnd settings_5020_limiteta.cmnd
    sed -i 's|^PDF:pSet*$|# PhaseSpace:etaMin = -0.9\n# PhaseSpace:etaMax = 0.9|' settings_5020_limiteta.cmnd
    
    cd ..
done