read -p "What directory is this going to? " workdir
echo $workdir

mkdir -p $workdir/STD-D
mkdir -p $workdir/SD-D
mkdir -p $workdir/WTA-D
mkdir -p $workdir/STD-SD
mkdir -p $workdir/STD-WTA
mkdir -p $workdir/WTA-SD

mv 'AnalysisResults_Name: hsparse_R0.4_STD-D Title: hsparse_R0.4_STD-D.root' $workdir/STD-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_SD-D Title: hsparse_R0.4_SD-D.root' $workdir/SD-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_WTA-D Title: hsparse_R0.4_WTA-D.root' $workdir/WTA-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_STD-SD Title: hsparse_R0.4_STD-SD.root' $workdir/STD-SD/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_STD-WTA Title: hsparse_R0.4_STD-WTA.root' $workdir/STD-WTA/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_WTA-SD Title: hsparse_R0.4_WTA-SD.root' $workdir/WTA-SD/AnalysisResults.root
