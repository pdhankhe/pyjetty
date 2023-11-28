#rerunning whole analysis, to be done after changing config file.
read -p "What directory?" dir
echo $dir

rm -r $dir/jet_axis
rm $dir/AnalysisResults_data.root
rm $dir/AnalysisResults_outputfd_1.root
rm $dir/AnalysisResults.root
rm $dir/fRoounfold_R0.4_1.root
rm $dir/STD-D/AnalysisResults.root

./generate_analysis_results_data_std-d.py -c configYE_std-d.yaml

mkdir -p $dir/STD-D
mkdir -p $dir/SD-D
mkdir -p $dir/WTA-D
mkdir -p $dir/STD-SD
mkdir -p $dir/STD-WTA
mkdir -p $dir/WTA-SD

mv 'AnalysisResults_Name: hsparse_R0.4_STD-D Title: hsparse_R0.4_STD-D.root' $dir/STD-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_SD-D Title: hsparse_R0.4_SD-D.root' $dir/SD-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_WTA-D Title: hsparse_R0.4_WTA-D.root' $dir/WTA-D/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_STD-SD Title: hsparse_R0.4_STD-SD.root' $dir/STD-SD/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_STD-WTA Title: hsparse_R0.4_STD-WTA.root' $dir/STD-WTA/AnalysisResults.root
mv 'AnalysisResults_Name: hsparse_R0.4_WTA-SD Title: hsparse_R0.4_WTA-SD.root' $dir/WTA-SD/AnalysisResults.root


cp $dir/STD-D/AnalysisResults.root $dir/AnalysisResults_data.root
#cp $dir/STD-WTA/AnalysisResults.root $dir/AnalysisResults_data.root
cp AnalysisResults_outputfd*.root $dir

./unfolding_and_run_analysis_std-d.py -c configYE_std-d.yaml

#next run the systematics
read -p "What would you like to name the systematics file?" sys_file_name
echo $sys_file_name

#cp $dir/jet_axis/main/fResult_R0.4_1.root $dir/systematics/$sys_file_name
cp $dir/jet_axis/main/Unfolded_obs/fResult_name_*R0.4_1.root $dir/systematics/$sys_file_name

#below is for feeddown variation only
#cp FixedDataBins25/STD-D/AnalysisResults.root FinalResult/Feeddown/AnalysisResults_$sys_file_name

echo "Done!"
