#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=1245121
OUTPUT_DIR=/rstorage/alice/AnalysisResults/preeti/ang/$JOB_ID

# Merge separate subsets, since otherwise it is too large for hadd
RUNLIST_LHC16d="000252235  000252271  000252317  000252322  000252326  000252332  000252368  000252371  000252375 000252248  000252310  000252319  000252325  000252330  000252336  000252370  000252374"
for RUN in $RUNLIST_LHC16d
do
   FILE_DIR=$OUTPUT_DIR/LHC16d/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16e="000252858 000252867 000253437 000253478 000253481 000253482 000253488 000253517 000253529 000253530 000253563 000253589 000253591"
for RUN in $RUNLIST_LHC16e
do
   FILE_DIR=$OUTPUT_DIR/LHC16e/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done
RUNLIST_LHC16g="000254128 000254147 000254149 000254174 000254175 000254178 000254193 000254199 000254204 000254205 000254293 000254302 000254303 000254304 000254330 000254331 000254332"

for RUN in $RUNLIST_LHC16g
do
   FILE_DIR=$OUTPUT_DIR/LHC16g/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16h="000254418 000254419 000254422 000254604 000254606 000254621 000254629 000254630 000254632 000254640 000254644 000254646 000254648 000254649 000254651 000254652 000254653 000254654 000254983 000254984 000255079 000255082 000255085 000255086 000255091 000255111 000255154 000255159 000255162 000255167 000255171 000255173 000255174 000255176 000255177 000255180 000255181 000255182 000255240 000255242 000255247 000255248 000255249 000255251 000255252 000255253 000255255 000255256 000255275 000255276 000255280 000255283 000255350 000255351 000255352 000255398 000255402 000255407 000255415 000255418 000255419 000255420 000255421 000255440 000255442 000255447 000255463 000255465 000255466 000255467 000255469"
for RUN in $RUNLIST_LHC16h
do
   FILE_DIR=$OUTPUT_DIR/LHC16h/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16j="000256219 000256223 000256225 000256227 000256228 000256231 000256281 000256282 000256283 000256284 000256287 000256289 000256290 000256292 000256295 000256297 000256299 000256302 000256307 000256309 000256311 000256356 000256361 000256362 000256363 000256364 000256365 000256366 000256368 000256371 000256372 000256373 000256415 000256417 000256418"
for RUN in $RUNLIST_LHC16j
do
   FILE_DIR=$OUTPUT_DIR/LHC16j/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16k="000256941 000256942 000256944 000257011 000257012 000257071 000257077 000257080 000257082 000257084 000257086 000257092 000257095 000257100 000257136 000257137 000257138 000257139 000257141 000257144 000257204 000257206 000257209 000257224 000257260 000257318 000257320 000257322 000257330 000257358 000257364 000257433 000257457 000257468 000257474 000257487 000257488 000257490 000257491 000257492 000257530 000257531 000257537 000257539 000257540 000257541 000257560 000257561 000257562 000257566 000257587 000257588 000257590 000257592 000257594 000257595 000257601 000257604 000257605 000257606 000257630 000257632 000257635 000257636 000257642 000257644 000257682 000257684 000257685 000257687 000257688 000257689 000257691 000257692 000257694 000257697 000257724 000257725 000257727 000257733 000257734 000257735 000257737 000257754 000257757 000257765 000257773 000257797 000257798 000257799 000257800 000257803 000257804 000257850 000257851 000257853 000257855 000257901 000257912 000257932 000257936 000257937 000257939 000257957 000257960 000257963 000257979 000257986 000257989 000257992 000258003 000258008 000258012 000258014 000258017 000258019 000258039 000258041 000258042 000258045 000258049 000258053 000258059 000258060 000258062 000258063 000258107 000258108 000258109 000258113 000258114 000258117 000258178 000258197 000258198 000258202 000258203 000258204 000258256 000258257 000258258 000258270 000258271 000258273 000258274 000258278 000258280 000258299 000258301 000258302 000258303 000258306 000258307 000258332 000258336 000258359 000258387 000258388 000258391 000258393 000258426 000258452 000258454 000258456 000258477 000258498 000258499 000258537"
for RUN in $RUNLIST_LHC16k
do
   FILE_DIR=$OUTPUT_DIR/LHC16k/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16l="000258962 000258964 000259086 000259088 000259090 000259091 000259096 000259099 000259117 000259118 000259162 000259164 000259204 000259257 000259261 000259263 000259264 000259269 000259270 000259271 000259272 000259273 000259274 000259302 000259303 000259305 000259307 000259334 000259336 000259339 000259340 000259341 000259342 000259378 000259381 000259382 000259388 000259389 000259394 000259395 000259396 000259473 000259477 000259747 000259748 000259750 000259751 000259752 000259756 000259781 000259788 000259789 000259822 000259841 000259842 000259860 000259866 000259867 000259868 000259888"
for RUN in $RUNLIST_LHC16l
do
   FILE_DIR=$OUTPUT_DIR/LHC16l/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16o="000262424 000262425 000262426 000262428 000262430 000262450 000262451 000262487 000262489 000262490 000262492 000262528 000262532 000262533 000262537 000262563 000262567 000262568 000262569 000262570 000262571 000262572 000262574 000262578 000262583 000262593 000262594 000262624 000262628 000262632 000262635 000262705 000262706 000262708 000262713 000262717 000262719 000262723 000262725 000262727 000262760 000262768 000262776 000262777 000262778 000262841 000262842 000262844 000262847 000262849 000262853 000262855 000262858 000263331 000263332 000263487 000263490 000263496 000263497 000263529 000263647 000263652 000263654 000263657 000263662 000263663 000263682 000263690 000263691 000263737 000263738 000263739 000263741 000263743 000263744 000263784 000263785 000263786 000263787 000263790 000263792 000263793 000263803 000263810 000263863 000263866 000263905 000263916 000263917 000263920 000263923 000263977 000263978 000263979 000263981 000263984 000263985 000264033 000264035"
for RUN in $RUNLIST_LHC16o
do
   FILE_DIR=$OUTPUT_DIR/LHC16o/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC16p="000264076 000264078 000264082 000264085 000264086 000264109 000264110 000264129 000264137 000264138 000264139 000264164 000264168 000264188 000264190 000264194 000264197 000264198 000264232 000264233 000264235 000264238 000264259 000264260 000264261 000264262 000264264 000264265 000264266 000264267 000264273 000264277 000264279 000264281 000264305 000264306 000264312 000264336 000264341 000264345 000264346 000264347"
for RUN in $RUNLIST_LHC16p
do
   FILE_DIR=$OUTPUT_DIR/LHC16p/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done



RUNLIST_LHC17e="000270822  000270824  000270827  000270828  000270830"
for RUN in $RUNLIST_LHC17e
do
   FILE_DIR=$OUTPUT_DIR/LHC17e/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17f="000270854  000270855  000270856  000270861  000270865"
for RUN in $RUNLIST_LHC17f
do
   FILE_DIR=$OUTPUT_DIR/LHC17f/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17h="000271870  000272036  000272100  000272155  000272388  000272413  000272468  000272607  000272747  000272783  000272870  000272933  000272983  000273101 000271871  000272038  000272101  000272156  000272389  000272414  000272469  000272608  000272749  000272784  000272871  000272934  000272985  000273103 000271873  000272039  000272123  000272194  000272394  000272417  000272521  000272610  000272760  000272828  000272873  000272935  000273009 000271880  000272040  000272151  000272335  000272395  000272461  000272574  000272620  000272762  000272829  000272880  000272939  000273010 000271886  000272042  000272152  000272340  000272399  000272462  000272575  000272690  000272763  000272833  000272903  000272947  000273077 000272018  000272075  000272153  000272359  000272400  000272463  000272577  000272691  000272764  000272834  000272905  000272949  000273099 000272020  000272076  000272154  000272360  000272411  000272466  000272585  000272712  000272782  000272836  000272932  000272976  000273100"
for RUN in $RUNLIST_LHC17h
do
   FILE_DIR=$OUTPUT_DIR/LHC17h/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17i="000273591  000273825  000273889  000273985  000274092  000274212  000274266  000274271  000274281  000274355  000274385  000274389 000273592  000273885  000273918  000273986  000274125  000274258  000274268  000274276  000274283  000274357  000274386  000274442 000273653  000273886  000273942  000274058  000274147  000274263  000274269  000274278  000274329  000274360  000274387 000273654  000273887  000273943  000274063  000274148  000274264  000274270  000274280  000274352  000274364  000274388"
for RUN in $RUNLIST_LHC17i
do
   FILE_DIR=$OUTPUT_DIR/LHC17i/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17j="000274593  000274594  000274595  000274596  000274601  000274653  000274657  000274667  000274669  000274671"
for RUN in $RUNLIST_LHC17j
do
   FILE_DIR=$OUTPUT_DIR/LHC17j/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17k="000274690  000274821  000274889  000275149  000275188  000275324  000275372  000275456  000275559  000275648  000276102  000276166  000276259  000276435 000274708  000274822  000274978  000275150  000275239  000275326  000275401  000275457  000275612  000275650  000276104  000276169  000276290  000276437 000274736  000274877  000274979  000275151  000275245  000275328  000275404  000275459  000275617  000275661  000276105  000276170  000276292  000276438 000274801  000274878  000275067  000275173  000275246  000275332  000275406  000275467  000275621  000275664  000276108  000276177  000276294  000276439 000274802  000274882  000275068  000275174  000275247  000275333  000275443  000275471  000275622  000275847  000276135  000276178  000276297  000276462 000274803  000274883  000275073  000275177  000275283  000275360  000275448  000275472  000275623  000276097  000276140  000276205  000276302  000276506 000274806  000274884  000275075  000275180  000275314  000275361  000275452  000275515  000275624  000276098  000276141  000276230  000276348  000276507 000274815  000274886  000275076  000275184  000275322  000275369  000275453  000275558  000275647  000276099  000276145  000276257  000276351  000276508 "
for RUN in $RUNLIST_LHC17k
do
   FILE_DIR=$OUTPUT_DIR/LHC17k/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17l="000276551  000276672  000276971  000277087  000277188  000277310  000277417  000277531  000277723  000277801  000277870  000277952  000278127  000278216 000276552  000276674  000276972  000277091  000277189  000277312  000277418  000277534  000277725  000277802  000277876  000277987  000278130 000276553  000276675  000277015  000277117  000277193  000277314  000277470  000277536  000277745  000277805  000277897  000277988  000278158 000276556  000276762  000277016  000277121  000277194  000277360  000277472  000277537  000277746  000277834  000277898  000277989  000278164 000276557  000276916  000277017  000277155  000277196  000277383  000277473  000277574  000277747  000277836  000277899  000277991  000278165 000276608  000276917  000277037  000277180  000277197  000277384  000277476  000277575  000277749  000277841  000277900  000277996  000278166 000276644  000276920  000277073  000277181  000277256  000277385  000277477  000277576  000277794  000277842  000277903  000278121  000278167 000276669  000276967  000277076  000277182  000277257  000277386  000277478  000277577  000277795  000277845  000277904  000278122  000278189 000276670  000276969  000277079  000277183  000277262  000277389  000277479  000277721  000277799  000277847  000277907  000278123  000278191 000276671  000276970  000277082  000277184  000277293  000277416  000277530  000277722  000277800  000277848  000277930  000278126  000278215"
for RUN in $RUNLIST_LHC17l
do
   FILE_DIR=$OUTPUT_DIR/LHC17l/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17m="000278914  000278964  000279041  000279106  000279157  000279238  000279274  000279354  000279487  000279642  000279688  000279826  000280052  000280131 000278915  000278999  000279043  000279107  000279199  000279242  000279309  000279355  000279488  000279676  000279689  000279827  000280066  000280134 000278936  000279000  000279044  000279117  000279201  000279264  000279310  000279391  000279491  000279677  000279715  000279830  000280107  000280135 000278939  000279005  000279068  000279118  000279207  000279265  000279312  000279410  000279550  000279679  000279718  000279853  000280108  000280140 000278941  000279007  000279069  000279122  000279208  000279267  000279342  000279435  000279559  000279682  000279719  000279854  000280111 000278959  000279008  000279073  000279123  000279232  000279268  000279344  000279439  000279630  000279683  000279747  000279855  000280114 000278960  000279035  000279074  000279130  000279234  000279270  000279348  000279441  000279632  000279684  000279749  000279879  000280118 000278963  000279036  000279075  000279155  000279235  000279273  000279349  000279483  000279641  000279687  000279773  000280051  000280126"
for RUN in $RUNLIST_LHC17m
do
   FILE_DIR=$OUTPUT_DIR/LHC17m/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17o="000280282  000280374  000280448  000280613  000280679  000280763  000280847  000280994  000281080  000281242  000281444  000281569  000281741  000281916 000280284  000280375  000280490  000280634  000280681  000280764  000280848  000280996  000281081  000281243  000281446  000281574  000281750  000281918 000280285  000280403  000280499  000280636  000280705  000280765  000280849  000280997  000281179  000281244  000281449  000281580  000281751  000281920 000280286  000280405  000280518  000280637  000280706  000280766  000280854  000280998  000281180  000281271  000281450  000281581  000281753  000281928 000280290  000280406  000280519  000280639  000280729  000280767  000280856  000280999  000281181  000281273  000281475  000281583  000281754  000281931 000280310  000280412  000280546  000280645  000280753  000280768  000280880  000281032  000281189  000281275  000281477  000281592  000281755  000281932 000280312  000280415  000280547  000280647  000280754  000280786  000280897  000281033  000281190  000281277  000281509  000281633  000281756  000281939 000280348  000280419  000280550  000280648  000280755  000280787  000280936  000281035  000281191  000281301  000281511  000281705  000281892  000281940 000280349  000280443  000280551  000280650  000280756  000280792  000280940  000281036  000281212  000281321  000281557  000281706  000281893  000281953 000280350  000280445  000280574  000280671  000280757  000280793  000280943  000281060  000281213  000281415  000281562  000281707  000281894  000281956 000280351  000280446  000280581  000280673  000280761  000280842  000280947  000281061  000281240  000281441  000281563  000281709  000281895  000281961 000280352  000280447  000280583  000280676  000280762  000280844  000280990  000281062  000281241  000281443  000281568  000281713  000281915"
for RUN in $RUNLIST_LHC17o
do
   FILE_DIR=$OUTPUT_DIR/LHC17o/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC17r="000282528  000282545  000282573  000282579  000282606  000282608  000282618  000282622  000282651  000282667  000282671  000282676  000282700  000282703 000282544  000282546  000282575  000282580  000282607  000282609  000282620  000282629  000282666  000282670  000282673  000282677  000282702  000282704"
for RUN in $RUNLIST_LHC17r
do
   FILE_DIR=$OUTPUT_DIR/LHC17r/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done



RUNLIST_LHC18b="000285009  000285012  000285014  000285064  000285066  000285108  000285127  000285200  000285203  000285224  000285328  000285364  000285396 000285011  000285013  000285015  000285065  000285106  000285125  000285165  000285202  000285222  000285327  000285347  000285365"
for RUN in $RUNLIST_LHC18b
do
   FILE_DIR=$OUTPUT_DIR/LHC18b/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC18d="000285978  000286025  000286064  000286130  000286201  000286230  000286258  000286284  000286308  000286312  000286337  000286348 000285979  000286027  000286124  000286159  000286202  000286231  000286261  000286287  000286309  000286313  000286340  000286349 000285980  000286028  000286127  000286198  000286203  000286254  000286263  000286288  000286310  000286314  000286341  000286350 000286014  000286030  000286129  000286199  000286229  000286257  000286282  000286289  000286311  000286336  000286345"
for RUN in $RUNLIST_LHC18d
do
   FILE_DIR=$OUTPUT_DIR/LHC18d/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18e="000286380  000286428  000286482  000286509  000286568  000286592  000286661  000286799  000286809  000286852  000286877  000286911  000286932  000286937 000286426  000286454  000286502  000286511  000286569  000286633  000286695  000286801  000286846  000286874  000286907  000286930  000286933 000286427  000286455  000286508  000286567  000286591  000286653  000286731  000286805  000286850  000286876  000286910  000286931  000286936"
for RUN in $RUNLIST_LHC18e
do
   FILE_DIR=$OUTPUT_DIR/LHC18e/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18f="000287000  000287072  000287202  000287249  000287324  000287349  000287381  000287451  000287516  000287578  000287784  000287912  000287977 000287021  000287077  000287203  000287250  000287325  000287353  000287385  000287480  000287517  000287654  000287876  000287913 000287063  000287137  000287204  000287251  000287343  000287355  000287387  000287481  000287518  000287656  000287877  000287915 000287064  000287155  000287208  000287254  000287344  000287356  000287388  000287484  000287521  000287657  000287884  000287923 000287066  000287185  000287209  000287283  000287346  000287360  000287389  000287486  000287524  000287658  000287885  000287941 000287071  000287201  000287248  000287323  000287347  000287380  000287413  000287513  000287575  000287783  000287911  000287975"
for RUN in $RUNLIST_LHC18f
do
   FILE_DIR=$OUTPUT_DIR/LHC18f/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18g="000288619  000288640  000288642  000288644  000288650  000288687  000288689  000288690  000288743  000288748  000288750"
for RUN in $RUNLIST_LHC18g
do
   FILE_DIR=$OUTPUT_DIR/LHC18g/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC18h="000288804  000288806"
for RUN in $RUNLIST_LHC18h
do
   FILE_DIR=$OUTPUT_DIR/LHC18h/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done

RUNLIST_LHC18i="000288861  000288862  000288863  000288864  000288868  000288897  000288902  000288903  000288908  000288909"
for RUN in $RUNLIST_LHC18i
do
   FILE_DIR=$OUTPUT_DIR/LHC18i/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18j="000288943"
for RUN in $RUNLIST_LHC18j
do
   FILE_DIR=$OUTPUT_DIR/LHC18j/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18k="000289166  000289167  000289169  000289172  000289175  000289176  000289177  000289198  000289199  000289200  000289201"
for RUN in $RUNLIST_LHC18k
do
   FILE_DIR=$OUTPUT_DIR/LHC18k/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18l="000289240  000289253  000289280  000289309  000289366  000289374  000289466  000289574  000289634  000289666  000289811  000289849  000289879  000289940 000289241  000289254  000289281  000289353  000289367  000289426  000289468  000289576  000289657  000289721  000289814  000289852  000289880  000289941 000289242  000289275  000289300  000289354  000289368  000289444  000289493  000289577  000289658  000289732  000289816  000289854  000289884  000289943 000289243  000289276  000289303  000289355  000289369  000289462  000289494  000289582  000289659  000289757  000289817  000289855  000289928  000289965 000289247  000289277  000289306  000289356  000289370  000289463  000289521  000289625  000289660  000289775  000289818  000289856  000289931  000289966 000289249  000289278  000289308  000289365  000289373  000289465  000289547  000289632  000289664  000289808  000289830  000289857  000289935"
for RUN in $RUNLIST_LHC18l
do
   FILE_DIR=$OUTPUT_DIR/LHC18l/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18m="000290323  000290467  000290627  000290887  000291005  000291257  000291400  000291618  000291769  000292062  000292167  000292434  000292563  000292809 000290327  000290469  000290632  000290888  000291006  000291262  000291402  000291622  000291795  000292067  000292168  000292456  000292584  000292810 000290350  000290499  000290645  000290894  000291035  000291263  000291416  000291624  000291796  000292075  000292192  000292457  000292586  000292811 000290374  000290500  000290660  000290895  000291037  000291265  000291417  000291626  000291803  000292080  000292218  000292460  000292693  000292831 000290375  000290501  000290665  000290932  000291041  000291266  000291420  000291657  000291942  000292081  000292240  000292461  000292695  000292832 000290376  000290538  000290687  000290935  000291065  000291282  000291424  000291661  000291943  000292106  000292241  000292495  000292696  000292834 000290399  000290539  000290689  000290941  000291066  000291284  000291447  000291665  000291944  000292107  000292242  000292496  000292698  000292836 000290401  000290540  000290766  000290943  000291069  000291285  000291451  000291690  000291945  000292108  000292265  000292497  000292701  000292839 000290404  000290544  000290787  000290944  000291093  000291286  000291453  000291697  000291946  000292109  000292273  000292500  000292704 000290411  000290549  000290790  000290948  000291100  000291360  000291456  000291698  000291948  000292114  000292298  000292521  000292737 000290412  000290550  000290841  000290974  000291101  000291361  000291457  000291706  000291953  000292115  000292397  000292523  000292739 000290423  000290553  000290843  000290975  000291110  000291362  000291481  000291729  000291976  000292140  000292398  000292524  000292744 000290425  000290588  000290846  000290976  000291111  000291363  000291482  000291755  000291977  000292160  000292405  000292526  000292747 000290426  000290590  000290848  000290979  000291116  000291373  000291484  000291756  000291982  000292161  000292406  000292553  000292748 000290427  000290612  000290853  000290980  000291143  000291375  000291485  000291760  000292012  000292162  000292428  000292554  000292750 000290456  000290613  000290860  000291002  000291188  000291377  000291590  000291762  000292040  000292163  000292429  000292557  000292752 000290458  000290614  000290862  000291003  000291209  000291397  000291614  000291766  000292060  000292164  000292430  000292559  000292803 000290459  000290615  000290886  000291004  000291240  000291399  000291615  000291768  000292061  000292166  000292432  000292560  000292804"
for RUN in $RUNLIST_LHC18m
do
   FILE_DIR=$OUTPUT_DIR/LHC18m/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18n="000293357  000293359"
for RUN in $RUNLIST_LHC18n
do
   FILE_DIR=$OUTPUT_DIR/LHC18n/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18o="000293368  000293424  000293496  000293578  000293587  000293695  000293741  000293776  000293807  000293831  000293893 000293386  000293474  000293570  000293579  000293588  000293696  000293770  000293802  000293809  000293856  000293896 000293392  000293475  000293571  000293582  000293691  000293698  000293773  000293805  000293829  000293886  000293898 000293413  000293494  000293573  000293583  000293692  000293740  000293774  000293806  000293830  000293891"
for RUN in $RUNLIST_LHC18o
do
   FILE_DIR=$OUTPUT_DIR/LHC18o/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


RUNLIST_LHC18p="000294009  000294131  000294205  000294242  000294524  000294530  000294562  000294591  000294653  000294722  000294746  000294775  000294816  000294880 000294010  000294152  000294208  000294305  000294525  000294531  000294563  000294593  000294703  000294741  000294747  000294805  000294817  000294883 000294011  000294154  000294210  000294307  000294526  000294553  000294586  000294632  000294710  000294742  000294749  000294809  000294818  000294884 000294012  000294155  000294212  000294308  000294527  000294556  000294587  000294633  000294715  000294744  000294769  000294813  000294852  000294916 000294013  000294199  000294241  000294310  000294529  000294558  000294590  000294636  000294718  000294745  000294772  000294815  000294875  000294925"
for RUN in $RUNLIST_LHC18p
do
   FILE_DIR=$OUTPUT_DIR/LHC18p/$RUN
   FILES=$( find "$FILE_DIR" -name "*.root" )
   hadd -f -j 20 $FILE_DIR/AnalysisResultsIntermediate.root $FILES
done


FILES=$( find $OUTPUT_DIR/LHC*/ -name "AnalysisResultsIntermediate.root" )
hadd -f $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
