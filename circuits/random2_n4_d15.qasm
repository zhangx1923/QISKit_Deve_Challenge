OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.43648765560089,1.04308654553329,-0.0346763093672431) q[3];
u3(1.79272781379865,-0.469340225562004,-4.30895143078715) q[0];
cx q[0],q[3];
u1(0.807126042367562) q[3];
u3(-1.41157909903435,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.357083910211695,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.84707995568474,1.27058697821318,-2.61723921552426) q[3];
u3(0.778846747965133,-1.51794368873874,-4.64845679784575) q[0];
u3(0.282350478838013,1.48674346777729,-1.72175657066643) q[2];
u3(0.796396694117834,0.927007738408398,-1.49814969274823) q[1];
cx q[1],q[2];
u1(0.689245003685467) q[2];
u3(-0.118816343671922,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.10928393103511,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.75858427410040,-0.0914458386246630,2.13190614272103) q[2];
u3(0.683207299856402,3.30953407098702,-2.25877038164811) q[1];
u3(2.18747922969139,0.330700291178790,-3.03069050409564) q[3];
u3(1.24021559152504,-2.91420293088227,3.06466988694822) q[1];
cx q[1],q[3];
u1(1.03665276566753) q[3];
u3(-3.33447098650324,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.86995853587272,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.57800848513107,2.53404014221114,-2.94886897478692) q[3];
u3(2.19669340249006,0.736074664423855,-2.10729151855850) q[1];
u3(2.58691673532544,1.43103888253325,1.58728193605967) q[0];
u3(1.30622563677918,-1.78693809146591,-1.99910581097280) q[2];
cx q[2],q[0];
u1(1.57830184925402) q[0];
u3(0.217377551204775,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.870738753398738,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.47405200846690,1.42019526956850,3.04567009938378) q[0];
u3(1.21285986809516,2.19240301505786,1.56841124180427) q[2];
u3(0.853338729863032,-2.69777336515898,3.52737515163611) q[3];
u3(1.51533996415222,0.338620161180784,-0.772574960924362) q[2];
cx q[2],q[3];
u1(1.49981270601618) q[3];
u3(-3.36981049791214,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.80642162184498,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.70100553060926,0.619979299617800,-2.23964603305248) q[3];
u3(2.25332717377737,-1.32189948615223,0.694700197409135) q[2];
u3(2.29433977054874,0.346857845428249,1.92514439372270) q[0];
u3(1.41363932307165,-0.760690429866063,-0.967447278522731) q[1];
cx q[1],q[0];
u1(3.76296166185487) q[0];
u3(-1.57293283042003,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.92309209762198,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.230824475553494,-0.533710833664831,1.89851885447012) q[0];
u3(1.61123141935146,-2.15649804296908,-1.88019082352330) q[1];
u3(2.63115375563757,1.36360790401595,1.10369311263311) q[2];
u3(1.31485849652460,-5.15241000366698,-0.0959663043751413) q[3];
cx q[3],q[2];
u1(0.634804388083371) q[2];
u3(-0.167466867369260,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.25358942827290,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.24248121752134,-3.07234320843303,1.48120382419640) q[2];
u3(2.02956224237840,-4.46393646411280,0.664197900498436) q[3];
u3(1.01933932208765,2.29865503091890,-0.252634070126404) q[0];
u3(1.11809643761879,1.53508597581452,-1.38232100798420) q[1];
cx q[1],q[0];
u1(1.90651275373293) q[0];
u3(0.929800924467778,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.50989718967277,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.85493045584862,3.33869946051487,-2.73771159461304) q[0];
u3(1.80037696863406,-1.87035726505058,-1.07631707022387) q[1];
u3(1.41565782617174,2.56697849810870,-2.76115909183858) q[0];
u3(1.07925178575965,1.35933474452584,-1.89843210354535) q[1];
cx q[1],q[0];
u1(-1.11239563745613) q[0];
u3(0.580566077150669,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.39220047443392,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.39209857586958,-1.00494416180901,-1.99973474965643) q[0];
u3(1.54827331524403,0.157327183074905,-3.47691057073673) q[1];
u3(1.63595750182911,-2.14196833455886,-0.682871468331611) q[2];
u3(1.75146182721431,-4.29618578256176,-0.539602573759441) q[3];
cx q[3],q[2];
u1(1.60796374010852) q[2];
u3(0.483247601534672,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.02821130241258,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.00218952242770,-1.61537867707290,2.82958979073715) q[2];
u3(0.268967774868600,-3.00843896301419,1.68396805743862) q[3];
u3(1.12883472559226,1.62303896338350,-1.97423177163439) q[3];
u3(0.0974496474999723,1.45134068297776,-3.90623931234668) q[0];
cx q[0],q[3];
u1(0.0488856997284102) q[3];
u3(1.46198125029260,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.99179098271383,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.49173025354070,-0.730483827487904,3.00687024651748) q[3];
u3(1.83906191358810,4.82456582948098,-0.246407224120496) q[0];
u3(0.601048529356623,-2.10450849576587,1.28438731481173) q[1];
u3(0.368076248078571,1.30412323116647,-3.35195743316881) q[2];
cx q[2],q[1];
u1(1.32212177033028) q[1];
u3(-0.942456810561492,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.551031431380008,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.35928948470344,-0.937222433789440,2.44679089057531) q[1];
u3(1.29455537762805,-0.749700035761657,0.846764028925580) q[2];
u3(1.68200961890753,2.07146810279509,0.0772625635175112) q[0];
u3(1.42325135967875,0.167351615327258,-4.47950060331580) q[3];
cx q[3],q[0];
u1(3.66688321133934) q[0];
u3(-3.96311154344095,0.0,0.0) q[3];
cx q[0],q[3];
u3(-1.11311051905676,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.74724512451973,0.909020312515167,-0.433767365748761) q[0];
u3(1.15904973392107,0.837510026311455,1.22474580820267) q[3];
u3(2.15552880785065,1.33856874912962,-1.22244806791489) q[2];
u3(2.03167999472918,1.81384168791203,-3.85460061979433) q[1];
cx q[1],q[2];
u1(3.34646161196013) q[2];
u3(-1.44319666023479,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.44883305693892,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.06798044837966,-0.112323373812411,0.895194496996776) q[2];
u3(1.15680989005080,-0.782459128731495,0.809206384123998) q[1];
u3(1.47039588608306,1.46300840671209,-2.44482037035239) q[3];
u3(0.773198702390972,2.43201867633984,-3.66737669258872) q[1];
cx q[1],q[3];
u1(3.16127058628530) q[3];
u3(-1.55350008461833,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.42608233299325,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.54535378006286,0.500766521459940,0.705168322642519) q[3];
u3(1.76572780941081,-1.89663997617040,-1.76791346781275) q[1];
u3(2.78289627025153,0.642332145052744,2.38773656443779) q[2];
u3(1.55981701719284,3.65990685319712,2.58237340032084) q[0];
cx q[0],q[2];
u1(0.478482486818556) q[2];
u3(-1.07041428068099,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.47351078232581,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.75922187146367,-0.555728104102114,4.32598297991673) q[2];
u3(0.866223275822000,-3.84000075071824,0.113387824393470) q[0];
u3(0.565109263697115,2.44076974337190,-2.80488330015857) q[1];
u3(0.681913625902051,0.501849146889451,-1.22039553346351) q[3];
cx q[3],q[1];
u1(1.71615341130245) q[1];
u3(-2.32177608416764,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.32486918776437,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.86005140956648,2.24607344067887,-0.116421730437771) q[1];
u3(1.12855731560478,-2.63478522778026,-3.07005087759846) q[3];
u3(1.17187244500579,1.38991939032619,-3.64126594618362) q[0];
u3(1.71505324879667,2.79149082611421,-2.82824281875162) q[2];
cx q[2],q[0];
u1(1.63293063757167) q[0];
u3(0.101671956771790,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.03484934727739,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.22329745480349,-0.641264484164529,1.82472612019675) q[0];
u3(1.88873804393595,-2.04438339414478,3.95445135548327) q[2];
u3(2.59283119685718,1.60940989623654,1.26582632908668) q[2];
u3(0.725576358485949,-4.84918355115160,-0.402727734346947) q[0];
cx q[0],q[2];
u1(4.45085716835124) q[2];
u3(-3.80440174639893,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.901883968666777,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.55375158636835,-1.71222372829876,2.63212269350724) q[2];
u3(1.53709093354942,2.27436652340238,2.10902139729597) q[0];
u3(1.74686152806991,-0.0863038934511951,0.673052827676468) q[3];
u3(2.84639408971917,-0.555225320069605,-1.37750974298829) q[1];
cx q[1],q[3];
u1(2.25957627305463) q[3];
u3(-1.61158648008245,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.413508298802820,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.10049092473916,0.0810447132600634,-2.83467470096502) q[3];
u3(0.985552048397292,0.248619943528658,-1.29384201497866) q[1];
u3(2.32435366132321,0.435987719776777,-2.14261924244206) q[3];
u3(2.99320596827620,4.12923583177183,0.222355976386588) q[1];
cx q[1],q[3];
u1(1.90848385790547) q[3];
u3(-3.12756985755037,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.683947026611847,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.77308048275306,-0.237460570377725,3.92873995908668) q[3];
u3(2.43094355983356,-0.759578240232429,-4.10121055209072) q[1];
u3(1.27155084893435,2.23109689465615,0.760027854805262) q[0];
u3(0.280818939513258,-0.404660890841621,-2.79089880378194) q[2];
cx q[2],q[0];
u1(0.707240147187241) q[0];
u3(-1.21617731041664,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.158734426075972,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.55154641367746,-0.802672341192610,0.925624601594618) q[0];
u3(2.39845679337181,-3.31428878357691,2.34775435641960) q[2];
u3(0.681287601865158,-0.127956683709583,0.625340915065711) q[2];
u3(1.19769004105515,0.117220169702402,-1.16568186385585) q[3];
cx q[3],q[2];
u1(1.73481078159851) q[2];
u3(-2.49065284577154,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.0767846210611665,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.23678229581975,1.82574343926760,1.64734670348324) q[2];
u3(1.17837921690191,2.52803396176859,-0.0174817407294314) q[3];
u3(0.789326966533533,-3.38024298887324,1.37539407299943) q[1];
u3(2.04306686294892,-5.67743984374995,-0.560628355343779) q[0];
cx q[0],q[1];
u1(2.13738274749837) q[1];
u3(-2.56358201974333,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.755864169043466,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.27891771749305,-2.78181633804672,-0.108621486039592) q[1];
u3(2.27324476400678,4.35659260314657,0.960578922897595) q[0];
u3(1.82618116217524,-3.73672577646646,0.644770511493868) q[1];
u3(1.61065165674423,0.228073927186023,2.62593237409969) q[2];
cx q[2],q[1];
u1(3.23003686142025) q[1];
u3(-1.97994998197935,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.19194728261130,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.76346547694907,0.771279475730960,-1.77204976069979) q[1];
u3(1.69539575153291,-1.95284884764523,3.69845222712610) q[2];
u3(0.195034964029936,0.658187088395147,-0.691846844484767) q[3];
u3(0.793882158972571,-3.41017699531730,1.14707657422418) q[0];
cx q[0],q[3];
u1(1.96878281088123) q[3];
u3(-3.04916339557719,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.56164656155306,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.37573517777007,-1.12688713694873,4.90325858228524) q[3];
u3(1.45302749792507,0.498732086420289,-3.44982212898499) q[0];
u3(1.02189294235964,-0.699689666192566,1.93772143279530) q[2];
u3(0.673655205160617,-0.596709206935527,-1.30321288692037) q[1];
cx q[1],q[2];
u1(3.05996584858771) q[2];
u3(-1.69424455494680,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.63881197234098,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.645554648979986,2.33565603458187,-2.73478091771765) q[2];
u3(0.942807105322894,0.330789592774302,-5.47360637603049) q[1];
u3(0.435396634033158,0.419097403680764,-0.366530147229707) q[0];
u3(0.135701642273073,-2.66933044639886,-0.0551810279013938) q[3];
cx q[3],q[0];
u1(1.80829704024869) q[0];
u3(-2.17415515253706,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.68447513549860,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.21386671512401,2.86382401755290,-0.335578369420175) q[0];
u3(1.36896902590972,2.92316489530154,-2.25482873632468) q[3];
u3(0.580845850370061,-1.13830569672375,0.123264306778335) q[0];
u3(1.92910392445560,-2.84679501316371,0.498301547375841) q[3];
cx q[3],q[0];
u1(0.867046435688347) q[0];
u3(-1.30857294664019,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.48321945414386,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.03866656831912,-0.770076130335757,-0.920444958645798) q[0];
u3(1.11337765894137,0.667689463007903,1.65744030314613) q[3];
u3(1.96800647045591,1.11256948039885,0.525878416691979) q[1];
u3(0.550306848473861,-3.21536248740125,-0.859022859168135) q[2];
cx q[2],q[1];
u1(2.82945311648530) q[1];
u3(-1.57934610558665,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.31556493983325,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.75938653025483,1.50388928964035,-0.516796967217240) q[1];
u3(2.42050716222483,2.58265223867827,-1.42751944124550) q[2];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];