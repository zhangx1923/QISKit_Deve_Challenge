OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.72450403341951,1.96484631526122,0.365124301544484) q[2];
u3(2.51787540379101,0.700382060446488,-2.27166572185734) q[1];
cx q[1],q[2];
u1(3.04698172602739) q[2];
u3(-1.85769841224477,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.415710912816834,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.56605430673410,0.420043706654100,1.90011351245570) q[2];
u3(0.714874871671902,-4.04224494386314,-1.13404146498832) q[1];
u3(0.474534156043167,-2.72798530803052,2.68178267765189) q[6];
u3(0.983071692660739,-3.36442374884542,2.72281357370355) q[5];
cx q[5],q[6];
u1(0.896549772486353) q[6];
u3(-1.17280373633045,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.125706392817915,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.16208458847222,-0.892142455031846,4.53070942094252) q[6];
u3(1.01141387307122,-0.881341323454698,-0.977523048713226) q[5];
u3(2.14212971859839,-3.75922089804978,0.642928158424189) q[4];
u3(1.40939843740114,-0.252148226026578,4.05266247047053) q[12];
cx q[12],q[4];
u1(3.34571738070435) q[4];
u3(-1.43046447556079,0.0,0.0) q[12];
cx q[4],q[12];
u3(2.04007699424559,0.0,0.0) q[12];
cx q[12],q[4];
u3(2.74964907914770,-0.339101288432585,0.890588296625428) q[4];
u3(1.49173507513920,-2.21436928776669,0.492392627329938) q[12];
u3(2.34491494615768,-1.11863469310735,-1.50504916278157) q[7];
u3(1.40902547129515,-1.90221835408452,-3.03079986719439) q[11];
cx q[11],q[7];
u1(-0.134471253314691) q[7];
u3(0.250866232666424,0.0,0.0) q[11];
cx q[7],q[11];
u3(3.98507865572844,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.72447412016431,-1.46855543385543,3.86311682295452) q[7];
u3(2.64512997683537,0.209517433642496,2.20124720726130) q[11];
u3(1.57427760898106,0.628089860831575,-2.25200582491059) q[9];
u3(1.06582345547884,-4.23309547532454,1.47248274276768) q[10];
cx q[10],q[9];
u1(3.48121900501859) q[9];
u3(-4.28631682804127,0.0,0.0) q[10];
cx q[9],q[10];
u3(-0.430789187718660,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.359737385603240,-2.82944408937515,2.07187718323608) q[9];
u3(1.15027822380653,2.77919890850019,-3.10055440373885) q[10];
u3(2.07661134290604,0.507943759699585,-0.877282561442725) q[0];
u3(1.25029394414290,0.858967394397668,-4.55733161338739) q[8];
cx q[8],q[0];
u1(3.68778636491032) q[0];
u3(-1.06163195045901,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.58405205591255,0.0,0.0) q[8];
cx q[8],q[0];
u3(3.04012169770938,1.99411384875115,-1.42144421357350) q[0];
u3(0.208023215861157,0.884971717040342,-4.35510203406678) q[8];
u3(0.655410504112586,-2.18143940476110,4.07618498392909) q[2];
u3(1.99149658163063,2.08735419136541,-1.78322813435414) q[5];
cx q[5],q[2];
u1(4.19287934407347) q[2];
u3(-3.60573633969753,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.220162414094120,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.678096966488662,-2.63911137707428,-0.263783169653081) q[2];
u3(2.69337670720775,1.13958424303978,-0.165034087359970) q[5];
u3(2.38516420440626,-3.89114212677980,1.60647253948034) q[8];
u3(2.19417246280881,-0.901772388795683,2.77627049491842) q[4];
cx q[4],q[8];
u1(2.89022947615971) q[8];
u3(-1.99450560903948,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.664580239903469,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.899151631319052,3.11084049443014,0.00758531926642658) q[8];
u3(1.82623735529526,-0.602878671748195,2.54424778208934) q[4];
u3(0.288749238744017,-1.13329521534389,1.33729784133396) q[11];
u3(1.29140771887488,-1.01896980663380,-2.10564756839608) q[7];
cx q[7],q[11];
u1(1.22623894557713) q[11];
u3(-0.744471419045781,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.133806333804807,0.0,0.0) q[7];
cx q[7],q[11];
u3(0.746291395269636,-1.98849635539633,-1.30409369357064) q[11];
u3(1.59314450314381,3.65250668472850,2.14834426436528) q[7];
u3(1.13727203534133,-1.98853274854276,0.0422631707760144) q[0];
u3(2.44724476315931,-3.69660534216167,1.15722373995255) q[10];
cx q[10],q[0];
u1(3.48713444495443) q[0];
u3(-1.35675538998701,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.11280507156142,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.84263235078775,-0.819718695430475,-2.19280713337178) q[0];
u3(0.842175350295318,-0.543753683930468,4.58708859449749) q[10];
u3(1.96440687038563,2.13130329536114,0.346858975722995) q[12];
u3(1.83614165919752,1.11252673043599,-3.70939225922891) q[9];
cx q[9],q[12];
u1(1.28714824949688) q[12];
u3(-0.552610709786110,0.0,0.0) q[9];
cx q[12],q[9];
u3(2.06123906067202,0.0,0.0) q[9];
cx q[9],q[12];
u3(2.22098665567279,-1.77787614604852,-1.72702466539570) q[12];
u3(1.59010967671895,2.50954151272698,2.67003239772021) q[9];
u3(1.85576730860844,-0.646810910719354,2.33731504975280) q[1];
u3(1.40498625847025,-1.32872420256079,-1.56828305459252) q[3];
cx q[3],q[1];
u1(2.05224234393384) q[1];
u3(-2.80163444720643,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.749321992078105,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.27765358765586,-2.05800066623167,0.597784570731382) q[1];
u3(1.99571814586879,-1.94662920887726,-1.93227446451646) q[3];
u3(1.42272637240288,1.93994213346202,-3.10194906325911) q[1];
u3(2.17900072313575,2.42174801248501,-3.78522697189474) q[12];
cx q[12],q[1];
u1(2.33072823251463) q[1];
u3(0.0889202662347905,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.21292391718941,0.0,0.0) q[12];
cx q[12],q[1];
u3(2.98652074493436,-0.469409573130123,-0.374368637862461) q[1];
u3(1.78608405758679,-0.400313265202218,1.46530232963216) q[12];
u3(2.49154286133042,2.58204870975659,-0.0385929706626797) q[5];
u3(2.35538041969807,1.78533446502340,-4.10342214920740) q[9];
cx q[9],q[5];
u1(-1.05379121302032) q[5];
u3(0.0761414500656303,0.0,0.0) q[9];
cx q[5],q[9];
u3(3.41982907131732,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.83739541464831,1.56462050033178,-3.63323306458320) q[5];
u3(0.168971476910273,4.93649138381726,-1.02036796685022) q[9];
u3(1.12319459592576,0.875177193315303,-2.38388682409145) q[2];
u3(0.536512229460895,-3.67217834042772,1.81643340590137) q[7];
cx q[7],q[2];
u1(2.88151284185090) q[2];
u3(-2.18398080620900,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.66536078545352,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.05148850538134,1.79015753819476,-3.37043328327489) q[2];
u3(1.21774477791451,-0.957705191387269,1.90693576927189) q[7];
u3(1.03329492822832,-0.567907159912436,0.823887228247585) q[0];
u3(1.40635864532005,-1.13042919043166,-1.51952259661468) q[6];
cx q[6],q[0];
u1(2.53052627117796) q[0];
u3(0.282839578902245,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.27284520905255,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.992881295338489,3.31573721665085,-1.88708941452292) q[0];
u3(0.904795052132831,2.23318936177329,0.726223934881025) q[6];
u3(2.84554579937323,1.53098947303117,-1.23304586099313) q[11];
u3(1.71600253964211,0.447067510419114,-3.39786439821391) q[8];
cx q[8],q[11];
u1(3.50653921498059) q[11];
u3(-1.36307155851556,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.06085117720280,0.0,0.0) q[8];
cx q[8],q[11];
u3(0.983813951625952,1.88833242850315,-3.56907240286082) q[11];
u3(0.676575131678028,-1.18320906198407,2.85137307544616) q[8];
u3(2.25458534711560,-2.65300873152563,0.777233807637177) q[10];
u3(2.29697892460412,2.54667490388195,3.12449268662271) q[3];
cx q[3],q[10];
u1(2.57486835135551) q[10];
u3(-2.22842195133515,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.40905254113591,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.16070700823301,2.76377245886339,-2.45706907103854) q[10];
u3(0.925028507079684,2.50856757952548,-2.56201039382589) q[3];
u3(1.82742372134823,0.0355052522257902,2.32836222043784) q[7];
u3(1.60796857574104,-2.40380705336781,-2.46433056439173) q[11];
cx q[11],q[7];
u1(3.10972480943285) q[7];
u3(-1.50415744387311,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.45529805019258,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.45954387011828,2.98999220536428,-1.76513345365890) q[7];
u3(0.397778719994198,2.99377856918234,1.94437499252830) q[11];
u3(1.07755574746557,-1.87165391179753,3.68162342612338) q[8];
u3(2.08707222696390,2.16500066081417,-1.74695114433760) q[12];
cx q[12],q[8];
u1(1.22394041285536) q[8];
u3(-3.38496146840211,0.0,0.0) q[12];
cx q[8],q[12];
u3(2.07702764810416,0.0,0.0) q[12];
cx q[12],q[8];
u3(0.412738425557106,-0.192405014256822,1.40085580889953) q[8];
u3(2.07427687296697,-5.03022872677057,-0.310091711625908) q[12];
u3(2.11204724753920,-2.75703111293990,0.583916600113541) q[3];
u3(2.45042612283481,0.707257838723423,1.85911259045619) q[9];
cx q[9],q[3];
u1(1.21014142222494) q[3];
u3(-3.47918020847587,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.05548573559837,0.0,0.0) q[9];
cx q[9],q[3];
u3(2.46837281308615,-1.73041710872439,1.52300695469818) q[3];
u3(1.27991823480073,0.141462290743595,-4.13598392087616) q[9];
u3(0.805669695401110,2.12976017921022,-2.62302507660729) q[6];
u3(0.591273983680329,2.24725030386731,-2.64433755150872) q[10];
cx q[10],q[6];
u1(0.0644130256855795) q[6];
u3(-0.833192549729364,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.81227430178889,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.05613552483837,-0.205564256142599,1.45914627050458) q[6];
u3(1.45459435581872,-3.91305557112830,-2.12423654036595) q[10];
u3(0.105452980864452,-0.430012926406503,0.582448246106436) q[4];
u3(0.473417544344697,0.174660917474456,-2.20143895628119) q[0];
cx q[0],q[4];
u1(-0.123039191391088) q[4];
u3(-0.815308198974075,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.81239732354508,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.773263021500612,-0.513237946937686,-3.62785815267757) q[4];
u3(2.05174762291244,-4.40053470683522,-0.897192713423100) q[0];
u3(2.04514554632850,2.19672475419986,-3.31740889447821) q[1];
u3(1.03777311386306,-1.80388407104840,2.30152737122774) q[5];
cx q[5],q[1];
u1(3.30938320745356) q[1];
u3(-0.796777512416207,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.02092336760790,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.69113021139545,1.29509488416575,-2.24724144895141) q[1];
u3(2.65434149240452,3.61518626788442,-1.49592117425994) q[5];
u3(2.37574952110695,-0.460581905496079,1.71595299312923) q[4];
u3(1.95300743859207,-0.749822022707337,-0.897053982324150) q[12];
cx q[12],q[4];
u1(1.51810196557878) q[4];
u3(-3.41772795325892,0.0,0.0) q[12];
cx q[4],q[12];
u3(2.44349802189849,0.0,0.0) q[12];
cx q[12],q[4];
u3(0.793665323644644,-2.85271537895158,2.06148140824770) q[4];
u3(2.25756244833293,1.83672278986472,-0.694456896175217) q[12];
u3(0.723328526454941,2.97043079119100,-3.08577562021279) q[8];
u3(1.06549370199048,-2.87857625168919,1.58586783879871) q[0];
cx q[0],q[8];
u1(0.00417482829436477) q[8];
u3(-1.11633656062722,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.70388717538606,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.52828534246031,0.908514102653698,1.01329966199160) q[8];
u3(2.01348575466934,-1.76505808490984,1.79815602561559) q[0];
u3(2.10469720992470,3.86353837001765,-2.08802282867478) q[10];
u3(1.11421613358028,2.17586620979419,-2.29403910655882) q[2];
cx q[2],q[10];
u1(2.86940019248199) q[10];
u3(-2.07680842275002,0.0,0.0) q[2];
cx q[10],q[2];
u3(0.648986515715802,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.94720831656231,1.25491051944112,-2.12622823468515) q[10];
u3(0.854336202080594,2.66957838881906,3.46203036605375) q[2];
u3(0.382267840793313,2.54400618684488,-2.00964714742306) q[9];
u3(0.644217411328680,1.41713269390895,-1.91473996121416) q[5];
cx q[5],q[9];
u1(1.82828870568607) q[9];
u3(-2.08121017097896,0.0,0.0) q[5];
cx q[9],q[5];
u3(3.42689171559927,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.79317915723373,-0.704567206419407,2.87494410336277) q[9];
u3(0.774486064733588,3.02866401699560,-2.81683162441631) q[5];
u3(1.52815995860276,-0.653180845617568,-1.97206100780425) q[7];
u3(1.01050310367460,-3.80748224844837,0.623065148443998) q[1];
cx q[1],q[7];
u1(1.28643774385179) q[7];
u3(-0.660373110413082,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.73400133957721,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.06394855346834,0.647980876644152,1.44772345665032) q[7];
u3(1.36692417748544,2.88607776162096,-0.771319174751253) q[1];
u3(0.156643631208304,0.743790001177860,-0.756810387721567) q[3];
u3(0.959175584521559,0.983508226263707,-1.22297188614750) q[11];
cx q[11],q[3];
u1(2.20169631662192) q[3];
u3(0.190369485297835,0.0,0.0) q[11];
cx q[3],q[11];
u3(1.34914753679136,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.21371498563469,-1.56256146138682,3.87707270649730) q[3];
u3(1.54328749512592,-4.27778538154037,1.18252528370478) q[11];
u3(1.35646768202308,2.12743023666125,-3.60386635398435) q[5];
u3(0.995575864150197,-1.37595329084636,2.26044819138371) q[11];
cx q[11],q[5];
u1(0.665212592090595) q[5];
u3(-1.46529479515955,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.91853249203980,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.48742761429972,1.42173048036161,1.01991036371674) q[5];
u3(1.46187778228246,-1.56773107772714,0.196927557500632) q[11];
u3(0.670763147752169,1.45717001261611,-2.71765653386607) q[3];
u3(1.37792546750302,2.42587938217154,-3.24446102986173) q[9];
cx q[9],q[3];
u1(3.03618647737185) q[3];
u3(-2.49862109851021,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.15483366405161,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.86781619510075,-2.63155345985511,1.71109330814865) q[3];
u3(0.737821415146992,4.76177892358026,-1.13345155727700) q[9];
u3(2.93011869155243,0.244227031533161,-2.15653108682244) q[8];
u3(2.58175542249211,3.68824738533017,-1.29599630213085) q[2];
cx q[2],q[8];
u1(2.73145019782145) q[8];
u3(-1.76689395597790,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.644141687935208,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.54892946849420,-1.34018473801646,-1.53732394610251) q[8];
u3(2.71157285128501,1.89059108773686,2.93975396050995) q[2];
u3(2.61401422879718,-1.89065852842095,0.876882358943973) q[0];
u3(1.48620300378272,-2.08553032213218,-0.808681813995765) q[6];
cx q[6],q[0];
u1(0.887808255867380) q[0];
u3(-3.23240824183957,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.87877550298013,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.36603116444191,-3.34453243630176,1.76186463769416) q[0];
u3(0.238342001672565,1.07229631049480,-3.64394025648600) q[6];
u3(1.36347658960621,-1.57508084168530,1.00404192338838) q[10];
u3(2.11011691582322,-4.51852685913100,-0.647936556962486) q[4];
cx q[4],q[10];
u1(-1.01960138834088) q[10];
u3(-0.0483449607835504,0.0,0.0) q[4];
cx q[10],q[4];
u3(3.63553800778157,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.51842392511086,0.349660098509831,2.23947534609494) q[10];
u3(1.52454534516349,-1.92976609842698,-0.412787975002082) q[4];
u3(2.20608819709796,-2.87283629814206,3.36218562391245) q[1];
u3(0.610190212490079,0.729881946201937,0.947772719661195) q[7];
cx q[7],q[1];
u1(0.704102131169192) q[1];
u3(-1.16048398556733,0.0,0.0) q[7];
cx q[1],q[7];
u3(-0.139088547164785,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.02273345188885,-0.321217920413877,0.703647114489998) q[1];
u3(1.76031397184942,4.33242126832981,1.35776389540383) q[7];
u3(1.19316768341484,1.82514591405919,-3.85691559206378) q[3];
u3(2.30918326190555,2.36109799129767,-3.25695986485866) q[10];
cx q[10],q[3];
u1(2.10748042911836) q[3];
u3(-3.09395867757287,0.0,0.0) q[10];
cx q[3],q[10];
u3(1.77963674553069,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.83461434693238,2.50423836396027,-3.02518517065971) q[3];
u3(1.66067341192194,0.655604032135078,-3.50181437951462) q[10];
u3(1.77295750522696,-0.158209251915726,0.395181167545647) q[12];
u3(2.20126579434072,-0.794422907244387,-1.43117896710389) q[7];
cx q[7],q[12];
u1(0.0229556626240177) q[12];
u3(1.07567580805436,0.0,0.0) q[7];
cx q[12],q[7];
u3(3.45130894374091,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.820360691850368,-4.15995670876586,1.30872940098479) q[12];
u3(1.95039627697579,-0.00186236262870298,1.19439281679379) q[7];
u3(1.67550328195815,-0.698782497418479,-0.0881877352962357) q[11];
u3(2.42450358448092,-2.89589718059453,-0.665832907886121) q[1];
cx q[1],q[11];
u1(2.61270792168026) q[11];
u3(-1.60270316014364,0.0,0.0) q[1];
cx q[11],q[1];
u3(0.437676834764647,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.462069277995055,-3.23831547543782,-0.418235470901027) q[11];
u3(0.495607833506508,-3.09811037888591,1.51800128127031) q[1];
u3(2.32435625843035,2.75285339177575,-3.31270126065975) q[0];
u3(1.30382069854556,3.43142732839147,-2.70371813662600) q[5];
cx q[5],q[0];
u1(1.13196298000413) q[0];
u3(-0.467380962360090,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.43757472378267,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.50881309465610,0.701907067625619,-2.59411102138553) q[0];
u3(1.20832813311539,-0.834962651245555,0.00224131591302063) q[5];
u3(1.04965734785599,3.66063593261673,-0.864745852780737) q[6];
u3(2.13404495734571,2.63496742133453,-1.11510246083008) q[2];
cx q[2],q[6];
u1(3.26302044687884) q[6];
u3(-1.08244678755844,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.80690150118117,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.56990885138109,-0.782548314426136,0.0519344057767958) q[6];
u3(0.341379104732977,1.35507429846121,-0.721212676866117) q[2];
u3(1.43466510853418,-0.248407297901508,0.616484368385934) q[4];
u3(1.28166480423738,-0.652354617339185,-2.00227893068383) q[9];
cx q[9],q[4];
u1(2.98895696652559) q[4];
u3(-2.06586876423616,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.21066833707986,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.71074830689307,-0.620292881688671,1.01669762103636) q[4];
u3(1.70588990656780,4.01698862635214,0.847667002790282) q[9];
u3(2.36560149624072,-0.00326937845235498,2.80691924010329) q[7];
u3(1.97483423589516,-1.88580164202038,-1.60182650688641) q[6];
cx q[6],q[7];
u1(3.04561071561921) q[7];
u3(-2.79347192538965,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.13188904991235,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.304747976355741,-2.31350021398982,-0.641621811445198) q[7];
u3(1.91839779749383,-1.98789735289859,-2.54219505183479) q[6];
u3(2.27581438008180,-3.49291622625757,2.36410297436640) q[9];
u3(0.728074190516360,3.94945755590973,-2.06466162705526) q[1];
cx q[1],q[9];
u1(1.98578089716571) q[9];
u3(-2.42294349608523,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.28706814653489,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.00751952971102,-3.76535552126437,-0.210676998409133) q[9];
u3(1.49173039719960,-5.80514043508157,0.272664888869192) q[1];
u3(1.68883622077167,2.73111457724970,-3.08310534039450) q[2];
u3(1.33569713926169,3.48902432135381,-2.68818013128992) q[3];
cx q[3],q[2];
u1(2.35606983247271) q[2];
u3(0.0978430992368173,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.21704297940360,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.24549261796805,-0.469861383164361,2.36497553677451) q[2];
u3(2.14545974917755,-2.08578378022971,0.886384702430944) q[3];
u3(2.05023884708288,0.777025126066474,-2.35175077108272) q[4];
u3(2.52654334600636,-0.756368422823148,-5.24869941098968) q[5];
cx q[5],q[4];
u1(2.51177906222121) q[4];
u3(-2.93961806551797,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.19113761470440,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.97831796892786,-1.71293483391027,2.94153155937081) q[4];
u3(1.11525828347978,-2.83392859563084,-2.86962689604487) q[5];
u3(1.98320021646937,-0.516519902765480,-2.46192663001110) q[11];
u3(1.53961452518521,-4.89480475111779,1.07777955241703) q[12];
cx q[12],q[11];
u1(-0.165113586020382) q[11];
u3(0.902288777115512,0.0,0.0) q[12];
cx q[11],q[12];
u3(3.64391504436931,0.0,0.0) q[12];
cx q[12],q[11];
u3(2.53055777200191,0.152806195185089,-0.565521787653398) q[11];
u3(1.72246477820189,1.38737474012681,-1.59643952622952) q[12];
u3(1.93115230380544,2.32117179360091,-1.63097540752912) q[8];
u3(1.73185716624904,0.913495164199345,-0.435559297582088) q[0];
cx q[0],q[8];
u1(-0.428497778853634) q[8];
u3(1.05655264272555,0.0,0.0) q[0];
cx q[8],q[0];
u3(3.46491879600220,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.96375406823439,-3.23795873705531,2.74672556616060) q[8];
u3(1.25197969426906,-3.16995338732947,-1.41520543130224) q[0];
u3(1.03093760608964,1.14573316805246,-1.14384525829657) q[1];
u3(1.06020206263829,-4.42125871921884,0.993747827667409) q[7];
cx q[7],q[1];
u1(-0.241987266677530) q[1];
u3(-1.49334032313997,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.01585160782799,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.60657195726058,1.47251613252197,0.449994445728124) q[1];
u3(0.873493403127164,-1.78852532518305,1.42696939305402) q[7];
u3(1.82928810232856,-0.224667815195957,0.686226318998151) q[12];
u3(1.37442378572499,-0.694795607683125,-1.37418179685692) q[9];
cx q[9],q[12];
u1(2.09413717499770) q[12];
u3(-3.22078779997395,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.31350802209460,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.85146091870783,-2.72903462080577,1.19170595284887) q[12];
u3(1.27456924436730,3.53226479909893,-1.56290367964539) q[9];
u3(2.42991907402563,-0.244583019897893,-1.41159142551502) q[4];
u3(1.33186652546169,-3.55215223455597,0.744836725243042) q[3];
cx q[3],q[4];
u1(0.317035155180645) q[4];
u3(-0.791343519990512,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.88200915258654,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.37561016372740,-0.526677233238540,2.05121158509270) q[4];
u3(1.20522329192657,-2.00511595533843,-2.16064199229852) q[3];
u3(0.435973286652882,0.981982707405148,-1.35849118845073) q[5];
u3(0.575174470707341,0.600410473777208,-2.31954105973865) q[0];
cx q[0],q[5];
u1(2.33394985815594) q[5];
u3(-1.84014264322602,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.58976978576640,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.54317398993084,-4.04197825917793,1.38719236024835) q[5];
u3(0.522696313901333,2.82666689189403,-3.33115604202412) q[0];
u3(1.63327052868996,1.77573679218140,-2.48572255580891) q[10];
u3(1.18825389734151,-2.73183237972788,3.02858893241851) q[11];
cx q[11],q[10];
u1(1.74715590927779) q[10];
u3(-0.449595715250194,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.53359792504388,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.56532112913497,2.43822192350412,-0.505311541067511) q[10];
u3(1.91113637970962,5.85995544532676,-0.367880854419524) q[11];
u3(0.780571037238984,-0.715632118100680,1.94457311559088) q[2];
u3(0.501362901232111,-1.30317928773733,-0.384779027073064) q[8];
cx q[8],q[2];
u1(1.42453554275500) q[2];
u3(-0.700277969595543,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.90300583735964,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.04833025131736,-3.95195437135377,1.76381212788931) q[2];
u3(2.07305197662979,2.34517341001651,-2.74632558032108) q[8];
u3(1.43582063620012,-3.03612081914252,0.132022878798631) q[3];
u3(1.58250322190734,-3.14572412462279,-0.913943628926619) q[12];
cx q[12],q[3];
u1(1.49577280010053) q[3];
u3(-0.534811652389124,0.0,0.0) q[12];
cx q[3],q[12];
u3(3.32309997204218,0.0,0.0) q[12];
cx q[12],q[3];
u3(0.820144269005247,0.778196677439078,-3.25196608212357) q[3];
u3(2.80383449994587,0.802862220064978,3.66003552928680) q[12];
u3(2.16842958979287,0.390712612811539,-2.61937611689994) q[5];
u3(2.11048504492038,0.660611520686774,-5.23307439375746) q[1];
cx q[1],q[5];
u1(2.38952279811658) q[5];
u3(-2.90685286009826,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.416959436725707,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.93595590704127,-2.54701737459572,1.37307516339473) q[5];
u3(1.27647066694064,-2.35107635056553,2.19814732315527) q[1];
u3(1.30274668010984,-1.08642806397157,1.66533526645316) q[6];
u3(1.69742640578431,-1.15182051450100,-1.65098445625540) q[0];
cx q[0],q[6];
u1(2.79825295403682) q[6];
u3(-1.70978189739008,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.573898346133049,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.46739581779166,-4.84898365275597,0.787206878056737) q[6];
u3(1.88208338692430,-0.444388063043687,-5.61180310728841) q[0];
u3(1.50742587795634,2.22196840706653,-3.40753097585256) q[8];
u3(0.910438511646071,2.60680145925195,-2.32528934600101) q[2];
cx q[2],q[8];
u1(2.97360447928343) q[8];
u3(-2.49110853572362,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.91750301281470,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.65076202209334,1.37341972857449,0.718986800398436) q[8];
u3(2.02558174726887,1.84151374905835,0.369789718404618) q[2];
u3(2.11767631111815,-3.60506691242885,1.02495142334549) q[7];
u3(0.988711903311698,-0.352587702611634,3.55532177460440) q[11];
cx q[11],q[7];
u1(0.670803054685154) q[7];
u3(-1.18605478012022,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.0544703658699890,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.25602056103313,-4.07010172163548,0.667798250192554) q[7];
u3(1.90717287969056,1.00680441986794,0.547786752616040) q[11];
u3(1.82284319630947,0.315105487363637,-1.40237196192151) q[4];
u3(0.601066787002379,-4.97336581256212,1.29034851700245) q[10];
cx q[10],q[4];
u1(0.955332280965803) q[4];
u3(-0.755542146724541,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.22675921767232,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.21551868650835,0.460062249418429,2.02419372567813) q[4];
u3(1.73282930077326,1.74952593575814,3.64535570080814) q[10];
u3(1.39367173206332,-1.78969660117943,-0.886224660154337) q[12];
u3(1.89266081727121,-2.81940604777761,0.236071870603482) q[9];
cx q[9],q[12];
u1(0.920902555728012) q[12];
u3(-3.63767634324008,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.90296030530956,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.82357976956669,1.19522276322006,1.04996070583938) q[12];
u3(2.77919332802232,0.389491501831501,-3.76287511919691) q[9];
u3(0.649494983862997,0.00203784116510963,-2.39272268999398) q[7];
u3(1.57308417924666,2.99692623938442,-2.71021514814631) q[0];
cx q[0],q[7];
u1(2.45282060571566) q[7];
u3(-1.83770652857419,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.124596044789934,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.78196069230153,-2.21681738606305,2.35842677746651) q[7];
u3(1.59196729783444,0.712945100254869,-1.41934485697678) q[0];
u3(1.51897835255147,-0.135754216336497,2.78130326030083) q[3];
u3(1.34183063104068,-1.30670470689770,-1.66163998007180) q[5];
cx q[5],q[3];
u1(1.35773463272886) q[3];
u3(-0.681504840602858,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.92249516568510,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.07083079840651,3.16273105725876,-0.605101696356582) q[3];
u3(0.864766449914551,0.0130778842532506,3.00174948676324) q[5];
u3(0.825061709682867,1.97802229885851,0.0314606868178822) q[4];
u3(1.50829023250757,0.486831294525041,-4.69453386363739) q[6];
cx q[6],q[4];
u1(3.41460455523285) q[4];
u3(-0.805731737686207,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.71875544997979,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.44795535145153,1.28912561784020,-0.377581902415854) q[4];
u3(0.800348752825580,4.66461023256577,0.425632909456466) q[6];
u3(0.487109454277189,1.51721265175719,-2.56211902117032) q[8];
u3(1.32883927635421,1.18459118457698,-4.12496461974800) q[11];
cx q[11],q[8];
u1(2.43383096049624) q[8];
u3(-1.85467994419671,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.49333792203526,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.10304154188402,0.970745950166509,-3.40927066979313) q[8];
u3(0.259801881849098,2.33177471676949,-2.74110808106144) q[11];
u3(1.44275322045010,-2.48566067751990,0.569340113986006) q[1];
u3(1.37555395227926,-2.76645694605471,-0.132672401480935) q[10];
cx q[10],q[1];
u1(0.598273027712275) q[1];
u3(-1.33754668725571,0.0,0.0) q[10];
cx q[1],q[10];
u3(-0.347934027480329,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.902030267753811,3.58317756633590,-0.782660165656359) q[1];
u3(0.561248799911696,0.759367154888162,-5.07833187368994) q[10];
u3(1.14076127161828,1.27996842233218,1.10727546796482) q[7];
u3(1.73908396079133,-1.50490357313669,-0.922728381630403) q[4];
cx q[4],q[7];
u1(0.340052948727644) q[7];
u3(-1.50611804941048,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.36553933089499,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.606805236565307,-2.48814493181836,-0.886222654621599) q[7];
u3(0.287510706920937,0.771881684295336,5.07300673890240) q[4];
u3(1.71073327415724,0.426885590491110,0.709075123362600) q[2];
u3(2.21803965444847,-1.14337674505937,-1.77639420867818) q[11];
cx q[11],q[2];
u1(0.848185385893479) q[2];
u3(-1.23984538538042,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.89927248245200,0.0,0.0) q[11];
cx q[11],q[2];
u3(2.82012175634921,1.53437632104241,0.369393620952723) q[2];
u3(1.41960166941561,-1.32051896495119,-1.33667525350537) q[11];
u3(1.80259691251136,-1.15903465769345,-0.836954411735473) q[9];
u3(1.70292918393068,-2.68379989873981,0.505260398153962) q[6];
cx q[6],q[9];
u1(2.25437537240425) q[9];
u3(-1.55405292975584,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.30864844820512,0.0,0.0) q[6];
cx q[6],q[9];
u3(2.72553742229208,3.19431798420345,0.664966332348350) q[9];
u3(1.24779410021575,-3.86562353034650,1.01028280150457) q[6];
u3(2.79794960843876,0.418589769031331,1.30457342722984) q[1];
u3(1.80175011271483,-1.27152807528798,-2.25139623457614) q[5];
cx q[5],q[1];
u1(1.84682136459653) q[1];
u3(-3.19448860359028,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.932937537139484,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.46160184166645,-1.31955388697594,2.75400820052357) q[1];
u3(1.13298122663996,-2.05516257138140,0.301251738319118) q[5];
u3(0.624835078736803,1.06516929298668,1.53416567616722) q[10];
u3(1.33985301455317,-0.526134535147732,-3.01960341486787) q[3];
cx q[3],q[10];
u1(0.0558868177671166) q[10];
u3(-0.797476540329170,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.26410234842445,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.44058566764812,-1.28126195608519,-0.614653436777337) q[10];
u3(1.07886416638181,-2.39222394815232,3.46484019971618) q[3];
u3(0.378851198448513,-2.31411026838155,1.86660991412986) q[0];
u3(1.23965221164123,-2.79590351871400,0.643100892061049) q[8];
cx q[8],q[0];
u1(3.28975978122115) q[0];
u3(-1.04722055612125,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.65302581253492,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.961748821238384,-0.505708774797359,1.37420846223202) q[0];
u3(0.847960176029162,-3.90168136951823,1.17365152736844) q[8];
u3(1.10325597673361,-0.122227970014913,1.12808609885286) q[6];
u3(1.20132951723786,-1.63708928143556,-1.17619036437956) q[7];
cx q[7],q[6];
u1(-0.281568775215819) q[6];
u3(1.29909969945244,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.64954959421847,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.501573029968932,1.03105676655073,1.72960633412531) q[6];
u3(1.27978082631691,2.55504064219892,-2.81003407248783) q[7];
u3(1.99036143858899,1.07913728433935,0.232621333545710) q[12];
u3(1.08760808671167,0.947026995114301,-5.25532934450468) q[9];
cx q[9],q[12];
u1(1.82249180646488) q[12];
u3(-2.91729819036957,0.0,0.0) q[9];
cx q[12],q[9];
u3(0.946028009934549,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.17500269765008,-0.0731246775529454,-2.04993904065947) q[12];
u3(0.556669505401539,1.27174705349345,-2.67547314270027) q[9];
u3(1.44297089424591,-0.950764392823100,-1.88709380595358) q[8];
u3(1.81694038258202,1.70445031961129,-4.03093561125636) q[2];
cx q[2],q[8];
u1(1.45403414754991) q[8];
u3(-3.25197499719696,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.28006100342500,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.40642060565947,2.55002072722285,-3.43629944624317) q[8];
u3(1.46976376114537,4.04255346377032,-1.41393078536110) q[2];
u3(1.89438918475807,1.44754151302669,-3.69336519330866) q[3];
u3(2.50823018199209,2.93982634182952,-2.57888244431099) q[11];
cx q[11],q[3];
u1(1.13097531204076) q[3];
u3(-0.572123508804369,0.0,0.0) q[11];
cx q[3],q[11];
u3(-0.0533551500508544,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.25249616215622,1.81165401657462,-2.23229269939286) q[3];
u3(0.896680081665276,-0.255981097619919,3.35181402281904) q[11];
u3(1.10319348076051,-0.998297435345757,-1.18457548035842) q[10];
u3(1.33276997780740,-3.77751983400990,-0.352557737128600) q[0];
cx q[0],q[10];
u1(1.77828981494496) q[10];
u3(0.110402602685066,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.663835706325651,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.872511858300298,-2.69521924652739,0.574444186601685) q[10];
u3(1.22454885623344,0.714092556091511,3.32122279263662) q[0];
u3(1.77525162713225,-2.12062373834431,-0.330406204385255) q[5];
u3(1.65232218276584,-4.56009756281568,-1.20850131671690) q[4];
cx q[4],q[5];
u1(1.52468209867956) q[5];
u3(-1.03836501060205,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.36196678004852,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.87543711262437,3.50398628727038,-1.54237682626213) q[5];
u3(1.89085603174287,1.86553463267671,0.271230732186930) q[4];
u3(1.79301187995767,2.64720736744455,-2.98037638124901) q[7];
u3(0.744897882917691,2.75011222508814,-2.99833373884402) q[0];
cx q[0],q[7];
u1(1.12755558973414) q[7];
u3(-0.766275159710456,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.61676994496836,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.99672960769963,0.630475288859280,-0.282957769311074) q[7];
u3(1.22933561295053,3.81247536643630,-0.209200291752858) q[0];
u3(1.59757426744936,0.568315087973341,-1.66654463571815) q[1];
u3(1.19317175352327,1.15738339792746,-4.44870712894240) q[2];
cx q[2],q[1];
u1(1.13506229005742) q[1];
u3(-0.927958378779953,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.0571912425310990,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.03361734365239,-2.46002288518837,3.52659593208018) q[1];
u3(2.86825616447617,-1.37296755648879,-2.50671851427411) q[2];
u3(2.43597534500369,-4.10161437707595,1.44942857407616) q[3];
u3(1.13560167069071,1.91631026549472,0.192373690647608) q[4];
cx q[4],q[3];
u1(1.58817133816682) q[3];
u3(-0.469954906409733,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.189132195099698,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.47106968819883,-2.67617359406595,2.79257822931660) q[3];
u3(2.32164608572008,-1.02759063152823,4.52499778204662) q[4];
u3(1.82045969796116,-0.732049144823606,1.79730951683748) q[6];
u3(1.55016359478390,-2.18298846890143,-2.26519679575689) q[10];
cx q[10],q[6];
u1(1.44882035360113) q[6];
u3(-3.53175854137741,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.42454408901503,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.518330812588449,0.0934867477249639,0.312714698316807) q[6];
u3(0.866849238950963,4.64790350961348,0.888526970611257) q[10];
u3(2.41777103995163,0.0350539447310329,0.831604618510259) q[5];
u3(1.59029832637485,-2.51706111992759,-1.05832351250169) q[12];
cx q[12],q[5];
u1(2.80433054455809) q[5];
u3(-1.93361586267507,0.0,0.0) q[12];
cx q[5],q[12];
u3(0.546980996505202,0.0,0.0) q[12];
cx q[12],q[5];
u3(0.109353194359719,3.99508192918240,-0.183253631821543) q[5];
u3(1.53242311991203,-4.59413588247317,-1.43642050248873) q[12];
u3(1.56184381440303,1.92863679159945,-2.87403406980111) q[9];
u3(0.813059244017183,1.92891259007865,-2.71006210383359) q[8];
cx q[8],q[9];
u1(0.622578319058092) q[9];
u3(-3.12842093672578,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.75707767265496,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.08855724324803,-0.515277075796686,0.530793310827430) q[9];
u3(2.57978880057461,-0.280361232061819,-2.59755053703953) q[8];
u3(0.827131647355072,-2.30258728874931,1.84985395327463) q[6];
u3(0.194499023956860,-1.73498286689068,-0.00539210865499240) q[12];
cx q[12],q[6];
u1(0.0904456006930316) q[6];
u3(-1.84335272648867,0.0,0.0) q[12];
cx q[6],q[12];
u3(2.77930492817740,0.0,0.0) q[12];
cx q[12],q[6];
u3(0.733213043592009,-0.283284086014426,0.983935555747121) q[6];
u3(0.481966436869591,-0.481051203085530,-2.68056976829334) q[12];
u3(1.91571304503136,2.08404050969618,-3.34767307610165) q[10];
u3(1.24981474210518,2.95216671354815,-2.96110310074131) q[2];
cx q[2],q[10];
u1(-0.0431718835402570) q[10];
u3(-0.479872397110327,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.85439177194433,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.39935377389539,0.624021741568736,-2.65216294260652) q[10];
u3(2.17543343056281,1.84643476992378,4.28499739770043) q[2];
u3(0.855437507651943,2.44437697892699,-1.12202963390430) q[0];
u3(0.175895145289743,-2.78080340983703,1.29944535284598) q[11];
cx q[11],q[0];
u1(2.58526985044181) q[0];
u3(-2.77910667615440,0.0,0.0) q[11];
cx q[0],q[11];
u3(1.24170596095701,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.67818948672812,0.740454568585731,-3.74878604611697) q[0];
u3(0.472238817669164,1.93418919424537,-3.41237125428794) q[11];
u3(2.48732664380554,-1.35930723306642,2.27666167837413) q[8];
u3(2.17710436404152,1.45570077348659,3.70225647374400) q[7];
cx q[7],q[8];
u1(-0.787606868655112) q[8];
u3(0.618496422483647,0.0,0.0) q[7];
cx q[8],q[7];
u3(3.66818337062009,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.50157826966515,1.69965166701089,0.315516504156022) q[8];
u3(0.632955386103372,3.89860584328014,-0.0146376172586158) q[7];
u3(0.926606724777040,-1.05190470302815,4.06137559043324) q[9];
u3(1.16989597786345,-0.981125655324787,0.626386124963601) q[5];
cx q[5],q[9];
u1(1.60247349244479) q[9];
u3(-0.0756224821091562,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.23765595967333,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.354294179501342,0.781187543283249,-0.515167672860470) q[9];
u3(2.74605682976405,-0.0644755333172109,-2.26991564872527) q[5];
u3(2.41312568263544,-1.61262402161717,0.560061343583792) q[1];
u3(2.75925695441751,-3.36828361698297,-1.67250734425575) q[3];
cx q[3],q[1];
u1(1.32796672900436) q[1];
u3(-3.11271224844299,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.23434655192658,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.62755342395345,-1.13672228787946,4.26443968032071) q[1];
u3(0.667986883218991,-5.42130465356781,-0.805812476372800) q[3];
u3(2.52226199272031,2.70517036121705,0.291327506244548) q[12];
u3(2.70140243035166,2.48251405509331,-3.13342420270473) q[6];
cx q[6],q[12];
u1(1.94216259742163) q[12];
u3(0.653515810761736,0.0,0.0) q[6];
cx q[12],q[6];
u3(1.06246115116135,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.55835721667539,-1.08176041618022,4.22563773117835) q[12];
u3(1.25153974291221,-2.64297574574149,-1.27323428691311) q[6];
u3(1.49608998760074,2.07303923579297,0.0515254332319068) q[11];
u3(0.293743229138493,0.623098523071897,-4.56725009417947) q[5];
cx q[5],q[11];
u1(1.15233396119607) q[11];
u3(-0.517687192895899,0.0,0.0) q[5];
cx q[11],q[5];
u3(2.85793516133813,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.10286731642709,1.28406855486776,-1.04719014067571) q[11];
u3(1.49033160841597,1.14950719733477,0.196950926534819) q[5];
u3(1.21545159783976,1.23178594930713,0.120911109458576) q[4];
u3(0.953298189649677,-0.311789703413395,-1.64953061213716) q[9];
cx q[9],q[4];
u1(2.05390958540242) q[4];
u3(-2.90901920400571,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.53100009001010,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.15049811780638,-0.132767354578483,4.21809059429953) q[4];
u3(1.67971965182634,-0.303211458534732,-0.816227121599629) q[9];
u3(1.50091513288114,1.93695774469073,0.297704726684500) q[7];
u3(2.69324177472355,-0.171878935596657,-4.38602242429869) q[3];
cx q[3],q[7];
u1(2.93110029721591) q[7];
u3(-2.44216847394082,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.30080867021385,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.12858716340023,3.83329227496722,-0.989004784990396) q[7];
u3(0.499851909271389,-3.26881578632999,0.160426429889340) q[3];
u3(1.59106708596160,2.25581093272696,-1.96405182010889) q[0];
u3(0.132599583277429,-3.22896454805912,1.76510493787984) q[10];
cx q[10],q[0];
u1(-0.181485727144501) q[0];
u3(-1.16750275903144,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.81343190802608,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.45280522733663,1.63933576483913,-1.71053259591558) q[0];
u3(1.68125026602333,3.92648860408130,-1.64121938632823) q[10];
u3(2.40843636615233,1.75695389714330,-1.14388563686413) q[2];
u3(2.88530947785664,-0.800638308575764,-4.32843708120641) q[8];
cx q[8],q[2];
u1(1.52165566502243) q[2];
u3(-0.991333901367211,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.72479412444853,0.0,0.0) q[8];
cx q[8],q[2];
u3(0.472405967591789,-0.835472749120511,-0.457725631695987) q[2];
u3(1.09086021628357,-3.84778296610202,1.18042116238119) q[8];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
