OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.05056118491397,1.25135600332963,1.80954222057402) q[4];
u3(1.57772253573982,-1.58618797376286,-0.692848350100978) q[5];
cx q[5],q[4];
u1(0.642110882454130) q[4];
u3(-1.08942433981271,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.135543765252338,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.92346963964873,1.72782992273004,-0.340481524279833) q[4];
u3(1.53303465024499,-1.96567385011297,-1.96579941593294) q[5];
u3(1.32542906376742,0.377980486329102,2.33453558163445) q[9];
u3(1.23581981866409,-0.930430684850171,-0.879696853324735) q[6];
cx q[6],q[9];
u1(2.17001901779687) q[9];
u3(-2.66404037301022,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.158428972691132,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.29329098956157,-4.64670531146215,1.36450425170386) q[9];
u3(2.01151617080736,1.79610145669044,-2.39472902721284) q[6];
u3(2.47713054324385,4.45795691087832,-1.42859848243887) q[2];
u3(0.971381764538238,1.45155103862418,-0.0176971783500576) q[3];
cx q[3],q[2];
u1(2.50577188783722) q[2];
u3(-1.91914843899551,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.61831111900169,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.57170150907916,-1.91224777436228,0.994323261560205) q[2];
u3(1.05352674464643,-4.35493256649009,-1.20583878806727) q[3];
u3(2.10572526343988,2.14680258558528,-1.91427650778183) q[1];
u3(1.76967034085097,1.30993821383111,-2.89988413868974) q[0];
cx q[0],q[1];
u1(0.612137876193992) q[1];
u3(-1.23567956075137,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.176573031250824,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.40723433390710,-3.70488292796461,0.426096394793015) q[1];
u3(2.49400730434676,1.45772210300398,1.83695644889627) q[0];
u3(1.47960123639367,2.74815581492630,-2.57174040694707) q[7];
u3(2.30431039822317,-3.08838071279335,3.18920620045176) q[8];
cx q[8],q[7];
u1(2.38551736965670) q[7];
u3(-3.18207528131387,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.08832111594146,0.0,0.0) q[8];
cx q[8],q[7];
u3(3.05196831167230,-1.45602806097772,-0.224321857553514) q[7];
u3(2.41897283406136,-3.35848403769395,1.05791692545186) q[8];
u3(0.711104566025819,2.01012950669397,-2.60031183704333) q[6];
u3(1.07928122600598,-3.30775050986642,2.53296873615370) q[9];
cx q[9],q[6];
u1(1.53485050828650) q[6];
u3(-0.480062881323903,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.04509370266762,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.967910614830505,2.67669057648071,-2.84330054306640) q[6];
u3(2.57428890154561,-1.53135537138363,-1.88278285060514) q[9];
u3(2.45615958082238,2.61629435466685,-0.564865410536488) q[2];
u3(2.80110066708225,5.57590713869395,0.691208461268997) q[0];
cx q[0],q[2];
u1(1.03319939166050) q[2];
u3(-0.0938177444930717,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.37593704873043,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.910927807624380,-1.64291065006059,4.61240485245918) q[2];
u3(0.319456868930609,0.606284590236492,-5.09163360118717) q[0];
u3(2.13881106975252,0.537185530202546,1.45734874208495) q[1];
u3(1.88185318487630,-2.45064693289191,-1.66468179643753) q[4];
cx q[4],q[1];
u1(-0.908961872041698) q[1];
u3(0.455282768986644,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.14278956043586,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.66007396158897,1.97808540448204,-2.79692914591589) q[1];
u3(0.730151632765246,-2.74097272507578,-1.97603535745546) q[4];
u3(0.638755270682974,2.22782570403113,-2.07404263445418) q[3];
u3(0.477858283146716,0.623043666767678,-2.79832288748579) q[5];
cx q[5],q[3];
u1(1.36946677675297) q[3];
u3(-3.38244246718479,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.10996100020869,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.42360262083180,-1.68334191170633,2.42679455301183) q[3];
u3(2.61041837359863,-1.48156311119177,-3.79223143275543) q[5];
u3(1.07840475946752,2.11248977225811,-2.95510921506794) q[8];
u3(1.40576511128619,-2.52251673849086,3.25402553222548) q[7];
cx q[7],q[8];
u1(1.61445308407546) q[8];
u3(-2.15561517402113,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.95056614695956,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.86005288538186,2.22024837061938,-2.68395876281422) q[8];
u3(1.60023749172872,-0.295014188529402,3.09849805795277) q[7];
u3(1.76611290337770,2.44237635050822,-1.24148647851794) q[0];
u3(2.31982969550319,2.28479549096918,-0.550633457873919) q[4];
cx q[4],q[0];
u1(0.0628687505423930) q[0];
u3(-0.868276448340076,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.63257637199580,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.60195252606369,2.68095778198348,0.0136920242031686) q[0];
u3(1.34654649101395,5.86749095143268,0.320303592925165) q[4];
u3(0.873951416926698,3.48102893362421,-1.63368660956757) q[1];
u3(1.95187761790634,1.26473985590283,-2.63106940469346) q[2];
cx q[2],q[1];
u1(2.34618858052928) q[1];
u3(-2.65010896918436,0.0,0.0) q[2];
cx q[1],q[2];
u3(-1.24205921998564,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.689859322621980,3.18533022614136,-2.45935751613715) q[1];
u3(1.43063043831639,1.34823613786183,-0.0327140714034266) q[2];
u3(2.60581473455892,-0.823149939287611,0.372740687287581) q[5];
u3(1.40606876902877,-3.17638077982463,-1.39670133656250) q[8];
cx q[8],q[5];
u1(2.17964503634134) q[5];
u3(-1.75318749556701,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.259144011737567,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.50538947662032,-1.90594401598106,2.48736448141672) q[5];
u3(1.57457833706182,-1.13151206315506,-4.32426531041681) q[8];
u3(0.269980542426850,-0.507718247287279,0.725476168358591) q[7];
u3(0.348594941752679,2.76379543816461,-2.97817947255167) q[6];
cx q[6],q[7];
u1(3.15591478010066) q[7];
u3(-1.52423705658127,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.54118029257881,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.33664949482807,2.71994950914561,-3.25382137472842) q[7];
u3(1.50521704783543,4.38305019432363,-0.609804316272740) q[6];
u3(1.73634410902862,1.01547593869059,-0.205914985306062) q[3];
u3(1.21300548646564,1.01225791115288,-4.32962676991165) q[9];
cx q[9],q[3];
u1(0.183476830536605) q[3];
u3(-1.34276706546691,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.22128215643856,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.555619284641214,-0.144013852633715,-1.88603411977541) q[3];
u3(2.89919865593157,-4.76900740418492,-0.212716715640776) q[9];
u3(2.41014841552260,-4.23670409630197,1.45993198479491) q[3];
u3(0.659186938372213,0.939677627383656,1.01393925214596) q[9];
cx q[9],q[3];
u1(2.90076147274217) q[3];
u3(-2.46907498248502,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.55066763062209,0.0,0.0) q[9];
cx q[9],q[3];
u3(3.07327776004144,-0.638971532131700,2.81654747173107) q[3];
u3(0.682564033454223,1.75274810749675,-3.37556519336287) q[9];
u3(0.396964197358250,0.593632265451538,-0.425410848399817) q[7];
u3(0.213467242372219,-0.447953239405028,-1.81491184717739) q[4];
cx q[4],q[7];
u1(0.933104182548665) q[7];
u3(-1.38692969567229,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.94290127959274,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.856628569326142,-3.08816055793468,2.62812292358808) q[7];
u3(0.903339921126386,-0.836604673550257,2.09600395651210) q[4];
u3(1.46900115985207,-1.61003128415221,1.16137937713178) q[6];
u3(2.07599893733548,-3.77301768264230,-0.214081265929436) q[2];
cx q[2],q[6];
u1(2.34831683155779) q[6];
u3(-2.86916135565280,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.25312647945996,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.46648463739596,-0.179675569519373,-0.0960385090455715) q[6];
u3(2.26322033177049,-2.14420658212430,-2.76491988803833) q[2];
u3(2.31196057493056,0.319696235744406,0.552312537339827) q[0];
u3(1.72444076896318,-2.01225005866124,-1.75949774849506) q[1];
cx q[1],q[0];
u1(0.197249926540674) q[0];
u3(-1.53479577666310,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.473379523967696,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.10730442700758,-1.80706031061538,-1.54156201504606) q[0];
u3(0.959661256885752,2.31385879927763,0.440966727668534) q[1];
u3(1.00669044013670,2.64614682905754,-2.40534386961077) q[5];
u3(1.01255880852386,1.44990452351576,-2.23581983300009) q[8];
cx q[8],q[5];
u1(1.15317972530003) q[5];
u3(-0.358858539751014,0.0,0.0) q[8];
cx q[5],q[8];
u3(3.12425821814382,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.28846871438647,-2.32175016108277,1.08270320101534) q[5];
u3(2.62355755344544,-1.84231131886145,-0.607266224795607) q[8];
u3(2.52262841238306,-3.06776169200947,2.59355987880875) q[0];
u3(0.755200803085782,-0.930387112759528,3.04809689872833) q[7];
cx q[7],q[0];
u1(0.400207785843126) q[0];
u3(-0.638412823118463,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.81311567656097,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.738970434004724,-1.19203812757822,-0.348621149609460) q[0];
u3(1.84060593388404,-0.728953275940911,-1.27069062664533) q[7];
u3(2.81986984745087,3.92079093660074,-1.74118253012289) q[1];
u3(1.27260974960718,-0.304766154393737,1.77023618716627) q[6];
cx q[6],q[1];
u1(0.673411924082701) q[1];
u3(-3.41550697952981,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.67779774302456,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.02065537864205,0.676493030893129,-3.89364924103014) q[1];
u3(1.18652986032638,-2.59917898925759,0.736116801453091) q[6];
u3(1.65336918795225,0.0567686229873104,-1.06617329431724) q[2];
u3(2.54600513628841,-3.11533505626392,2.37484063539250) q[9];
cx q[9],q[2];
u1(1.48139606455632) q[2];
u3(-0.974114452668550,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.63130779010266,0.0,0.0) q[9];
cx q[9],q[2];
u3(0.341117794713852,0.747489557929951,-1.39106433605183) q[2];
u3(1.48021924793390,3.64260410814036,0.292258451081645) q[9];
u3(0.223978731141997,3.31698323365487,-1.75548200924086) q[4];
u3(1.52554624140067,2.45501331296165,-1.24705826205845) q[8];
cx q[8],q[4];
u1(2.17765995862164) q[4];
u3(-2.85936896637771,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.909046076618092,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.58921054224456,-1.09373245865322,0.250007224550354) q[4];
u3(2.78410214745576,-1.35494500265587,1.59045450166775) q[8];
u3(0.659759506744302,-2.47442654060247,2.62656092328696) q[3];
u3(0.976583600536542,-3.27823506290948,2.08580729154440) q[5];
cx q[5],q[3];
u1(2.03735929200319) q[3];
u3(-2.56834632606512,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.0294457378036475,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.34930686836333,-2.10518870556818,1.35179675298027) q[3];
u3(1.90481431945769,-4.95453442372263,1.10371156806841) q[5];
u3(0.425667887208830,-0.994853645890123,2.18066646305222) q[1];
u3(0.709038778846782,-3.24994143659662,1.76486022301853) q[6];
cx q[6],q[1];
u1(3.77811064407648) q[1];
u3(-4.44834505150242,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.894963808907575,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.55973286644018,1.54070641487716,0.358387368560688) q[1];
u3(1.92983058522757,4.71759410125360,0.0156012919302433) q[6];
u3(1.92618309492353,-1.06038323669880,2.47093056936929) q[9];
u3(1.29701060076498,-1.32415596139854,-1.45224066102196) q[5];
cx q[5],q[9];
u1(-0.149126600431503) q[9];
u3(-1.59304384780574,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.24498390320491,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.948102313876563,0.183657416237950,1.02981613359854) q[9];
u3(0.878670392894638,-0.0607578105297288,-3.47483153955505) q[5];
u3(2.30387152633214,-0.995495466583322,0.000189957306781063) q[0];
u3(0.640752617344643,-5.06076938262712,0.0220350378274077) q[2];
cx q[2],q[0];
u1(0.932518533385410) q[0];
u3(-1.35716404864598,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.04879703897211,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.928317616277004,-1.98058126044477,0.148076611176637) q[0];
u3(1.72740731569097,2.44475567287984,-0.890325772465664) q[2];
u3(2.15062238859859,-2.82933242086554,2.93826264239713) q[7];
u3(0.793356394381022,-0.471403076450043,2.32788293521533) q[3];
cx q[3],q[7];
u1(3.54734241374732) q[7];
u3(-0.742717285728249,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.83524683063048,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.07192998273915,-2.59556186773501,3.39594926092427) q[7];
u3(1.25169942969876,0.649517638161521,-0.0990600009765155) q[3];
u3(1.59238701861029,0.451272499838212,2.63744350204976) q[4];
u3(2.19330020521357,-2.40751302335198,-2.49000935919948) q[8];
cx q[8],q[4];
u1(2.34070166581251) q[4];
u3(-1.69870609086610,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.136421572407797,0.0,0.0) q[8];
cx q[8],q[4];
u3(0.330455947464760,-4.09075585466264,0.0109609110347546) q[4];
u3(1.81158360556745,2.99672462408687,-0.378725828335984) q[8];
u3(2.48926007888021,0.769905631186466,-1.54252604893484) q[1];
u3(1.28435628205327,-4.03981344798941,1.26202338038502) q[8];
cx q[8],q[1];
u1(0.131602355976593) q[1];
u3(-1.37885185773190,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.50992719569413,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.85684224545998,-0.111192191996151,-2.79058093791041) q[1];
u3(1.97838398444641,1.29486025479227,2.36493896401934) q[8];
u3(1.41628702379812,1.26286393630267,-3.65472503265965) q[7];
u3(2.63339025715121,2.44118472184299,-3.00822161709556) q[3];
cx q[3],q[7];
u1(3.12630810909935) q[7];
u3(-1.37213288792798,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.64049822936489,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.10900427724296,2.58626037767468,-0.456591127301750) q[7];
u3(1.71052643353944,-0.116768876873258,2.51502453518312) q[3];
u3(1.15733231625871,1.12225280768546,-1.94180408475256) q[6];
u3(0.101561032083791,0.00336888006326674,-1.38734669390729) q[2];
cx q[2],q[6];
u1(2.81257320657604) q[6];
u3(-2.32402970870990,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.06946418483169,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.63287354945967,-1.94917199055070,2.31584104279179) q[6];
u3(1.16820104468350,-3.12203878439551,-0.567148656143889) q[2];
u3(0.611078587212303,1.05339022823063,-1.18033545097297) q[5];
u3(0.861404912365307,-0.276232040824958,-0.214585880115229) q[4];
cx q[4],q[5];
u1(0.714124242742778) q[5];
u3(-1.33744716141195,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.70179363629047,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.73084139806036,0.0842670176511359,-1.59324507559389) q[5];
u3(1.07505747649673,-3.08514609319090,2.38318210436769) q[4];
u3(1.39731381972484,2.57410678214908,-2.63395641361316) q[9];
u3(2.19292886091802,-2.87969618301071,2.70406001346450) q[0];
cx q[0],q[9];
u1(1.95920405243837) q[9];
u3(-2.75505380417379,0.0,0.0) q[0];
cx q[9],q[0];
u3(0.348507504579946,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.54182660410195,-0.389358278398668,1.56758712618332) q[9];
u3(2.22322153963937,-2.94771159100733,0.296326702813768) q[0];
u3(1.46909175998829,0.516243923675190,-2.56628794882947) q[3];
u3(1.81953748584630,-2.61215227862153,2.70061036422572) q[8];
cx q[8],q[3];
u1(1.24439279345167) q[3];
u3(-0.714141980960911,0.0,0.0) q[8];
cx q[3],q[8];
u3(-0.151419716641898,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.24494223384825,3.10650203658048,-0.427076873530627) q[3];
u3(1.89959344241180,0.491824927194773,3.74008804772127) q[8];
u3(2.09360337954646,-0.690047734835369,-0.136354623652213) q[2];
u3(1.60842835629450,-2.78173644980491,-0.928019672220782) q[0];
cx q[0],q[2];
u1(0.793689761452567) q[2];
u3(-3.13713607748129,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.64208454824380,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.56780470158928,-1.01355267519077,2.07362682701566) q[2];
u3(1.39473993567035,4.22202906403692,-0.306719347681936) q[0];
u3(1.17706960443933,0.866713151408011,1.78981387660769) q[1];
u3(1.05770973101606,-1.88501925568616,-0.587148629450214) q[6];
cx q[6],q[1];
u1(1.38654342256085) q[1];
u3(-0.758698190854460,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.573747204888688,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.875848508326160,-1.06128309589807,1.38554815083328) q[1];
u3(2.11346644130940,-0.689658730158131,4.93984284345690) q[6];
u3(2.87935235529619,2.46027303385726,-2.47566834037796) q[9];
u3(1.59816388610035,-3.42275385644869,2.51887381937954) q[5];
cx q[5],q[9];
u1(0.195473877593617) q[9];
u3(-0.902905489457235,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.77701785564391,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.79796550929259,-1.19955258820529,-0.0756328550092989) q[9];
u3(1.98984802919523,2.97570511738931,-3.13301436277954) q[5];
u3(0.427873336179560,0.254855079325330,-0.521811317382763) q[4];
u3(1.09276504931623,-3.50450472055622,1.46187329531266) q[7];
cx q[7],q[4];
u1(-1.02659190978705) q[4];
u3(0.301715520695649,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.85040473888128,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.52683857369642,-2.63331581702193,0.462581100824650) q[4];
u3(1.61425491547317,-3.38810235202762,-2.20765607073323) q[7];
u3(1.15886823763602,0.465347013954404,-1.81023675931509) q[7];
u3(1.87880631769106,-3.98707160570035,0.537582827488564) q[2];
cx q[2],q[7];
u1(0.924386870880093) q[7];
u3(-1.49524436204816,0.0,0.0) q[2];
cx q[7],q[2];
u3(-0.695661746134202,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.13886032917205,3.00933107256228,-2.70446957249663) q[7];
u3(0.358833065776682,-1.55754853160246,1.31553989998936) q[2];
u3(0.605608132885706,1.33156477405010,-0.390716693719713) q[1];
u3(1.18459031995786,0.914060079345013,-2.38383046487571) q[8];
cx q[8],q[1];
u1(1.43667854037259) q[1];
u3(-3.58408385140579,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.31888387489207,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.54728075486219,-1.87938820606242,2.52955134825493) q[1];
u3(1.34418241370567,-3.81506232788481,1.66020168575303) q[8];
u3(2.59741880651147,1.45585296254229,-2.37222013459199) q[0];
u3(2.33166475236903,-0.666728210527790,-5.34140706047546) q[5];
cx q[5],q[0];
u1(1.14002779989484) q[0];
u3(-0.757949746958026,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.33708966915001,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.31927373043908,1.27635136128443,0.909619565543912) q[0];
u3(1.12314705388273,4.07240681097479,1.95295781603321) q[5];
u3(0.967410275227488,2.40596989625265,-0.686938283393939) q[4];
u3(1.08438276033532,0.789063515101261,-2.60191020432402) q[6];
cx q[6],q[4];
u1(2.25610181757981) q[4];
u3(0.464595882099167,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.70373359766081,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.21716724563556,-1.06487672484428,1.72387686559280) q[4];
u3(0.872590093164946,0.365970276125266,4.98251762823338) q[6];
u3(1.37734327322208,-0.966602645241332,2.64108797227734) q[3];
u3(1.25765478306595,-2.15498833750043,-2.06080113208587) q[9];
cx q[9],q[3];
u1(0.338552022092177) q[3];
u3(-1.34105400259806,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.23334461028414,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.04172161479149,1.31939748663795,0.00234475691953329) q[3];
u3(1.32886583530916,-3.81225690905815,-0.0464791702637124) q[9];
u3(0.432123237524603,-0.746030281825440,1.15803601653033) q[4];
u3(0.227303642353204,0.492987551110612,-1.99942676429293) q[8];
cx q[8],q[4];
u1(0.196090385211262) q[4];
u3(-1.79813986724876,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.62167096778791,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.14745666898483,-0.683474127235526,0.570772511702878) q[4];
u3(1.16970184320948,-1.52724178790135,2.87986815588676) q[8];
u3(2.68814381051563,3.62774186488939,-0.684132046747620) q[3];
u3(1.75070654781930,2.21968095724058,-2.58294896134399) q[7];
cx q[7],q[3];
u1(2.90096780057561) q[3];
u3(-2.09704261739249,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.318874424921706,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.321281270344358,-1.65053324708141,1.78035052837814) q[3];
u3(0.141547007015154,0.0192267696628443,-6.17791917609320) q[7];
u3(1.33204665301788,1.59272945234743,-3.36241743888657) q[6];
u3(2.64355145752390,1.66616717397176,-3.41753847098529) q[9];
cx q[9],q[6];
u1(1.60703730344477) q[6];
u3(-2.91704966074303,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.859214669313956,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.03315841535494,-3.24615930040010,2.16017308224874) q[6];
u3(0.243747179547141,-0.0764084426672143,-0.312849537853929) q[9];
u3(1.65134327167442,-0.712405305501017,-1.22178670846949) q[2];
u3(1.00788966356267,-3.80876383140686,1.01826518903682) q[1];
cx q[1],q[2];
u1(1.88758773218244) q[2];
u3(-2.39766922154193,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.643253870926163,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.23863086842696,0.298344900920746,-1.43862727201393) q[2];
u3(2.20547199081721,1.11175844787701,1.61194451050373) q[1];
u3(1.55618165121185,-1.83230386938800,-0.309967439711341) q[5];
u3(1.26235140141027,-4.43116860681558,-1.17066332279376) q[0];
cx q[0],q[5];
u1(1.80371743837157) q[5];
u3(0.306562892281427,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.57130517182679,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.24925700175921,2.50170841506035,-0.667199797353088) q[5];
u3(2.24007535577205,-1.70666306029263,1.04206130670783) q[0];
u3(2.82994417181377,1.55795941509172,-4.64979313992356) q[4];
u3(0.824904895131024,0.0397511639306385,2.09557790687945) q[0];
cx q[0],q[4];
u1(1.25037741681262) q[4];
u3(-0.176607667751539,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.49231059031790,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.98682719375615,-1.96368664527887,-1.70540622718623) q[4];
u3(1.26683876053926,-0.601465613984877,5.06309881086643) q[0];
u3(1.11930943954141,0.278613372361584,-1.00178101805682) q[3];
u3(1.78464627077455,-5.44111380151098,0.771002359774163) q[8];
cx q[8],q[3];
u1(2.51802752301175) q[3];
u3(-2.77547472617565,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.35180909677948,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.475325289303326,0.955576773080019,3.53835542583120) q[3];
u3(2.68898514739248,-2.49710706071034,1.07811387556142) q[8];
u3(2.16559438042116,2.48837834835218,-2.52097917561039) q[9];
u3(0.914090265196650,-2.72109960258997,3.40430723331787) q[2];
cx q[2],q[9];
u1(3.13263959860464) q[9];
u3(-0.796648381157626,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.59104724739305,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.30707675363234,-1.07297787056875,0.580578263387848) q[9];
u3(0.575474767447324,-0.260323141718979,-3.18975842757868) q[2];
u3(2.28955029785064,2.73788049112520,-1.35651696623292) q[5];
u3(1.96431973785978,1.88565948418824,-0.477385795989051) q[7];
cx q[7],q[5];
u1(1.48368738656168) q[5];
u3(0.0350144882528285,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.286927826027520,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.997327892201732,-2.51635970581825,2.25514697222126) q[5];
u3(1.98136109529361,1.69109862956728,-3.71225343435682) q[7];
u3(2.49052916216983,-0.0237328177028981,-2.92470110246044) q[6];
u3(2.58518728662247,4.16984801385583,0.0380401671157364) q[1];
cx q[1],q[6];
u1(-0.540555696862791) q[6];
u3(0.349792806609092,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.91322235871049,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.99754699301085,1.96784365420009,-2.40403354178893) q[6];
u3(0.486081718227822,-1.86076397169560,1.28055092299059) q[1];
u3(1.12990532519266,3.26753838505955,-0.721315003163361) q[1];
u3(2.05157977715319,2.03423078241949,-1.63135545952862) q[4];
cx q[4],q[1];
u1(1.72640582683808) q[1];
u3(-3.49114122298701,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.14919701407863,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.21137924589078,1.33448078820998,0.495268354214106) q[1];
u3(0.737494585042278,-0.342299200961556,-2.03879704391859) q[4];
u3(1.26845314176661,0.700892172153812,1.86336308267582) q[9];
u3(1.51493951210241,-1.00116505717607,-2.78336798841534) q[3];
cx q[3],q[9];
u1(1.74662495734535) q[9];
u3(0.237249238858853,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.585892155748951,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.13370960286507,1.89090704685763,-2.84481104434460) q[9];
u3(1.79953405334533,-2.03295862139182,-4.15959199265513) q[3];
u3(1.04354517617675,-1.04560483404824,-1.03850341492144) q[6];
u3(2.02830467721445,-4.22080004770533,0.946883262145678) q[0];
cx q[0],q[6];
u1(1.11479760003420) q[6];
u3(-3.34283730864298,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.95031819489776,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.94082974583424,-2.12741266214308,0.461309869096054) q[6];
u3(2.27836129341560,-2.69539887124098,2.43547019352861) q[0];
u3(0.962767470835007,1.60331817660019,-1.85036303065314) q[5];
u3(0.183839590013085,0.0610884630228004,-1.09095838963741) q[2];
cx q[2],q[5];
u1(2.75266390394661) q[5];
u3(-1.48978917161233,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.71021715394102,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.74067118021796,-4.39731277551720,1.88187198243625) q[5];
u3(0.914173300364521,0.190550746490001,1.98955763313426) q[2];
u3(2.32657672974846,2.53917592672270,-2.53907439895805) q[7];
u3(1.75837104909083,2.09764140918008,-2.42356596684706) q[8];
cx q[8],q[7];
u1(-0.0679295040123158) q[7];
u3(-1.18089387779429,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.43551165145164,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.299481938447104,-1.26007822507566,1.16106824809321) q[7];
u3(1.76335548250146,2.00528454149696,-1.62401603559200) q[8];
u3(2.86274135338115,2.56338966630527,-0.394421032924650) q[4];
u3(2.54317979755860,0.139512173685147,-3.44891001646029) q[6];
cx q[6],q[4];
u1(0.113906234261210) q[4];
u3(-0.683855474294309,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.11390017586673,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.61768162654922,-1.20949735793315,1.20209655847575) q[4];
u3(2.08896223459816,-0.0381360571208527,-6.03101587219724) q[6];
u3(1.24109054095727,1.42749705456455,-2.19969535684360) q[0];
u3(1.31349227787554,1.63503383660228,-2.78452204013001) q[7];
cx q[7],q[0];
u1(2.10548799892979) q[0];
u3(-3.05212436558362,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.44285119253859,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.33772872682153,1.50282440486985,-2.22202827925067) q[0];
u3(1.47228426318596,1.87537804664150,-3.20219585570134) q[7];
u3(0.919311766901849,-0.308594445108063,1.08296349896441) q[2];
u3(0.208009016964506,1.78405690579720,-2.78518537272836) q[3];
cx q[3],q[2];
u1(1.33460713563123) q[2];
u3(-0.155931028643549,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.12953082832420,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.81186252069003,3.43038564675326,-2.29394441697452) q[2];
u3(2.35467615540848,-0.0474390591561238,-5.57936373076084) q[3];
u3(1.15569572849844,-2.37038243215845,2.08176748056297) q[9];
u3(0.516607282990266,0.923475339441135,-3.51276699905616) q[5];
cx q[5],q[9];
u1(1.03833125322856) q[9];
u3(0.00586267884869818,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.96305107775781,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.975137665863726,1.30487224881783,-2.25450038475140) q[9];
u3(1.27787778846022,-4.22272808677116,0.786628532082562) q[5];
u3(0.542571177161378,1.34452601036857,0.764578896605576) q[1];
u3(1.31546473463052,0.814208283283763,-3.06807943841348) q[8];
cx q[8],q[1];
u1(3.51310904353515) q[1];
u3(-4.05369000027582,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.486001771465947,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.68539371577953,-2.16678358214107,1.77743239763707) q[1];
u3(1.85905627776365,1.42095485131533,-0.702627981789992) q[8];
u3(1.75939860366548,-2.90279703391240,2.35096023400599) q[9];
u3(2.93575571802966,-1.23316855115694,2.82683862511314) q[1];
cx q[1],q[9];
u1(2.63703546198250) q[9];
u3(-2.01725404893515,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.799985750562948,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.356870525626979,2.16604475821126,-0.238523503152776) q[9];
u3(2.53791357878768,-0.373431654322756,-1.99639502105996) q[1];
u3(1.50965530075963,0.184530160800789,-1.94824605552862) q[8];
u3(0.311555052022979,-4.62090620192195,0.634121175299828) q[6];
cx q[6],q[8];
u1(1.31358914430533) q[8];
u3(-0.484753861045104,0.0,0.0) q[6];
cx q[8],q[6];
u3(2.47382828858710,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.02840486095772,1.67729025756167,-3.36910420949718) q[8];
u3(1.62410509443270,3.07698230570355,-1.00209318090227) q[6];
u3(0.600971037408466,-3.19720585897864,2.64477406028196) q[0];
u3(1.44076116097220,0.563003179735838,-1.58556357371989) q[3];
cx q[3],q[0];
u1(-0.0810164557362292) q[0];
u3(-2.18643262586554,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.51573696495629,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.28017656038224,0.363640986592290,-1.89292495312401) q[0];
u3(0.516789319602716,1.63517303609234,1.45446204547793) q[3];
u3(0.424924736925445,-1.57176220943977,1.43555713741806) q[7];
u3(0.542173038208786,-2.83354161263969,1.12182022065974) q[4];
cx q[4],q[7];
u1(2.98198539210986) q[7];
u3(-1.10448218930514,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.27123394596222,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.77639871989501,-0.713550833887443,-0.0572725273810694) q[7];
u3(1.97085869180129,-0.566133195251314,-0.759517567486458) q[4];
u3(2.32450500747645,2.66697099967845,-3.10515346946509) q[2];
u3(0.566713600560578,2.82536874133973,-2.01957170332837) q[5];
cx q[5],q[2];
u1(1.04078404125198) q[2];
u3(-3.62349340358463,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.50160531721656,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.901465267278631,-3.18351227449955,0.801244728904242) q[2];
u3(1.73439620701983,-0.136234590932209,-5.82541591460883) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
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
