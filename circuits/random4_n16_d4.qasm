OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(1.63007209755483,1.09386630330633,1.28652398403986) q[12];
u3(0.410181810491915,-1.14876254915688,-3.02760604164088) q[13];
cx q[13],q[12];
u1(1.50815942603288) q[12];
u3(-2.75253695806123,0.0,0.0) q[13];
cx q[12],q[13];
u3(0.919193959328542,0.0,0.0) q[13];
cx q[13],q[12];
u3(2.85167439063361,2.21270524775370,-0.568291894997227) q[12];
u3(1.91701347445321,0.542037106413681,4.69970058014164) q[13];
u3(0.976193979133303,2.90097761016389,-1.17784625305523) q[4];
u3(1.25137499538636,1.57988597018147,-1.31958872779801) q[10];
cx q[10],q[4];
u1(3.78354801252453) q[4];
u3(-4.39543676817652,0.0,0.0) q[10];
cx q[4],q[10];
u3(-0.384859197947357,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.10362356074993,-3.08676816803116,1.70627059153017) q[4];
u3(2.41131523971670,4.31054369866910,1.27386620785219) q[10];
u3(1.79771612675041,2.73058810994408,-0.113368357408841) q[8];
u3(2.56359218101571,1.61232815851113,-1.78014627524473) q[0];
cx q[0],q[8];
u1(2.02228719760039) q[8];
u3(-2.60875807604264,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.277998792403254,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.28137699864247,1.29322749606653,0.192974474202182) q[8];
u3(0.912251187896830,1.89833939633091,-3.75341759784860) q[0];
u3(1.78223912045754,1.42419968726305,-0.676658763459870) q[1];
u3(1.36503689131245,0.701297613477297,-4.09723832597149) q[5];
cx q[5],q[1];
u1(0.959466451030928) q[1];
u3(-0.218351215666595,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.91382523287635,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.90034334547014,0.258304040515415,-3.31284653207406) q[1];
u3(1.62041879074103,-4.79690546588392,0.671840355854586) q[5];
u3(1.65533430022482,0.833438019017453,-0.0375535204017648) q[11];
u3(1.25037270147901,-0.303895898123797,-4.47761571407316) q[9];
cx q[9],q[11];
u1(1.31131253235396) q[11];
u3(-0.981750152632115,0.0,0.0) q[9];
cx q[11],q[9];
u3(2.60810783558566,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.22772874798868,1.46220306758002,2.34162566156959) q[11];
u3(2.11938535041753,4.37814360813175,-1.56869704937669) q[9];
u3(1.77238856762708,0.247254949546506,-2.07639529322372) q[7];
u3(2.89930643920535,5.55964721607070,0.659895506958490) q[3];
cx q[3],q[7];
u1(3.01613983393453) q[7];
u3(-1.85065619485256,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.903910255033622,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.57379579366646,2.26607053777712,-1.11084809624516) q[7];
u3(1.42437035395585,-2.31970187301342,-0.0225872901077353) q[3];
u3(0.354672960115139,-0.489427637056676,0.0990406121459106) q[15];
u3(0.897476452320735,-2.52386470177844,0.809486258431187) q[2];
cx q[2],q[15];
u1(1.52138365335726) q[15];
u3(-2.96507272069520,0.0,0.0) q[2];
cx q[15],q[2];
u3(0.651929126748265,0.0,0.0) q[2];
cx q[2],q[15];
u3(1.08893065442491,0.713477019489872,-0.837387944874672) q[15];
u3(0.822314191736388,1.41766020102512,4.01227896574416) q[2];
u3(1.17636031434344,2.42394587437862,-3.01768133784611) q[6];
u3(1.00452144973403,3.45493016166908,-2.61082126990353) q[14];
cx q[14],q[6];
u1(0.536381503057702) q[6];
u3(-1.85106287932557,0.0,0.0) q[14];
cx q[6],q[14];
u3(2.77099206410835,0.0,0.0) q[14];
cx q[14],q[6];
u3(1.28040803743056,3.15207401742470,0.500913764264165) q[6];
u3(1.09032602087021,1.92805975930650,-1.62138834095054) q[14];
u3(2.28053964169795,1.45240155497901,-4.04453317201633) q[9];
u3(1.03879669241584,-1.11850816177858,3.23946825710546) q[2];
cx q[2],q[9];
u1(1.52841047019090) q[9];
u3(-3.42213177712574,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.05994702768211,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.34319492028917,-4.18720004976604,1.24583443021656) q[9];
u3(2.73862932966197,-4.13415179876787,-0.110833064065381) q[2];
u3(1.89971377042258,1.76863729290226,-3.57173139888142) q[12];
u3(0.993813570942452,-2.46132056782478,2.99908223768438) q[7];
cx q[7],q[12];
u1(1.44742632128341) q[12];
u3(-0.0558915587479054,0.0,0.0) q[7];
cx q[12],q[7];
u3(1.01670009998179,0.0,0.0) q[7];
cx q[7],q[12];
u3(2.02675405246878,-3.85790812718185,1.37844684226735) q[12];
u3(0.439887304286493,-3.56224726675995,-0.648848496542324) q[7];
u3(0.949976579834070,0.641468498483467,-2.78139593271771) q[5];
u3(2.07324392878824,2.17303836209922,-3.83445211781778) q[6];
cx q[6],q[5];
u1(4.22227298196898) q[5];
u3(-3.44328793683887,0.0,0.0) q[6];
cx q[5],q[6];
u3(-0.258005868069862,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.38401977749160,-0.0689334831016784,-0.663350470994436) q[5];
u3(2.18297545186226,-4.78292260117816,-1.32916267636848) q[6];
u3(1.39263565857541,-1.91387308646534,1.27402624335349) q[3];
u3(0.295684271180271,1.22520387886837,-2.58519584451982) q[11];
cx q[11],q[3];
u1(3.54518702299407) q[3];
u3(-4.38790653769669,0.0,0.0) q[11];
cx q[3],q[11];
u3(-0.266495167017257,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.29782495195273,1.04889751273955,-3.49447842320918) q[3];
u3(1.66520521135731,0.133757610426780,-1.80197675046836) q[11];
u3(2.86903631240000,0.132950992933447,0.0351466105765395) q[15];
u3(0.994315469376790,-0.0428974383076972,-4.49725292806955) q[10];
cx q[10],q[15];
u1(3.30911807776835) q[15];
u3(-1.16213425253254,0.0,0.0) q[10];
cx q[15],q[10];
u3(1.80391339130546,0.0,0.0) q[10];
cx q[10],q[15];
u3(1.25573107441217,-2.82904102479641,1.81395823291992) q[15];
u3(1.94100719134358,-0.439817994175397,3.71158160797054) q[10];
u3(1.21947200693320,-1.33945266546556,1.53115457053764) q[8];
u3(0.194038993910674,2.46636745924511,-3.29663134961307) q[0];
cx q[0],q[8];
u1(0.274859088557278) q[8];
u3(-0.641771394690854,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.75805140073191,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.68288321156954,-3.04819726110351,0.904335120174314) q[8];
u3(1.49394987075505,-0.823762865064921,-0.937236793162228) q[0];
u3(2.16993013520974,-2.96188663980347,-0.144674092670671) q[13];
u3(2.97955211604575,-1.88468402733627,-1.27266632611180) q[14];
cx q[14],q[13];
u1(1.56903525135840) q[13];
u3(-3.18460774946056,0.0,0.0) q[14];
cx q[13],q[14];
u3(1.29946284952807,0.0,0.0) q[14];
cx q[14],q[13];
u3(0.651931800494611,-1.64890211644449,-1.37287861849799) q[13];
u3(1.22375268397890,1.44924108955318,1.82725354920114) q[14];
u3(2.02541023732383,1.05159727236592,-3.73533440030245) q[1];
u3(1.47197296984606,2.82966951495884,-3.03337835620197) q[4];
cx q[4],q[1];
u1(1.74816822873594) q[1];
u3(-0.548002010172462,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.23576035431773,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.03257143959432,-2.38881298968728,0.0251651116749418) q[1];
u3(2.03361422710212,-2.11824475672817,3.25941219336018) q[4];
u3(1.22533999699588,1.26113189676330,0.645866594588852) q[15];
u3(0.295933899437141,0.180862675017203,-4.26979051325142) q[8];
cx q[8],q[15];
u1(1.43892847842179) q[15];
u3(-3.06710928530923,0.0,0.0) q[8];
cx q[15],q[8];
u3(0.484160088012417,0.0,0.0) q[8];
cx q[8],q[15];
u3(1.15417755918260,0.268475658291157,0.377812040353340) q[15];
u3(1.62267911670752,2.32939207999202,-1.75921963721493) q[8];
u3(2.09035472447295,0.482717123151389,-3.15456370416713) q[2];
u3(2.22049294465474,-1.66888484525552,4.40157314509722) q[4];
cx q[4],q[2];
u1(1.28962766944037) q[2];
u3(-3.10445336820667,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.47595840571605,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.76209934053747,-1.99091821797149,-0.162768668791763) q[2];
u3(2.48054100201867,2.98292476708435,-0.668000366424627) q[4];
u3(0.765893489037314,-0.0403602834278423,-1.87282553756960) q[11];
u3(1.32120904009902,2.52359424443034,-3.22947244059862) q[0];
cx q[0],q[11];
u1(0.232999729377952) q[11];
u3(-0.866930986020076,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.94602357433391,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.74572153109459,-1.13691136098154,3.31989668858760) q[11];
u3(0.468617083410779,1.38385366449805,2.90885881164343) q[0];
u3(0.972443729674468,-2.46921430417431,-0.400349103844946) q[12];
u3(0.962668907208180,-3.53894081087991,0.710615510304388) q[1];
cx q[1],q[12];
u1(2.57633204060709) q[12];
u3(-2.98196209032552,0.0,0.0) q[1];
cx q[12],q[1];
u3(1.72717138766676,0.0,0.0) q[1];
cx q[1],q[12];
u3(1.76153427723686,2.66899755242925,0.935648107586647) q[12];
u3(2.40832491401453,1.71482414649431,-1.35366439632269) q[1];
u3(1.48079220898336,0.142494891965352,1.78669093051338) q[7];
u3(1.82724587686551,-2.05171153965857,-1.37642444681834) q[14];
cx q[14],q[7];
u1(2.48773771644606) q[7];
u3(-2.17292942301219,0.0,0.0) q[14];
cx q[7],q[14];
u3(1.33304806232057,0.0,0.0) q[14];
cx q[14],q[7];
u3(2.36630508230626,3.60758691470387,-0.875967105072088) q[7];
u3(0.705837181823276,0.868082148206671,3.84960651205416) q[14];
u3(1.07491570323731,1.01069721911834,-0.973811138995608) q[3];
u3(0.971495350716230,-0.843293166343925,-1.04106425880122) q[5];
cx q[5],q[3];
u1(3.24573844046192) q[3];
u3(-1.39088774041934,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.32858305676425,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.07379145281264,3.14340539905134,-1.48327023056413) q[3];
u3(0.701789859392092,3.77810101254191,2.13482843023695) q[5];
u3(2.73556051534306,-1.51596734003447,1.51742257675641) q[6];
u3(2.82863954404841,-2.09956480618740,0.968814717182062) q[13];
cx q[13],q[6];
u1(0.770088242136247) q[6];
u3(-1.35280194896467,0.0,0.0) q[13];
cx q[6],q[13];
u3(2.91706900236224,0.0,0.0) q[13];
cx q[13],q[6];
u3(1.01636844483789,-2.90595862806125,1.92924289523511) q[6];
u3(2.76858307089844,2.85921845843994,-0.799661209542660) q[13];
u3(2.00386273347303,-1.77510622080091,-0.0380887001066146) q[9];
u3(1.23295841141838,-3.96588199674486,-0.478584103645686) q[10];
cx q[10],q[9];
u1(1.85514908312910) q[9];
u3(-2.66603921288196,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.941358034521319,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.93358595350602,1.63133478932917,-3.49414443275515) q[9];
u3(0.998722112586237,2.70081337713242,-1.53368043014200) q[10];
u3(0.653734219833148,2.79578569575267,-1.54824320512798) q[13];
u3(1.88184672048463,1.39425221446104,-2.37473693243484) q[3];
cx q[3],q[13];
u1(2.06901214136706) q[13];
u3(-2.75310719271882,0.0,0.0) q[3];
cx q[13],q[3];
u3(0.705908745189778,0.0,0.0) q[3];
cx q[3],q[13];
u3(2.46122034286859,-0.117304952308381,-1.85627450867178) q[13];
u3(1.70821895173009,-4.16909789625866,1.99716856442566) q[3];
u3(1.28951259690838,0.769333506026561,-2.23796097956327) q[7];
u3(1.04832721022560,-3.93541537490507,2.26792207782865) q[0];
cx q[0],q[7];
u1(-1.34615300512997) q[7];
u3(0.622220676641370,0.0,0.0) q[0];
cx q[7],q[0];
u3(3.63075161363708,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.55655080763159,1.36605964132299,0.325336618936587) q[7];
u3(2.07353486671687,-2.29125160807449,3.12350776409937) q[0];
u3(1.85413876604575,2.36819616086856,-2.97838130983602) q[9];
u3(1.79144351064320,3.06585141882373,-3.20431862749289) q[10];
cx q[10],q[9];
u1(-0.0369787392524288) q[9];
u3(-1.26101092212908,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.777372984002982,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.01058728642292,0.998430071119041,-1.16277493895528) q[9];
u3(0.990435018476309,-4.09251038332658,-0.858154296655790) q[10];
u3(1.25473000764831,0.786382895828864,1.17167691369849) q[2];
u3(1.43002861735998,-1.41520912498867,-2.25113945485669) q[6];
cx q[6],q[2];
u1(1.46164978323027) q[2];
u3(-3.54315072073702,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.01727643692600,0.0,0.0) q[6];
cx q[6],q[2];
u3(3.07804041872037,0.664785991607506,-0.779705794652947) q[2];
u3(0.431014315373040,-4.27094587386250,-0.989034407910216) q[6];
u3(1.76468792235722,-1.08322236518813,0.445283489674196) q[8];
u3(1.78200214730597,-2.35688468799179,-1.13878424135477) q[15];
cx q[15],q[8];
u1(0.0994106005058861) q[8];
u3(-1.08494425581589,0.0,0.0) q[15];
cx q[8],q[15];
u3(2.59433271796996,0.0,0.0) q[15];
cx q[15],q[8];
u3(0.862831665241094,-0.519532099371203,-0.356743354929914) q[8];
u3(1.96205584561886,0.339014086970631,-0.409254027090743) q[15];
u3(0.970285066034389,2.10741351284037,-3.24136699387715) q[5];
u3(0.529720797408161,1.22307394558456,-2.38656556929833) q[1];
cx q[1],q[5];
u1(2.85947878662474) q[5];
u3(-2.24690860414315,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.28811574716654,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.91794169665016,1.15682813002624,-4.34669900974923) q[5];
u3(1.23426182926576,-2.25773757412435,1.40498220227714) q[1];
u3(2.74572786472500,-1.29597627856006,1.86537250399681) q[12];
u3(1.73693900429087,1.87177109615471,4.12819098775082) q[4];
cx q[4],q[12];
u1(1.33735485938124) q[12];
u3(-0.319857874919560,0.0,0.0) q[4];
cx q[12],q[4];
u3(2.23474545274249,0.0,0.0) q[4];
cx q[4],q[12];
u3(0.201802256947935,-3.18211024537135,1.65618239371045) q[12];
u3(0.551220086763147,-0.380989098075927,5.51778993254214) q[4];
u3(0.779584616941814,-1.63885854010116,2.03048077943628) q[11];
u3(0.558165757770405,-0.363279202265872,-1.70372162013812) q[14];
cx q[14],q[11];
u1(1.40465445029402) q[11];
u3(-0.899638303791963,0.0,0.0) q[14];
cx q[11],q[14];
u3(2.57042054961513,0.0,0.0) q[14];
cx q[14],q[11];
u3(1.74987832631497,2.34011727782905,-2.70802903159065) q[11];
u3(1.94034182271968,-0.312803617606566,-3.91276226178750) q[14];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
