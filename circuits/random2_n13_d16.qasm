OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.53965736748969,1.39842587971411,-1.90305344534313) q[1];
u3(1.88920221374429,-5.08260967767266,0.513129779008034) q[11];
cx q[11],q[1];
u1(0.809095624455915) q[1];
u3(-1.50878151077898,0.0,0.0) q[11];
cx q[1],q[11];
u3(2.84686834838683,0.0,0.0) q[11];
cx q[11],q[1];
u3(2.08034788671684,-3.88725811455846,2.33291988536865) q[1];
u3(2.48350979629963,-1.24218939809537,-4.54587798337936) q[11];
u3(2.61523936351674,-1.18781311162281,4.14360160215136) q[10];
u3(0.568684562098995,-0.803568024707636,2.48309485803589) q[0];
cx q[0],q[10];
u1(-0.422488726480266) q[10];
u3(0.988246832232616,0.0,0.0) q[0];
cx q[10],q[0];
u3(3.98423643296781,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.07213831226742,-0.622934374220813,-0.913583925002825) q[10];
u3(1.18638928587119,4.32250196885549,0.303615001241097) q[0];
u3(0.380173607329188,-2.05738689604652,2.95693572004988) q[4];
u3(1.29758338531477,1.53100251399440,-2.55279193971365) q[9];
cx q[9],q[4];
u1(1.53098631864968) q[4];
u3(-2.90092506704007,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.760711502980985,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.46309113723534,0.593053128844276,-4.82980050743530) q[4];
u3(1.52210817595013,2.38115841163819,-0.401789427105151) q[9];
u3(1.74829867007183,2.41832013524757,-3.15487626262714) q[5];
u3(0.748895076914223,-2.22517211602482,2.74829452068951) q[2];
cx q[2],q[5];
u1(1.81613449442320) q[5];
u3(-3.18880806601337,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.603760529728548,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.46000139259068,-3.36785125577320,2.43724090618652) q[5];
u3(1.00554949003570,0.255440122897689,1.09738277767828) q[2];
u3(1.30729832234732,0.451047113630313,2.39589475201503) q[3];
u3(1.86359490479285,-2.58199905319468,-2.00980513441575) q[6];
cx q[6],q[3];
u1(1.57770396138980) q[3];
u3(-0.319202020813707,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.33683788511602,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.12107096489576,-1.18170447186508,-0.532939568072493) q[3];
u3(2.36753724120355,-0.395035120276463,-4.11653178804894) q[6];
u3(1.79184858891117,-0.0766732057506467,1.60520148191637) q[8];
u3(1.04029955548869,-0.316854198557969,-1.26985070293176) q[12];
cx q[12],q[8];
u1(3.35497877228388) q[8];
u3(-0.985428520132385,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.53712563183835,0.0,0.0) q[12];
cx q[12],q[8];
u3(1.19802835870085,-0.829510875405434,3.14187182919785) q[8];
u3(1.30052194780903,-1.41875364607684,2.35995757551888) q[12];
u3(0.910044254022830,3.11546468427601,-2.95583516538579) q[0];
u3(1.06965664830170,1.13860585150172,-1.38684364751364) q[9];
cx q[9],q[0];
u1(2.49978956420767) q[0];
u3(-2.02948211015633,0.0,0.0) q[9];
cx q[0],q[9];
u3(-0.204896158456201,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.97286872696578,2.11320417773139,-3.27407969130211) q[0];
u3(0.576976505934791,-3.41318523135688,-2.61761120468041) q[9];
u3(2.22240714369689,0.202520042153991,-1.15554159237444) q[11];
u3(1.78593724384972,-3.80014145133973,0.997444957646093) q[2];
cx q[2],q[11];
u1(-0.535252651692629) q[11];
u3(0.0876346264892849,0.0,0.0) q[2];
cx q[11],q[2];
u3(4.01007974499184,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.47876694086933,-0.0853875852797472,1.25466950694278) q[11];
u3(2.58857349014204,-4.69178098187072,0.212751375900083) q[2];
u3(1.71817855120408,-3.60401557181606,2.61657203919871) q[4];
u3(0.297234469514840,-2.09978306898151,2.96801582784371) q[8];
cx q[8],q[4];
u1(1.21348543164654) q[4];
u3(-3.25623739760796,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.31916077223115,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.71751372780962,2.75746871805613,0.146837406078920) q[4];
u3(0.619375537289606,0.190592435523417,4.17301922242019) q[8];
u3(2.59527769021072,-1.82886991659922,0.651901614223509) q[5];
u3(2.48488313712729,2.15955735981670,3.42906477384540) q[10];
cx q[10],q[5];
u1(2.21573331174339) q[5];
u3(-2.63341508994814,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.24803085566920,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.91030388707824,1.24825324278965,-1.15305526714813) q[5];
u3(1.64004967885174,-3.36987586934319,-1.62420055029397) q[10];
u3(0.800663519982291,1.45351146054002,-3.19707834942498) q[1];
u3(2.35209498438848,2.80559584308066,-2.99364122564269) q[6];
cx q[6],q[1];
u1(1.26122306705174) q[1];
u3(-3.20954352672724,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.20919444447294,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.96502477679259,1.02793082155292,-0.887131559876297) q[1];
u3(2.61668366805696,-0.100098439003514,-3.38903303141986) q[6];
u3(1.19292950215986,0.385649154126781,-1.79735670428198) q[3];
u3(1.09745159635932,-3.38608182817584,1.14252801928978) q[7];
cx q[7],q[3];
u1(1.43765050723939) q[3];
u3(-0.566575820139669,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.99289024054462,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.53055315130494,-1.27784774639155,2.03773781043227) q[3];
u3(0.982801879771895,0.0702875532560472,-0.932478953511043) q[7];
u3(1.44099958904451,-1.81909918762715,0.451660974686870) q[6];
u3(1.89010611973838,-2.22355508974714,0.929484465602053) q[10];
cx q[10],q[6];
u1(1.70144549226796) q[6];
u3(-2.80332600216569,0.0,0.0) q[10];
cx q[6],q[10];
u3(0.922605563559028,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.33138483177235,-0.956573256118135,2.08149919596186) q[6];
u3(2.51652500632688,-2.59223714298712,3.38210620508444) q[10];
u3(1.57703359828731,-0.815584818052547,2.17220000392558) q[7];
u3(2.15171324584358,-2.00015700213347,-0.920104908373842) q[11];
cx q[11],q[7];
u1(1.22910565002636) q[7];
u3(-3.36104400316408,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.19708147550663,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.30811054379777,-3.38970000439067,1.40327478965774) q[7];
u3(0.649994287580410,-0.840845221610070,5.42401969225795) q[11];
u3(1.22767861617671,0.241497631059144,2.31261578456575) q[9];
u3(1.29149877683659,-0.426992591535118,-1.51535617174157) q[12];
cx q[12],q[9];
u1(0.455937457394798) q[9];
u3(-0.944814876465442,0.0,0.0) q[12];
cx q[9],q[12];
u3(1.67910966255447,0.0,0.0) q[12];
cx q[12],q[9];
u3(1.19703202197969,2.08400972235035,-2.74286628088546) q[9];
u3(2.07792684068892,-2.68601215465170,3.40107288471907) q[12];
u3(1.93053167483049,2.33429361317296,-0.603015392227778) q[3];
u3(1.71540823872887,0.556931985741443,-3.56446236308322) q[2];
cx q[2],q[3];
u1(1.98765098907190) q[3];
u3(0.582799870683008,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.41162217465892,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.38964510733212,1.21802471286002,0.0384613844645146) q[3];
u3(2.73663379882722,2.86946921627234,3.02615230986251) q[2];
u3(1.31043451693375,1.26753840116862,-0.378606871450861) q[5];
u3(2.42744788592499,-1.12374763154016,-4.67653383252094) q[8];
cx q[8],q[5];
u1(4.49784202751635) q[5];
u3(-3.70426037280259,0.0,0.0) q[8];
cx q[5],q[8];
u3(-0.613377217231905,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.43126018334280,-2.19676188173711,-0.856754125779387) q[5];
u3(1.05312307848198,1.05009911886143,-3.48190223059565) q[8];
u3(0.193392141603180,-0.0529108945211537,-0.768301205878030) q[4];
u3(1.44272366820338,-3.70450719947876,0.913172580884219) q[0];
cx q[0],q[4];
u1(0.957936658620863) q[4];
u3(-1.36876168860732,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.307618382091015,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.467474935833016,0.836067060931354,2.90421041152748) q[4];
u3(1.60526430825925,-5.01151679730376,1.15624589942231) q[0];
u3(1.26354249064329,2.45732696557114,-0.163225070748909) q[7];
u3(0.573362236137781,0.285870997330232,-3.06046799786550) q[8];
cx q[8],q[7];
u1(1.88013063480940) q[7];
u3(-3.05791722435918,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.750601621065299,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.55744209743808,1.23892317349143,-1.32994601683916) q[7];
u3(0.0713678856473080,2.91511112436781,1.89671608125633) q[8];
u3(1.16187288053242,-2.25299880760215,0.842096286985036) q[10];
u3(1.42001740975014,-2.11383622003154,0.158330037700444) q[6];
cx q[6],q[10];
u1(-0.978826338796395) q[10];
u3(0.131963763338881,0.0,0.0) q[6];
cx q[10],q[6];
u3(3.84135518389912,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.15815562541663,-1.61531468387898,-0.127565871016905) q[10];
u3(0.935667461172466,-1.09772209527064,3.14161681310091) q[6];
u3(0.559620036257775,2.56964828707431,-2.95813998229036) q[1];
u3(1.05464622656139,-2.92854428737748,1.84623063974930) q[4];
cx q[4],q[1];
u1(1.70263819047631) q[1];
u3(-2.64692824799430,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.428860481048116,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.03828017502363,3.08772442161025,-2.09498853610771) q[1];
u3(1.61555304257658,-0.631123152583325,5.06703956621028) q[4];
u3(0.882559339643053,-0.354651803607700,1.91571061856835) q[0];
u3(0.771139641076954,-0.600890704504334,-0.975333168462742) q[12];
cx q[12],q[0];
u1(1.17400234023845) q[0];
u3(-2.93621913063879,0.0,0.0) q[12];
cx q[0],q[12];
u3(1.59799791500139,0.0,0.0) q[12];
cx q[12],q[0];
u3(1.41193123835433,-1.12028677675380,-1.69522827873631) q[0];
u3(2.61134536440189,2.40660001333117,1.02584210378568) q[12];
u3(3.06914394797564,2.20040344051785,0.358335267323444) q[5];
u3(1.94845597088305,0.898911502387753,-3.51388438632862) q[3];
cx q[3],q[5];
u1(1.77897385505799) q[5];
u3(0.206233233736574,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.23511442869891,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.44867681768287,-0.00937210481815975,-1.14378853169023) q[5];
u3(0.621685562853627,1.02580343269478,-4.68437874613120) q[3];
u3(1.87522258053963,-0.143410824414131,-2.59801118109901) q[2];
u3(1.45655683311524,-2.99278767320191,1.97995183306418) q[9];
cx q[9],q[2];
u1(3.50545163307363) q[2];
u3(-1.68678757112001,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.56233018436162,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.44985285147649,-1.64857574384330,-2.01833830139343) q[2];
u3(0.783944316483106,-0.421996263393464,5.10665830317672) q[9];
u3(1.10916993114259,0.859804901159514,-3.25488866653489) q[9];
u3(2.25486956223387,-4.16226516582858,2.04984203665861) q[10];
cx q[10],q[9];
u1(1.52256609997966) q[9];
u3(-3.57689729332562,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.60887798629566,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.70630112182171,-2.76222062997694,-0.290900311348102) q[9];
u3(2.62033145946458,-2.97432769029411,-2.61659907287357) q[10];
u3(1.28748486902856,2.48615859868102,-1.04483021857241) q[12];
u3(2.42039965067736,0.684008826768680,-1.99978887587101) q[8];
cx q[8],q[12];
u1(0.353436623282127) q[12];
u3(-0.747822780857239,0.0,0.0) q[8];
cx q[12],q[8];
u3(1.98719420966979,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.86221524754088,-0.606027440709988,-0.0745006928835832) q[12];
u3(2.07364184298808,0.434370495927637,3.11739684381932) q[8];
u3(0.496202809855575,-0.0849783297376979,0.294974172156572) q[0];
u3(1.38757123952506,-0.210415105917889,-1.17593462914016) q[7];
cx q[7],q[0];
u1(2.41229090112517) q[0];
u3(-2.98353758238205,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.20966301575213,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.703156814256459,2.27550991605156,-1.67989465977909) q[0];
u3(1.77112648615783,4.92474769069178,0.328984220841607) q[7];
u3(1.44263029371109,3.01148923817700,-3.02924956730460) q[1];
u3(1.67071789439419,-3.12131172211744,2.79341469868770) q[6];
cx q[6],q[1];
u1(1.56997404026435) q[1];
u3(-3.51170823057566,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.774506647953500,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.794837982986119,-0.246536754038126,-2.38079672505498) q[1];
u3(2.83206342162480,-2.65602247364568,0.924138290785046) q[6];
u3(2.43615570533999,2.04321398656014,-1.20442703839657) q[3];
u3(2.48105563892470,0.520821997847222,-5.15872407565987) q[5];
cx q[5],q[3];
u1(1.26041970799080) q[3];
u3(-3.25854167417664,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.02116044152159,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.712238663623939,-2.60176691773011,0.638832522982313) q[3];
u3(1.95636989094623,2.00355930971521,3.67591516001289) q[5];
u3(0.432613024425418,0.478713195912274,-1.51439646468408) q[4];
u3(0.207849070736195,-2.59669396921587,0.554757941130190) q[11];
cx q[11],q[4];
u1(2.21906047909808) q[4];
u3(-0.0690194322517808,0.0,0.0) q[11];
cx q[4],q[11];
u3(1.21178553340992,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.28244602463976,-1.65899206201726,-0.690824971667480) q[4];
u3(1.78895807600799,2.60117113233221,1.21529411500903) q[11];
u3(0.837413231604207,-1.88231364033417,-0.171229468658882) q[5];
u3(0.880226781905913,-3.29975093734190,0.628135240467975) q[1];
cx q[1],q[5];
u1(1.14841787799020) q[5];
u3(-0.810448800915481,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.153176628809561,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.61957328672471,-3.59731616281327,1.17235900284293) q[5];
u3(2.34524105610114,1.74442733224514,2.95139365007851) q[1];
u3(0.422596228029140,0.174202692949843,-1.37260220060006) q[9];
u3(1.20684993057523,1.09608246705373,-4.97051615540170) q[11];
cx q[11],q[9];
u1(2.09231829669338) q[9];
u3(-2.35837341841314,0.0,0.0) q[11];
cx q[9],q[11];
u3(0.106119762736913,0.0,0.0) q[11];
cx q[11],q[9];
u3(0.685805593643809,-2.74782575913794,-0.384789086069020) q[9];
u3(1.43337865062615,2.98724698780270,0.578062875757570) q[11];
u3(1.38030208007794,0.383062343921574,0.342830407838568) q[8];
u3(1.00972874834224,-0.0543582357246415,-3.69732763008936) q[12];
cx q[12],q[8];
u1(-1.37715603092823) q[8];
u3(0.644789283730980,0.0,0.0) q[12];
cx q[8],q[12];
u3(3.97056073573750,0.0,0.0) q[12];
cx q[12],q[8];
u3(1.25802486417937,0.458873390719442,-2.89125615224805) q[8];
u3(0.164651549360935,1.77270934813885,-3.04369421677872) q[12];
u3(1.51372675119486,-2.24340820869612,0.00279129877368756) q[4];
u3(1.41317965070353,-3.77555034440132,-0.630152903137584) q[10];
cx q[10],q[4];
u1(3.35590614547643) q[4];
u3(-1.50829140394659,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.01803412838406,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.69573724335548,-0.737002494653694,-1.25543090989368) q[4];
u3(0.814395492371488,4.77822738089317,0.109772334350529) q[10];
u3(1.67725835775928,0.0987653375122921,2.40934898335370) q[6];
u3(0.859354810872365,-3.11640790634311,-2.21742615039730) q[7];
cx q[7],q[6];
u1(-0.353774296462666) q[6];
u3(-1.58450315944951,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.897889919567182,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.91735327787940,4.11213468341204,-0.677489998546308) q[6];
u3(3.12260799220393,1.03387126508859,2.76336311388991) q[7];
u3(2.29416827506327,0.322397339135629,-1.20860459904192) q[2];
u3(1.35237339384881,0.200924926630819,-2.91263258175328) q[0];
cx q[0],q[2];
u1(-0.0483137464476686) q[2];
u3(-1.68750447063704,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.38652960008539,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.65138581690013,-0.860889404052303,4.78293836488397) q[2];
u3(2.94098712532116,-0.477554230075558,-3.42386687324359) q[0];
u3(1.40985421946000,-1.06294996471606,-1.65139487678179) q[2];
u3(1.57893157148714,1.48965628162969,-3.36025858065730) q[1];
cx q[1],q[2];
u1(2.51830021730520) q[2];
u3(-1.63030017113996,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.316977977999351,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.30866117761603,-0.132086205444476,-0.831711527420739) q[2];
u3(1.91981997722092,2.19970955310654,-3.43819621197822) q[1];
u3(2.60522240037115,0.252234058270248,-2.68244134481413) q[3];
u3(1.90828710029564,-3.49373192774090,2.49623035674619) q[10];
cx q[10],q[3];
u1(0.563618299300245) q[3];
u3(-0.880094863897180,0.0,0.0) q[10];
cx q[3],q[10];
u3(1.94345951091946,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.85302405142740,0.558065891502216,0.841380769587265) q[3];
u3(0.654547574304272,3.15541765084127,0.302352517427401) q[10];
u3(1.64579429477673,-1.28589973544251,0.601348227359004) q[8];
u3(0.674848368977628,-1.88337713895824,-0.868298838676988) q[9];
cx q[9],q[8];
u1(1.23055434033926) q[8];
u3(0.192876109637122,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.77086127830906,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.93995717572410,0.333624397011764,-3.19810780099673) q[8];
u3(0.652393131515116,2.67669930658048,-3.49263039178947) q[9];
u3(1.53426485222921,0.768974708490415,-3.75029310555220) q[12];
u3(2.42287506902525,-1.15250419178175,4.51480483503293) q[0];
cx q[0],q[12];
u1(-0.610584351716713) q[12];
u3(-1.79861180525153,0.0,0.0) q[0];
cx q[12],q[0];
u3(0.973441938074971,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.21031900528614,-1.74875655525433,1.49681413417574) q[12];
u3(1.99480227623047,1.74196077425866,1.90783180046823) q[0];
u3(1.80518345390189,-0.390270485024914,-0.729131431092042) q[7];
u3(1.43771602444677,-2.79205939979927,0.503775358752529) q[6];
cx q[6],q[7];
u1(2.90506673797572) q[7];
u3(-1.63663532900174,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.840479611976313,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.66236312206703,-2.07042961884272,4.15637564323321) q[7];
u3(2.65504722581101,0.503507076231394,4.69548912519858) q[6];
u3(1.35483523397897,-1.34008869238113,-0.571798508002172) q[4];
u3(0.869580780795293,-4.23142088034905,0.349625477967732) q[11];
cx q[11],q[4];
u1(2.22694484090339) q[4];
u3(-0.0517825435511110,0.0,0.0) q[11];
cx q[4],q[11];
u3(1.42704319289194,0.0,0.0) q[11];
cx q[11],q[4];
u3(2.32769371072423,1.08540047963113,-3.16312061306047) q[4];
u3(2.14996141577897,-4.67036885157120,-0.352942218786954) q[11];
u3(1.74731446191099,-3.03255556471784,0.473479072769939) q[5];
u3(2.44685833071633,-3.01003344428549,-0.554458261796499) q[0];
cx q[0],q[5];
u1(1.28516191935888) q[5];
u3(-0.636952379832113,0.0,0.0) q[0];
cx q[5],q[0];
u3(-0.286216906837345,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.44978925341859,-3.75568646324562,2.34504410019425) q[5];
u3(1.67398499548462,1.91639549817889,-3.90472555480085) q[0];
u3(1.23248353278316,1.92365460401630,-1.70457611322241) q[9];
u3(0.142583473281424,-1.80308127469283,1.18816082399847) q[6];
cx q[6],q[9];
u1(1.37855126397689) q[9];
u3(-1.08084851663095,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.571347225699408,0.0,0.0) q[6];
cx q[6],q[9];
u3(2.67044670152763,-4.39805965734025,1.17412365969673) q[9];
u3(1.99424500314377,1.97836856418483,2.69416038666200) q[6];
u3(2.35213985251651,-0.568682079567818,-0.345514200029840) q[10];
u3(0.969051436141819,-0.126331900622231,-5.19128686536812) q[1];
cx q[1],q[10];
u1(1.99045664509728) q[10];
u3(-1.76596400875435,0.0,0.0) q[1];
cx q[10],q[1];
u3(-0.448593135135600,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.29209795007207,0.907591479007220,0.428050323417448) q[10];
u3(0.482866208711696,0.476510733453400,1.71894150588032) q[1];
u3(1.21974040825036,-1.60779412061096,3.04670447992893) q[11];
u3(1.31384308707835,0.635767617662785,0.0716133012876061) q[12];
cx q[12],q[11];
u1(1.18166802085307) q[11];
u3(-3.56776367332319,0.0,0.0) q[12];
cx q[11],q[12];
u3(2.31975271188336,0.0,0.0) q[12];
cx q[12],q[11];
u3(1.79863782705971,-1.18307514535015,-1.99929236880410) q[11];
u3(0.505346521273508,-2.24324220883196,1.76584358427961) q[12];
u3(1.90661648723033,3.40417758627812,-0.914624272272646) q[7];
u3(2.45667912600858,1.70198349499537,-0.536281732816905) q[3];
cx q[3],q[7];
u1(2.20808959495215) q[7];
u3(-2.74552982465034,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.20900478221728,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.268753393423274,-1.86519552922991,2.56613878119952) q[7];
u3(3.00053693317204,4.71008691829029,1.26711474083873) q[3];
u3(2.85344904845208,4.20151982851377,-1.58852713978880) q[4];
u3(0.778572243346101,-0.314065351670896,1.20127208992334) q[8];
cx q[8],q[4];
u1(2.17809213109671) q[4];
u3(0.409287203719149,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.26875105538676,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.92062224799464,1.26677564920405,-1.41543201108443) q[4];
u3(1.38683460091592,1.96550863152583,2.68458896163441) q[8];
u3(1.86184194979866,-0.605365899488160,2.31406994910580) q[1];
u3(1.56940975933746,-2.09905426013043,-2.04103001372523) q[4];
cx q[4],q[1];
u1(-0.00443601581091113) q[1];
u3(-1.49819687868569,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.728941566823856,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.14172958929628,1.25461750241953,2.41556542498575) q[1];
u3(2.68417736765275,1.05275228784006,4.07665128008795) q[4];
u3(0.761095478440698,2.81066963051648,-1.53535962802280) q[0];
u3(1.08704006039269,1.44059884317213,-2.14288795845731) q[10];
cx q[10],q[0];
u1(1.60110410765159) q[0];
u3(0.573642466522826,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.14984401276707,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.87867321227029,0.833790283859738,0.968346586876081) q[0];
u3(2.05573715826371,-1.31457648055525,2.28489752872627) q[10];
u3(0.916259314337295,2.69897036273618,-1.44024802509467) q[7];
u3(1.02454338404516,1.21224806620368,-3.04422741737708) q[12];
cx q[12],q[7];
u1(-0.0662787238150353) q[7];
u3(-2.67452657558503,0.0,0.0) q[12];
cx q[7],q[12];
u3(1.57538676821716,0.0,0.0) q[12];
cx q[12],q[7];
u3(1.99020101192371,-3.28277810878761,-0.728830671148645) q[7];
u3(1.27819336135719,1.90722611070275,-1.95029590942652) q[12];
u3(1.69359102795523,1.03327753265756,-3.83483821474102) q[11];
u3(1.81170747201040,3.20469247944417,-2.45011523040027) q[3];
cx q[3],q[11];
u1(0.848538136006757) q[11];
u3(-3.23468887387252,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.93679108805956,0.0,0.0) q[3];
cx q[3],q[11];
u3(0.843550284186973,1.54926977354108,-2.47440124364518) q[11];
u3(2.68424463394483,2.20845711477290,-0.275452855163859) q[3];
u3(1.27070178671238,-3.53169787187351,1.60551205269684) q[5];
u3(2.03655164563865,-5.65143184234818,-0.293274636600247) q[8];
cx q[8],q[5];
u1(0.358623438530369) q[5];
u3(-0.567794506314571,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.82158963228638,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.18570049875400,1.88318193063648,0.978879247707856) q[5];
u3(2.04743334373944,-3.05673570352560,1.74170152823890) q[8];
u3(1.34078615383724,0.660238524808064,-1.65028193150668) q[2];
u3(1.57513631545215,1.49790315751679,-4.75116484379974) q[6];
cx q[6],q[2];
u1(3.79385252008382) q[2];
u3(-1.22147703072287,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.65518386013084,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.16076722955025,2.40038001686022,-0.600267638847051) q[2];
u3(1.90762785042941,-0.877882974327446,5.16596439943238) q[6];
u3(0.300559866532087,-0.958026102118697,1.38019188278942) q[8];
u3(0.572302196877334,-3.62010181400445,1.19575502377721) q[9];
cx q[9],q[8];
u1(-0.0431718113969517) q[8];
u3(0.506333760243669,0.0,0.0) q[9];
cx q[8],q[9];
u3(4.16411281493250,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.11056129560006,-2.44711575690214,3.65008771054748) q[8];
u3(1.64690277375968,1.69878262016724,-0.789663917418272) q[9];
u3(2.66231862500431,1.35237473253719,-0.119301114023907) q[0];
u3(2.06513941356358,0.992072769448363,-3.20342728718558) q[3];
cx q[3],q[0];
u1(1.75763161294135) q[0];
u3(-2.53865344552149,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.24280420254181,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.779177928692089,1.73813998235105,-2.08579024905496) q[0];
u3(1.15981049551068,0.773796597093807,-3.35429120311212) q[3];
u3(2.44802215027119,0.945076717893011,1.07678084607830) q[4];
u3(0.719526826656851,-3.39135023650118,-1.16989366939158) q[5];
cx q[5],q[4];
u1(2.30815657440226) q[4];
u3(-3.03115977281233,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.26292393415474,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.615493728943061,1.24258689966624,2.95252025585949) q[4];
u3(1.61064427264054,-2.40658585502806,-1.54704234774134) q[5];
u3(1.56078001577861,-0.199559623187302,-1.08846456202353) q[6];
u3(1.80817099637460,0.967917015750038,-5.05029993922793) q[12];
cx q[12],q[6];
u1(3.65713860165769) q[6];
u3(-3.26743216530802,0.0,0.0) q[12];
cx q[6],q[12];
u3(-1.09746495310186,0.0,0.0) q[12];
cx q[12],q[6];
u3(2.22070341469410,-2.07900393568300,1.28742356163735) q[6];
u3(1.72150896244472,0.603030937558058,5.65050616700573) q[12];
u3(2.07750771753527,2.11547934222943,0.506373364274799) q[11];
u3(1.62863061678378,0.311324956797383,-2.16143377000812) q[1];
cx q[1],q[11];
u1(3.38324475105080) q[11];
u3(-1.37581602725547,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.68477131015144,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.76717543403633,-0.914297562260152,1.11403360588785) q[11];
u3(1.41734396335261,4.28031283486149,1.11990223193498) q[1];
u3(2.15185954868484,2.01327505080068,-3.40173404377608) q[10];
u3(1.17085028242612,-2.07453105785576,2.54287166703407) q[7];
cx q[7],q[10];
u1(1.19695763127882) q[10];
u3(-0.0154586221363853,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.18133631761884,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.28022466478057,2.21130423741496,-1.49717770656892) q[10];
u3(1.00472784634804,-2.39211945792295,1.39278945954791) q[7];
u3(0.602553936934623,-2.49637264322440,1.80864976685798) q[1];
u3(0.965384564796250,1.87936258102299,-3.75368918709742) q[7];
cx q[7],q[1];
u1(0.608437485307218) q[1];
u3(-0.0870769633000357,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.66449103314880,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.753616952324528,-0.514862269156241,2.12247013754871) q[1];
u3(1.28786874381048,-1.97961818385920,-0.149969514474705) q[7];
u3(3.08859752036623,2.10440408612623,-0.443190290708600) q[12];
u3(2.97077491428134,4.90224902695703,0.332714043512186) q[2];
cx q[2],q[12];
u1(1.74670717751999) q[12];
u3(-2.74134381228235,0.0,0.0) q[2];
cx q[12],q[2];
u3(1.07780648597306,0.0,0.0) q[2];
cx q[2],q[12];
u3(0.395465115174783,-3.19689324070039,0.363313156850210) q[12];
u3(1.30864219951056,-0.906347326242248,2.69024184296963) q[2];
u3(2.34837328234825,-1.00019444258763,1.15654478663135) q[3];
u3(1.87231824428645,-2.79024951820947,-0.652581029841350) q[6];
cx q[6],q[3];
u1(3.79145948891254) q[3];
u3(-4.30675144804574,0.0,0.0) q[6];
cx q[3],q[6];
u3(-0.854337551066163,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.32736680610441,-3.58921952839772,2.50800290352667) q[3];
u3(0.844392694801133,2.44268044556280,-2.82983877477395) q[6];
u3(1.40159612761935,0.198657173113144,-1.10769801930745) q[4];
u3(0.487867121781021,-3.34732698626320,1.07207355211122) q[11];
cx q[11],q[4];
u1(1.13917948577972) q[4];
u3(-0.634170796091646,0.0,0.0) q[11];
cx q[4],q[11];
u3(3.04897979586493,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.15013528412300,1.93603323816493,0.250749153072600) q[4];
u3(2.25115020458788,-1.75330623455205,-2.82071548967407) q[11];
u3(1.47477771561015,-0.613800653749591,-0.930893921892402) q[9];
u3(2.07237611055901,1.75334765559406,-3.78229841405086) q[10];
cx q[10],q[9];
u1(1.36487664676267) q[9];
u3(-3.52726798644811,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.22027138841718,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.981657253093215,3.00667774761604,-2.60692282412606) q[9];
u3(2.02470080526601,2.65550381560715,-1.54901808062766) q[10];
u3(1.39735345803297,-1.67586944089547,-0.249897608618396) q[0];
u3(0.171623775084269,-3.72113270681424,-0.719920362394691) q[8];
cx q[8],q[0];
u1(2.26103636436043) q[0];
u3(-3.00756222206231,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.977333282886944,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.905721417949847,0.0396468761278190,-2.25344446399960) q[0];
u3(1.59208535808031,-5.84067374018205,0.0383999863632214) q[8];
u3(2.18831017944070,1.65381006971409,0.887072159456604) q[2];
u3(2.00108209185504,-0.183741584039897,-3.74932112231943) q[4];
cx q[4],q[2];
u1(1.27504432501750) q[2];
u3(-0.862101889533802,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.67486937157434,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.63184941277951,-2.27150923829076,2.07686625367167) q[2];
u3(2.17332279994712,-0.924049739163030,2.65604776370276) q[4];
u3(1.06189728197734,-2.59397553973783,1.72497924248081) q[3];
u3(0.420089651624436,-3.79512406326745,2.01895312841113) q[5];
cx q[5],q[3];
u1(0.741281106506388) q[3];
u3(-1.30443024253285,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.308181715438397,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.811535754583093,-0.723791402399901,1.96693853461908) q[3];
u3(2.46629055352297,-2.46351349409546,-1.95931022154573) q[5];
u3(2.31283894322769,-2.42206360022679,3.24860430268206) q[6];
u3(0.975260931555588,2.57548592859829,-1.81794117657983) q[11];
cx q[11],q[6];
u1(-0.0172191318608530) q[6];
u3(-1.66458234779939,0.0,0.0) q[11];
cx q[6],q[11];
u3(0.887358811307944,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.83825605899770,-3.01728693616483,2.39286684121486) q[6];
u3(2.22051100065766,0.701638736418687,-2.64261366773050) q[11];
u3(1.53450846656673,3.09030615705919,-2.67867062805407) q[9];
u3(2.51765977038368,0.591394141212458,-2.54108137257259) q[12];
cx q[12],q[9];
u1(2.73552158961350) q[9];
u3(-2.29632918377984,0.0,0.0) q[12];
cx q[9],q[12];
u3(1.27943151408504,0.0,0.0) q[12];
cx q[12],q[9];
u3(2.49622050127567,-1.34125275432809,0.911501599557175) q[9];
u3(0.779422969635166,-1.80648783760839,3.16630500868294) q[12];
u3(1.70289309495740,2.04363024575953,-3.23529056148363) q[10];
u3(2.02671714450191,-3.01756809214409,2.93098391530342) q[0];
cx q[0],q[10];
u1(-0.0720980862516674) q[10];
u3(-2.14266681321332,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.73015101088579,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.00747837764946,0.428480436055694,-2.90066012942135) q[10];
u3(2.15283601902492,3.00350431370486,0.496679881693381) q[0];
u3(0.666139693858443,2.28745124307530,-3.81358482444666) q[8];
u3(1.41413383704389,3.30884664129339,-2.23603167147279) q[7];
cx q[7],q[8];
u1(1.65183030941154) q[8];
u3(-0.0489857560348195,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.506801880137538,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.53844767818895,2.62240631549259,-2.01475644363783) q[8];
u3(0.164861013462798,3.11392228435731,-1.09145680421521) q[7];
u3(1.77738550579528,0.978322929322617,-2.83220270175057) q[7];
u3(1.22368500779312,-2.26206595998549,2.51592756248049) q[0];
cx q[0],q[7];
u1(1.47481442885450) q[7];
u3(-3.02293058042468,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.52137126656155,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.56312548677119,1.55103529825050,-3.18890373335784) q[7];
u3(2.21922719312212,-5.39161955803476,0.870200653199063) q[0];
u3(2.82232613464804,-3.00223306669309,2.35309582255776) q[12];
u3(1.31195478137230,3.57712190730273,-1.57550129733869) q[9];
cx q[9],q[12];
u1(-0.0177928771049527) q[12];
u3(0.464633511282300,0.0,0.0) q[9];
cx q[12],q[9];
u3(4.14223925027102,0.0,0.0) q[9];
cx q[9],q[12];
u3(2.71903633967206,-3.03531144075036,2.66576753996292) q[12];
u3(1.26848931517034,-0.370733167254879,5.39068258303503) q[9];
u3(2.45742936117510,2.37920342352670,-0.179765114158108) q[1];
u3(2.15178238840754,0.518134615859191,-3.04694940634305) q[4];
cx q[4],q[1];
u1(2.49591363292145) q[1];
u3(-2.97927330582446,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.43659230894403,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.02142177651698,-0.955904942681041,0.603648020612864) q[1];
u3(0.763686546653503,-0.676384017340991,-1.84968677381209) q[4];
u3(0.823276836739625,1.42375870976612,-2.33560609055559) q[10];
u3(1.75139812047768,-2.42347225586844,3.31389303588295) q[3];
cx q[3],q[10];
u1(3.30008784126084) q[10];
u3(-1.52608036264603,0.0,0.0) q[3];
cx q[10],q[3];
u3(0.955838006713247,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.98873683455923,-0.931360719240255,2.40068954712093) q[10];
u3(2.55930675610008,-0.101115910670373,4.37711551344784) q[3];
u3(0.782007210763492,2.18355678679398,-1.02455886240780) q[6];
u3(1.02581966859341,0.360393627229032,-2.79161421548531) q[5];
cx q[5],q[6];
u1(0.120049865374867) q[6];
u3(-1.33703733017822,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.94457709022543,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.02734858944034,1.17398740311436,-0.294614653745805) q[6];
u3(1.88829349226243,-4.18193436508678,1.06286523503210) q[5];
u3(0.995531426494304,-1.41116616722768,-1.31616907650168) q[2];
u3(1.42948730294566,1.50475616072223,-4.51519796750242) q[8];
cx q[8],q[2];
u1(2.22387108000192) q[2];
u3(-2.50541606476710,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.76736936557357,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.02595125079182,-2.81556025132426,-0.912614843806858) q[2];
u3(1.32912922148220,3.45811910387770,0.632268863346895) q[8];
u3(2.99339525131162,-3.46386967742161,0.604401921915985) q[0];
u3(2.71243698073306,0.938844217103268,2.56288273029672) q[12];
cx q[12],q[0];
u1(3.55152346150261) q[0];
u3(-3.69172427665342,0.0,0.0) q[12];
cx q[0],q[12];
u3(-0.707761110568435,0.0,0.0) q[12];
cx q[12],q[0];
u3(2.57957832741265,0.233232449909807,-3.16714702355044) q[0];
u3(0.426819748223506,4.49580821024005,0.134429785269671) q[12];
u3(1.54872977667209,1.43860172968970,-3.44487333313591) q[1];
u3(0.939663452383280,-1.67679182254676,2.17868347950085) q[2];
cx q[2],q[1];
u1(1.71188472533701) q[1];
u3(-0.0571302013424395,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.07263484527497,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.58655724162429,0.117131394918018,2.47510527439099) q[1];
u3(1.97309120360156,2.95408742273888,1.98594929173933) q[2];
u3(2.08820879763877,-1.25325866296377,0.0961888227568591) q[3];
u3(0.898482958783511,-2.69439734886635,1.18775583277993) q[6];
cx q[6],q[3];
u1(1.68680645009783) q[3];
u3(-2.93872613493672,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.09046083444652,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.97778448516639,-1.78956600057274,3.07551731113459) q[3];
u3(1.70014604716556,4.22195688235252,1.99822898652163) q[6];
u3(1.55373415230230,0.777846807148104,0.381638085098034) q[4];
u3(0.556635044420346,0.0479657230973913,-4.22221460032529) q[7];
cx q[7],q[4];
u1(1.06283088440666) q[4];
u3(-0.747343636640773,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.55313207902494,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.807734446812841,0.268053066562292,-2.54963490725953) q[4];
u3(0.396062379876555,2.02174468441285,1.72366446091285) q[7];
u3(1.69895948443007,0.619042746481294,1.13700072512500) q[10];
u3(0.907225860246755,-2.72763134501282,-1.09098422691187) q[8];
cx q[8],q[10];
u1(2.65812726583093) q[10];
u3(-1.63547227690313,0.0,0.0) q[8];
cx q[10],q[8];
u3(3.20773263162819,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.10070625674236,0.464536097247707,-1.14068600178967) q[10];
u3(0.601254239028557,-0.121086989068028,5.52884477394572) q[8];
u3(0.757964868992210,0.903850693902084,-2.95742076668892) q[9];
u3(1.15434872697049,-4.92955633754282,1.24730955887817) q[5];
cx q[5],q[9];
u1(2.87039556567476) q[9];
u3(-1.71626802717549,0.0,0.0) q[5];
cx q[9],q[5];
u3(0.887899024875892,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.988365324277992,1.70508529323939,-2.22284025553657) q[9];
u3(0.303383373003306,-2.64502732592280,2.73001150397319) q[5];
u3(1.33382082770959,-0.850868466344285,0.350851004797104) q[12];
u3(1.15696967336821,-1.33220240934767,-1.42402340448077) q[6];
cx q[6],q[12];
u1(1.68587286635864) q[12];
u3(-2.59580191993593,0.0,0.0) q[6];
cx q[12],q[6];
u3(0.0149972917860131,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.95979391911015,-1.70629733555296,3.16846101694637) q[12];
u3(1.44524229005441,-0.0711760822953442,6.08149924934022) q[6];
u3(1.92165570537236,0.856501841938230,-3.10935912726896) q[3];
u3(1.73694905181908,2.80893745022463,-2.88623938052260) q[8];
cx q[8],q[3];
u1(-0.442870825374171) q[3];
u3(1.21803498789852,0.0,0.0) q[8];
cx q[3],q[8];
u3(3.68910565490105,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.21542238433829,-1.70313883985802,-0.0239066330714778) q[3];
u3(1.34384154827318,-4.47808984617663,0.0529959253575556) q[8];
u3(1.01616659236579,2.90418013742824,-2.21415790853097) q[5];
u3(0.998056701351446,1.41188218314454,-1.31550406986700) q[2];
cx q[2],q[5];
u1(-1.02694419183033) q[5];
u3(0.546145794949470,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.39067266830715,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.20743324371021,0.153688152360358,1.41845608835082) q[5];
u3(2.54429378519184,0.893168577724983,-0.783551552257213) q[2];
u3(0.659181895925789,0.676432703904488,-0.572739370400958) q[9];
u3(1.04702804530512,-0.955315406535375,-0.834367847226947) q[0];
cx q[0],q[9];
u1(-0.194969374539024) q[9];
u3(-1.94822213193015,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.37261870237655,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.49244849028262,4.54174126140503,-1.63134210124232) q[9];
u3(1.50578615184548,0.767748927842111,1.28963816866896) q[0];
u3(2.03575996313003,0.900163207654956,0.543562523040273) q[11];
u3(2.69795428633782,-0.467905537172149,-3.34085279281290) q[1];
cx q[1],q[11];
u1(0.147126791856792) q[11];
u3(-1.34744666852483,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.15304919338426,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.601030052968139,5.02848927694221,-0.714872014757824) q[11];
u3(1.16822273902313,0.232581903628176,5.45878008013373) q[1];
u3(1.73928504177167,1.15805839477653,-1.21816718436191) q[4];
u3(1.98051125578371,-0.749140307629975,-3.09564570468840) q[10];
cx q[10],q[4];
u1(-0.746530188278147) q[4];
u3(1.37638989808391,0.0,0.0) q[10];
cx q[4],q[10];
u3(3.69395071818108,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.23443546750961,3.06308748766543,-2.19244952100881) q[4];
u3(2.24122594052175,-0.980343254583523,4.82005242714018) q[10];
u3(1.46599793950154,-1.37412613748622,0.732729247916818) q[5];
u3(1.34890918332122,-3.28516590367908,0.00727915034416537) q[2];
cx q[2],q[5];
u1(2.11741946779302) q[5];
u3(-2.38568372275235,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.625196487271611,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.27702746632518,-3.71452846487198,1.82918345063622) q[5];
u3(1.92344656438401,-5.41824533246854,0.787154653687496) q[2];
u3(1.18392178278872,2.80424088866875,-0.938748301030973) q[8];
u3(1.94354240458528,0.289750020551089,-3.23622546881562) q[9];
cx q[9],q[8];
u1(2.80639872302329) q[8];
u3(-1.82967669263927,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.909576173152499,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.71079918182688,1.87486725835434,-2.51574071698699) q[8];
u3(1.47682364787694,-3.73670954852320,-0.0379453070637250) q[9];
u3(1.92514622862895,-0.131571443072263,-1.13790129059955) q[12];
u3(1.64063167223996,-4.96847172369380,0.979574044351056) q[11];
cx q[11],q[12];
u1(1.87930111988652) q[12];
u3(-2.96524436609321,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.670443624320424,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.81281189625154,2.53210248585912,-3.50423570769980) q[12];
u3(2.83642188726078,-3.20525707533849,1.25256434683629) q[11];
u3(1.24719955298837,0.934305425046525,-2.53482762882024) q[7];
u3(1.78625522749945,-2.65502819575770,2.85532871206912) q[4];
cx q[4],q[7];
u1(3.41068026326018) q[7];
u3(-3.22964313804854,0.0,0.0) q[4];
cx q[7],q[4];
u3(0.791967284574660,0.0,0.0) q[4];
cx q[4],q[7];
u3(3.07596102605489,2.06081415762671,-0.504341930316881) q[7];
u3(1.47217367320196,-1.19800563131481,-3.84517895724134) q[4];
u3(1.19589971039268,1.84226880677068,0.123350646703275) q[1];
u3(1.28646986059176,0.432057858037118,-3.39382133133477) q[3];
cx q[3],q[1];
u1(1.32733021030667) q[1];
u3(-0.759716916406402,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.414454219748583,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.34160389729799,-0.992153882716802,-2.14115381761889) q[1];
u3(2.32772507522580,3.22618698944756,-2.90078481874480) q[3];
u3(0.846211410036812,0.885269530960760,-3.51342588865505) q[0];
u3(1.70238223787697,2.91845607933750,-2.42333387267822) q[6];
cx q[6],q[0];
u1(-0.320027430673957) q[0];
u3(-2.44045805547802,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.44692775517894,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.29057431240877,3.05844716084985,-2.68684545473684) q[0];
u3(1.21886584119509,0.548473544466499,0.938119711360945) q[6];
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