OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.07807922345503,1.44715164232173,-3.19639575352387) q[3];
u3(1.90723928859745,1.69257041935973,-2.92127753356560) q[8];
cx q[8],q[3];
u1(3.10456222673516) q[3];
u3(-2.05233261501281,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.71430609648087,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.56809046754128,2.61474843769480,-2.85905125653962) q[3];
u3(2.69490328162940,0.116019955332967,-5.78954887025594) q[8];
u3(0.879482419328366,-2.10286715353983,0.397502906661153) q[5];
u3(1.74801580263868,-4.41983866360096,0.540565641387478) q[7];
cx q[7],q[5];
u1(2.51200144419839) q[5];
u3(0.217572317369367,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.60408363866397,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.57292355962259,1.08218196813706,-0.149597229176932) q[5];
u3(2.16267012917020,-0.199407885118558,-0.489100398660746) q[7];
u3(2.88957746295071,1.05521023192444,-3.12562311382300) q[10];
u3(2.22697377305302,1.87542095719991,-2.66298250198644) q[9];
cx q[9],q[10];
u1(0.720040265908202) q[10];
u3(-1.28163240516691,0.0,0.0) q[9];
cx q[10],q[9];
u3(2.97714838826267,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.66801446444347,1.23702372108479,-2.42472427431646) q[10];
u3(1.78035598414121,0.462084083883148,-1.27761982899137) q[9];
u3(1.57552517584722,-1.01039681526683,-0.101753771574846) q[0];
u3(1.35568857844672,-2.82112748073901,-0.552599896649634) q[1];
cx q[1],q[0];
u1(-0.216332715420403) q[0];
u3(-2.27312538050415,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.26742503150398,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.61353841232448,4.50208574844336,-0.597222889663521) q[0];
u3(1.11366044484934,0.626404038338130,2.09365949507435) q[1];
u3(1.46533193183710,2.04658092826713,0.329372367311579) q[6];
u3(0.398507946549320,0.816817916134860,-3.46898946470847) q[2];
cx q[2],q[6];
u1(1.91936453192625) q[6];
u3(-1.69482374035286,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.648520052396971,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.25238604712064,3.14544249450897,-1.61873752083881) q[6];
u3(1.59805713525634,3.81035691011801,1.51845695876371) q[2];
u3(2.86474296037697,-0.845847027850335,-0.225370360552336) q[8];
u3(0.367530358274859,-5.21079268889253,0.617722122761272) q[2];
cx q[2],q[8];
u1(0.365619865725827) q[8];
u3(-0.942110548696957,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.45475126261927,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.16076419411483,-0.188507245447053,1.77400462435207) q[8];
u3(2.68734243903534,1.79929032647871,-0.0116222871560302) q[2];
u3(0.808200666955079,1.29011740065920,0.623663462903265) q[3];
u3(1.36607739880975,0.115898943038080,-4.00245707874755) q[1];
cx q[1],q[3];
u1(3.56381562330148) q[3];
u3(-1.36741483132815,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.29691972607012,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.750635667771059,0.930890841577523,-3.13150191572532) q[3];
u3(1.30877131509295,5.19951690906493,-0.0317616918873997) q[1];
u3(2.44848542478474,-1.76833381526521,3.99732989705446) q[7];
u3(0.223989231574777,-1.15601288170063,3.20439945418397) q[9];
cx q[9],q[7];
u1(2.51837072385935) q[7];
u3(-1.78203140345535,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.647707475728834,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.550745139196055,2.70762139083186,0.528686978728338) q[7];
u3(1.79816688074115,-3.19672899594031,-1.24216062476275) q[9];
u3(2.65007975902614,0.419862922465489,1.30380308928901) q[0];
u3(1.58089618043562,-2.39971344039719,-2.37444599865216) q[5];
cx q[5],q[0];
u1(1.43667261713353) q[0];
u3(-0.515203591442778,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.00483794534488,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.80153137676308,-3.01308167847435,1.16250833099134) q[0];
u3(2.99032378666648,0.622933622800697,-0.556570492921234) q[5];
u3(2.92429148160343,-0.451821502135870,2.72058769188054) q[6];
u3(2.42970900944097,-1.50681213098118,0.988882401188339) q[4];
cx q[4],q[6];
u1(4.12462953040043) q[6];
u3(-3.24857974117615,0.0,0.0) q[4];
cx q[6],q[4];
u3(-0.476976663006396,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.55439487080259,-4.38756696068032,1.19813437989377) q[6];
u3(0.456220082745023,-0.0387318410864081,-2.06577538142957) q[4];
u3(1.38754746872187,-2.04887294304453,1.48467268264951) q[0];
u3(2.09891107717143,-3.49314915019560,-0.0873121543902704) q[10];
cx q[10],q[0];
u1(1.69145751121330) q[0];
u3(-3.06799068501904,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.02746338555604,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.77115562048601,-1.82726262636572,3.76856871012163) q[0];
u3(2.27367860989080,3.47377962590504,0.270689921477501) q[10];
u3(1.68728410630044,1.87800890936003,-2.92060483238700) q[6];
u3(2.17811151684003,-2.55450252194570,2.84256325312164) q[2];
cx q[2],q[6];
u1(1.14349275852508) q[6];
u3(-0.610939760101030,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.37485052855215,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.000822897831423372,0.644054588906891,-0.0748209714361933) q[6];
u3(1.46435523111105,2.02519975127658,-2.20607373469183) q[2];
u3(0.954384737147574,-0.246102154207549,1.32881516115235) q[9];
u3(1.21825900390982,-2.65183338776493,0.0392552800499704) q[7];
cx q[7],q[9];
u1(4.14939590277218) q[9];
u3(-3.25532658191858,0.0,0.0) q[7];
cx q[9],q[7];
u3(-0.542047196077103,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.69612993403197,-0.191564330242616,1.48233832988792) q[9];
u3(1.40955222401587,-0.689519709277733,-0.0779042168090381) q[7];
u3(1.23704907307316,1.68608876049898,-1.83836103183035) q[3];
u3(0.684393342304569,0.882703200992226,-2.98740911093306) q[8];
cx q[8],q[3];
u1(2.13185606570874) q[3];
u3(-2.74405441023268,0.0,0.0) q[8];
cx q[3],q[8];
u3(-0.0195328238462946,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.101166947088820,-2.07131858357793,1.73407889895795) q[3];
u3(0.818682861582403,-2.49619923962576,-1.42533876544204) q[8];
u3(1.43370346785875,2.25944567667107,-2.93111341678180) q[1];
u3(1.17224752403866,2.62762088072197,-3.05998906115363) q[5];
cx q[5],q[1];
u1(1.63894283887961) q[1];
u3(-0.891895136195677,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.20676510626983,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.71702797480455,-1.20005594098647,0.697002733649885) q[1];
u3(2.25736051061501,-1.63159318345438,1.47860139229038) q[5];
u3(1.98333636008898,0.302961209602103,0.981231585873881) q[0];
u3(1.58149443787410,-2.57577936381580,-1.80452006401169) q[3];
cx q[3],q[0];
u1(2.89352575926544) q[0];
u3(-2.09137982705393,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.19094281929519,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.53056876448045,-2.26573827669047,3.43243042302187) q[0];
u3(2.61830886999910,2.38872550642156,-0.410966770487311) q[3];
u3(0.262682766216027,1.71330054334164,-2.13780370469774) q[2];
u3(0.922085645569025,-4.03133426800297,2.00408839057586) q[1];
cx q[1],q[2];
u1(2.34655299706956) q[2];
u3(-3.05684945952636,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.37218470164570,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.09660672710013,-1.19612676543081,4.07851841217246) q[2];
u3(1.37734791997670,4.05052908756043,0.925838761246290) q[1];
u3(0.437663870141616,-2.31083035083177,2.82722805351430) q[5];
u3(1.19986988032720,-2.67740495622443,1.80196878922605) q[4];
cx q[4],q[5];
u1(2.75993389040968) q[5];
u3(-1.64689258309011,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.493557632539932,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.40922704644658,0.253992612600447,3.77453753768517) q[5];
u3(1.77421604376233,0.196511649126921,-1.10932814357257) q[4];
u3(2.24973869690755,-0.0292259035751005,1.62577984856718) q[10];
u3(2.28807946400670,-1.73612241784835,-1.04495314900122) q[6];
cx q[6],q[10];
u1(1.15401975130192) q[10];
u3(-0.399136922281662,0.0,0.0) q[6];
cx q[10],q[6];
u3(3.22662007570659,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.872220334744256,1.84477358701133,-0.384504827506012) q[10];
u3(0.600841226805138,-1.68456312796158,-1.14215363982757) q[6];
u3(0.932569026281955,2.88913641542690,-1.98015555905382) q[8];
u3(0.928646836125466,2.09210358513821,-1.98590053137786) q[7];
cx q[7],q[8];
u1(-0.870783471247075) q[8];
u3(0.553712503650017,0.0,0.0) q[7];
cx q[8],q[7];
u3(3.73830782028988,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.55963423201090,0.137287077966072,-2.93564442763602) q[8];
u3(1.91205266901389,3.72106285226685,0.280736151207471) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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