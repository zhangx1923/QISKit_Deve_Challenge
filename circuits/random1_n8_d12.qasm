OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.60450443578975,1.02930610904539,-3.44034698881355) q[2];
u3(2.08835043962235,3.25674384513931,-2.67528806501547) q[6];
cx q[6],q[2];
u1(3.49421036064838) q[2];
u3(-3.68709899963900,0.0,0.0) q[6];
cx q[2],q[6];
u3(-1.22491494498164,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.41236734824494,-0.812810732320788,-2.16153908121607) q[2];
u3(1.17686337985840,1.82721660185379,2.82976212990518) q[6];
u3(2.05602944245824,1.81355169855145,-2.25051265880974) q[0];
u3(1.86049114805846,1.11067726776313,-1.25354365849235) q[1];
cx q[1],q[0];
u1(-0.0346689880435274) q[0];
u3(-1.91025807261151,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.667506763479480,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.694392768017168,2.00912395105407,1.08605420224007) q[0];
u3(0.814975155893419,1.54979315750787,4.10498075438444) q[1];
u3(1.25859689094588,-2.14007352554185,1.15020401772245) q[3];
u3(1.92450271016076,-3.52261910090879,-0.0496293441109339) q[7];
cx q[7],q[3];
u1(1.37067962959334) q[3];
u3(-0.525041960499737,0.0,0.0) q[7];
cx q[3],q[7];
u3(-0.230696659769804,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.14977089022159,0.568140686686336,0.143655588529228) q[3];
u3(0.468843778012453,-0.177835239149770,2.95692599626687) q[7];
u3(2.47569911435295,-2.60754448083544,0.592045382990912) q[4];
u3(2.65124332541900,-0.538376218750941,0.102039720967505) q[5];
cx q[5],q[4];
u1(0.892204244137902) q[4];
u3(-0.329922015618464,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.54866140829067,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.64736320523955,1.82736141604560,-0.105199364542973) q[4];
u3(2.12078372884754,3.55910532574707,-1.86917658520344) q[5];
u3(1.47181763562252,-0.0417077339538401,0.221902447663174) q[7];
u3(1.29150653357804,-3.12642254412586,-0.484928922534547) q[0];
cx q[0],q[7];
u1(1.38586227169551) q[7];
u3(-0.920788038282560,0.0,0.0) q[0];
cx q[7],q[0];
u3(-0.333292957340362,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.28256198792704,-2.73952929107571,0.484525291668925) q[7];
u3(1.60054438786983,0.179396663253333,3.55823815707504) q[0];
u3(2.49497093460349,0.755314838004557,1.05299359985912) q[6];
u3(1.83373457847067,-1.75742065052129,-1.36650422971714) q[5];
cx q[5],q[6];
u1(2.78466001497645) q[6];
u3(-1.52552397560602,0.0,0.0) q[5];
cx q[6],q[5];
u3(-0.0406404101526758,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.769408031255386,1.12438894126678,-0.262544063596663) q[6];
u3(2.58348001233906,-0.0195293333785496,-4.20085927608245) q[5];
u3(0.525536620338778,1.75368066283938,-1.49596740854174) q[2];
u3(0.726019306915123,-2.75921846180264,1.12990067625231) q[1];
cx q[1],q[2];
u1(1.44300246319329) q[2];
u3(-2.29374036806372,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.0179692376466336,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.29731720983601,-2.35420572694653,1.65155619391603) q[2];
u3(0.817641946414384,-4.37931371331668,-1.67504128625286) q[1];
u3(1.37307926966801,1.14953990442268,0.349096818864426) q[4];
u3(1.98686363331420,-1.18705673449559,-2.43077367721450) q[3];
cx q[3],q[4];
u1(1.48285294575474) q[4];
u3(-3.67093320491507,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.46922916878126,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.43979503681061,2.02375659981056,2.30466580912430) q[4];
u3(0.802241347332433,1.16467902320881,-1.70217592727724) q[3];
u3(1.82355745658197,2.57958367396373,-2.65165835338208) q[3];
u3(0.823211817592125,2.63475017301571,-2.40020703239651) q[1];
cx q[1],q[3];
u1(0.177541675846472) q[3];
u3(-0.565331888790364,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.98054316085387,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.03423424347553,0.0241988379760660,1.51304691393829) q[3];
u3(2.07674290844323,-1.67266502005874,-0.572939343370721) q[1];
u3(0.447783962777568,-1.15824292100143,0.914991717672356) q[4];
u3(0.902155877052750,-2.62006038596098,-0.335935341883528) q[6];
cx q[6],q[4];
u1(2.03745835835919) q[4];
u3(-3.13799884537678,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.24760930395079,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.218512924038041,2.17809848448355,-1.81135072043600) q[4];
u3(2.90310748674208,1.49308653266110,-1.83392116298692) q[6];
u3(1.26687006094086,0.902277868549147,-0.389322312722738) q[2];
u3(1.09345384943877,-0.662368694393298,-3.77178991347278) q[7];
cx q[7],q[2];
u1(2.23494031692709) q[2];
u3(-2.54493032363054,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.187748217942981,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.93902524100647,2.31156447206766,-3.79568498558815) q[2];
u3(1.74853985499800,0.913969619185820,0.324320267380795) q[7];
u3(2.95951454807384,1.13585976017566,-1.56820739692674) q[5];
u3(1.64934579157267,3.94745946379981,0.273105027069666) q[0];
cx q[0],q[5];
u1(2.36543097934886) q[5];
u3(-3.05815862201474,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.37745298114276,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.786363006213835,-2.61572774503556,1.52534322875622) q[5];
u3(2.62735622833240,-0.576022561845314,0.551875415463205) q[0];
u3(1.00726532057738,-1.94372069266110,3.17633629157427) q[3];
u3(0.234249176012098,-2.81661681798027,1.52575067212042) q[6];
cx q[6],q[3];
u1(2.70130352926813) q[3];
u3(-1.57574628521318,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.19670965292332,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.81612460407707,2.49048397883966,-1.40406608605773) q[3];
u3(1.27826053546673,-0.293113523709282,3.40665906940026) q[6];
u3(1.76644098019731,-1.17895601550387,0.812395495749214) q[1];
u3(2.12816828653639,-1.28575175617776,-1.47442312852208) q[4];
cx q[4],q[1];
u1(-0.253643064586247) q[1];
u3(-1.56473583852282,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.08011486230125,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.52406377809647,-1.18000634848878,2.29219797035777) q[1];
u3(1.06924758556231,2.45864607382405,-0.0111943655992728) q[4];
u3(2.16008714459478,-2.34792755092330,2.44735139883278) q[2];
u3(2.39446082050546,-2.09288233777164,1.89808041981553) q[0];
cx q[0],q[2];
u1(2.90886550870406) q[2];
u3(-1.62245404498605,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.560461678642104,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.34096128587775,1.14740184408807,-0.671967802169743) q[2];
u3(0.858270903106303,2.77440876222295,3.14554932721631) q[0];
u3(0.402132082031181,-0.847361200420655,1.02452797703612) q[7];
u3(0.390931101561637,-1.54372767225354,-0.557688323187784) q[5];
cx q[5],q[7];
u1(-0.214354387643469) q[7];
u3(-0.703210442840941,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.58022025826652,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.93892577519787,-0.536514949482384,0.0944220850587044) q[7];
u3(1.76756125427814,-1.55642777556790,3.23311757226214) q[5];
u3(1.89825073232420,0.580114399687138,1.06248712427430) q[7];
u3(0.303717822634992,-4.74646727458180,-0.124920105017701) q[6];
cx q[6],q[7];
u1(0.581307669412324) q[7];
u3(-1.46928596997295,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.13322593398735,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.473907315405289,-1.21465692337463,-0.282498118780111) q[7];
u3(1.31874758982468,-0.0107466380832899,5.53586016620446) q[6];
u3(1.41344071370900,-4.30430147535112,1.57862394956432) q[4];
u3(2.02556707350223,-5.05381001241744,-0.0195272183988839) q[2];
cx q[2],q[4];
u1(1.82045831240944) q[4];
u3(-2.52441369028255,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.14999335822774,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.72133066973185,1.50357048320112,-3.67771805705281) q[4];
u3(2.92550642522603,0.744307031603389,-5.08617200091367) q[2];
u3(2.73741642544102,4.03851233914176,-1.83501099856638) q[0];
u3(1.33128465385482,1.60212850365355,-1.79267655272136) q[5];
cx q[5],q[0];
u1(1.60240548077159) q[0];
u3(-2.66417514078179,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.33748370681988,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.906359749167037,-2.49970060093903,-0.491284474835433) q[0];
u3(1.92159266891122,0.345730557789774,2.46071023533047) q[5];
u3(0.706706934377278,-3.98174048865753,1.90612823077928) q[1];
u3(1.37615368310235,2.69558464487760,-2.74381230899406) q[3];
cx q[3],q[1];
u1(1.54084972623470) q[1];
u3(-2.43939788490633,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.36378598880204,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.74952968813210,-2.18090569235470,0.0278601102448277) q[1];
u3(0.566439039106253,-1.23768192187511,-4.48916980312412) q[3];
u3(0.786848495837463,0.444416927384769,-1.60115258896426) q[1];
u3(0.912913089492388,-4.55348775843204,1.72814939669268) q[6];
cx q[6],q[1];
u1(2.01637382955690) q[1];
u3(-3.42067596657158,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.37625446650488,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.89909109529393,-1.67101220427687,3.82633281580939) q[1];
u3(0.789551273781293,1.68243479918484,-1.89335056011627) q[6];
u3(2.43093508157666,-0.510307744519065,-1.09649178884739) q[5];
u3(0.233406611272859,0.541031505734597,-5.51279182180022) q[3];
cx q[3],q[5];
u1(1.78027967740891) q[5];
u3(-2.71129337674031,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.557352901972526,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.306717865268021,-3.75329207383653,0.849393809096891) q[5];
u3(2.44927734493115,4.31542869264138,-1.41392055803093) q[3];
u3(1.94920611056082,2.91323587592878,-3.02322649144300) q[2];
u3(0.703218360864796,3.20222461712162,-1.88947026807779) q[7];
cx q[7],q[2];
u1(3.06674521353051) q[2];
u3(-1.79250413845022,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.703625181010105,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.693047253872550,-0.659841515837191,3.40862624335582) q[2];
u3(1.45879968186872,-3.25079511528047,2.58192478944766) q[7];
u3(1.76596998397123,1.86067947238481,1.25415239801437) q[4];
u3(2.44009800522845,0.529679516595748,-2.59752765539656) q[0];
cx q[0],q[4];
u1(-0.403552435590415) q[4];
u3(-2.23088140488875,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.57854281131210,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.93764161144242,-3.69363299681718,0.112983298963972) q[4];
u3(0.709472015928755,-0.935220152233780,-1.71490440855345) q[0];
u3(1.20067386326370,-1.61999826101868,1.77635403262318) q[6];
u3(0.774210784032875,1.62105654614151,-2.19587580937057) q[7];
cx q[7],q[6];
u1(1.43117239250647) q[6];
u3(-3.59816051158460,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.42518203332453,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.30042138438299,4.10813893325599,-1.23660774486335) q[6];
u3(1.67905365660356,0.258703425267190,0.626704358750960) q[7];
u3(1.84572076319697,-1.76797077699052,0.398181295350581) q[4];
u3(1.34637633214532,-3.87042418551419,-0.728195564947705) q[3];
cx q[3],q[4];
u1(2.67566997717672) q[4];
u3(-1.83214001046124,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.430866408568318,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.590675200163828,-0.559468722480909,4.82890390449835) q[4];
u3(0.378913339216233,-0.190740593589537,-5.87363150259998) q[3];
u3(0.492371828856711,0.685665483759755,-1.19573815119371) q[2];
u3(0.656991874527479,-1.24603826409645,0.652555732683038) q[5];
cx q[5],q[2];
u1(1.77333343851701) q[2];
u3(-0.712238392616559,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.31312975657196,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.12352389817273,2.98484314720447,-2.73920588320603) q[2];
u3(1.05343073513744,5.58149264710277,0.290628513919029) q[5];
u3(1.88978543594327,-0.0617062312468839,2.53670592666670) q[1];
u3(2.26722907112050,-3.27113680747259,-1.80959374428172) q[0];
cx q[0],q[1];
u1(1.39050049678034) q[1];
u3(-0.297992415463943,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.04843141206138,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.32617918983552,0.689551972209881,1.24417447860179) q[1];
u3(2.63833936377276,-1.40897784098802,-3.21232784243170) q[0];
u3(2.32790349369903,3.22768039404225,-2.72418136687525) q[0];
u3(0.487774544903641,0.888723696214413,1.29009645836806) q[6];
cx q[6],q[0];
u1(0.341733848735553) q[0];
u3(-0.999212636745519,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.80899089069888,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.12055475948816,2.24022102889987,-2.70965342406147) q[0];
u3(0.793200104900826,-2.01278039988024,2.75098215318385) q[6];
u3(0.581579038144823,-2.11480189934749,1.18335221055559) q[3];
u3(0.295882506280146,2.01073933256625,-4.12603296791075) q[5];
cx q[5],q[3];
u1(1.64277729683604) q[3];
u3(-2.22401462154844,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.438452167466367,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.403317213299961,-1.14670822167669,-0.596093027206045) q[3];
u3(1.47305649964616,1.02330496513885,-2.92047888488502) q[5];
u3(2.96502928812633,0.617025283150539,0.178272883260330) q[2];
u3(1.36938543267143,0.692165566803856,-4.64523028114920) q[7];
cx q[7],q[2];
u1(0.113033969040307) q[2];
u3(-1.67171533749519,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.697957491444101,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.42747498957080,1.49302001659237,0.664850312015061) q[2];
u3(2.27955399937839,-0.376199339617357,-5.40980935240294) q[7];
u3(1.29381681537117,1.43671193647041,-3.94623002408519) q[1];
u3(0.687798220020655,2.26628225848673,-1.93745384713255) q[4];
cx q[4],q[1];
u1(1.53254055996543) q[1];
u3(-3.47260814739525,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.15666685103103,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.29734669945141,-0.443444506161404,-0.346401686100994) q[1];
u3(0.640422441515180,4.59507252234273,0.0356683292340128) q[4];
u3(3.07947770715168,-0.906641952843423,-0.771773101640609) q[7];
u3(1.12768738375240,-0.200050271503729,-4.68428085707064) q[6];
cx q[6],q[7];
u1(0.131186884931217) q[7];
u3(-0.982961887537933,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.68586298322229,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.862156600111988,1.76154842736843,2.16673157873632) q[7];
u3(1.96204547343407,2.69464663232376,-2.67363126924846) q[6];
u3(2.30313111253953,2.14613631265420,-1.91841197663193) q[5];
u3(2.20280228121573,2.32710055867568,-3.16240984611081) q[0];
cx q[0],q[5];
u1(0.201222403437754) q[5];
u3(-0.831014079101222,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.83792989970830,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.197130556859436,-1.26153279283233,-1.76742225497309) q[5];
u3(1.90466735110988,-1.46876181800321,3.00479189378619) q[0];
u3(1.22871234904114,2.01915721250815,-3.02435237156562) q[2];
u3(0.642241645193975,2.88996229473236,-2.16495132449195) q[4];
cx q[4],q[2];
u1(1.56597770160174) q[2];
u3(-0.195401288883962,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.68707611796419,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.470234345503177,2.41083975458973,-2.32960104923262) q[2];
u3(1.92391265112642,0.253378831464710,3.59850994537827) q[4];
u3(2.17423557961247,-3.61812834829353,1.87062114109097) q[1];
u3(0.427605776646897,3.83870792918741,-1.38368669230099) q[3];
cx q[3],q[1];
u1(3.49094200747531) q[1];
u3(-0.513644108561886,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.43864819403330,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.19139683769898,0.343759385199556,-2.12323675685680) q[1];
u3(2.55092303610236,1.47570306769259,-3.15510377907268) q[3];
u3(2.41467807517386,0.195639408286726,1.99466123963603) q[7];
u3(1.70691771681385,-2.79678015279452,-2.62741253359298) q[0];
cx q[0],q[7];
u1(0.120060929245950) q[7];
u3(-2.28620814652222,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.15085863229600,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.39156057407384,-0.923286291653480,2.62863529756990) q[7];
u3(0.572310241535683,-1.24905797590754,-4.74102207586913) q[0];
u3(2.29312508415057,0.239262655084472,2.76750999161826) q[4];
u3(2.16044630738315,-1.98772056752519,-2.34292366826648) q[2];
cx q[2],q[4];
u1(3.22117150668057) q[4];
u3(-2.17987390172266,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.62274756067123,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.922162481106654,-0.664625911623847,-1.95299907838044) q[4];
u3(1.68281058284268,-2.74338324230780,-0.305030384433710) q[2];
u3(1.53823291354560,1.22268480910231,1.01112021741857) q[6];
u3(2.22258083073403,-0.770177334861818,-3.71271183953517) q[1];
cx q[1],q[6];
u1(-0.106965095367048) q[6];
u3(-2.06071601299142,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.823327296103790,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.20636793223286,-1.90244043634544,1.71675913555644) q[6];
u3(0.856998134872434,0.588805751670935,1.92159713568928) q[1];
u3(1.32552215760726,1.32294966582740,0.447592848947795) q[5];
u3(1.03025841628308,0.0566305831488703,-2.05892004125850) q[3];
cx q[3],q[5];
u1(3.06559335633548) q[5];
u3(-2.59665792329484,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.709584955111142,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.45824331580795,2.41822179460034,0.584278284510494) q[5];
u3(0.981305454560602,1.24658960615010,-1.96395230805115) q[3];
u3(2.38291211266603,2.49072577783976,-3.51317003800961) q[7];
u3(0.741790951878935,3.36655539240315,-2.57862142349546) q[1];
cx q[1],q[7];
u1(-0.0366283853133442) q[7];
u3(-2.25380508854357,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.55573518020526,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.27220694308446,-2.31688340705822,3.81220518556079) q[7];
u3(0.573005422095620,-0.534598666707475,2.03609176659823) q[1];
u3(1.84570090719217,1.36831107791683,0.228621030446445) q[3];
u3(0.541139509456860,-0.773331226065175,-3.78390669384379) q[5];
cx q[5],q[3];
u1(2.90889357397503) q[3];
u3(-2.25083131985641,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.948131698881541,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.33524357162759,-2.82855585334100,1.70626850366997) q[3];
u3(1.91419259566467,1.84064103892149,-4.43876014485088) q[5];
u3(2.01014658996097,3.26329177095970,-2.88112167043293) q[6];
u3(2.40553258745736,2.01580482540864,-1.15628527629979) q[0];
cx q[0],q[6];
u1(1.53686267314429) q[6];
u3(-0.0330148401562429,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.474442897899837,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.801399282285698,2.16875141489180,-2.08672384021147) q[6];
u3(2.68409179190048,-0.842836674523606,1.72614312783225) q[0];
u3(1.10172642670407,0.734870354618691,-1.83432567802963) q[2];
u3(0.771339681245419,0.847664606012054,-3.76519492554743) q[4];
cx q[4],q[2];
u1(3.40565146179586) q[2];
u3(-1.38905425828901,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.58966257084836,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.37734263343486,0.747463234239292,1.76742217227934) q[2];
u3(0.702984527793608,-1.58795828299385,-4.16853214921586) q[4];
u3(1.15908663185048,-0.375988126953301,1.44269744563445) q[5];
u3(1.63176461097090,-0.821657651120568,-2.50643944422914) q[2];
cx q[2],q[5];
u1(0.271384751848867) q[5];
u3(-2.09646382344270,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.47860120545380,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.87989722000918,1.24111527218173,-1.85667681518445) q[5];
u3(1.27075479983350,5.18576546255288,0.494308002685643) q[2];
u3(0.863756842261461,-1.33727262368113,1.08350917969849) q[4];
u3(0.564438846664007,-2.19855528664924,0.181662191594011) q[1];
cx q[1],q[4];
u1(1.61947759645160) q[4];
u3(0.00454189045387521,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.983508814199890,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.12713898888191,-2.02208373786461,4.16665316496856) q[4];
u3(1.71471677950026,-0.358218523935435,4.84351564146604) q[1];
u3(0.131419071325160,0.606132325731096,-0.279500158003966) q[7];
u3(0.927212627944151,0.646920728201821,-1.24834898079397) q[3];
cx q[3],q[7];
u1(1.63146001305579) q[7];
u3(-2.17080638823331,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.33439506551095,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.04614902859511,0.905106948015880,0.528735196744979) q[7];
u3(1.02788483645737,1.17174114397949,-1.27221030168675) q[3];
u3(2.26340086670357,2.13864479182702,-1.89554030270345) q[0];
u3(1.13123920789955,1.63981962583544,-2.62476609376553) q[6];
cx q[6],q[0];
u1(1.44905395233831) q[0];
u3(-0.650648086482044,0.0,0.0) q[6];
cx q[0],q[6];
u3(-0.173436372969155,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.97003193939199,-4.03237836880382,0.248991359273502) q[0];
u3(2.85317800109049,-2.87247704394174,-0.659409823888165) q[6];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
