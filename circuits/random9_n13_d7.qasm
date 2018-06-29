OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(0.898234099202226,0.671617000701877,-2.06675437188561) q[3];
u3(1.01862798075880,1.00408335433022,-5.14481377372266) q[4];
cx q[4],q[3];
u1(0.776238460182794) q[3];
u3(-1.46128243928027,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.67559442335984,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.719019968593518,-1.96933423638116,-1.27879144008340) q[3];
u3(1.04898746120123,-1.93386987868340,-1.36763843095150) q[4];
u3(2.37917392311283,0.458195153339690,1.71376285549244) q[10];
u3(2.06333834976618,-1.83601091352758,-2.21889011721337) q[8];
cx q[8],q[10];
u1(1.69038833882370) q[10];
u3(-1.94945720424955,0.0,0.0) q[8];
cx q[10],q[8];
u3(3.35201713846060,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.12486639332434,-0.0139270416736588,-1.13545059086211) q[10];
u3(2.79420495051341,0.718315404132638,1.58136045183025) q[8];
u3(1.07357357856859,0.0759712374164807,1.10796771835314) q[11];
u3(0.894104491209223,-1.90176090186100,-1.38101532130046) q[9];
cx q[9],q[11];
u1(2.96399540073325) q[11];
u3(-0.956898998341830,0.0,0.0) q[9];
cx q[11],q[9];
u3(2.48153785234160,0.0,0.0) q[9];
cx q[9],q[11];
u3(0.620486697146508,-0.771263944608761,3.22819023932639) q[11];
u3(1.26699670711679,-0.767926175984099,-1.61312348167577) q[9];
u3(0.898429544000317,-0.473493686598175,1.29508462314908) q[1];
u3(0.613337191621644,-2.76616202084988,1.47743246026311) q[12];
cx q[12],q[1];
u1(-0.120172447807150) q[1];
u3(-2.27589432499590,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.14059878994085,0.0,0.0) q[12];
cx q[12],q[1];
u3(2.05455963478819,-2.36549356334006,0.481477360137948) q[1];
u3(1.12068565255345,1.08433126731672,-4.68124894698574) q[12];
u3(1.84720158119206,1.01465769222630,-2.89044684610703) q[6];
u3(1.84240738806548,0.909554331446142,-4.79750748212768) q[2];
cx q[2],q[6];
u1(0.234089847432633) q[6];
u3(-1.55478316850742,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.07904771471291,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.03745032556083,-0.952278607917790,1.69841764749420) q[6];
u3(1.17850163899209,3.49310434934584,2.76337046874723) q[2];
u3(0.760123220769922,1.39304761452230,-1.83426044082168) q[5];
u3(0.311891275787087,0.258678238167680,-2.34246469531797) q[0];
cx q[0],q[5];
u1(-0.119222900955138) q[5];
u3(-2.52563670271956,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.31300161080922,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.92334499695831,1.65590240976051,-1.75115334599532) q[5];
u3(0.880129711966485,-1.73498233500730,-2.98142628384974) q[0];
u3(1.20785313867279,0.816075995391911,-3.78991933715122) q[10];
u3(1.62473934313629,3.40681281279840,-2.47233012839349) q[11];
cx q[11],q[10];
u1(2.00027917545079) q[10];
u3(-2.62208865035237,0.0,0.0) q[11];
cx q[10],q[11];
u3(3.05486619788813,0.0,0.0) q[11];
cx q[11],q[10];
u3(2.27825559947351,2.49428048946025,-0.379271471618845) q[10];
u3(1.21537439010235,-3.56806084242004,1.53821627551456) q[11];
u3(1.19934735474926,0.431485664466617,0.673320260103718) q[6];
u3(1.69638964508143,-0.416303057381554,-1.48022398967387) q[7];
cx q[7],q[6];
u1(0.613375083399674) q[6];
u3(-1.49256746049302,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.86694019298316,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.14159736561320,0.211934364564697,1.62003618842420) q[6];
u3(1.08128137289485,2.83704043194456,-0.484229243040670) q[7];
u3(2.15776711972600,-1.17265226181107,-1.30089699602206) q[8];
u3(0.898951420031936,-4.72381709372684,0.815745231382537) q[3];
cx q[3],q[8];
u1(1.13923762462665) q[8];
u3(-0.232432175563326,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.64354678842215,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.28264125648107,-0.0698423662823731,-2.64350564308830) q[8];
u3(1.96199384558753,3.56642483925611,-0.397149499392657) q[3];
u3(1.10984403812803,3.40849537483618,-0.729351323155010) q[1];
u3(1.29261919811917,1.48615531714323,-1.52713068239647) q[5];
cx q[5],q[1];
u1(-0.161690853324288) q[1];
u3(-2.53294911830791,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.12867943279828,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.688157430266343,-1.66268030753682,3.71325466021409) q[1];
u3(2.08805504540288,-5.10561865278814,0.259575932866356) q[5];
u3(2.40534840128476,0.379287598631471,2.72765677166082) q[9];
u3(1.95334912450418,3.18271585549224,2.89999411557029) q[12];
cx q[12],q[9];
u1(1.32032397908009) q[9];
u3(-3.19471748762282,0.0,0.0) q[12];
cx q[9],q[12];
u3(0.341300991501401,0.0,0.0) q[12];
cx q[12],q[9];
u3(2.05648679866858,-0.437434793613373,0.915825174878131) q[9];
u3(1.98604102987632,-0.651100846508898,-0.692514101229106) q[12];
u3(2.17024382978816,-0.159834085563135,2.32871618932305) q[4];
u3(1.74729299765933,-1.81358512948542,-0.822592193171042) q[0];
cx q[0],q[4];
u1(1.52201365489869) q[4];
u3(0.550432948128745,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.16827944944046,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.01751652645286,-0.877853404205616,-2.29245801498837) q[4];
u3(0.534629147544645,1.24361839299513,-3.61939309218229) q[0];
u3(0.479820055449167,-1.84145501558003,2.22791842387922) q[4];
u3(1.14907003062116,0.338355909872563,-1.14705894813426) q[7];
cx q[7],q[4];
u1(2.93072368323925) q[4];
u3(-1.92690298856909,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.972408649148491,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.06794243424654,1.54261307877256,-2.68896076890295) q[4];
u3(1.99895913190474,-1.05738291036404,0.891689633123009) q[7];
u3(2.33605258715338,-1.27009219779054,4.25941584238670) q[9];
u3(1.70168643420067,1.52682033250856,1.58237625130969) q[8];
cx q[8],q[9];
u1(3.46421021061573) q[9];
u3(-0.674574059543767,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.84810813185497,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.24584093079113,3.09483544168804,0.0731559719045884) q[9];
u3(0.792241651072622,3.41756484181300,1.46197335343677) q[8];
u3(0.855549389699550,3.98980785799158,-1.82855139596020) q[6];
u3(1.51761319736480,2.02117950252708,-2.29295671947761) q[0];
cx q[0],q[6];
u1(2.04750767045518) q[6];
u3(-2.78067974126367,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.217196916961054,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.73720266778300,2.07972247464256,-3.91257202647265) q[6];
u3(1.00652289279059,-2.11552834671927,1.62486625091507) q[0];
u3(2.50621313577805,-2.42768023765029,1.27969999013407) q[12];
u3(2.19210650687996,-1.46474713358787,0.0805541429192094) q[11];
cx q[11],q[12];
u1(1.84993694537155) q[12];
u3(-2.86410837978451,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.788647926257609,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.39673596258903,-3.18678079689685,1.76452291600963) q[12];
u3(1.93111299537322,-2.42026470501475,0.951556783512861) q[11];
u3(0.896734526365027,2.13443520679120,-3.22237042984850) q[3];
u3(1.63475306093800,-2.47866641937285,3.16305277383681) q[10];
cx q[10],q[3];
u1(0.734473210446206) q[3];
u3(-3.42579652519654,0.0,0.0) q[10];
cx q[3],q[10];
u3(1.46403064602944,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.503409362945654,-3.62621630890949,2.21152458013899) q[3];
u3(2.82315196311387,1.75891806077810,1.40214319070126) q[10];
u3(2.37905754239333,0.817167116594892,-3.29943578600184) q[2];
u3(1.89418949356130,-2.67447392101200,3.46831803884846) q[5];
cx q[5],q[2];
u1(0.829999804337030) q[2];
u3(-3.41200652381675,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.96526966688327,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.764099539682443,-2.52772063764766,2.13580098489458) q[2];
u3(2.62922037395573,-0.844552850198477,3.22752221473221) q[5];
u3(0.994275264774238,1.02974846159028,1.08822823988409) q[1];
u3(1.18575497051356,-0.274621034170347,-2.36890779523267) q[7];
cx q[7],q[1];
u1(-0.0159330645803755) q[1];
u3(-0.944387120650013,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.74055103712277,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.00548102636409,1.48540099360207,1.14007619985702) q[1];
u3(2.75367540863579,0.139827162899001,3.44748114946942) q[7];
u3(2.27039460218850,0.853721253695355,1.00209060473333) q[5];
u3(1.42424781798425,-5.56216282820664,0.354152566389681) q[3];
cx q[3],q[5];
u1(1.94710757244586) q[5];
u3(-2.75387568670793,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.637909202986560,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.28025646728858,3.38095965535733,-1.98196316725778) q[5];
u3(0.575824554873283,-2.54006179080299,-3.42041965328936) q[3];
u3(1.62559684203880,2.13935854084689,-3.37623998559402) q[0];
u3(1.94592049713358,2.57099419526214,-3.05376849804505) q[9];
cx q[9],q[0];
u1(2.97596016429917) q[0];
u3(-2.04817411470350,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.625992773249112,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.22026325103569,0.174406816070944,-3.06376013387815) q[0];
u3(1.11336977911268,3.77461361096067,-0.577774003404542) q[9];
u3(2.12716217877327,0.0956304888533031,-2.12613604516133) q[10];
u3(1.75639311588890,0.662711470429012,-3.77898345250850) q[2];
cx q[2],q[10];
u1(1.46737421919987) q[10];
u3(-0.956782787787027,0.0,0.0) q[2];
cx q[10],q[2];
u3(0.0908319744346389,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.930434552994723,2.95542225611624,0.341930070836280) q[10];
u3(0.897993370085240,0.218927129177879,-5.01256980717863) q[2];
u3(1.38025108521152,1.24883663266107,-1.63803782140160) q[8];
u3(1.63628661877613,-4.98715600732518,0.283251941268302) q[11];
cx q[11],q[8];
u1(2.27365425576803) q[8];
u3(-2.76243572726256,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.27936849428956,0.0,0.0) q[11];
cx q[11],q[8];
u3(0.597770244972558,2.00426696900568,-3.96476688099020) q[8];
u3(0.335733143545845,-0.475259665926605,5.45343488965350) q[11];
u3(1.70152377937138,-0.546485692709268,2.06534182193024) q[4];
u3(2.92131619963264,0.908430050568048,3.53770524912559) q[12];
cx q[12],q[4];
u1(2.34762296834967) q[4];
u3(-2.65127487556963,0.0,0.0) q[12];
cx q[4],q[12];
u3(1.13236856677270,0.0,0.0) q[12];
cx q[12],q[4];
u3(2.31306555310056,-1.16338220445984,1.27667304285186) q[4];
u3(1.03732415585837,-0.387865935922507,5.75298585146505) q[12];
u3(1.17770686267088,2.16233259927748,-3.19984257081808) q[10];
u3(1.35662394885422,-2.67277574185124,3.09272419490569) q[5];
cx q[5],q[10];
u1(0.104406488518794) q[10];
u3(-1.16396282994970,0.0,0.0) q[5];
cx q[10],q[5];
u3(0.799386227834733,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.955088626310305,2.91980573347745,1.04690932387833) q[10];
u3(2.51260178437476,-4.73324974298769,-1.32457673930423) q[5];
u3(1.25369965033437,-0.605041435556619,-1.29999095210770) q[0];
u3(1.00983605205927,-4.39209239555627,1.05872042154660) q[4];
cx q[4],q[0];
u1(2.82660672409290) q[0];
u3(-1.86140086925876,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.0424457324530678,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.42162154817199,0.972204559830421,-1.04816260424326) q[0];
u3(2.38309273468526,2.18980117681043,0.554189648098339) q[4];
u3(0.819201278313093,-2.51229925609730,1.69689784583903) q[9];
u3(0.218747412155113,-2.31021922759520,0.628956239373191) q[6];
cx q[6],q[9];
u1(3.72164248577705) q[9];
u3(-4.41102105935204,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.701486520391027,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.37573090700006,0.180453033046524,0.156914123707612) q[9];
u3(1.57835004811543,-2.38495786179493,2.52458547193742) q[6];
u3(1.66023850260523,-0.125631327668997,2.04618380030494) q[1];
u3(2.07026985307648,-2.63548736991390,-1.23605805297730) q[12];
cx q[12],q[1];
u1(2.75359243369265) q[1];
u3(-3.05740171905751,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.61090193732720,0.0,0.0) q[12];
cx q[12],q[1];
u3(0.594691338024389,1.33222959850920,-2.04694927829417) q[1];
u3(2.20295874435054,-2.22242822667402,1.27479429248113) q[12];
u3(1.53894122551313,1.64811549626316,0.248510985036832) q[11];
u3(0.680950205932844,0.783024830769219,-4.72288912305594) q[3];
cx q[3],q[11];
u1(0.835203531379231) q[11];
u3(-1.10289125676736,0.0,0.0) q[3];
cx q[11],q[3];
u3(0.104125025261282,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.24006475994765,1.58457817217644,2.05904879616957) q[11];
u3(0.732133512002309,0.819276955521740,3.89812419683109) q[3];
u3(1.19900689624350,-1.47096654370376,0.709301741862479) q[2];
u3(1.80007057525438,-1.96549781494643,1.08344750130385) q[7];
cx q[7],q[2];
u1(1.67291371130245) q[2];
u3(-2.72603112611183,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.351189528968921,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.965139097962579,-0.423557230893299,2.38024082565426) q[2];
u3(1.69675057578032,0.405820129905946,1.52447743838105) q[7];
u3(1.31040770239429,0.189664488191661,1.02693624774215) q[7];
u3(1.76494605387533,-0.889672914095452,-1.33014278588117) q[4];
cx q[4],q[7];
u1(1.76706938177661) q[7];
u3(-2.69212621466494,0.0,0.0) q[4];
cx q[7],q[4];
u3(3.16212265681027,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.614944279099468,1.28500351735870,-0.0851625254893180) q[7];
u3(2.65569083094848,5.95608467024407,-0.291228080982441) q[4];
u3(1.31730080723561,3.08980096971083,-1.87676178365157) q[8];
u3(1.05092274990196,1.39331834067824,-2.40525334434586) q[9];
cx q[9],q[8];
u1(0.508856918749119) q[8];
u3(-1.15864340604278,0.0,0.0) q[9];
cx q[8],q[9];
u3(3.17436342862304,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.53944421338794,-0.620194207301184,-0.180067320685514) q[8];
u3(1.23005013688241,-4.04980982509486,-0.596585228144884) q[9];
u3(1.19253978578334,0.585638495077954,0.949790358409156) q[11];
u3(1.42307603467585,-2.11661737955852,-1.93241583182050) q[3];
cx q[3],q[11];
u1(0.00549816596074026) q[11];
u3(-2.29655714239140,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.21988560941507,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.82919408013580,0.0596424672397483,0.301064637988101) q[11];
u3(1.74524640520230,0.421095638327186,3.00740001074384) q[3];
u3(1.79501613265593,-1.48110702702273,0.749876260225174) q[1];
u3(1.52074537397943,-2.50192126188385,-0.343295548498251) q[12];
cx q[12],q[1];
u1(1.89669297895537) q[1];
u3(-2.99091935608023,0.0,0.0) q[12];
cx q[1],q[12];
u3(0.504505588073808,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.21560389837010,0.376223007685175,-2.95689437152029) q[1];
u3(2.38934939155478,-0.119609473084202,-4.89970174691649) q[12];
u3(1.56831046016675,0.532396415429975,-3.47209528428176) q[5];
u3(0.836977389467577,-0.907406206398501,4.96905789884829) q[6];
cx q[6],q[5];
u1(0.136696352840903) q[5];
u3(-1.41756656928014,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.67396703740354,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.17289384905419,0.420184945306218,1.83134086691221) q[5];
u3(2.12447633103047,3.44765270103688,-2.04678403619929) q[6];
u3(0.556907549663082,1.34775259836585,-1.33032809079312) q[2];
u3(0.636795865648494,0.110303402311964,-0.397251853564937) q[0];
cx q[0],q[2];
u1(1.25432513770886) q[2];
u3(-0.536614394591681,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.86640918610853,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.784603531227802,0.0912171786311803,-0.439416374332034) q[2];
u3(0.991024488063022,0.853159711959536,4.88515555458128) q[0];
u3(1.26052691345378,-0.180659614323723,0.893232847603558) q[9];
u3(1.15765643428482,-2.09329890990583,-1.60819967285748) q[0];
cx q[0],q[9];
u1(2.53198786480962) q[9];
u3(-2.99529533560183,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.81907853101531,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.908427518828771,1.27408844079788,1.34579365851682) q[9];
u3(2.00094578064912,0.313909359496470,2.67052458944246) q[0];
u3(1.12996161257822,1.14570322515721,1.11333765266804) q[2];
u3(1.55191347292187,-1.56256079456966,-0.699890442759533) q[1];
cx q[1],q[2];
u1(2.48734786716416) q[2];
u3(-1.57868233614650,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.12363439866284,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.453885479191635,-3.18427728713719,0.205824603875150) q[2];
u3(2.55086095898357,-0.128702862251539,-1.12654323661829) q[1];
u3(1.66527117359151,2.86562601338452,0.187016829997422) q[6];
u3(1.06400529784717,-0.205977642552927,-2.27874559876186) q[10];
cx q[10],q[6];
u1(0.966584139409155) q[6];
u3(-0.593628340679468,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.62585064603143,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.203311281277386,2.26357858727817,-2.69583025356798) q[6];
u3(1.53835871975079,-5.01310446763437,0.760630967085764) q[10];
u3(2.51753817670488,-1.68104580451604,3.65712524401587) q[3];
u3(1.65052482931926,1.24202151574044,1.18033543247985) q[5];
cx q[5],q[3];
u1(-0.178244477488588) q[3];
u3(-1.65473902600964,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.427896768398098,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.10982492844479,0.409550413708971,2.51303438631495) q[3];
u3(1.55019558327755,0.239669435108134,5.49940837579276) q[5];
u3(0.903712522787983,0.0688529429835391,1.51911312931182) q[12];
u3(1.74040112573597,-2.87926757298377,-1.06423723752259) q[11];
cx q[11],q[12];
u1(1.64591656240062) q[12];
u3(-2.20143074351754,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.433817251601880,0.0,0.0) q[11];
cx q[11],q[12];
u3(0.905988826728011,-0.831677434465631,0.489935771240213) q[12];
u3(2.23674342506911,-2.56224200876812,-3.29340528245589) q[11];
u3(1.74624453659996,-1.14914823484374,0.00227453873465322) q[7];
u3(1.81550190246152,-2.18546406993880,-0.257744499160401) q[8];
cx q[8],q[7];
u1(2.55469605103093) q[7];
u3(0.0959878682305029,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.39527550669612,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.85167378096102,-2.44725579817647,-1.23496334832266) q[7];
u3(1.22658338973875,-0.684634813535842,2.48380101809463) q[8];
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