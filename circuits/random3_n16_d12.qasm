OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(2.37613094997230,-2.47358084289731,-0.0491070636496445) q[15];
u3(2.03340343004189,-3.94725428104620,-0.870588208611913) q[2];
cx q[2],q[15];
u1(2.45317679537390) q[15];
u3(-1.82065833324773,0.0,0.0) q[2];
cx q[15],q[2];
u3(0.274561612928055,0.0,0.0) q[2];
cx q[2],q[15];
u3(1.03685219242113,-2.04202739074524,2.03957814597473) q[15];
u3(2.81190794927652,-1.19181310017798,-1.38340502151484) q[2];
u3(1.46257148367820,-0.368802631365395,0.280871618329640) q[11];
u3(1.19808417735605,-2.78328606254308,-0.655919530175473) q[3];
cx q[3],q[11];
u1(3.18783625190226) q[11];
u3(-0.700160862191518,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.91162283417007,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.82689491174513,0.828827027905824,-2.65896654093843) q[11];
u3(0.944237698488225,-1.88087448985654,2.99011373408159) q[3];
u3(1.91271118288810,-3.37170345148617,2.87498533336625) q[9];
u3(1.30746846882205,3.13563628936814,-2.66104879348583) q[14];
cx q[14],q[9];
u1(2.15013922005849) q[9];
u3(0.208462049540414,0.0,0.0) q[14];
cx q[9],q[14];
u3(1.32466580583222,0.0,0.0) q[14];
cx q[14],q[9];
u3(1.17663760641949,1.19721395283008,-4.40267953739927) q[9];
u3(2.47298006695390,3.11565083909848,-0.811556807590420) q[14];
u3(2.40375579941815,-1.51084851549492,4.19367660625227) q[1];
u3(0.942723259015449,1.22046808092569,0.724081267558510) q[4];
cx q[4],q[1];
u1(1.97970794057890) q[1];
u3(-2.32753072763278,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.489286738609213,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.791882091926162,2.74412677994246,-2.74203934065195) q[1];
u3(1.69521172629280,-1.88524396155113,-0.186097896561644) q[4];
u3(2.38941480507382,0.0955095483864186,2.51090222325672) q[0];
u3(1.61248959824730,-1.90426430613108,-1.45519100238447) q[5];
cx q[5],q[0];
u1(0.0609483543223146) q[0];
u3(-2.27081645014266,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.40286702012992,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.27146504758851,0.0855634996868609,-1.43966200757951) q[0];
u3(0.265543760829283,-0.490408334430554,5.40066081554968) q[5];
u3(1.96898464419670,-1.74716486909660,-0.838282480876027) q[7];
u3(0.877229358318118,-3.63489937537888,0.305672925026693) q[13];
cx q[13],q[7];
u1(2.76615836580140) q[7];
u3(-1.54375059092819,0.0,0.0) q[13];
cx q[7],q[13];
u3(0.978619602069058,0.0,0.0) q[13];
cx q[13],q[7];
u3(2.12438231091768,3.47434419719512,0.309254391211948) q[7];
u3(0.196718154577320,0.604150668179451,-1.51868104555705) q[13];
u3(2.38461394055832,-3.79095864591630,1.85347237218130) q[8];
u3(0.0232345234780548,-1.94774583639351,3.67006640831649) q[12];
cx q[12],q[8];
u1(3.44289043239092) q[8];
u3(-1.34079374705881,0.0,0.0) q[12];
cx q[8],q[12];
u3(2.38938251272960,0.0,0.0) q[12];
cx q[12],q[8];
u3(0.987427091778511,-1.37213175771563,4.47256698507214) q[8];
u3(1.00265099525637,0.188377055386500,5.12748942224210) q[12];
u3(1.96398050669598,-3.81300704806361,0.865969440122633) q[10];
u3(1.14941084578963,0.0805693703597592,3.25260022396361) q[6];
cx q[6],q[10];
u1(1.89231649122512) q[10];
u3(-2.53246591708486,0.0,0.0) q[6];
cx q[10],q[6];
u3(3.37172179005727,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.16971499609540,2.16924756729414,-1.96098183041389) q[10];
u3(2.51982582119305,2.44200885291660,-0.328981834721559) q[6];
u3(1.13119981914707,-2.45596989322753,0.503307685175818) q[7];
u3(1.41235102048316,-3.60834151395013,-0.615710915462099) q[9];
cx q[9],q[7];
u1(1.22024723190450) q[7];
u3(-3.27820207641380,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.44913286858303,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.37787679368260,-0.489746620950321,4.50045049140121) q[7];
u3(1.66125259461805,-3.42661018841263,2.02357228988016) q[9];
u3(2.34410615201993,4.49628313182451,-1.59732555793714) q[14];
u3(0.0371583616542635,2.91101356583708,-2.20203467335383) q[13];
cx q[13],q[14];
u1(3.13786572642253) q[14];
u3(-0.901509071083387,0.0,0.0) q[13];
cx q[14],q[13];
u3(1.73241544520654,0.0,0.0) q[13];
cx q[13],q[14];
u3(1.94153848971552,-1.82057792686619,3.36928351579826) q[14];
u3(1.10033346030818,-2.99222280069290,2.03034511329187) q[13];
u3(2.60085904301395,1.84357076810819,-2.54025817209426) q[15];
u3(1.70634666076340,2.27753946106376,-2.66185443183506) q[12];
cx q[12],q[15];
u1(0.649736126279630) q[15];
u3(-0.995669140816503,0.0,0.0) q[12];
cx q[15],q[12];
u3(3.08924390518058,0.0,0.0) q[12];
cx q[12],q[15];
u3(1.86836503463380,-2.23157698107180,-0.840897476146953) q[15];
u3(1.58624941504977,0.518543882883003,0.886566319607261) q[12];
u3(0.321275130121597,-0.463837797263651,-0.455792986669788) q[6];
u3(1.19869421923587,-3.77110425674439,1.65984401737778) q[11];
cx q[11],q[6];
u1(2.53510278169978) q[6];
u3(-2.07254912118730,0.0,0.0) q[11];
cx q[6],q[11];
u3(0.230587911621151,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.75326576444631,-3.66918525893206,1.90796397004055) q[6];
u3(1.88344839030427,0.300704566888901,0.175033198361141) q[11];
u3(1.88076244873131,1.61539683491195,-2.83342425657462) q[5];
u3(2.19069719205262,-3.01079729517751,2.84118491399413) q[2];
cx q[2],q[5];
u1(2.08591636509982) q[5];
u3(-3.08672360629339,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.32693872142211,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.59748435905916,-1.19995431530204,2.47747752238058) q[5];
u3(0.103163969524823,-3.00564441456625,2.67003562797473) q[2];
u3(1.42318116068670,1.37265876326831,-3.07404974487654) q[0];
u3(2.17193048970552,1.96577278424506,-3.53424911396638) q[10];
cx q[10],q[0];
u1(2.50836437638428) q[0];
u3(-1.88762132366207,0.0,0.0) q[10];
cx q[0],q[10];
u3(-0.120840147169551,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.19791923442804,-0.618791541475661,-1.70855601419393) q[0];
u3(2.80547127196291,-0.455192857510184,-5.26347541136077) q[10];
u3(0.954655403505892,0.781561357868020,-3.92195439267691) q[1];
u3(1.57926247204354,-1.10925443691257,4.75656354666734) q[4];
cx q[4],q[1];
u1(1.78127950468727) q[1];
u3(-2.20799539015074,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.52383052107691,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.19879028108603,-2.79916848473177,1.47317360920699) q[1];
u3(1.65867789705993,-1.24380511722502,4.50447867951208) q[4];
u3(2.92893891026886,-0.990060480785262,-1.12889751083119) q[3];
u3(1.23888662861993,0.612775609910881,-4.71533881903238) q[8];
cx q[8],q[3];
u1(1.10911121250484) q[3];
u3(-3.13108926063381,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.46760045726505,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.33145316739277,-2.35804958872040,2.30735392679716) q[3];
u3(0.184878360718447,-3.42222765134747,0.508370953549770) q[8];
u3(1.59345824576637,1.93324253923032,-4.34862694703576) q[12];
u3(0.572422874720187,-1.58574317966618,3.00795597021575) q[0];
cx q[0],q[12];
u1(3.08986459603377) q[12];
u3(-2.59788722567824,0.0,0.0) q[0];
cx q[12],q[0];
u3(0.822440061518462,0.0,0.0) q[0];
cx q[0],q[12];
u3(2.17748586529202,-0.181943129259684,-0.972786675937389) q[12];
u3(2.35548183955005,-1.75081680655701,4.30621015006134) q[0];
u3(2.47953729469232,0.552546390606084,-1.78858064049953) q[1];
u3(2.30766637961182,4.83721148175447,0.556712997330689) q[9];
cx q[9],q[1];
u1(0.319368518467896) q[1];
u3(-1.53351139306790,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.00207295827516,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.85852183208965,0.0338566744933957,-2.16310105013351) q[1];
u3(1.51167017998223,-0.382196196922564,-4.33627326757805) q[9];
u3(1.90291205609000,3.19698305421570,-2.37687094516718) q[15];
u3(0.378421066969900,1.50339103819417,-0.220808022781721) q[6];
cx q[6],q[15];
u1(1.28112797399495) q[15];
u3(-1.09121529202484,0.0,0.0) q[6];
cx q[15],q[6];
u3(-0.528494821613932,0.0,0.0) q[6];
cx q[6],q[15];
u3(1.30420704214244,-1.19380644652271,0.264093476375142) q[15];
u3(0.885508979737157,-2.09659871232394,-2.82897192705707) q[6];
u3(1.94337830584873,-2.58512135817357,0.687814487997784) q[11];
u3(1.72366100415697,-3.26513224175393,-0.503144691670401) q[2];
cx q[2],q[11];
u1(0.620170767594894) q[11];
u3(-1.33032583510677,0.0,0.0) q[2];
cx q[11],q[2];
u3(-0.252371965730113,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.937071376596354,-1.46404263285479,3.77777625301144) q[11];
u3(2.46607239396127,0.510324730419375,5.24025877439627) q[2];
u3(1.79178376834679,-1.63750270838945,-1.34456253186744) q[4];
u3(0.343613103639623,-5.19015824181054,0.916180627402096) q[8];
cx q[8],q[4];
u1(1.14525667359348) q[4];
u3(-3.69307041636463,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.54950319201814,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.04880336283906,-2.76261445164069,0.197913609265838) q[4];
u3(1.17630374470416,0.722879849067148,-3.30882501163726) q[8];
u3(2.02416097313247,-3.53195432059826,2.57610742479622) q[13];
u3(2.25544297849358,2.89949582326185,-3.14576406401249) q[10];
cx q[10],q[13];
u1(0.540514784205587) q[13];
u3(-1.46582489837971,0.0,0.0) q[10];
cx q[13],q[10];
u3(2.08582117697344,0.0,0.0) q[10];
cx q[10],q[13];
u3(2.38869544298070,0.392241133614700,-2.39178383980225) q[13];
u3(2.16361001296663,-0.556052447368125,-5.42827689812239) q[10];
u3(1.59850358793714,2.18655025277296,-2.30269897026641) q[14];
u3(0.279423978578162,-3.05577693020434,2.85846155075958) q[5];
cx q[5],q[14];
u1(1.82774254334253) q[14];
u3(-2.17524828854136,0.0,0.0) q[5];
cx q[14],q[5];
u3(-0.332469315383461,0.0,0.0) q[5];
cx q[5],q[14];
u3(1.49513498007458,-1.50867892450024,1.72771698254135) q[14];
u3(0.735352325705226,-0.0426975115391661,-2.89477734610097) q[5];
u3(2.01379568040771,0.197689066139911,-1.55879909204735) q[3];
u3(1.14421968917325,0.801287846760222,-4.43399599313062) q[7];
cx q[7],q[3];
u1(1.28860885657051) q[3];
u3(-0.808231949987727,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.0374326696230576,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.65030102753344,-0.00955188447881172,-1.14740590935312) q[3];
u3(1.25574548512883,-0.189053166817656,3.00441871872094) q[7];
u3(2.22609130913389,1.88350365846867,-0.838320726815409) q[15];
u3(1.43524933891099,-0.0752562990260834,-2.63368625021946) q[13];
cx q[13],q[15];
u1(3.20648756544839) q[15];
u3(-0.994561804472093,0.0,0.0) q[13];
cx q[15],q[13];
u3(1.79092971483641,0.0,0.0) q[13];
cx q[13],q[15];
u3(1.56524936222327,2.11990661927310,-2.01158636642464) q[15];
u3(1.21977281644302,-0.168252475614022,-3.79142239355736) q[13];
u3(2.67909573480704,-0.362517268148422,-1.17077039807538) q[9];
u3(1.39503514414802,-4.06524045446988,0.955957710055763) q[11];
cx q[11],q[9];
u1(0.303413958782917) q[9];
u3(-1.25060011245366,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.18835168756820,0.0,0.0) q[11];
cx q[11],q[9];
u3(3.04777893830758,-4.65360300862855,0.629610483362130) q[9];
u3(1.46981297512348,-1.35651029778385,-2.49560278073443) q[11];
u3(1.82823995622937,-0.473597756656135,1.33465902139006) q[6];
u3(1.48952994096782,-0.638439731213866,-0.911361009780675) q[8];
cx q[8],q[6];
u1(0.309512427998646) q[6];
u3(-1.79626465252351,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.71002177663385,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.25852689934504,3.34110927699972,-2.05791971333467) q[6];
u3(1.61685292035847,0.371618319739954,-4.37535317513773) q[8];
u3(2.53746600285865,0.570058233778103,-1.57989998406169) q[12];
u3(1.68954425611927,1.88994099490709,-4.22391735893846) q[2];
cx q[2],q[12];
u1(2.60361266927957) q[12];
u3(-2.02476400397013,0.0,0.0) q[2];
cx q[12],q[2];
u3(3.22178883038763,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.48190458869829,2.88070822726546,0.294652266768084) q[12];
u3(1.36872959722728,-3.41348943383962,-0.812110258229664) q[2];
u3(0.723150779974453,0.355533182690006,-0.717960023524262) q[10];
u3(0.806098822194065,-2.73695732622495,1.35460487366084) q[4];
cx q[4],q[10];
u1(2.05712183396991) q[10];
u3(-3.19504091456498,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.68953703503114,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.953913801963038,2.77367203668774,-3.06956782435506) q[10];
u3(2.17349380244436,-0.702695564464322,2.18951867575406) q[4];
u3(1.96419407324917,-0.162757163458807,-0.991246593444142) q[5];
u3(1.47838038970934,-4.09232298848880,0.909594616442515) q[1];
cx q[1],q[5];
u1(2.47871123871423) q[5];
u3(-2.38490256076503,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.851128020361487,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.463858709919396,-0.677105324560479,1.88901881855715) q[5];
u3(1.14740741448430,-2.59445579570954,-1.16651494473435) q[1];
u3(2.20423647720818,-1.74376458399559,-0.765356002339697) q[7];
u3(1.67825864558816,-3.76912671929111,-0.130585722593529) q[14];
cx q[14],q[7];
u1(0.940448190262017) q[7];
u3(-0.603041139394486,0.0,0.0) q[14];
cx q[7],q[14];
u3(1.98720155549302,0.0,0.0) q[14];
cx q[14],q[7];
u3(1.57728234325341,1.13035058670624,-2.30202432425425) q[7];
u3(1.14417399511769,1.08827990594700,3.45279509018625) q[14];
u3(0.973326383992139,1.90524641809168,-1.13395785245387) q[0];
u3(1.13527039604535,1.33304017986903,-1.65711008600391) q[3];
cx q[3],q[0];
u1(-0.384499849295280) q[0];
u3(1.02501795570788,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.93034026586098,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.08183282824467,1.53208782017071,-1.73475729167844) q[0];
u3(1.38977492340906,0.670705928569292,1.39889716444822) q[3];
u3(1.80332665989467,-2.62412727614102,0.0382013022995744) q[6];
u3(0.900953598236490,-2.92637948183975,0.369384561761062) q[2];
cx q[2],q[6];
u1(1.04423720589894) q[6];
u3(-3.70281367424631,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.99110815586145,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.498014108232555,1.19490128651236,-4.07687448151942) q[6];
u3(2.18786432871876,0.597190033257815,-2.31083392880222) q[2];
u3(1.44064007226196,1.75827245574396,-0.412514144282003) q[0];
u3(2.66369938002222,0.108023356491971,-2.70053400531425) q[9];
cx q[9],q[0];
u1(2.94239256146178) q[0];
u3(-1.41111664056482,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.59860197592638,0.0,0.0) q[9];
cx q[9],q[0];
u3(0.968199819801119,4.12168596566968,-0.821772654213957) q[0];
u3(1.59694021491886,2.79577130054194,-3.21297495489268) q[9];
u3(0.754108820951758,-2.83858519682434,-0.234440026532779) q[14];
u3(0.876314029178943,-2.19128606465510,-0.921444265181268) q[4];
cx q[4],q[14];
u1(1.80559814610212) q[14];
u3(-2.68335240738640,0.0,0.0) q[4];
cx q[14],q[4];
u3(0.965884425669417,0.0,0.0) q[4];
cx q[4],q[14];
u3(1.46746099994312,-0.284545908539718,-0.797006124753047) q[14];
u3(2.82653085943835,-0.195299779512438,-1.85183469818527) q[4];
u3(2.37348124670027,1.68544976620298,-0.0926312255745825) q[1];
u3(1.73350610494513,0.187490057431335,-1.85878851055434) q[12];
cx q[12],q[1];
u1(2.26238275343387) q[1];
u3(0.436238867838796,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.18908055862934,0.0,0.0) q[12];
cx q[12],q[1];
u3(0.921259913443855,3.95117942383993,-1.14872918253964) q[1];
u3(2.05050330118951,-5.16942018830539,-0.801788573474096) q[12];
u3(0.515636233803445,0.129940373060320,-0.476365255046000) q[10];
u3(0.504740904743028,-2.38740138450872,0.758032783843757) q[11];
cx q[11],q[10];
u1(0.814523944137585) q[10];
u3(-3.20104984447835,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.55125752355866,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.41006678824906,-3.13127069864402,-0.896406146091047) q[10];
u3(1.58482582912251,-0.246975884460218,-0.964946828275828) q[11];
u3(1.69927922682381,0.217381009309823,1.06427440278202) q[8];
u3(1.63609167519388,-1.23821251255430,-1.40497836712771) q[7];
cx q[7],q[8];
u1(3.59218271380366) q[8];
u3(-1.04665357411567,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.45019541311510,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.17763971666380,-1.44101447592314,-0.367388283466743) q[8];
u3(0.828967403397636,1.18913038522522,-0.405839416070014) q[7];
u3(0.676421973429739,-1.86831175850131,0.937416792237993) q[5];
u3(0.591462897876276,-1.85596095424991,-0.0145123568051907) q[3];
cx q[3],q[5];
u1(0.297462686004764) q[5];
u3(-1.45121292106622,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.64594302826838,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.61753553074612,-2.21099185123096,1.98908044161362) q[5];
u3(1.27440056349876,2.82900311775229,-0.189302658321794) q[3];
u3(0.605004716274528,-1.16655803975164,0.937346303511904) q[15];
u3(1.19679104473408,-1.00346185928787,-1.30812117013545) q[13];
cx q[13],q[15];
u1(1.50456039907437) q[15];
u3(-0.944647142607542,0.0,0.0) q[13];
cx q[15],q[13];
u3(0.233861687157628,0.0,0.0) q[13];
cx q[13],q[15];
u3(1.28009013537134,0.235300831371942,0.749761177209630) q[15];
u3(2.72206330875883,1.33299300229712,4.60753893973903) q[13];
u3(1.05127193574370,1.48150830968857,-3.84458340168579) q[5];
u3(1.86999780941380,3.09196767286804,-2.41008811866450) q[15];
cx q[15],q[5];
u1(0.686766683560052) q[5];
u3(-3.28749693700873,0.0,0.0) q[15];
cx q[5],q[15];
u3(1.21806595542194,0.0,0.0) q[15];
cx q[15],q[5];
u3(2.37919159065502,-2.94481835222867,3.24742328363618) q[5];
u3(0.206883808602203,4.03860875985259,0.327631088937643) q[15];
u3(1.62717634163342,0.450317864249636,1.68059684159291) q[3];
u3(1.49671200427148,-1.53722131018072,-0.193532369858445) q[9];
cx q[9],q[3];
u1(1.45662181955815) q[3];
u3(0.333415637804190,0.0,0.0) q[9];
cx q[3],q[9];
u3(0.867531139254984,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.36967003355137,-0.476039876546939,2.38136227028364) q[3];
u3(0.858379380183186,2.98321795029302,2.31639775440008) q[9];
u3(0.924585445580843,1.26641892793106,-3.46476733818586) q[14];
u3(1.62132614065280,2.32331797045621,-3.55401116065482) q[0];
cx q[0],q[14];
u1(2.88225072381256) q[14];
u3(-1.80792723058711,0.0,0.0) q[0];
cx q[14],q[0];
u3(0.750440771945565,0.0,0.0) q[0];
cx q[0],q[14];
u3(1.81772892677553,-0.751899732676007,-1.31972473667476) q[14];
u3(2.55370486822547,1.05275688184308,-4.52330315302906) q[0];
u3(3.05896111029110,0.758593457308520,-2.58702283760467) q[2];
u3(2.68078006606487,1.69468901256669,-1.79089056439087) q[4];
cx q[4],q[2];
u1(1.79230910745892) q[2];
u3(-3.06530391403401,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.50869625900861,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.347639687864453,-3.50830286689561,1.16258072792463) q[2];
u3(2.33287674638832,-5.15281224690247,0.826443264364934) q[4];
u3(1.53464955876030,0.649165259203429,1.57867528124075) q[11];
u3(1.67125149757172,-0.733343603350666,-2.76036306077353) q[12];
cx q[12],q[11];
u1(1.99232271945266) q[11];
u3(-2.95629019452430,0.0,0.0) q[12];
cx q[11],q[12];
u3(0.901609807040205,0.0,0.0) q[12];
cx q[12],q[11];
u3(0.204229941427519,-1.65191253299128,3.72379940622477) q[11];
u3(2.09989063047844,1.37857248529275,-1.81916826570585) q[12];
u3(2.62790353795105,-3.60668059377329,1.24643634960495) q[13];
u3(1.88563955253012,-0.332171799416833,3.25063995208821) q[1];
cx q[1],q[13];
u1(1.64198528290761) q[13];
u3(-3.59623482208927,0.0,0.0) q[1];
cx q[13],q[1];
u3(2.05705269078852,0.0,0.0) q[1];
cx q[1],q[13];
u3(2.63219904482550,2.74344674855188,-0.524817232761405) q[13];
u3(2.03618873641201,-0.664949098186039,1.42556823883987) q[1];
u3(2.16970545689834,-0.877789142887734,-1.37500789896815) q[6];
u3(0.538901435445653,-0.281933820679351,-4.16500401601566) q[7];
cx q[7],q[6];
u1(0.711107951447778) q[6];
u3(-1.09915135641184,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.08850819960953,0.0,0.0) q[7];
cx q[7],q[6];
u3(3.11795638388067,-1.70887997243113,-0.150727619938108) q[6];
u3(2.31514551520592,1.14285138247621,1.00034419799190) q[7];
u3(1.99564551771978,-0.0224693185157211,1.40259643089473) q[10];
u3(1.24802371760421,-1.21448043033667,-2.92555374986116) q[8];
cx q[8],q[10];
u1(1.49052026095869) q[10];
u3(-2.67288878652993,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.348118030719832,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.24478602079530,-2.41432641993470,2.44530165376655) q[10];
u3(0.438481932037752,0.392935843660278,-1.76127184239782) q[8];
u3(0.993019405798932,1.70204359037855,-3.19153098463414) q[6];
u3(1.86293178635342,1.84872770306587,-3.91903807677104) q[12];
cx q[12],q[6];
u1(0.213174981745739) q[6];
u3(-1.24082763663683,0.0,0.0) q[12];
cx q[6],q[12];
u3(2.20671755475301,0.0,0.0) q[12];
cx q[12],q[6];
u3(0.949902470167894,-2.58430296373594,0.0954044704456873) q[6];
u3(1.66279971061593,2.17697536954953,2.82388807661460) q[12];
u3(1.38332220468902,3.10591397676613,-0.963759142848334) q[14];
u3(1.87604626569445,1.80662495579277,-1.27590437680595) q[10];
cx q[10],q[14];
u1(0.249398059354076) q[14];
u3(-1.06959521838443,0.0,0.0) q[10];
cx q[14],q[10];
u3(1.64813838458146,0.0,0.0) q[10];
cx q[10],q[14];
u3(1.75270266438182,2.87711139801966,-1.13997838193433) q[14];
u3(0.892024163057202,-0.298892309732090,3.61748944580163) q[10];
u3(0.893775295432638,0.590559471341453,-2.55359979498101) q[2];
u3(1.21971508949908,2.37864689805220,-2.49394297584953) q[15];
cx q[15],q[2];
u1(0.797312866491328) q[2];
u3(-1.53056345986637,0.0,0.0) q[15];
cx q[2],q[15];
u3(-0.517748607797684,0.0,0.0) q[15];
cx q[15],q[2];
u3(2.43981420760292,2.70321097354763,-1.20136970327478) q[2];
u3(0.748540842681428,0.352419103500713,-4.61670772407589) q[15];
u3(1.05130624707667,0.154212550991137,1.93051326601510) q[4];
u3(1.27442877684183,-2.09072955078146,-1.72886403431482) q[11];
cx q[11],q[4];
u1(2.06262987737344) q[4];
u3(-2.80000398261126,0.0,0.0) q[11];
cx q[4],q[11];
u3(0.603986808759091,0.0,0.0) q[11];
cx q[11],q[4];
u3(2.37454366355472,-2.55403969271132,3.03103525083462) q[4];
u3(1.49084132112453,-2.84786722164804,0.794640129431344) q[11];
u3(0.841673009037480,1.79669845450326,-2.59185166734833) q[0];
u3(1.60561314011016,-3.39152146104535,2.19815137112505) q[7];
cx q[7],q[0];
u1(0.0298683190878357) q[0];
u3(1.00340959583945,0.0,0.0) q[7];
cx q[0],q[7];
u3(3.60288179716721,0.0,0.0) q[7];
cx q[7],q[0];
u3(3.03214029674288,0.778636595924471,-1.59092066683963) q[0];
u3(1.44378426229003,-0.762647474492511,2.21878881062768) q[7];
u3(2.18290510726896,-0.366411450208509,0.774090670804753) q[8];
u3(2.60025048842494,-0.679311043253053,-1.94704617709845) q[13];
cx q[13],q[8];
u1(0.617653067800433) q[8];
u3(-3.11888063608824,0.0,0.0) q[13];
cx q[8],q[13];
u3(2.17099171680770,0.0,0.0) q[13];
cx q[13],q[8];
u3(2.05837360444002,-2.89247390699934,2.42591351360100) q[8];
u3(0.548252693443112,2.28187359637117,0.994027486667801) q[13];
u3(1.59157454427505,-0.751062383651176,-1.54929880781882) q[9];
u3(2.78816915393799,0.805525362877726,-4.92309634744227) q[5];
cx q[5],q[9];
u1(0.528671196132119) q[9];
u3(-3.31544275211723,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.74085199750562,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.89903605195239,1.52487821718584,-0.0362176314725995) q[9];
u3(2.53017184357337,-2.07718283450559,-1.69493771743937) q[5];
u3(0.234060935042885,-1.11335602749643,1.30532360260143) q[1];
u3(0.526252705817666,-0.594758661237515,-1.58672059689824) q[3];
cx q[3],q[1];
u1(1.43404872174547) q[1];
u3(-0.614344176504410,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.323052903493134,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.60209554837262,-2.73123622564606,-0.618395256952440) q[1];
u3(1.26944802065464,0.819588579431203,-0.402512726068233) q[3];
u3(1.23848704359519,0.743911411713313,0.138646525847622) q[4];
u3(2.18264556433589,-1.23132035758138,-4.28901282855577) q[12];
cx q[12],q[4];
u1(4.17283330667489) q[4];
u3(-3.42647886284316,0.0,0.0) q[12];
cx q[4],q[12];
u3(-0.582925125294415,0.0,0.0) q[12];
cx q[12],q[4];
u3(1.29081896992794,-2.26954036614666,2.67884648954902) q[4];
u3(1.66387860872626,0.800208780736382,0.448533915832949) q[12];
u3(2.99917783431685,-0.439646318967474,2.13089463641289) q[3];
u3(2.70492515146586,1.83358896645568,3.13321504005281) q[1];
cx q[1],q[3];
u1(4.19358327864283) q[3];
u3(-3.05873616259113,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.400981725691906,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.80648542421849,-1.98672418017912,1.85618563893710) q[3];
u3(1.59952922534765,1.19356639566386,-0.494983205011666) q[1];
u3(1.61822786037638,-1.18873302725380,-0.0405493434395161) q[9];
u3(1.37887328935739,-3.79276127631124,-0.833960594618339) q[11];
cx q[11],q[9];
u1(3.80612267994951) q[9];
u3(-3.48603922500217,0.0,0.0) q[11];
cx q[9],q[11];
u3(-1.06903972339137,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.58138480223952,1.32536538834347,-1.44987693160475) q[9];
u3(2.54738901779157,-0.693900348650975,-1.05054641080371) q[11];
u3(0.746004135357403,0.0730407165625301,0.350738873878768) q[8];
u3(0.625505911875459,-0.973822354819997,0.101753619262156) q[2];
cx q[2],q[8];
u1(1.16277740335640) q[8];
u3(-0.783703534838785,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.14282532805532,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.51181607852107,-3.77841699171223,1.82391460349727) q[8];
u3(1.77516114150803,-5.28828963626373,0.0725446597798407) q[2];
u3(1.93844462651490,0.835039726066618,-1.34976104907339) q[5];
u3(1.46135427829164,-4.12357841931977,1.63847267418965) q[15];
cx q[15],q[5];
u1(0.596255355781603) q[5];
u3(-1.34218471655628,0.0,0.0) q[15];
cx q[5],q[15];
u3(-0.144603924243243,0.0,0.0) q[15];
cx q[15],q[5];
u3(0.961268529260614,-4.19848072042358,1.84584402801822) q[5];
u3(1.82976844570073,3.44088601831722,-2.64034703932355) q[15];
u3(1.17473931530990,-0.949748768406531,1.38520770238287) q[13];
u3(1.43112515426418,-1.37527403681909,-1.81383046779089) q[14];
cx q[14],q[13];
u1(1.14575717571153) q[13];
u3(-0.713591388826846,0.0,0.0) q[14];
cx q[13],q[14];
u3(-0.157993192533268,0.0,0.0) q[14];
cx q[14],q[13];
u3(0.869819753004736,3.14456227361774,-1.37068850591871) q[13];
u3(2.18845483133335,3.90124803717763,2.04638734910071) q[14];
u3(1.89372514358932,2.12360761947392,-3.59028526620400) q[7];
u3(1.07140375185622,3.30031958534869,-2.13654123262319) q[6];
cx q[6],q[7];
u1(2.52313192395171) q[7];
u3(-1.78354302554250,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.13144383572270,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.62152440375357,3.87215547158318,-2.37026681266539) q[7];
u3(1.74088036831868,-0.405747188080127,-5.38663187716689) q[6];
u3(1.90922988881457,1.33745964926643,-0.171441063929697) q[0];
u3(0.823278576858719,-0.506247774060977,-2.49220132900135) q[10];
cx q[10],q[0];
u1(-0.571671645791314) q[0];
u3(0.423428172695182,0.0,0.0) q[10];
cx q[0],q[10];
u3(4.38199898967631,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.04647166508111,-0.458090156962755,0.712954217728170) q[0];
u3(1.39921343405768,1.56424137467898,-3.30512592379750) q[10];
u3(1.49484544779482,-2.40481246272443,3.18127105956177) q[11];
u3(1.95908656738791,1.53617394165318,-1.70882064553372) q[13];
cx q[13],q[11];
u1(2.22736175128184) q[11];
u3(-1.76260486867000,0.0,0.0) q[13];
cx q[11],q[13];
u3(0.336128985959474,0.0,0.0) q[13];
cx q[13],q[11];
u3(2.14140229543809,0.00539282360392690,-2.59100474170404) q[11];
u3(1.08483274056891,-2.61294732696705,1.09594487326024) q[13];
u3(2.01150099221248,0.576292668372247,2.53981187667829) q[9];
u3(1.23027008659911,-2.65472463198292,-3.26558882764399) q[10];
cx q[10],q[9];
u1(-1.36384891896097) q[9];
u3(-0.335724704201385,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.63900307168637,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.03503927701557,-1.91681281117528,-0.568932003947164) q[9];
u3(1.19367021001838,-1.52480335102175,2.19783880128633) q[10];
u3(1.88041177134217,0.0193755583277465,0.194520050094557) q[1];
u3(0.261781474582110,-0.167095266694035,-4.93242094424521) q[0];
cx q[0],q[1];
u1(3.68395599390763) q[1];
u3(-3.76910156683602,0.0,0.0) q[0];
cx q[1],q[0];
u3(-1.03193311217298,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.33748595284749,-0.660256620616147,-0.953310110403150) q[1];
u3(1.34104060874760,-2.73216678888576,2.12927392646379) q[0];
u3(2.04595725165355,0.450057305026595,1.47347175998419) q[12];
u3(2.14031532425425,-1.47838382383487,-1.24349934409882) q[15];
cx q[15],q[12];
u1(1.16007202667260) q[12];
u3(-1.35925699872605,0.0,0.0) q[15];
cx q[12],q[15];
u3(2.61769298328442,0.0,0.0) q[15];
cx q[15],q[12];
u3(1.15633129443694,2.65665460012416,0.332474503332933) q[12];
u3(0.808282372175785,-0.706072386498322,-3.63501688728310) q[15];
u3(0.494487101082773,-0.921768785514586,0.222846582595495) q[4];
u3(0.248519229530159,-1.98010759671297,1.27194702415147) q[2];
cx q[2],q[4];
u1(1.92479437856955) q[4];
u3(-0.0455116192383083,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.25689559510371,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.95692016342653,-0.829631266063873,-2.45259958695663) q[4];
u3(1.30621855214239,0.285567707284736,5.79525362165715) q[2];
u3(1.15173392277185,2.13545969242406,-2.71629021815855) q[8];
u3(0.946747439931720,2.64524288106099,-3.47722902940671) q[14];
cx q[14],q[8];
u1(2.04724044622075) q[8];
u3(-2.72911891407070,0.0,0.0) q[14];
cx q[8],q[14];
u3(-0.0779827848649171,0.0,0.0) q[14];
cx q[14],q[8];
u3(0.216984112000009,-2.42195552713252,1.82908089639267) q[8];
u3(0.988712172545240,2.05851510771493,-2.54278950074553) q[14];
u3(0.460114034452986,2.75155822815106,-3.25852920425389) q[5];
u3(0.916579346890974,-3.20832840752674,2.60891269072613) q[6];
cx q[6],q[5];
u1(2.24004596827932) q[5];
u3(0.0895948921282581,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.25852457240162,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.75713235292833,0.549374703512765,2.30745068827971) q[5];
u3(1.27966186133109,-0.529725337547244,-4.37253266124829) q[6];
u3(2.52992489203562,2.07198267392062,-3.82529338462783) q[3];
u3(0.954066344600001,1.76005889517730,-0.470362624270200) q[7];
cx q[7],q[3];
u1(1.29782954224741) q[3];
u3(-3.32050771912900,0.0,0.0) q[7];
cx q[3],q[7];
u3(2.13940106795360,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.412161113431303,-0.0536540371011407,2.34342990338663) q[3];
u3(0.830553564435555,2.90858931181199,-0.473161128084824) q[7];
u3(0.997499127562010,0.199747755878330,-1.57673185211385) q[0];
u3(1.58553677628309,0.283216577974861,-5.08370715240237) q[7];
cx q[7],q[0];
u1(1.83781694436713) q[0];
u3(-2.22292991310326,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.93348506263901,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.160648488539442,3.35179544383579,-2.23771528485863) q[0];
u3(2.77229726072471,-1.49577262802484,-2.02064584576220) q[7];
u3(1.71213365596347,-2.25037929005781,1.40877284666696) q[13];
u3(1.77509074836961,-2.85751070646913,0.102793354175956) q[14];
cx q[14],q[13];
u1(1.38138181490923) q[13];
u3(-3.31609466812875,0.0,0.0) q[14];
cx q[13],q[14];
u3(2.49521984379016,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.12675486969789,-0.771715694194001,-1.66272702312648) q[13];
u3(1.03574464493104,1.70575935033567,0.755671485999243) q[14];
u3(1.97935934398001,-1.26567586851946,-1.44141920911236) q[9];
u3(0.440182959908296,-4.36066100177637,-0.343727491830998) q[3];
cx q[3],q[9];
u1(3.18969500054954) q[9];
u3(-0.889511021853107,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.64644812159220,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.455247268430678,-0.926571689528832,2.28687202932265) q[9];
u3(1.34212089727073,0.954384119244279,0.696874520302316) q[3];
u3(1.52263746857750,-1.02505142099077,0.299165982832789) q[6];
u3(1.40680531277263,-2.88431987704351,1.16603372801915) q[4];
cx q[4],q[6];
u1(-0.248090625142301) q[6];
u3(1.07715390104823,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.32110392700512,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.22391669660345,-2.39058473792873,2.18643113518622) q[6];
u3(0.577111010004401,-0.0599820448193933,-2.99504227405839) q[4];
u3(2.18401262825711,2.14120009340528,-2.25294143029056) q[5];
u3(1.23980498570327,-2.62564890072046,2.15588066734077) q[12];
cx q[12],q[5];
u1(3.12522174275415) q[5];
u3(-1.46636270151754,0.0,0.0) q[12];
cx q[5],q[12];
u3(0.585800667329370,0.0,0.0) q[12];
cx q[12],q[5];
u3(1.71427371815224,1.34262918986493,-2.09714030829783) q[5];
u3(2.39096955650844,-5.75042150050597,-0.502250631996818) q[12];
u3(1.84065141917114,-0.458028753365354,0.680189179701918) q[15];
u3(2.62434689280264,-1.44337921234321,-1.95111356750985) q[2];
cx q[2],q[15];
u1(0.778244364900943) q[15];
u3(0.0205827881240845,0.0,0.0) q[2];
cx q[15],q[2];
u3(1.92741147415512,0.0,0.0) q[2];
cx q[2],q[15];
u3(0.630466147638651,0.541592223388858,-2.08328847363502) q[15];
u3(0.734535456004742,0.867506103834546,3.51089276426539) q[2];
u3(2.85275682617959,3.14143387361713,-2.69266485631690) q[10];
u3(0.637873025653953,2.21456493962445,-0.908343145415191) q[8];
cx q[8],q[10];
u1(4.55008541784255) q[10];
u3(-3.27010223333784,0.0,0.0) q[8];
cx q[10],q[8];
u3(-0.166677126623588,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.66084034471566,-1.17013569513619,-0.844884095795207) q[10];
u3(1.61545593105227,2.12241346457826,-3.63596717688972) q[8];
u3(1.24473602513642,1.73903544830669,-2.38735717117138) q[1];
u3(2.00148681249339,-2.53119440817947,2.88979930482506) q[11];
cx q[11],q[1];
u1(2.60810923090763) q[1];
u3(-1.62158561415674,0.0,0.0) q[11];
cx q[1],q[11];
u3(3.47339062067023,0.0,0.0) q[11];
cx q[11],q[1];
u3(3.04236147136119,-2.62588386705152,-0.612239432120309) q[1];
u3(1.56384569150178,-2.36749227742042,0.620507888229269) q[11];
u3(1.20291599858992,1.86821055866574,-3.44303628879980) q[13];
u3(1.91800905228416,-2.47545567421774,2.72833292819243) q[8];
cx q[8],q[13];
u1(2.37989784217826) q[13];
u3(0.284168220835294,0.0,0.0) q[8];
cx q[13],q[8];
u3(1.10500401624646,0.0,0.0) q[8];
cx q[8],q[13];
u3(1.50091753856860,-2.75414032084763,0.436973706439145) q[13];
u3(2.58588481294661,3.67154571905477,1.28032478682330) q[8];
u3(2.38125599577778,2.00297528186735,-1.06181685377509) q[12];
u3(2.46393178551655,2.00630290057959,-3.31101219641156) q[15];
cx q[15],q[12];
u1(1.80179691709306) q[12];
u3(0.0298663132838626,0.0,0.0) q[15];
cx q[12],q[15];
u3(0.416976772815753,0.0,0.0) q[15];
cx q[15],q[12];
u3(2.54324899881761,0.248176646681783,3.23380846648636) q[12];
u3(2.32580267669117,-1.80984155684124,-2.74599885378464) q[15];
u3(1.98616337574951,0.630919184625698,2.32391362720433) q[2];
u3(1.27390402690120,2.51389464021538,3.25464546221297) q[9];
cx q[9],q[2];
u1(1.54721542490203) q[2];
u3(0.215461329376838,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.32953561580479,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.68561720913214,3.70131384763960,-2.22985610558740) q[2];
u3(0.868833351367806,-1.37942842312875,-0.946376491267171) q[9];
u3(2.54257118566609,1.76646730586058,-1.12504128097852) q[14];
u3(2.17977126291533,-0.479700440818815,-5.73726661014798) q[3];
cx q[3],q[14];
u1(-0.804572496175449) q[14];
u3(-1.77188498040897,0.0,0.0) q[3];
cx q[14],q[3];
u3(1.08663890007761,0.0,0.0) q[3];
cx q[3],q[14];
u3(0.845936232806836,2.09329344926745,-3.76168257799628) q[14];
u3(1.25684119512444,0.0653900632978266,-2.88564514249274) q[3];
u3(1.91590876056038,3.56145084830008,-1.49058450259695) q[6];
u3(2.17259440552073,0.508836910811753,-2.55119559841598) q[4];
cx q[4],q[6];
u1(1.18222477577650) q[6];
u3(-0.310602321158195,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.70808807063569,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.288925917945831,0.510604742592393,-2.82443543635832) q[6];
u3(2.14271563132013,4.31136856669147,0.154284726401217) q[4];
u3(1.06936605676721,-0.0908808237874836,1.46439550226178) q[7];
u3(0.663480606886850,-1.48098017412804,-2.36491577131372) q[10];
cx q[10],q[7];
u1(2.57800366189479) q[7];
u3(-2.94993457399707,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.795647377636789,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.443322194958632,3.87757621286259,-0.607136781036138) q[7];
u3(0.824171089159814,-3.08925431781373,0.0498579145004001) q[10];
u3(2.21579338297817,2.47424339316462,-2.26616020160753) q[11];
u3(0.636516960422127,-2.01277083877009,3.54480875918633) q[0];
cx q[0],q[11];
u1(1.19884549044918) q[11];
u3(-3.24611750160379,0.0,0.0) q[0];
cx q[11],q[0];
u3(2.24867224953583,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.447668659304265,-0.848782293788506,-1.06581577670245) q[11];
u3(1.86565299568992,3.90229115088046,-1.48683524705008) q[0];
u3(0.654110126416106,0.221771362418472,-1.56640880496423) q[5];
u3(1.44428124487785,0.303473051558829,-4.75664316872710) q[1];
cx q[1],q[5];
u1(1.41140983614931) q[5];
u3(-0.533266182910596,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.08487022608487,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.17701634714713,-4.32938976620701,0.401163845657448) q[5];
u3(1.53276007277538,-2.97699720762644,3.08836809390927) q[1];
u3(0.842645051605474,-2.15435028583630,1.64557981272565) q[10];
u3(0.377679634536129,-1.33070270995712,-0.915989582236700) q[12];
cx q[12],q[10];
u1(1.48824151688588) q[10];
u3(-0.730466448858406,0.0,0.0) q[12];
cx q[10],q[12];
u3(-0.461903183958146,0.0,0.0) q[12];
cx q[12],q[10];
u3(2.80062223190879,0.528169705134324,1.41763043918355) q[10];
u3(1.71494729894224,0.868018080901934,-4.01017617951249) q[12];
u3(1.35466259949497,0.810223273560136,-3.53467792263539) q[2];
u3(0.509245748523035,-2.47511218037592,2.15064288171254) q[3];
cx q[3],q[2];
u1(3.17975992123169) q[2];
u3(-1.33688492524228,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.72232258539499,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.91851409441831,-2.97345028981218,1.73039422132350) q[2];
u3(0.461340695116242,-4.15384066403033,0.0737534434397391) q[3];
u3(0.402632219055145,-1.99270059886227,1.97312226267044) q[13];
u3(1.01275160283850,-2.44872087780053,2.30739717913412) q[9];
cx q[9],q[13];
u1(0.588439551474466) q[13];
u3(-1.21205936392064,0.0,0.0) q[9];
cx q[13],q[9];
u3(0.182328664058131,0.0,0.0) q[9];
cx q[9],q[13];
u3(0.850501973492179,1.56404325525154,-1.65133634852776) q[13];
u3(1.85822658643294,0.780704952148269,-5.44214917417641) q[9];
u3(1.67434707541518,0.962337746135363,0.0835875566744182) q[11];
u3(1.52012170496683,-0.467309667025456,-4.55456452984374) q[8];
cx q[8],q[11];
u1(2.68387291165971) q[11];
u3(-2.61131722310624,0.0,0.0) q[8];
cx q[11],q[8];
u3(-1.38566434724549,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.09859099026928,4.46605973073048,-0.341101770723792) q[11];
u3(2.19448229753893,-1.64618481740807,1.14806894435836) q[8];
u3(1.01710412251099,1.36063028315444,-3.47677975132377) q[14];
u3(1.65382959161932,-2.23819223978431,3.68843934997768) q[0];
cx q[0],q[14];
u1(0.503849625212401) q[14];
u3(0.00284748481519848,0.0,0.0) q[0];
cx q[14],q[0];
u3(2.11736001133297,0.0,0.0) q[0];
cx q[0],q[14];
u3(0.760207634276368,1.30042584944079,-3.57437381492469) q[14];
u3(1.67255787497411,-0.474347666746884,-1.29232547279239) q[0];
u3(2.01300424757903,2.96674082330242,-0.455257972592748) q[7];
u3(1.49798624050072,1.72035043421760,-1.59272682296528) q[5];
cx q[5],q[7];
u1(1.72141467587842) q[7];
u3(-0.222557325606111,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.63310646388854,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.81302997703718,1.87326242354501,0.867695545963187) q[7];
u3(2.49930289573905,-1.52245279573876,-1.80915864084564) q[5];
u3(1.88316296979813,-1.36263482671003,-0.148092143727237) q[1];
u3(1.14821003664754,-2.02982486869046,0.391612517860898) q[4];
cx q[4],q[1];
u1(-0.335163232323403) q[1];
u3(-1.83894069786519,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.732457375457405,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.24347085719413,2.27476030424281,-1.99703976195922) q[1];
u3(2.29964304910511,1.01071561795094,-0.925615676899288) q[4];
u3(1.86861957754712,-0.487157842874340,0.696544401383927) q[15];
u3(1.04358332103914,-3.40511705637097,0.611392631129788) q[6];
cx q[6],q[15];
u1(1.82720717577445) q[15];
u3(0.180640697766911,0.0,0.0) q[6];
cx q[15],q[6];
u3(0.969180753460943,0.0,0.0) q[6];
cx q[6],q[15];
u3(0.370465387244043,2.70932291256023,-3.20789177290573) q[15];
u3(0.979052250976372,2.45258447751610,-3.32858545343426) q[6];
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
