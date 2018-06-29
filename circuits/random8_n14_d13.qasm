OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(2.43482497286478,-0.0490484399407058,2.95743281498904) q[10];
u3(2.59972196212245,-2.45704953242765,-1.62461132973324) q[13];
cx q[13],q[10];
u1(-0.0845971334315618) q[10];
u3(-1.94400431386484,0.0,0.0) q[13];
cx q[10],q[13];
u3(0.800885602222359,0.0,0.0) q[13];
cx q[13],q[10];
u3(0.737653190795995,1.02749484627873,0.720199978174121) q[10];
u3(1.69086411784923,0.970626618204531,1.54180932793614) q[13];
u3(1.44373989118672,3.24738703250517,-0.806624815235934) q[11];
u3(1.27034260954036,1.64941525676073,-1.12978159255755) q[12];
cx q[12],q[11];
u1(2.43325348277588) q[11];
u3(-1.76466505856516,0.0,0.0) q[12];
cx q[11],q[12];
u3(0.213756875961852,0.0,0.0) q[12];
cx q[12],q[11];
u3(2.48956008070409,-2.56173366181719,2.67918473195190) q[11];
u3(2.45399172910783,-1.19133928368087,-0.796252770154614) q[12];
u3(1.13263689816178,0.511753113611649,-2.01079508473431) q[2];
u3(1.79826779512449,-4.01887406109524,1.94899147207971) q[7];
cx q[7],q[2];
u1(3.35752478446171) q[2];
u3(-3.57280363824830,0.0,0.0) q[7];
cx q[2],q[7];
u3(-1.18310603213863,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.91467637723834,-0.149565675147614,-3.48042481874483) q[2];
u3(2.70295449242031,-0.466928149036476,-1.17941156750522) q[7];
u3(1.66082286989481,-0.166921549371738,1.89899942512116) q[8];
u3(1.31855364772163,-0.891110540091748,-0.820043407863828) q[1];
cx q[1],q[8];
u1(0.104083854540500) q[8];
u3(-0.718568820465679,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.68408954754258,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.74730777289126,0.873517596286367,0.902141638677697) q[8];
u3(1.09376927273837,-3.41312018026752,2.29714249338820) q[1];
u3(1.97815759325022,-0.366491484787997,2.23705641512654) q[0];
u3(1.41821942436703,-2.30641106724675,-2.33847848484626) q[3];
cx q[3],q[0];
u1(2.56869371257292) q[0];
u3(-1.60402341123480,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.46289605262683,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.910585507600179,-0.885156389766513,0.127185137335426) q[0];
u3(1.48148050299725,-2.74949403590101,-0.339253185145190) q[3];
u3(2.13179739824513,2.25868073988691,-2.57466507546673) q[6];
u3(2.09542945347980,-3.04724939135051,2.57077063889112) q[4];
cx q[4],q[6];
u1(1.90401695532644) q[6];
u3(0.0670523331855351,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.35110165982517,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.99107229113644,1.15445843556464,-1.84458013110307) q[6];
u3(2.30696701211364,-0.534668498104163,0.494144580990162) q[4];
u3(2.24620809045665,1.52499792505711,-3.25007538465022) q[5];
u3(1.59952672654040,2.25667152477437,-3.47696141455071) q[9];
cx q[9],q[5];
u1(2.41688071744436) q[5];
u3(-1.70886740820599,0.0,0.0) q[9];
cx q[5],q[9];
u3(3.39401125098886,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.08012937153120,-2.07677103119369,3.26800500101216) q[5];
u3(0.950129973366262,0.434715206202286,-2.58482913961990) q[9];
u3(2.73720874484125,2.40811266786658,-0.936056593028386) q[1];
u3(2.52014850902964,4.71911785654340,0.919545801419078) q[3];
cx q[3],q[1];
u1(2.28400866038017) q[1];
u3(0.400644629675597,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.63134452500760,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.56182265440758,2.12172335219558,-0.925884277213780) q[1];
u3(0.408725247609865,-1.84019605361806,3.91957614706576) q[3];
u3(1.90632044485607,1.43501834267777,0.292439944872258) q[12];
u3(1.38246989241243,-1.23615561356602,-2.57005507388847) q[2];
cx q[2],q[12];
u1(3.28902984750217) q[12];
u3(-4.14747447731479,0.0,0.0) q[2];
cx q[12],q[2];
u3(-0.700501344651273,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.21885554360749,0.482725874059318,0.760047449037454) q[12];
u3(1.01512958401963,4.52796065391435,-0.344694771954953) q[2];
u3(1.45451240262913,-0.898198522193419,2.66428668181812) q[11];
u3(1.18657213647699,-0.672023427309624,-1.25503267360574) q[9];
cx q[9],q[11];
u1(1.38597591871376) q[11];
u3(-0.532757965108375,0.0,0.0) q[9];
cx q[11],q[9];
u3(-0.247725016768329,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.77809331576903,3.51933886517939,-0.739727253051001) q[11];
u3(2.99153688641807,-1.33282783103478,-3.15696657390598) q[9];
u3(2.68856664297811,-3.94716552633736,2.25575470814493) q[7];
u3(0.968169622262945,4.50177677770172,-1.74532240488178) q[10];
cx q[10],q[7];
u1(0.254018504603125) q[7];
u3(-1.09142082015364,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.55643484838361,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.901891377594165,1.94226028378791,0.0215844054809121) q[7];
u3(1.37674924705323,0.147821456587640,1.05132689611536) q[10];
u3(2.47851123323613,-2.33997514256962,-0.482085106938092) q[0];
u3(2.51876813324115,-2.89297947963640,-1.37085587310809) q[5];
cx q[5],q[0];
u1(0.732873602126344) q[0];
u3(-1.24885251436306,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.0614444469358826,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.80136704984787,1.09849321784029,-3.33541302537815) q[0];
u3(0.629254504093054,2.07690669059532,2.69586775369638) q[5];
u3(2.62263461207368,1.96706268974810,-1.25712189418379) q[8];
u3(2.52240643039581,0.0433339679364764,-4.41847919620310) q[13];
cx q[13],q[8];
u1(0.427154592824805) q[8];
u3(-1.15071513704698,0.0,0.0) q[13];
cx q[8],q[13];
u3(2.37162015299567,0.0,0.0) q[13];
cx q[13],q[8];
u3(2.40407137857465,-0.116142221264161,-1.33833562297324) q[8];
u3(1.47851240778331,0.775513457275175,-4.60291793664506) q[13];
u3(1.69683680927698,4.54350161926879,-1.67885155324189) q[6];
u3(0.153489698471633,0.530158085117980,0.298170056411072) q[4];
cx q[4],q[6];
u1(0.265745230394362) q[6];
u3(-1.05566366476278,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.61271892486541,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.05582292709414,2.29533265764099,-2.50662537570745) q[6];
u3(1.77275934054066,-1.87610838044570,1.05605009226187) q[4];
u3(2.30517537627050,0.133479063426154,-2.27448618598307) q[1];
u3(1.29819216718567,-3.79362581528921,1.79608830058130) q[6];
cx q[6],q[1];
u1(2.29572209646285) q[1];
u3(-2.00255922631636,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.06295019266329,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.90354522821021,-0.410465672882479,0.716828490109536) q[1];
u3(1.72814690411513,-0.0346391652464884,4.84816794395490) q[6];
u3(1.44678158605001,-0.965671383115220,0.759183608302803) q[7];
u3(1.67347133496638,-1.41840499574586,-1.18481935750063) q[3];
cx q[3],q[7];
u1(1.58129629470616) q[7];
u3(-2.46466183175816,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.53564199869425,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.66559041049103,0.970965172906128,1.49959947462700) q[7];
u3(1.08051491357381,-3.64092240779286,-2.39746910191289) q[3];
u3(0.543822389694132,1.00327084020305,-1.70368660054356) q[10];
u3(0.857897192742744,-0.934879455116343,-0.335219237451776) q[13];
cx q[13],q[10];
u1(1.37820760934353) q[10];
u3(0.458415966681438,0.0,0.0) q[13];
cx q[10],q[13];
u3(0.653499563743111,0.0,0.0) q[13];
cx q[13],q[10];
u3(1.82930050563445,-2.08743204605657,0.977654979360574) q[10];
u3(1.71637599374950,0.901281808178383,0.495307105437045) q[13];
u3(1.87323308368746,0.621996731184173,1.45633038246666) q[2];
u3(2.15470724451466,-1.32856285309979,-0.790085547547106) q[5];
cx q[5],q[2];
u1(0.930702934903032) q[2];
u3(-1.54511020664378,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.50989517874355,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.610564600889438,-0.741742961914243,-1.28392590380874) q[2];
u3(0.782348339921213,-2.79551371793004,0.216761025995736) q[5];
u3(0.647460044277765,-0.0716607292015884,0.594881168909978) q[0];
u3(0.887128520723413,-2.41501818615784,0.204360651875903) q[11];
cx q[11],q[0];
u1(1.73538271461005) q[0];
u3(-3.21935262348144,0.0,0.0) q[11];
cx q[0],q[11];
u3(2.19937928894705,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.11254066404372,-1.06003487749415,0.220980497203961) q[0];
u3(1.46095169854000,-4.54956291389376,0.0668103699949660) q[11];
u3(2.20997038662322,-0.0505897247526219,1.23099344753030) q[12];
u3(1.43980826807821,-2.79287316500087,-2.05532221660115) q[4];
cx q[4],q[12];
u1(2.51636031068510) q[12];
u3(-1.63855327645859,0.0,0.0) q[4];
cx q[12],q[4];
u3(3.22500455481690,0.0,0.0) q[4];
cx q[4],q[12];
u3(0.353932097507471,-1.09893542220701,-0.270927509973825) q[12];
u3(2.40966084866004,-0.157698528211570,0.957762807487792) q[4];
u3(2.21122156065226,2.19063396751015,-2.55253745140030) q[9];
u3(1.35496286102620,1.61518154377511,-2.58670038744772) q[8];
cx q[8],q[9];
u1(2.24251251245607) q[9];
u3(-1.85338685580953,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.74828458632028,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.48680242432067,2.94215708743959,-1.34502815068522) q[9];
u3(1.88938440900951,1.08670104479126,1.14893513934478) q[8];
u3(0.710951605106401,-0.655265530835805,0.559122172839364) q[13];
u3(1.49840490184528,-2.69088051056525,0.783718298100929) q[2];
cx q[2],q[13];
u1(1.73313138465653) q[13];
u3(-2.06897549259815,0.0,0.0) q[2];
cx q[13],q[2];
u3(3.91753859587443,0.0,0.0) q[2];
cx q[2],q[13];
u3(1.96455172789748,-1.89133248460579,1.00992807747715) q[13];
u3(2.84765758427067,-0.580535357290541,-0.0247573691737208) q[2];
u3(2.15854044641500,2.11604984306135,-2.81986142840980) q[3];
u3(0.980791444009532,2.36235326077727,-2.76017650979365) q[1];
cx q[1],q[3];
u1(2.87562719234252) q[3];
u3(-1.95508789948901,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.984964851995365,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.48346262776020,4.69126034993830,-1.05484610601293) q[3];
u3(0.451763946213156,-0.0826880487655068,-4.04688676030822) q[1];
u3(0.456723117598337,2.10269438370336,-4.07135546076935) q[0];
u3(1.76423835454640,2.94765546653218,-2.40938237864916) q[7];
cx q[7],q[0];
u1(1.57855904512060) q[0];
u3(-2.10352331957490,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.544888035163325,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.59377262573364,1.71638144249119,0.0559773461590803) q[0];
u3(1.41482183028806,0.967560727536643,-3.96663878158111) q[7];
u3(0.602022518150156,-0.241418407116808,-0.118354813627365) q[9];
u3(0.822955488534219,-1.10997894444498,-1.51671328218304) q[11];
cx q[11],q[9];
u1(1.29744552756032) q[9];
u3(-0.806745619379693,0.0,0.0) q[11];
cx q[9],q[11];
u3(-0.371747216015968,0.0,0.0) q[11];
cx q[11],q[9];
u3(2.54491931171354,-2.96011739369230,2.64258614309475) q[9];
u3(0.812783209567275,-0.695602346349069,-1.46839270813689) q[11];
u3(1.28283541002786,0.566745683860369,0.814445611111681) q[8];
u3(1.98637591482048,-0.603981619848359,-2.05434054898031) q[4];
cx q[4],q[8];
u1(1.14771987323605) q[8];
u3(-0.520713513998160,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.73310570007761,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.13997397989648,1.86956309795909,-4.02506340947748) q[8];
u3(2.10114822386479,1.64028921767174,-2.88442884128638) q[4];
u3(1.40849980962292,1.83387186309967,-2.50786159871065) q[10];
u3(1.20250597571371,-2.86795177823549,2.55999503799824) q[5];
cx q[5],q[10];
u1(3.34318801297604) q[10];
u3(-3.90101292915373,0.0,0.0) q[5];
cx q[10],q[5];
u3(-0.922593608240006,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.34242087245455,-2.45476738531369,3.54217965736321) q[10];
u3(1.70451738205553,-2.59741411052463,-0.750065333134505) q[5];
u3(1.13792034835207,0.987200341614560,-1.07009927834247) q[12];
u3(0.597183437788492,1.72472190001796,-4.54802135398980) q[6];
cx q[6],q[12];
u1(1.76311534344195) q[12];
u3(-2.38812288967204,0.0,0.0) q[6];
cx q[12],q[6];
u3(0.0221518148339146,0.0,0.0) q[6];
cx q[6],q[12];
u3(0.656845433008960,-2.11000302296658,1.98889555976583) q[12];
u3(1.48003917240526,2.62873131954136,0.525283054794442) q[6];
u3(1.90536161832379,3.91975039718701,-1.11336554871497) q[9];
u3(1.48192948677440,2.45071819968581,-0.504291418489629) q[12];
cx q[12],q[9];
u1(0.0955847492775641) q[9];
u3(-0.956091221167351,0.0,0.0) q[12];
cx q[9],q[12];
u3(1.72920492526899,0.0,0.0) q[12];
cx q[12],q[9];
u3(1.06287502840654,3.73394739217664,-1.88921706615714) q[9];
u3(0.719980577113018,-0.747520490663903,2.45180385256434) q[12];
u3(1.53330998303716,1.77305829489985,0.109041478054460) q[4];
u3(2.65799212255782,0.125634003658366,-2.07842809616393) q[11];
cx q[11],q[4];
u1(0.554811008105641) q[4];
u3(-1.34025060923357,0.0,0.0) q[11];
cx q[4],q[11];
u3(-0.286832278024369,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.05170548971151,-0.348791218473034,-1.03026963057983) q[4];
u3(2.77656908724942,-1.28930763541592,3.17340319392278) q[11];
u3(2.35321357865309,1.49632511199933,0.0202982974061462) q[0];
u3(1.51283087752677,-1.12454045827394,-2.76120415552378) q[13];
cx q[13],q[0];
u1(-1.18670344511153) q[0];
u3(0.711948652625773,0.0,0.0) q[13];
cx q[0],q[13];
u3(3.87236841794897,0.0,0.0) q[13];
cx q[13],q[0];
u3(1.97381311302504,-3.92924330983561,1.79463190961596) q[0];
u3(2.14377190907858,-1.71263345865649,-1.28391665380371) q[13];
u3(1.79875536777040,1.48106077889309,-1.92811248377311) q[5];
u3(2.82821947682666,2.08224403855537,-4.16569364662352) q[1];
cx q[1],q[5];
u1(3.24901536085950) q[5];
u3(-1.57149866209286,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.43963620913110,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.57316496289298,2.38063962962092,-1.78201437326096) q[5];
u3(2.13934412847176,-0.565346351237644,-3.23382954736886) q[1];
u3(1.98704138698402,1.83756896066263,-3.07349797920103) q[10];
u3(1.61521107864797,-2.55306607733720,2.61870158323734) q[3];
cx q[3],q[10];
u1(1.53613320807031) q[10];
u3(-0.825545872571290,0.0,0.0) q[3];
cx q[10],q[3];
u3(-0.291961738563047,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.03134765883385,4.49941738329159,-0.947678818012475) q[10];
u3(0.610135611314784,3.48047407938585,2.07497450006371) q[3];
u3(0.0336636692303671,0.512883826970884,-1.03523791333248) q[6];
u3(1.01527132899778,-0.936281654338219,-1.32155872997747) q[7];
cx q[7],q[6];
u1(1.47125595545395) q[6];
u3(-3.36557268560461,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.66033328315106,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.32686678446907,-0.780336246679264,1.34222184044124) q[6];
u3(0.618338561463814,-0.524696997369487,-1.91811256197797) q[7];
u3(1.95602262041886,0.807199080853525,1.75951709732832) q[8];
u3(1.69042762708234,-1.17330125019949,-0.416889643273651) q[2];
cx q[2],q[8];
u1(0.150921910252241) q[8];
u3(-1.27544615900618,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.08131541852287,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.31327642500608,-0.116621605817367,2.44567503601547) q[8];
u3(2.03877483276960,-2.77292392363840,-0.206256526093738) q[2];
u3(2.80272168896506,-1.97984204825757,1.83086937901036) q[2];
u3(2.74077169919035,-2.71858699427984,-1.28565491465439) q[1];
cx q[1],q[2];
u1(-0.515273753807574) q[2];
u3(1.19820784597333,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.43117085082869,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.35411707710259,-2.00297694676126,1.86533738533327) q[2];
u3(1.70783387192140,2.22913986023525,2.94492450901067) q[1];
u3(0.443752781858380,1.90394669083463,-3.66023557980272) q[11];
u3(1.66183267044828,-2.56901526411689,3.61583537154580) q[5];
cx q[5],q[11];
u1(-0.0538687442109438) q[11];
u3(-0.904976851209210,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.58824909186460,0.0,0.0) q[5];
cx q[5],q[11];
u3(0.806312484446569,-0.460322421505087,-1.87627722530051) q[11];
u3(1.27748340828277,4.34244262283858,-0.749755542391163) q[5];
u3(1.19765900201881,-1.49321796217757,1.14167634344676) q[4];
u3(0.907432085598008,-2.09510461279919,-0.638529839139623) q[12];
cx q[12],q[4];
u1(1.13272704199027) q[4];
u3(-0.598588187053889,0.0,0.0) q[12];
cx q[4],q[12];
u3(2.85756750475211,0.0,0.0) q[12];
cx q[12],q[4];
u3(2.58761737801280,-1.05326461539952,-0.279579122635527) q[4];
u3(0.159799761188668,-3.30911631117235,2.44812788547406) q[12];
u3(2.95688330147253,-1.10329303190176,1.92112313013761) q[7];
u3(2.20299207522356,1.27233014417938,4.07606328112046) q[9];
cx q[9],q[7];
u1(1.09627965142202) q[7];
u3(-3.21834311338585,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.55004099576230,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.286988627002288,-4.30420654043468,1.61853052143244) q[7];
u3(1.05338191181465,-1.69478228005429,-4.45440913615537) q[9];
u3(2.31105529618901,0.413762670143759,-2.96182349237653) q[0];
u3(1.97670627055957,-0.456827626614559,-5.01107926870733) q[13];
cx q[13],q[0];
u1(3.23378660629427) q[0];
u3(-4.23216647137709,0.0,0.0) q[13];
cx q[0],q[13];
u3(-0.352833831629224,0.0,0.0) q[13];
cx q[13],q[0];
u3(2.00507109202714,-1.52801368774826,1.88533197299978) q[0];
u3(2.53639705947642,2.37408152761518,3.12434553908560) q[13];
u3(1.26754910938489,0.871387767962975,1.93921812641882) q[6];
u3(0.863128282275519,-1.21571563649092,-2.04419649088267) q[3];
cx q[3],q[6];
u1(2.25432008082417) q[6];
u3(-1.67853109603265,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.566262947335324,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.987757450931778,0.492395841210947,0.840972170235740) q[6];
u3(1.92740111570030,-2.20275094846479,-3.74225203700480) q[3];
u3(1.08648288480890,-3.46454467094500,2.44888051408822) q[8];
u3(2.01494435033502,-2.66744619269127,3.22136817706735) q[10];
cx q[10],q[8];
u1(-0.0206031873210879) q[8];
u3(0.867604556883619,0.0,0.0) q[10];
cx q[8],q[10];
u3(3.66595923971028,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.67425512856033,2.10741651115423,-1.29712704085970) q[8];
u3(1.71817794633757,-5.17384097200846,-0.267426035435991) q[10];
u3(2.62691392796581,2.61642524000769,-1.25295569404152) q[12];
u3(1.41813555943864,1.62488478060316,-2.24059711296949) q[1];
cx q[1],q[12];
u1(-0.332428920598350) q[12];
u3(0.645474549102830,0.0,0.0) q[1];
cx q[12],q[1];
u3(4.11266858029627,0.0,0.0) q[1];
cx q[1],q[12];
u3(2.11478188720829,3.23612947810929,-2.42619021833065) q[12];
u3(0.728566192887262,-3.62920770408047,-1.56309386198494) q[1];
u3(1.77786448570762,-0.176127037701142,-1.44806219498206) q[11];
u3(1.95273911395307,-3.55773704344871,1.60990587816431) q[7];
cx q[7],q[11];
u1(1.97898763873666) q[11];
u3(-3.27797565236436,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.442226298143729,0.0,0.0) q[7];
cx q[7],q[11];
u3(0.293153861112313,-2.71324052868095,3.27495805228338) q[11];
u3(0.559015893366884,-3.11874675844461,-1.63788341311316) q[7];
u3(1.26058732932150,-1.11604678333061,1.99714883658695) q[4];
u3(0.878108781672354,-2.03322752785053,-1.72519845783689) q[6];
cx q[6],q[4];
u1(-0.0579267839347168) q[4];
u3(-2.32386716007282,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.68225568720373,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.70097191594647,-2.01228022309208,2.48449433875542) q[4];
u3(2.12053167212297,4.32029437648168,-0.202305561588561) q[6];
u3(0.535645126310377,-3.02856736955781,3.19865074953102) q[3];
u3(1.11441620471492,0.324480106157958,-2.76579885924314) q[9];
cx q[9],q[3];
u1(-0.172344605052928) q[3];
u3(-1.10569599214866,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.39836623702473,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.69214954486830,-3.12144258977780,-0.0518652545248797) q[3];
u3(0.472577476338651,0.0693084649940237,-3.95833878907204) q[9];
u3(1.52310908545047,-0.277727310243980,0.857293304649031) q[2];
u3(1.53471678399465,-2.60628513366934,-1.03568362071685) q[10];
cx q[10],q[2];
u1(3.40221169715121) q[2];
u3(-1.57347988860806,0.0,0.0) q[10];
cx q[2],q[10];
u3(2.39145414533105,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.30970148331624,-1.65194188159561,3.99975588139175) q[2];
u3(1.00914001879373,1.96060573261689,3.89390099627551) q[10];
u3(2.41257983148106,2.32618315055410,0.523074400533600) q[8];
u3(1.67928713010915,-0.763169246850090,-2.63916506980981) q[5];
cx q[5],q[8];
u1(2.94478295207906) q[8];
u3(-1.77418923723174,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.447082863907014,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.553446322341934,-1.37869267944128,3.91352623661670) q[8];
u3(0.286680203213180,0.782585727638293,4.31968543049457) q[5];
u3(0.526567186078555,1.94114580135364,-0.322353817757341) q[13];
u3(1.73160243176627,-0.919436639036279,-3.54716125714586) q[0];
cx q[0],q[13];
u1(-0.113160129959270) q[13];
u3(0.802965484346911,0.0,0.0) q[0];
cx q[13],q[0];
u3(3.55732916707702,0.0,0.0) q[0];
cx q[0],q[13];
u3(0.808406176269324,1.82081993222957,0.827661763996419) q[13];
u3(1.21441089487721,1.61259071875402,3.98037274217985) q[0];
u3(2.09742914459994,-0.994044140260204,3.45870450957919) q[1];
u3(2.94544337123839,-2.52596107106254,-1.10377711003939) q[8];
cx q[8],q[1];
u1(2.09716939840863) q[1];
u3(-2.53152382840205,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.0687979637864389,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.439971006207255,-3.90791291925571,1.69023970222328) q[1];
u3(2.58379336879873,-3.80239511720955,-1.81447827535308) q[8];
u3(2.03396266554683,-0.357253871991320,-0.644277125554078) q[7];
u3(1.86574516465168,-2.79833271127092,0.537787189096832) q[4];
cx q[4],q[7];
u1(2.59254695556288) q[7];
u3(-1.61856214506493,0.0,0.0) q[4];
cx q[7],q[4];
u3(3.25368612865379,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.48277227440965,-3.13894917322947,0.118039485482829) q[7];
u3(1.96419055310246,-4.21116709581505,0.916613984049640) q[4];
u3(2.61820629386293,1.84373333754827,-1.39180868842883) q[11];
u3(2.39557084227730,4.29549904398163,0.698457371616419) q[6];
cx q[6],q[11];
u1(1.67848249113626) q[11];
u3(0.437391023401233,0.0,0.0) q[6];
cx q[11],q[6];
u3(0.781827045203633,0.0,0.0) q[6];
cx q[6],q[11];
u3(0.946921223792786,0.995798094092815,0.615551666203808) q[11];
u3(1.36171011711770,-4.62340655207886,-0.993876198464759) q[6];
u3(2.51885839871434,-1.25246464518300,-1.30758668845094) q[2];
u3(0.405992446537812,-1.42545675133782,-3.40274335921837) q[13];
cx q[13],q[2];
u1(4.35728668515740) q[2];
u3(-3.60338860912554,0.0,0.0) q[13];
cx q[2],q[13];
u3(-0.134980656048776,0.0,0.0) q[13];
cx q[13],q[2];
u3(0.592405884437257,-2.48888211661732,0.213747932193742) q[2];
u3(0.640817119788670,-0.324826940884762,5.37010988632773) q[13];
u3(1.69027288978234,-0.536502257098306,0.255693877114780) q[12];
u3(2.13336320544491,-0.663532778517662,-1.23254241309058) q[3];
cx q[3],q[12];
u1(-0.0374571867109947) q[12];
u3(-1.85203794879181,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.528062095794211,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.72880200762782,-0.533286753510776,1.78827797989173) q[12];
u3(2.10941961799793,-1.56513130759284,1.18641299787776) q[3];
u3(2.09294831452514,0.970313459737694,0.165802914914714) q[10];
u3(0.655598925559319,1.01561748169539,-3.95699632465986) q[0];
cx q[0],q[10];
u1(1.62773179752860) q[10];
u3(0.324958336386060,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.456127704635833,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.69187448643015,-3.61446068687866,1.86204381712322) q[10];
u3(1.55277705928120,-0.440351092125131,2.58046714764860) q[0];
u3(0.153616761610364,2.32888890648648,-1.66583996227521) q[5];
u3(0.277177801921173,-2.84699585411772,1.32122642080863) q[9];
cx q[9],q[5];
u1(3.21454705597332) q[5];
u3(-1.98331710491170,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.60391150155239,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.48377954452126,2.75858934493113,-1.77799729949218) q[5];
u3(1.29277942251397,1.39692500016199,1.02353730505369) q[9];
u3(1.64111533923307,2.95250246382476,-0.864546822380269) q[11];
u3(2.32417949062682,1.17313982481290,-1.85247633892057) q[3];
cx q[3],q[11];
u1(1.44959599882589) q[11];
u3(-0.845927796421164,0.0,0.0) q[3];
cx q[11],q[3];
u3(-0.357340501622769,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.72338924628184,2.47099735132241,0.343340933290279) q[11];
u3(0.347813669200919,-2.55875593158800,-2.45649437666403) q[3];
u3(1.56460555924397,0.185476180376384,2.70877949603361) q[12];
u3(1.83438929164680,-2.11571260578315,-1.93223669771896) q[8];
cx q[8],q[12];
u1(2.32561138276632) q[12];
u3(-3.12645759540336,0.0,0.0) q[8];
cx q[12],q[8];
u3(1.79713153106819,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.11617798861954,-0.691750464578161,0.398363898408527) q[12];
u3(0.861309721739066,-0.0871396756189287,-3.27429402661807) q[8];
u3(2.12090462177349,3.13207493793236,-2.40422050457348) q[10];
u3(1.41938839703934,2.92841781995945,-2.53697671867299) q[2];
cx q[2],q[10];
u1(1.58704730051956) q[10];
u3(-0.545917849919403,0.0,0.0) q[2];
cx q[10],q[2];
u3(2.09240782190042,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.54606834461214,-1.56636443535160,4.37171369202596) q[10];
u3(0.810422318621638,-3.29288912136545,-1.87020462278554) q[2];
u3(0.161431907493123,1.99825518707543,-1.42843178506661) q[6];
u3(0.708220264381079,0.0766797867065238,-0.822338273131133) q[4];
cx q[4],q[6];
u1(1.45807274124112) q[6];
u3(-0.718772375153891,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.91275921037134,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.44195599830740,-0.766752708134925,0.568129774658644) q[6];
u3(1.76086337826091,0.468063002608075,4.27077746223680) q[4];
u3(1.85452724527522,2.70682954825091,-2.95283983537764) q[5];
u3(1.86240677804389,0.245778860527453,-1.71406866586954) q[7];
cx q[7],q[5];
u1(2.77099390254438) q[5];
u3(-2.24141223805402,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.807609110757033,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.59461327466003,2.49463728900303,-3.41944382739083) q[5];
u3(1.05014415595574,-0.388242223470494,-4.91722648036562) q[7];
u3(1.22978078749772,-3.53441510962994,2.03787048577903) q[1];
u3(1.41699829608557,2.92296388817108,-2.50323449807434) q[9];
cx q[9],q[1];
u1(1.51279147583375) q[1];
u3(-0.846361084715633,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.33903901835161,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.64679268162139,-3.49834985698633,0.0260743040221398) q[1];
u3(1.28203992407056,1.04107643629917,-2.65993003838625) q[9];
u3(2.30236364365582,-0.425253024660496,0.968367027699534) q[13];
u3(2.03076530426483,-2.42248200566561,-1.59750734431841) q[0];
cx q[0],q[13];
u1(0.841206159533828) q[13];
u3(-0.246475226271955,0.0,0.0) q[0];
cx q[13],q[0];
u3(1.74028585647038,0.0,0.0) q[0];
cx q[0],q[13];
u3(1.75722203273644,1.09499697840673,-0.795635241798395) q[13];
u3(2.96876831280620,1.28547149040088,0.883275322253915) q[0];
u3(2.53338900140112,-0.0896620295796339,0.147545854015334) q[11];
u3(1.28024444316189,-2.76838630711766,-1.91161662650893) q[10];
cx q[10],q[11];
u1(3.69513677088902) q[11];
u3(-1.42466768819828,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.34549959528043,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.32427766493803,-1.42039818368205,4.53860750144493) q[11];
u3(2.72714548667817,2.64853746812878,2.18000188837837) q[10];
u3(1.97309256435901,-1.49443765806586,-1.27445543725206) q[0];
u3(0.569809487594116,-5.30926815361424,0.667147082597948) q[2];
cx q[2],q[0];
u1(0.127708509171057) q[0];
u3(-2.41502467052858,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.24169646459841,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.31369098879698,3.11125338093704,0.0196037287309585) q[0];
u3(0.949506565259650,-2.56424537000274,3.29312878384263) q[2];
u3(1.77411400152512,1.43090097807692,-2.70681217241477) q[12];
u3(1.78043611599539,2.29127596211387,-3.80340401370480) q[8];
cx q[8],q[12];
u1(1.93388922590320) q[12];
u3(-2.86622429747151,0.0,0.0) q[8];
cx q[12],q[8];
u3(1.01894520680408,0.0,0.0) q[8];
cx q[8],q[12];
u3(2.62347181326817,3.67697814257154,-1.75952789238200) q[12];
u3(2.21409865125132,0.0425713900071005,5.56055955570718) q[8];
u3(0.572406733165233,-2.70689025785940,3.26069066320199) q[1];
u3(1.16892166714422,1.56789623390401,-1.80862307043955) q[6];
cx q[6],q[1];
u1(1.72339288413059) q[1];
u3(0.0301370836813608,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.487410919655003,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.513574769276250,-1.88525579667147,2.51014813430593) q[1];
u3(2.68838666255539,0.237062534639712,4.00684324359004) q[6];
u3(1.61641382213688,3.13203084328364,-1.96893181341209) q[9];
u3(2.50797351648298,1.03916128426973,-0.858294714082180) q[4];
cx q[4],q[9];
u1(0.215375807363757) q[9];
u3(-0.558202559116196,0.0,0.0) q[4];
cx q[9],q[4];
u3(2.09240463105703,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.688249131595061,-1.47777468995921,0.101994020633155) q[9];
u3(1.43404946366560,-0.410941789945874,5.66599475704237) q[4];
u3(1.51366512863167,-0.183751357805593,1.96620332420582) q[7];
u3(1.53931010808441,-2.16270422201589,-1.02744259770832) q[5];
cx q[5],q[7];
u1(4.28906066110211) q[7];
u3(-3.25800350857592,0.0,0.0) q[5];
cx q[7],q[5];
u3(-0.516141670640930,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.37093030597480,1.08346449414618,-0.0722187746904746) q[7];
u3(1.44733227177694,-0.751180158446391,0.621598722498438) q[5];
u3(0.865664878395870,1.94487072863942,-2.39548966162697) q[3];
u3(0.895094659510715,-4.44137994973851,1.58641772563254) q[13];
cx q[13],q[3];
u1(1.70317223077008) q[3];
u3(0.263543063557502,0.0,0.0) q[13];
cx q[3],q[13];
u3(0.730630480961198,0.0,0.0) q[13];
cx q[13],q[3];
u3(1.47346829212060,2.17733806733250,-2.66138923637602) q[3];
u3(1.77265861975823,0.953136895214053,4.47895031628681) q[13];
u3(0.803026380378771,-0.0476551438655735,-2.12340482208120) q[8];
u3(1.68107941164267,0.586152334565512,-4.96807158998691) q[13];
cx q[13],q[8];
u1(-0.412260693534717) q[8];
u3(-1.97389996615593,0.0,0.0) q[13];
cx q[8],q[13];
u3(1.68294296070341,0.0,0.0) q[13];
cx q[13],q[8];
u3(0.721905067467346,0.422303252086973,1.07057307146256) q[8];
u3(1.82968243947902,-3.01657612864551,-0.150929644706588) q[13];
u3(1.72846838409037,0.595783206733409,1.41454455137131) q[10];
u3(0.989438307541196,-1.38680313033733,-2.72495260324895) q[5];
cx q[5],q[10];
u1(-0.902290416462957) q[10];
u3(-1.90123021765333,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.29051098938995,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.74465370347944,1.33481429017634,-2.67414737909629) q[10];
u3(1.75167600290589,2.66328176750317,-1.74169208311857) q[5];
u3(0.382064736761819,1.08821902943695,-1.32802024485636) q[6];
u3(0.865496724905535,0.105372322211445,-0.521099014481826) q[7];
cx q[7],q[6];
u1(1.75341651253776) q[6];
u3(-2.20123515919585,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.601495620221740,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.48081743704167,-0.100695862041306,-0.377619719326317) q[6];
u3(1.48451156325767,-3.16379869598744,-2.45495108508957) q[7];
u3(0.590410667059871,-0.0689920364582020,0.892584888136311) q[2];
u3(1.15340410185945,-0.860391443308351,-1.25940025431106) q[0];
cx q[0],q[2];
u1(2.28984994030556) q[2];
u3(-0.00697290891619429,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.28654624642220,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.10092934805105,1.95383165587366,-4.09822254939239) q[2];
u3(0.774589717124863,0.953788282895078,0.610185343583188) q[0];
u3(1.93766085897740,-1.34025814846717,2.55417118001369) q[11];
u3(2.12089380264500,-1.29423306519009,0.712546979152680) q[1];
cx q[1],q[11];
u1(-1.17612421569729) q[11];
u3(0.247032525996362,0.0,0.0) q[1];
cx q[11],q[1];
u3(3.52097702279728,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.79498349535341,2.56611246782629,-1.37110776283678) q[11];
u3(1.77520016445903,-0.995591312568758,3.20895487246154) q[1];
u3(1.74036124654571,-0.149967279683438,1.86238749159680) q[12];
u3(1.41730547891887,-2.58055335570531,-2.42975307533606) q[3];
cx q[3],q[12];
u1(0.803401519972753) q[12];
u3(-0.506111610343070,0.0,0.0) q[3];
cx q[12],q[3];
u3(1.22100551850711,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.25411895205150,2.64910847439447,-2.32662825631689) q[12];
u3(1.34329151943138,-1.25881880449347,-1.36925501795272) q[3];
u3(1.18402427668755,1.49805005051793,0.990583245935621) q[4];
u3(0.979682486385522,-1.12832899612029,-1.98909382236462) q[9];
cx q[9],q[4];
u1(1.00104562158475) q[4];
u3(-0.232457384075456,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.65696308979052,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.49653786744086,1.70892879389187,-2.71055801659525) q[4];
u3(1.80420704242826,1.74628656110220,3.92167527030943) q[9];
u3(1.05510304267707,0.864096473403173,-1.58170488676697) q[9];
u3(0.529129070205153,-3.23986529900527,0.856662369323833) q[0];
cx q[0],q[9];
u1(1.41734011415204) q[9];
u3(-0.145351225865495,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.40415045673793,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.30821702251744,0.879913394531702,-0.0617908280501142) q[9];
u3(0.885588819691899,-3.65011628520661,0.401283339891092) q[0];
u3(2.05746441179456,-0.0678871826854717,-1.21935464941793) q[7];
u3(1.43181590848272,-3.56637755829605,0.739844222758013) q[1];
cx q[1],q[7];
u1(0.765181573151784) q[7];
u3(-3.33246103191907,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.69532713234415,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.24932371562581,-0.693032247176663,0.0204806022510546) q[7];
u3(1.82179905057438,-1.96135470921912,-2.26728293089172) q[1];
u3(0.446433018237275,1.64723781591533,-3.37983322372519) q[8];
u3(1.56777954743848,2.98362135465698,-2.99774802559812) q[4];
cx q[4],q[8];
u1(0.954774972117809) q[8];
u3(-0.0591515945072143,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.40618973744745,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.65989654676963,-0.0831157288891009,-1.28007061758073) q[8];
u3(0.491871837678752,1.41627970543470,0.487999836597600) q[4];
u3(2.06332593676464,1.50715921969515,-2.73806007988476) q[3];
u3(0.770676315177051,-2.41435969590894,1.99748953527958) q[13];
cx q[13],q[3];
u1(0.143805363085886) q[3];
u3(-1.40292144915769,0.0,0.0) q[13];
cx q[3],q[13];
u3(2.01298557620524,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.43752955942808,1.95160631948488,0.0727782546045457) q[3];
u3(2.03851835209827,-5.13611989078888,0.777791012403719) q[13];
u3(2.46274735164768,-3.54821394466961,1.74475203973510) q[5];
u3(0.884436689045162,3.51888400829472,-1.24265574242027) q[11];
cx q[11],q[5];
u1(4.33640181675442) q[5];
u3(-3.87368931895164,0.0,0.0) q[11];
cx q[5],q[11];
u3(-0.434663756587793,0.0,0.0) q[11];
cx q[11],q[5];
u3(0.259530873214116,1.03235400135900,1.63811011479449) q[5];
u3(1.45550902605726,2.10165400876273,2.76261285245300) q[11];
u3(1.35109581413253,-2.99877915939213,0.610306652476790) q[2];
u3(0.758785249545465,-3.04161077188002,0.714579009053065) q[12];
cx q[12],q[2];
u1(3.94528520674024) q[2];
u3(-1.21249981067156,0.0,0.0) q[12];
cx q[2],q[12];
u3(1.65164099216132,0.0,0.0) q[12];
cx q[12],q[2];
u3(1.56395720146324,-4.31445808187632,1.41702906282882) q[2];
u3(0.746648653985354,2.00839382083890,-3.49343302502086) q[12];
u3(0.791152832816254,0.539758445550967,0.611660650642301) q[6];
u3(1.85209704754799,-0.585529263849286,-1.70355723074110) q[10];
cx q[10],q[6];
u1(2.83510391976181) q[6];
u3(-2.49477633773817,0.0,0.0) q[10];
cx q[6],q[10];
u3(0.964652850073072,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.801731702412580,-0.866391364508862,1.61533216689381) q[6];
u3(1.11129554320078,0.675375552559433,1.66854336605742) q[10];
u3(1.95044544183146,1.69704792345462,-0.0981421953782716) q[11];
u3(0.883216104362397,0.516687123220111,-4.80738874211440) q[10];
cx q[10],q[11];
u1(2.00159289725572) q[11];
u3(-2.64012076164055,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.35189699214426,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.44340517586073,-3.29907197266343,2.80215249952111) q[11];
u3(0.803206816248323,-3.16655394445688,1.11000708751436) q[10];
u3(1.50901653532027,3.27355608253424,-0.370050500148565) q[8];
u3(1.06120673718374,1.36559645131975,-1.53340624774026) q[2];
cx q[2],q[8];
u1(2.51515069054090) q[8];
u3(-1.65547367552622,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.93787751345541,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.11289631369154,-1.56442772462497,1.22816173820052) q[8];
u3(2.39179204372888,-1.81464872153738,-1.48339771817036) q[2];
u3(2.57888431532325,2.26391548044418,0.537668523809111) q[5];
u3(1.59234736914382,0.129932489298457,-4.24087982152687) q[9];
cx q[9],q[5];
u1(1.45890348815048) q[5];
u3(-0.159115793995264,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.36213130040883,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.51980916791601,0.865875876638423,-3.22993092036021) q[5];
u3(2.82210162267480,-0.560330213405311,1.05418309701154) q[9];
u3(2.45954018667592,1.21432539987866,-2.61344579206274) q[3];
u3(1.68319776315790,-2.05062571893990,1.82114114186424) q[12];
cx q[12],q[3];
u1(2.85188856099278) q[3];
u3(-1.57744583164952,0.0,0.0) q[12];
cx q[3],q[12];
u3(-0.0468310551827995,0.0,0.0) q[12];
cx q[12],q[3];
u3(0.450631675194846,-0.392169199682571,-0.919375451954195) q[3];
u3(1.73522204015098,1.82975966309342,3.86274707054216) q[12];
u3(1.11153970932665,1.22908981601794,1.15097655093808) q[13];
u3(1.25885484894017,-0.751294402441223,-3.15443905644929) q[7];
cx q[7],q[13];
u1(3.57208346926216) q[13];
u3(-1.43461965535983,0.0,0.0) q[7];
cx q[13],q[7];
u3(1.90326672659890,0.0,0.0) q[7];
cx q[7],q[13];
u3(1.64300965005736,-3.29295658349295,1.91049779437143) q[13];
u3(2.04470869444857,-1.73571493778094,-1.34084610669691) q[7];
u3(2.19586557493007,1.27021699334788,-1.14506188486991) q[6];
u3(1.76651201996868,-0.442930949380863,-3.50220547126838) q[1];
cx q[1],q[6];
u1(3.15018238361175) q[6];
u3(-1.76562312468150,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.815149507434225,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.59386194509117,0.166878531833680,1.93925583404613) q[6];
u3(0.732125340110507,-3.75425573629718,2.39597044596912) q[1];
u3(0.849959737349773,1.81228670507138,-1.68365862041878) q[4];
u3(0.815599585254556,0.898017810569248,-3.13438037258905) q[0];
cx q[0],q[4];
u1(0.537057689723444) q[4];
u3(-1.14471311986132,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.0528562758218274,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.707572833065231,-1.29341256165381,-1.19003194621077) q[4];
u3(1.25794002309458,-2.54539939597350,1.64097655572422) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13];
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
