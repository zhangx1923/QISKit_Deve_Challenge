OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.69180329095720,-1.03776547644061,-0.978417870943987) q[8];
u3(2.37509290181375,1.95426833668931,-4.02160823772686) q[3];
cx q[3],q[8];
u1(1.04850749473339) q[8];
u3(-1.40463726605385,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.90886528746788,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.38847149064494,2.91280265932587,0.867970680049503) q[8];
u3(2.45697003859764,-2.44459429719999,-2.25778371408184) q[3];
u3(0.776169672958378,0.514566259732120,-2.95520243494528) q[6];
u3(1.10273445918679,2.76508114900862,-3.06711562698039) q[0];
cx q[0],q[6];
u1(3.00641440260857) q[6];
u3(-1.80312583867650,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.734166744032595,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.57039167902039,-3.62183684686405,2.39002871005466) q[6];
u3(1.52997874764150,-0.456672451714918,-2.49790162448924) q[0];
u3(0.706677592259499,-1.21467040468156,2.23600407546909) q[7];
u3(0.0813446370338149,1.02739459641168,-2.92579190351335) q[5];
cx q[5],q[7];
u1(1.63634704262005) q[7];
u3(-0.232819942170260,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.31718031956776,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.09937077473390,-0.501052969459049,1.32363786007267) q[7];
u3(1.90409501819163,-0.0192756672026995,1.89674253199464) q[5];
u3(1.78706667342521,2.84601048760083,-0.857209430133010) q[1];
u3(1.54316516795519,1.06222785224024,-1.41111447526301) q[4];
cx q[4],q[1];
u1(1.48835057922967) q[1];
u3(0.477915198127038,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.719447636443569,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.63044634441094,0.584466008937323,-1.56816551514604) q[1];
u3(0.647849573100787,-1.20396388008880,-3.59732864016327) q[4];
u3(1.66890570592200,-0.0918417459834312,1.12224864595264) q[2];
u3(1.70488172217803,-1.35403654743706,-2.12655123635388) q[9];
cx q[9],q[2];
u1(1.55571859860695) q[2];
u3(-2.88548411168132,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.03587478205609,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.24667191641716,-2.40845816988023,2.93349724920066) q[2];
u3(2.36480771223389,-1.88396490045477,4.27123014230622) q[9];
u3(1.67711492255647,-1.18089874365121,1.43273576740086) q[6];
u3(0.766344790446919,-3.67822721600572,0.145743816492815) q[9];
cx q[9],q[6];
u1(1.79025305683988) q[6];
u3(-2.78732829724270,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.929927164013186,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.991959323510421,1.17808303897423,0.474929146627838) q[6];
u3(2.10955371965054,-0.456057024280175,-3.56910378768029) q[9];
u3(2.67223664260414,2.35343691677978,-3.66306131523797) q[7];
u3(1.24276255371431,0.458642243228465,0.764670574142259) q[5];
cx q[5],q[7];
u1(1.89212588593087) q[7];
u3(-2.69338181086930,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.18998329687494,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.925058102043254,-1.16074326123476,-1.40850451646713) q[7];
u3(1.04398380230723,-0.0607542556802199,-3.84105506796891) q[5];
u3(1.29609985795948,2.73038756614316,-3.52331321817297) q[3];
u3(1.76620606838288,-2.81752324617945,2.78548127096257) q[8];
cx q[8],q[3];
u1(1.74435634203781) q[3];
u3(-2.49274387640682,0.0,0.0) q[8];
cx q[3],q[8];
u3(3.52750519714156,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.67523978689659,1.06919180767914,-3.09351709365140) q[3];
u3(1.40861366202522,1.33901445908053,-3.16955701265144) q[8];
u3(2.69295980679857,-2.73319301101426,-0.262690534018060) q[0];
u3(2.44486955948718,-2.75557287669739,-1.73294123945135) q[1];
cx q[1],q[0];
u1(0.497885468602200) q[0];
u3(-1.20801295282627,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.55933733703301,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.74867183283485,1.80229477739546,-1.64086674346076) q[0];
u3(2.89301335105153,-2.73902957357405,-2.72416881488631) q[1];
u3(0.525997511889967,1.11089157316145,-0.278379769387967) q[2];
u3(1.18682477039712,-0.472272679447710,-2.41302812746129) q[4];
cx q[4],q[2];
u1(1.49485725332175) q[2];
u3(-2.56135806132475,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.35303889662203,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.44801963691642,2.05702796360323,1.92989863802602) q[2];
u3(2.25289744459425,0.237299620116358,-2.75929798347233) q[4];
u3(0.766870324688791,3.21594789476058,-2.83462027604203) q[5];
u3(1.22697431757348,-0.0176756340580648,-2.11934955285520) q[9];
cx q[9],q[5];
u1(1.49677502279637) q[5];
u3(-3.06380632902836,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.653122590561263,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.419183135248751,-0.964694652684437,0.508363335987346) q[5];
u3(2.11273390399712,1.09481237367783,-4.60764297832429) q[9];
u3(1.20538032538054,-1.44992265535222,-1.03574522404683) q[2];
u3(2.31922355903737,-2.54534805143718,-0.0968990160112575) q[8];
cx q[8],q[2];
u1(0.484415286878419) q[2];
u3(-0.886708301870023,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.41779335000967,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.23749554220357,0.638182794629326,-0.863449021732985) q[2];
u3(0.872808961777994,-2.28014255227465,-3.92225862722143) q[8];
u3(1.26947269408694,1.06217968043038,-1.14299870965403) q[7];
u3(1.79478142109867,-4.91779974107875,0.875973639126105) q[6];
cx q[6],q[7];
u1(0.817581135584692) q[7];
u3(-1.36016474782656,0.0,0.0) q[6];
cx q[7],q[6];
u3(3.12007379491861,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.67564174568621,1.68297721914552,-2.03597970723296) q[7];
u3(2.38063125880218,0.106428505723990,4.07558032375518) q[6];
u3(0.797162590371687,3.54530476686777,-2.14242168014436) q[4];
u3(2.28518683492206,2.40761227664071,-1.67381642254969) q[0];
cx q[0],q[4];
u1(1.55050446113990) q[4];
u3(-0.392911500073501,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.42453710802490,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.22490856103698,-2.24249624146667,2.49374612781889) q[4];
u3(1.38211264621205,-2.09852576087777,2.29109628418228) q[0];
u3(1.61631517454640,0.351345542104269,-1.10645064483881) q[3];
u3(0.327788327110233,-3.90389234108503,0.618335968819673) q[1];
cx q[1],q[3];
u1(3.83829067753807) q[3];
u3(-3.34849334375626,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.637157348035267,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.16601887719902,-2.82347404442936,0.443317615912189) q[3];
u3(0.495611068085938,0.124434625822704,6.11775707045554) q[1];
u3(2.45282478808570,1.54930385387609,1.47592604419722) q[8];
u3(0.835507208397441,-3.63862673054103,-1.04634285539393) q[3];
cx q[3],q[8];
u1(2.28945768981910) q[8];
u3(-2.63462600078566,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.91749008717866,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.566206130991415,1.78552054672274,0.0582469785128341) q[8];
u3(1.79375554359888,1.63406772629687,-1.89416268461308) q[3];
u3(2.66018771946965,-4.01972716366762,0.957099102812107) q[0];
u3(0.347393764546757,2.88853293152015,-1.41269849952449) q[9];
cx q[9],q[0];
u1(-0.922833860835465) q[0];
u3(0.325308369769101,0.0,0.0) q[9];
cx q[0],q[9];
u3(3.51455638280065,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.98588916550821,-2.15668918325574,2.71995007068395) q[0];
u3(2.69117338310415,-0.367828783305646,-5.74441902020391) q[9];
u3(1.76260270644513,0.737678089976038,0.449006036349785) q[7];
u3(2.12340931545110,-1.50929871127145,-1.14669586514409) q[4];
cx q[4],q[7];
u1(1.31601747250281) q[7];
u3(-0.348242591323686,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.32267392273595,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.84464817401391,-0.505069266816431,1.19033228927232) q[7];
u3(1.57166969334367,-4.34923679718928,1.70477427430199) q[4];
u3(1.69668902212689,-0.638164174586716,-1.01976529483435) q[2];
u3(1.61884658632184,-3.81263411603106,0.899875284990490) q[6];
cx q[6],q[2];
u1(2.65752550200171) q[2];
u3(-2.88928910799131,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.43644482662860,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.11270956159787,-0.237422379237873,2.21162557388771) q[2];
u3(1.02921188769406,-1.58244047662143,3.55574394682588) q[6];
u3(2.37872527022179,1.41755675417387,-1.21930316652542) q[5];
u3(1.58256596008998,-4.39405948124928,1.33054505124034) q[1];
cx q[1],q[5];
u1(0.821398119206425) q[5];
u3(-3.47323367975310,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.55219528259292,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.85075142426539,-0.594331230517417,-0.863893041378603) q[5];
u3(2.39310151991748,1.62244895465238,2.26406833274071) q[1];
u3(1.48399319729530,-0.229853676833656,-1.39842022758363) q[0];
u3(1.84548898628954,-4.32372582543228,0.685878347241807) q[7];
cx q[7],q[0];
u1(2.33118594022122) q[0];
u3(-3.00460945945197,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.916750448539392,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.96445750931488,3.56901967581310,-2.31724939458579) q[0];
u3(0.398071791821836,5.36015174143125,0.849682188937300) q[7];
u3(1.86714119051467,-0.325570419324411,-2.03998186701986) q[4];
u3(1.76847154264044,-0.199655234768070,-5.53486128030891) q[9];
cx q[9],q[4];
u1(2.76538915583769) q[4];
u3(-2.10044332004477,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.231546660304147,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.49238831197340,-3.43448657307485,2.74424020145817) q[4];
u3(1.58767663339151,3.53216298057514,-2.45913549828545) q[9];
u3(0.425088281697852,-2.74704639569707,2.27458249435976) q[3];
u3(1.04372622660875,0.423757326718019,-1.90716650654439) q[8];
cx q[8],q[3];
u1(1.55866144428273) q[3];
u3(0.126626610064753,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.15649961738557,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.392915849144728,-4.02936646053376,0.321883552877770) q[3];
u3(2.18176620520341,4.34954365344813,1.03218634723201) q[8];
u3(1.82982948880640,-2.27045163817935,3.94907650801206) q[2];
u3(0.327678787655122,2.18633570106108,-0.849582514215651) q[1];
cx q[1],q[2];
u1(2.25018744784455) q[2];
u3(-2.98941819015860,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.709251345373847,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.608613546780712,-3.00273899653130,-0.0373711516645843) q[2];
u3(1.29686350156203,-0.815939777740073,4.04060475111512) q[1];
u3(2.28252382911530,4.19136396726772,-1.23029235010229) q[6];
u3(1.35560051383534,3.05506769101889,0.377009098065791) q[5];
cx q[5],q[6];
u1(-0.183545901408134) q[6];
u3(-1.88240778663414,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.710728241446495,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.88407257853641,-2.35434226193516,1.04767616719830) q[6];
u3(1.25823419382170,2.96979466131793,2.21270709436040) q[5];
u3(1.85905075872216,-0.767164449001038,-2.00250827334780) q[0];
u3(1.46223609922039,-4.91722028242059,0.807144236122743) q[2];
cx q[2],q[0];
u1(-0.375238700169555) q[0];
u3(-1.66606355432530,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.789862316146891,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.75781031157388,1.21649256090565,1.26810051781167) q[0];
u3(0.867586665682038,0.114630616087072,3.60034593024345) q[2];
u3(2.19897505676788,0.880931224013387,1.89020295529159) q[9];
u3(1.92198542861585,-1.12241980444707,-1.05610239957235) q[6];
cx q[6],q[9];
u1(3.51374452593313) q[9];
u3(-1.58744305611686,0.0,0.0) q[6];
cx q[9],q[6];
u3(2.35491676414923,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.61249432354434,-2.39734951664633,-0.873742776171790) q[9];
u3(2.09116290378438,4.32024551658570,-0.478935574216398) q[6];
u3(1.37599210248634,-0.456106463195278,-2.19040276959291) q[3];
u3(0.989050710595736,1.14777412908541,-4.13773223739611) q[1];
cx q[1],q[3];
u1(0.901199287373176) q[3];
u3(-1.43750721814721,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.0495284222928662,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.30880011892021,1.81550454943145,-2.75247722542161) q[3];
u3(1.49864715114747,1.54339942875892,-1.64518519656845) q[1];
u3(0.852403198914644,1.13478750638660,-0.373307278799293) q[7];
u3(2.06394135658164,-1.27941466120653,-3.92430628626834) q[5];
cx q[5],q[7];
u1(0.342844977352297) q[7];
u3(-1.49464984583602,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.73248122884434,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.14241884630815,1.37581355071638,1.78268632179269) q[7];
u3(1.35411497865941,-0.660461572579184,3.53594480707613) q[5];
u3(1.11907560932549,-1.12686404419704,-1.24488242497447) q[4];
u3(1.63621764210971,-4.34679562576201,1.19435029914851) q[8];
cx q[8],q[4];
u1(0.181144205157718) q[4];
u3(-1.54498541732537,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.47385478602458,0.0,0.0) q[8];
cx q[8],q[4];
u3(0.376239367262805,-2.25013785971980,3.21132985893830) q[4];
u3(0.983904128963852,1.70614225833885,-1.45902512336007) q[8];
u3(0.890842473887848,1.01125989190629,-0.169319647233459) q[6];
u3(0.988656474306825,-0.826842228902146,-0.720738775523334) q[0];
cx q[0],q[6];
u1(1.38388497194888) q[6];
u3(-3.53608803982622,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.98338606021728,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.82842076681823,0.00786344119311577,2.60371212413565) q[6];
u3(1.26502997624151,3.22211271840336,-0.521227245367762) q[0];
u3(2.16221329562042,3.75538456138172,-1.90654086333986) q[9];
u3(0.894769696016178,-1.38279146579785,2.94115865308744) q[2];
cx q[2],q[9];
u1(0.200560109587177) q[9];
u3(-1.28779301531585,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.23419362176656,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.08103511755798,0.890932640903919,-3.16407533622777) q[9];
u3(1.47390667416723,-4.33059095403097,1.16436236570917) q[2];
u3(1.52472461002465,1.11297624575006,0.679084013805781) q[4];
u3(0.482477307781069,-0.549054072810958,-2.55407918623031) q[1];
cx q[1],q[4];
u1(2.49038517668700) q[4];
u3(-1.82514547704882,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.990142877161575,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.674095099651838,-1.15840541584841,2.80675645280060) q[4];
u3(2.56489575132470,-2.02659482070554,1.10697803366616) q[1];
u3(1.30542260735014,1.97067480227706,-1.81634250566886) q[7];
u3(1.01962865517367,1.32951572795023,-1.67957150743679) q[8];
cx q[8],q[7];
u1(-0.0157655943394317) q[7];
u3(-1.77558366942988,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.01095672672848,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.618943524888117,-4.15973718925865,1.83230309652063) q[7];
u3(1.70371502701387,5.61193099852644,-0.0679691214543383) q[8];
u3(2.49326938163323,-0.366571821918773,-0.704708985413319) q[5];
u3(0.721679100742240,-3.04095327031967,-1.06265796021682) q[3];
cx q[3],q[5];
u1(-0.292262270589166) q[5];
u3(-1.83161518371988,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.961609583730749,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.733352471547140,0.285841923341404,0.418582230587130) q[5];
u3(1.75389539821078,1.83394333794036,-3.04388666293633) q[3];
u3(0.792660928506895,1.71561218190292,-2.87341288818037) q[3];
u3(2.15261392471428,-2.92987109403685,2.88966772815213) q[8];
cx q[8],q[3];
u1(2.48835359385632) q[3];
u3(-2.70813274302671,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.72092087666502,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.63287046491073,2.43123365050485,-0.428961758761860) q[3];
u3(1.40813764729280,1.47641818691939,-3.13276718076051) q[8];
u3(0.715252739492527,-0.506551125014026,0.390929128350505) q[6];
u3(0.612185468133351,-2.85851772698641,0.260064093590968) q[5];
cx q[5],q[6];
u1(-0.108992399968818) q[6];
u3(-2.24107437177587,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.38055714389321,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.37525199659503,-1.16459857229396,-0.718058568055197) q[6];
u3(0.896393920691387,-1.09171012769647,-4.63153437604960) q[5];
u3(0.240012858201827,2.70699802255394,-2.69219576866143) q[2];
u3(1.13243915816039,-3.41335140628558,2.07291345431909) q[1];
cx q[1],q[2];
u1(1.32189915230237) q[2];
u3(-3.05431308242891,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.47322118539844,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.30573350462764,2.72091606849535,-3.24366548717447) q[2];
u3(2.08662477413370,0.786627094804945,-0.628705932797316) q[1];
u3(1.47296801101041,2.20793348860868,-0.735365053172012) q[7];
u3(2.69936085756648,0.747919068222699,-3.01666205453524) q[9];
cx q[9],q[7];
u1(-1.09756685575822) q[7];
u3(0.544777139839765,0.0,0.0) q[9];
cx q[7],q[9];
u3(3.37143718847435,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.42703168212256,-3.02821856900181,1.81986819314196) q[7];
u3(2.62287673882558,-1.35201690394237,-1.66761900559301) q[9];
u3(0.824251890867846,0.169061822181069,-1.51393273115574) q[4];
u3(0.830204780684977,-1.26759857438604,-0.0915364139176197) q[0];
cx q[0],q[4];
u1(1.48629652500518) q[4];
u3(-3.21582419781673,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.858388612480145,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.47684763342993,0.301933140319255,-1.64883415559222) q[4];
u3(1.02965300047690,-4.47523835291876,1.09844907674748) q[0];
u3(3.03543789588795,1.99867190760733,0.915476268342748) q[6];
u3(1.30935411763003,-3.88476433787969,-0.457144559406048) q[3];
cx q[3],q[6];
u1(2.89880709833457) q[6];
u3(-2.12802299254939,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.13818914217692,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.53889505244393,0.169531966155038,0.607680364381442) q[6];
u3(0.407891669130315,1.29758703611261,2.01701665374616) q[3];
u3(1.41627571194177,-1.00015776270747,2.73862647524652) q[5];
u3(1.24670174427896,-1.28440366597300,-2.28722677184625) q[8];
cx q[8],q[5];
u1(1.51630798287375) q[5];
u3(-0.851119901628324,0.0,0.0) q[8];
cx q[5],q[8];
u3(2.90658772054316,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.40859571244917,-3.80850691477803,2.16075351293673) q[5];
u3(0.192739059654005,3.35184209216849,1.86669517817272) q[8];
u3(1.01968590621477,0.639649997845395,1.01430998257052) q[4];
u3(1.55929009844251,-0.598470582553562,-2.46523363881207) q[7];
cx q[7],q[4];
u1(2.30474924140970) q[4];
u3(0.0428715119567438,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.40962911969473,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.587938193047683,1.28566505941635,-3.51148481506316) q[4];
u3(1.39513757985986,2.99640592236806,-1.00592136236661) q[7];
u3(1.18327196791724,2.84564417812414,-0.971335862537701) q[0];
u3(1.04277905641826,1.18634739807941,-1.47042677876837) q[9];
cx q[9],q[0];
u1(1.32096694347238) q[0];
u3(-0.819633542836333,0.0,0.0) q[9];
cx q[0],q[9];
u3(-0.400424106574407,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.39013001902939,-1.76648688262464,3.36130038721750) q[0];
u3(2.10464076034352,0.479590076887185,-5.02468397691872) q[9];
u3(1.78027750558499,1.89950074078775,0.0331566927771663) q[1];
u3(2.58378910387422,0.470925645915618,-3.09810965442889) q[2];
cx q[2],q[1];
u1(2.21264018671422) q[1];
u3(0.364951744798208,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.28158171577405,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.16056082392403,1.89271123359906,1.26557961124423) q[1];
u3(2.28965765302680,-0.0597257432415176,-0.508566799618353) q[2];
u3(0.335667776703175,2.45951050047416,-3.44910853482387) q[4];
u3(0.993138965527681,0.153401354688237,-1.76902392918145) q[2];
cx q[2],q[4];
u1(2.07200885774718) q[4];
u3(-2.71581779241324,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.220445918303877,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.92039400640225,-3.99479196056020,1.40008497514797) q[4];
u3(0.622222522911926,-1.24653895425287,-3.01347996002854) q[2];
u3(0.303311073588936,-2.79365788008388,2.84455045442695) q[3];
u3(0.581266206857636,2.23488553097280,-3.25380206069262) q[6];
cx q[6],q[3];
u1(-1.04954808773997) q[3];
u3(-0.0612311547010695,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.59531235387035,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.871080895559287,2.72942895900596,-2.78011170751115) q[3];
u3(0.983888616980098,-3.53194414197231,2.17587948970775) q[6];
u3(2.56772272852926,0.618750194660141,-1.78193874802425) q[9];
u3(2.25678062734562,5.17844662818339,0.493504392393957) q[8];
cx q[8],q[9];
u1(1.22876084286665) q[9];
u3(-3.42259413724947,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.77568238820105,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.29558560490872,-1.03474247778365,2.68444841844881) q[9];
u3(1.27577315612179,-0.993200585380783,4.09179136650012) q[8];
u3(1.37954754665245,1.14524671399707,-0.829386712020451) q[5];
u3(0.456911006836285,-1.47952234075962,-0.903100406426429) q[7];
cx q[7],q[5];
u1(3.15358429877293) q[5];
u3(-1.79632943891702,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.732724138926298,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.334745498195044,0.709196002529370,-3.91574334433348) q[5];
u3(2.39066630284459,2.28024371342609,0.603804207574638) q[7];
u3(1.61208103756919,1.37539049564757,-0.941304082387052) q[1];
u3(1.55920187855416,-4.60256750937416,1.12199924501746) q[0];
cx q[0],q[1];
u1(3.06983062342831) q[1];
u3(-2.14439285344566,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.25052683473605,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.970452110815483,-0.140968038424566,1.78627826374282) q[1];
u3(2.20288401882967,1.21147642808813,2.16203646877555) q[0];
u3(0.946353828696350,-1.55469191216085,0.0118347283258661) q[2];
u3(1.39887698642628,-3.47624598114985,-0.401171231093193) q[5];
cx q[5],q[2];
u1(2.32160789900862) q[2];
u3(-1.56754266314843,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.57563879747273,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.03544875991875,-0.997592947770406,2.66041033879984) q[2];
u3(1.61252976911970,-0.767692202242150,-4.61167359783740) q[5];
u3(1.26215974393186,3.85102041650836,-0.791976006816523) q[8];
u3(0.841227786435409,1.56718394004480,-1.25645398680893) q[3];
cx q[3],q[8];
u1(0.161309614771063) q[8];
u3(-1.98061087524878,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.985277562180876,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.19096032989030,-1.48588165666940,0.485582269775749) q[8];
u3(1.02786188060110,-4.60259021461848,-0.401917896200392) q[3];
u3(1.75519349830536,0.999639347092712,-1.38816554611920) q[7];
u3(1.54414231122795,-1.11512158211253,-3.25053081995096) q[9];
cx q[9],q[7];
u1(-0.0283650238799915) q[7];
u3(1.14479206861539,0.0,0.0) q[9];
cx q[7],q[9];
u3(3.45859088800306,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.865901087550346,0.687720663294402,0.952130856551958) q[7];
u3(1.61585568260202,3.40514352699937,-1.21252716552271) q[9];
u3(0.876122308866985,0.826455732921783,-0.660846647109292) q[0];
u3(1.01954798618795,-0.624812989602212,-0.238826558993248) q[1];
cx q[1],q[0];
u1(2.06561516565283) q[0];
u3(-2.44645953989892,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.0913577190250037,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.880784725121972,1.34907025447482,-0.775394961655397) q[0];
u3(2.01223425897716,-0.848214391920093,-5.22346645012866) q[1];
u3(1.65358008596556,1.86869354420183,0.298832767248985) q[4];
u3(0.566599376142184,-1.54515389067159,-1.69164643653196) q[6];
cx q[6],q[4];
u1(1.48381738146250) q[4];
u3(-3.13414941125572,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.38613836419522,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.02904730888244,-1.41771881136035,-3.20054725424687) q[4];
u3(2.21847124912937,4.95262870541322,0.946667082268989) q[6];
u3(2.81507511221108,2.69696585370814,-0.268131658265194) q[4];
u3(2.54680916613058,1.07961238566909,-2.23316493371722) q[7];
cx q[7],q[4];
u1(0.958196245250114) q[4];
u3(-0.354291997065222,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.63940374330936,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.00093388626024,-3.38534525449707,1.30920069493068) q[4];
u3(0.820184938658647,-3.89885365614218,-1.52748101455893) q[7];
u3(1.74964194230882,1.88139690534685,-0.0104376695842849) q[0];
u3(0.715181542621281,0.883841812359529,-4.49579890263073) q[3];
cx q[3],q[0];
u1(-0.171267836783928) q[0];
u3(-1.56850633862080,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.19716063417163,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.07073022280227,3.74412776893076,-0.780319681529537) q[0];
u3(2.75627577254252,2.12055081460659,-2.07275195335565) q[3];
u3(0.584125713465160,1.84867531084705,-0.0474299989227325) q[8];
u3(1.19046913984873,0.280091054488036,-3.77104680632243) q[6];
cx q[6],q[8];
u1(0.00416608779105321) q[8];
u3(-1.19033009192448,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.56567573797210,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.87385593843984,-2.17823613966634,3.54817843313683) q[8];
u3(2.61416473182765,-1.92766525677635,2.65939984008387) q[6];
u3(1.79269135646670,-0.142591048161434,-1.26793626435668) q[2];
u3(1.19458891552969,0.229637540480989,-3.53242609864579) q[5];
cx q[5],q[2];
u1(3.51372140254869) q[2];
u3(-1.19361167335060,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.22536200331411,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.788281876801006,1.51357273713940,-3.75272206610720) q[2];
u3(1.65127088223874,0.366165899501207,0.122695926187351) q[5];
u3(2.01033402980165,1.51793923489215,-3.16152249353722) q[9];
u3(1.41319074750590,-2.12498629140358,2.01122868642105) q[1];
cx q[1],q[9];
u1(0.643703218318670) q[9];
u3(-1.56477698127607,0.0,0.0) q[1];
cx q[9],q[1];
u3(-0.414483549385495,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.59740515339276,-0.0972907441734345,-0.00338156267885412) q[9];
u3(0.657007735047451,-4.13615638253139,-1.73306003167865) q[1];
u3(1.60964894724809,3.12525441485698,-1.32747426697193) q[2];
u3(1.05533589096998,1.75599042913525,-2.11572456777253) q[6];
cx q[6],q[2];
u1(3.73452201018319) q[2];
u3(-1.21177563345218,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.87940697712591,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.65412346666414,2.30590102816348,-0.272165238642434) q[2];
u3(1.89424107212599,-0.0556187761244638,4.38316211206035) q[6];
u3(0.230192984829154,2.36990863048709,-1.70914171285554) q[5];
u3(0.466949408405959,-0.00970483029573244,-1.07919963382204) q[4];
cx q[4],q[5];
u1(0.925762875181853) q[5];
u3(-3.14795777565378,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.87714510127755,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.19446561099436,0.389836440357899,-1.17316392726809) q[5];
u3(2.38485587459430,-1.30569826639521,-2.13824929537261) q[4];
u3(1.24138204017715,1.41943409148620,0.192811357751891) q[1];
u3(2.02217338898831,0.814625845147629,-2.56967919247651) q[7];
cx q[7],q[1];
u1(1.50889211864523) q[1];
u3(-0.489217946539724,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.50537568722602,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.04777926178838,0.165742597472416,-0.980849141191067) q[1];
u3(1.57369840551617,-2.64110503794134,2.35589308545151) q[7];
u3(1.92609631772435,1.14533880001388,-3.50660486067600) q[0];
u3(1.32182376624327,-1.92502138419617,2.47648956370086) q[8];
cx q[8],q[0];
u1(2.67507005223443) q[0];
u3(-2.22819573830601,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.0881049030606889,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.96932027878864,0.496269681843942,1.37286087216272) q[0];
u3(1.22423110165848,-1.71699905099129,-3.13386003490128) q[8];
u3(2.16481189518152,-1.74189608742467,-0.669227805676159) q[9];
u3(2.00244578596681,-4.07923984815926,-0.666193560648897) q[3];
cx q[3],q[9];
u1(2.41716042802329) q[9];
u3(-1.90307440615779,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.48650824146699,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.45326178162923,0.935172713507196,2.09036107170478) q[9];
u3(2.75872248603871,4.50366514137894,1.50631267246922) q[3];
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
