OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(0.534704696819449,2.55517421823062,-2.62991754371664) q[2];
u3(0.638030423432675,-0.202472908758145,-1.87999712079278) q[0];
cx q[0],q[2];
u1(0.0123389086452808) q[2];
u3(-1.86859373801974,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.56871401257557,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.461093667834015,1.58454748098247,-2.72975225140421) q[2];
u3(1.31787992005462,-1.14245462728261,-3.14961080763034) q[0];
u3(2.36585539796282,0.0349979709394596,0.622321319229548) q[3];
u3(1.55483774334348,-2.57983449773272,-0.636513327249879) q[1];
cx q[1],q[3];
u1(2.07778357899113) q[3];
u3(-1.51999499621401,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.83364811905948,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.15303912896734,2.91247536234335,-1.43956566513852) q[3];
u3(1.81277274444178,0.723159413567048,-1.37797308288779) q[1];
u3(0.845681063748727,0.0525264200773330,-2.19235522028609) q[2];
u3(1.01951531371602,0.599072786916418,-4.92535648285192) q[1];
cx q[1],q[2];
u1(1.86468942062664) q[2];
u3(-2.81359867032253,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.802560033790128,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.60998882773964,-0.391999089208051,0.824315906549448) q[2];
u3(0.827003103411047,4.49788365122073,-1.07766161402682) q[1];
u3(1.43397510382902,0.0291201238091378,-0.522652645595955) q[3];
u3(1.63481147062138,-3.10504601348332,0.990075476056618) q[0];
cx q[0],q[3];
u1(3.18059973823124) q[3];
u3(-4.02840419776310,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.587666998986824,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.91824866165249,-3.17459079573642,2.25574947952922) q[3];
u3(2.91478564529894,0.915498609665656,3.28272607923497) q[0];
u3(0.280805531320914,-1.54811921847650,1.22670265597147) q[0];
u3(0.826519107541062,-2.64621680720722,2.10308533667929) q[3];
cx q[3],q[0];
u1(1.50195725931669) q[0];
u3(-0.985448101482420,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.545275840571761,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.68435710796558,0.237256484866072,-2.17212921052420) q[0];
u3(0.531449565212384,-3.67196679231292,0.633568186315589) q[3];
u3(1.45127257800669,-1.99062711821121,-0.501905177846252) q[1];
u3(0.490987517950954,-3.41892621822970,-0.143879751053128) q[2];
cx q[2],q[1];
u1(2.71985928368363) q[1];
u3(-1.91517678909851,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.647273006878351,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.62966707513755,3.15116626144033,-2.55219490561765) q[1];
u3(0.971894177145932,1.60160880539789,-1.36064121275872) q[2];
u3(2.03706971702393,2.98996649996264,-1.15569519720230) q[1];
u3(2.06442476418396,0.821044814160965,-0.934146793245374) q[3];
cx q[3],q[1];
u1(1.92394215084325) q[1];
u3(-2.17687648175248,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.913921475608981,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.01964141868390,-2.45529733393066,0.208138315095427) q[1];
u3(0.625577362095166,3.88694781225197,1.86413789187734) q[3];
u3(2.03407159761334,2.59935394904106,-2.24095274208457) q[0];
u3(2.41289683745639,1.24724883827848,-1.75521353174292) q[2];
cx q[2],q[0];
u1(2.70435455399198) q[0];
u3(-2.23608130546608,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.44388830983033,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.40272322163762,-2.55235659743531,-1.55286788522635) q[0];
u3(2.00875525952157,3.35805062075495,-1.88963236741370) q[2];
u3(2.82719364006839,0.703064837857088,2.24442311163738) q[0];
u3(1.88553541970675,3.40737223397273,2.60090586499007) q[2];
cx q[2],q[0];
u1(3.30612556494173) q[0];
u3(-1.56968160897036,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.07566769232255,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.69348746444834,-0.705759126848139,0.746372263345020) q[0];
u3(2.26989973292485,2.67386873580437,-3.09841257821063) q[2];
u3(2.03731997534839,2.40218334367597,-1.85061166946292) q[1];
u3(1.41051763680540,1.54304167899816,-2.78511278553227) q[3];
cx q[3],q[1];
u1(0.932198085780615) q[1];
u3(-3.44026911317653,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.32162690605207,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.72806960916320,4.04840232753536,-1.15873329631422) q[1];
u3(1.29750675059442,-4.19516051952397,0.436262296731713) q[3];
u3(1.46594951609045,2.29251690827128,-3.80796839820587) q[2];
u3(1.15216548861746,-2.30343987050741,2.88497683400325) q[1];
cx q[1],q[2];
u1(0.701622785032135) q[2];
u3(-1.45021987943248,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.383977611407164,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.00326472810872,1.31447293402518,-4.66203502980588) q[2];
u3(2.01317834753865,3.43121431839542,-0.708809985345776) q[1];
u3(2.42923916512564,-1.73503266178126,0.970058288592283) q[3];
u3(2.97785772971381,1.06675503500010,1.74715893890734) q[0];
cx q[0],q[3];
u1(0.916356375157647) q[3];
u3(-3.17850159242692,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.80633753162443,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.739159617539421,-0.565068721162384,2.41866850977684) q[3];
u3(0.556064040504416,0.759287274279278,2.95583817476589) q[0];
u3(1.86311087046883,-0.419868326668319,0.696043667992712) q[0];
u3(2.00433358701463,-0.666484592346222,-2.06924280190780) q[2];
cx q[2],q[0];
u1(3.01539376531651) q[0];
u3(-1.63325623150671,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.29670198663284,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.56227443660542,2.22414901028965,-2.75827356086485) q[0];
u3(0.507357866560242,-3.95079700038164,1.61145753352100) q[2];
u3(0.954159614104838,1.48983739440767,-0.979575391698214) q[1];
u3(1.98758825363752,-1.30803634076232,-3.90179596399844) q[3];
cx q[3],q[1];
u1(0.562531974411469) q[1];
u3(-3.45661519237407,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.67455379382662,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.40375004264046,-0.854880239977576,1.96655046740298) q[1];
u3(1.18160461889775,5.39529465141605,-0.224400239588903) q[3];
u3(1.75011968455545,1.26874588326295,-0.833582594540370) q[0];
u3(1.15995948531335,-1.03833357283237,-3.24863945967479) q[3];
cx q[3],q[0];
u1(0.422169291484008) q[0];
u3(-1.31940668327403,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.56454941257079,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.759296858018670,1.38573502332916,-4.65920575734898) q[0];
u3(0.800410511466963,2.98539782969640,3.15159107610339) q[3];
u3(1.69193774786230,1.90948887052613,-3.72330728790977) q[1];
u3(2.78344303543793,2.96252959268439,-2.60275310600300) q[2];
cx q[2],q[1];
u1(3.37570582305065) q[1];
u3(-1.81098036848897,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.997947597187389,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.65894326349102,-1.32363478370947,2.85074994412795) q[1];
u3(0.560273756774173,-1.93373843188197,-0.894194276663102) q[2];
u3(0.341662320037372,2.06069317938931,-0.498230821327426) q[2];
u3(1.76448513669503,0.139335435616247,-3.83189719499320) q[0];
cx q[0],q[2];
u1(1.94727622707417) q[2];
u3(-3.01889978450051,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.03617830566670,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.186556359123543,-1.65967914280651,2.61700165596972) q[2];
u3(1.98582170419912,-0.351716521123238,4.50256200738449) q[0];
u3(2.71160595666464,-2.49423513593274,0.536529234873280) q[3];
u3(2.02364146463080,-2.50888083199646,0.593976770484123) q[1];
cx q[1],q[3];
u1(1.28680750290809) q[3];
u3(-0.812602739446826,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.19434046748939,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.01694543327892,-1.09802970253721,-0.136681713241285) q[3];
u3(1.39608552776099,2.84197435525522,2.89737493486472) q[1];
u3(1.69379715026603,-0.0251403182151849,0.964669591379639) q[1];
u3(1.79552461907686,-0.176722114941263,-1.24069872752614) q[0];
cx q[0],q[1];
u1(1.66582835133505) q[1];
u3(-3.12515702145666,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.58772364434639,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.584352363461091,1.22070468801912,-3.17731760338483) q[1];
u3(2.88439141046316,-2.34726209864721,1.58456488710995) q[0];
u3(2.17720894015833,1.73529811839780,0.205040971001833) q[3];
u3(1.90476340608992,-0.130391682409815,-1.76468300901604) q[2];
cx q[2],q[3];
u1(1.29583570600450) q[3];
u3(-3.95421329925566,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.51310734557515,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.888492703839095,1.50024971161041,-1.62877779094904) q[3];
u3(0.531225865389534,-1.73991270049365,2.10428141243648) q[2];
u3(1.55091274853006,0.480503736694342,0.600527545977584) q[2];
u3(1.03844374106152,-2.52558666941060,-1.48287554471359) q[3];
cx q[3],q[2];
u1(1.59452510296639) q[2];
u3(-3.02253780068053,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.540522340199930,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.72870749069102,3.07367802384371,-0.884799570648063) q[2];
u3(1.75504638998247,-2.98383186200282,1.57184002931502) q[3];
u3(1.94276399039923,1.12333169063978,1.34846238884500) q[1];
u3(1.41992405090078,-1.60739219492806,-1.91399137552566) q[0];
cx q[0],q[1];
u1(1.64850535196466) q[1];
u3(-0.307069613573358,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.00998291642421911,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.74747476201077,-0.684905198106371,0.367827837647501) q[1];
u3(1.18219702464128,-3.34649666309140,-1.11214300150329) q[0];
u3(1.57354731759599,2.03433659694240,-3.18989360739327) q[1];
u3(2.06283641514311,1.94860936524705,-3.54401039855261) q[3];
cx q[3],q[1];
u1(3.39475191363407) q[1];
u3(-0.980988460724172,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.35989182813091,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.25982869819374,3.12100303276869,-1.82939072748913) q[1];
u3(2.18476882441195,-3.03675769314494,-1.17095110819122) q[3];
u3(1.55216543179898,2.62131922626766,-0.921432578898004) q[2];
u3(2.52292600719342,2.34619069623666,-1.33138444451851) q[0];
cx q[0],q[2];
u1(2.65687737277675) q[2];
u3(-1.93261010786320,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.229028574868573,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.33365081664327,-0.00247861510305647,-1.50856490910856) q[2];
u3(1.08680823155865,2.51960701737532,0.0122555760421288) q[0];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
