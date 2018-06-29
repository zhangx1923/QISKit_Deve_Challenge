OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.25440353443541,-0.0694330961007453,2.00075567960249) q[5];
u3(0.865531685044529,-1.31125108149613,-0.944672248280583) q[0];
cx q[0],q[5];
u1(1.64257315100992) q[5];
u3(-2.11886450213747,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.32939435170550,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.01032952853426,3.09545785117469,-0.822034959167628) q[5];
u3(2.13604133202315,-4.20025534177229,-1.37041491334181) q[0];
u3(1.75697408652080,-0.237167307839303,2.09346934730465) q[1];
u3(1.72170263635624,-2.24770121156573,-1.96013141933421) q[4];
cx q[4],q[1];
u1(1.92389826446911) q[1];
u3(-3.21044749221914,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.846245055894461,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.43236473365870,-2.34179988237947,2.38292520760908) q[1];
u3(1.87576911510220,-0.307895129902210,-1.76887754526026) q[4];
u3(0.986623990523901,-0.178745490929110,2.32959245972719) q[3];
u3(0.754673986176179,-2.12834633450021,-1.83006821509023) q[2];
cx q[2],q[3];
u1(1.80823162026872) q[3];
u3(-3.18150348500718,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.952835428983231,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.302106475896344,-3.54851798510716,0.293679084712590) q[3];
u3(1.57715871396407,3.06174106702078,-2.07424515137373) q[2];
u3(1.62167892395674,1.71369194484629,-2.90441962859103) q[0];
u3(1.94448006951969,-1.90720116417070,2.82108252644397) q[2];
cx q[2],q[0];
u1(1.64723384484007) q[0];
u3(0.0912170528909821,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.950403521038925,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.17354291439228,0.211121720091606,-0.918150474965824) q[0];
u3(2.70190236129924,-0.408151400353125,4.53241354476639) q[2];
u3(1.93775735263066,1.23851079714473,-2.67034290409000) q[5];
u3(2.22559407136665,-3.17627989671470,2.37072553366660) q[4];
cx q[4],q[5];
u1(0.699525269762919) q[5];
u3(-1.81223039193548,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.282019527386339,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.31123269201993,2.92204807793066,-2.21026728934508) q[5];
u3(2.38459839270491,2.48659629069923,-2.85602943918549) q[4];
u3(2.15342665231507,-2.59114938797797,0.0992003500982392) q[3];
u3(2.55670978367620,-4.13887922802351,-2.01146375649236) q[1];
cx q[1],q[3];
u1(2.30170199622682) q[3];
u3(0.198151005208434,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.33208815281916,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.69000813354171,-1.06830633166781,-1.70368792922488) q[3];
u3(1.35945744525741,-5.43135586513111,0.454807438404419) q[1];
u3(1.66860943065218,-0.534109376660884,1.41204316944286) q[4];
u3(1.76286894517175,-1.46457748098948,-0.872980206297133) q[1];
cx q[1],q[4];
u1(0.870713878405933) q[4];
u3(-3.39922366514956,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.00979424299477,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.59751124539404,-4.32843379446861,1.86630830252813) q[4];
u3(1.47605765589622,2.37119000781501,-3.42252218781319) q[1];
u3(2.21730544899502,-3.89549493268117,1.77888461636268) q[0];
u3(0.128141507657568,3.35389167984694,-1.13063776955726) q[3];
cx q[3],q[0];
u1(2.83830429295731) q[0];
u3(-1.98125351335709,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.547083287134124,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.864931722774473,-0.0584775200817773,-3.94937395617937) q[0];
u3(1.61631252663824,3.84136626747236,2.19036719811171) q[3];
u3(1.40206172859568,0.530332816776884,1.70212788917337) q[2];
u3(2.17703014394607,-2.38288564273912,-1.26630125166748) q[5];
cx q[5],q[2];
u1(2.45774331876159) q[2];
u3(-1.67860808919859,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.0713229316377182,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.79738378954456,-2.37519917386939,-0.288550044983845) q[2];
u3(1.68054391828203,0.648299756042816,-2.57661139725520) q[5];
u3(0.866658236907101,1.93653020872104,-2.93950410286163) q[4];
u3(0.802848937061365,-3.10292535425289,3.01250828935496) q[0];
cx q[0],q[4];
u1(2.65962793557049) q[4];
u3(-1.82297610394002,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.532629325380638,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.22662374910696,-0.156391449929406,-0.0395752467330530) q[4];
u3(0.498811834394991,-2.43355453127218,0.639904224633172) q[0];
u3(1.93833038071914,0.111606157072227,2.64039229342711) q[1];
u3(1.10144050085474,-0.984455484609684,-1.67649252723079) q[5];
cx q[5],q[1];
u1(-1.20544094359877) q[1];
u3(0.336170556248724,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.70547194381488,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.59394911720216,-0.167170440719130,0.00469697030949176) q[1];
u3(2.99157999188560,-1.05999353284444,-2.95718120181636) q[5];
u3(0.990132234569116,1.89869143939419,-3.21290615467365) q[2];
u3(2.30625839016612,-2.82380786835840,3.29085426091545) q[3];
cx q[3],q[2];
u1(1.70791504495694) q[2];
u3(-3.06344273909635,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.36385408262278,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.19559271698711,-1.89659012372656,2.72353057266218) q[2];
u3(2.18219231351167,-5.84089835409720,-0.0461011663420554) q[3];
u3(2.30434872773941,-1.73927653752270,1.67788947083096) q[5];
u3(2.30971999092268,-3.80257721653722,-2.16616210802213) q[1];
cx q[1],q[5];
u1(1.37925726681485) q[5];
u3(-1.07786400745447,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.04976824777054,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.93365796824811,2.75319508451006,1.30845167212037) q[5];
u3(1.88096141095566,-0.794775656035177,2.20801107904845) q[1];
u3(1.56990113261705,-3.65793242282082,1.34517125313617) q[3];
u3(1.99607647112720,0.325819599998886,2.79338294781763) q[2];
cx q[2],q[3];
u1(0.442387983220812) q[3];
u3(-1.21219767298970,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.08298701070566,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.76621698515043,-2.72026681634819,1.87144065453825) q[3];
u3(2.39579782249362,-1.86664881848483,-3.08240252922920) q[2];
u3(2.04322766515538,1.45219808844287,1.13250228197236) q[0];
u3(2.05831840425914,-0.156400955176910,-3.40671187743853) q[4];
cx q[4],q[0];
u1(2.57808227708722) q[0];
u3(-2.39570679656647,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.53973644198185,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.90121784275978,1.45617065736339,0.00903576545349805) q[0];
u3(1.44520012687281,1.32212772995085,-1.68587264138060) q[4];
u3(1.94047656047775,-1.13526330183924,1.62349768552268) q[2];
u3(1.58382878209443,-1.35557278434363,-0.951050279809684) q[1];
cx q[1],q[2];
u1(3.18759875517653) q[2];
u3(-1.46309867177095,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.490670414438147,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.91058277343057,2.89406227443916,-2.63715920405894) q[2];
u3(0.534687315408658,-2.84938143634324,-0.0904577611076749) q[1];
u3(1.85745682135310,0.972265650637430,-3.33191124424004) q[5];
u3(2.45097625558020,4.00128282610303,-1.80595387614140) q[3];
cx q[3],q[5];
u1(1.02363142387026) q[5];
u3(-0.266436518548678,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.60200889561626,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.53708948605565,0.313117412586952,-1.51547259664578) q[5];
u3(0.790167932414416,-0.481587164607783,-5.27550930846069) q[3];
u3(1.65091997192609,0.637579917299637,0.169708876246662) q[4];
u3(1.58434133880864,-1.02485346189266,-1.55682296709460) q[0];
cx q[0],q[4];
u1(0.236249033540457) q[4];
u3(-0.954654800919380,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.65307081654898,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.760786536983315,2.86514123570658,-1.28904462277544) q[4];
u3(2.34201214261270,3.33967542983422,-1.38831888389054) q[0];
u3(1.68300854593646,-0.808895604700824,-0.955080713780866) q[1];
u3(1.94545455888152,-4.85617101693386,1.17178600343497) q[0];
cx q[0],q[1];
u1(4.16548173231734) q[1];
u3(-3.64777530109525,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.824562929097658,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.56872456665754,-3.65183790074655,1.76749295279889) q[1];
u3(1.47367205780714,-2.43391897483049,-3.18127151246503) q[0];
u3(2.04221422357020,2.25125093363953,-2.98530617430332) q[4];
u3(2.33565653648997,2.51326032404095,-3.73247728018859) q[2];
cx q[2],q[4];
u1(0.893600808234103) q[4];
u3(-1.67539888452232,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.83513007868195,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.29387890867920,-1.35958531856856,3.32014816473464) q[4];
u3(1.33612685002503,-5.57425412802434,-0.182801507882643) q[2];
u3(0.948233667440802,1.48990892683470,-2.93209277894211) q[3];
u3(0.399310292686776,2.01224296954474,-3.41314469000038) q[5];
cx q[5],q[3];
u1(3.10217557306150) q[3];
u3(-1.59302089642546,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.567384568071386,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.47762652041339,-2.53855973116857,2.87763200478060) q[3];
u3(0.696411418830334,0.124672575819894,3.46715916337204) q[5];
u3(1.78067074327297,1.65496146774436,-4.49297919437704) q[2];
u3(0.226612574683439,-1.10146313259038,2.93126848794523) q[1];
cx q[1],q[2];
u1(2.28061110444911) q[2];
u3(0.335847853306451,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.08535070326532,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.701092510914833,2.39265737964020,-0.472832220393118) q[2];
u3(1.39847741667895,0.623373456547290,-2.01400953247742) q[1];
u3(1.86486238611522,2.13772778994386,-2.25606480095682) q[4];
u3(1.41329530312910,-2.95993416152949,2.52682992191411) q[3];
cx q[3],q[4];
u1(0.863354959581427) q[4];
u3(-0.0161378677769610,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.86760703093914,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.61044740324309,-2.40985380857379,0.0116431193053126) q[4];
u3(2.06484104743241,3.00939139163886,-2.63213465392993) q[3];
u3(1.81470869394180,-0.625614308594404,-1.63684221497380) q[0];
u3(2.27768103747372,-0.160608760080918,-5.86308752860455) q[5];
cx q[5],q[0];
u1(2.02327870433994) q[0];
u3(0.00965617095400617,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.60224545116835,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.17359991464550,2.85369015146765,-3.15809153226273) q[0];
u3(1.96541141779854,2.64771626684092,2.16234973169511) q[5];
u3(1.10143274362424,0.668491504559044,0.978566744754872) q[3];
u3(2.15508790999670,-0.831594812937670,-3.08342043187185) q[0];
cx q[0],q[3];
u1(4.51533407396352) q[3];
u3(-3.93642712661796,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.667035550943274,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.36871413488439,-2.36634604941559,-0.371547103753306) q[3];
u3(0.265798972677070,0.438209335709562,0.589952727051961) q[0];
u3(1.63793666557447,-1.63434951019062,-0.251709748559291) q[4];
u3(1.23243175992198,-4.58332092364426,0.181904501961977) q[2];
cx q[2],q[4];
u1(1.87080352863792) q[4];
u3(-0.114060670710763,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.36021059945308,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.37143246405530,1.32611996349694,-0.501145050104398) q[4];
u3(1.38856264077716,0.640735437478718,4.30226384623116) q[2];
u3(1.92402185640711,4.18572828316347,-1.56111002342131) q[5];
u3(0.922070826992779,1.59752967282652,-0.0110444122630982) q[1];
cx q[1],q[5];
u1(-0.179557129191436) q[5];
u3(0.0821735183901511,0.0,0.0) q[1];
cx q[5],q[1];
u3(4.12083123705342,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.94044343281524,3.39713785391014,-1.57487485564290) q[5];
u3(1.53803518222981,-2.47101712578926,2.65124933530263) q[1];
u3(0.872912900666522,0.975177285037772,-0.0219692887572518) q[1];
u3(1.08000197776219,-1.07940980137453,-0.782982163325023) q[0];
cx q[0],q[1];
u1(1.70672162150433) q[1];
u3(-2.49283414307016,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.91530554980012,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.17556875547671,0.0144706398771790,-3.14570319165700) q[1];
u3(0.835852725274048,3.47265181745014,-1.93174544720147) q[0];
u3(1.36019936375957,0.773155445673348,1.17190523860069) q[3];
u3(1.03813038835470,-1.22314346638797,-2.73086375044630) q[5];
cx q[5],q[3];
u1(1.96463275108240) q[3];
u3(0.210827818605121,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.656539601447329,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.768363195499546,2.05332217007172,-3.02827719486931) q[3];
u3(2.00687682322038,1.49026862504354,-3.08772918071689) q[5];
u3(2.49816682001531,-1.09880611354162,3.86170988506209) q[2];
u3(1.57828886193791,1.76926032602582,2.13328923112737) q[4];
cx q[4],q[2];
u1(0.228527890867696) q[2];
u3(-1.33436927702487,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.37618326702668,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.09960618658005,-2.02604595218349,0.337272946852202) q[2];
u3(1.38872178991049,1.11987875927344,-3.92968695530132) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
