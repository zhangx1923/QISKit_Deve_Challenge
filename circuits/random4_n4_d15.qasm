OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(2.40972100704112,1.13688777965431,-2.79871270847769) q[2];
u3(1.01700467070136,-2.08140455385856,2.70489525884371) q[0];
cx q[0],q[2];
u1(-0.386115155132372) q[2];
u3(-1.98772091763519,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.826347559619727,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.652546167902342,-0.751958649322446,-0.384435861354734) q[2];
u3(2.51293940817066,1.28680179302702,-4.65153745099345) q[0];
u3(2.12097932959054,1.10599313741290,-2.41401550711233) q[1];
u3(2.78450648081862,5.48013326984016,0.622393137605716) q[3];
cx q[3],q[1];
u1(1.75006981847258) q[1];
u3(0.299955141954279,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.06669114501215,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.91137854429628,-2.45981114375604,2.19060755639354) q[1];
u3(2.04281234602594,-5.13635368286559,-1.01784287532649) q[3];
u3(1.68773763460979,-2.21545768464537,-0.578432204804010) q[1];
u3(1.21435947630105,-2.92964354013158,-0.768901365829140) q[3];
cx q[3],q[1];
u1(3.44972886545986) q[1];
u3(-1.13119514987181,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.90149626971320,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.943595341592619,-0.468970405683564,-3.54717870353416) q[1];
u3(1.94673268994621,2.16628525330504,0.943005552727967) q[3];
u3(1.10286980021463,2.31717246899895,-1.85849090349306) q[0];
u3(1.23825252581049,1.10023010365920,-2.70676904152147) q[2];
cx q[2],q[0];
u1(-0.531578040263669) q[0];
u3(0.979176081182905,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.88568427595871,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.07300136237639,0.691954529838716,0.484214066929960) q[0];
u3(2.20028408325089,-1.79817885741269,4.31611162806845) q[2];
u3(1.53686198000317,-0.789004466136812,-0.681354352122236) q[1];
u3(1.12687493429159,-2.85940409931043,-0.0374482828064580) q[2];
cx q[2],q[1];
u1(0.582985769923920) q[1];
u3(-1.27726414104048,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.15222075722584,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.62430456499932,0.922352638671956,1.97440308846832) q[1];
u3(2.52910008004463,-1.51580473132649,-0.989607088065995) q[2];
u3(1.87522739820846,1.02811595918417,1.70907116464293) q[0];
u3(1.20590539014430,-0.819918372542457,-2.73592835345838) q[3];
cx q[3],q[0];
u1(1.89759695830929) q[0];
u3(-2.58838577181611,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.0667431226326305,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.13265797278413,3.30870138608873,-0.138329805888602) q[0];
u3(2.27920975160031,-1.70457378431171,-3.21012943868774) q[3];
u3(2.66221947330563,1.22607843638399,0.461437861196050) q[0];
u3(1.61917115688366,-0.810491568705842,-3.30352645549641) q[3];
cx q[3],q[0];
u1(0.857614139520181) q[0];
u3(-0.379337734997078,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.42694311102315,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.54397432112720,-1.45793068721060,0.446395217868219) q[0];
u3(0.985565280475458,0.373752751006921,-2.57351990148289) q[3];
u3(1.77205528899787,2.28935518607397,-3.89895888140242) q[1];
u3(0.760717426024029,2.05048528904194,0.214453179771414) q[2];
cx q[2],q[1];
u1(1.31994947994285) q[1];
u3(-0.414031109870336,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.83737631005319,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.15950590692374,0.143642581617194,-4.18588364096499) q[1];
u3(1.00716956397254,-3.54306901719642,0.547967120643999) q[2];
u3(1.65764158514272,3.32209585345938,-1.37239515827162) q[0];
u3(1.11424829595797,1.41368203172901,-0.907587542614172) q[2];
cx q[2],q[0];
u1(0.783176289302468) q[0];
u3(-1.30477195831711,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.210281801078132,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.861745662144724,-4.05490423840092,0.484326399103849) q[0];
u3(1.82179589583778,0.582771346202906,4.62159077187254) q[2];
u3(2.56272693195376,1.41966389366264,-3.76363374534318) q[1];
u3(1.59385400142203,-1.93202581717827,3.38745912521175) q[3];
cx q[3],q[1];
u1(3.25844288240251) q[1];
u3(-4.00255585367573,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.526785559814464,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.08036464470449,-1.51796583540161,-0.156713199113881) q[1];
u3(2.14797312554238,-4.84163596059331,0.446526665223649) q[3];
u3(2.08512656317141,-0.840831350672039,0.931270575373609) q[0];
u3(2.41130401619606,-2.88484960788041,-0.0748135764582272) q[3];
cx q[3],q[0];
u1(1.14271052462876) q[0];
u3(-0.596380421510666,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.14178446147506,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.385628528275546,0.434720991485471,-2.57905561080054) q[0];
u3(0.328454700660290,-1.86145044173356,-2.67114076935663) q[3];
u3(2.66405517346427,1.91912334728011,-2.53104616974082) q[2];
u3(1.70534993068911,-3.18658719909206,2.49109299816081) q[1];
cx q[1],q[2];
u1(-0.411818805378770) q[2];
u3(1.26577858941256,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.38441190385481,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.17939628365733,1.98271959918611,1.49868766518226) q[2];
u3(1.49190872764393,-3.76045023632176,-2.04520257211300) q[1];
u3(1.06457300350058,3.44745900161363,-1.84837653149130) q[0];
u3(1.55608570335982,1.25399748984650,-1.90396856554587) q[3];
cx q[3],q[0];
u1(0.836175920446453) q[0];
u3(-1.18107749128780,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.77287072043137,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.32377128723269,0.636604376861506,2.46772901096845) q[0];
u3(0.607674649644486,-0.0368391772770180,0.446180201676358) q[3];
u3(1.62895810739406,2.01363160738993,0.376911281314509) q[1];
u3(1.89659605980218,0.269724703588932,-1.90508968687220) q[2];
cx q[2],q[1];
u1(4.23495999057916) q[1];
u3(-3.52792147283598,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.805286830881102,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.613889741320078,-1.64093361383120,2.92511096184262) q[1];
u3(2.28244901824282,3.62975782581569,1.02512205575390) q[2];
u3(1.19547813513122,2.00760835216930,-3.62622664656401) q[3];
u3(2.06692649798674,3.59225896703938,-2.37640364730957) q[0];
cx q[0],q[3];
u1(1.60642987967913) q[3];
u3(0.0805446988965530,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.804813503706988,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.39270825387204,-1.19260013531315,4.06111519656688) q[3];
u3(0.747518149521754,2.17489536634996,4.09898927721324) q[0];
u3(0.717604135907448,2.56353660925416,-0.457436186821132) q[1];
u3(1.52142319781121,0.515751867241433,-3.76322248291275) q[2];
cx q[2],q[1];
u1(2.96383731330220) q[1];
u3(-2.27074865408992,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.43471200983248,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.09569071919599,-3.44925076012603,1.38531346388189) q[1];
u3(0.693325912918254,-4.75231581628801,-0.431348134220426) q[2];
u3(2.51882472267373,-1.35786299126684,2.15713031313465) q[0];
u3(1.95148509445249,-2.29391527515548,-0.792515256862680) q[1];
cx q[1],q[0];
u1(-0.803250876672359) q[0];
u3(0.674258373900647,0.0,0.0) q[1];
cx q[0],q[1];
u3(4.45165703129059,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.23071927662455,2.20793569193375,-0.925156281110749) q[0];
u3(1.31975007790940,-1.96012814624624,-1.45231764769522) q[1];
u3(1.46365214456659,-0.691473819129376,-1.04400132374970) q[2];
u3(0.651210756709750,1.48946239972020,-4.77541093112421) q[3];
cx q[3],q[2];
u1(3.55444980291232) q[2];
u3(-1.28027339281293,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.10056962056205,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.26797256230501,-3.26176532191003,2.04063716633121) q[2];
u3(1.96450445578558,3.19630912994970,-0.985430237525680) q[3];
u3(2.06813029726975,2.34154518724759,-3.56101753488165) q[3];
u3(2.42240555755807,-3.54883651794402,2.71879997897042) q[0];
cx q[0],q[3];
u1(3.31408947498147) q[3];
u3(-4.57659780333715,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.244789341991021,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.08110603109519,0.336735391563021,2.26895604569212) q[3];
u3(1.28069395498560,-0.364955752212907,0.875439128820800) q[0];
u3(0.853243598522039,-0.183893324173000,-2.25508505443097) q[2];
u3(2.28763468363732,-2.78213016724657,2.71636899686736) q[1];
cx q[1],q[2];
u1(3.71118722800818) q[2];
u3(-1.42437648637997,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.87303880332541,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.64576129027275,0.410897766259935,-2.99389717513922) q[2];
u3(1.17137549337986,4.95898085659971,1.00922136141297) q[1];
u3(2.10987117566676,-1.79080143260444,-0.217978646039784) q[3];
u3(2.26661448794343,-2.65168556002400,-0.909748661155933) q[1];
cx q[1],q[3];
u1(2.96128146385277) q[3];
u3(-2.64153955053458,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.53037426922489,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.443755274077206,-1.38181044911100,1.77676107008039) q[3];
u3(1.15613439375846,-3.79266274620069,2.01469897483661) q[1];
u3(0.608443172520103,-2.29127365789627,2.65933527132210) q[0];
u3(0.683001871298021,-3.71660035886617,2.35539460800849) q[2];
cx q[2],q[0];
u1(1.78343250311422) q[0];
u3(-3.17562678412946,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.368500980818790,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.35381014947385,2.52434909349686,-0.389162057744955) q[0];
u3(2.55734980332470,-2.68534394431393,1.36586673278443) q[2];
u3(1.03534925436382,-1.32615438575543,-0.835305819644869) q[0];
u3(2.27978639619560,-3.16942397414753,-0.586982079669479) q[2];
cx q[2],q[0];
u1(0.394862116147507) q[0];
u3(-1.32227811524033,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.14445330510334,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.36939406293768,2.12543730487844,-2.07464653475876) q[0];
u3(1.50181395814226,0.268613316862746,-1.44472975558679) q[2];
u3(0.334370805242926,-1.95056308724816,2.63171940888905) q[1];
u3(1.24946580136431,0.821692311581153,-1.30308966051700) q[3];
cx q[3],q[1];
u1(1.58734654568620) q[1];
u3(-3.06216678340811,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.268997896723907,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.45721533782825,3.54152270060660,-0.601056715816779) q[1];
u3(0.643090464021546,2.37113965364008,2.94870258776637) q[3];
u3(1.79330722488372,-0.607382962392846,0.621979303362174) q[2];
u3(0.978864130827324,-2.66317611883785,-0.995777701262968) q[1];
cx q[1],q[2];
u1(2.76700533800048) q[2];
u3(-1.53386908391103,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.03729958985544,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.779664362261217,4.67091020661358,-1.28168881638979) q[2];
u3(2.42830874711919,0.216437191166503,3.77587430666354) q[1];
u3(1.14546635154202,-0.354196595093503,2.96889509309966) q[3];
u3(0.975219885101147,-1.30252071468179,-1.80529085031271) q[0];
cx q[0],q[3];
u1(4.53692276674680) q[3];
u3(-3.58072358485176,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.552038431981698,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.51782854361820,-1.25979857378849,2.98282250619273) q[3];
u3(1.09420266149412,4.40239450406303,-0.650602652668159) q[0];
u3(2.07490192266257,-2.09922050742950,3.45782206178565) q[0];
u3(1.16092480238923,1.12268031017220,1.28384662949453) q[1];
cx q[1],q[0];
u1(1.72851820002677) q[0];
u3(0.453936077857656,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.942567290929892,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.49414831114661,-0.914081004886986,-1.32668907318013) q[0];
u3(2.05481755530075,3.51959739729377,0.499407722239806) q[1];
u3(2.57240108657001,-1.38862548511271,-0.840075515218876) q[3];
u3(0.556124371855197,1.01336057209485,-4.93239991239750) q[2];
cx q[2],q[3];
u1(3.21459606168032) q[3];
u3(-3.84917423623793,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.904620643994952,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.371137994282333,1.69802584909806,-0.272652992607001) q[3];
u3(1.92346812027839,3.84149737085941,-0.789381687539906) q[2];
u3(1.63003543056145,-0.629083187422324,-1.55869766595192) q[1];
u3(0.757603576567670,-4.37519577958844,0.704591206630623) q[0];
cx q[0],q[1];
u1(1.96646941275254) q[1];
u3(-2.87897747072093,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.695540337468778,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.37382193377207,0.294424086386034,0.289082403465412) q[1];
u3(2.01412946679573,-2.78311810135542,1.16692675039678) q[0];
u3(1.50346517818848,-4.23486315064812,1.94976446378046) q[3];
u3(2.58920857452072,-2.49968370410112,3.52026938221511) q[2];
cx q[2],q[3];
u1(1.74465501643645) q[3];
u3(0.225225355933571,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.629172678012174,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.78017709774896,2.02476704626798,-4.01342322179698) q[3];
u3(1.62436676365918,3.31878819309136,-0.794574626127706) q[2];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
