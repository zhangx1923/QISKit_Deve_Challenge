OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(2.00474295481797,-0.807606145880985,-1.34891043277612) q[1];
u3(0.259785848568386,-4.18406406858652,0.504720449352299) q[4];
cx q[4],q[1];
u1(0.902918297052976) q[1];
u3(-1.32091013384052,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.54309396790708,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.81415458030390,-2.06156609326268,-1.36197538082169) q[1];
u3(1.92100589138699,-3.39246895935694,1.15892081023275) q[4];
u3(1.33601780866818,0.00957425188703931,1.30639941556461) q[5];
u3(1.49223034174271,-0.237553062500607,-1.41517232681483) q[3];
cx q[3],q[5];
u1(-0.273601444751716) q[5];
u3(-2.24282417235177,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.97015795577738,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.48548819414571,-3.48453289346502,0.366813920077612) q[5];
u3(2.01851304438495,2.43701320097546,-3.43126271915964) q[3];
u3(1.77122015589984,-0.632603071877087,-1.12538579078169) q[2];
u3(1.86252139064043,1.21532499841229,-4.37011133666828) q[0];
cx q[0],q[2];
u1(2.79343755861158) q[2];
u3(-1.39010383826083,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.191079991105453,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.54730880247065,-1.35701195197172,-0.638624183894675) q[2];
u3(2.19186971665718,4.54027442907244,1.62807340930530) q[0];
u3(2.90846286829247,-1.59506190687755,0.275406697200804) q[5];
u3(2.79023673733969,-1.68127547567133,0.150032826033663) q[3];
cx q[3],q[5];
u1(0.942693555542478) q[5];
u3(-1.14472183802305,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.12771869056082,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.04242025498660,2.87217437847217,-2.98074777450298) q[5];
u3(1.19333548153217,-3.23294143045747,-2.83194258650543) q[3];
u3(2.52277354464978,-0.325153839928703,-0.279890681320027) q[1];
u3(0.542549970325057,-0.0170554031495218,-5.26180666784677) q[0];
cx q[0],q[1];
u1(2.47746463264663) q[1];
u3(-1.92718009269705,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.0574348110309180,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.64384405852006,4.77947909512285,-1.33774976312876) q[1];
u3(1.59162323778914,-1.93215173451032,-0.285234295110719) q[0];
u3(0.156226157505988,0.747990521092523,-0.714518225979578) q[2];
u3(0.893691635842317,1.22396071570114,-1.80912338040951) q[4];
cx q[4],q[2];
u1(1.57472366137475) q[2];
u3(-2.31206421203939,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.24700454327643,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.76669767155861,0.354091370753557,0.000569707669070163) q[2];
u3(1.75307913606834,-0.0610604060542141,-4.11904021647192) q[4];
u3(2.74168887561390,0.286643158569761,-0.339789321983480) q[3];
u3(1.15087794940873,-2.59521466996618,-2.12903161608230) q[1];
cx q[1],q[3];
u1(1.79250881915764) q[3];
u3(-2.97705123625108,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.500669904602143,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.32778976867773,1.02320362572306,-0.722429670247224) q[3];
u3(2.24687842399670,-0.689834529003795,0.416147770522450) q[1];
u3(2.07232554842602,1.31784257658412,0.828112445543563) q[4];
u3(2.35495570526281,-0.105511404886162,-3.67282887560480) q[0];
cx q[0],q[4];
u1(2.87882363862686) q[4];
u3(-2.44253586589582,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.88693490269095,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.98485341817767,1.07195511078784,2.39868083357792) q[4];
u3(0.613212877383860,-2.88752058246099,3.15140433722243) q[0];
u3(0.455172973687600,-0.544487224318616,-0.319989717512114) q[2];
u3(0.713237873371929,-2.97239284781699,1.69936352491266) q[5];
cx q[5],q[2];
u1(0.899190637779410) q[2];
u3(-0.192853316657014,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.88065491300367,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.66053683998672,-0.622935305156407,3.11436918398111) q[2];
u3(0.501638080464967,-4.56984881874072,0.0803449914628476) q[5];
u3(2.02060208418293,3.32714179015099,-2.52148019622744) q[3];
u3(2.34149128094663,1.06450332291811,-1.77156043758114) q[1];
cx q[1],q[3];
u1(2.31813244792217) q[3];
u3(-2.71146365893771,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.998396510183103,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.585747874191760,-2.65072413591691,2.98723981857006) q[3];
u3(2.03670370678414,-0.0206054358869929,3.41473738542821) q[1];
u3(2.21819935738440,1.35759185315941,-0.833964557558989) q[4];
u3(1.96231666993820,-0.245765749392621,-2.99124012646315) q[0];
cx q[0],q[4];
u1(2.55352598416140) q[4];
u3(-2.07462512324778,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.832496566589910,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.16112747712458,2.29831623179215,-3.44442218777213) q[4];
u3(1.79728326337655,-2.40475160606095,1.91012860549296) q[0];
u3(1.97809547469156,-0.0122170837966050,-0.183057704207302) q[5];
u3(1.13411544536232,-0.0253329649574039,-4.95629166229645) q[2];
cx q[2],q[5];
u1(2.95693309278519) q[5];
u3(-1.88729699742181,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.555370993252575,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.14690812974791,1.52017047721428,0.782298262372390) q[5];
u3(1.89910149281452,2.54197413627692,-2.05042683084339) q[2];
u3(1.60942769539380,0.879643085806067,2.19454728266819) q[2];
u3(1.63088080816518,-0.913548690066021,-1.00837227695165) q[3];
cx q[3],q[2];
u1(2.34575124513682) q[2];
u3(-1.65200461437050,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.36463824469608,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.86775945990209,0.925002713928572,-0.118646818152986) q[2];
u3(2.21468673628902,2.95206703758639,-2.86416295800125) q[3];
u3(0.796161735572611,3.34866424267071,-1.22987191345311) q[4];
u3(1.56532894030577,1.74216255921824,-1.15392483888671) q[0];
cx q[0],q[4];
u1(0.815989650567180) q[4];
u3(-1.43501358651684,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.278283972198824,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.68934114498229,-2.65670676537061,2.60827309131735) q[4];
u3(2.15025556857564,-3.98799089220167,1.91132408861873) q[0];
u3(0.168251118073259,1.47250927481126,-2.69986005362286) q[5];
u3(1.22456690989668,-3.51903763475094,2.20766285728552) q[1];
cx q[1],q[5];
u1(0.512614952820221) q[5];
u3(-1.23941783509944,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.67895152113315,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.630891756337663,0.338948099477047,-0.510094217683087) q[5];
u3(0.979011603106156,4.27742209328785,1.47665139945317) q[1];
u3(1.82052979566297,-0.608707865176568,-0.0343726222844567) q[2];
u3(1.44965821860179,-3.51426481348974,0.671628533495167) q[0];
cx q[0],q[2];
u1(0.646793790751124) q[2];
u3(-1.58020788969150,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.20762884374711,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.20308016409154,1.75269278217476,-3.08514310344059) q[2];
u3(2.04928666674492,2.34225671349822,-3.45895848794814) q[0];
u3(2.69549923574247,-2.85532205379232,2.43530290023267) q[3];
u3(1.49026464472873,3.07990285623701,-1.42305354290336) q[1];
cx q[1],q[3];
u1(-0.462640341783210) q[3];
u3(0.307631979313338,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.93688219579695,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.31332495165857,1.48542265610976,0.713921764244404) q[3];
u3(2.31551483545886,2.74760513015954,-2.52210654588048) q[1];
u3(1.70488070389212,1.00025572319683,-3.33678737197155) q[5];
u3(2.93928282604528,4.29348637098086,-0.964764588080827) q[4];
cx q[4],q[5];
u1(3.07045436743727) q[5];
u3(-1.64450442652149,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.97771663722981,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.86803538138616,0.877537203267402,-0.460149850674465) q[5];
u3(2.14092062365605,3.15961868836652,-2.18716391066265) q[4];
u3(0.888814980761867,2.19980144585296,0.777589183343597) q[0];
u3(0.947225824465906,1.17144114675825,-4.23874433980524) q[3];
cx q[3],q[0];
u1(2.52285838632967) q[0];
u3(-1.74127359151724,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.163281004532438,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.12384474449792,3.40831134362719,-1.65422694469748) q[0];
u3(2.59040638691882,-4.24455574634902,-1.73717859832093) q[3];
u3(1.21362321478183,2.28105369271227,-0.836533016847553) q[5];
u3(1.64613864464165,-0.783925152830994,-3.14088397539616) q[2];
cx q[2],q[5];
u1(3.26304288576479) q[5];
u3(-3.86619852700555,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.550863368192850,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.41145562175356,-2.82710691571444,3.26793346541480) q[5];
u3(1.45247366078456,2.08260560638538,2.42499905364777) q[2];
u3(1.35998486369072,3.65878348964299,-1.85664364038012) q[4];
u3(1.22429361460485,2.32971521140906,-2.24770296519691) q[1];
cx q[1],q[4];
u1(0.812221314435191) q[4];
u3(-1.23422560538268,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.471852894587284,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.74775982520898,2.01214159526151,-3.07105090264740) q[4];
u3(1.01760109550823,-4.32953970227404,-0.803005276705454) q[1];
u3(0.445235817849258,-2.14329561071244,0.649925878507871) q[1];
u3(1.51745754370400,-4.40930884739519,-0.0806622363752743) q[0];
cx q[0],q[1];
u1(2.20811282436412) q[1];
u3(0.208348467064193,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.15732362340906,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.48609811329079,-1.59918641370461,3.22318957761493) q[1];
u3(1.01553610197909,3.36104602249019,1.37854656526853) q[0];
u3(1.30893425935971,1.92023179769697,-3.13407519649275) q[3];
u3(0.766354316450912,-2.28746495972945,2.67607972391847) q[5];
cx q[5],q[3];
u1(-1.33100858469199) q[3];
u3(0.362062687382488,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.59922893345670,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.520597180989089,-1.37441958543883,0.638982917786798) q[3];
u3(1.26674268785703,-3.57482545803193,-2.27918303641561) q[5];
u3(2.77928073012854,1.46187342094015,-3.75573333624481) q[2];
u3(1.77891616697278,-2.71505507152372,2.99631477163409) q[4];
cx q[4],q[2];
u1(2.36644004470100) q[2];
u3(-1.68169328753220,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.08815673655742,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.489071870777013,-2.24323634889489,1.64665375587495) q[2];
u3(0.365111377752813,0.811790094153513,-3.53959095176822) q[4];
u3(1.73548241180055,-0.782879454733252,1.63424891945877) q[1];
u3(1.12281904385534,-1.20875475071141,-0.619977921360911) q[0];
cx q[0],q[1];
u1(2.32527906308426) q[1];
u3(0.312492047652739,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.67823708199541,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.41929386474580,-2.13384640442806,-0.336194944185016) q[1];
u3(2.58014367221981,1.84413565369352,3.76597428464493) q[0];
u3(2.54444731264029,2.68008999343705,0.270036163181861) q[5];
u3(2.10266748428528,0.236463983786493,-3.41157182094138) q[2];
cx q[2],q[5];
u1(3.43544582284140) q[5];
u3(-1.35883504233503,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.18098493663349,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.61549808482828,-2.20910066786961,0.0133067582273241) q[5];
u3(0.766178799753715,4.63310551277175,0.256574754109477) q[2];
u3(1.89662555832266,-1.28124897748901,-1.14983755134660) q[4];
u3(0.828241827803181,-1.97426130974539,-3.14298024118321) q[3];
cx q[3],q[4];
u1(0.795958980933220) q[4];
u3(-0.287421276684721,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.58764168261898,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.921406979445182,-1.67825038523643,1.32065679054859) q[4];
u3(0.986981134607949,-0.324189918858070,-5.88083863389736) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];