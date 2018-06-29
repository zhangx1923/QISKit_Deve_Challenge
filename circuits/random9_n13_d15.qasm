OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(0.418733293972544,0.405539609510214,-0.725088031329875) q[5];
u3(0.780872474856782,-2.67394739228598,1.47632754335159) q[7];
cx q[7],q[5];
u1(1.88505918050644) q[5];
u3(0.112382031161289,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.36776504718212,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.11640132524287,-1.54177541105715,0.301658280002049) q[5];
u3(1.22178677619849,-0.319919449262389,3.94034052633488) q[7];
u3(1.26982078739100,1.67406298218604,-1.79292236233717) q[12];
u3(0.0981323393455333,-2.50383174026671,1.72791722271985) q[6];
cx q[6],q[12];
u1(0.740882821738720) q[12];
u3(-1.12819132436629,0.0,0.0) q[6];
cx q[12],q[6];
u3(-0.179933445759071,0.0,0.0) q[6];
cx q[6],q[12];
u3(2.57607958423496,1.06414360960000,1.96110210556754) q[12];
u3(1.95029650215627,3.21542504828314,-2.39915467460249) q[6];
u3(2.47131692238418,2.84621312424018,-2.52052730717832) q[10];
u3(0.282138342139835,1.44140712250251,0.139613213717110) q[0];
cx q[0],q[10];
u1(-0.679048392546307) q[10];
u3(-1.73810457784659,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.09644156060417,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.33778450312887,1.67891120816688,-0.395032716339882) q[10];
u3(0.488665045918526,2.71536783801480,-2.57706875057207) q[0];
u3(1.98055367198852,1.25104144292928,-4.03163954125491) q[3];
u3(2.11183514410817,-2.32102294807937,2.73784832314508) q[2];
cx q[2],q[3];
u1(0.778411903005780) q[3];
u3(-1.52322299841942,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.63825178351938,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.43596225045548,-0.672164299309637,-0.676401620055777) q[3];
u3(1.03854032117210,-1.72137473652952,-0.0968071352454605) q[2];
u3(1.77925644338755,-0.901887018709224,0.428559821386444) q[4];
u3(1.97262615156455,-3.17241569760420,0.241408353275815) q[9];
cx q[9],q[4];
u1(1.36140790219073) q[4];
u3(-3.42841469911136,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.36192039160776,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.61362590991624,-2.94366295868798,1.91648842933874) q[4];
u3(2.24343695773254,-3.71973528032960,-1.26361978326423) q[9];
u3(1.20334512157762,-0.322627625278741,1.71460931399630) q[8];
u3(1.63674927516188,-1.29936002059054,-2.02433287474338) q[1];
cx q[1],q[8];
u1(1.68246226256264) q[8];
u3(-2.61531885834807,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.729031054377785,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.17840365386544,-3.20154634604761,2.59580658369080) q[8];
u3(0.908972497974437,-3.63795835020247,1.89261378404869) q[1];
u3(1.67052401628244,-1.52493448599562,-1.35363541825653) q[12];
u3(0.573161617987116,-4.62776561510270,0.524921404702711) q[9];
cx q[9],q[12];
u1(1.55664704610150) q[12];
u3(-2.38007784591286,0.0,0.0) q[9];
cx q[12],q[9];
u3(-0.0118029818553758,0.0,0.0) q[9];
cx q[9],q[12];
u3(0.247100382147635,2.06436147663507,-0.378720705251194) q[12];
u3(2.01100190355373,-3.83105172014790,-1.82831298045774) q[9];
u3(0.273712429487242,2.09148962377544,-2.87802135943544) q[2];
u3(0.417479163582622,1.24429940449916,-2.45397569000255) q[0];
cx q[0],q[2];
u1(1.99265503477000) q[2];
u3(-3.11265767439135,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.594738342405779,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.00548811080051,-1.35800253257876,1.85982979150266) q[2];
u3(1.96492605762388,-0.240885859911736,1.55225380626496) q[0];
u3(1.55892967832055,-1.83834592291281,-0.551006898527121) q[7];
u3(1.74827376833749,-4.52013688373396,-0.720046211585272) q[5];
cx q[5],q[7];
u1(3.42329954900281) q[7];
u3(-1.15747984936970,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.28318322143345,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.46916249839644,-1.86558241252910,3.38169087862989) q[7];
u3(2.25821201083947,1.62720167144872,1.47833894579263) q[5];
u3(1.34501599486097,0.542111152279203,-0.913046665885768) q[1];
u3(0.940499914986585,0.868842014432878,-4.52475929576330) q[10];
cx q[10],q[1];
u1(0.804523306569332) q[1];
u3(-0.255477183178168,0.0,0.0) q[10];
cx q[1],q[10];
u3(2.73053539121396,0.0,0.0) q[10];
cx q[10],q[1];
u3(2.67172890398847,-2.58875173786634,-0.921157235476930) q[1];
u3(1.96093938099039,0.926626446423357,-2.54394375024238) q[10];
u3(0.706689036466443,-1.24788980199476,1.28298900102981) q[11];
u3(0.751490512279489,-1.57848661033885,-0.441851932620690) q[3];
cx q[3],q[11];
u1(0.435227255994987) q[11];
u3(-1.61981809138032,0.0,0.0) q[3];
cx q[11],q[3];
u3(2.93289899554286,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.20397331886973,-0.175251939537596,2.18792123567717) q[11];
u3(2.85837173834851,-2.24475380700522,-3.39389480062250) q[3];
u3(1.00450236192992,-0.0714706632630813,2.56857940242759) q[8];
u3(1.32471347095810,-2.41134553856719,-1.66441694289298) q[6];
cx q[6],q[8];
u1(0.345739278752449) q[8];
u3(-1.60655730231056,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.08084351489883,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.04473090366724,-2.23029680617626,-1.75861750824706) q[8];
u3(2.18964913160712,-2.52542313049078,3.32138449260173) q[6];
u3(2.48960489097390,2.29107381450014,-2.70578014815512) q[8];
u3(1.81120783783766,-2.74216644225331,2.86195373027217) q[1];
cx q[1],q[8];
u1(0.958733078330040) q[8];
u3(-1.45450883012629,0.0,0.0) q[1];
cx q[8],q[1];
u3(3.44808329509344,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.479856925753684,0.597068882018198,-1.16593501338825) q[8];
u3(1.84854002190521,-0.818235450327594,2.07035119656567) q[1];
u3(1.72572418481656,-1.21227456597471,4.13005916974357) q[5];
u3(1.69354438213135,1.67699623474883,1.92183476339546) q[11];
cx q[11],q[5];
u1(3.20572849291533) q[5];
u3(-1.00331470582948,0.0,0.0) q[11];
cx q[5],q[11];
u3(2.26699564821355,0.0,0.0) q[11];
cx q[11],q[5];
u3(0.373027117676741,-0.296763820936680,1.52966448185334) q[5];
u3(2.40402892108246,1.65462782043427,-0.237237665187724) q[11];
u3(0.889792487966941,-0.221174504701031,-1.75278660433024) q[7];
u3(2.22117371688393,0.0967795571119510,-4.77439983527387) q[4];
cx q[4],q[7];
u1(0.925604675134035) q[7];
u3(-3.31520032765220,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.71135444887205,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.32902204400965,-3.12273981722052,-0.130195607203664) q[7];
u3(1.01457317582231,-4.41634631974039,0.604134745983083) q[4];
u3(1.10768947515551,1.87271667627040,1.12246456717248) q[10];
u3(0.494679450129377,0.117524612280500,-3.18069407779022) q[6];
cx q[6],q[10];
u1(3.04545110504564) q[10];
u3(-1.76516111474900,0.0,0.0) q[6];
cx q[10],q[6];
u3(0.475594988704985,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.737298308118832,3.05515553733209,-2.69962243599738) q[10];
u3(2.29266423133653,2.23454644542116,0.398707049293794) q[6];
u3(2.21378968827093,1.52572421620719,-2.94792303659416) q[12];
u3(1.39263998811331,2.16781820559585,-2.55912305041897) q[3];
cx q[3],q[12];
u1(1.70513326445976) q[12];
u3(-2.65559590324038,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.864853171475661,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.69615551832743,1.52159709839392,-0.0986044200842189) q[12];
u3(0.661407486741247,0.897997041158977,4.25301061602653) q[3];
u3(0.579656340636831,-2.49490167856375,1.89931397675304) q[9];
u3(0.778509003126821,-2.81263721594884,2.26216606388440) q[0];
cx q[0],q[9];
u1(2.90382498848894) q[9];
u3(-1.90286976049470,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.37749500596583,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.533815290844843,-1.58678852382562,2.34024787040378) q[9];
u3(1.75053091284169,1.11882157492771,-3.48842434227090) q[0];
u3(1.63743473535144,0.872007142744967,0.695155887166625) q[1];
u3(1.07402483806039,-1.02301704198383,-2.78738695064529) q[0];
cx q[0],q[1];
u1(2.57921066138363) q[1];
u3(-2.31139136932044,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.527413828078891,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.78165263538223,2.58298668718898,-2.52731491319834) q[1];
u3(1.35119462609416,-1.11226119537745,-2.39430957665924) q[0];
u3(1.28920736526318,0.368955350826172,-2.15306747832696) q[3];
u3(1.67592292263293,-5.17898544313358,0.664037599796173) q[2];
cx q[2],q[3];
u1(1.93247645801832) q[3];
u3(-2.84526887164381,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.04903889485266,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.849504597448991,0.273993386941036,-0.398981290799188) q[3];
u3(2.48442649383844,-2.81317573581869,-1.31791644870678) q[2];
u3(2.18508228896812,-0.477384704724189,-1.56925151231745) q[12];
u3(1.87062925379509,-3.89571134406737,0.767038311246839) q[6];
cx q[6],q[12];
u1(1.07963443639802) q[12];
u3(-3.23416450774202,0.0,0.0) q[6];
cx q[12],q[6];
u3(2.37909900760011,0.0,0.0) q[6];
cx q[6],q[12];
u3(0.981303456714355,-0.918423883328017,-0.582171437988211) q[12];
u3(1.38356056579196,2.20217101993097,1.10256702320301) q[6];
u3(2.14750906402931,-0.0802211403132873,-1.19532456794903) q[8];
u3(0.906132170393782,-4.14678611463714,0.541200550620559) q[9];
cx q[9],q[8];
u1(1.43441214753791) q[8];
u3(-0.524200602755017,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.98511212217405,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.91126167533085,-0.685555316743096,-0.499627332521394) q[8];
u3(0.828294182888972,1.00514368268453,0.889622433526965) q[9];
u3(1.58332164417104,0.105764091274127,1.42022961100607) q[10];
u3(2.57552844181190,-0.818985488525619,-2.12948114837985) q[4];
cx q[4],q[10];
u1(3.36421094744022) q[10];
u3(-0.959800666937888,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.55883632476637,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.85893766977472,-1.00683655139848,-0.591289599190429) q[10];
u3(1.39129383194751,-4.56404901806786,1.59951818918226) q[4];
u3(1.28473958179257,-1.85664981191205,3.57293042434921) q[7];
u3(1.94057319539467,2.06962010184336,-1.82670899959493) q[11];
cx q[11],q[7];
u1(0.878113772688849) q[7];
u3(-3.37647168273041,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.94310900336113,0.0,0.0) q[11];
cx q[11],q[7];
u3(0.215331847718049,2.77087756730374,-1.09017158829636) q[7];
u3(1.21996556709149,-2.79470451158400,1.35469148806261) q[11];
u3(1.63191363333885,2.87600049822507,-0.450407922671802) q[1];
u3(1.87651070882889,2.66126651062933,-1.50827805790700) q[12];
cx q[12],q[1];
u1(1.35049386169229) q[1];
u3(-0.521601402856670,0.0,0.0) q[12];
cx q[1],q[12];
u3(2.71328446526277,0.0,0.0) q[12];
cx q[12],q[1];
u3(2.03701598680019,0.325373699125753,3.39829090674753) q[1];
u3(0.195057582252716,-3.04500141470539,1.73609888436698) q[12];
u3(1.65694742198573,3.95982755458413,-1.26745182678039) q[7];
u3(1.22208375845273,1.90804244244218,-0.929198390058955) q[11];
cx q[11],q[7];
u1(1.56744613687624) q[7];
u3(-2.71528331223352,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.16060568475584,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.23668954151901,-3.77834396724507,1.63387342201238) q[7];
u3(2.43099065175525,3.21652744006546,-1.88326212847556) q[11];
u3(1.12835950052599,-0.959497253726801,1.09207042887404) q[6];
u3(0.291333031962851,-2.26601182750815,1.21070887330788) q[5];
cx q[5],q[6];
u1(-0.315413198083311) q[6];
u3(-2.02828557083280,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.993197558552828,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.88820085088036,0.749051227129306,-0.816749170886942) q[6];
u3(2.29053134521900,-3.71854960516105,-2.33239531418029) q[5];
u3(1.96624087300198,3.05498992518320,-0.264345397625674) q[2];
u3(1.37216673146782,-0.128364019076322,-2.30213196687929) q[8];
cx q[8],q[2];
u1(1.50522768378181) q[2];
u3(-2.77249367781500,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.271493903684040,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.32230986470709,-0.791011590081317,-0.146352497650091) q[2];
u3(1.15770312223556,-0.661919886635100,1.39772223929903) q[8];
u3(2.16493569577361,1.97051779354250,-2.11957932916367) q[9];
u3(1.91018905760460,1.45274977040277,-2.96860502300326) q[10];
cx q[10],q[9];
u1(1.67987309404382) q[9];
u3(-3.11366918544550,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.14320857602408,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.48062571879769,2.26348885020547,-1.54302479836161) q[9];
u3(0.940896649843761,0.820456850356390,-0.416401601086467) q[10];
u3(3.09960702037514,-2.50941697293082,0.950536345908941) q[3];
u3(2.07019463339731,1.97579442699942,2.85160468839365) q[0];
cx q[0],q[3];
u1(0.827012252960382) q[3];
u3(-1.26560774178643,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.75360697017512,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.66239183223075,2.01304439953389,-2.16586534476826) q[3];
u3(0.554223284555301,4.72822565474509,-0.611973541979345) q[0];
u3(1.33008968660069,-1.37292331178764,0.454423666564792) q[9];
u3(1.74859444603213,-2.62667009212763,0.162943623399173) q[4];
cx q[4],q[9];
u1(0.769581236325131) q[9];
u3(-1.29987008733710,0.0,0.0) q[4];
cx q[9],q[4];
u3(-0.426919191525287,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.506239331484489,1.12415418677263,-1.93310437853554) q[9];
u3(2.95783762652383,-1.63036654261038,-4.60465267973476) q[4];
u3(1.66669764895372,-1.58576790523421,-0.0848395481366976) q[11];
u3(2.21556958680115,-3.96472617139380,1.02014701206249) q[1];
cx q[1],q[11];
u1(-0.679111898228320) q[11];
u3(0.216339531408495,0.0,0.0) q[1];
cx q[11],q[1];
u3(4.16706250790148,0.0,0.0) q[1];
cx q[1],q[11];
u3(2.92424205476566,-0.628895476683135,0.466683461038918) q[11];
u3(1.97429847073065,-1.55686676606471,-2.25214990066170) q[1];
u3(1.78890082013625,-1.60527390503532,-0.377327818517406) q[0];
u3(1.97194207632149,-3.18149512674927,-1.04292916473126) q[12];
cx q[12],q[0];
u1(1.34050775744886) q[0];
u3(-3.48417872004016,0.0,0.0) q[12];
cx q[0],q[12];
u3(2.43013482725958,0.0,0.0) q[12];
cx q[12],q[0];
u3(0.686739121488259,-1.80261542771755,-1.02689900074570) q[0];
u3(1.72141968199134,2.60652373278407,1.03433498115170) q[12];
u3(1.95705018431665,-0.558200024386882,1.53187971050532) q[10];
u3(1.93704757183972,-1.98508775216205,-0.525617977618680) q[3];
cx q[3],q[10];
u1(1.34868598779206) q[10];
u3(-0.722588201561747,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.95684811209345,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.81699503876121,3.07630016601943,-2.50749678131341) q[10];
u3(1.31828858793474,3.33478762856326,1.10500442135487) q[3];
u3(2.15473628949779,-3.03084765136968,2.48182746340829) q[8];
u3(1.02729876486449,2.94930894572711,-1.02039697321338) q[5];
cx q[5],q[8];
u1(1.56213250139502) q[8];
u3(-0.876744572315015,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.595714584446684,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.55741594582870,2.07485025706345,-3.28834574423428) q[8];
u3(1.40104472855097,-1.80618234781652,-2.45470486655978) q[5];
u3(1.97137423563263,-1.42766228999643,-0.569584673758890) q[6];
u3(1.32656832662669,-3.70495032047730,-0.442778495706022) q[7];
cx q[7],q[6];
u1(4.29348150703214) q[6];
u3(-3.78105928969741,0.0,0.0) q[7];
cx q[6],q[7];
u3(-0.411033172615734,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.15811836545362,1.10104802390082,-2.40990143088884) q[6];
u3(1.01313531892681,2.98051320704429,1.59543546772371) q[7];
u3(1.48757036554396,0.396758466204184,1.68467227089558) q[9];
u3(1.96338070273895,-2.07743482477612,-0.371794652281650) q[1];
cx q[1],q[9];
u1(1.33628294409511) q[9];
u3(-2.34826279712195,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.212585314328675,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.44518403468737,2.08769989401654,-0.495318523964332) q[9];
u3(1.85014290319165,2.91443930480894,2.98207366759896) q[1];
u3(1.47264105322600,0.229886825055616,-2.16433649751854) q[7];
u3(2.18578549580544,0.486643701301626,-4.96557838558191) q[5];
cx q[5],q[7];
u1(2.47039384625713) q[7];
u3(-1.86486367970416,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.206152465250566,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.28546871283012,-1.85376468131967,-1.64267785625637) q[7];
u3(1.96798645362080,1.69965476148348,1.44110197595503) q[5];
u3(0.923464909350496,2.05472472454318,-2.94832586792105) q[10];
u3(2.25843268603545,-2.59424043022255,2.82570279476012) q[0];
cx q[0],q[10];
u1(3.31456943956797) q[10];
u3(-0.559524954281335,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.62818709167032,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.28044652805400,0.838907593727279,-0.205214989729199) q[10];
u3(0.514059877806452,-1.89139509888507,-4.30013979025855) q[0];
u3(0.587096278858133,-2.95782028915109,2.26203815760902) q[6];
u3(0.757168234571972,2.26937868688785,-3.77628372972475) q[4];
cx q[4],q[6];
u1(3.57824732464996) q[6];
u3(-0.816511307016620,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.73852235678600,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.11829807911106,-1.81488021944271,0.410746010253373) q[6];
u3(0.327765176741573,-0.540972982436793,-4.52989926872012) q[4];
u3(2.10650833580467,0.0728793468916764,-2.31827151234300) q[3];
u3(1.60439930708148,1.20515208891863,-4.07758406800360) q[8];
cx q[8],q[3];
u1(-1.14131451847378) q[3];
u3(0.125620047585312,0.0,0.0) q[8];
cx q[3],q[8];
u3(3.65405673159221,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.57672347225866,2.09529300067192,-1.39122963094651) q[3];
u3(1.34569215400902,0.158124005502572,0.291325650574748) q[8];
u3(2.15743930204735,1.76204814964576,-3.46523592881397) q[2];
u3(1.04205991053106,2.90138583670639,-2.58573917007273) q[12];
cx q[12],q[2];
u1(1.00663367242018) q[2];
u3(-1.40130841907874,0.0,0.0) q[12];
cx q[2],q[12];
u3(3.39373531938610,0.0,0.0) q[12];
cx q[12],q[2];
u3(2.11014573206613,1.75250868073556,-1.54248110029557) q[2];
u3(1.44518876999695,-4.52934806912044,0.213988650923927) q[12];
u3(1.15932754568130,-0.479066026294936,-1.43932981147021) q[12];
u3(1.49189382626509,-3.66460640567371,0.520162800336077) q[0];
cx q[0],q[12];
u1(3.43601116165818) q[12];
u3(-4.50041415049423,0.0,0.0) q[0];
cx q[12],q[0];
u3(-0.265351433237969,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.06573886841561,0.910916082034766,0.691778287031258) q[12];
u3(1.65240765672466,-5.68804235224082,0.240754886468874) q[0];
u3(2.66964512036092,-0.637251054385594,0.0144961300116295) q[8];
u3(1.87508167027682,-2.45413920153840,1.35630333925820) q[11];
cx q[11],q[8];
u1(3.57690015717180) q[8];
u3(-4.41964090504860,0.0,0.0) q[11];
cx q[8],q[11];
u3(-0.613128513178427,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.11623154463148,-2.87695795452803,1.46504347640616) q[8];
u3(1.42564531028629,-4.96763299397635,0.746023547493403) q[11];
u3(0.781586998838929,1.29273304139089,-3.41451918407684) q[2];
u3(1.79892449491345,-1.67128095076399,4.16473468880754) q[1];
cx q[1],q[2];
u1(0.868052216677997) q[2];
u3(-1.25012551650607,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.90925045209393,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.769746675475855,1.23462352422086,-3.91376302428880) q[2];
u3(1.15467267439447,-0.250060185282297,-0.863036332771376) q[1];
u3(1.68540586166035,1.57962960865525,-4.58266627287160) q[7];
u3(0.330742921067396,0.800468858420121,-0.0762309279775977) q[6];
cx q[6],q[7];
u1(2.13914586599859) q[7];
u3(-2.72006660884529,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.10628030302134,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.46794449503912,3.11323553531419,0.603160721285847) q[7];
u3(2.49790707135608,4.96676604405923,0.00439080200367448) q[6];
u3(0.832195781794298,0.608838933802775,-3.00509344534002) q[5];
u3(1.86240263991176,2.88565170789585,-2.50710597770313) q[4];
cx q[4],q[5];
u1(1.95773825066209) q[5];
u3(-2.27889154921690,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.44649838705960,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.98930336708968,0.883406429981504,-2.15057034105812) q[5];
u3(2.53978398610256,4.87512034333941,1.03947646811051) q[4];
u3(0.968863397909856,-1.04187268876176,1.09754532026105) q[10];
u3(1.00785518428022,-2.00875816815435,-0.0533395094329852) q[3];
cx q[3],q[10];
u1(0.598515857641180) q[10];
u3(-0.177372825965161,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.37021194243900,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.712013250761683,2.07907088130937,-3.98590620907946) q[10];
u3(2.94731315684426,-2.45211910087122,3.02477112225115) q[3];
u3(2.30275361317968,3.29895817980598,-1.18572832511614) q[5];
u3(1.88006644208403,2.14082375218171,-2.48602155932919) q[8];
cx q[8],q[5];
u1(0.829925624840419) q[5];
u3(-3.36372101277161,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.61084240186596,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.589566123718292,1.52735468662516,-2.10745871699589) q[5];
u3(2.35207566946506,-3.23302623833092,-2.62303202020700) q[8];
u3(2.20096677166603,2.78418283862576,-0.632944059977155) q[6];
u3(1.58134873301548,1.81375635984660,-2.04749436850705) q[7];
cx q[7],q[6];
u1(2.05490651461937) q[6];
u3(-2.80869044825597,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.800800620155968,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.00915396753065,-1.31587768306276,3.51968665054663) q[6];
u3(0.481644115944301,-1.68195102564935,2.16239982926872) q[7];
u3(2.09747329435349,-1.98052025607395,1.03525997978290) q[11];
u3(2.07187154985944,-2.89995906570253,0.385081495455726) q[4];
cx q[4],q[11];
u1(-0.641129029274983) q[11];
u3(0.256120988619889,0.0,0.0) q[4];
cx q[11],q[4];
u3(4.18053755432441,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.08775798703728,1.83370577808774,-0.454079078856392) q[11];
u3(2.83721386163481,2.46282248443337,-3.64075315324343) q[4];
u3(1.18645707775327,-1.33593518970106,0.763506611438884) q[0];
u3(1.42055309968589,-4.08172344139693,-0.545652130156472) q[12];
cx q[12],q[0];
u1(1.79550807354086) q[0];
u3(-2.37779347576636,0.0,0.0) q[12];
cx q[0],q[12];
u3(3.23966380272453,0.0,0.0) q[12];
cx q[12],q[0];
u3(1.65590365523967,3.90361205275682,0.303596609036237) q[0];
u3(2.61360135967439,4.17907592225588,-1.06694783655857) q[12];
u3(0.455189047630164,2.33281040978787,-2.10190437029851) q[10];
u3(0.608291674244088,0.922445363082417,-1.18653625761875) q[2];
cx q[2],q[10];
u1(3.46317668153940) q[10];
u3(-0.901957318465137,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.70250631023946,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.63875777940186,-0.741096985318540,0.445919045752354) q[10];
u3(2.51830557017971,-2.53691952665770,-0.881026020265260) q[2];
u3(1.69077859539731,0.603042912846159,1.89897680981324) q[3];
u3(1.49638796698279,-1.10821900885947,-0.523947742668967) q[1];
cx q[1],q[3];
u1(1.52885710494505) q[3];
u3(-0.954545701582284,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.193013442199591,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.19999095722691,-0.382204081672097,3.07739313281331) q[3];
u3(2.03053017876987,0.714376194972640,-1.12875949023120) q[1];
u3(1.13733457865713,3.50756417159996,-2.04361191403633) q[7];
u3(1.89205934918018,1.71522159165948,-2.27255057533432) q[9];
cx q[9],q[7];
u1(-0.205982214972184) q[7];
u3(-1.69322434256768,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.809293333518940,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.53679725911060,0.758184773800215,-3.36584451623798) q[7];
u3(2.03091818119614,2.27478088523732,-1.45999801725047) q[9];
u3(1.97102265964665,1.81883102656760,0.542614630177308) q[3];
u3(1.49021068556702,0.676487077692598,-3.19540531508531) q[2];
cx q[2],q[3];
u1(1.61533727850034) q[3];
u3(-2.43575344615319,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.67836653147363,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.54947911381680,1.60232433064446,-4.06953619255389) q[3];
u3(1.76481094675101,1.01692597759672,5.14184142560721) q[2];
u3(2.39741802239010,0.935887589177230,-1.86651749904521) q[5];
u3(1.08561382830346,1.12599342969920,-4.07837350247283) q[4];
cx q[4],q[5];
u1(1.71425635687797) q[5];
u3(0.170708941315352,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.916122675475216,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.16583704849881,1.99029038970715,1.08012867242262) q[5];
u3(2.57664785201718,0.469308488164563,-1.17622730900575) q[4];
u3(0.710827698155055,0.332005662666379,-0.829436475606125) q[11];
u3(1.19072372856367,-3.91556860396386,0.639335576349616) q[10];
cx q[10],q[11];
u1(2.86737255843618) q[11];
u3(-1.81761747754211,0.0,0.0) q[10];
cx q[11],q[10];
u3(0.815115308338134,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.92878096627855,0.823109292675385,2.29309467266621) q[11];
u3(1.04059220177225,-4.11900604986612,2.15252410232596) q[10];
u3(1.46408890905874,0.416845174492829,0.191910694643093) q[0];
u3(1.20014543636466,-0.00739427055153785,-1.48173059341173) q[8];
cx q[8],q[0];
u1(-0.174223315296755) q[0];
u3(0.966811319508805,0.0,0.0) q[8];
cx q[0],q[8];
u3(3.67660589440544,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.41791620910223,-3.36513761850484,0.0640651913230768) q[0];
u3(1.69375526045985,3.51710464289448,-1.92622292300735) q[8];
u3(1.07380623548932,2.36622936203501,-2.19525595553849) q[1];
u3(0.337323026767561,-2.71677742363628,1.86798016537487) q[12];
cx q[12],q[1];
u1(0.684924316096173) q[1];
u3(-0.197120277834846,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.57054763901891,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.72190288997002,-1.64997104026372,3.24052175656036) q[1];
u3(1.23872650732091,-3.04782937534174,-1.49667876305515) q[12];
u3(2.54810430900153,1.22845165435500,-1.48825058920572) q[7];
u3(2.07378856583336,4.02439158941296,-0.562006066431765) q[3];
cx q[3],q[7];
u1(3.41449494988724) q[7];
u3(-1.48270168193157,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.68319950063087,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.62296549221916,-0.405428936097865,-0.486304936079619) q[7];
u3(0.893189161838108,-1.52735647722486,-2.15511778139796) q[3];
u3(1.91095518510531,-0.323416657349299,1.39315364509859) q[12];
u3(2.02392161118922,-0.209148834493813,-1.51342266217829) q[2];
cx q[2],q[12];
u1(2.98181138695627) q[12];
u3(-0.640748164247723,0.0,0.0) q[2];
cx q[12],q[2];
u3(1.43662088679970,0.0,0.0) q[2];
cx q[2],q[12];
u3(0.921405282138111,-3.09384985539225,2.13907389659010) q[12];
u3(1.28918478544461,-1.19418847359453,-1.05214699808474) q[2];
u3(2.89464239606599,-0.928617163968676,-1.59884998455235) q[11];
u3(1.46977075380911,-4.84623113562524,0.461567668843064) q[9];
cx q[9],q[11];
u1(2.49400860843263) q[11];
u3(-2.00829525382266,0.0,0.0) q[9];
cx q[11],q[9];
u3(1.33630951874328,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.85746905772279,2.54074850532810,-1.28217943115504) q[11];
u3(1.33839770915236,0.822813793967909,2.11594441659108) q[9];
u3(0.0603030525075279,1.34730736803830,-0.797658433644040) q[5];
u3(0.867933848906857,-2.44818442665600,1.51802256619696) q[8];
cx q[8],q[5];
u1(2.36435842143207) q[5];
u3(0.00143505646524300,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.02153921024675,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.13061203492662,1.06349112401363,-1.23206266766178) q[5];
u3(0.563992119360485,-0.0935028723502627,1.53674621886867) q[8];
u3(0.978391133058682,-1.92063636249196,1.09238438676525) q[4];
u3(0.163625134873533,-2.68068712270873,1.31839758981520) q[10];
cx q[10],q[4];
u1(0.994563238779299) q[4];
u3(-0.104097185191869,0.0,0.0) q[10];
cx q[4],q[10];
u3(1.43940270836081,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.752894829583402,-2.51651124188646,3.43637209325362) q[4];
u3(1.73036547206897,1.77999519176152,-2.80416857140213) q[10];
u3(2.10325162163593,0.281014186770340,1.54150779005578) q[1];
u3(1.97898333079677,-0.887631030310228,-1.09979090109411) q[6];
cx q[6],q[1];
u1(2.73625528146938) q[1];
u3(-2.84742599899784,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.85944371167599,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.40037921784437,0.493190253873897,0.353219458046034) q[1];
u3(1.58723764016414,2.38621745859715,-3.40354160289382) q[6];
u3(2.71524735594162,3.42529045297700,-2.47332562819956) q[11];
u3(0.638237853114424,-0.0509363455754519,2.18260429633006) q[9];
cx q[9],q[11];
u1(1.55867303055411) q[11];
u3(-0.808514950421616,0.0,0.0) q[9];
cx q[11],q[9];
u3(-0.427481484907037,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.37600083725963,-0.237682294057226,-1.80673047344357) q[11];
u3(1.72293597685141,1.93841379093238,1.30192117508963) q[9];
u3(3.08577991201664,1.92036399371268,-0.893787149925937) q[4];
u3(1.80986571183093,-1.07541744037989,-4.71538852283585) q[1];
cx q[1],q[4];
u1(0.976844224405918) q[4];
u3(-0.514778991371158,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.00247679142490,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.29220309262836,0.128594533066354,-0.0906996748616059) q[4];
u3(2.44613457451461,-0.389476988842040,-0.238113642213768) q[1];
u3(1.12833442018762,0.219035301983685,1.83320529990423) q[8];
u3(1.75052289963298,-1.80812456202303,-0.652440794477110) q[7];
cx q[7],q[8];
u1(1.91165448865818) q[8];
u3(-0.439007473397922,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.56809298995050,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.539859075711939,-1.08210240534833,0.503381028600867) q[8];
u3(1.13015646465816,-0.314136734178013,-0.736067242851747) q[7];
u3(1.70103216594005,3.03648330125591,-1.32913682721448) q[12];
u3(2.16369930104309,0.979171308636553,-2.75316804847355) q[2];
cx q[2],q[12];
u1(1.80994425038303) q[12];
u3(-2.38917136382482,0.0,0.0) q[2];
cx q[12],q[2];
u3(3.40884240605323,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.41728153692883,0.473954527602653,-2.32534468814585) q[12];
u3(1.24381164457009,-1.69540284920740,-1.76555169334061) q[2];
u3(0.800544771988614,0.104934925654221,0.997265958479029) q[10];
u3(0.984404017199908,-0.270982635112437,-2.01943528207364) q[5];
cx q[5],q[10];
u1(-0.592823160195465) q[10];
u3(-2.01193038086565,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.61080069174858,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.611090732683708,0.822762526376702,-0.880968744964254) q[10];
u3(2.38025710931071,0.882943462209826,3.11332685808359) q[5];
u3(0.870940527288904,1.05164520741347,-0.910370012185108) q[6];
u3(0.356882379487137,-3.17789879991346,0.516562044336605) q[0];
cx q[0],q[6];
u1(0.817546249154591) q[6];
u3(-3.28706875235376,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.80719249157249,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.28253391155668,1.66213499848987,1.58809024118669) q[6];
u3(1.27672607199465,-0.470049973363027,1.05472331848031) q[0];
u3(1.63638852518591,0.601639664231427,1.67959200251040) q[2];
u3(1.52363507051013,-0.760462557877260,-2.05518015096735) q[1];
cx q[1],q[2];
u1(1.68468898127351) q[2];
u3(0.483737001674849,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.931787633400904,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.30263563165980,-2.59268559309749,0.238688132742815) q[2];
u3(2.15840093804448,-0.357470076004667,-3.90709888849952) q[1];
u3(1.43018149029035,0.570772802579902,1.00195511094620) q[4];
u3(1.20513598986337,-0.306053720162607,-2.23174687682913) q[3];
cx q[3],q[4];
u1(0.857327464066782) q[4];
u3(-3.29198871130306,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.47261790827340,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.80727585105969,-0.753183542066326,3.03833207908740) q[4];
u3(1.90558208322912,2.25717644766070,2.73429000620888) q[3];
u3(2.63218539815234,1.23871305145573,1.37438901603499) q[0];
u3(1.10093448035774,-4.88711904906433,-0.631016265301533) q[7];
cx q[7],q[0];
u1(0.421032201891166) q[0];
u3(-1.45793461879199,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.07694643287025,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.608626065096295,-0.137449336061661,1.20864368906520) q[0];
u3(0.820377873129545,-0.600609641616924,2.29894651509023) q[7];
u3(1.19675519567196,-0.408880228870519,1.19129926012800) q[6];
u3(1.42755525407620,-1.74364478930174,-2.30977000662331) q[9];
cx q[9],q[6];
u1(0.751113015771586) q[6];
u3(-1.34635295030128,0.0,0.0) q[9];
cx q[6],q[9];
u3(-0.0440255106874927,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.70483023008938,0.368912709947679,-3.75837397943463) q[6];
u3(2.56806046336220,2.25720968728827,3.74703411124512) q[9];
u3(1.25907725776272,0.629371061590173,-2.04574649039823) q[10];
u3(0.984017182206336,-4.11508827888050,1.51527017887362) q[8];
cx q[8],q[10];
u1(3.57955254828522) q[10];
u3(-1.54177221086596,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.31154055137676,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.79404138468610,-0.915269912070645,1.11710781597523) q[10];
u3(1.73425949399849,-5.05260848191643,0.288211757980051) q[8];
u3(0.908656887795031,3.64277592233185,-1.86816768324744) q[12];
u3(1.57348335861679,1.18686687258845,-2.72725077273283) q[11];
cx q[11],q[12];
u1(1.48367382535012) q[12];
u3(-2.93092550832572,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.957924235138716,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.99551316711831,-1.70839826320419,1.29671293055672) q[12];
u3(0.846685410786895,-2.08065495418502,1.49810164946580) q[11];
u3(1.25305389615671,2.64628945719650,-1.38791331686309) q[0];
u3(2.51472802727111,1.42784739617629,-1.10331319611638) q[11];
cx q[11],q[0];
u1(2.57432772267745) q[0];
u3(-1.98264333383743,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.842611603099441,0.0,0.0) q[11];
cx q[11],q[0];
u3(2.70319363189680,-1.80304155376007,3.42632761190823) q[0];
u3(1.48583062179610,2.42406373414579,2.09827644227466) q[11];
u3(1.16551931168863,-1.46227238764094,0.986344938861406) q[3];
u3(0.524905525070964,-3.74247794798232,1.51133348328706) q[6];
cx q[6],q[3];
u1(3.17621517525046) q[3];
u3(-1.66363814861575,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.785404949311253,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.262655139729384,3.67721682178459,-1.55569680943695) q[3];
u3(1.89218640690174,0.121683242545219,2.65190854899148) q[6];
u3(2.54956574898519,0.466549705393398,0.554042298969426) q[2];
u3(1.09219170782920,-2.01261724020499,-2.42972861585489) q[12];
cx q[12],q[2];
u1(2.54408312718414) q[2];
u3(-3.04971769231796,0.0,0.0) q[12];
cx q[2],q[12];
u3(1.39998937181483,0.0,0.0) q[12];
cx q[12],q[2];
u3(0.696354510294119,4.22066067934904,-0.163547259188033) q[2];
u3(1.97680973944411,-0.464963691131681,2.14162742346996) q[12];
u3(1.05299614619977,-1.59925253997974,-0.818459248840262) q[10];
u3(1.35241111830928,-4.11331989228267,-0.177151827268486) q[5];
cx q[5],q[10];
u1(1.82197312301179) q[10];
u3(0.0744699433248792,0.0,0.0) q[5];
cx q[10],q[5];
u3(0.944165155525656,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.419838833497887,-2.76659143338683,1.73388371610086) q[10];
u3(1.51301635713089,1.34371086182578,4.75022173756883) q[5];
u3(0.747294541425637,-0.418914345154208,0.797150733868873) q[7];
u3(1.33088736153743,-0.774907240514674,-2.13439848283566) q[1];
cx q[1],q[7];
u1(2.60638100887574) q[7];
u3(-1.81495995829532,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.131326645170106,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.92753797295808,2.27067194004964,-3.89285639358359) q[7];
u3(1.88618547307052,2.60819102110848,-1.32078860321849) q[1];
u3(2.37119267201119,1.73773496413063,-3.42120110248768) q[4];
u3(1.05806707473506,3.27389180631758,-2.92508059426650) q[9];
cx q[9],q[4];
u1(0.248461678267829) q[4];
u3(-1.04139374927146,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.32696668999792,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.339632255113602,-3.52839896835354,-0.131450063437236) q[4];
u3(2.69008368724247,3.42363377711856,2.30263456594232) q[9];
u3(1.53192608695161,3.85138363606931,-2.05895181088265) q[2];
u3(0.146026881560993,-0.547031693372013,2.53115167409242) q[4];
cx q[4],q[2];
u1(3.80884304271587) q[2];
u3(-1.35552565952505,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.21545605584631,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.965776183849980,-1.20486783244505,0.644615081434832) q[2];
u3(2.05824452311089,-0.710128574142261,-0.726051788089835) q[4];
u3(1.75959524758146,-0.220772727233829,2.23751397741259) q[8];
u3(2.02918413467384,-2.59199400081902,-1.67240973225131) q[3];
cx q[3],q[8];
u1(1.40238399301235) q[8];
u3(-0.390440204443777,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.236494615263733,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.55708955077691,1.13173752497871,-1.93099151424687) q[8];
u3(2.27983619542435,0.787869757053919,0.525952336098966) q[3];
u3(0.589129343904975,-1.11752651340629,1.22312123681897) q[10];
u3(0.171025870366465,-3.29020498666291,2.73119835553992) q[5];
cx q[5],q[10];
u1(1.83381316771616) q[10];
u3(-3.18046966123219,0.0,0.0) q[5];
cx q[10],q[5];
u3(0.998441648225225,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.85283671961092,1.19462102350625,0.358421317247298) q[10];
u3(1.11087491925221,-3.02547286943754,2.95691746929924) q[5];
u3(2.36529722878650,-1.36247848913187,-0.198837373878435) q[9];
u3(2.19174741527869,-2.73744962184822,0.723723299913900) q[7];
cx q[7],q[9];
u1(2.67916146936116) q[9];
u3(-1.80376648766908,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.723083253463999,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.49364625868370,-1.15313490639377,4.81354790386926) q[9];
u3(2.40997601653847,2.06040836413262,0.126500381432446) q[7];
u3(0.941124799893081,-2.31264700096438,0.938312151235777) q[11];
u3(1.04198722407216,-2.86757465446184,-0.562498091306870) q[6];
cx q[6],q[11];
u1(2.15177415171188) q[11];
u3(-2.78792185613966,0.0,0.0) q[6];
cx q[11],q[6];
u3(0.721255287474237,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.24899029632411,-1.52099869434871,1.27779899430701) q[11];
u3(1.14041554958965,-1.08579446740957,-2.22247755206105) q[6];
u3(1.13532165893095,-0.598042008605218,0.824035725976204) q[0];
u3(1.07180592258111,-1.74490499222194,-1.19306397464308) q[1];
cx q[1],q[0];
u1(3.37878944517532) q[0];
u3(-0.674264605338273,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.93628916881845,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.13002672733294,-2.57306932642949,1.17052024022129) q[0];
u3(2.15482343308941,-2.68659526261989,1.92214320094660) q[1];
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
