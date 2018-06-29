OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.66012516613248,0.504675389665770,-3.12903138511215) q[1];
u3(0.353604802839387,-2.73020107727956,2.82674405749080) q[2];
cx q[2],q[1];
u1(0.100759256783784) q[1];
u3(-1.41579738948508,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.66604124442718,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.0942616007786818,3.05665057687769,-1.55801719698263) q[1];
u3(1.30248497812937,2.33340175917485,-1.77293060594026) q[2];
u3(1.92568738689889,-4.47079853144333,1.67452071335121) q[0];
u3(0.244755415792866,-2.38191405431890,2.44921489813977) q[3];
cx q[3],q[0];
u1(3.14724452166969) q[0];
u3(-2.22590337510435,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.29399769102300,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.50260079314666,0.248620108135051,2.00101753097099) q[0];
u3(1.70669708026538,1.66380839818454,3.29746939687885) q[3];
u3(0.933338580983794,0.551757021705529,-0.879114768684588) q[2];
u3(0.380797768506054,1.16331934270411,-4.17522451066281) q[1];
cx q[1],q[2];
u1(2.32706600171152) q[2];
u3(-2.75140024647865,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.35870831575213,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.98168794037108,0.330370458187312,-1.77263688197714) q[2];
u3(1.59219965275311,0.794534501457387,-3.69883890875390) q[1];
u3(1.74060054777460,1.28312630673256,-2.44678159831307) q[0];
u3(1.07068245220653,-3.29634816877098,2.47225924622154) q[3];
cx q[3],q[0];
u1(-0.213274285204622) q[0];
u3(-2.22986139892264,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.988465863414198,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.50216948264995,-4.13891577542466,2.13391931893810) q[0];
u3(0.684343232873377,-1.87401285493681,-0.718118238431515) q[3];
u3(1.20255637876673,0.805162342889456,-3.79930571516784) q[2];
u3(1.74331355658916,4.83307904287056,-1.13718774625956) q[0];
cx q[0],q[2];
u1(3.21454837383989) q[2];
u3(-2.48257507966590,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.16023144816936,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.72991394832790,-0.104367828836181,-1.94055266827130) q[2];
u3(0.724087417419215,-4.06701962028557,-1.54180576152317) q[0];
u3(1.27244907367085,1.68685817010542,-3.41170374321946) q[3];
u3(1.77662872485459,1.81489009865370,-3.46381541001189) q[1];
cx q[1],q[3];
u1(2.15702577270687) q[3];
u3(-1.75674235314104,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.33944347169143,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.918261466045649,3.29734422523348,-1.12341684709128) q[3];
u3(1.48714972982515,1.41212205328297,3.91854306738008) q[1];
u3(0.809654591461539,1.40774920414662,0.0971239169243915) q[3];
u3(1.13770107949845,-0.525433150762590,-2.35701208803412) q[2];
cx q[2],q[3];
u1(1.82167444039131) q[3];
u3(0.108556198471127,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.710540989907041,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.89349138135758,-2.40894042338340,1.39078094588211) q[3];
u3(2.21894823879956,-0.515795977131509,2.35874173325481) q[2];
u3(1.79278727990581,0.506737448050887,1.89304286295015) q[0];
u3(1.33212082750388,-0.487005259049795,-1.28878182333382) q[1];
cx q[1],q[0];
u1(1.38315933823334) q[0];
u3(-0.732867371159567,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.94338237196182,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.74380097967012,1.28151853878023,-4.32854551467332) q[0];
u3(1.38250618134465,0.261161059877417,-2.49980417431800) q[1];
u3(1.13089793859835,2.36379085302078,-2.74675143250659) q[3];
u3(0.425134387685266,-2.37304024585772,2.17863698227871) q[0];
cx q[0],q[3];
u1(1.45661521805125) q[3];
u3(-3.17362597768054,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.97836467284640,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.12386470313663,-0.183183365516381,-1.19885930849522) q[3];
u3(1.79383960959308,4.37764342469293,-1.57328749859362) q[0];
u3(1.71983429750222,-1.00704773237366,-0.650277172077876) q[2];
u3(1.01671227319074,-4.29285848370426,0.329270827703301) q[1];
cx q[1],q[2];
u1(2.96860574452818) q[2];
u3(-2.31331654340726,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.32614545799983,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.26317039456496,-2.67091628040685,1.86484917316391) q[2];
u3(1.67186110154203,-1.22337445357019,3.93908894687927) q[1];
u3(2.57545491787610,-3.02557123202924,3.12356188650725) q[1];
u3(1.15937196020928,3.37220943415513,-2.18088452477105) q[3];
cx q[3],q[1];
u1(-0.402880520283393) q[1];
u3(-1.91345141541824,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.861011677239727,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.51273904801067,-2.04001778591757,-0.643124399259648) q[1];
u3(0.734458909188141,-0.306240603429308,0.292261811963133) q[3];
u3(2.17330748360269,0.319524849109126,-0.832387383479144) q[0];
u3(1.42607691121305,0.683132061040458,-4.51464902117139) q[2];
cx q[2],q[0];
u1(4.08850166581996) q[0];
u3(-3.26185915019182,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.612966272103638,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.685370553487339,0.0696166963237059,-3.49078934233324) q[0];
u3(0.486691488732447,-4.37019393184746,0.386997771230606) q[2];
u3(1.50985300066999,1.32055790861213,-0.815453412743335) q[1];
u3(0.444110139434397,1.42961368517637,-4.49558927061671) q[3];
cx q[3],q[1];
u1(-0.0549823741492945) q[1];
u3(-1.49778234718753,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.62933841269344,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.51869751064242,-2.33244351160966,2.50556139608716) q[1];
u3(2.61249120681738,0.284874286353037,-3.07081213517274) q[3];
u3(1.41124398891345,-0.754142781032543,-1.35403672951493) q[2];
u3(2.65202069131818,1.47481868666454,-3.66204693903097) q[0];
cx q[0],q[2];
u1(1.37891826570333) q[2];
u3(-3.31849366057531,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.99230723182927,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.29360514098820,3.51651706836091,-1.34450894561215) q[2];
u3(0.786473653763886,2.25418103965698,1.00424565160632) q[0];
u3(0.155942994142579,1.11497039547618,-1.95058900977687) q[3];
u3(0.734324220262104,-1.06902316049033,-0.122755032368314) q[2];
cx q[2],q[3];
u1(-1.18432150902742) q[3];
u3(0.120608492025492,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.33208400297061,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.966101646463989,-0.151623841441289,-0.554917421676167) q[3];
u3(2.42424025245579,0.170582213210578,-1.37852457979828) q[2];
u3(2.44550214180411,-0.562595330679488,-0.158256966956218) q[0];
u3(0.990537697061112,0.576186903339965,-5.67202735081973) q[1];
cx q[1],q[0];
u1(1.22785948583061) q[0];
u3(-3.19117859439823,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.48260100981923,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.01881780195133,-1.79501404802242,1.11606290694298) q[0];
u3(1.10377139942512,0.259339286885321,-0.805619004535451) q[1];
u3(0.731353225964229,0.0947488662067937,0.0415361595541010) q[1];
u3(0.887465497159234,-2.98515791000597,0.763743931811696) q[3];
cx q[3],q[1];
u1(0.831341422166618) q[1];
u3(0.140353293923037,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.96890008932536,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.990733021763644,2.53407797756511,-0.849746379069826) q[1];
u3(2.52790391028215,-1.31143213203761,-2.52442780199847) q[3];
u3(0.822536103905558,0.343014523619533,-1.02853591793899) q[2];
u3(0.186159037691187,-2.08703237900848,0.132515220820858) q[0];
cx q[0],q[2];
u1(2.72703557867981) q[2];
u3(-1.49920691149945,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.0432290186061555,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.91472228350989,-1.02649487603791,0.311118184360564) q[2];
u3(1.65755581031223,1.82663673432142,2.47770661347999) q[0];
u3(1.94304595872374,1.80029948573536,-3.35363574767905) q[1];
u3(2.28844812175739,-2.04832480198857,3.76910099832340) q[2];
cx q[2],q[1];
u1(1.96883803000991) q[1];
u3(0.308171012479562,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.46079890949849,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.30399029457242,-1.05148251879845,0.0806938027731995) q[1];
u3(2.96694401569018,-1.77783744461622,-1.12050676443263) q[2];
u3(2.76945175577387,1.32338991992759,-1.42903525879915) q[0];
u3(2.16045454862393,2.27107344706342,-2.21489231106028) q[3];
cx q[3],q[0];
u1(3.64332724454290) q[0];
u3(-4.41419173588356,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.285238970646799,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.26910179191577,4.06858524224228,-1.62291942345761) q[0];
u3(1.02113127828712,-0.277070636509283,-2.73709828380626) q[3];
u3(0.909866163437018,1.22607082074022,-0.321786853880972) q[1];
u3(0.148107420230863,-3.44981284032906,2.07008507353252) q[3];
cx q[3],q[1];
u1(0.880006057942137) q[1];
u3(-1.48655265691055,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.62491405584910,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.73787749249280,0.0781854992817352,0.315816366651802) q[1];
u3(1.21170349978900,1.23904079233187,2.05071425371048) q[3];
u3(0.589961169215762,-0.281425754722971,-0.741746215050411) q[0];
u3(1.21108802657815,-2.75352359481402,1.34254599095521) q[2];
cx q[2],q[0];
u1(1.45820384124028) q[0];
u3(-0.362847350746472,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.115225811488981,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.39761998740678,-1.68509729284601,1.20001730287191) q[0];
u3(1.31315579357987,-1.64429407951362,-3.23095565881163) q[2];
u3(2.27012598160295,3.10032772064901,-1.76959017792292) q[2];
u3(1.36520895052263,1.71571340579149,-2.54236566766970) q[0];
cx q[0],q[2];
u1(1.66293366896413) q[2];
u3(-2.46523141444716,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.10988493540618,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.94095286149350,0.0617222655448365,-1.17878099832388) q[2];
u3(1.75206331738772,0.551136144848219,3.96306338238848) q[0];
u3(1.98963714513514,-1.20993558356204,1.44013525130643) q[3];
u3(1.32789128197736,-1.50502382326713,-0.876275903998362) q[1];
cx q[1],q[3];
u1(1.56819820031156) q[3];
u3(-3.54968740818821,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.64677821786377,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.16373632470740,-0.0791975857864649,-0.263189879487193) q[3];
u3(1.66802632037813,-1.28466085934265,-3.53788259155199) q[1];
u3(1.99972270933632,-1.92516784333328,2.89199443168149) q[2];
u3(2.60387947033723,-2.80354804118882,1.79274921912746) q[0];
cx q[0],q[2];
u1(1.30516726203134) q[2];
u3(-1.45135546215507,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.736897539112418,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.88001733079389,-0.590804793673175,3.43661927168749) q[2];
u3(0.882874371564491,-0.386057545937581,3.64403717251770) q[0];
u3(1.55161114417611,-0.294963457735836,2.15012490486078) q[1];
u3(1.55056711114190,-1.78476309561939,-1.57098110291711) q[3];
cx q[3],q[1];
u1(0.490008056624036) q[1];
u3(-0.548314788259637,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.86196045253910,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.27351637447696,-2.11503099120763,-2.00555169913531) q[1];
u3(1.23542616199132,1.78767091279307,-2.90214550287846) q[3];
u3(2.33863875361699,0.359070457404732,-3.35613776308147) q[0];
u3(2.46978244004101,-0.686762675013017,-4.60220263060772) q[2];
cx q[2],q[0];
u1(1.13997050465665) q[0];
u3(-0.667053875068739,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.173352408955315,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.34571682755221,-1.83417719361984,3.06858608709033) q[0];
u3(0.902741500859274,2.49530958776080,0.751488220639218) q[2];
u3(0.361846557880510,0.780266582404981,-0.661032736353594) q[3];
u3(0.943233086907051,0.832427381270346,-1.77848733077205) q[1];
cx q[1],q[3];
u1(1.57005461266550) q[3];
u3(-2.53888993785001,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.164238556469919,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.21168674057223,-1.66682895917727,3.32459666199340) q[3];
u3(1.20768070381115,2.29364092484116,-1.62440335368244) q[1];
u3(1.89673230096938,4.00765621148683,-1.40959416378049) q[0];
u3(1.77483278793270,2.20960990007110,-0.164890307447519) q[2];
cx q[2],q[0];
u1(-0.558212251998061) q[0];
u3(-2.04074676151146,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.63717109265695,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.900393613299430,2.21349420964170,-0.937448758157400) q[0];
u3(1.81116936809490,-2.48866100384797,-0.530275261226668) q[2];
u3(1.90662900867998,-0.786476747736273,1.10538829067550) q[1];
u3(2.44253757814413,-1.55400161330830,-2.97743225785771) q[3];
cx q[3],q[1];
u1(3.69946258257294) q[1];
u3(-1.57306181058721,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.10620582035809,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.41956768781177,-3.49601334616743,2.20544136227430) q[1];
u3(1.01243860678775,-0.557493329309906,0.721383165661609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
