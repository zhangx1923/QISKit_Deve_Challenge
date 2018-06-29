OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.56430747000454,0.945169219505827,-3.45189613446739) q[5];
u3(2.31766979543561,2.60027983188985,-2.59514493222576) q[4];
cx q[4],q[5];
u1(2.53294666652877) q[5];
u3(0.00445061878392194,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.64284823163762,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.68612370456636,-1.20726042482094,-1.71595739753176) q[5];
u3(0.821159272645567,4.02995191557461,-1.04424804412512) q[4];
u3(1.90884889897634,0.921869491906997,-0.805759022369143) q[0];
u3(2.07871999406946,0.687031041627943,-4.33995296873697) q[8];
cx q[8],q[0];
u1(1.14343308955915) q[0];
u3(-3.27890225530114,0.0,0.0) q[8];
cx q[0],q[8];
u3(2.34112512227689,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.95440735048311,0.0715990954520907,0.476355159025297) q[0];
u3(1.61700770075364,-0.0678572403838799,1.25013149538430) q[8];
u3(0.261926115532475,-0.194270616216869,0.557153836719779) q[7];
u3(0.472069703066350,-0.514917229211175,-1.93208305675083) q[2];
cx q[2],q[7];
u1(1.35120129135568) q[7];
u3(-2.97864054302650,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.82429391212733,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.39278634582477,2.97723009386590,-3.16316667936513) q[7];
u3(1.86143964363601,-0.610482347841009,5.61961741865047) q[2];
u3(2.23897421785697,1.15620040181117,-0.263209062192972) q[9];
u3(2.28736519554435,-0.452406796946714,-3.19929432564209) q[3];
cx q[3],q[9];
u1(2.10977428145135) q[9];
u3(0.652637473047318,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.65054546597921,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.34978550771618,-0.999977658267637,-0.870515258898941) q[9];
u3(2.07821681669106,3.44137800118251,2.19685567209525) q[3];
u3(0.681421959898446,-1.50122469107061,4.22913799775412) q[6];
u3(1.42841458985975,-0.734081009817987,0.693362698582291) q[1];
cx q[1],q[6];
u1(1.66185755796048) q[6];
u3(-2.85501263946475,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.921360299879532,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.62117937521105,-1.55458971712894,3.31946576857812) q[6];
u3(1.62579200939087,-2.14786254963921,0.820417724654084) q[1];
u3(1.22000692034086,-1.95018741083761,0.734124869592434) q[5];
u3(1.80339247545613,-3.99443651033880,0.528165231513388) q[1];
cx q[1],q[5];
u1(4.03475743334709) q[5];
u3(-3.39168767555960,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.590991886108169,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.914529298922581,-1.70685546750233,2.57531716945935) q[5];
u3(1.56107303912850,-3.09992806395124,-3.04471304532185) q[1];
u3(1.58978275371947,2.72858821966444,-2.01707422754891) q[0];
u3(2.06128733801516,2.52973386658002,-0.214488253410985) q[7];
cx q[7],q[0];
u1(1.67974747089242) q[0];
u3(0.529751214531038,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.911504131872564,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.88150982114824,-2.13846370213685,-0.148397951589675) q[0];
u3(2.35286073908936,2.92888336704373,3.18201036083798) q[7];
u3(0.236019917692715,-1.84755711327246,2.10081528998046) q[4];
u3(0.841654827603955,-0.258584993278341,-2.17985641489390) q[8];
cx q[8],q[4];
u1(-0.352112815251225) q[4];
u3(-1.51393380660783,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.93972239020804,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.88273252223618,-1.28224013823400,2.47338457089141) q[4];
u3(1.24673937333658,2.91098359185114,0.533676271292234) q[8];
u3(2.04852593021842,2.57021918773058,-0.339271587812146) q[3];
u3(2.42328257837026,0.813114246953619,-2.37124800495607) q[9];
cx q[9],q[3];
u1(4.32183741406773) q[3];
u3(-3.48280276597226,0.0,0.0) q[9];
cx q[3],q[9];
u3(-0.642367956603510,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.64320598591497,-1.22698204862188,-2.50842973442675) q[3];
u3(1.29599233169983,2.81024755972087,-2.10578167724041) q[9];
u3(0.186547927947264,-0.650619233191280,1.92137422130056) q[6];
u3(0.521809969419443,0.183685262812161,-1.92925813846404) q[2];
cx q[2],q[6];
u1(2.38093456215445) q[6];
u3(-2.94650421480560,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.64923070475204,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.85144711729628,-0.270167097921754,-2.22038086113262) q[6];
u3(2.31830096007269,-4.28270810371855,-1.22804605912005) q[2];
u3(1.70037278994652,-0.768828173828300,1.78429999046288) q[2];
u3(1.94058787154061,-1.63068763044352,-2.19906795934628) q[5];
cx q[5],q[2];
u1(0.778430418428817) q[2];
u3(-3.11686269667429,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.73761956235428,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.68423949833294,0.894940473308901,0.700314205330192) q[2];
u3(0.473392692242398,1.56537469580891,1.67718799282387) q[5];
u3(1.17889476954988,0.802252004058017,-0.915814001355692) q[4];
u3(0.771007321875100,-0.987989115612141,0.0663850143726281) q[9];
cx q[9],q[4];
u1(2.07006913609866) q[4];
u3(-1.73425995638289,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.0793422216376616,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.30216333113022,2.77706884670365,-2.02657863825994) q[4];
u3(1.79260475427581,-0.171486807803699,-1.61804227618349) q[9];
u3(0.914262823365285,0.450174000553149,2.17009175721577) q[6];
u3(1.20130215450395,-1.19755899453325,-1.64698618354964) q[1];
cx q[1],q[6];
u1(1.60272816845443) q[6];
u3(-2.64131098106919,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.849015322819376,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.734955302703305,0.781939234036581,3.29298501190819) q[6];
u3(1.33149268539426,-3.48545955476717,0.407211355009544) q[1];
u3(1.32876361431239,1.00255341444925,-2.67837732758123) q[0];
u3(1.63224179953499,1.98526605461595,-4.07387003029541) q[8];
cx q[8],q[0];
u1(3.01089047748337) q[0];
u3(-1.93894178873877,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.60063586732379,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.883944811389226,-0.346063768931118,-0.942894518216621) q[0];
u3(0.908273738174633,1.40662466713298,-3.69676871529607) q[8];
u3(1.67925302546972,3.30275780637048,-2.86721315932797) q[3];
u3(0.164492483437836,0.154737285234263,1.49288003102924) q[7];
cx q[7],q[3];
u1(3.70835837870716) q[3];
u3(-3.36550166197884,0.0,0.0) q[7];
cx q[3],q[7];
u3(-1.04808159290981,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.02091344586911,-3.62725172873178,1.70842981867759) q[3];
u3(1.05215387649884,5.06128828505306,-0.157689285973170) q[7];
u3(1.34592264801951,2.29001786618267,-0.840500978020300) q[9];
u3(0.976523136069219,1.38742942242375,-0.577986974530765) q[8];
cx q[8],q[9];
u1(0.449098099565034) q[9];
u3(-1.26945656554449,0.0,0.0) q[8];
cx q[9],q[8];
u3(3.12463249310863,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.11259920233931,-1.57638105852999,-2.21252068509775) q[9];
u3(2.80901644920653,0.837160454931689,3.09930286750792) q[8];
u3(0.639372255031338,-1.05023772111705,-0.617023589996939) q[5];
u3(1.70393766338774,1.20841880236456,-4.51786058385603) q[1];
cx q[1],q[5];
u1(1.58081370489404) q[5];
u3(0.252543664628282,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.08827568691748,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.892907304195795,3.91389511553978,-2.07146970630525) q[5];
u3(1.46400598676682,-1.54515260506674,-1.21417801925674) q[1];
u3(1.84621834215284,-1.33828861427157,-0.860295726683706) q[4];
u3(0.408533506642453,-3.07670879506026,-1.33426931425456) q[6];
cx q[6],q[4];
u1(0.591820693664259) q[4];
u3(-1.03191346354636,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.96235286179778,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.58588292077944,1.06619240433549,2.71980114562581) q[4];
u3(1.40158424365322,-2.96134842522329,0.681415213346099) q[6];
u3(1.64895372540019,3.82186931547771,-2.08008699887759) q[2];
u3(0.660881360806128,0.832984525353287,0.156996889163666) q[3];
cx q[3],q[2];
u1(1.42610619869926) q[2];
u3(-0.702970031211301,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.06625096356210,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.31474403066125,1.88330017951978,-3.54295710522002) q[2];
u3(0.208045073797831,2.27552140478524,3.61728809724816) q[3];
u3(0.943112226545694,-2.29828283792531,-0.398116245142775) q[7];
u3(2.13275207672300,-3.57747067544020,-0.0523824693263517) q[0];
cx q[0],q[7];
u1(1.38779402968271) q[7];
u3(-3.15173598186334,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.00999226239012,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.54583580364822,-3.02464694761095,-1.55021751879045) q[7];
u3(1.00309851116246,-1.85598873106174,0.863256677937551) q[0];
u3(2.94419926167567,1.79171191695901,-0.422863087201465) q[1];
u3(2.97363331811671,-0.242434847165284,-5.18734468194487) q[4];
cx q[4],q[1];
u1(3.63486921046794) q[1];
u3(-0.903758838187178,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.41238718859510,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.30777226238014,-1.78537826372912,-0.224928623928874) q[1];
u3(0.847888465281930,-1.44265153237628,-1.09740970135473) q[4];
u3(2.10227724797792,-2.47224348460185,-0.412490256717321) q[2];
u3(1.77093314103391,-3.38003948668646,-0.452853288919921) q[6];
cx q[6],q[2];
u1(1.28021844248041) q[2];
u3(-3.06354871889736,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.51429144251508,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.38620177828758,0.143548241184647,4.03567776902657) q[2];
u3(0.886891988100571,3.54580769665799,-1.27211458394699) q[6];
u3(1.73352728703211,0.218353774609065,1.90503781631130) q[5];
u3(1.29634026405021,-2.78646336900632,-2.70443210027542) q[7];
cx q[7],q[5];
u1(2.86443513827339) q[5];
u3(-2.24381131408499,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.05177012196524,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.44768187327837,3.15505152261914,-1.11997497164120) q[5];
u3(1.15378949418675,-1.10814984523870,4.21803956679914) q[7];
u3(1.02412455729889,-0.518384534491690,1.49571020055713) q[9];
u3(0.128196548689682,-2.31428432643232,0.668907409330324) q[8];
cx q[8],q[9];
u1(1.47208950724616) q[9];
u3(-3.20022864509008,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.83444454192174,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.22309817348658,0.397098236595922,-2.30152262971576) q[9];
u3(1.53934767426383,-1.75523362916576,-0.214675896408037) q[8];
u3(2.36567041503053,0.840642254354599,-2.68715650280804) q[0];
u3(1.94108532264056,-3.10542452313324,2.73317264008642) q[3];
cx q[3],q[0];
u1(2.40972455613136) q[0];
u3(-2.89274924732610,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.03291117089848,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.15515021205970,0.374482002237181,-3.50855562600361) q[0];
u3(1.32018040115682,-3.82237374733236,1.26229706956071) q[3];
u3(0.669052007228243,-2.18058856648136,1.08150626140669) q[6];
u3(0.753117075713274,-0.712725317033751,-0.905306466561869) q[8];
cx q[8],q[6];
u1(-0.414459033363202) q[6];
u3(-1.57608234916245,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.793648239899526,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.63989631858539,2.58492391380140,-0.917096874799483) q[6];
u3(0.669438342667222,2.76843454633020,-2.64640572359555) q[8];
u3(1.73282034865869,1.06192158498283,-2.01562593125313) q[4];
u3(1.37013660310646,1.74865993097296,-4.39908296254786) q[2];
cx q[2],q[4];
u1(0.596940845201932) q[4];
u3(-0.0143058508716902,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.82788705246647,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.458087109402923,-1.19047298494370,0.923881125848746) q[4];
u3(1.21893828534377,-1.21568141744047,4.27223112649264) q[2];
u3(1.56500933753968,3.09786684129055,-1.40418163614686) q[9];
u3(1.73542142281811,1.25517660747720,-2.45832091524341) q[1];
cx q[1],q[9];
u1(-0.0760991036264020) q[9];
u3(-2.15888686426509,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.32654156295713,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.667626800617616,-2.02748067303922,-0.108714169581150) q[9];
u3(1.45904240772054,0.437552174417433,2.78133439034065) q[1];
u3(0.894796386897218,0.405847175826207,1.21324714858971) q[0];
u3(0.912795404742482,-1.89714234477981,-1.28236547449165) q[7];
cx q[7],q[0];
u1(1.61425643294380) q[0];
u3(-2.96120362044093,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.243308786404057,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.33983020263098,2.07480877786460,-1.64809992747050) q[0];
u3(1.15500749285118,3.02690082202737,0.708847517616639) q[7];
u3(2.12095949361406,-1.18044617312674,-0.920747605654129) q[5];
u3(2.05436173504570,-2.18165944932412,0.124815634234608) q[3];
cx q[3],q[5];
u1(1.31486808962249) q[5];
u3(-0.495544516001700,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.0271207127231046,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.212068395347863,-3.69585452098330,1.31343075338994) q[5];
u3(1.29757007568533,4.09137655838677,-2.08679885113237) q[3];
u3(1.96730915935950,0.241821111369714,0.604994825692411) q[5];
u3(1.47445243073028,-2.21445605911059,-1.42604960569722) q[0];
cx q[0],q[5];
u1(1.74335500917080) q[5];
u3(-0.305041797665467,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.72366479109225,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.14513798862815,-3.35101486674163,0.163669476060295) q[5];
u3(1.03116267899091,-4.45926056556767,1.09485566634209) q[0];
u3(1.86739130892009,0.941350759055986,-4.00740735446881) q[1];
u3(2.13144603285510,4.18227853137743,-2.09690411011245) q[9];
cx q[9],q[1];
u1(1.33151569573465) q[1];
u3(-0.862973165743364,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.107906873659432,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.05862455928853,-2.92617668884861,1.04324552438942) q[1];
u3(1.03352885368121,2.14318507784256,0.528111663264106) q[9];
u3(1.50718484448088,1.23035499403906,-1.13595306633076) q[3];
u3(0.865804180295085,-1.32025032534512,-0.00897228716969911) q[7];
cx q[7],q[3];
u1(0.566336823271941) q[3];
u3(-3.60708661887469,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.59839388030036,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.33464362311021,-1.66893839528755,-0.752354750301168) q[3];
u3(2.78506719641032,-2.27830164716072,-2.85385078887344) q[7];
u3(0.873780505252712,1.61426615296174,-2.09765645540675) q[6];
u3(1.24349108995862,-2.81739109119912,2.86415519897265) q[4];
cx q[4],q[6];
u1(0.380082162397504) q[6];
u3(-1.88551976325755,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.87262898919721,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.29037621915325,-0.250800500608989,0.271359360781421) q[6];
u3(2.29418059763007,-1.65595028676565,-2.57520331612581) q[4];
u3(0.478150348154756,1.09851316287770,-0.278805919927509) q[8];
u3(0.304804948269007,-0.0817128278430744,-1.58204267391546) q[2];
cx q[2],q[8];
u1(0.624510569155887) q[8];
u3(-1.35134176324630,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.27833807504392,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.13729909325618,-1.41590198666581,-0.00454119460550440) q[8];
u3(2.43904290408398,-2.05740961871266,-3.58640418974254) q[2];
u3(1.29891787464921,-1.14406657715898,2.72310466835222) q[6];
u3(1.05670619640043,-1.45124858668558,-1.70086538961304) q[7];
cx q[7],q[6];
u1(1.78102810460673) q[6];
u3(-0.741439936192658,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.06154719497951,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.56611322823041,-2.32134819430847,1.58056588665374) q[6];
u3(1.95022567213911,0.427937514256176,3.44488110833160) q[7];
u3(0.933738754457336,0.301231907859031,-0.441318145034645) q[1];
u3(1.37759475426114,-3.92867488367261,1.13051329222920) q[3];
cx q[3],q[1];
u1(1.21991999439615) q[1];
u3(-0.636738347070657,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.214639244025903,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.57442852659216,-1.70777495111631,3.96282661021627) q[1];
u3(0.759293264913763,2.18383628388809,-1.29939492200191) q[3];
u3(0.923498826366145,1.19095756019521,-2.47345761285403) q[2];
u3(1.77213526643644,-3.51385141181908,2.46628220299623) q[8];
cx q[8],q[2];
u1(1.54863602593125) q[2];
u3(-2.62937422380089,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.957827579504121,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.62587229224330,2.74508953474814,0.221835234691049) q[2];
u3(0.599313406706464,0.253746213178018,3.29386219480278) q[8];
u3(1.71282321899003,-0.448452065617147,1.16222556357531) q[5];
u3(1.73939408239291,-2.12084168598328,-0.960489430278166) q[4];
cx q[4],q[5];
u1(1.10713836336152) q[5];
u3(-0.574442139048612,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.441846950226976,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.58338440324230,-0.865152409539679,5.41192068008610) q[5];
u3(1.07712067211610,3.01855521574516,0.0645196013223290) q[4];
u3(2.24665833513531,2.79134923844710,-0.0776778350391856) q[0];
u3(2.82948294149736,0.716236458862383,-5.22310904244834) q[9];
cx q[9],q[0];
u1(2.39774231057073) q[0];
u3(0.0328894106582058,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.38766867049386,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.70262046890383,-1.62349733254682,4.12627517942218) q[0];
u3(0.997281565590237,-0.811200420147022,-0.798645391591486) q[9];
u3(2.25612932276720,4.17737754383800,-2.00024797500107) q[5];
u3(0.435758862133991,-2.47413768843768,3.57297757754774) q[9];
cx q[9],q[5];
u1(3.69108292246477) q[5];
u3(-0.983883780963947,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.60613309225441,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.868814524408780,0.256693190133902,1.79584347083502) q[5];
u3(1.84165705668480,3.43945047361869,1.85442408157233) q[9];
u3(2.90258390089373,2.24864538318539,-0.593907842373599) q[3];
u3(1.88047408858906,3.58965782726996,-0.0276953484934126) q[1];
cx q[1],q[3];
u1(-0.297474863232272) q[3];
u3(-2.13852774881659,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.18565472597373,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.849030813391768,3.87813309328547,-1.02459317364378) q[3];
u3(1.50245924834400,1.76854731751079,2.20503537018412) q[1];
u3(0.611347104283308,0.170189853808376,-1.89315137710371) q[2];
u3(1.37515381864598,0.449413486921376,-5.17495791357012) q[4];
cx q[4],q[2];
u1(0.854540072242738) q[2];
u3(-3.11208455065090,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.38556164335326,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.71128828136184,-1.76224953558074,2.39324208999033) q[2];
u3(1.80254467549436,-2.83282162571396,-1.15310738397084) q[4];
u3(1.49598682755932,-1.66372938058813,1.08095115305018) q[7];
u3(1.23059798853742,-4.26129702053354,-0.221094684110972) q[0];
cx q[0],q[7];
u1(2.34357774935150) q[7];
u3(-1.68487106613688,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.298245724362752,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.46667196206077,-1.78294158919129,-1.49335185152724) q[7];
u3(1.91794532096386,0.498845614324444,-1.65551839900513) q[0];
u3(0.623825469785842,-2.28719307417078,-0.671464586759793) q[6];
u3(2.13953218127543,-3.27537128549801,0.251623953589716) q[8];
cx q[8],q[6];
u1(3.32501187030457) q[6];
u3(-1.33393692624780,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.54974963433073,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.38091526010143,3.51695012891605,-1.96294792417831) q[6];
u3(0.637441530667874,-1.72486586402218,-0.556837582912707) q[8];
u3(1.28713883634258,-2.71605399202837,1.38909278531013) q[9];
u3(1.46630654772394,-3.61692560556974,-0.151195763713601) q[2];
cx q[2],q[9];
u1(-0.731370587399693) q[9];
u3(0.450884016256255,0.0,0.0) q[2];
cx q[9],q[2];
u3(4.25671627866054,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.49040275260555,0.933623012391570,3.13130005616501) q[9];
u3(1.06408432033098,1.07830158818440,-2.95404594658652) q[2];
u3(2.77170146742052,-2.46408785342803,0.684579126294994) q[5];
u3(2.41098631635677,-1.25505532337593,0.921750302233469) q[7];
cx q[7],q[5];
u1(3.72754358645375) q[5];
u3(-1.51786059568983,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.16491908663012,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.62746486567423,-0.638820359136118,2.51053814089868) q[5];
u3(1.32087614799028,-1.46339966352750,-2.86625334300419) q[7];
u3(2.82732567062826,1.39014579572533,-0.533891237442195) q[0];
u3(2.41730153930030,5.06950351936518,0.571729834940181) q[1];
cx q[1],q[0];
u1(1.37608064803862) q[0];
u3(-0.383236598435873,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.68788245259324,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.01396982790303,0.198216993077727,1.25635744777866) q[0];
u3(2.91282620347474,-1.92114861378479,-2.46295357776577) q[1];
u3(2.05498057632458,-3.33138063830656,0.212315080674711) q[4];
u3(2.61110284966981,-2.46833453339956,-0.287010978833669) q[6];
cx q[6],q[4];
u1(0.823404193067276) q[4];
u3(-3.41191605987669,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.53485250527782,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.91429293515110,1.83954331083192,-2.27947349018383) q[4];
u3(0.939494180827163,-0.578993121705469,-5.58551329022829) q[6];
u3(2.47795506922984,3.14153456112417,-1.48054950933711) q[8];
u3(1.68370241845494,2.22008986457564,-2.72463015051089) q[3];
cx q[3],q[8];
u1(0.965789892474930) q[8];
u3(-1.36171417923175,0.0,0.0) q[3];
cx q[8],q[3];
u3(3.13192200466491,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.37034551683359,-1.97833308575103,-0.163946661268346) q[8];
u3(1.23810067229858,-5.43327629464196,0.200616390516263) q[3];
u3(1.54243188862502,3.14629407369989,-0.617332746859748) q[5];
u3(1.28582520491375,0.435387569076425,-1.35202306892921) q[1];
cx q[1],q[5];
u1(-0.214343771171186) q[5];
u3(-1.87854155568266,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.668017222041689,0.0,0.0) q[1];
cx q[1],q[5];
u3(3.10238022025708,-2.73488362702057,1.25342163521505) q[5];
u3(1.89704730826471,2.14997323984710,-4.09275045286669) q[1];
u3(1.65998400313556,-0.0778445328834994,-0.950450444583567) q[3];
u3(0.841764814536150,-4.31092976211119,0.874010898721182) q[6];
cx q[6],q[3];
u1(-1.23360568369690) q[3];
u3(0.635359741250806,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.44131227568597,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.08931863810618,1.09057924976119,-1.65662814873325) q[3];
u3(1.19754120317016,1.55386628584032,-2.08632234291246) q[6];
u3(2.19948470116225,0.461615686339753,-1.26686452045815) q[0];
u3(1.09959734112539,0.698459865574929,-3.99321395666800) q[9];
cx q[9],q[0];
u1(3.05726665431752) q[0];
u3(-2.05756718502182,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.329295010798392,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.68114533331678,3.67877947108427,-0.338633339489457) q[0];
u3(1.55672121400822,-1.66862754000327,-2.11745507924226) q[9];
u3(2.10360176972547,0.0115525660553163,-1.31919795880768) q[2];
u3(1.68337627574361,1.27658875958147,-4.93621767654218) q[4];
cx q[4],q[2];
u1(2.45191058176926) q[2];
u3(-2.11571207647454,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.18865975152130,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.807347098866026,-3.29477210165428,1.55684803179688) q[2];
u3(1.17429290496377,-1.48793625693866,-1.81125466859307) q[4];
u3(0.509775099491262,1.38089575291128,0.833940356979375) q[8];
u3(1.41676775503319,-0.665274184228284,-3.57071845288271) q[7];
cx q[7],q[8];
u1(1.40606123813004) q[8];
u3(-0.335818756672602,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.67612655466160,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.91442663050098,1.58224645691136,-1.58873966030586) q[8];
u3(1.62871540702321,0.561175725700455,2.38216579211036) q[7];
u3(1.07705978160853,-0.561826877653975,0.890473584587684) q[4];
u3(0.752747404100579,-1.05168590260050,-1.97251163377556) q[7];
cx q[7],q[4];
u1(0.951121779902163) q[4];
u3(-0.0659550612155371,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.62464657349747,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.89847386584523,1.53174914320517,-4.21594484640895) q[4];
u3(1.17662153070557,-1.44223172387857,-2.24968701538278) q[7];
u3(2.59224441687107,1.71040850502987,-3.55784650538110) q[9];
u3(1.64446472346643,-2.42617298341507,3.57785412476934) q[5];
cx q[5],q[9];
u1(0.728157857996363) q[9];
u3(-3.22696771627499,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.75395860253337,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.95001036700717,1.51570998868677,-4.31291140817905) q[9];
u3(1.21442846145038,2.05154177932774,3.00933051792766) q[5];
u3(2.19819414488946,-4.09779579409033,2.16268396751076) q[0];
u3(1.25531280041600,-2.21365802825814,2.75812891410496) q[6];
cx q[6],q[0];
u1(-0.196503350435816) q[0];
u3(-1.91418486380148,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.721642399385624,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.933541582495709,-1.19635064809994,1.35243516792896) q[0];
u3(1.33746099798780,0.130445251364581,-3.21315378375712) q[6];
u3(0.660869756842737,-1.62413211886580,3.99404101346098) q[8];
u3(1.86409955122064,-1.04939715570114,0.312687439816918) q[1];
cx q[1],q[8];
u1(3.57089100424743) q[8];
u3(-1.15931134492683,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.21879689589706,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.81698711874010,1.93567534571537,-1.68809422972338) q[8];
u3(2.25699673938400,2.63942999493233,0.384680644654861) q[1];
u3(1.88591488162952,-0.286669603747993,0.804766617255170) q[3];
u3(2.34740098411784,-0.531790913825640,-1.37841056566215) q[2];
cx q[2],q[3];
u1(0.266763531578270) q[3];
u3(-0.693843227302164,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.63006989407795,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.50102403571212,2.78229207999250,-3.45434130224815) q[3];
u3(1.52775547481005,2.26280735162165,3.71299719736245) q[2];
u3(1.04802650966524,-0.121557165462706,2.10957464161655) q[7];
u3(1.09751419363248,-1.32534457274260,-2.68668814819329) q[3];
cx q[3],q[7];
u1(-0.149613504053551) q[7];
u3(-2.39477526587433,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.37297840969539,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.82513328562437,0.482962858523410,1.90993303264579) q[7];
u3(1.15633517624714,-2.63724833925497,-2.42788597038504) q[3];
u3(1.98165112858769,0.195252090818770,2.81435578931068) q[6];
u3(2.95847079906265,-2.10523671225761,-1.22658263495269) q[9];
cx q[9],q[6];
u1(2.00261805733196) q[6];
u3(-2.93416813943997,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.623525308661869,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.90010742265914,-2.46527405636189,3.29874739639176) q[6];
u3(3.06665356344139,-1.64370183520761,-4.50847307540598) q[9];
u3(1.27802584629605,-0.124751442107170,-1.61236653927358) q[2];
u3(2.12192146170668,-3.41087889269431,2.12374689300619) q[1];
cx q[1],q[2];
u1(3.20689454172939) q[2];
u3(-4.61653698280065,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.0358147376968774,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.76169824227959,1.54824423495572,1.25902797619009) q[2];
u3(1.38944922340566,-4.67172041065945,-0.0422423576678215) q[1];
u3(1.33716096078922,0.254115415288290,1.37620139602607) q[4];
u3(2.44055340266632,-2.86954613025853,-1.21286083188877) q[5];
cx q[5],q[4];
u1(1.44145431592006) q[4];
u3(-2.98034328730200,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.662874251750395,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.456265639832163,3.94411367424350,-1.06497034144472) q[4];
u3(1.57476097265021,-1.96627877142153,-3.77923163209245) q[5];
u3(1.58985599613696,-1.22366149270552,-0.866142526364075) q[0];
u3(0.215236503021797,-1.54333699063736,-3.36633526613606) q[8];
cx q[8],q[0];
u1(-0.110728380076351) q[0];
u3(-1.30879439887964,0.0,0.0) q[8];
cx q[0],q[8];
u3(2.61939086410235,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.92849581767378,0.803727165903367,-3.62091991697673) q[0];
u3(1.19565117687014,-1.83876352019333,-0.619511136475625) q[8];
u3(2.36550935988652,0.421980020276284,-3.27976420716114) q[8];
u3(2.69987601409758,3.34770201192771,-2.11795378641153) q[1];
cx q[1],q[8];
u1(1.69158213064106) q[8];
u3(-2.56899236416497,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.19981608051225,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.78616301209386,-0.278206092883397,0.289311614760504) q[8];
u3(1.32555989655770,-5.36240294381758,0.190043908210697) q[1];
u3(2.20584952703576,1.87998076413787,0.104012740604867) q[5];
u3(2.23012604287688,-0.650649957106459,-4.98224051605825) q[3];
cx q[3],q[5];
u1(0.0402397440029649) q[5];
u3(-1.55940193872674,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.370052139574887,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.902330767320422,-2.94027129304999,2.21709453289521) q[5];
u3(1.46005611290884,4.13625294240241,1.82809432036293) q[3];
u3(1.88896970695267,-0.701563916063437,2.62329192265710) q[7];
u3(2.32238315532124,-2.50832533902735,-2.04481300251003) q[9];
cx q[9],q[7];
u1(1.24468192662846) q[7];
u3(-3.10682463266624,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.07091727669771,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.388185403112741,-1.46201292159324,4.01696265824973) q[7];
u3(1.22112507882410,-4.77551292021988,-1.25054811337021) q[9];
u3(0.346153966519757,-2.58068893506787,2.40721539486966) q[0];
u3(0.448593633692844,0.985965877626175,-1.96162839119768) q[4];
cx q[4],q[0];
u1(1.61740163931121) q[0];
u3(-3.01924874038054,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.68389728628788,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.71371487323954,-1.55053221106589,1.97558184939221) q[0];
u3(2.17298542947668,0.801620011452214,3.45162854747632) q[4];
u3(1.28404334808289,1.31695892139943,-2.84632320396038) q[2];
u3(0.663529999292461,-3.02964304466433,2.20258483695208) q[6];
cx q[6],q[2];
u1(1.56689581507130) q[2];
u3(-0.626515065422745,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.328404366566216,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.17718934853146,1.22244840943037,2.31978205829301) q[2];
u3(1.09382804483480,2.11727170727274,-1.65533054929949) q[6];
u3(0.213242164342139,1.19809090009790,-0.976682553126418) q[7];
u3(1.32094011459224,0.628365336940617,-1.59169074814055) q[8];
cx q[8],q[7];
u1(3.17874699987536) q[7];
u3(-2.11053573718452,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.486998754387371,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.27649629438012,2.60887268983217,-2.51052090583594) q[7];
u3(2.42738569386321,0.695762009182276,1.46051622335724) q[8];
u3(1.79049145086265,3.73419246798898,-1.85634855628118) q[1];
u3(2.10422957941948,2.22203746528435,-2.02255007608584) q[9];
cx q[9],q[1];
u1(0.355899734133765) q[1];
u3(-1.21807649278025,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.41267227885485,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.16258308931436,3.18954393279592,0.303651230828256) q[1];
u3(2.79475803432118,4.58183525184171,-0.172434182361233) q[9];
u3(1.58519378677737,1.35881411484190,-2.91376125825733) q[5];
u3(1.29663062981920,-2.61124666122876,2.30861629290092) q[6];
cx q[6],q[5];
u1(3.42878846990345) q[5];
u3(-1.97895972799752,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.50969116887405,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.62149844758149,1.23928753257900,1.08068203624231) q[5];
u3(1.60358368737440,-1.23524034958656,-3.43854329472732) q[6];
u3(1.86828375522682,-0.500901836359650,0.489931215848770) q[4];
u3(2.85599846709369,-0.601366768598276,-1.26422163407627) q[3];
cx q[3],q[4];
u1(3.11588226636346) q[4];
u3(-0.654969015474862,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.68305625294749,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.31242781265608,-1.73985518514939,-0.223254716567172) q[4];
u3(2.16549654112811,-2.18340448437066,3.31868972660984) q[3];
u3(2.36256836151276,-3.47381497927968,0.442526334051754) q[2];
u3(2.45354322544247,-0.774424169909904,1.24189227388267) q[0];
cx q[0],q[2];
u1(-0.0700531906360622) q[2];
u3(-1.93337176158850,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.647548392604340,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.23010462836625,-1.07350687254655,2.56725663620793) q[2];
u3(2.46221694295761,-0.125555062766713,4.30068210398749) q[0];
u3(0.369746129508389,1.07510255001601,-2.23551546101039) q[1];
u3(1.63346388888289,-2.30167292287368,3.43931521399985) q[6];
cx q[6],q[1];
u1(0.00992550358911326) q[1];
u3(-1.68577933708152,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.763896314476462,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.618418137514040,-2.66729757057215,3.57999892273803) q[1];
u3(2.43746863346330,0.856029139934881,1.81901300115631) q[6];
u3(1.12689835968421,-0.921627275228099,-1.38967961538558) q[3];
u3(2.18042560237610,1.21215016269843,-4.79376618029688) q[4];
cx q[4],q[3];
u1(0.852970782977051) q[3];
u3(-0.606929534248428,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.58378121543651,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.829017067612302,-0.217560360321009,-0.617054535175099) q[3];
u3(2.62466385301618,3.14968728651127,-2.75673007142447) q[4];
u3(1.84859053038854,1.37772078059100,-0.427419783948053) q[2];
u3(0.426461678830188,-0.00941225874400664,-2.85938776558338) q[5];
cx q[5],q[2];
u1(3.04688071842311) q[2];
u3(-0.635902052259279,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.05248597036255,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.12690303633765,-0.234838739715597,0.102952722765977) q[2];
u3(1.12379158052818,-1.04514721213066,2.83475666660546) q[5];
u3(0.871386062803675,2.12082310571651,-3.05520778462142) q[9];
u3(0.833262471876548,-2.62919081481130,3.09922761341923) q[0];
cx q[0],q[9];
u1(2.90385099660075) q[9];
u3(-1.48046797858188,0.0,0.0) q[0];
cx q[9],q[0];
u3(0.0433275500601209,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.28330291843566,0.0764667652992939,2.87115182543652) q[9];
u3(1.15778988358419,0.664487258649661,2.75428145343216) q[0];
u3(2.38714592426454,1.70892558480696,0.275526937477485) q[7];
u3(0.456075161089902,0.194569175972114,-4.87075184725161) q[8];
cx q[8],q[7];
u1(0.351131209378642) q[7];
u3(-1.43213433506609,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.30192126514758,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.75105054543336,0.852063457656736,-2.88756143807436) q[7];
u3(1.55488226026760,-2.54107667756717,3.20740109991532) q[8];
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
