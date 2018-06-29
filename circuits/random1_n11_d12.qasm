OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.90345188986591,0.428574706858142,-0.522099302701443) q[1];
u3(1.19836140855463,0.532599130097258,-4.76180050781862) q[3];
cx q[3],q[1];
u1(1.28480214270039) q[1];
u3(-0.859519925226983,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.98683047917972,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.928130403386407,-3.53110136572857,0.818182700655595) q[1];
u3(2.54499017381182,2.06756796796004,2.40389209122688) q[3];
u3(0.977180210812892,-0.614950850375275,0.993261107403019) q[8];
u3(0.356886569185359,-1.43192618008566,0.699521602950532) q[5];
cx q[5],q[8];
u1(-1.16005986996167) q[8];
u3(0.406171711848976,0.0,0.0) q[5];
cx q[8],q[5];
u3(3.89217720310556,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.938488443716941,-2.10516897929217,0.257053677014570) q[8];
u3(1.37844957536144,-1.52168940300607,2.43618372416147) q[5];
u3(1.61672355509924,-1.98211818076360,-0.772604166189589) q[2];
u3(2.28518270821465,-4.73855293040909,-1.15223725334015) q[7];
cx q[7],q[2];
u1(1.51894307341274) q[2];
u3(-3.38059848648567,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.35614821017704,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.236281875597334,2.39516595463639,-0.668494470338549) q[2];
u3(1.46080266862750,-3.28010703359393,-0.145288859223985) q[7];
u3(0.608988109100187,-2.68480670314625,3.22963750999312) q[0];
u3(1.49931320345937,0.524440998377532,-2.15435359988051) q[6];
cx q[6],q[0];
u1(-1.21491363019845) q[0];
u3(0.691031599891732,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.97123461800832,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.98127925217480,1.91349509981737,2.14895850413985) q[0];
u3(2.52502552986731,1.55353008175987,-0.911280773434032) q[6];
u3(2.78852832870213,-1.57095256873746,0.404548243745061) q[10];
u3(2.59981741940040,-1.21398263812594,0.136098122967972) q[9];
cx q[9],q[10];
u1(0.958554864456212) q[10];
u3(-1.34541443900320,0.0,0.0) q[9];
cx q[10],q[9];
u3(-0.202645722281383,0.0,0.0) q[9];
cx q[9],q[10];
u3(0.997402383146640,-1.38390266904538,4.09127511553852) q[10];
u3(1.53901674427272,-0.431334519170298,-2.68140699795540) q[9];
u3(0.0814577225802480,2.69344418793854,-2.55091217907404) q[8];
u3(1.23580859772266,-2.63159508722520,2.06142769310042) q[3];
cx q[3],q[8];
u1(1.69907405718728) q[8];
u3(-3.11866520636322,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.21066855314308,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.612126479934004,1.80162116964633,-1.05922651077343) q[8];
u3(1.67133750643926,2.21190858305673,-3.11360861939202) q[3];
u3(1.65479866515431,2.34780433980538,0.132230787253828) q[7];
u3(2.84306564205766,-0.854717714428362,-5.05962976716674) q[9];
cx q[9],q[7];
u1(1.13133325979242) q[7];
u3(-3.28686856054507,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.31249660926157,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.00220454738228,-2.30116686051222,3.25383031077916) q[7];
u3(1.81311212287388,1.30043101964533,-0.520533021982012) q[9];
u3(1.54403444426460,3.20685338277088,-0.332338063628746) q[0];
u3(1.74806106777833,1.86214645505480,-1.54348004706854) q[5];
cx q[5],q[0];
u1(0.242624689521095) q[0];
u3(-0.477766020577349,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.05458510236050,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.71047475198909,0.378235556625848,2.65133963999991) q[0];
u3(1.98078070395523,0.377721991489508,3.06777999931050) q[5];
u3(1.44288196377758,-0.157265691755446,1.90038475202136) q[6];
u3(1.25733932599656,-0.508355199713498,-1.41648865030904) q[10];
cx q[10],q[6];
u1(-0.730266386074071) q[6];
u3(0.0810326206096730,0.0,0.0) q[10];
cx q[6],q[10];
u3(4.05459847102014,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.28993151141093,2.54833614678710,0.102942647102516) q[6];
u3(1.06695291361249,2.50025649476478,3.51070720688491) q[10];
u3(1.07322388582095,1.14957103123837,1.45997356834555) q[1];
u3(1.30594093634825,-0.236225221221142,-3.06983193588314) q[2];
cx q[2],q[1];
u1(3.51902306153001) q[1];
u3(-1.33058919559879,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.39816701021357,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.02501777982075,-0.392287602135345,-0.784663980800042) q[1];
u3(0.754245506647306,-0.184853924064115,1.53553050178993) q[2];
u3(1.26624007354148,1.64137310999587,-2.71925112743594) q[8];
u3(0.784118280443185,2.38155911114396,-3.57258889812914) q[1];
cx q[1],q[8];
u1(2.93090027695969) q[8];
u3(-1.93015620818667,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.781984643357581,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.832735062308443,-1.48321553461259,-0.656136949756806) q[8];
u3(0.765422004146485,3.25852554868291,2.06658239250053) q[1];
u3(2.20865202044390,-0.163648985215463,0.414553442856052) q[4];
u3(1.34950304329320,-2.79624988978696,-1.81968778737587) q[10];
cx q[10],q[4];
u1(1.41891903195044) q[4];
u3(-0.0891003750694748,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.422555137756962,0.0,0.0) q[10];
cx q[10],q[4];
u3(3.00807513549329,0.305243469381570,-1.82350193734195) q[4];
u3(1.48612299935271,-0.647923739844336,2.33289874714802) q[10];
u3(0.830576502353115,-1.81689605843310,2.55390973160223) q[0];
u3(1.17154539014598,1.07760102483669,-1.41766695799708) q[2];
cx q[2],q[0];
u1(3.98761741163589) q[0];
u3(-4.48586730316116,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.636209114268528,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.429736808329054,0.177348918120199,1.34329342499483) q[0];
u3(1.62895588559223,3.57852147208831,1.18785252545406) q[2];
u3(1.91086955695510,-0.649729017876231,0.618946918552353) q[5];
u3(2.12443449327603,-3.56516713858656,-0.417982178427703) q[9];
cx q[9],q[5];
u1(0.838278707870707) q[5];
u3(-1.41449881193253,0.0,0.0) q[9];
cx q[5],q[9];
u3(3.03538115242373,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.25860914354505,2.95666059111625,-1.05895662849257) q[5];
u3(0.576028982011110,-3.71707813212748,1.43628938026855) q[9];
u3(1.90682735186801,-1.22800987830442,-0.647883191939075) q[3];
u3(1.15271926025147,-2.25351302908101,0.816036388255270) q[7];
cx q[7],q[3];
u1(2.31587937760918) q[3];
u3(-3.20130563247941,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.26378487133447,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.24891713648716,2.52128307176380,-2.22443778146124) q[3];
u3(1.69965347040641,0.559106161138300,-2.98477612209049) q[7];
u3(1.87534405759224,-0.949342782759122,1.52293807343594) q[4];
u3(2.38141445144008,-2.29210060862382,0.0502748222710017) q[7];
cx q[7],q[4];
u1(1.50947606700797) q[4];
u3(-0.334156179126960,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.04523174441248,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.09471235455094,0.499039676296700,2.17235996122804) q[4];
u3(1.45638045707672,3.76739770592284,1.35604788294996) q[7];
u3(2.41484675115696,0.917728985659150,-1.93627310667988) q[5];
u3(1.81701529734885,1.62475954928218,-4.20590961740691) q[10];
cx q[10],q[5];
u1(2.46831263320462) q[5];
u3(-2.84120695015338,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.63631289913832,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.22224034024052,-2.55063755345874,3.16245016802667) q[5];
u3(2.58137679179780,-3.85905192323851,0.720199855780270) q[10];
u3(2.40233168423042,3.25841569694390,-0.334369913140270) q[2];
u3(2.68454015007957,5.93064033457795,0.162601755168791) q[0];
cx q[0],q[2];
u1(2.93870108253630) q[2];
u3(-1.80681589327859,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.365671134563420,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.04256491448247,-2.10769464970984,3.78799896894584) q[2];
u3(2.41019582531281,-1.77257231007640,1.28364440449165) q[0];
u3(1.93195965459574,-0.333971890548116,1.08457271178263) q[9];
u3(1.78119234001846,-1.83043904928463,-1.70857640246793) q[1];
cx q[1],q[9];
u1(0.145199399270285) q[9];
u3(-1.83967550378556,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.564316785000333,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.37290242797515,3.67384258827462,-2.30401279432356) q[9];
u3(2.18227925620440,-0.630144005224589,3.63449587022389) q[1];
u3(1.19477437620302,0.922486401217750,-2.09218863732019) q[3];
u3(1.16772295520365,1.62429528969163,-4.13291198490853) q[6];
cx q[6],q[3];
u1(0.879767694389125) q[3];
u3(-1.69814728815629,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.64560063656335,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.157978026968046,-1.08574432580965,0.812444402780078) q[3];
u3(0.633583609604983,-3.05094899367866,-3.20067073504790) q[6];
u3(1.07679717474581,4.05280239714406,-1.08946825234205) q[10];
u3(1.94903027107142,2.64385666351167,0.449724691003292) q[8];
cx q[8],q[10];
u1(2.78304484202875) q[10];
u3(-1.76633490895655,0.0,0.0) q[8];
cx q[10],q[8];
u3(1.01592775144820,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.34877927819844,0.949204332694668,0.820322526115213) q[10];
u3(1.58207523077890,2.96885681912508,0.868452212485397) q[8];
u3(0.472824792888713,0.527542870951569,-1.04464966336962) q[1];
u3(1.25458288654194,-3.95677119197020,1.51435519355677) q[2];
cx q[2],q[1];
u1(-0.0539683799999400) q[1];
u3(-2.49967710228221,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.38329142810243,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.53242702550762,3.79587063282919,-2.06471607352043) q[1];
u3(0.815882535783523,-2.54799797235257,0.200126268835266) q[2];
u3(2.32792436887231,3.18028614542337,-0.0619272942596094) q[0];
u3(2.88765329755335,2.88860141116195,-1.21004641148652) q[6];
cx q[6],q[0];
u1(1.13139850281985) q[0];
u3(-0.949140833264461,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.73234825898767,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.46331980412303,-2.34811400577232,1.99887460193971) q[0];
u3(2.44266457165488,0.435576226969821,-1.50524650853726) q[6];
u3(1.24073789678727,0.497169794963684,-0.935715738233244) q[3];
u3(1.52523992107392,-3.60632832140912,0.559066550951274) q[7];
cx q[7],q[3];
u1(0.799801207024082) q[3];
u3(0.193665564297578,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.89837966565547,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.60114275315207,1.61881904064328,0.531166509892437) q[3];
u3(0.731139679180259,-1.30010865158430,3.39632314638975) q[7];
u3(2.20299567335600,-3.27632469958934,0.238740092762811) q[9];
u3(2.93054466723469,0.702473875316705,2.97102273795972) q[5];
cx q[5],q[9];
u1(-0.371514902775909) q[9];
u3(1.27263925182217,0.0,0.0) q[5];
cx q[9],q[5];
u3(3.80649696625764,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.30001787041168,1.46339005803293,-2.16889402490889) q[9];
u3(2.29052817641411,5.00016339653812,-0.846119915948425) q[5];
u3(0.980832557632480,0.188996694292101,-2.39075444406012) q[1];
u3(1.13334632858220,-4.97855460959715,1.13160049800798) q[7];
cx q[7],q[1];
u1(2.72025274393287) q[1];
u3(-1.52641413317161,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.345116491105761,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.876885531090642,-2.60595593849264,-0.951881781864628) q[1];
u3(1.23041472767762,-1.81219574170128,-3.59842343947440) q[7];
u3(1.36352293517256,-1.29330868874112,0.232050368320597) q[2];
u3(2.14712705988853,-2.17277876070363,-0.453312698563033) q[3];
cx q[3],q[2];
u1(3.62855639863035) q[2];
u3(-1.53098292997758,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.50126121906197,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.82518472511416,-0.779717700224420,1.03838784935939) q[2];
u3(0.750621130133257,5.15560746667260,-0.822795216540293) q[3];
u3(1.02140237840721,0.793964752686115,-3.74675334065379) q[10];
u3(0.512412486091580,3.06528608593835,-2.91065140654543) q[9];
cx q[9],q[10];
u1(-0.0990999521157954) q[10];
u3(-2.52632022775987,0.0,0.0) q[9];
cx q[10],q[9];
u3(1.22100531104510,0.0,0.0) q[9];
cx q[9],q[10];
u3(0.532907416593690,0.415950601931546,-0.356066835069134) q[10];
u3(1.53708836643812,1.46976662650215,4.79418480464329) q[9];
u3(0.225960860956602,0.512257841996437,0.500031495800801) q[0];
u3(0.841623427740454,0.444183392375245,-2.31365091715907) q[4];
cx q[4],q[0];
u1(3.26589416415258) q[0];
u3(-0.949911918961332,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.62805190640359,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.38040812449001,-1.48353200696051,2.23407070983199) q[0];
u3(1.94429313466198,0.805615035308426,5.43901911283318) q[4];
u3(0.270205369041208,2.66109344852060,-1.69298530441580) q[5];
u3(1.13047159089722,-2.74817536046930,1.18339518224594) q[8];
cx q[8],q[5];
u1(1.18961576579915) q[5];
u3(-3.17403999370444,0.0,0.0) q[8];
cx q[5],q[8];
u3(2.11074834098076,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.36749339025025,-0.290795141882836,1.33192649137743) q[5];
u3(2.75111572527636,-2.52766872869611,-1.67741483889787) q[8];
u3(0.642311522341188,-1.61383456248495,3.16106852996806) q[1];
u3(0.923028951023446,-1.67943639924556,0.0904695619150456) q[0];
cx q[0],q[1];
u1(-0.534913256971296) q[1];
u3(-1.59386672215628,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.811208414624389,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.65836283135961,0.384900497891897,-3.97500228149080) q[1];
u3(0.926303039641664,0.112310576340951,-2.67901909237137) q[0];
u3(0.183690049555366,1.67807839443351,-0.567800695474283) q[2];
u3(0.839497871448051,-2.01809836235499,0.584626076279587) q[3];
cx q[3],q[2];
u1(1.68661538110723) q[2];
u3(-2.63673742731741,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.999638679842181,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.67661973832486,-1.56060996206981,3.52568768023424) q[2];
u3(2.53924179577861,-4.36424657242050,-0.424114046352945) q[3];
u3(1.08427144712210,-2.24609647549968,1.98291235517786) q[5];
u3(0.508262563327918,-1.76623104169817,-0.508405387573651) q[4];
cx q[4],q[5];
u1(-0.583004973632298) q[5];
u3(0.131731990976438,0.0,0.0) q[4];
cx q[5],q[4];
u3(4.25861323785328,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.69523621969558,4.15666500596850,-1.60612281691992) q[5];
u3(2.01050913222716,-3.18572035040852,-2.80242322750561) q[4];
u3(0.330475573549047,-3.18680723308506,2.46090331211400) q[6];
u3(0.824389995746014,0.418246228789780,-1.70091936085013) q[7];
cx q[7],q[6];
u1(0.941941831218274) q[6];
u3(-3.08658524102165,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.78644777567487,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.01853397785884,2.72178900361137,0.822794850224514) q[6];
u3(1.50451632356646,1.84851063227156,-0.0569560610160996) q[7];
u3(1.64881806164502,0.874814830008406,1.36992855864296) q[8];
u3(1.47448095265816,-0.838501160323111,-2.78051844147230) q[9];
cx q[9],q[8];
u1(2.91285633924484) q[8];
u3(-1.77268391716825,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.879698893649820,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.70675437618289,0.678042091358858,-5.00663633083459) q[8];
u3(1.35260352055149,-4.29145966887915,-0.0618072471378412) q[9];
u3(2.90904501940174,0.861198282159179,-2.26533558943868) q[0];
u3(2.31215824230458,2.76851610623977,-0.533123852324555) q[9];
cx q[9],q[0];
u1(0.509778111609455) q[0];
u3(-1.43885814562389,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.41241237976585,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.07378674642668,1.31821135501060,-3.53710331246080) q[0];
u3(2.26792340031502,1.88665223135777,-1.68841801056342) q[9];
u3(2.14487895639564,0.634800908416066,-1.14881102786073) q[7];
u3(1.15999631308074,-0.0568982086037277,-3.71369276359856) q[1];
cx q[1],q[7];
u1(3.35266039300478) q[7];
u3(-1.30103204091540,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.49745994008107,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.68796897624827,4.03575781658581,-0.282646210207999) q[7];
u3(1.50178861666619,-2.68436202518113,-3.59643571940306) q[1];
u3(0.629498276207459,-1.69240735920202,0.967700975872940) q[4];
u3(0.451237404025538,-3.22766745220106,1.95792274764320) q[6];
cx q[6],q[4];
u1(-0.0235390080751452) q[4];
u3(-2.08573079568032,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.693696208207426,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.35549714314941,-3.64251486644374,1.74067484122853) q[4];
u3(2.45313330807025,-1.01469820566946,1.61114963547679) q[6];
u3(1.61296995092068,0.544963027691611,-3.13803312862595) q[2];
u3(2.86518416570197,-2.51266156881810,3.44165743231944) q[10];
cx q[10],q[2];
u1(2.80746348925225) q[2];
u3(-1.70731624357348,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.25155530591375,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.18782705667398,-0.795658267262293,-0.766874912025186) q[2];
u3(1.79794208833901,1.03270045457938,3.43584905395324) q[10];
u3(1.88874094308582,0.697045306658538,-3.44821325466284) q[3];
u3(2.26011263320318,3.46310830672826,-2.42900676693245) q[5];
cx q[5],q[3];
u1(0.997705555110904) q[3];
u3(-3.55078685871710,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.56776738631548,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.38703740088341,0.757023817575325,-1.79825182947310) q[3];
u3(2.35184401010650,0.543796591421191,1.88559227801905) q[5];
u3(1.76199280030201,1.42235936493733,0.226826412571440) q[7];
u3(1.65009082118531,0.486069133123150,-3.84350537104278) q[3];
cx q[3],q[7];
u1(1.88639697942900) q[7];
u3(-2.66668970465178,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.659265891781846,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.48042305295352,-0.624826503029575,3.07779489439143) q[7];
u3(2.78459825857496,4.18911695885185,1.11447424730316) q[3];
u3(1.43852791362765,-1.27077709034141,-0.701976693594219) q[0];
u3(1.80992986621501,-2.11998409202855,-0.127945790536271) q[6];
cx q[6],q[0];
u1(1.40032912118468) q[0];
u3(0.00754466216107486,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.74282208126459,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.84246704989492,1.19379504049493,-1.59055582721580) q[0];
u3(2.33595482566835,0.210492979331608,-0.459091945759861) q[6];
u3(0.768033512390764,-0.921888621021826,0.913615528380354) q[2];
u3(0.675553526891567,-1.42164119325049,-0.455023890590040) q[5];
cx q[5],q[2];
u1(0.820096850626722) q[2];
u3(-1.48366296378839,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.411806162166830,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.62537784691265,-0.942465710301333,1.48073602705190) q[2];
u3(1.88281667953675,4.15522565008176,-1.30610709589337) q[5];
u3(1.86030253715914,3.49145610635830,-0.618366246819229) q[8];
u3(2.23805900115412,2.36479488062321,-1.04011893079090) q[1];
cx q[1],q[8];
u1(0.530429700479384) q[8];
u3(-0.615546626239390,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.92732650127529,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.97541104454611,0.400105870660328,-2.10704840110197) q[8];
u3(2.45117321643105,-0.507133520082857,-0.361698146263095) q[1];
u3(0.943616998771953,-2.31194164786574,3.23813890149530) q[10];
u3(1.40979619404393,1.86964777587402,-1.37199263116749) q[9];
cx q[9],q[10];
u1(2.04314090427726) q[10];
u3(-2.63170076515759,0.0,0.0) q[9];
cx q[10],q[9];
u3(-0.109859610578553,0.0,0.0) q[9];
cx q[9],q[10];
u3(2.20990446362519,-4.72257511970570,1.02572730861245) q[10];
u3(1.33087408635111,2.34619928965533,-3.74105756664371) q[9];
u3(0.815564269234549,-0.689309237733316,1.08251713438115) q[4];
u3(0.655116995117205,-1.47456673002455,-0.234399440481178) q[2];
cx q[2],q[4];
u1(2.26269886704995) q[4];
u3(-0.215083512527336,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.75953869515059,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.43981965643116,-0.617024030292790,-0.500495976001192) q[4];
u3(2.33921364971129,4.50322506641739,0.175833648242579) q[2];
u3(1.17022925088160,1.18607870280780,-2.56797057917251) q[10];
u3(1.81856103380047,-2.01362404883028,2.61085589469600) q[9];
cx q[9],q[10];
u1(3.20653340397029) q[10];
u3(-1.39821207955037,0.0,0.0) q[9];
cx q[10],q[9];
u3(2.19023907719658,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.27969794375373,1.29966156451118,1.70450126116323) q[10];
u3(1.17609133291405,-0.267479436139471,4.91618253441257) q[9];
u3(2.34426602988570,1.48830913289989,-0.419927577338953) q[0];
u3(2.29916074767325,0.812665903613308,-3.57568199625063) q[6];
cx q[6],q[0];
u1(1.70545849028620) q[0];
u3(-2.85035955436941,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.912841992804032,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.25700961409825,-0.0952500707947328,-2.51285003532533) q[0];
u3(1.43282498724568,-3.96585182398111,-0.427338308334331) q[6];
u3(1.56191657275004,0.846298087393648,1.77930026309183) q[1];
u3(2.23005577414564,-1.28999770113451,-0.877760509829990) q[8];
cx q[8],q[1];
u1(1.18350966041495) q[1];
u3(-0.0517570090781818,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.68119798722194,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.18211728345213,-0.941344123823261,-0.724761166009153) q[1];
u3(1.21555888744818,-1.50944829372531,-2.67382947259332) q[8];
u3(2.67149401119200,2.37148580422918,-1.21847006175825) q[5];
u3(2.24711256656682,2.21271656797594,-3.17227259857438) q[7];
cx q[7],q[5];
u1(-0.209653826299637) q[5];
u3(-2.40473038932351,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.47883911247603,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.832890594973465,3.54933074343906,-0.548983303922240) q[5];
u3(1.07766763084087,2.87340419173251,-0.675146645302738) q[7];
u3(0.647562164747550,-0.157961895443522,0.196450380294772) q[9];
u3(0.500621403767395,-1.50940468039874,-1.40005348964759) q[6];
cx q[6],q[9];
u1(2.12291217643930) q[9];
u3(-2.57687537355837,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.395683392589762,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.24867372737032,0.964671306740774,1.91698004168631) q[9];
u3(0.116960470699241,-1.48545193772190,-4.55059546451309) q[6];
u3(0.988884486441455,-0.937238698038066,-0.171945024593246) q[4];
u3(0.879422628287703,-2.89851708493956,1.33421410496173) q[10];
cx q[10],q[4];
u1(-0.432666985633683) q[4];
u3(-1.77389282804868,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.725810767421298,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.50333947300104,2.72230334092607,-1.21564490394939) q[4];
u3(2.22482975235118,-3.49905598249779,2.75560335126131) q[10];
u3(2.26997407055460,0.753486127762981,1.00936302234224) q[3];
u3(1.33218271073685,-1.43455325443715,-2.23054667060709) q[2];
cx q[2],q[3];
u1(1.55282803162746) q[3];
u3(-1.01734585527343,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.548026443053780,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.40610774276275,-0.0779792656947736,-1.03804974541323) q[3];
u3(0.671187887189328,-1.88505994446066,-3.40983787687170) q[2];
u3(1.11048527751456,-0.841058372623837,1.29067667761596) q[7];
u3(1.73831973813286,-1.47025617063543,-2.48610636301619) q[1];
cx q[1],q[7];
u1(-0.320479072094970) q[7];
u3(-1.67437794268278,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.586080208336306,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.899282900260052,0.459550757549456,-1.87228656283643) q[7];
u3(0.355365803628818,-0.554480572218702,-4.21333287889289) q[1];
u3(2.49477779632392,3.07644941368507,-3.18698984571415) q[8];
u3(0.986035508245504,0.524062082845129,2.20195798068453) q[0];
cx q[0],q[8];
u1(1.51025248694170) q[8];
u3(-0.929650736089741,0.0,0.0) q[0];
cx q[8],q[0];
u3(-0.308201851589708,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.90432468168979,1.81527251027910,-2.51425201440555) q[8];
u3(1.21916375418603,1.82746534678953,-4.32249740903867) q[0];
u3(1.04655874697557,1.13677186670011,0.965590316433788) q[7];
u3(0.992605236876624,-0.00267220232689613,-3.36064901433667) q[10];
cx q[10],q[7];
u1(0.812192949707702) q[7];
u3(-0.367000562413228,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.13256295925168,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.926396762179500,-1.24552912920399,2.11055802581062) q[7];
u3(2.17065908445985,-3.81137438724119,-0.179092664761687) q[10];
u3(1.68835555730419,-1.85364694297427,0.495931043341462) q[1];
u3(2.10253173935778,-1.89574959850490,0.337281857599254) q[9];
cx q[9],q[1];
u1(2.86581989005559) q[1];
u3(-2.29365992225949,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.368914764158232,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.44528628214226,-1.09292927488843,4.96293225818235) q[1];
u3(0.996916176240590,-4.76573952607352,1.26241442869181) q[9];
u3(2.47652415756794,1.48986190880196,1.08402032065967) q[6];
u3(1.01195090165938,-1.63103802551148,-1.65059068660599) q[3];
cx q[3],q[6];
u1(0.00163965009323830) q[6];
u3(-0.845944035268920,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.37966431595723,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.49654912269149,2.51874469215768,-0.781614119552715) q[6];
u3(2.85089085315839,5.24314581397172,-0.737142591741567) q[3];
u3(1.40567968344698,-1.31384943468617,1.36034782620880) q[5];
u3(1.62220139839891,-3.15670690822117,0.152189312052137) q[0];
cx q[0],q[5];
u1(-0.498094003584443) q[5];
u3(0.723675905140987,0.0,0.0) q[0];
cx q[5],q[0];
u3(4.43705496602017,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.672673473268100,-2.86754598739465,0.378307100366964) q[5];
u3(2.10151393602135,0.890546106722255,2.05769875691460) q[0];
u3(1.95264699807035,-0.740873976492421,1.06137236087381) q[2];
u3(1.94457937918510,-3.15256482098065,-0.402120789242239) q[8];
cx q[8],q[2];
u1(0.984547871555811) q[2];
u3(-0.581513311525804,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.0560821992018086,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.38010567490604,1.90130314802498,-1.81720542295254) q[2];
u3(1.27672887449858,2.73028637952497,-2.54416625267478) q[8];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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