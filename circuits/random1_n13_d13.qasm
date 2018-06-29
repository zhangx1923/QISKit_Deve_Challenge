OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.65982793642435,-0.788217723075266,0.996087560377781) q[7];
u3(0.982484211525206,-1.80693338734463,-1.52805201551574) q[8];
cx q[8],q[7];
u1(0.971470164100359) q[7];
u3(0.216963768863315,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.74741164362383,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.76855519191339,0.107770543657241,-1.15417801111759) q[7];
u3(1.42702996833338,-3.12495109656116,-1.53043594123032) q[8];
u3(1.71709568437608,-0.690806112873276,1.09787877343549) q[2];
u3(1.33713407484918,-2.84247990583520,-0.526240056425531) q[10];
cx q[10],q[2];
u1(1.60905813399121) q[2];
u3(0.279453428024004,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.702101758731640,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.12721165526504,-1.63693122244803,0.482746145495263) q[2];
u3(1.10395983721846,-0.772808162541995,3.08103435649796) q[10];
u3(2.22448937741847,-0.984245841300494,2.34484166398680) q[9];
u3(2.41528717629066,-1.00345091841737,-0.558592589748849) q[1];
cx q[1],q[9];
u1(2.03792533700739) q[9];
u3(0.185507606459490,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.894723784452978,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.606220193263087,0.854331255869687,-4.27372347228542) q[9];
u3(2.49734061083173,3.78494798979968,0.210006522719236) q[1];
u3(0.926294084153702,-0.648326321553647,1.77018930840688) q[12];
u3(1.63142939404808,-1.26586807963796,-2.86392942933014) q[0];
cx q[0],q[12];
u1(2.11547001584634) q[12];
u3(-1.96357901112038,0.0,0.0) q[0];
cx q[12],q[0];
u3(0.738670984074280,0.0,0.0) q[0];
cx q[0],q[12];
u3(0.432384927297108,4.60040686803153,-0.435757088330648) q[12];
u3(3.00736262713926,-4.40961965235820,0.302837820286174) q[0];
u3(1.95721036051400,1.86676644278981,-0.0179801138067444) q[11];
u3(1.21269351702905,0.464170452465497,-4.26793477904431) q[5];
cx q[5],q[11];
u1(2.87561744890904) q[11];
u3(-2.35086071974441,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.27373774667072,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.68655202449471,3.95306999283374,-1.44019562351259) q[11];
u3(1.82052106064711,1.52419003531681,4.21250331681983) q[5];
u3(1.05046739690540,0.840297788619485,-3.71646998053654) q[4];
u3(0.804376741667803,2.24991422638178,-2.61992743031807) q[3];
cx q[3],q[4];
u1(1.69681513463571) q[4];
u3(-2.35517059974710,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.27546551868178,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.793494848774065,-2.62739720641527,2.37732622048162) q[4];
u3(2.04040054099414,3.12202150098329,-1.09458915238297) q[3];
u3(0.973467826814515,-1.30540923183045,0.140655081919228) q[12];
u3(1.47907824047322,-3.94851088042388,-0.283011394270196) q[10];
cx q[10],q[12];
u1(-0.493101346725775) q[12];
u3(1.22157320037538,0.0,0.0) q[10];
cx q[12],q[10];
u3(3.53452945828920,0.0,0.0) q[10];
cx q[10],q[12];
u3(1.22513455638004,-0.417968889560454,0.900426125128859) q[12];
u3(2.15752614396305,0.465187167775169,-2.28057004045865) q[10];
u3(2.59282644282977,-0.731068439196996,2.35730066748726) q[8];
u3(2.78352585713351,-0.546106733588489,2.38035252208782) q[5];
cx q[5],q[8];
u1(0.428518061760844) q[8];
u3(-0.0554600076455825,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.06857402794159,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.93791217677416,-2.97888735888230,2.08733919263173) q[8];
u3(1.23374084620502,-2.93747764649598,0.716224263010178) q[5];
u3(2.42316174469390,-1.20659608841729,0.800856790135779) q[0];
u3(2.65124936628608,-0.232213625245730,1.69054767977254) q[9];
cx q[9],q[0];
u1(2.80955446656409) q[0];
u3(-4.49231322468820,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.193516327798241,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.14350192017251,-3.49640640050478,1.92412223546073) q[0];
u3(2.52276972774960,-4.53204986434501,-1.23021420569515) q[9];
u3(0.438785726756583,1.31152665604628,-0.346641380083460) q[4];
u3(0.862037945837586,1.18840720360268,-2.40164195979068) q[1];
cx q[1],q[4];
u1(1.89864349356513) q[4];
u3(-2.77702674503519,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.54792464285847,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.50595101842784,-1.75155942273828,-0.609669311190941) q[4];
u3(1.06106835390262,-1.37430545154896,2.39307451127755) q[1];
u3(2.03702815547234,1.62756334168282,1.12787052568136) q[7];
u3(0.157117024873420,-2.77363899583326,-0.992973605612209) q[11];
cx q[11],q[7];
u1(3.47066235040569) q[7];
u3(-1.64266472245891,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.24649623332513,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.27528256109668,0.966477447272322,0.199547630628169) q[7];
u3(1.18498837214898,1.83249994696603,2.66831746484881) q[11];
u3(0.775248668032167,-0.888055791867160,0.510464880026270) q[2];
u3(1.30685581719458,-2.03386929325976,-1.04262907542517) q[6];
cx q[6],q[2];
u1(0.562044756342960) q[2];
u3(-3.29045251627429,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.93645652113899,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.31040053946718,0.340037297150173,-1.19106555037658) q[2];
u3(0.169405224389105,-0.981753875991597,-0.567335701169022) q[6];
u3(1.41139409147975,-0.747856490372643,-0.673356893318179) q[10];
u3(1.61722109265419,-4.37150288254278,1.35608743952468) q[3];
cx q[3],q[10];
u1(4.12219483459038) q[10];
u3(-3.10383419709234,0.0,0.0) q[3];
cx q[10],q[3];
u3(-0.522196421523125,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.717733158426287,1.13934758856434,0.365024642455289) q[10];
u3(1.11404059848563,-0.925934714754718,-2.25229757146631) q[3];
u3(0.997152020382216,0.444897913347698,2.01168779606913) q[11];
u3(1.38789534994338,-1.86777660974745,-1.57150409914425) q[7];
cx q[7],q[11];
u1(-0.220411938535942) q[11];
u3(-1.47399579439630,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.778254472873267,0.0,0.0) q[7];
cx q[7],q[11];
u3(2.44865492505340,-4.30486089590811,0.633412891730659) q[11];
u3(1.11229250814123,0.240707935798531,1.96406969531131) q[7];
u3(1.27387381365388,-1.98045857880900,0.202120068884543) q[5];
u3(1.28792555933568,-3.54581140409710,-0.725834316405706) q[0];
cx q[0],q[5];
u1(1.84911488721706) q[5];
u3(0.163340249655928,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.48207759363734,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.845829638272960,-0.737703699779186,-0.146925576542865) q[5];
u3(1.86346330989383,0.261121987482120,-5.98443131838522) q[0];
u3(1.25277508982590,2.08262971161745,-0.0130001225070524) q[12];
u3(0.916291440447053,0.347210751354036,-4.56661691598445) q[4];
cx q[4],q[12];
u1(0.0926867522263999) q[12];
u3(-1.62264294274476,0.0,0.0) q[4];
cx q[12],q[4];
u3(0.536361825200951,0.0,0.0) q[4];
cx q[4],q[12];
u3(0.934070869853602,2.65975981690783,-1.67797922376203) q[12];
u3(2.10188790585211,-1.54351378869151,-4.17941594590976) q[4];
u3(2.13561887399592,-1.21001810734005,0.305480992269167) q[6];
u3(0.562057398896457,-2.39788572380808,-2.02832710668007) q[1];
cx q[1],q[6];
u1(2.76360797400409) q[6];
u3(-1.94218332828084,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.08440782438757,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.356039902277977,-2.63607743784300,3.16160589771059) q[6];
u3(1.43054400583425,4.27724915093256,1.43402294707755) q[1];
u3(1.94600437914507,-1.69077204542583,-0.655722272030313) q[9];
u3(2.19939115783892,-3.77273315586049,0.145487494807394) q[2];
cx q[2],q[9];
u1(3.47204177948078) q[9];
u3(-1.22130405067134,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.16572817404953,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.08939289311822,2.90320123523112,-0.699359225905069) q[9];
u3(2.01327157898907,0.428746252712891,-0.892920595224796) q[2];
u3(2.00079665645819,2.12602575612354,0.839716828598749) q[6];
u3(2.58621519122834,-0.750561915028630,-2.87429275847785) q[4];
cx q[4],q[6];
u1(2.34012413123589) q[6];
u3(-2.95625565304669,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.36202236707351,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.14916973199046,-0.237652243250795,0.947293408773292) q[6];
u3(0.276936850268831,-4.15945668555950,-0.871003035404268) q[4];
u3(2.56588788907429,-3.96532183151612,2.30896150225297) q[10];
u3(0.944649312339458,1.78624613684738,-0.735455551800054) q[12];
cx q[12],q[10];
u1(3.33737160683013) q[10];
u3(-1.00905738612203,0.0,0.0) q[12];
cx q[10],q[12];
u3(1.95468108748005,0.0,0.0) q[12];
cx q[12],q[10];
u3(1.37398183063466,-1.46272472672513,-2.26171268127367) q[10];
u3(1.89293417404205,-0.761509783509163,4.43216372958194) q[12];
u3(2.34028968680275,2.06529890117328,0.533579180036521) q[0];
u3(2.09376971285974,0.393292556990796,-3.91952002837537) q[1];
cx q[1],q[0];
u1(1.09029288006523) q[0];
u3(-3.42948239786850,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.55045272980020,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.85916202726605,-3.81164713153459,0.424270874412588) q[0];
u3(1.21928993054103,-0.714452865145971,-4.27138064798225) q[1];
u3(2.66476372045746,2.49721994342577,0.357871722312914) q[11];
u3(2.42104806870600,1.77070744021549,-2.69993834002671) q[9];
cx q[9],q[11];
u1(1.54394212013701) q[11];
u3(0.617470953073682,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.906724406983667,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.36707363675646,0.456104496790350,-3.04791898856884) q[11];
u3(2.13063445708053,1.76859053849725,2.71717061576821) q[9];
u3(1.81740498688774,-1.45976159794443,1.72007428666694) q[3];
u3(2.00340431990615,-1.75241539951578,-2.43938986771174) q[5];
cx q[5],q[3];
u1(-0.220115102186975) q[3];
u3(-1.54092951761107,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.919340948875798,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.336455929193436,3.31076827914104,-0.971384370943643) q[3];
u3(2.90909803369760,3.06143364024054,1.46531065430805) q[5];
u3(1.00936132506390,2.41870059329250,0.0955710380692920) q[8];
u3(1.12700002063334,-0.444357628777878,-2.31515915051200) q[2];
cx q[2],q[8];
u1(2.02880483375803) q[8];
u3(0.0635412657472585,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.669776880050608,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.826041921793921,-1.00187106475599,0.147105149543134) q[8];
u3(1.27163602992192,2.23734404275148,-0.817056623129375) q[2];
u3(2.11054675153639,1.06116991780948,2.00353490392301) q[2];
u3(1.36507369036281,-2.21638760568457,-2.83940942201928) q[3];
cx q[3],q[2];
u1(0.169005814487386) q[2];
u3(-1.64341118362184,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.514794196861942,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.905158416564298,-2.04580396349375,0.0894299770415762) q[2];
u3(0.720972774707418,1.50462349464479,-4.16696163292723) q[3];
u3(1.27251309506771,2.74429713667105,-1.98281272239959) q[1];
u3(0.947015252795891,1.28977094679185,-2.13956864896736) q[4];
cx q[4],q[1];
u1(1.26460358003211) q[1];
u3(-0.890204463096398,0.0,0.0) q[4];
cx q[1],q[4];
u3(-0.115773486966827,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.23823242596053,2.09001662499639,-0.914417472828241) q[1];
u3(2.19795166728109,-1.68871194116265,0.363131122754103) q[4];
u3(1.48191548235517,-0.686008412620465,1.25260085060500) q[7];
u3(1.05165954800355,-1.31829610331642,-0.358894419226286) q[5];
cx q[5],q[7];
u1(1.58636205193599) q[7];
u3(-3.31988675927054,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.01053233237156,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.37213125489456,-0.672574709070851,3.21913014004550) q[7];
u3(0.670609007643236,1.09241731358128,0.0783203782856694) q[5];
u3(0.912374767297580,-2.55492656302952,1.02223044233800) q[10];
u3(0.892411378052665,-2.16625097862601,-0.332203079323851) q[12];
cx q[12],q[10];
u1(1.71076127974229) q[10];
u3(-3.19804919826318,0.0,0.0) q[12];
cx q[10],q[12];
u3(1.09556074015991,0.0,0.0) q[12];
cx q[12],q[10];
u3(1.35003739246227,-2.96042779408205,2.03664070846499) q[10];
u3(1.72992689260915,0.0802320667871919,5.18917860970895) q[12];
u3(0.133250485055467,1.28339396083381,-1.36753592976785) q[11];
u3(0.916113104469579,-3.23651639257683,2.88429390403743) q[0];
cx q[0],q[11];
u1(1.75791432863395) q[11];
u3(-3.29590458432280,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.07778086789192,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.91930645950034,-0.897940835333634,-0.128853892916773) q[11];
u3(0.695430610495632,0.853232167192423,0.799407152603843) q[0];
u3(1.65763059060536,2.92090487212072,-2.69114603134852) q[6];
u3(1.84903552418013,1.02353103153596,-1.76056891205586) q[9];
cx q[9],q[6];
u1(0.413055719826657) q[6];
u3(-1.17755883424344,0.0,0.0) q[9];
cx q[6],q[9];
u3(1.61773718095073,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.88719965457944,0.182808515962808,1.73265254361425) q[6];
u3(2.02590640651849,2.93294743163444,-0.421031063608819) q[9];
u3(1.43996151868833,-3.56929179817130,2.46610546264632) q[1];
u3(0.460429526542743,3.08411899230654,-1.91041369471642) q[2];
cx q[2],q[1];
u1(1.54526089707757) q[1];
u3(-0.157083938143784,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.30170439236458,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.19092504478749,0.884255351817514,-2.24093349922698) q[1];
u3(2.42670119889019,-0.707510517379202,-0.609489517193035) q[2];
u3(1.88303194493842,3.25130370224078,-0.486994011475774) q[5];
u3(2.57659452058284,3.03226680475827,-0.359332466073153) q[8];
cx q[8],q[5];
u1(2.88287040082434) q[5];
u3(-2.23470021292383,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.50990146335435,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.20184015100119,2.18211112017855,-2.27374409623691) q[5];
u3(2.37968278264280,-0.364978264314705,-2.16939851092740) q[8];
u3(1.51371948004034,-0.649055918867540,-0.000841536380228325) q[12];
u3(1.87811255838640,-3.72985067907319,-1.49704415047395) q[3];
cx q[3],q[12];
u1(4.13651823162413) q[12];
u3(-3.21809738009884,0.0,0.0) q[3];
cx q[12],q[3];
u3(-0.618925020431438,0.0,0.0) q[3];
cx q[3],q[12];
u3(2.50929295996465,2.05577052053085,-0.862921634845802) q[12];
u3(0.782398094798574,1.87024826153568,1.07807203155461) q[3];
u3(2.36083539682217,-1.35363741474830,2.32814725624513) q[7];
u3(2.53513703706509,-3.07473703479022,-1.21490237869834) q[11];
cx q[11],q[7];
u1(2.62070319115919) q[7];
u3(-1.74970210417756,0.0,0.0) q[11];
cx q[7],q[11];
u3(0.709970296761197,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.65414508176957,-2.07907080032875,3.54391981633251) q[7];
u3(1.31713869220504,-0.401328607770012,-4.10700153812127) q[11];
u3(1.98584339610964,-0.657967245393510,1.45256716086570) q[4];
u3(1.67975157458009,-2.26566744502999,-1.08001798307937) q[0];
cx q[0],q[4];
u1(-0.402099060467890) q[4];
u3(-1.64801563814298,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.867753582335116,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.27833759411375,3.10712747977772,-1.94289096438012) q[4];
u3(2.09967956759156,-1.23199559175196,-0.616497191498937) q[0];
u3(1.12690028111926,-0.375294508423274,-2.31510337053278) q[6];
u3(2.02986659935425,2.51259166081550,-3.74289710188907) q[9];
cx q[9],q[6];
u1(2.41229725920497) q[6];
u3(-1.63979512258629,0.0,0.0) q[9];
cx q[6],q[9];
u3(3.57990432236762,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.27781722041457,1.56321157567758,-2.24301664571432) q[6];
u3(2.13228175562922,-4.60201494520045,1.62580660259227) q[9];
u3(2.26543938801295,-1.64351687346250,4.52378948352316) q[9];
u3(0.435109914468714,3.12788001161967,-1.59240300524771) q[7];
cx q[7],q[9];
u1(1.72812798761664) q[9];
u3(-2.16879227039282,0.0,0.0) q[7];
cx q[9],q[7];
u3(3.33982667315677,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.77199252025765,2.64040440174763,-1.42864854469918) q[9];
u3(2.61792669208812,-0.133746521347720,-0.946058925154368) q[7];
u3(2.04069766050865,-2.92134585087448,0.00107537029824534) q[5];
u3(2.84402705217603,1.98155129215759,2.95668088706442) q[10];
cx q[10],q[5];
u1(1.30997910451131) q[5];
u3(-4.00464464223733,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.99052543187331,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.982837503416172,1.18300699835269,0.452983229560127) q[5];
u3(1.59447593011907,-1.60462687508892,-1.78315712525083) q[10];
u3(1.91126279432598,-1.09809215173725,-0.551840922873524) q[0];
u3(1.70021753188261,-2.36724915112660,1.07904903732176) q[3];
cx q[3],q[0];
u1(0.569838430865538) q[0];
u3(-1.25887718386136,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.44223483884812,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.10096782467397,-1.79676287674271,-0.129300517858686) q[0];
u3(0.757622102039829,-1.86340850607029,4.05427548346563) q[3];
u3(1.87436484261506,1.35222971843223,0.698656566033449) q[11];
u3(2.50150843268961,1.20795353358536,-1.97406995711850) q[8];
cx q[8],q[11];
u1(1.18594305396361) q[11];
u3(-0.113477503715378,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.19860033800804,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.56415265800027,-2.73164092135336,0.462844390574045) q[11];
u3(1.66773239210382,-0.284884083102955,1.59776704505847) q[8];
u3(2.85212825473156,2.70236738321260,-1.47859278208966) q[2];
u3(2.42071259591172,2.61520450551181,-2.03160318774867) q[6];
cx q[6],q[2];
u1(2.31463095685661) q[2];
u3(-3.11622206949633,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.47153463173659,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.51922036972585,-0.0300317118436879,-2.91765915022737) q[2];
u3(1.11069385383898,-0.0103505504499243,1.35205893107902) q[6];
u3(1.82850470603596,0.433049214257787,1.73471300048678) q[12];
u3(0.686722531468172,-0.695877011770413,-0.660825245064295) q[4];
cx q[4],q[12];
u1(3.16004602460745) q[12];
u3(-1.52783479198561,0.0,0.0) q[4];
cx q[12],q[4];
u3(2.36110886673900,0.0,0.0) q[4];
cx q[4],q[12];
u3(1.32049538877222,-4.54328892766048,1.54353584617671) q[12];
u3(1.12511601661818,-1.22871191887024,4.58820724220097) q[4];
u3(0.409698849239466,-2.57527731553715,3.30887582414239) q[7];
u3(1.35184507130793,-3.19942933921803,1.34192466459780) q[3];
cx q[3],q[7];
u1(0.0916688299350825) q[7];
u3(-1.77249074172885,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.663249299450543,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.64004998127321,-0.131823851217959,0.802370192266999) q[7];
u3(1.34030541500562,-1.54900602657512,2.99559824949852) q[3];
u3(0.104664678576882,-1.92057721605259,1.09698545455801) q[2];
u3(1.27183576170550,-3.30152160295596,1.71426868118429) q[0];
cx q[0],q[2];
u1(3.11717674502477) q[2];
u3(-1.14880419875549,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.52954134742614,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.81254497377573,-0.353011008187191,1.59586237960637) q[2];
u3(2.22311978102792,0.366414703510357,4.34783580751088) q[0];
u3(1.90996790482148,0.224005911077330,-0.779883755107141) q[1];
u3(1.31323896137814,0.683860377513552,-4.60761267706098) q[8];
cx q[8],q[1];
u1(3.77168548761374) q[1];
u3(-4.23294361215045,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.154047297856812,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.408570635339577,-1.58791327170986,4.04925519037410) q[1];
u3(2.76288692847739,1.46123155591922,-0.667833609636772) q[8];
u3(1.94526565953349,1.84347816912328,-2.59018515876316) q[10];
u3(1.41863433961693,1.93744160272920,-2.22087827380509) q[6];
cx q[6],q[10];
u1(3.12100932054613) q[10];
u3(-1.75495898067211,0.0,0.0) q[6];
cx q[10],q[6];
u3(0.579196760974721,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.14682674940552,1.80457398538049,-2.98656268733896) q[10];
u3(0.567261503332715,-5.06372448591809,1.15213913149609) q[6];
u3(0.699900037259347,-0.526911395328048,-1.38434376819563) q[5];
u3(1.88068640959929,0.399878917834767,-5.37723198414732) q[12];
cx q[12],q[5];
u1(1.69243565922656) q[5];
u3(-2.10600186406516,0.0,0.0) q[12];
cx q[5],q[12];
u3(3.93502615917843,0.0,0.0) q[12];
cx q[12],q[5];
u3(2.38991873363570,0.322611104004155,-3.68245248257732) q[5];
u3(0.561777618342923,-4.69815278530510,-0.927707196137512) q[12];
u3(0.660063517380816,0.995947967012194,-0.200670125192629) q[11];
u3(0.469362854179472,-2.16165875979449,0.105838591313147) q[9];
cx q[9],q[11];
u1(2.93463661181323) q[11];
u3(-1.66937954075548,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.867356053690154,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.82692617754561,0.288251492459683,0.702160586714000) q[11];
u3(0.321345696863237,0.462831453561414,0.481775148329201) q[9];
u3(1.69436571838594,0.826378700224390,0.180430520360006) q[10];
u3(0.819611863671183,0.0722378248790594,-2.05925804112439) q[8];
cx q[8],q[10];
u1(1.68865131290894) q[10];
u3(-2.53410536344725,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.145336914049580,0.0,0.0) q[8];
cx q[8],q[10];
u3(0.972540336562016,-0.221521810012440,1.58498600349562) q[10];
u3(0.985707550921268,-4.17111680220359,-1.68739636925387) q[8];
u3(2.81341980214289,0.0886139333580240,3.04564250818110) q[9];
u3(1.82937589684976,-3.45903741737638,-1.05699699875980) q[4];
cx q[4],q[9];
u1(1.44136602211038) q[9];
u3(-3.59293056130467,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.75939641556867,0.0,0.0) q[4];
cx q[4],q[9];
u3(2.40637253076486,-4.58027052172178,0.915906790569132) q[9];
u3(0.940513442103993,-1.48935278721951,-4.39999257531867) q[4];
u3(2.67993628653187,2.05576248050098,-4.22725514409387) q[7];
u3(1.41164090712818,2.18642654972107,-0.951654503777503) q[1];
cx q[1],q[7];
u1(0.741008887487038) q[7];
u3(-1.38219536985889,0.0,0.0) q[1];
cx q[7],q[1];
u3(3.01931077376496,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.631150069546750,0.779362103994934,-1.18104451060548) q[7];
u3(0.527716935447473,3.69416284037571,-1.13582330675073) q[1];
u3(1.97656703174747,-0.673662984785770,-0.244123985211789) q[0];
u3(1.75449092385174,-3.08435457441214,0.902733602699572) q[12];
cx q[12],q[0];
u1(1.34294466333755) q[0];
u3(-0.475354235698454,0.0,0.0) q[12];
cx q[0],q[12];
u3(2.12379067248937,0.0,0.0) q[12];
cx q[12],q[0];
u3(2.03883881567941,0.451935242944211,1.32696421051279) q[0];
u3(1.75981960060101,-1.79005656003831,-3.84656187641274) q[12];
u3(2.57968613023387,-2.60154550111242,3.34089166901153) q[11];
u3(1.19642110678654,3.75306207062028,-1.94454299464656) q[2];
cx q[2],q[11];
u1(0.641293581316382) q[11];
u3(-1.03318457903262,0.0,0.0) q[2];
cx q[11],q[2];
u3(3.26388167968926,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.687345758386826,-1.28902974272625,4.19636468790875) q[11];
u3(1.55743799115111,-3.13193017627616,0.975497163493419) q[2];
u3(1.80488866507161,2.16178651762155,-3.60478259030568) q[6];
u3(1.54738219984352,3.22640291711708,-2.62404245650816) q[3];
cx q[3],q[6];
u1(-0.338065516409973) q[6];
u3(-1.09847250294626,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.53474594686196,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.53078317510834,-1.97635856529499,3.45921001831145) q[6];
u3(0.906256727483842,2.60215591728235,-1.61857913402292) q[3];
u3(0.842066223873914,2.57346911364328,-0.185439155737572) q[12];
u3(1.89869438458463,0.174225890560506,-4.00782017219328) q[7];
cx q[7],q[12];
u1(1.48891235649173) q[12];
u3(-0.652649500713387,0.0,0.0) q[7];
cx q[12],q[7];
u3(-0.527791397004765,0.0,0.0) q[7];
cx q[7],q[12];
u3(1.32290978151341,1.66937157279004,-0.0993079287210099) q[12];
u3(2.17974153845051,2.42389005489069,-0.607061305640980) q[7];
u3(2.39926989119535,-2.36899367179935,0.415945336716511) q[0];
u3(1.45830450515303,-2.90267312453661,0.537747351448698) q[4];
cx q[4],q[0];
u1(2.45045918668558) q[0];
u3(0.100929367973138,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.13703317184907,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.45910688842497,2.88017185503967,-0.759867231730577) q[0];
u3(1.02135272400186,5.76512116657925,-0.0668056330728075) q[4];
u3(0.806910759179685,1.38630752070842,-1.28975023682964) q[11];
u3(0.831349303198222,-1.85130942958839,0.949326767206730) q[2];
cx q[2],q[11];
u1(0.548222360420050) q[11];
u3(-1.58720237547038,0.0,0.0) q[2];
cx q[11],q[2];
u3(-0.168484433476994,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.811514754069135,-1.85782007141225,1.38954822563065) q[11];
u3(1.38302339800567,-2.70199359386276,-2.68944811565542) q[2];
u3(2.16195841001400,0.586347767631041,0.471024315878595) q[1];
u3(1.88315338577542,-1.91269050667548,-1.71172036817337) q[8];
cx q[8],q[1];
u1(-0.00446332033410646) q[1];
u3(-1.50313319477430,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.47935446408349,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.29585809747862,-2.11797626850163,3.26517001792665) q[1];
u3(2.05604718238221,0.512709611640201,-0.137208418348009) q[8];
u3(1.68807432289805,-1.66302391218465,-0.992543795043752) q[9];
u3(1.29128069169836,-4.29981045715993,0.425674819092208) q[10];
cx q[10],q[9];
u1(1.66778266877241) q[9];
u3(-2.29260379816819,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.446340483486540,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.19514979301994,-3.96276059641359,2.08017023033418) q[9];
u3(2.07617961280107,0.836846933230575,3.65193547788472) q[10];
u3(1.82777994716871,2.62236058199834,0.133446665722696) q[3];
u3(1.50620265554717,0.159927522797774,-2.57048418967426) q[5];
cx q[5],q[3];
u1(1.69573911000818) q[3];
u3(0.286706444780072,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.661685851765500,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.07302450621271,0.667531299970486,-1.78493889045574) q[3];
u3(2.38844996366634,2.09968795585450,-1.52243034373672) q[5];
u3(0.815104843375436,2.90541471597662,-2.96788665719544) q[7];
u3(2.01486933155956,3.49770086956151,-2.57138631511002) q[11];
cx q[11],q[7];
u1(1.05770165897988) q[7];
u3(-1.55933557651725,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.25330131112874,0.0,0.0) q[11];
cx q[11],q[7];
u3(0.905008939469337,-2.46228902498497,-2.01624915720636) q[7];
u3(0.365843636007910,-5.91586599740439,-0.00232661098532194) q[11];
u3(1.48689421016139,-0.164427710953296,1.63976801134979) q[6];
u3(1.25199216073484,-0.644272614596798,-1.16565288169186) q[5];
cx q[5],q[6];
u1(-0.0604932868631749) q[6];
u3(-0.962085339578057,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.43288645366246,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.29374685700779,3.73163958504142,-2.04972722546857) q[6];
u3(1.19424372422256,1.31733567774287,1.75867512928664) q[5];
u3(2.12197379817494,1.89901528722055,-2.13238544129718) q[8];
u3(0.548271615394734,2.14748908865375,-3.17048903099161) q[0];
cx q[0],q[8];
u1(0.967959001753217) q[8];
u3(-1.26331618853365,0.0,0.0) q[0];
cx q[8],q[0];
u3(-0.115077051560733,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.10580272476361,2.30274818878405,-2.36654890375769) q[8];
u3(1.43709284954433,3.03511803673794,-1.39472504317189) q[0];
u3(0.450710458188577,2.67232227246460,-2.39682177868496) q[2];
u3(0.656093397145387,-4.07264112518949,1.80495230452651) q[3];
cx q[3],q[2];
u1(0.892129325806464) q[2];
u3(0.286719077666858,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.61280711258658,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.71475932661029,-0.932075601959613,1.56758115011452) q[2];
u3(2.03249242907963,2.80624392534338,1.09625450517125) q[3];
u3(2.46059923877647,2.66723263932437,0.430461015464954) q[4];
u3(1.71318843830438,0.661473085606618,-3.39650877756707) q[12];
cx q[12],q[4];
u1(2.34737755344846) q[4];
u3(0.224024880764586,0.0,0.0) q[12];
cx q[4],q[12];
u3(1.52737521928600,0.0,0.0) q[12];
cx q[12],q[4];
u3(1.67992681375414,0.564477545635080,-1.89152854424262) q[4];
u3(1.80742370986418,-4.24375764010175,0.411553678694175) q[12];
u3(1.55834386220071,3.07728365689920,-1.27653081459530) q[10];
u3(2.76643257550532,2.38399927897041,-0.745139382015800) q[9];
cx q[9],q[10];
u1(1.13696965326454) q[10];
u3(-1.41323986803214,0.0,0.0) q[9];
cx q[10],q[9];
u3(-0.485756722258962,0.0,0.0) q[9];
cx q[9],q[10];
u3(2.74812392235747,1.17183495853957,-4.11437549180137) q[10];
u3(2.15205262501007,-3.46912170721399,2.28080174677937) q[9];
u3(1.88129524053414,-1.65813560209967,-0.787174296243392) q[5];
u3(1.26530953448550,-4.67209268742526,0.374829678036344) q[9];
cx q[9],q[5];
u1(1.62240252775702) q[5];
u3(0.362980752133212,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.697178598294624,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.52637258410835,-0.416055971714465,3.04033685081644) q[5];
u3(1.53493819968913,-2.21089567950668,0.0526489777210126) q[9];
u3(0.757208319191753,1.66083764834591,-2.51810286132628) q[8];
u3(0.447126034417174,0.756232177648564,-3.09204084200450) q[6];
cx q[6],q[8];
u1(0.307854979483841) q[8];
u3(-1.37443353954045,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.99873461226690,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.374265431323080,0.990646092973447,-1.11921366523444) q[8];
u3(1.29613914328508,0.490222529429586,2.42056405770935) q[6];
u3(1.01087499422037,2.96015714198044,-1.60977688461150) q[7];
u3(0.366998451834849,1.20562313458910,-2.53820410392297) q[0];
cx q[0],q[7];
u1(-0.683329063700346) q[7];
u3(1.17645571119101,0.0,0.0) q[0];
cx q[7],q[0];
u3(3.61346952954023,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.57515251087875,-0.0284128967042842,0.000211259669550080) q[7];
u3(1.78317632667273,0.799832701712857,-1.60780440722593) q[0];
u3(2.24047753794574,2.18194379859655,0.606569077189319) q[10];
u3(1.89449254360276,0.283293670145272,-2.75375065091748) q[12];
cx q[12],q[10];
u1(2.19744554213465) q[10];
u3(-1.66823446882215,0.0,0.0) q[12];
cx q[10],q[12];
u3(3.14857707319137,0.0,0.0) q[12];
cx q[12],q[10];
u3(0.743413745015649,1.71430135394828,-2.00039313630800) q[10];
u3(1.42727398238882,-0.290684837195566,-0.443042072226780) q[12];
u3(2.16091051533950,0.249431708795046,-1.13708685818930) q[3];
u3(1.32527985461433,-4.04960325288785,1.65942607243025) q[11];
cx q[11],q[3];
u1(2.24948198767375) q[3];
u3(-1.71539336415571,0.0,0.0) q[11];
cx q[3],q[11];
u3(1.06428107874468,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.39405949998565,-3.12108430340190,2.70206604254592) q[3];
u3(2.00791902785922,0.596350641152388,-3.67996886599598) q[11];
u3(0.944792162477306,1.19215584034621,-0.848873747407056) q[1];
u3(0.917627255628365,-2.78972404712808,0.820842682952813) q[4];
cx q[4],q[1];
u1(-0.203000474193290) q[1];
u3(-1.25094826816953,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.29993551152412,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.53538309443274,-2.51874345248055,-1.54242191450144) q[1];
u3(1.30675818544607,2.16762046297128,-0.397530448512490) q[4];
u3(0.930287231080929,0.647207158155878,-2.17229679134965) q[9];
u3(1.05675349919488,-3.55480483365771,2.26754332396318) q[1];
cx q[1],q[9];
u1(2.47237196173442) q[9];
u3(-1.98432144481718,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.45801257912313,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.33781500126012,1.07697525181467,-0.419733562640053) q[9];
u3(0.870508101568198,0.817752035671340,5.20869709739047) q[1];
u3(1.55826018278400,-3.41150996058142,0.832435016762465) q[0];
u3(0.275710416770336,-1.85057925706742,1.02473184668599) q[6];
cx q[6],q[0];
u1(2.50713073320843) q[0];
u3(-2.80873417531678,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.825154860452604,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.07186014329838,-1.56684581641056,-2.25433209648332) q[0];
u3(1.01105885855817,0.986926233110275,-1.01174288665667) q[6];
u3(0.675834777044956,-1.75074358566132,1.88454567858007) q[5];
u3(0.481211970649947,-2.75769461664802,-0.273368036547647) q[3];
cx q[3],q[5];
u1(2.07143131067997) q[5];
u3(-2.25571014671762,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.11015065771420,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.05238584600130,0.684067006900552,-0.870733186923934) q[5];
u3(0.435261859517680,0.731889896154731,4.43898752929938) q[3];
u3(2.60546928292105,1.77938328489247,-0.0350471382669841) q[11];
u3(1.74995471229857,0.0390021746484537,-3.31900529690786) q[4];
cx q[4],q[11];
u1(1.87503942805854) q[11];
u3(-2.75023210322093,0.0,0.0) q[4];
cx q[11],q[4];
u3(0.358319173196014,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.51005426835800,-3.03190310074822,1.32314152361462) q[11];
u3(1.78512506504618,-2.82486222703182,2.80047912833797) q[4];
u3(1.36739035775546,0.693006213776904,-2.86390735885898) q[10];
u3(1.88811347616534,2.60320935251859,-2.72707284400639) q[2];
cx q[2],q[10];
u1(2.51239572982129) q[10];
u3(-1.54622307199937,0.0,0.0) q[2];
cx q[10],q[2];
u3(0.176564128403326,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.54838494771976,0.914088443147232,-1.05226334070288) q[10];
u3(2.47681403683857,-1.58798494296457,3.73151654166656) q[2];
u3(1.54199003763445,1.14137454349384,-3.01532322601362) q[7];
u3(2.45311833435455,2.88906929280339,-3.23637110072473) q[8];
cx q[8],q[7];
u1(3.10240982418270) q[7];
u3(-2.19255922552210,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.12271057315612,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.09719588446916,0.193349892649277,1.23649883661476) q[7];
u3(1.08433166539079,1.09455138778335,-3.14016969192846) q[8];
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
