OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.24320063325283,2.13730884566055,0.401155033558822) q[5];
u3(2.21770671128584,-0.0458636757062147,-2.11721187751861) q[3];
cx q[3],q[5];
u1(3.39376089331946) q[5];
u3(-1.37005855190114,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.42953945026950,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.26198662980434,-3.48265248813188,2.52850466210294) q[5];
u3(1.32014324190271,2.15933396915821,2.57044670379292) q[3];
u3(2.08410835616337,2.10870728604024,-0.778823529685609) q[2];
u3(2.39092551290568,0.665891772301781,-2.91055408132676) q[4];
cx q[4],q[2];
u1(3.40171328496518) q[2];
u3(-3.50658218663766,0.0,0.0) q[4];
cx q[2],q[4];
u3(-1.26980622078861,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.41346038062175,-2.24087236220854,1.48611187382558) q[2];
u3(2.41627664783380,-0.550768802085186,-0.522940701364433) q[4];
u3(2.03230743047718,-1.34518782196604,-0.554751723959575) q[1];
u3(1.65225688952715,-3.47934506262063,-0.469285984201885) q[0];
cx q[0],q[1];
u1(0.573467461174186) q[1];
u3(-1.22966567795574,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.14388302257386,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.57487004686541,1.93005730830791,-1.86792812725945) q[1];
u3(2.10497664776334,5.25343279438318,-0.587960196912071) q[0];
u3(2.56373082119933,-1.14618512498889,-0.839557150196462) q[4];
u3(0.568735720727356,-5.37527061220355,0.419225388363439) q[1];
cx q[1],q[4];
u1(2.33603836877899) q[4];
u3(-2.94449890053786,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.07464931299541,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.46100467033549,2.63862125529267,1.24077284605018) q[4];
u3(1.89473180166150,-0.615400741610189,1.54515513356182) q[1];
u3(2.03130707969538,1.93463415861502,-0.323709474964831) q[5];
u3(2.24377346441991,-0.553085884508494,-4.02951434082394) q[3];
cx q[3],q[5];
u1(0.341991443115146) q[5];
u3(-0.630429627477020,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.03329444575319,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.51385809643606,-3.91583420775855,1.90642480957372) q[5];
u3(1.40569595856730,2.38873536380470,3.53817787346563) q[3];
u3(0.788803588314031,-2.68350418743319,2.08574874249424) q[0];
u3(0.806773819748027,1.37861231316844,-2.79490982245626) q[2];
cx q[2],q[0];
u1(1.49056673506945) q[0];
u3(0.0279996116240744,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.26594554107248,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.64584559037430,0.719666877476753,-2.92784947018438) q[0];
u3(2.68129019152937,1.03083029821854,3.91409114572133) q[2];
u3(0.618839808021160,2.28832470379254,-2.93493145886482) q[4];
u3(0.715681541095934,-0.325672929591491,-1.03935141509732) q[2];
cx q[2],q[4];
u1(2.47033961979994) q[4];
u3(-1.51388325191082,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.0183496495655515,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.890113121461472,-2.68661756093660,1.44121455507100) q[4];
u3(1.79524794113950,-2.86689916215636,-2.50007169659229) q[2];
u3(2.74183952057265,-3.08880798631387,2.95226203541149) q[0];
u3(0.916852717912897,-2.21961555184865,3.65311694538243) q[5];
cx q[5],q[0];
u1(-0.100239066771680) q[0];
u3(0.932230623971850,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.78781117136015,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.81743126840208,-2.32587890254282,0.447334603332405) q[0];
u3(1.05197827541593,-1.36178742816088,-4.15743014690562) q[5];
u3(1.83175737526773,1.19591876543021,1.69528812398104) q[3];
u3(1.80944550234093,-1.47183968805076,-2.00397356374169) q[1];
cx q[1],q[3];
u1(-0.164476771643249) q[3];
u3(-1.72341931388333,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.23253877392702,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.42766024085248,-0.938314009612439,-0.411433552341781) q[3];
u3(0.952068488699818,0.476562653193308,0.134750962721259) q[1];
u3(2.00573244245596,-0.0765832797818495,3.09030431532139) q[5];
u3(2.59177350000466,-3.50731575897996,-2.56801699875174) q[4];
cx q[4],q[5];
u1(2.42025635609784) q[5];
u3(-2.00316817348256,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.25435998688563,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.50367117255496,4.54716212932853,-0.889197363809399) q[5];
u3(1.59575785703586,3.49138852157505,-2.76540894585728) q[4];
u3(0.769614317519047,-0.518005502990670,0.595297767560055) q[2];
u3(0.749154768401236,-2.77601897374731,2.52852960863661) q[1];
cx q[1],q[2];
u1(0.645948600843885) q[2];
u3(-1.42669220828930,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.0461359414768963,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.44094604361879,-2.72292511595695,1.93884129348034) q[2];
u3(0.337271540306634,1.95878983645389,-2.83420013836008) q[1];
u3(0.338367743680170,1.94491746619264,-3.64622694556924) q[0];
u3(1.49996453381737,3.03326032271217,-2.17596001265184) q[3];
cx q[3],q[0];
u1(2.80106405996298) q[0];
u3(-1.83474529284490,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.843229805038235,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.26350050988901,2.21096189310093,-2.39209053945331) q[0];
u3(0.978450992101295,0.611884703603733,-1.56412229486387) q[3];
u3(1.67203343886792,4.19179359016925,-1.92072735553560) q[5];
u3(1.57478650410172,1.91362092098082,-2.63013861627346) q[3];
cx q[3],q[5];
u1(1.21309815815334) q[5];
u3(-0.197834498510775,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.76208221551435,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.666150553121621,0.529135706488246,-1.10813229610551) q[5];
u3(2.08899384194504,-0.144340783769155,3.31024891639666) q[3];
u3(2.62310391392134,1.36331930054635,-1.28833271373913) q[1];
u3(1.92763607589030,4.72383948858288,0.643161254363596) q[2];
cx q[2],q[1];
u1(-1.11837442038404) q[1];
u3(0.898788803283670,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.67522842485636,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.08636436707909,4.20543896436269,-0.671620928690675) q[1];
u3(1.61565087409699,-0.506117825162327,3.04813415556716) q[2];
u3(2.36022819353988,-2.96315190503878,3.00117426082618) q[0];
u3(0.439141399517647,-2.51093148847161,3.24782907737333) q[4];
cx q[4],q[0];
u1(0.467128147947872) q[0];
u3(-1.24232083253741,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.06094035827934,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.67438115474448,-0.230574241055104,1.54480332497975) q[0];
u3(2.13760793332300,-0.791382259285315,-4.63667510498813) q[4];
u3(1.57325064462329,0.475767766523671,2.24571595059527) q[1];
u3(1.90113110069683,-0.226952781189182,-0.909212287429329) q[5];
cx q[5],q[1];
u1(1.07373966304813) q[1];
u3(-3.55244694961621,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.68338106705594,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.196830418539992,-0.325617312975452,-1.71606730443279) q[1];
u3(1.25526779490783,0.529466303919708,-1.88257940370852) q[5];
u3(1.77849382397004,0.517300689814465,0.968223886766280) q[4];
u3(0.145192235647300,-4.27082418073419,-0.890103285587663) q[0];
cx q[0],q[4];
u1(0.207739502340271) q[4];
u3(-0.914650871674536,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.98940286808462,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.56181105956367,2.34727788734340,-0.0667640179977746) q[4];
u3(1.87896836646951,-0.890072493277529,-0.150197467886918) q[0];
u3(1.94408166405849,3.03750462545891,-1.75904265429487) q[3];
u3(2.57574325018864,1.84979822632106,-0.663530802580288) q[2];
cx q[2],q[3];
u1(2.58511897448890) q[3];
u3(-1.86531812588477,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.238601879817800,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.53315908334197,0.207832707744679,2.73298122025094) q[3];
u3(1.72167910946699,0.909072194807336,-2.65690827146759) q[2];
u3(1.68002328446619,0.480692931374397,1.51049943761637) q[2];
u3(0.616437227585260,-5.21931536906457,0.200507849549755) q[4];
cx q[4],q[2];
u1(0.136372448872674) q[2];
u3(-1.33242402110774,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.59680136686726,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.77930362591215,2.50434997451952,-1.21548934516461) q[2];
u3(1.54950165359603,-1.14972758369863,5.06459022282210) q[4];
u3(0.800761479830714,1.72998674264035,-3.11378561656317) q[1];
u3(0.689225035937924,1.51468698259730,-2.89505816412572) q[0];
cx q[0],q[1];
u1(3.10028915042069) q[1];
u3(-0.571839065410080,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.33341316492101,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.35548688011608,-2.21318373328362,0.00893880284855952) q[1];
u3(0.211380308285250,-1.76233336141282,-3.70536308493591) q[0];
u3(2.68162347038536,0.285842017277345,-1.53745028888641) q[3];
u3(2.19316917260245,4.01364435881003,-0.636565166049698) q[5];
cx q[5],q[3];
u1(1.82423015967996) q[3];
u3(-3.18109813111452,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.522325934371473,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.84705855410433,1.56177683587534,-0.465319993047661) q[3];
u3(1.93029150053740,-2.02819651829099,-3.88780742206027) q[5];
u3(0.740299377856793,-1.14347929874921,0.860802155484482) q[5];
u3(1.13134588781443,-2.22946498495054,-0.116319462416016) q[1];
cx q[1],q[5];
u1(0.288655032726167) q[5];
u3(-1.66455067686480,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.44598667806352,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.331112956319785,-0.393944832437754,3.64772610627680) q[5];
u3(1.10198384354440,-1.03088834018881,0.305926752068021) q[1];
u3(2.69941201686659,-0.395895732618359,-0.0572327100559613) q[0];
u3(0.921191355121779,-3.36246356102366,-0.250550516231632) q[3];
cx q[3],q[0];
u1(1.92546950118446) q[0];
u3(0.122417199375465,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.743183976974072,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.857775743892321,-2.99288688987077,2.35602020757064) q[0];
u3(0.818715875011074,-3.08188636984360,-0.218822396729729) q[3];
u3(2.11300598651026,-0.976765223681322,3.60017002283910) q[4];
u3(1.99222137668628,1.39438415260304,1.34430188429840) q[2];
cx q[2],q[4];
u1(1.97723676680891) q[4];
u3(-2.96059052434683,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.36361635158769,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.506536242899720,2.75865653198471,-2.89090269860125) q[4];
u3(2.07269219081603,-2.13115771339270,-1.64068790131993) q[2];
u3(0.672938096875287,3.61222636813722,-1.74881376687013) q[3];
u3(2.24985773051646,2.11441070839713,-1.94868089205670) q[5];
cx q[5],q[3];
u1(2.63947211820360) q[3];
u3(-1.82309775560090,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.749988567113540,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.517691703577344,1.40243545385365,1.41499705318294) q[3];
u3(0.795117552690987,-2.64923177988716,0.818561387795197) q[5];
u3(2.43547935001531,0.351288108012488,0.999066936678435) q[1];
u3(1.72178918214473,-2.01484976454774,-1.64590677138162) q[4];
cx q[4],q[1];
u1(1.83403683723571) q[1];
u3(-3.37712301918807,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.40913984424830,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.69403829063242,-1.52026185278379,2.87040946973320) q[1];
u3(1.73269365944654,0.837722692788252,4.38659845758357) q[4];
u3(1.61861152359942,-1.52837155405209,0.481117302286469) q[2];
u3(2.27342311291676,-1.50841001676405,1.41755427430755) q[0];
cx q[0],q[2];
u1(4.12961430371288) q[2];
u3(-3.42969588613181,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.739660006543048,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.829004838641941,-3.14781961245141,2.71410599036277) q[2];
u3(2.40353458955398,1.95226324550697,-2.23753782094675) q[0];
u3(0.565110247681952,0.383714872933891,-0.976317297371792) q[1];
u3(1.85827345906375,-4.29771035948656,0.902238354940651) q[2];
cx q[2],q[1];
u1(2.57646369591478) q[1];
u3(-2.77032036965359,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.07452349041917,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.44703315458317,4.34781712435285,-1.58835584282938) q[1];
u3(1.14610138990507,-2.23102944481552,-2.73869345250216) q[2];
u3(1.17283952035795,0.801535085262191,-1.21498706911725) q[4];
u3(0.384821450069906,-2.85127001487290,1.23955076529996) q[0];
cx q[0],q[4];
u1(2.72528544041545) q[4];
u3(-1.79086594047231,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.04812132316790,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.764503031758146,2.85553941904281,0.977046172169609) q[4];
u3(1.89157917912396,-1.68183771185378,3.84043263597948) q[0];
u3(0.597534870034585,2.71955295521719,-1.96573641829617) q[5];
u3(0.997322478607750,-3.43171762520849,1.95718270093129) q[3];
cx q[3],q[5];
u1(1.87447712448513) q[5];
u3(-3.09370028835316,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.49055265317601,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.57883015340470,-0.155179011550779,-3.49361946551044) q[5];
u3(0.387491306550673,-3.64429693308448,1.36739171915985) q[3];
u3(1.93289303179541,1.51995774778179,0.766821198203261) q[1];
u3(0.482109040310405,-4.61972783285680,0.123817452048885) q[2];
cx q[2],q[1];
u1(0.522819809582297) q[1];
u3(-1.41770373550580,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.31834947702916,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.06970971589646,2.72648058221328,-3.44671620038183) q[1];
u3(2.32171410511142,-4.05728660402543,1.88569947948052) q[2];
u3(1.87421698043724,-3.19775766609810,2.18782613087120) q[3];
u3(1.69704737655378,-2.84002877324886,1.91554358698838) q[0];
cx q[0],q[3];
u1(2.25744382343809) q[3];
u3(-2.93077412774225,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.556591742617526,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.19739323545979,-0.444090321410083,-0.882770678730349) q[3];
u3(1.68169627566458,-4.87998414019947,-0.844364985626272) q[0];
u3(0.805898719362482,2.19482801585352,-1.21164829152129) q[4];
u3(1.55686509535567,0.668092393981982,-3.47976355896868) q[5];
cx q[5],q[4];
u1(1.04203087362296) q[4];
u3(-0.722455593286290,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.62583309792291,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.35720868722455,0.316662326145591,1.63296520206511) q[4];
u3(2.32739876255508,-0.937936242586521,4.75314422792645) q[5];
u3(1.17323749142617,-0.275145889234691,1.29417095533473) q[2];
u3(1.38808049074430,-1.47113002564687,-1.60030853310137) q[5];
cx q[5],q[2];
u1(0.814068092250723) q[2];
u3(-0.0373244210156185,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.41704301299425,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.50060351237599,-3.10202410432656,0.177798255717891) q[2];
u3(0.735296779418224,4.95063849550900,-1.00607634271612) q[5];
u3(1.20601663711749,0.724699467013817,-2.47064281327136) q[1];
u3(1.58959046308367,2.13390132593150,-4.07163727391198) q[4];
cx q[4],q[1];
u1(2.54426930953024) q[1];
u3(-2.04806566333762,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.364745397969942,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.75550779128282,-0.873411083907929,4.08879794192228) q[1];
u3(0.971504232256825,-1.21190912308279,4.61634513426624) q[4];
u3(1.80345651603791,3.48115139851053,-1.73541113053612) q[3];
u3(1.69562931419267,1.39188879872334,-2.43973320577142) q[0];
cx q[0],q[3];
u1(-1.26692724563735) q[3];
u3(0.400383115449318,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.82660933729269,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.932265183747792,-0.510016050156908,3.57028038533402) q[3];
u3(1.08640925920517,0.00773662092780447,-2.84346259447765) q[0];
u3(0.845397933383516,-1.33771052781883,1.13328590956168) q[2];
u3(0.800137027230252,-2.09491232379200,-0.215755516065953) q[3];
cx q[3],q[2];
u1(1.93414459516265) q[2];
u3(0.212405711182158,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.908822449568384,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.14129224831965,-2.57517209893249,1.77368830296468) q[2];
u3(2.10445322778122,3.82971111519878,0.398723945759857) q[3];
u3(1.51092619223072,1.22363513431627,1.16865911956598) q[0];
u3(0.965972176737346,0.175730120361487,-3.53624145715521) q[1];
cx q[1],q[0];
u1(1.90912045052417) q[0];
u3(-2.37566299330683,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.496415654843369,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.32079727334109,-2.06812047641603,-0.916319541691299) q[0];
u3(1.33377006245614,5.05735265488896,1.08223058383606) q[1];
u3(2.66348902794009,-0.190300276191822,1.01834832778417) q[5];
u3(1.82799512687600,-1.96035785318909,-2.96762061550349) q[4];
cx q[4],q[5];
u1(1.03728889846969) q[5];
u3(-0.758231351927686,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.83143647661609,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.32582419101770,-2.90833366461659,1.79979221463814) q[5];
u3(2.53845429003670,5.01838424731659,-1.02930266635250) q[4];
u3(2.16262191996779,2.22268131662427,-3.62387650845286) q[1];
u3(2.45098070542195,3.30905561245983,-2.78887046699458) q[3];
cx q[3],q[1];
u1(1.25751876697809) q[1];
u3(-1.60853974110476,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.66921092961416,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.64150981322015,3.27064848229910,-1.00960069651210) q[1];
u3(0.579023533735000,1.54780622362530,-1.21199281357269) q[3];
u3(1.70841360515653,1.67957226791250,-2.96421329761010) q[0];
u3(0.847366894349479,-2.01268526303405,1.62920665030168) q[2];
cx q[2],q[0];
u1(1.43302080920208) q[0];
u3(-0.629990490632764,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.29942174452548,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.711977699994638,1.11537489531964,0.221552693314551) q[0];
u3(1.17266128139439,-0.840463901565526,-1.80829430503601) q[2];
u3(1.45667424671731,1.51968636007891,0.356506733386215) q[4];
u3(1.03301561128945,0.295295468117025,-3.33534917844377) q[5];
cx q[5],q[4];
u1(3.14736582533751) q[4];
u3(-2.51139865429239,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.987580563683647,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.08219689999537,3.24041505105432,0.609767285018672) q[4];
u3(0.747576425360983,-4.75433573115200,1.15584548205327) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
