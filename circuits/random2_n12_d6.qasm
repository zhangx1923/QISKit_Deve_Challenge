OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.31924694980461,1.04808561233884,-0.964561114187344) q[9];
u3(1.71896716138120,-0.381556497144868,-2.52080870197288) q[3];
cx q[3],q[9];
u1(2.41682865975497) q[9];
u3(-1.85056701176823,0.0,0.0) q[3];
cx q[9],q[3];
u3(3.20416646108534,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.84131307826385,-2.75102377594459,3.42862292231698) q[9];
u3(1.14024657497826,-1.16741183648064,2.37468881255823) q[3];
u3(0.347295698428966,-1.87329921290526,3.05304823138768) q[5];
u3(1.00825024442661,-3.65211364554945,1.86689430934822) q[11];
cx q[11],q[5];
u1(-0.0540345586837980) q[5];
u3(-1.10635573175329,0.0,0.0) q[11];
cx q[5],q[11];
u3(2.27240626406791,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.90529679354575,-2.51988827754223,2.36066734219653) q[5];
u3(2.69879385220142,2.20935459910452,-0.796695445067908) q[11];
u3(1.93656397863093,-1.07655741483847,1.21501091876416) q[6];
u3(1.60620536521348,-2.15443426062645,-0.570385360216114) q[2];
cx q[2],q[6];
u1(1.89735738452653) q[6];
u3(-2.10873835680162,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.99014378037476,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.00181072606535,-2.67022717189854,3.58096767911269) q[6];
u3(2.38562370461225,-3.34634807211072,-0.457971596600066) q[2];
u3(1.11314497888912,-1.10672817568201,-0.0212027597996779) q[1];
u3(0.830471080033911,-3.10299728067046,1.11785724553402) q[4];
cx q[4],q[1];
u1(2.06438553043971) q[1];
u3(-1.83134405724711,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.78208065557773,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.929864243710882,0.320617686563696,-2.76782649096246) q[1];
u3(2.27164629197849,1.80685891782052,-3.92146646681807) q[4];
u3(0.659776073285565,3.86683918008287,-2.03169998037561) q[0];
u3(1.44685257780413,1.42831481891187,-0.728904514061840) q[7];
cx q[7],q[0];
u1(1.41758804992649) q[0];
u3(-0.538920622646832,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.96330378371331,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.03172078562293,3.67702521891877,-0.268042813464432) q[0];
u3(1.85973989345217,-2.61204628452500,2.35388526081536) q[7];
u3(1.02905515919430,2.40561847703451,-0.250051399677778) q[10];
u3(1.27931675649705,1.57693181301985,-1.24942939266209) q[8];
cx q[8],q[10];
u1(3.55035530348337) q[10];
u3(-3.86156057540497,0.0,0.0) q[8];
cx q[10],q[8];
u3(-0.956925149767624,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.03309074112510,1.06072345847858,-2.88392300221885) q[10];
u3(2.93880586922706,3.85847621839492,-0.316089580998044) q[8];
u3(1.80475985225077,-1.64835008466311,1.09019848205500) q[0];
u3(2.03741217512031,-3.50020453204386,0.345217496181736) q[6];
cx q[6],q[0];
u1(2.48668714196741) q[0];
u3(-2.07067822339874,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.13441769393815,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.995991411208699,-2.85616540961792,2.63917343223725) q[0];
u3(2.36703850523412,-1.81767953707519,3.36453130679186) q[6];
u3(2.44028368307541,3.13603381065590,-2.03035266131189) q[1];
u3(1.43121153117759,1.86517782273029,-2.16882067636268) q[7];
cx q[7],q[1];
u1(1.49458227853495) q[1];
u3(-0.848501234248408,0.0,0.0) q[7];
cx q[1],q[7];
u3(3.06955560371985,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.777394078054256,1.92876164223315,-0.990393137483747) q[1];
u3(1.06655812899618,-1.66186695957809,-1.72413830852675) q[7];
u3(2.49249700110617,-2.44577390611919,1.63857067379545) q[4];
u3(2.12540803146607,2.05517649366785,3.63926433039340) q[10];
cx q[10],q[4];
u1(1.32741582386828) q[4];
u3(-0.793463091954381,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.0670392498170376,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.38573832383048,2.35692158393868,1.35811896400019) q[4];
u3(1.71877211807797,2.51690167437742,-2.69862389483574) q[10];
u3(2.37231189365876,0.105139628778349,-3.24184899847657) q[5];
u3(1.46225739332781,-3.06517939685794,3.15331234359138) q[8];
cx q[8],q[5];
u1(2.33402283579675) q[5];
u3(-1.71057453213132,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.215183660237489,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.838599606773347,1.43746133956632,-3.92516875802435) q[5];
u3(0.216545812309016,4.41735730706621,1.23160781463904) q[8];
u3(1.86258235523449,-1.30505434517141,4.39328569024635) q[11];
u3(0.903622854614100,2.07868166135400,0.886152467699632) q[9];
cx q[9],q[11];
u1(1.33425695705978) q[11];
u3(-0.547632860385336,0.0,0.0) q[9];
cx q[11],q[9];
u3(2.11694769923802,0.0,0.0) q[9];
cx q[9],q[11];
u3(0.612719064458050,1.88011928939183,-2.10549973339602) q[11];
u3(1.84072173420263,2.01208832702445,-0.162630720394129) q[9];
u3(1.07222011388933,-1.56233809671186,-0.368754361571709) q[3];
u3(0.634601028281964,-2.61096767978841,-1.45356489953982) q[2];
cx q[2],q[3];
u1(2.16704803712073) q[3];
u3(0.109703062774692,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.08153295061899,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.728109646634861,1.05819471953998,-1.64496611884078) q[3];
u3(1.01703585227105,0.189905467057662,1.39606865750584) q[2];
u3(1.66483788517831,1.82589628091312,-2.49383366421906) q[11];
u3(2.04321396456250,-2.80957579869062,2.43327227710529) q[6];
cx q[6],q[11];
u1(1.57781511096115) q[11];
u3(-2.38053815005077,0.0,0.0) q[6];
cx q[11],q[6];
u3(-0.191817142436626,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.36667701419051,-0.560598938802060,2.04484223069753) q[11];
u3(2.55652075161051,0.420388193063562,-4.91604286582562) q[6];
u3(1.14130864733588,0.294638861113846,1.43789620904907) q[7];
u3(1.05986742486466,-1.17315474739896,-0.243916400556830) q[8];
cx q[8],q[7];
u1(3.11880342869270) q[7];
u3(-1.70177931571893,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.16606376922797,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.74728702580361,-1.26832556574742,-1.04938354802915) q[7];
u3(2.50107843112249,3.01488904225270,0.142544985633694) q[8];
u3(2.41574446300472,0.885563566366213,-0.765520158176234) q[0];
u3(1.67306274329892,-0.254319184444200,-3.38886047399266) q[10];
cx q[10],q[0];
u1(2.30399376367580) q[0];
u3(-1.65851576934796,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.97965959802759,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.58009830350189,3.25630061147679,0.0995803524249341) q[0];
u3(1.86685712290240,-1.26430242633350,1.09195134281809) q[10];
u3(1.45490395236180,1.50667261058258,1.00027465530336) q[1];
u3(0.760073463010133,0.642812945728801,-3.89873066185631) q[4];
cx q[4],q[1];
u1(-0.570558828052588) q[1];
u3(-1.86598760266252,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.00892357804729,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.948333995993085,0.0739345524950187,-0.361296597884492) q[1];
u3(1.14664529104853,-0.299652231148665,-3.55572153299733) q[4];
u3(1.27587911860613,-2.42412215922045,3.63100672017324) q[5];
u3(2.03865963438191,2.11760694187200,-1.80855662183218) q[9];
cx q[9],q[5];
u1(0.769424369700090) q[5];
u3(-1.30693651087102,0.0,0.0) q[9];
cx q[5],q[9];
u3(-0.214970297293988,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.67443364902506,-1.54846326428992,1.20419462051015) q[5];
u3(0.888337423079331,-2.91967746503058,0.420751617344342) q[9];
u3(0.858001027155497,-0.398023476972542,0.0707189211931001) q[2];
u3(1.06081314149082,-0.288292607398271,-1.08764188835445) q[3];
cx q[3],q[2];
u1(-0.331539456913891) q[2];
u3(-1.05593975191023,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.70842212007789,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.837658796945600,-1.48668235136984,0.733732292013620) q[2];
u3(0.939408623399962,1.77298091541916,-3.45175940336941) q[3];
u3(2.48415640185118,0.595370147907978,-0.875972548693697) q[0];
u3(1.43316008926947,0.330306139280888,-3.86928489467442) q[3];
cx q[3],q[0];
u1(1.82884697667453) q[0];
u3(-2.44028314050887,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.28537675911284,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.913008559972994,3.71131801745271,-2.21617067068886) q[0];
u3(1.80007141984755,-2.34316768454603,3.29160540610917) q[3];
u3(1.81672331287881,1.83792599890551,-0.147494443681509) q[5];
u3(2.41594352772708,0.609934597444471,-3.71368532690442) q[1];
cx q[1],q[5];
u1(2.94105957548003) q[5];
u3(-2.12273653087910,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.54631765509233,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.663147904842232,2.68294026511106,0.777721299466040) q[5];
u3(2.22841439210558,-3.17336963411933,-1.92202317526053) q[1];
u3(2.34697090334598,-2.77172516080890,0.788484893641803) q[2];
u3(2.27395060863780,2.08949539400429,4.14516882140861) q[8];
cx q[8],q[2];
u1(2.74323018098204) q[2];
u3(-1.87102129229540,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.68202154190104,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.50061580269760,1.20816925159870,-3.52073491752146) q[2];
u3(1.56817881184591,2.80588537161592,-1.18931594186427) q[8];
u3(0.442157374900782,-1.77598634540343,0.698368064105340) q[9];
u3(0.278812307891273,-2.17509359010605,0.132288430377719) q[10];
cx q[10],q[9];
u1(0.695392561064712) q[9];
u3(-1.23994207690032,0.0,0.0) q[10];
cx q[9],q[10];
u3(3.15914612104914,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.29807197208828,-2.43589750085218,1.42801386380246) q[9];
u3(0.392852736677948,-1.43844671193147,-2.48442204950395) q[10];
u3(2.45948587667323,2.81349710222501,-0.0871324589428357) q[4];
u3(2.90900675042777,1.00821870258653,-2.67345968789761) q[6];
cx q[6],q[4];
u1(3.48375116308897) q[4];
u3(-1.42594940614559,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.13307606003991,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.05271647999176,-0.346946285030727,-0.553622471528418) q[4];
u3(0.431530345112449,-0.728417846214103,-4.35318935863618) q[6];
u3(2.35481546435510,1.17872914911610,1.16228735327325) q[7];
u3(1.05539024033706,-4.10062163513704,-0.558299553396068) q[11];
cx q[11],q[7];
u1(-0.425394717726213) q[7];
u3(-1.60856750194646,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.42411834017885,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.63530364826898,-3.27548753113693,2.56537810426337) q[7];
u3(1.82899550030786,-0.0175156959198048,-4.04963257762506) q[11];
u3(1.30806457336768,-0.319637609316919,1.74787463460237) q[5];
u3(0.596782575679787,-0.695884633297478,-1.69520637566880) q[4];
cx q[4],q[5];
u1(0.598010528086464) q[5];
u3(-1.43769676422032,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.197340056627866,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.24381854006643,1.19342992761169,-2.67606060322995) q[5];
u3(0.884099894418625,2.53015445591152,2.00302476965765) q[4];
u3(2.34901086717511,-4.22544734448748,1.36989352498782) q[7];
u3(0.608278659737494,1.91935421294796,0.232760822266442) q[2];
cx q[2],q[7];
u1(-0.0744527579943821) q[7];
u3(-1.53921590911600,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.81854418093557,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.46534619028088,-3.04365591662724,-0.896457815513039) q[7];
u3(0.439351699164391,-1.57901669907535,3.57144156029953) q[2];
u3(2.56285529020972,-0.585647658509045,2.21703397937181) q[9];
u3(2.12633254620913,-3.36536363506113,-2.31856312042072) q[3];
cx q[3],q[9];
u1(3.10773603301515) q[9];
u3(-1.63200076041430,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.102700007714012,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.94465733243539,3.12153457186581,-0.109501650109514) q[9];
u3(0.319950102265543,0.373778269073862,-1.50984142186411) q[3];
u3(0.503768788117196,-2.23510617277262,1.78696482298521) q[8];
u3(0.262477035501308,1.62009634799456,-3.74613736148082) q[10];
cx q[10],q[8];
u1(0.291096246582503) q[8];
u3(-1.29790251971563,0.0,0.0) q[10];
cx q[8],q[10];
u3(2.50421412493417,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.952583224875436,1.83558871234355,-0.0291318771498296) q[8];
u3(1.56925511182093,1.72573024174940,-2.20955981152578) q[10];
u3(1.86519516325774,-0.454044301560502,1.90969807319515) q[6];
u3(1.57844468287984,-1.27964984989011,-0.642227250500751) q[1];
cx q[1],q[6];
u1(0.0703467628289760) q[6];
u3(-0.563668458258610,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.77223742416659,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.663749204516813,-0.447417322023162,-0.436317837232570) q[6];
u3(1.24188579709976,0.576056246852009,-5.28782898727084) q[1];
u3(0.704050021166812,-1.45633384025527,-0.0191579870116528) q[11];
u3(1.65030368724140,1.43187152181655,-4.44240495698384) q[0];
cx q[0],q[11];
u1(1.75372970204802) q[11];
u3(0.415941820139179,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.762569096669058,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.519112307168692,-0.385916407260735,1.14367596226042) q[11];
u3(1.98989151674290,-0.423783328370073,3.68772564729578) q[0];
u3(1.61219121181162,2.58566659964622,-3.31776067312593) q[5];
u3(2.70525967504078,3.60661428498178,-2.48308290605121) q[4];
cx q[4],q[5];
u1(1.08495863452989) q[5];
u3(-0.420880319235737,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.88717392620099,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.18287174351872,2.38596812166535,-2.23982140730458) q[5];
u3(1.15563889936738,2.92267335271541,-0.350990147711026) q[4];
u3(2.35363669072091,3.28018056758080,-2.99187839511609) q[9];
u3(1.68019754091193,2.04536322084276,-1.76233816417755) q[1];
cx q[1],q[9];
u1(3.02181720519414) q[9];
u3(-1.18963977804899,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.43636911581044,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.19921314336120,3.40465834879573,-0.485897187129669) q[9];
u3(0.733694615104945,-2.51190852957814,-3.41244323952104) q[1];
u3(2.67301233803184,1.86977694960636,-3.71955878097280) q[10];
u3(0.733691325047549,-1.57828414672838,3.67086569531340) q[7];
cx q[7],q[10];
u1(3.09042279769103) q[10];
u3(-1.76430494797618,0.0,0.0) q[7];
cx q[10],q[7];
u3(0.749183739416261,0.0,0.0) q[7];
cx q[7],q[10];
u3(2.41172630780260,1.23523876623685,0.424204442037843) q[10];
u3(1.68275888223017,2.34107039103577,-1.53326497122409) q[7];
u3(1.38475567321082,1.66608749319492,0.485567871308298) q[8];
u3(0.790182950345474,-0.824440929791750,-1.98720819031891) q[0];
cx q[0],q[8];
u1(0.417123133203450) q[8];
u3(-0.684091686940607,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.11866994678252,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.51746397279378,-0.428376406461291,2.75598737505565) q[8];
u3(2.09204690450131,-0.713817846067784,-4.36978260094051) q[0];
u3(2.49144926210087,1.59842750449653,-3.04231698242451) q[3];
u3(1.75213777444123,-2.94481458768863,2.67407396494201) q[6];
cx q[6],q[3];
u1(0.667143644202976) q[3];
u3(0.0348973579021610,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.86178426932576,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.49139824399770,2.03864637535263,-0.0100053690190975) q[3];
u3(2.95824547180625,-1.01771031143954,-0.312268546946747) q[6];
u3(1.87746319388695,0.0230396514304118,1.73039860324561) q[11];
u3(1.52743147858362,-0.532879992321567,-1.25397399375867) q[2];
cx q[2],q[11];
u1(3.03880558711544) q[11];
u3(-1.29553503505726,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.72934149614791,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.987627150739721,1.08172251806961,-2.66964499701953) q[11];
u3(1.25126985974606,-4.57770588199355,-0.507437995515765) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
