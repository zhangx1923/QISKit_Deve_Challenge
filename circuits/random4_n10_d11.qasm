OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.39220629294875,0.280230869317779,0.858174296046971) q[7];
u3(2.18148010894745,-1.98016337044444,-1.55568954781279) q[3];
cx q[3],q[7];
u1(2.47558515186801) q[7];
u3(-1.30392370412549,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.24809070433201,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.63469566944136,2.36667265730085,-0.911973762371320) q[7];
u3(1.41408044335697,-5.72913302291493,-0.543613763290320) q[3];
u3(2.61269161758309,-1.25475653426905,4.18593134681573) q[0];
u3(1.55235336355331,1.68605836233983,1.75236784285498) q[8];
cx q[8],q[0];
u1(1.10772579022956) q[0];
u3(-3.58939581770109,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.57556049661311,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.51024431486278,-0.786690590092083,-0.634377198826614) q[0];
u3(1.40121709601499,4.88187882048386,-0.611776796159258) q[8];
u3(1.84866055536440,-2.32374491464297,-0.303689555834041) q[6];
u3(0.898627451513718,1.62334474522451,4.63570388141994) q[1];
cx q[1],q[6];
u1(1.48362472290842) q[6];
u3(-3.29963485073089,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.17045950257234,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.58820012688684,-1.55085346646861,-0.935417809930035) q[6];
u3(0.367265041015587,1.00585229310184,3.65546556458572) q[1];
u3(1.98527729329439,3.30816770641247,-0.793634017982358) q[5];
u3(1.32502364057498,1.73266342431698,-0.896364533918406) q[9];
cx q[9],q[5];
u1(1.35720123226557) q[5];
u3(-0.583477777924679,0.0,0.0) q[9];
cx q[5],q[9];
u3(-0.326333982061923,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.50879874393982,-0.395745347916387,2.88143898327801) q[5];
u3(1.56020399186676,1.90717762353972,2.66737255939428) q[9];
u3(2.01352663963494,-0.382385840198438,-0.375429575785644) q[2];
u3(0.537138423988488,-2.88728900661775,-1.11840029338734) q[4];
cx q[4],q[2];
u1(3.64617797311478) q[2];
u3(-1.21234562037627,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.25327605776454,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.53153353603539,-3.77810123079151,0.338398130501600) q[2];
u3(2.18921769045934,-4.71638035305467,1.25358562143049) q[4];
u3(2.55882521067920,3.55951206057301,-2.37924217278006) q[3];
u3(1.49191159395929,-0.576814457491384,2.32134243010042) q[5];
cx q[5],q[3];
u1(0.952132236235249) q[3];
u3(-3.24925623864245,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.55695470710510,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.651867260397801,-1.21905927414261,-1.19046772232422) q[3];
u3(2.69718493204711,-2.39798306179115,0.432404013617397) q[5];
u3(2.12789111686139,0.614435022158510,0.898856426338539) q[6];
u3(1.59282755074370,-1.76190822584957,-2.00158903573063) q[2];
cx q[2],q[6];
u1(1.67398842455378) q[6];
u3(-2.19779640651516,0.0,0.0) q[2];
cx q[6],q[2];
u3(3.69831125767142,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.43709889275070,-1.16392324838600,2.18151069373466) q[6];
u3(1.07666717552171,3.14446430502200,0.0226691751975951) q[2];
u3(1.74930343964173,-0.193258376390663,0.726065524195996) q[7];
u3(1.15155576403751,-2.90213300335320,-1.68565203495574) q[1];
cx q[1],q[7];
u1(0.958873878391945) q[7];
u3(-0.385233156564234,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.58997398034590,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.78395776439095,-1.35489913497991,1.27432092290022) q[7];
u3(1.97732237471610,-3.48843834207134,-0.410315787446550) q[1];
u3(2.48552665748682,0.650423927294891,1.18227761933805) q[8];
u3(1.15821382899974,-1.16679247446165,-3.67660955221894) q[9];
cx q[9],q[8];
u1(1.34501946514649) q[8];
u3(-0.373353558284640,0.0,0.0) q[9];
cx q[8],q[9];
u3(2.48482601001159,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.11868065670578,1.94102742064319,1.26916642155250) q[8];
u3(1.30680125062362,-0.388248427805954,2.02413160692181) q[9];
u3(2.34780115914331,3.51742688439794,-2.60874527775553) q[0];
u3(0.994297710129388,3.15373428754544,-2.33121288479105) q[4];
cx q[4],q[0];
u1(0.570056684757920) q[0];
u3(-1.39402030347568,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.85679118527766,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.85551208617965,-2.18597451734528,0.251435315583449) q[0];
u3(2.60527959906461,-2.78576437939435,3.34834301600762) q[4];
u3(1.21770494327339,-0.0294231931292628,-2.27380991360803) q[9];
u3(2.14841644229679,0.219857487269916,-4.90680316796086) q[0];
cx q[0],q[9];
u1(-0.0979282492075153) q[9];
u3(-2.38888316367932,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.23364360155718,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.919224504365598,-1.65493291126782,1.16733238914841) q[9];
u3(0.740348713736694,0.901738706448083,3.94996322756551) q[0];
u3(2.46234360625663,0.419139500913065,-0.213623207889049) q[7];
u3(0.985965407086573,-0.173090292807456,-4.87162548819042) q[5];
cx q[5],q[7];
u1(1.42852925813323) q[7];
u3(-3.78137868956345,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.08672456843814,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.99325008490592,2.16070169081312,-2.43406662301381) q[7];
u3(0.624180675633214,-2.32388800090799,-0.249060514817197) q[5];
u3(2.08135034654328,-0.952443877425915,1.54922334065685) q[4];
u3(2.07436750575903,-1.56417733195902,-0.976617488086691) q[6];
cx q[6],q[4];
u1(1.42705315898972) q[4];
u3(-3.82795591853439,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.19441423727529,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.06173013605031,1.89550214833603,-2.36322661934523) q[4];
u3(0.765636546192111,-0.213964751135693,-5.38277254500275) q[6];
u3(2.62500460447609,0.0534158875560096,-2.44017059251772) q[8];
u3(1.68332744454588,3.51040766716949,-0.912239952126394) q[1];
cx q[1],q[8];
u1(1.32249122845566) q[8];
u3(-0.195985688009570,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.23018413655662,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.35986871577329,2.98923660454931,-0.848396554185468) q[8];
u3(1.86228058448820,0.726763469345211,4.59109552241742) q[1];
u3(1.72191227604634,3.29368002628984,-1.28427275133399) q[3];
u3(2.16824518851581,1.66564672583218,-0.799838572518127) q[2];
cx q[2],q[3];
u1(2.47674974656614) q[3];
u3(0.242271201086397,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.68405364583124,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.33815912109574,-1.76782626404118,3.90025227288691) q[3];
u3(0.463383403109934,1.34300808589088,-0.227187328155060) q[2];
u3(2.44734882516488,2.06847666635559,-1.77490501802195) q[1];
u3(2.15990073365727,2.37174179957545,-3.09807338324137) q[5];
cx q[5],q[1];
u1(3.32384349938346) q[1];
u3(-4.38424811669286,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.120545968698258,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.05783762289317,-1.10375076569175,1.51151029243060) q[1];
u3(2.44047407409708,-2.29877590743726,1.71751727794365) q[5];
u3(1.82317916662896,4.69079362296612,-1.58945945542170) q[3];
u3(0.664524181471083,1.70649933704870,0.00343527531614829) q[4];
cx q[4],q[3];
u1(-0.000953523591412209) q[3];
u3(-1.86504763502904,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.374159386290904,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.484283902434012,-0.700539739993712,2.51702294212223) q[3];
u3(2.33008784135633,3.47184690160645,0.921663207399512) q[4];
u3(1.46146066250260,1.83049620608080,0.508944471716595) q[8];
u3(1.96352166439453,-0.165097917537123,-4.07749603451685) q[2];
cx q[2],q[8];
u1(1.69815148172980) q[8];
u3(-2.39277513504412,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.49925580102722,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.762681519643425,-4.58195038604135,1.51387745129794) q[8];
u3(1.53859444521167,0.897224006064930,0.527805593140094) q[2];
u3(2.33797445913403,1.55116840513704,-3.74248797548376) q[7];
u3(1.14458667015361,-2.56864371246006,3.52999735525377) q[0];
cx q[0],q[7];
u1(-0.992848156237722) q[7];
u3(0.0811056794896143,0.0,0.0) q[0];
cx q[7],q[0];
u3(3.67090825546199,0.0,0.0) q[0];
cx q[0],q[7];
u3(3.05686417069042,2.97638972422300,-2.32357372888948) q[7];
u3(1.11920803147726,-2.41519964796734,-2.48387002330547) q[0];
u3(0.598006886914498,-0.192837502060574,-1.67696363276416) q[9];
u3(0.864180762326761,0.717474505301043,-4.86240885689794) q[6];
cx q[6],q[9];
u1(0.120904883181297) q[9];
u3(-0.973877329600620,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.53725783759079,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.12733571758421,-1.58877157814280,2.94615122723080) q[9];
u3(1.04874708609840,-2.23715129586626,3.58365298632517) q[6];
u3(1.12434739701825,1.54854020213080,-2.76360945816481) q[9];
u3(2.42042432771199,2.48768761362579,-3.74087580868200) q[5];
cx q[5],q[9];
u1(-0.511747026073842) q[9];
u3(0.220019092720653,0.0,0.0) q[5];
cx q[9],q[5];
u3(4.39621240617500,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.235721797788397,2.05513103417684,-2.41923325165446) q[9];
u3(0.753117435680582,3.19833517211874,-1.06233937114208) q[5];
u3(2.34377614501647,-1.32667954417327,0.886585829280051) q[7];
u3(1.94726479945125,-1.39389200581117,0.442500473512864) q[3];
cx q[3],q[7];
u1(1.22198478229232) q[7];
u3(-0.439271411362780,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.37396741741077,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.37471653790062,-3.36784626443359,1.98569455467829) q[7];
u3(2.50553531305565,-3.42817348872345,0.641612689787920) q[3];
u3(1.77066411454342,1.60764622861255,-3.08522691333561) q[2];
u3(1.12343590898472,-2.19209365781260,2.85786730111169) q[0];
cx q[0],q[2];
u1(1.29328050380354) q[2];
u3(-2.98897282992556,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.13538764377034,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.95896042730905,-4.70759577100121,1.12890997040139) q[2];
u3(2.02170257383507,3.77689735615015,-1.66003503848232) q[0];
u3(1.04535976838354,0.545222646486987,1.64611761583997) q[8];
u3(0.723246284417724,-1.86107707393961,-2.44055136876998) q[6];
cx q[6],q[8];
u1(1.85441176028624) q[8];
u3(-0.172260614837288,0.0,0.0) q[6];
cx q[8],q[6];
u3(0.945743529417595,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.56613358830388,-1.93188370624636,0.788454083009323) q[8];
u3(2.09221718842908,-1.64718375600861,3.28077241535310) q[6];
u3(2.43040228462541,-0.259561630054554,3.17228114790407) q[4];
u3(2.70741788763334,0.124815140851711,2.67770556238476) q[1];
cx q[1],q[4];
u1(3.96238879803836) q[4];
u3(-3.61589318143902,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.972190392840898,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.40314536281651,0.471491138922495,0.779318870656117) q[4];
u3(1.54561522542476,-0.877757302529614,-2.18605097384743) q[1];
u3(1.85448171978997,-2.33686371529448,-0.284959270614801) q[2];
u3(1.10883164554748,-3.56039887246046,1.15991979893162) q[9];
cx q[9],q[2];
u1(1.91315596973736) q[2];
u3(-2.97844836126731,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.622151773845013,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.56438321064466,2.69785161851781,-0.332664263697833) q[2];
u3(0.947572664038586,2.05128506965876,-2.16093406759232) q[9];
u3(2.39705672707164,3.01962380118372,-0.816565481443998) q[3];
u3(1.82471717952480,1.60264945707948,-1.44367334258929) q[5];
cx q[5],q[3];
u1(1.42030544216730) q[3];
u3(-0.678758019805762,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.000138637241245521,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.33898414526516,1.99413466450792,-3.59448471101787) q[3];
u3(0.236353242513693,-0.411336535193937,3.11394376604993) q[5];
u3(2.41586655981366,-3.13403338874960,3.12317192802776) q[0];
u3(1.11153616561726,0.0349920102840597,2.05153902426517) q[6];
cx q[6],q[0];
u1(0.783198535537719) q[0];
u3(-0.352950327257922,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.29731736783942,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.30359433275306,-1.75646505764651,2.46318424481613) q[0];
u3(2.15137082920121,2.40479665948284,1.95985874316682) q[6];
u3(2.23559134730550,-2.44541898551221,0.513203005773957) q[7];
u3(1.79913902441430,-2.38941474993168,0.622551354276094) q[1];
cx q[1],q[7];
u1(0.338812454833711) q[7];
u3(-1.32689901553520,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.92914547932546,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.253908787016290,1.47083089838008,-3.52305616628520) q[7];
u3(2.49294703700540,0.0452075292352386,-1.12438990484927) q[1];
u3(1.66802877220110,0.392286833182587,-1.30169323216249) q[8];
u3(2.34454068430983,-4.25574822215820,1.49298133012499) q[4];
cx q[4],q[8];
u1(2.59130814594612) q[8];
u3(-2.41072166856112,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.22550549510565,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.61681462240370,-1.72898255913958,3.64633057865778) q[8];
u3(1.56443862989663,0.427727551280844,3.95035838248266) q[4];
u3(1.94349723816015,1.47550150807611,0.760732469625501) q[8];
u3(0.221636287406098,-3.12416923051830,-1.06690582788607) q[9];
cx q[9],q[8];
u1(0.0670946478355683) q[8];
u3(-2.55911020027871,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.24921068749319,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.63217663781570,1.61648474782036,-4.65307354007832) q[8];
u3(1.86347634688986,-0.220255011103468,2.89840460647318) q[9];
u3(1.37327236763067,-2.33480875545410,1.22059882857636) q[7];
u3(1.47219757321806,-3.41364656007262,-0.368148690905908) q[3];
cx q[3],q[7];
u1(2.83684341370995) q[7];
u3(-1.65468594308627,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.27033781623883,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.715001587456997,2.84355937539416,-0.251355125745707) q[7];
u3(2.15234727411349,0.0398980744231059,-4.78951166144864) q[3];
u3(2.82896336156606,3.42099877062558,-1.66766677299547) q[6];
u3(1.71810734614591,2.08711577876150,-1.63804547995770) q[1];
cx q[1],q[6];
u1(3.66838170196543) q[6];
u3(-4.52222673694764,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.441268088320319,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.74991915827871,1.40985089647640,2.51332904042428) q[6];
u3(1.89939030242377,5.13378300854953,1.05072257551931) q[1];
u3(1.24823693239078,0.249361687749572,-1.95129593405947) q[0];
u3(2.32598703684929,0.871127711639234,-3.84442549526381) q[5];
cx q[5],q[0];
u1(1.49804836940286) q[0];
u3(-0.601604508446141,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.12985900248482,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.85256960930139,-1.14615729524669,2.34452915974335) q[0];
u3(0.255965709928214,4.36780929600897,1.14017319492387) q[5];
u3(1.31897608561519,3.52197854159290,-0.604268608399668) q[2];
u3(1.30419165248149,2.05732450752215,-2.10112970771236) q[4];
cx q[4],q[2];
u1(0.00155685261526761) q[2];
u3(-2.41211005951681,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.15584172032497,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.23186471605887,-1.97278992126651,3.12475985031165) q[2];
u3(0.653329743435476,-0.800688050785844,-2.68615354532659) q[4];
u3(0.473536957855895,2.67249450360148,-2.43039467339099) q[6];
u3(0.655223131617199,0.732684105107534,-1.81203037431624) q[9];
cx q[9],q[6];
u1(1.47423589033363) q[6];
u3(0.125318255374543,0.0,0.0) q[9];
cx q[6],q[9];
u3(1.07514617996693,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.60691404683885,0.525152900754443,-1.24919572129718) q[6];
u3(1.79252158666986,-0.775060336694444,-1.74759570337272) q[9];
u3(1.23794184838100,2.05500425505301,-4.13045595941417) q[0];
u3(1.71368807855646,2.34465179780015,-2.76724845561985) q[8];
cx q[8],q[0];
u1(2.49889595316762) q[0];
u3(-1.55323813058868,0.0,0.0) q[8];
cx q[0],q[8];
u3(3.12953945141253,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.63490461114065,5.21376652686831,-1.05705155577506) q[0];
u3(1.99719198002883,2.69016957830892,-3.07381450575915) q[8];
u3(0.966811460741262,-1.64860681266598,0.190168300294952) q[7];
u3(1.49263212184266,-3.33240224097508,-1.41734423132148) q[4];
cx q[4],q[7];
u1(1.13740160769715) q[7];
u3(-0.355709018412120,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.62191799559829,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.92070501076694,0.179072361239410,0.795869395578789) q[7];
u3(2.64818442672372,2.94086392087691,-3.12313240993298) q[4];
u3(0.821741315546319,0.750096990267548,1.22214546472602) q[3];
u3(1.35278035239059,-1.93697896357965,-1.04018001314924) q[5];
cx q[5],q[3];
u1(2.55576796098876) q[3];
u3(-1.93905576440962,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.701582915191453,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.74845634027429,-1.76013635484381,1.97487736067353) q[3];
u3(0.461933455484718,-1.02236519039786,-0.882447641857643) q[5];
u3(0.636372383918721,3.01230063895112,-2.93971424292906) q[2];
u3(1.18964285013163,-0.292505915513273,-2.01120345644259) q[1];
cx q[1],q[2];
u1(0.805359616079480) q[2];
u3(-1.51420158588190,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.64029955088445,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.02536479490950,4.38381454423801,-1.57740617809058) q[2];
u3(2.53817761439414,2.91760433782030,-1.36337477315648) q[1];
u3(1.78818844393775,1.36223956635744,0.854064113118052) q[5];
u3(1.51416958509062,-1.80447549292353,-0.975957236601273) q[3];
cx q[3],q[5];
u1(1.31105237574503) q[5];
u3(0.233179163310530,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.51469581210738,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.69556718068908,-0.346487006421097,4.28296237229091) q[5];
u3(2.11261309911059,-1.55579552252489,1.95376702744493) q[3];
u3(1.28464942053131,-0.176691375459365,2.09809375860838) q[7];
u3(1.23374415060409,-1.56551512961550,-1.17714712497024) q[0];
cx q[0],q[7];
u1(2.24666368197711) q[7];
u3(-1.61545278217334,0.0,0.0) q[0];
cx q[7],q[0];
u3(-0.224766123008260,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.865260538688610,-1.59873172920475,3.56249061841016) q[7];
u3(1.13811214104805,0.752523500621641,-0.919092322148187) q[0];
u3(2.48212561556532,-4.09386798449138,1.59981075312154) q[9];
u3(1.63281617354889,0.429461721627439,3.67219353406277) q[1];
cx q[1],q[9];
u1(1.86384457108543) q[9];
u3(0.232381060372230,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.681952746563951,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.57750439295524,0.216583355420569,1.70333021101349) q[9];
u3(1.51951599066322,-2.66325693374387,1.83333041342137) q[1];
u3(2.39478557302147,-0.557611998206232,2.00031781453337) q[8];
u3(1.81280450813118,-0.608815965465191,-1.47390273008396) q[2];
cx q[2],q[8];
u1(1.78358290811917) q[8];
u3(-3.07610287653326,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.653678313053145,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.38624871636046,-0.747043130429453,-1.38630205996476) q[8];
u3(2.73010121733748,-1.07131894633474,-1.54711089294503) q[2];
u3(2.65264594265136,2.83528860769847,-0.732404497513882) q[4];
u3(2.77332247688629,1.57910717451465,-4.27844133029099) q[6];
cx q[6],q[4];
u1(0.0668009701466286) q[4];
u3(-0.762655203700935,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.98929501295104,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.28127614683165,3.45949607474433,-2.34770359810029) q[4];
u3(1.94827541254407,-1.42831898835125,-4.76733504851251) q[6];
u3(2.49344182897777,0.133687460429742,-0.502511175843656) q[8];
u3(1.27949953916578,0.120838730970412,-3.95084987700531) q[4];
cx q[4],q[8];
u1(1.72952826736508) q[8];
u3(0.465458344041839,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.00377415473580,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.13042374026538,-1.10468967714965,1.03137428219962) q[8];
u3(0.364098064231423,1.29845029840309,1.45978988369916) q[4];
u3(0.329043586622454,0.639692499925906,-1.72263702104613) q[1];
u3(1.35773205829444,-3.83036095101771,1.51773881165200) q[6];
cx q[6],q[1];
u1(2.62516989212987) q[1];
u3(0.224023439670109,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.39319234231286,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.21811545486339,-0.375424521520455,0.119603184358391) q[1];
u3(0.426195493914201,3.77425056760517,-0.597164074789108) q[6];
u3(2.18398934615753,2.26505010978608,-3.61364423762122) q[2];
u3(2.60718201028849,-3.15586580416502,2.64059796406249) q[3];
cx q[3],q[2];
u1(-0.0112705070127850) q[2];
u3(-1.50504006610194,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.96036650208466,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.24495465323561,-2.54249718082681,1.57311041646241) q[2];
u3(0.711543255900658,-3.78524630823251,0.623210761588248) q[3];
u3(0.599761907992881,-0.931006763528926,-0.0295531660885577) q[7];
u3(1.11376170473444,-3.55110562714167,1.01153847359857) q[9];
cx q[9],q[7];
u1(0.621646626364382) q[7];
u3(-1.21572285319345,0.0,0.0) q[9];
cx q[7],q[9];
u3(-0.219955631380821,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.89922645655342,-0.845994749022572,2.57813299901464) q[7];
u3(2.70795785984186,-3.68205654724789,-2.18445109379405) q[9];
u3(0.782339869209879,-0.466989162469116,-0.444439018828819) q[5];
u3(1.21741638489345,-3.07764235374554,1.41120943949015) q[0];
cx q[0],q[5];
u1(4.16947517694015) q[5];
u3(-3.22591482105529,0.0,0.0) q[0];
cx q[5],q[0];
u3(-0.519285907337697,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.26781477356561,1.76698853206333,-0.987848197522193) q[5];
u3(2.69104244319171,3.18851000484621,0.280087441236809) q[0];
u3(1.78665630601681,2.95341621754468,-2.42175971885963) q[4];
u3(1.53131680130484,2.40304028739539,-1.89555882195444) q[6];
cx q[6],q[4];
u1(4.63499942315675) q[4];
u3(-3.64517825725950,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.429353196459853,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.95891527323152,-3.96789376778676,1.98755044397722) q[4];
u3(1.96227658913661,-4.97843042784694,0.382654287972759) q[6];
u3(3.03267209532458,-1.66201592423307,4.22908173156269) q[1];
u3(1.31013940036008,3.21821697983002,-2.58159651234872) q[9];
cx q[9],q[1];
u1(0.892028391836071) q[1];
u3(-0.135913849691938,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.76523233355352,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.643542928699200,-1.59977937087700,1.20992042860105) q[1];
u3(0.927505900415261,-1.53610002393752,4.53816772335234) q[9];
u3(0.689270765427319,-1.11149982651247,-1.04128521301818) q[3];
u3(1.33457468251533,1.17574413231398,-4.93198870650086) q[0];
cx q[0],q[3];
u1(0.890376522904778) q[3];
u3(-3.50244949491547,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.56878519953472,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.16533831486442,1.84220933873378,-1.37177353780911) q[3];
u3(2.23151203320346,1.46417316854538,-2.00964595092987) q[0];
u3(2.16175557373416,-1.16122869295392,-1.95990766194916) q[2];
u3(1.61965411125469,1.33610456210053,-4.83507353331113) q[5];
cx q[5],q[2];
u1(1.23367695399365) q[2];
u3(-0.957723682916478,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.551353518932657,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.01831557509159,1.89302201851866,2.71163798367844) q[2];
u3(1.61033229705215,2.51083081837648,-1.37801453937286) q[5];
u3(2.40081401639432,2.99504430608472,0.00986620630157531) q[8];
u3(2.88663951173785,3.95289878065894,0.429381094721808) q[7];
cx q[7],q[8];
u1(2.41970122475972) q[8];
u3(-3.08636075884646,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.02181625735502,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.44401948129391,1.74466117219837,0.552240460668776) q[8];
u3(2.36516832605804,-0.649622843928906,-4.26251093587701) q[7];
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
