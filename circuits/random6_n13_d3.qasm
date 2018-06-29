OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.00243281405555,-0.988617608704430,2.39970228307906) q[4];
u3(0.904440203872517,-2.64164703853047,-1.10107963910452) q[0];
cx q[0],q[4];
u1(0.938044571095180) q[4];
u3(-3.28954948990432,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.65195203892421,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.31229046627758,-0.651764557731746,3.61151663099906) q[4];
u3(1.08648461331613,-2.01527120111727,-2.06242695451491) q[0];
u3(2.13715707037431,-0.332266288444132,2.91945105006021) q[10];
u3(1.94706167038161,-3.33691970913014,-1.88260085611142) q[1];
cx q[1],q[10];
u1(2.22783559770428) q[10];
u3(-2.95134340676850,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.12921604729049,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.48360477822184,-0.345667024918613,1.23857159170003) q[10];
u3(1.79738076655772,3.19550393892361,0.815498629230470) q[1];
u3(1.53330952234716,3.10876868419178,-1.76362276921587) q[12];
u3(0.725234420911001,0.721838984859603,0.0847330616086615) q[8];
cx q[8],q[12];
u1(2.07429626885571) q[12];
u3(-1.89981144017964,0.0,0.0) q[8];
cx q[12],q[8];
u3(-0.330867285663124,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.36779236901534,0.744071842341103,-3.85752316083366) q[12];
u3(1.19059237370580,2.91900519652956,0.963501002234169) q[8];
u3(1.95506074949799,1.92707993989117,-4.35395911236108) q[9];
u3(0.598139894445367,-0.703884625881478,2.27076579969597) q[6];
cx q[6],q[9];
u1(1.92524876864290) q[9];
u3(0.730444251850729,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.36478145682731,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.720033616174101,-0.103473100131535,0.0344675343705642) q[9];
u3(1.21021426674690,1.20371625709612,-2.44367932920967) q[6];
u3(1.39353601100537,1.78237581738737,-0.811310196663235) q[2];
u3(2.52424808194761,-0.318857208438062,-3.30311864797286) q[7];
cx q[7],q[2];
u1(1.67719951421709) q[2];
u3(-2.89867399287824,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.916874433157344,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.09710341482646,-1.04887601598493,-0.929591701463786) q[2];
u3(2.64872017591631,-1.55478173480467,-4.00124042625086) q[7];
u3(0.971586926178713,2.49388061381923,-1.71851667821613) q[11];
u3(1.37692273601823,2.03844807498249,-2.36927344815086) q[5];
cx q[5],q[11];
u1(-0.00810385943448177) q[11];
u3(-0.586426244783662,0.0,0.0) q[5];
cx q[11],q[5];
u3(2.40848962095154,0.0,0.0) q[5];
cx q[5],q[11];
u3(0.295891775984177,1.16222948804120,-3.08204594714327) q[11];
u3(2.22396227864130,0.114922484393298,2.76516615792743) q[5];
u3(1.31228986824549,1.09359213014620,-0.365756135626194) q[2];
u3(0.505013766837196,-0.175666717965091,-3.51072788452507) q[1];
cx q[1],q[2];
u1(2.25899554873968) q[2];
u3(-2.50200582036866,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.97511778999213,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.33324002190882,1.54863893137601,-2.29557516165963) q[2];
u3(0.359193516934997,1.11958572952798,-3.34261075521056) q[1];
u3(2.03422196402222,0.262270490376179,2.63760584457522) q[12];
u3(2.73873290095419,-1.93490345890043,-1.10876085201780) q[6];
cx q[6],q[12];
u1(3.25576971189129) q[12];
u3(-1.30053908034644,0.0,0.0) q[6];
cx q[12],q[6];
u3(2.57121264922286,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.18403102896532,-0.555592786584226,-0.331387852107673) q[12];
u3(0.991841891520151,-0.563645423711361,-3.62971156603020) q[6];
u3(1.56773972126557,-0.766535419904551,0.185809806165542) q[9];
u3(1.22871413364506,-3.08205400427630,-0.148643786955556) q[4];
cx q[4],q[9];
u1(1.93436848444483) q[9];
u3(-3.05721548029846,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.25884290419565,0.0,0.0) q[4];
cx q[4],q[9];
u3(2.37210556051236,1.72994182312555,0.567721079709481) q[9];
u3(2.02388340089768,1.16884056682191,1.83622742553726) q[4];
u3(1.53784528648151,1.17972360949832,0.0630520623916772) q[10];
u3(0.445579069105853,-0.158833968596754,-3.46793643187450) q[7];
cx q[7],q[10];
u1(2.98281922312314) q[10];
u3(-2.82140449013274,0.0,0.0) q[7];
cx q[10],q[7];
u3(1.35495883816283,0.0,0.0) q[7];
cx q[7],q[10];
u3(2.89772112897428,-2.82086257038762,-0.564789191962842) q[10];
u3(2.28371396592291,-3.85868889579788,1.66297336149937) q[7];
u3(2.19433934894593,1.60262013215576,-2.87487357714930) q[11];
u3(2.44900969775707,2.64383651653085,-3.25768306794537) q[5];
cx q[5],q[11];
u1(1.59507202944414) q[11];
u3(-2.61377306282698,0.0,0.0) q[5];
cx q[11],q[5];
u3(3.54967313216180,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.02556931409036,-3.92681174743607,0.805463569841242) q[11];
u3(1.93987287774208,2.21651507329768,0.349325126032976) q[5];
u3(0.760635759379036,-1.55737445829518,1.08515262720950) q[3];
u3(0.722760491048457,-2.23793567725835,-0.494219647353781) q[0];
cx q[0],q[3];
u1(0.138994981447555) q[3];
u3(-1.36984818508497,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.36870576864664,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.12026299398191,1.40760197549362,-2.05605524021806) q[3];
u3(0.769792759035709,0.852222147827027,4.06466360513106) q[0];
u3(1.70398512786876,-1.10713487678365,0.708544571658032) q[0];
u3(1.59359548966329,-2.02071088146159,-1.58218698512440) q[9];
cx q[9],q[0];
u1(3.19677630301593) q[0];
u3(-1.14912262620765,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.37500475479437,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.74980445445670,1.31381057575175,2.59483778571403) q[0];
u3(1.54063535363144,-1.05980999924554,-4.93383709371291) q[9];
u3(1.16231037413187,-0.0382536630083274,1.14872494947732) q[12];
u3(2.11550830497833,-0.666989555268112,-2.32144244307475) q[8];
cx q[8],q[12];
u1(0.583972305827031) q[12];
u3(-1.47493922461678,0.0,0.0) q[8];
cx q[12],q[8];
u3(-0.264839289276242,0.0,0.0) q[8];
cx q[8],q[12];
u3(2.61918154719601,1.82768464511754,-0.334230625312196) q[12];
u3(2.38957680436366,-3.57027142527281,2.28250785550711) q[8];
u3(1.30544717415454,3.38099115658680,-2.35027539548849) q[10];
u3(2.39327837811227,1.46934464722425,-1.82506172178186) q[1];
cx q[1],q[10];
u1(0.145416325540605) q[10];
u3(-1.83072761713749,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.549511996336801,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.17770702383922,-1.50061472788398,-0.0875228324789894) q[10];
u3(1.40011976505527,-5.32849329845592,0.657057739281339) q[1];
u3(1.52475751957176,3.81158955058064,-2.40571771224477) q[11];
u3(0.168798973118255,-0.320675163814872,2.36235530643012) q[4];
cx q[4],q[11];
u1(3.17130078741587) q[11];
u3(-1.35741067394891,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.17187161726657,0.0,0.0) q[4];
cx q[4],q[11];
u3(0.403713489163680,-1.27064550569146,1.25655126774984) q[11];
u3(1.47372987917148,3.62151115265737,-2.00693472377529) q[4];
u3(1.99450233264139,1.48910341167670,-0.372772064264637) q[6];
u3(1.21432647847605,0.513404333249092,-3.94976447238265) q[3];
cx q[3],q[6];
u1(2.44142601062485) q[6];
u3(-2.08463628529461,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.38541627006932,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.32178471169314,-0.0604937246633306,3.77861296781885) q[6];
u3(1.94215792976026,-4.56036808758666,-0.736529465280628) q[3];
u3(1.48214986404494,-0.577681970749933,0.0651669130145691) q[7];
u3(1.30394962333969,-2.60756626115921,-0.970763343105185) q[2];
cx q[2],q[7];
u1(3.00899458191615) q[7];
u3(-1.74667469102263,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.361473868417777,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.459444452931990,-0.173958470353360,0.604910826414915) q[7];
u3(0.833059027375731,0.880763272754658,2.21347462502064) q[2];
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
