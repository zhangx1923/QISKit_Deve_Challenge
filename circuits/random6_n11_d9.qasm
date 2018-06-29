OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.47669982881068,1.37346023201904,-0.0788960609840607) q[3];
u3(1.61313065489892,0.603547580613885,-4.19710376981301) q[1];
cx q[1],q[3];
u1(1.76556162503335) q[3];
u3(-2.13867522849526,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.46307304281282,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.46774809519100,1.43643910847661,-2.89381816929602) q[3];
u3(1.73766457212080,3.53263427815081,2.33660997347165) q[1];
u3(0.328658334164166,3.09060457233858,-2.32031522427821) q[0];
u3(0.662750958151757,-3.63680850675616,1.58786153576341) q[5];
cx q[5],q[0];
u1(1.88245975883769) q[0];
u3(-1.62711142907534,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.93089865318002,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.09284820902471,-1.85984302328555,3.72685319941110) q[0];
u3(2.77379826307146,1.28760642660164,-1.37973552073894) q[5];
u3(1.17543748774963,3.34622819798325,-0.397574496312338) q[8];
u3(0.863633598665231,1.68682274452535,-1.91794685999073) q[6];
cx q[6],q[8];
u1(3.38385722419517) q[8];
u3(-0.999049839517037,0.0,0.0) q[6];
cx q[8],q[6];
u3(2.05647708457530,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.19053178757571,0.757863687056889,-3.31850968894397) q[8];
u3(1.76079674344071,-0.173965746957147,3.73355950639031) q[6];
u3(1.77885287510528,-0.855992272855612,-1.49415807777225) q[4];
u3(1.51426715994548,1.27048494786549,-4.69604896683088) q[7];
cx q[7],q[4];
u1(3.04157619575898) q[4];
u3(-1.39817168041445,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.950518145387854,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.83434264345611,1.64688220811480,-3.01952791569374) q[4];
u3(0.913759490599166,-2.97143328871616,1.06344636224036) q[7];
u3(0.514250678533388,1.89876147910385,-2.24339319538450) q[2];
u3(0.306884293766161,-1.74496383496038,-0.306313686454843) q[9];
cx q[9],q[2];
u1(2.28150580099873) q[2];
u3(-3.24446225197020,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.36346127197249,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.66716832443411,1.05257729770026,0.423085014998736) q[2];
u3(1.15593230398104,-1.39561473244315,1.08172924191421) q[9];
u3(0.506876517500145,-0.765705892325157,0.956955496783557) q[0];
u3(1.03829679038332,-1.51826427606237,-0.962624353979779) q[4];
cx q[4],q[0];
u1(2.93001405941917) q[0];
u3(-4.21039674331120,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.376542276297253,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.92407497282449,3.12781281741682,-3.09780482253230) q[0];
u3(0.318729662031343,-4.71479589624199,-1.05534303634814) q[4];
u3(1.23976531352120,0.138991506434611,1.36425170614441) q[1];
u3(2.34027914765316,-1.61252216665289,-2.54822949957955) q[9];
cx q[9],q[1];
u1(0.101338885736959) q[1];
u3(-0.573434750372756,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.33737933021785,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.83578461887314,1.83221039218098,-2.97348375365542) q[1];
u3(3.00699646495360,2.37311338578393,-1.85824278959744) q[9];
u3(0.888744069809274,0.701756124181232,-1.79660205619808) q[6];
u3(0.488647205691489,-3.25846022346869,1.44153370364629) q[5];
cx q[5],q[6];
u1(1.59099162842506) q[6];
u3(-3.24367164839507,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.521806022435931,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.707310031518622,-1.00952459145365,3.00606940492192) q[6];
u3(2.18145692544029,2.56757143350897,-1.52617004768664) q[5];
u3(2.11363769447703,-0.962329060185040,-0.809958058629400) q[10];
u3(1.66614828205334,-2.63675599473781,0.255641733025269) q[8];
cx q[8],q[10];
u1(1.46452399505037) q[10];
u3(-3.12344321182636,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.62402526511382,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.42802599579410,-1.08196504521665,5.14061786477953) q[10];
u3(1.83830674704298,4.84332493041518,0.453340914393457) q[8];
u3(1.07411697461882,2.83241152593597,-1.23021686120933) q[2];
u3(2.37739453346427,1.34812229783838,-2.54514850315220) q[3];
cx q[3],q[2];
u1(2.83984179491095) q[2];
u3(-2.13071287814027,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.604592991883462,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.25866823507277,-0.797933944690438,1.83972171081290) q[2];
u3(1.30224788552647,-3.82676421801309,0.504995242887591) q[3];
u3(2.03748368287024,1.81141305029847,-3.20941041455258) q[8];
u3(2.58554805160619,2.79682759134909,-3.29853553589539) q[0];
cx q[0],q[8];
u1(1.31681382067862) q[8];
u3(-0.0465101492008897,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.38187083394868,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.74426172952540,-1.07453996687371,2.83775641294203) q[8];
u3(2.97526034071320,0.242692374129326,2.70937497547862) q[0];
u3(2.15468154848197,1.87402065947520,-4.17983373927930) q[5];
u3(0.599171648916493,3.57782488759900,-2.33379332116186) q[4];
cx q[4],q[5];
u1(0.177420373068494) q[5];
u3(-0.794624276103011,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.93845354060745,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.72744970376160,0.503944562146563,-2.76136834546840) q[5];
u3(2.86829216875540,1.00364102048692,-0.511705814127326) q[4];
u3(0.258723779239818,3.16266364124563,-0.436370824042212) q[1];
u3(1.17345699432052,0.545482605567333,-3.96302043766894) q[9];
cx q[9],q[1];
u1(2.65737313220438) q[1];
u3(-1.83330611209778,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.24807082335705,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.05313170115947,1.22226813811156,-4.02739681070806) q[1];
u3(0.465767737296074,0.00791315761719336,-0.405305312745077) q[9];
u3(1.92776724274334,2.13542396985067,-3.42088279632437) q[7];
u3(1.61125405318966,3.19265579660625,-2.65433099244185) q[6];
cx q[6],q[7];
u1(1.45426092604654) q[7];
u3(-3.28075418179837,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.51158457595404,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.17862938232406,1.37438741361630,0.331971138111504) q[7];
u3(1.75961608428357,1.25334388962971,-0.558753886341824) q[6];
u3(2.45336515724891,1.95428818878024,-3.61414647260378) q[10];
u3(1.00921950229238,-2.36788899612674,2.94753129439049) q[3];
cx q[3],q[10];
u1(1.39683316231662) q[10];
u3(-0.0290769057867413,0.0,0.0) q[3];
cx q[10],q[3];
u3(0.324543570198810,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.93394047804460,0.239099361719771,-1.83488631107371) q[10];
u3(1.71035434345248,-3.44595940351873,2.46864587360156) q[3];
u3(1.24403092249177,0.181315587668673,1.64421031391756) q[2];
u3(1.13416358700346,-1.26287223055345,-1.81050440289224) q[7];
cx q[7],q[2];
u1(1.77697328164155) q[2];
u3(-3.01195894142277,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.04578639403999,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.640186442744933,0.100928092964120,-3.76735498753060) q[2];
u3(1.93672905954628,-5.30435183038151,0.688760198314432) q[7];
u3(2.57854112445098,0.494440979193710,-0.629995247728896) q[4];
u3(1.74199153229472,0.471936407339722,-4.23123228209059) q[10];
cx q[10],q[4];
u1(2.35856669693162) q[4];
u3(-2.76346390597233,0.0,0.0) q[10];
cx q[4],q[10];
u3(1.51508647508730,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.45244292445265,1.43607394092206,-0.0418976528907996) q[4];
u3(1.23157014021474,-0.769111427189630,4.41407634266426) q[10];
u3(0.699770045042041,0.542396401162767,-0.449177304148876) q[8];
u3(1.12432546454633,-0.558344474207937,-0.628571971589382) q[0];
cx q[0],q[8];
u1(3.31914733360326) q[8];
u3(-0.716887673628366,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.53999484391479,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.967319197322884,-2.62441812418015,0.954359560931611) q[8];
u3(1.18154214515622,1.74216613068054,-1.03444145297001) q[0];
u3(1.97068959650234,-2.27997640058181,3.67294954692896) q[3];
u3(1.18188824179989,2.49018613555131,-2.09181773938018) q[6];
cx q[6],q[3];
u1(1.51151420746607) q[3];
u3(-0.219173389350376,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.54197571488192,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.70631597335355,-3.57049683494880,0.180282977688472) q[3];
u3(0.456316595564368,0.657885241653175,-3.43331511383795) q[6];
u3(2.30080216207890,1.68750576773923,-4.56905409179231) q[1];
u3(0.904083663959314,-2.47503224872419,3.57782685258922) q[5];
cx q[5],q[1];
u1(1.37351058183011) q[1];
u3(-2.87853735755182,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.19018179206738,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.611500771590590,2.09896772620415,-1.01644561482437) q[1];
u3(0.749792784904121,0.404257976616882,-5.13171147527043) q[5];
u3(1.30900039555660,-0.600462679408968,0.851516668681395) q[8];
u3(1.80874922734946,-3.33226626253449,-0.198256167047274) q[9];
cx q[9],q[8];
u1(0.633158060683618) q[8];
u3(-3.24798278385086,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.99861678418016,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.09831221692019,3.66356637792794,-0.188649189150370) q[8];
u3(1.41206051928293,-1.73809901370301,3.28790860925175) q[9];
u3(1.40235918055402,-1.62902219502498,-0.796897586506257) q[4];
u3(1.19025858999570,-3.79610105432006,0.634042065201256) q[0];
cx q[0],q[4];
u1(0.224395096684127) q[4];
u3(-1.41722989341889,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.35606394646731,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.75205057792972,-0.279101966627496,-0.777068679921195) q[4];
u3(1.29262243862186,-3.67849984435162,1.88026906165766) q[0];
u3(0.718020137425612,-0.698478335351825,-0.555710370056975) q[7];
u3(1.44847730301878,-3.74323202053873,-0.579386471663231) q[6];
cx q[6],q[7];
u1(2.50189106498162) q[7];
u3(-2.12686060325666,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.703910240424167,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.85606004248109,3.10487134823341,-1.79023542332650) q[7];
u3(2.09958001972112,-1.20173180294785,-2.21682951242927) q[6];
u3(1.50175922653358,-3.90267399971331,2.09375005971920) q[1];
u3(2.46080119615395,-2.45080495307903,3.51096200908745) q[10];
cx q[10],q[1];
u1(2.82036434893733) q[1];
u3(-2.16018494134409,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.13318586473923,0.0,0.0) q[10];
cx q[10],q[1];
u3(2.04904163666474,0.410845668108651,-1.00204601254453) q[1];
u3(1.53557498460232,-1.11351830111505,3.60447787559606) q[10];
u3(0.935320260098961,2.45391476396634,-0.428221021105905) q[5];
u3(1.44657372365707,0.280273831447996,-3.87954346960946) q[2];
cx q[2],q[5];
u1(1.62157920645108) q[5];
u3(-0.982882472560572,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.45354248267092,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.45607855637811,-2.42194826943123,2.17187791386595) q[5];
u3(2.27119019970966,-3.95359729803356,-0.980432260978317) q[2];
u3(1.80744869768971,0.642414753071093,-0.376090337829891) q[2];
u3(0.415851068574740,0.575664220575214,-3.62548401126411) q[1];
cx q[1],q[2];
u1(-0.122134588138747) q[2];
u3(1.10309152055388,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.59646316892993,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.718333156874865,0.166700499193199,-1.85783438317531) q[2];
u3(1.33433241273336,2.44210860478656,1.33763327773720) q[1];
u3(1.10382710988034,1.27846547463681,-1.05910344720545) q[3];
u3(0.644183061966825,-4.21182894806945,1.06234759233270) q[5];
cx q[5],q[3];
u1(3.02478857464596) q[3];
u3(-1.50874835423996,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.885137072714853,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.97189052113490,0.856100070778622,0.488185777435498) q[3];
u3(1.43295130259321,-0.895006998317993,-2.15210416764566) q[5];
u3(1.16462673908648,1.07937198783386,1.38722999952601) q[0];
u3(0.757572103017700,-1.51273725213024,-1.94738362035955) q[9];
cx q[9],q[0];
u1(2.41193410224043) q[0];
u3(-1.85757601818474,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.381974016589498,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.16935361784705,-1.81426505776859,3.58026362784619) q[0];
u3(2.74264861222891,-4.04973476320119,-0.823546047579378) q[9];
u3(0.180156989579144,1.20183065786554,-2.09245627800373) q[8];
u3(1.15586798715649,-3.46522795527525,1.62187179519620) q[4];
cx q[4],q[8];
u1(3.22561920133003) q[8];
u3(-0.782128586296342,0.0,0.0) q[4];
cx q[8],q[4];
u3(2.09798640282663,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.56235633690178,-0.305555652764772,1.53891442496158) q[8];
u3(1.99544813216990,-3.51639289639264,0.366304860654824) q[4];
u3(1.66522323776311,0.508470326613524,-3.50655371000656) q[10];
u3(1.21588020946921,-1.05015151689727,4.65434282426270) q[7];
cx q[7],q[10];
u1(3.28493312799663) q[10];
u3(-0.846023659824504,0.0,0.0) q[7];
cx q[10],q[7];
u3(1.75315749842179,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.353424547159337,-0.722411134090980,0.400534343866718) q[10];
u3(2.20125170064289,0.958316368212399,-2.17166201943689) q[7];
u3(1.78990986591063,1.89221737194473,-3.24752630027278) q[9];
u3(1.25545701179327,2.16263613059576,-2.12976535164277) q[6];
cx q[6],q[9];
u1(3.25524398036513) q[9];
u3(-1.95813852315571,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.31349461115168,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.800353317608232,1.31153588286037,-1.74177954886926) q[9];
u3(2.39663816464464,-2.42231408391761,2.53816633603697) q[6];
u3(1.09987097254504,-0.460866100984047,1.60088399482198) q[3];
u3(0.991592241522868,-0.946032179065743,-0.992696717670579) q[4];
cx q[4],q[3];
u1(0.897831663095909) q[3];
u3(-1.45715286615611,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.131661767099392,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.66843387306801,-0.795097386513528,-0.370445766813627) q[3];
u3(1.23685425984060,-2.10374341831372,-3.30264072267989) q[4];
u3(1.43204820498363,2.82600433620215,-3.12693944720279) q[7];
u3(1.60484346980726,2.93933320087009,-2.29497686758030) q[5];
cx q[5],q[7];
u1(2.24708521121125) q[7];
u3(-1.48910806145207,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.191267642666318,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.01972620848674,1.84999129277002,-0.443138162075316) q[7];
u3(1.92595497864019,-2.44601288184330,3.11439042987172) q[5];
u3(0.437424646238624,0.862287095104287,-1.80529252308618) q[0];
u3(0.670631936744812,-1.04839788031803,-0.0358866245426652) q[10];
cx q[10],q[0];
u1(0.854536885217997) q[0];
u3(-3.16030355148997,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.47061594838723,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.69195534243856,-0.279192844885403,-3.33777833421241) q[0];
u3(2.72075483375087,0.914529560658775,2.51132825575224) q[10];
u3(1.57970073005792,0.198406238976232,-1.76405052231864) q[2];
u3(1.68875796033707,3.07186324324062,-3.14092743418889) q[1];
cx q[1],q[2];
u1(0.976651671444843) q[2];
u3(-3.12945077028043,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.81978566124700,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.18785958879145,-0.645636113900836,1.88500730256289) q[2];
u3(1.54819573085203,1.66398884748064,-2.22050469961728) q[1];
u3(1.66724606367840,0.428524177396865,0.366449101839498) q[4];
u3(0.808039506412224,-2.33865044271926,-1.44591702094888) q[6];
cx q[6],q[4];
u1(2.65584528911777) q[4];
u3(-2.83938072171471,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.35153893121486,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.42302438426922,0.769135032921610,-0.916222416818524) q[4];
u3(0.580807994855807,2.13019174529026,1.21611126418375) q[6];
u3(0.880261854086383,-0.306793136910303,-1.26408091777123) q[10];
u3(1.94471289495436,0.712202558472632,-4.93330185691026) q[8];
cx q[8],q[10];
u1(-0.100433808367190) q[10];
u3(-1.90085757878232,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.490467717869288,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.67482739851175,3.67810190671238,-1.97778272025538) q[10];
u3(0.655563494485357,-3.29950920976561,1.39321647625620) q[8];
u3(0.644936423375434,1.65273902817338,-3.09635350717892) q[0];
u3(0.849980568990346,-3.46164956721835,2.34672343243871) q[5];
cx q[5],q[0];
u1(1.24860040774514) q[0];
u3(-0.457352318821392,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.01865001905741,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.88838799474746,-0.276559533856837,2.31732614799551) q[0];
u3(2.62592439989661,0.292420267140899,2.26413580186717) q[5];
u3(1.27195430488628,-0.719279284546478,-0.729005300712639) q[3];
u3(1.12508134721752,-3.07881047529682,0.640896570109743) q[1];
cx q[1],q[3];
u1(2.16288311175228) q[3];
u3(-1.66517653463372,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.258003561422828,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.65064324537649,-2.77742137573846,-0.406761254601466) q[3];
u3(1.68472814678497,-1.28725861085297,-2.15483764728206) q[1];
u3(1.43074876466676,-0.168052412868441,1.45876148407312) q[9];
u3(1.55369693685107,-2.82819191098117,-2.61157866321966) q[7];
cx q[7],q[9];
u1(1.34853092624138) q[9];
u3(-3.41602867536551,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.48640044484200,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.945066839233221,-2.35283793831893,-0.297765262200665) q[9];
u3(2.25308375592929,3.89214765527431,0.815387115691613) q[7];
u3(0.619931111921812,0.133160074752702,0.141248101202205) q[3];
u3(1.15366736735405,-1.04277886060410,-1.68942757704616) q[4];
cx q[4],q[3];
u1(3.26848952261653) q[3];
u3(-0.599716254587089,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.51106189392607,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.87102478779152,1.46164198200465,1.41479272145877) q[3];
u3(1.18892308364401,-0.102000300902745,-1.23119538411256) q[4];
u3(0.205821328226916,-2.05418288006920,1.82197248227791) q[8];
u3(0.961778959471441,-4.21808000404750,1.66698403058455) q[6];
cx q[6],q[8];
u1(4.34979586743075) q[8];
u3(-3.94972152353761,0.0,0.0) q[6];
cx q[8],q[6];
u3(-0.507946461938384,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.715655194030994,1.74946039170295,-2.52580750736913) q[8];
u3(1.80616156360667,-0.696676748693017,-2.93137194853980) q[6];
u3(1.77198537989266,-1.24810741354348,0.656084580524136) q[10];
u3(1.54665179003562,-3.06009182448012,0.177777766695159) q[0];
cx q[0],q[10];
u1(2.16746569833323) q[10];
u3(-2.58790429049265,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.34638127126078,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.26009603456239,0.103946532700060,-0.168244283389551) q[10];
u3(2.40512989794945,0.868933074683263,1.93167828207344) q[0];
u3(0.509856469940781,-2.31836022331242,2.51446930676983) q[1];
u3(0.526771296824739,0.955301594216986,-3.29530343381788) q[5];
cx q[5],q[1];
u1(2.41206729381369) q[1];
u3(0.127276276693487,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.16965002968183,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.23108811829387,-1.42979345052419,2.15698091754093) q[1];
u3(1.61679042020125,1.04645721309736,-2.18135659570979) q[5];
u3(0.394959992009974,3.44544636768538,-2.47239678903254) q[2];
u3(1.06183553073908,-2.54831587272888,1.20250962675810) q[9];
cx q[9],q[2];
u1(0.0660415175520201) q[2];
u3(-1.32367992962456,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.63858239290401,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.93478441052311,0.227332165667847,-1.55787984163888) q[2];
u3(1.82495475880652,-4.92942376490302,-0.306458359345445) q[9];
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