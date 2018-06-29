OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(0.0808627800754111,3.11281636568190,-2.58781826791941) q[7];
u3(1.45147523159901,2.12155543721359,-1.81293606409080) q[3];
cx q[3],q[7];
u1(2.88018942394915) q[7];
u3(-1.55963461367189,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.776772788681099,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.44357399310089,1.80878959282650,-0.311793663596952) q[7];
u3(0.858739321661816,0.730824964841032,0.505572254821749) q[3];
u3(0.280528592485687,-2.27371471675758,1.41577480407316) q[9];
u3(1.02066380318497,-0.854422357744390,-1.07129204956066) q[5];
cx q[5],q[9];
u1(3.23376857579592) q[9];
u3(-1.48683851722745,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.23978230825643,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.38393660813910,0.702040725161784,1.43196327623742) q[9];
u3(1.66087570629478,-0.373392645892544,0.100504437840019) q[5];
u3(0.546810381418521,2.65351407696602,-0.316977907810649) q[10];
u3(1.46658854894871,0.713675728629283,-3.27612720689147) q[8];
cx q[8],q[10];
u1(1.76224020592179) q[10];
u3(-2.19308725437700,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.304876434511011,0.0,0.0) q[8];
cx q[8],q[10];
u3(0.804329710153432,-1.75661141567993,1.23692926782336) q[10];
u3(0.647221981785226,0.507118539203789,-1.28118841727866) q[8];
u3(0.927088523474289,1.48903075459644,-3.11844622220853) q[2];
u3(1.77391282209356,2.15083968129311,-3.25671217923116) q[0];
cx q[0],q[2];
u1(0.203196065519953) q[2];
u3(-1.53481489234213,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.32516141846981,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.75307043981021,-1.48525473836598,3.07061856532928) q[2];
u3(2.43245086757933,0.310526135402873,0.618343982387624) q[0];
u3(1.85422355430245,-2.15634257337669,1.01372766086101) q[4];
u3(2.46189332673647,-3.71495154853971,0.175919288874502) q[1];
cx q[1],q[4];
u1(1.78957840357845) q[4];
u3(-2.58365593054422,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.0555919714823476,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.47089023583822,2.35035244776865,-2.04998991432971) q[4];
u3(2.64145705831014,2.87142367288996,0.0493349088788264) q[1];
u3(2.04840735504937,0.528811154003805,-0.268720767354414) q[8];
u3(0.411552548813204,-2.93795241842454,-1.75340265709329) q[2];
cx q[2],q[8];
u1(1.29372370029810) q[8];
u3(-3.36609961440475,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.51509087873221,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.457746012470724,0.103972721537061,-3.27424229539525) q[8];
u3(0.432681528320893,1.41099156233933,2.44497914056264) q[2];
u3(1.18391974313771,3.71299298794469,-1.11619419947844) q[5];
u3(0.495705882245206,0.618483808152835,-0.565591666999062) q[3];
cx q[3],q[5];
u1(-1.00416801099801) q[5];
u3(0.320025755574509,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.11428640103030,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.34792048186005,1.86097699037340,-1.04709436720769) q[5];
u3(2.05189402040705,2.02976152392262,-0.0680751002163600) q[3];
u3(2.26765126282088,2.53049925380914,0.328813483196027) q[7];
u3(2.28881613141400,0.248105497544993,-3.52318297724177) q[10];
cx q[10],q[7];
u1(1.37532917413790) q[7];
u3(-3.37423397776243,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.39778808719115,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.51221856750573,-4.01310946393167,1.74350851359065) q[7];
u3(1.45708300374785,-3.28861461177985,1.50741151526104) q[10];
u3(0.661852020755036,-1.16024364252433,-0.229750245119597) q[1];
u3(1.29147275423617,-2.82106886442867,0.565328950140654) q[4];
cx q[4],q[1];
u1(2.00760906590292) q[1];
u3(-2.87902151347774,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.782035388669759,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.771844031064765,1.89486612827032,-1.33771436708493) q[1];
u3(2.69700707484486,5.51404198699320,0.384945971681347) q[4];
u3(2.34709897525771,-2.94259923484926,1.02500615858778) q[9];
u3(2.68083304007680,-2.52654450401264,-1.14484253861461) q[0];
cx q[0],q[9];
u1(3.72206506980208) q[9];
u3(-1.33030675596804,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.09496476546325,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.705394744349992,0.956050100352396,0.514352628811467) q[9];
u3(1.56014181199371,-0.338568440557733,-3.06669967911337) q[0];
u3(2.40787865678968,-1.83951422970828,0.797147037253452) q[5];
u3(2.52568049476134,2.28271290620373,3.15742514139955) q[4];
cx q[4],q[5];
u1(1.94717822775964) q[5];
u3(-2.26554303293850,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.40997968533924,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.73446524014526,-2.50048222019798,0.916710401888085) q[5];
u3(1.69664677392736,-1.13149691446894,-4.89136958702645) q[4];
u3(2.39821301436617,3.05305053130509,-0.145082195899341) q[7];
u3(2.96640085700529,1.97766215084386,-3.92928404946052) q[1];
cx q[1],q[7];
u1(1.76705851795465) q[7];
u3(-0.0902612341051998,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.571885623777652,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.28613965836662,-0.494733215000582,1.95312517303507) q[7];
u3(2.52992407448575,2.65833981429581,-2.46822555170596) q[1];
u3(2.32144601157422,-3.79831794389797,1.30720917647956) q[9];
u3(1.53944718005217,0.241558511299256,3.29399067911436) q[6];
cx q[6],q[9];
u1(2.18664558288293) q[9];
u3(-2.60277180984715,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.113031987688159,0.0,0.0) q[6];
cx q[6],q[9];
u3(2.22188146159343,1.26709005851000,2.30819156600030) q[9];
u3(2.21821013963943,1.41473289472735,2.23428806639650) q[6];
u3(2.45940091206918,2.19557184981184,0.0412317695636488) q[10];
u3(1.83357418302870,0.617617996996647,-2.81504064445601) q[8];
cx q[8],q[10];
u1(4.03010184213781) q[10];
u3(-3.57518261288355,0.0,0.0) q[8];
cx q[10],q[8];
u3(-0.104797368669566,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.06741755062751,-0.102345132451992,1.78215465688069) q[10];
u3(1.08964203585562,-0.148151764816242,-6.01389552868369) q[8];
u3(2.70395483789069,0.465909779667416,0.685345212750171) q[0];
u3(1.18108197864369,-5.40337569855606,0.127126224956291) q[3];
cx q[3],q[0];
u1(2.65471160126175) q[0];
u3(-1.93524892528816,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.839167311026416,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.01884057301917,-0.115113643340461,2.50263419422787) q[0];
u3(1.94126524795093,-1.99816056451863,-4.08170107954729) q[3];
u3(1.85777641006312,1.75987418495526,-2.59051071982235) q[3];
u3(1.45726368734368,-3.17810600555088,2.28143860045492) q[9];
cx q[9],q[3];
u1(-0.0908734609358024) q[3];
u3(-2.39255030688161,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.21876393673396,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.91585488363529,-0.864098717081748,3.41535883232365) q[3];
u3(2.12324463944696,-2.75353465746938,1.53780159459816) q[9];
u3(0.937461337491117,-0.0503860404635854,2.09603600868784) q[6];
u3(0.901173746089232,-1.62305490502346,-2.32545385111110) q[7];
cx q[7],q[6];
u1(1.83204115488021) q[6];
u3(-0.201270091878671,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.59486468519276,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.00910472346226,-2.28901571401607,2.48769619774640) q[6];
u3(2.21112072043220,2.59189415047714,2.90196038170168) q[7];
u3(2.17321810256270,-0.0917114772772135,-0.00674537685297257) q[2];
u3(2.44809730283210,-2.32431566006240,0.869790481659742) q[0];
cx q[0],q[2];
u1(1.73636055256719) q[2];
u3(0.0117328771814733,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.863182157336798,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.56905165872903,2.48327737424128,-1.94242288840512) q[2];
u3(2.70367544665710,-0.356382778049076,2.26109249073619) q[0];
u3(1.86726315829350,-0.680355669307929,1.52843629322193) q[5];
u3(1.18006770340799,-2.32917475517388,-0.442985232633431) q[10];
cx q[10],q[5];
u1(0.109031311718771) q[5];
u3(-0.996858599119429,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.61241366713305,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.01229116369168,0.815897463557528,0.260596595221647) q[5];
u3(1.68110540068887,4.50503463115749,-0.320423698427499) q[10];
u3(2.48842479086313,-0.342422044693434,-1.29542233718148) q[8];
u3(1.22380094142778,1.36364004286659,-4.00821103587117) q[1];
cx q[1],q[8];
u1(-0.196041472052845) q[8];
u3(0.594773319291171,0.0,0.0) q[1];
cx q[8],q[1];
u3(3.37658047975262,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.54416801457239,-1.67771509261234,-1.46682810319552) q[8];
u3(1.36183248485820,2.78128986698039,1.44790878118988) q[1];
u3(1.80153484973069,-0.334146774534835,-0.726112761398219) q[0];
u3(1.67113585814538,0.788435085829188,-5.37458563459852) q[9];
cx q[9],q[0];
u1(1.38422419205420) q[0];
u3(-0.225488847088470,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.39668476779321,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.64415882185094,-1.47368045192891,1.86243085235458) q[0];
u3(2.37567524006765,0.921190975159373,-0.764879151099922) q[9];
u3(1.63345986729374,1.21206249128009,-0.814412071381385) q[1];
u3(1.59998469736524,0.0468440037852571,-3.61918964983633) q[5];
cx q[5],q[1];
u1(1.86805443575483) q[1];
u3(0.275357375259477,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.851863837073569,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.670518834349679,0.147816497049714,-1.13230938812102) q[1];
u3(0.672942322988848,-0.506020855079198,-4.88389416010879) q[5];
u3(1.17923424480602,-0.530221153049641,-1.87659327611958) q[4];
u3(1.06465381766462,-4.18748057787581,0.934874711295209) q[7];
cx q[7],q[4];
u1(1.84436780007785) q[4];
u3(-2.64959230792522,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.997992118710447,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.743631870944117,-0.0407497133111143,0.514648245114074) q[4];
u3(1.37930112599409,2.20301948664776,0.190464671615356) q[7];
u3(1.74160340933942,0.927623848917581,1.05301700785223) q[3];
u3(1.73028984561813,-1.46296222263406,-2.21807035447093) q[10];
cx q[10],q[3];
u1(1.68596859606954) q[3];
u3(0.348394100415828,0.0,0.0) q[10];
cx q[3],q[10];
u3(0.638188896382125,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.990570964065975,-2.09359578477414,0.560663513223888) q[3];
u3(1.45919488780902,1.07113187313074,-1.99482021962245) q[10];
u3(1.62346537385487,1.02528559719624,-2.68534120487077) q[6];
u3(2.83076055012425,4.07018260916187,-0.404781038532655) q[2];
cx q[2],q[6];
u1(0.249201772757574) q[6];
u3(-1.08767997932925,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.38673056213113,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.52822226008456,0.636619482744888,-1.77235013162246) q[6];
u3(0.183083784411147,-1.60727574967724,-1.28883053297810) q[2];
u3(2.59735525147172,-0.393654272885066,3.02870410904167) q[5];
u3(2.97540983023111,-1.21081779019262,1.53447437990488) q[7];
cx q[7],q[5];
u1(0.530694905078766) q[5];
u3(-0.770559552618159,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.23014328369863,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.869264220696681,2.89728504025365,-1.71249518363168) q[5];
u3(1.27227572044026,-2.61364348960970,0.606706908557240) q[7];
u3(1.90366104053468,1.87728343612273,0.0782598487759683) q[2];
u3(1.02789807735560,1.02055225908327,-4.33598208567027) q[1];
cx q[1],q[2];
u1(0.611601341323997) q[2];
u3(-3.41510096769077,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.34773634480258,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.50759413176974,2.38279673831540,-3.65533637675555) q[2];
u3(2.55167143314091,0.721781890028853,5.26826666482822) q[1];
u3(2.29776956849538,2.37126241036274,-2.38481868195696) q[0];
u3(1.27549958610288,1.58004051213523,-1.71059298057663) q[3];
cx q[3],q[0];
u1(0.0931876661418278) q[0];
u3(-1.05766284802925,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.70502378225329,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.766202936931167,1.65124064639902,-0.731117482916485) q[0];
u3(1.72580966573839,-0.0903251389239319,-4.88031780027550) q[3];
u3(1.41275608904603,2.30331601871184,-3.47195543964694) q[10];
u3(1.31285166831806,-2.68932228559205,3.54215081979724) q[6];
cx q[6],q[10];
u1(1.78977589813882) q[10];
u3(0.420884836440506,0.0,0.0) q[6];
cx q[10],q[6];
u3(0.725753347300827,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.20533271126976,2.64134440422481,-1.05777877878611) q[10];
u3(1.33026959326100,1.46608857096327,-1.02310838074753) q[6];
u3(1.53016395629091,-2.28484504157451,-0.822276435365194) q[8];
u3(2.06832195324603,-3.03582049634671,0.351201126602294) q[9];
cx q[9],q[8];
u1(1.26594600608272) q[8];
u3(-0.691864756062827,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.0525164498380195,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.482525151348257,2.66731173795275,-2.81101351256815) q[8];
u3(2.75696908773390,-1.53605045706565,-4.04484916686421) q[9];
u3(1.59839058159933,2.01711820751053,-3.93071248853160) q[1];
u3(2.08397325267490,3.00652077911649,-2.56799831175892) q[7];
cx q[7],q[1];
u1(2.20034056193821) q[1];
u3(-1.43947124058741,0.0,0.0) q[7];
cx q[1],q[7];
u3(3.37118001779524,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.76198253397733,3.29328170774005,-0.789519624811190) q[1];
u3(0.319071587609370,-4.27735964250069,0.857384609409134) q[7];
u3(0.825894696580739,-1.91147592967668,4.03349210946690) q[10];
u3(0.944577027908678,-0.463755285901235,1.19635996484715) q[2];
cx q[2],q[10];
u1(1.89451886893475) q[10];
u3(-2.94592674392659,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.20592195750811,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.505241429139114,2.67549790140883,-0.687024618281119) q[10];
u3(1.69678066146406,-4.12989378549357,1.50375074332878) q[2];
u3(1.69458875805476,0.117640914516753,2.04285593480243) q[6];
u3(2.06548290125483,-2.62879412122907,-1.49078573555255) q[5];
cx q[5],q[6];
u1(1.82171157449713) q[6];
u3(-2.96287243739789,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.695195213972458,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.53338895664293,-2.45902985641196,1.61858258219035) q[6];
u3(2.21191981622173,1.83428163252531,0.0968310844253131) q[5];
u3(2.36893093634268,-0.264625913064371,2.01961630713215) q[4];
u3(2.15753099425145,-1.88059174055168,-1.24542134995354) q[8];
cx q[8],q[4];
u1(1.19315236252526) q[4];
u3(-0.210783151926810,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.386770020153419,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.01900564263853,0.364794487521700,-3.45166997255581) q[4];
u3(0.538607088080821,-3.32163793529216,-1.49215109829193) q[8];
u3(0.461959062052637,1.83756027563636,-1.60902315040872) q[3];
u3(0.827114335342329,-0.341402359819032,-0.642722559807022) q[9];
cx q[9],q[3];
u1(1.60491580629575) q[3];
u3(0.167004542321961,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.01808048454407,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.610095673291818,3.08406087244168,0.895656140990919) q[3];
u3(0.902262394165417,-1.27560663084342,-1.23412404608520) q[9];
u3(1.07729629307313,2.71928224830232,-2.97583542806793) q[6];
u3(2.41813602361123,-2.85328985385221,2.76898406580047) q[2];
cx q[2],q[6];
u1(3.23406535973892) q[6];
u3(-0.775263066106353,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.50734637612731,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.57938444128103,-0.832698974043055,2.07510045323296) q[6];
u3(1.06852888001657,-0.0885709294177716,-1.23159614494927) q[2];
u3(1.90348438803630,-2.59578600771370,3.40604545365197) q[10];
u3(0.797722283949039,0.136615833780895,0.836076300823030) q[3];
cx q[3],q[10];
u1(2.07968215314939) q[10];
u3(0.120186922188687,0.0,0.0) q[3];
cx q[10],q[3];
u3(0.690496460935276,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.57996805672195,-2.82889887786271,1.93329350044956) q[10];
u3(1.48514493176725,3.94667119633604,0.103768698978438) q[3];
u3(1.78030395109102,2.80239980970466,-0.0791006903017095) q[1];
u3(1.76739769623226,0.0806790513448039,-4.58169878819156) q[9];
cx q[9],q[1];
u1(3.40165944343588) q[1];
u3(-1.18945351845303,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.29255345942215,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.60817192277995,1.72899035617136,-3.57609467008656) q[1];
u3(1.22318589414106,-3.20863725284535,-2.93963233799555) q[9];
u3(1.06693610649713,1.92251833262052,-0.902829548123233) q[0];
u3(0.338935864566591,1.61027836677497,-3.22645035257217) q[5];
cx q[5],q[0];
u1(-0.110971182423354) q[0];
u3(-2.04091194288156,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.916220720006269,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.54096119479455,-0.343718181508132,2.20732238806018) q[0];
u3(2.18546981052310,-4.68615245464231,-0.215618722684991) q[5];
u3(1.06375411796546,-2.47882150813174,0.691471818696743) q[4];
u3(1.99685107273722,-3.44630766855388,0.275044023144188) q[7];
cx q[7],q[4];
u1(1.40614199792615) q[4];
u3(-3.68899616360015,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.96048539605526,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.00590493744646,0.754917031059509,-3.09706516384623) q[4];
u3(1.49857963687698,-2.47361875903394,-2.12489214452649) q[7];
u3(1.42907100256364,-0.827644342894179,1.97629005119559) q[3];
u3(1.42950341596680,-1.94011949591609,-1.62256002393901) q[5];
cx q[5],q[3];
u1(0.241438166989906) q[3];
u3(-1.12798692975637,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.53422860134421,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.12374127958845,0.138875688868632,1.93320567212881) q[3];
u3(2.77486759848138,-3.23968070448141,-2.08163561376035) q[5];
u3(2.57758975303992,-2.13689182239495,1.59034013453275) q[1];
u3(2.58703472473844,1.39397649649114,3.44737269162180) q[4];
cx q[4],q[1];
u1(-0.591958592462477) q[1];
u3(1.33180258677889,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.81988362917074,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.69747511051478,-0.944787300145147,-0.431690339285170) q[1];
u3(2.47678893417854,-3.50394870693635,-0.960803852550364) q[4];
u3(1.29764702030687,2.82159916918975,-0.999835104244144) q[10];
u3(1.32286118700247,0.476283532934155,-0.375263241580160) q[0];
cx q[0],q[10];
u1(1.34587725885173) q[10];
u3(-3.10847348004270,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.70544362199259,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.74363340819555,-1.70955254402930,-1.75020958600107) q[10];
u3(0.152442882703058,4.08646692407903,-0.319830895003862) q[0];
u3(0.586993870981775,2.57218268792697,-2.18550727980016) q[9];
u3(1.05279960058909,-2.42218400786481,1.36395199628434) q[7];
cx q[7],q[9];
u1(4.04784641999927) q[9];
u3(-1.55227237852959,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.02810924813748,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.34274211778358,2.98884192279363,-3.07674511549983) q[9];
u3(1.57782234196556,-1.67934438865668,-2.82501719220856) q[7];
u3(2.09630212660680,2.33450607072913,-1.42520071836592) q[6];
u3(1.87504503071439,1.32075841490176,-2.70007019118970) q[8];
cx q[8],q[6];
u1(-0.529924555979421) q[6];
u3(-2.30377725965070,0.0,0.0) q[8];
cx q[6],q[8];
u3(1.27829936798339,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.64104103848355,-3.76867844313274,0.775004496362601) q[6];
u3(1.00611065128926,-5.18645121363323,-0.220113647554479) q[8];
u3(2.92230257468866,-2.25433485592545,0.0673782428357352) q[0];
u3(2.19574284460416,-0.432870654240661,0.964595214585142) q[1];
cx q[1],q[0];
u1(1.74147090075202) q[0];
u3(-3.08919806594350,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.609738724640716,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.14755252035848,-1.15882501520807,-2.61624266085272) q[0];
u3(1.49232018360652,-0.331444762955722,0.799338010698178) q[1];
u3(0.557357854721514,-2.38958389756140,1.32435058582030) q[4];
u3(0.826929813861577,1.84425962259385,-3.34892169340870) q[8];
cx q[8],q[4];
u1(0.397424844732903) q[4];
u3(-0.934706363366202,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.85487181812515,0.0,0.0) q[8];
cx q[8],q[4];
u3(0.996236363910586,-2.86566703509131,3.39740357919506) q[4];
u3(1.32147162322249,0.824195228109909,3.26252069941162) q[8];
u3(1.33431920681778,2.95686528555703,-1.84096231859386) q[9];
u3(0.515509344849835,1.40423438710413,-2.14777831644587) q[10];
cx q[10],q[9];
u1(0.904508294591187) q[9];
u3(-1.42885206044114,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.89016814717905,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.516271436962939,-2.26242513843129,2.62515414857880) q[9];
u3(0.979744634409530,0.360926804852351,4.99866877666368) q[10];
u3(2.46498584670521,-1.77021966724677,-0.378032816470526) q[3];
u3(1.65890338168093,-2.76349220423367,-0.648910888104381) q[2];
cx q[2],q[3];
u1(3.37811121421359) q[3];
u3(-1.43909858762124,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.35607510413181,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.637414777498690,0.342249065961526,2.29981430390489) q[3];
u3(1.80303801488552,-2.16765298316517,2.31126867370211) q[2];
u3(2.60254028040784,-2.94136885857307,2.65568043767027) q[7];
u3(0.833565363534844,2.89563815206804,-1.08100392884133) q[5];
cx q[5],q[7];
u1(0.0569118623022595) q[7];
u3(-2.22286220404523,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.23155939548764,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.56515691114903,-1.46851226275728,4.49886888513685) q[7];
u3(1.97190057863921,-0.657479355791220,4.33515151550827) q[5];
u3(1.83768082102776,3.30598411944966,-1.27919529177620) q[0];
u3(1.81078969904919,3.06735728881691,-0.0638556481927457) q[10];
cx q[10],q[0];
u1(3.26352211263204) q[0];
u3(-1.05611706292642,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.24695341878091,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.65204183643599,0.179626249422292,-2.66526864785372) q[0];
u3(0.793668160782091,-2.17958247596742,-2.12625873516843) q[10];
u3(0.759916618027087,-1.72391902743766,-0.925955970070075) q[6];
u3(1.01519041407081,-2.31069274859262,-0.851074031976184) q[5];
cx q[5],q[6];
u1(2.95301410254991) q[6];
u3(-1.79577509169534,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.828936580913318,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.67765647526158,2.70889873435489,-3.21625340359310) q[6];
u3(1.16915034863212,3.11682440488150,2.88793615064898) q[5];
u3(1.30414972313730,0.780744738762056,-3.48228603280717) q[3];
u3(2.24425069328722,3.76511673490582,-2.37594744777384) q[9];
cx q[9],q[3];
u1(-0.305770723453059) q[3];
u3(-1.27201700615696,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.02110532140114,0.0,0.0) q[9];
cx q[9],q[3];
u3(2.38193924049838,2.52526016126736,-0.690796362733544) q[3];
u3(1.97860718491987,-0.425292226589708,4.39473198887471) q[9];
u3(2.89471616094593,3.36713292004298,-2.46693669570748) q[1];
u3(1.28157181254629,-1.25914335742863,2.64558746575463) q[7];
cx q[7],q[1];
u1(3.59387241778246) q[1];
u3(-4.00210563464378,0.0,0.0) q[7];
cx q[1],q[7];
u3(-0.785286950291044,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.29446097088188,0.410284450355262,-0.900411232999594) q[1];
u3(2.74258580517853,-0.586050549937903,3.48731718381625) q[7];
u3(1.92318657805134,-3.99290077549763,2.20479217368864) q[2];
u3(0.361300584110476,2.32169180193626,-1.53505292779977) q[4];
cx q[4],q[2];
u1(1.32642615060763) q[2];
u3(-0.294760586736295,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.61208505303187,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.848594238854037,-1.10551942955020,1.09163279954065) q[2];
u3(1.33368729292120,0.316857678779041,2.41799163754353) q[4];
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
