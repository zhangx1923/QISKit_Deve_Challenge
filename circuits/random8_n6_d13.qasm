OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(2.13451418152018,1.48652078599925,0.715550118054298) q[1];
u3(0.849851077442944,-5.06564355524911,0.539084329864032) q[2];
cx q[2],q[1];
u1(2.52118480790919) q[1];
u3(-1.77788652308717,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.05314248332667,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.16159776235610,1.56550381110905,-3.88210291551247) q[1];
u3(1.82592117720262,-1.64777338977652,2.72332377172715) q[2];
u3(2.16050601499006,3.20825852515926,-2.07506742064647) q[0];
u3(1.44803278437684,2.17300538961932,-1.45313151763912) q[5];
cx q[5],q[0];
u1(-0.619900014243652) q[0];
u3(1.04920042208386,0.0,0.0) q[5];
cx q[0],q[5];
u3(4.05405122572003,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.629536350570102,1.65339337603130,-1.96211475134146) q[0];
u3(1.57813523187103,0.802892007598369,-2.53563638040736) q[5];
u3(0.922036804162014,1.89494105308894,-0.272254333840783) q[4];
u3(1.96363782335575,0.521005601860037,-3.32894157868864) q[3];
cx q[3],q[4];
u1(1.58100762044517) q[4];
u3(-0.277127981103813,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.40308223006646,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.57129087667871,3.47349227917091,-2.56293997435953) q[4];
u3(2.51081519790185,-0.690222950224683,-0.535689183495259) q[3];
u3(2.01927324507065,2.06432613746727,-4.11912170092376) q[0];
u3(0.0630220232249601,0.531578398109256,0.925955347026083) q[1];
cx q[1],q[0];
u1(1.19722965982088) q[0];
u3(-1.48630327630785,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.426203190094401,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.15116639117465,-2.57461317959992,1.58316786839627) q[0];
u3(1.10381580880651,-0.425806729962697,-2.51143103378338) q[1];
u3(2.25412573729555,3.95329054891634,-2.29936602236314) q[2];
u3(0.905209252787150,1.01119846329052,0.844153051188216) q[5];
cx q[5],q[2];
u1(1.66525552235074) q[2];
u3(-2.13563091626747,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.155497728512618,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.623833313605927,2.92179343737931,-2.00284815883256) q[2];
u3(1.91724574675384,0.455032129983226,-4.25866601750311) q[5];
u3(2.09316690535374,-2.36526161023628,3.76229428883605) q[4];
u3(1.09786033338542,-0.785390619225323,1.86669939333059) q[3];
cx q[3],q[4];
u1(1.17040646387090) q[4];
u3(-0.235662959422123,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.40105394464970,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.66852467066116,1.93027549980847,-3.75088898271269) q[4];
u3(1.55614924903436,-2.67906114819769,2.83686095355652) q[3];
u3(1.55832911270390,1.91849805283295,-3.85385740745861) q[3];
u3(1.83963982724218,3.25270605853990,-2.55050581222231) q[5];
cx q[5],q[3];
u1(1.59619278554761) q[3];
u3(-3.41291438506327,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.31791107571426,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.10380896535049,0.276947459128817,-0.830563433370500) q[3];
u3(2.05684980222576,0.389531610474981,-2.16216102332330) q[5];
u3(2.60946435818119,-0.146902084097636,2.48822567441979) q[1];
u3(2.32462223621814,-1.35270058050480,1.13532550738684) q[4];
cx q[4],q[1];
u1(1.16944778724850) q[1];
u3(-0.266392069108536,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.53401400985054,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.976070661127407,-0.913752027639406,2.60786659326094) q[1];
u3(2.79873878914012,0.742220494239944,-1.61547251510537) q[4];
u3(1.31975676337401,2.89712085197245,-2.12467471838458) q[2];
u3(0.722170724883555,-3.49185605402932,2.74218907746243) q[0];
cx q[0],q[2];
u1(2.84977030887645) q[2];
u3(-1.97999377646878,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.33108056821763,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.66821401081009,-2.26698085338158,1.13327149499006) q[2];
u3(2.83302568845706,5.59454761976881,0.264495853724291) q[0];
u3(2.25806801111383,2.12684717970787,-2.77305883932438) q[2];
u3(1.07081778835645,-2.55807560057822,3.02800143392452) q[3];
cx q[3],q[2];
u1(4.08933299447838) q[2];
u3(-3.50745352844860,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.0978574591871992,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.03244789442086,-1.99768906823564,3.03885059258600) q[2];
u3(1.37069546264710,-1.76483181869512,-0.536614395441689) q[3];
u3(2.22888473958036,-0.491821280541492,-1.40873226970720) q[5];
u3(1.36210340590078,-5.89965630055192,0.380948833855866) q[4];
cx q[4],q[5];
u1(0.793547281353795) q[5];
u3(-1.31869013162588,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.49206675179754,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.95321880761304,-1.50817624262051,2.26569350599650) q[5];
u3(2.12473014909019,3.03530822279865,0.139848782448174) q[4];
u3(1.16965988076861,0.555488843046069,-3.38763687316414) q[0];
u3(0.546866862398584,-2.93766568782282,2.88370119895178) q[1];
cx q[1],q[0];
u1(0.815034409972431) q[0];
u3(-3.37667553511895,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.52293799184946,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.61323559225962,-0.425138298430430,-0.804182350045199) q[0];
u3(2.71831709150889,-1.66966785652800,0.203132806526525) q[1];
u3(2.72678947854700,1.38472938744027,-4.43756485296514) q[5];
u3(1.03405482995603,-1.00511563973217,1.75511518869379) q[3];
cx q[3],q[5];
u1(2.77721745430552) q[5];
u3(-1.88650455647857,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.791937120523730,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.168998277239764,3.44913290720419,-1.33110911432843) q[5];
u3(1.96080914038559,1.34622185464634,-4.90308868012991) q[3];
u3(1.08259291369304,1.43448971038687,-2.99252345348787) q[1];
u3(1.61516079747163,2.06301032997522,-4.00201784035471) q[4];
cx q[4],q[1];
u1(-0.460281629918367) q[1];
u3(-2.33963406043385,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.23112348869465,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.86737935098597,0.269927827779700,0.994851784596304) q[1];
u3(1.90531441126826,1.38324850971258,-4.66714634005869) q[4];
u3(2.64050522709552,-2.35618051901402,2.94651828851872) q[0];
u3(1.20895333232604,0.571808063063362,0.712665936965881) q[2];
cx q[2],q[0];
u1(2.48878513213118) q[0];
u3(-1.81157326468895,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.08480062369354,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.05231403147877,-1.39493286995109,2.63419338972758) q[0];
u3(0.962643613278623,0.367490432546122,-4.63656553723493) q[2];
u3(2.09121787618954,2.37270470387423,-3.46809718405547) q[5];
u3(0.660203654229819,-1.98805837412538,3.71274500944653) q[4];
cx q[4],q[5];
u1(0.185072569548610) q[5];
u3(-1.47073092382359,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.31726449433726,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.96266093004990,-0.151467804556035,-2.55228266610481) q[5];
u3(0.661648940529238,-2.18219884644332,1.69642420220293) q[4];
u3(1.47387841208062,0.234784470610057,1.10206577527211) q[0];
u3(1.06673701086243,-2.41339319131520,-1.83044975663881) q[2];
cx q[2],q[0];
u1(2.94081569147255) q[0];
u3(-2.67751026745418,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.628350502740100,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.25397190778702,-2.85679089637814,0.774162041613668) q[0];
u3(0.331077725238854,-3.49703006426407,-1.50930197300070) q[2];
u3(1.98554505119777,-2.06135138921858,-0.187425676018183) q[1];
u3(1.03925366764967,-3.87125610093839,-0.00802330430564435) q[3];
cx q[3],q[1];
u1(3.40070213501587) q[1];
u3(-4.32563459770520,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.218399047138321,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.24808408438214,0.752380416484250,-0.343150843607534) q[1];
u3(2.29240191262263,-1.15397184437174,-4.34008213223353) q[3];
u3(2.45843718349294,2.29473424909155,-3.27466283329536) q[3];
u3(0.984959672363613,-0.276143979747670,1.21688342337625) q[1];
cx q[1],q[3];
u1(0.184677883551954) q[3];
u3(-0.652308314715185,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.69493452280250,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.43548603886329,-3.35455182879421,-0.228429189361286) q[3];
u3(1.85967139789426,1.21702393540210,-4.07738919761536) q[1];
u3(1.77320556907951,1.99103499661468,0.476982920720640) q[4];
u3(2.01954460205251,0.915876199873784,-2.24967597593245) q[2];
cx q[2],q[4];
u1(1.86473289217679) q[4];
u3(0.363534085149200,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.05086569585120,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.34947175185189,-1.45316300057458,3.09541724800466) q[4];
u3(2.35832299545079,-0.0907914296277600,1.19378714524627) q[2];
u3(1.31572480950542,-0.213699865625787,-1.69978046452140) q[5];
u3(1.13050951688585,0.781983109723637,-4.59351871272190) q[0];
cx q[0],q[5];
u1(1.33559252845746) q[5];
u3(-3.03308277625418,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.41979348891105,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.459595700365699,-0.0680699996588707,1.25288383084932) q[5];
u3(2.27237111099413,1.30985917589161,-1.86062740970987) q[0];
u3(2.25277782479696,1.17635313544810,0.919347124504430) q[5];
u3(0.797382508103724,-1.85673652000712,-2.88734973513410) q[4];
cx q[4],q[5];
u1(2.82710369827510) q[5];
u3(-1.90475031968210,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.26699726522302,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.12768322194691,-0.173477495998094,3.07335307721894) q[5];
u3(0.608371466909194,3.55871774196013,-0.980330512377154) q[4];
u3(0.945557930327888,-1.83274108632535,0.0665084312533996) q[0];
u3(0.793948295259227,-2.93265963236574,-1.30374292641211) q[1];
cx q[1],q[0];
u1(0.673862747156691) q[0];
u3(-1.37670093850794,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.00644219202893,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.07681251899518,-0.658361287702854,4.24045677699757) q[0];
u3(1.35774601541025,-2.66046632047231,2.34502935598988) q[1];
u3(1.23982688229709,-0.911509550650295,2.01179010090353) q[2];
u3(0.688343165859302,-1.08727261058690,-1.18903916390741) q[3];
cx q[3],q[2];
u1(2.61851105475754) q[2];
u3(-0.0553461845654766,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.98314433573539,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.11323257027824,1.02410296315348,-3.70009811565129) q[2];
u3(2.06525338281558,-4.55863198896096,-1.63039085318250) q[3];
u3(1.01678511823753,2.26164759732647,-3.27358036639528) q[0];
u3(1.23345640528701,2.36437005373690,-3.18334459963944) q[4];
cx q[4],q[0];
u1(0.194817670862263) q[0];
u3(-0.531502510839208,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.36122169254770,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.02863716076738,-1.68689224042287,4.36323181906647) q[0];
u3(0.872986444497216,0.760665618129463,1.83290607655396) q[4];
u3(1.91235654475514,2.30424937655627,-0.399829848505185) q[2];
u3(2.52317045207381,0.0923195817078404,-3.84067848239294) q[5];
cx q[5],q[2];
u1(0.823939161642284) q[2];
u3(-3.36454621323185,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.75587696828351,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.49268794133325,1.98450233573779,-1.65031758666297) q[2];
u3(1.15411752155327,-0.598823253740201,0.690100736523358) q[5];
u3(1.61480538880058,0.126861600553574,1.07665700398550) q[1];
u3(1.44532449696851,-2.11097905776886,-1.37897819768701) q[3];
cx q[3],q[1];
u1(0.589234765006735) q[1];
u3(-1.01861603305470,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.50923150878235,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.45913998181728,-3.08890479402197,2.13323082105061) q[1];
u3(2.87263455956021,1.97967569420866,-4.30233740135957) q[3];
u3(0.689398014054473,-1.94140956270462,2.68122849968526) q[1];
u3(0.624301071804908,1.04085058524485,-2.82935089737756) q[4];
cx q[4],q[1];
u1(2.76695850044067) q[1];
u3(-1.47700587217506,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.15846291954093,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.00255291739108,-1.46621950956714,2.62903347983158) q[1];
u3(1.09814096912136,-2.74421332449700,-2.04965010608013) q[4];
u3(1.05009547252375,1.16235370640208,-0.689241002832905) q[3];
u3(0.457711674412647,1.00085626530000,-1.93090694685027) q[0];
cx q[0],q[3];
u1(0.570013976985269) q[3];
u3(-0.139383365852362,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.49428692350081,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.936050769375569,-1.45366588352933,3.14514824106261) q[3];
u3(1.61377308763155,-1.91718776773790,2.92604317907176) q[0];
u3(1.13315356971906,-0.286642504571006,1.55639616071050) q[2];
u3(1.11663758119833,-2.13058170155002,-1.80352895087245) q[5];
cx q[5],q[2];
u1(1.58572390969455) q[2];
u3(0.121124280879654,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.31554349824726,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.76342126853924,-1.53167514025068,1.91553015463983) q[2];
u3(0.291616065382606,-3.71898902355809,-2.24844501372263) q[5];
u3(1.08985761579193,3.46154966179279,-0.618012306322715) q[5];
u3(1.08865926206930,1.03687807844388,-1.38389254511755) q[1];
cx q[1],q[5];
u1(2.08931741106311) q[5];
u3(-1.77252309991064,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.833735411288259,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.801686750467190,1.54376765839691,-0.262349841404955) q[5];
u3(2.20694513586840,-4.41777337786154,-1.18849798996342) q[1];
u3(1.19650299785903,0.644535787850426,-2.22916970379986) q[0];
u3(0.966073928883259,-3.56401724747828,2.14859089117966) q[3];
cx q[3],q[0];
u1(2.05010251490526) q[0];
u3(-1.71683237696659,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.14797728169968,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.38582126971287,0.0294183100120865,1.75851768338177) q[0];
u3(2.44671856860625,-2.14868855513882,-1.53499538569664) q[3];
u3(2.52625986328197,0.0581069928767233,-2.24013839913869) q[2];
u3(2.14390002175174,0.384454756704869,-5.22530120099434) q[4];
cx q[4],q[2];
u1(2.93632872769024) q[2];
u3(-1.77191554728756,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.153573015582498,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.906363506041003,1.72416186013953,-0.148941500974704) q[2];
u3(2.32730424527324,4.45592963479048,-0.625467987395712) q[4];
u3(1.33246964975989,0.521042076812047,-2.83209700176098) q[2];
u3(1.10642441170766,-2.78895051678332,2.69355847056655) q[5];
cx q[5],q[2];
u1(1.69776534062181) q[2];
u3(-2.23454580044009,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.96138748646449,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.34523124146394,3.04074602422382,-3.17780674071361) q[2];
u3(1.70220027708409,0.887605124829151,3.06183146599519) q[5];
u3(1.67620180458861,-2.23848130404443,-0.218674451793396) q[1];
u3(0.194003975065594,1.15609870222433,4.15089491928626) q[3];
cx q[3],q[1];
u1(2.52205674665745) q[1];
u3(-1.76200358013135,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.09262349955736,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.04933783194635,-2.89043175669401,0.843306990662144) q[1];
u3(2.66200077007474,1.30554244005380,1.28247296850023) q[3];
u3(2.18455132736259,2.11868090045916,-2.97124854706891) q[0];
u3(1.81208559074109,0.471025269212005,-1.68296283152181) q[4];
cx q[4],q[0];
u1(0.653296634198108) q[0];
u3(-0.215499050225030,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.03008127444213,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.17116816890308,-1.16937065215362,-1.58119254997981) q[0];
u3(1.28092945495442,0.266754280925279,-4.15543048083554) q[4];
u3(0.749391127553780,-0.412685188819060,1.55043472285199) q[4];
u3(1.17375736655823,-0.558938955456782,-1.73578129922703) q[0];
cx q[0],q[4];
u1(1.52873379866764) q[4];
u3(-2.80402695814763,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.0355747963882758,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.43420898268262,-0.538159812997183,-1.11881363362980) q[4];
u3(0.882032825411264,-2.25087934149532,0.227324166428065) q[0];
u3(2.30466249043053,0.790201323734295,2.03697688297398) q[1];
u3(1.39873815553915,-2.38512439429107,-2.86204424117082) q[5];
cx q[5],q[1];
u1(1.12257405836872) q[1];
u3(-1.31999816792500,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.0699757427240151,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.14591932590562,0.957456705515231,-1.10677391879235) q[1];
u3(2.75275194234201,-3.64039159926788,-0.470085612723498) q[5];
u3(2.01031791394301,-1.39842379246121,-0.283228441078706) q[3];
u3(2.01472025946567,-1.94014997686738,1.49401706123235) q[2];
cx q[2],q[3];
u1(0.673044671708670) q[3];
u3(-1.43446798875822,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.09084851754014,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.281731008655072,2.18286138438299,-2.28367544477546) q[3];
u3(1.90995942267214,3.96121242132208,-0.846594837673313) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];