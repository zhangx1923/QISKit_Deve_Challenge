OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.08877350067311,0.280263494838391,1.88305884114119) q[3];
u3(1.32197886649094,-0.299801523979480,-2.86071394771219) q[5];
cx q[5],q[3];
u1(0.185518269946416) q[3];
u3(-0.409837947605342,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.46958421569578,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.61055323237899,1.75496016719482,-1.95688781627892) q[3];
u3(2.07789502395277,-3.25798292235470,-1.21273185590848) q[5];
u3(1.62800410406123,1.00649475657896,-3.42786322676316) q[8];
u3(1.11908736228411,3.04203233713673,-2.85189807142753) q[1];
cx q[1],q[8];
u1(0.587912383832725) q[8];
u3(-1.13655679467545,0.0,0.0) q[1];
cx q[8],q[1];
u3(-0.0688306172638957,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.40403033606556,-2.11734907859272,2.38908013789750) q[8];
u3(1.94520660913061,-1.05239123326117,4.36511829751450) q[1];
u3(2.33834920906282,-1.38329738713924,0.645802155317440) q[7];
u3(2.20570710183263,-2.25607131127213,-0.000469160850110617) q[9];
cx q[9],q[7];
u1(2.40971908293211) q[7];
u3(-2.65217945305123,0.0,0.0) q[9];
cx q[7],q[9];
u3(-1.26197546471356,0.0,0.0) q[9];
cx q[9],q[7];
u3(3.01413113816657,2.24474383763405,-1.95227702759604) q[7];
u3(2.22841331499044,-3.63900472481695,-2.43334638003755) q[9];
u3(1.50317272983912,-0.112752721464367,2.40887112602517) q[6];
u3(1.31000830286311,-2.24151381457790,-0.830113443479460) q[2];
cx q[2],q[6];
u1(3.67290138305188) q[6];
u3(-4.28055836801386,0.0,0.0) q[2];
cx q[6],q[2];
u3(-0.920866164850527,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.637482125265413,-2.00442221735700,0.675235007823166) q[6];
u3(0.561921442064786,-0.770743628289236,-1.78052793763735) q[2];
u3(1.56653910496670,1.82163751301032,-3.71689644252963) q[4];
u3(2.37189924090957,4.05390834979486,-2.18946612699129) q[0];
cx q[0],q[4];
u1(0.0627128240490298) q[4];
u3(-0.948880875711739,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.51204649847232,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.09984016265461,-1.24200884712497,2.41495732981585) q[4];
u3(2.50311150435131,-0.475852683813848,-4.88900765469063) q[0];
u3(1.42467776602272,-3.36221750246929,2.83799967603439) q[8];
u3(1.67033871332158,-3.26010919029414,3.01250446413427) q[2];
cx q[2],q[8];
u1(1.19735978817924) q[8];
u3(-3.14000516746554,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.67210689598511,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.445773327834713,-1.60319622957467,3.09555518999724) q[8];
u3(0.733178638090325,2.94539964661828,-2.51707874632435) q[2];
u3(1.09433724869975,-2.13176771448516,2.59501631580180) q[6];
u3(0.576933451466699,0.420195083191363,-2.77896237090266) q[1];
cx q[1],q[6];
u1(0.351880623795163) q[6];
u3(-1.52365861890797,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.93951687546352,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.16711508871201,1.66495223819399,-0.375133931205472) q[6];
u3(2.45800064742094,0.0441105866168960,-0.175751281992431) q[1];
u3(1.79759031675706,-1.59395486654427,-1.20177248962279) q[5];
u3(0.126021807996095,-5.19679923196444,0.827221103665894) q[9];
cx q[9],q[5];
u1(1.57603205252887) q[5];
u3(-2.87465210992476,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.523838650778136,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.24467675537793,0.0693375825597036,2.38601214864827) q[5];
u3(2.63429912478571,-3.37059517123037,-0.309291429989551) q[9];
u3(1.95354407139848,-1.55723554241224,0.101302139689782) q[4];
u3(1.72042407959484,-1.70707193089248,0.0963588657208279) q[3];
cx q[3],q[4];
u1(1.14979117680519) q[4];
u3(-1.00839085741069,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.72847455119304,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.973580209344573,0.923879977973692,-2.51021299139563) q[4];
u3(2.37477687931108,-0.231251626800923,3.39747805806909) q[3];
u3(0.815625833324908,-0.813568403760995,3.95449807988255) q[0];
u3(1.51913025970316,1.73194408436425,1.91414033139765) q[7];
cx q[7],q[0];
u1(2.10075968998584) q[0];
u3(-2.68105785674869,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.22651880156484,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.05149676921153,-1.66328614669691,-1.22399646382672) q[0];
u3(1.69063063539123,3.36714938553188,-1.90685650558108) q[7];
u3(1.37319621750376,0.912676513337389,1.55566177278035) q[7];
u3(1.80504544286787,-1.12236260179146,-0.967459418283127) q[8];
cx q[8],q[7];
u1(0.841774663308253) q[7];
u3(-1.34958376688181,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.91834164290643,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.18422621877655,-0.711347321467176,0.465816963663541) q[7];
u3(1.89358992768735,-1.41336580670829,-1.47132126989437) q[8];
u3(1.43336078449027,-1.23546514387968,1.27736990834254) q[6];
u3(0.796794827366140,-1.48566858625259,0.0750239012836008) q[5];
cx q[5],q[6];
u1(2.87243541170607) q[6];
u3(-1.61069308953214,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.62170961436573,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.01405112205692,1.92278351819306,-1.03469661222934) q[6];
u3(0.286421908118995,-3.82673981062052,-0.344496027975796) q[5];
u3(2.34602386287198,-0.395245198554173,2.22307762518857) q[3];
u3(2.00794054394712,-2.35842844052436,-1.54008506636495) q[1];
cx q[1],q[3];
u1(0.203455716918538) q[3];
u3(-1.36004902760815,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.64626881471148,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.99468875413139,-0.974196554078280,0.102519937145840) q[3];
u3(2.56538099920679,1.96233814041558,-3.27881826620143) q[1];
u3(1.74650600509191,2.33900866925187,-2.47477821485108) q[4];
u3(1.24492031514243,-2.72341446043502,2.59479545113110) q[9];
cx q[9],q[4];
u1(3.51389800886420) q[4];
u3(-1.35470026120439,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.36124866467942,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.04436793728899,3.91187682777794,-1.55683011477225) q[4];
u3(1.64727135331238,-1.02330032382583,-0.117447476310831) q[9];
u3(1.13368949276815,0.693900851270184,-0.477933707362759) q[0];
u3(0.365285283343716,0.328673781430291,-3.38863726790068) q[2];
cx q[2],q[0];
u1(4.08922856155107) q[0];
u3(-3.35333878914501,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.788645126777027,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.93921865214367,-4.23388543923393,0.468816522832193) q[0];
u3(1.86743071025126,2.62382160366626,-0.832694360212109) q[2];
u3(1.03549671601351,0.884594490190744,-3.13624551289850) q[7];
u3(1.47810917280109,-2.30971517065083,3.27511450029714) q[3];
cx q[3],q[7];
u1(-0.0504423818401534) q[7];
u3(-1.80635879985571,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.22561610968473,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.290673449622232,-2.14846723445564,3.75092638179013) q[7];
u3(2.00706823504042,0.413251682553294,-5.84530291499717) q[3];
u3(0.483443767992611,-1.44126687207691,1.29427050983838) q[5];
u3(0.736385293107325,-2.82638336673066,2.26090148554086) q[1];
cx q[1],q[5];
u1(1.07143643795228) q[5];
u3(-1.49596170197822,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.82985172947471,0.0,0.0) q[1];
cx q[1],q[5];
u3(3.03395366770856,1.48151778933303,0.900637967083797) q[5];
u3(0.847590100905832,-1.00593485033423,3.36717432342598) q[1];
u3(1.72603426373332,1.05283400214179,-4.05549611215484) q[8];
u3(1.44293163164421,-2.03765748459326,3.33960660742940) q[2];
cx q[2],q[8];
u1(0.446240919508589) q[8];
u3(-1.71098342258169,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.19523795332416,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.81735577264742,0.324751763758005,0.251013476191012) q[8];
u3(0.902267154386594,-3.22598354480669,3.05136341259412) q[2];
u3(0.477262925437796,1.46862061882579,-0.543039205958559) q[9];
u3(0.868111908593056,0.667414097304501,-2.15668229094628) q[6];
cx q[6],q[9];
u1(2.68730616337277) q[9];
u3(-1.74725373937584,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.725513425556333,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.50628717761825,-2.47729484805714,0.627357693120015) q[9];
u3(0.681479326749832,3.65910918444517,0.857612180752324) q[6];
u3(1.62930117864400,3.73146379834431,-1.75313128357819) q[4];
u3(1.04210279082315,2.48523403996770,-2.36896403386800) q[0];
cx q[0],q[4];
u1(1.28960494551227) q[4];
u3(-3.49710011614692,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.38895998696804,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.633181857956004,-2.91368970824427,3.24548857704324) q[4];
u3(1.93053713303403,0.270023556355452,0.386240983176921) q[0];
u3(2.11618364969820,1.51807720792119,0.825507595265962) q[3];
u3(0.762585370293782,-1.59515849354979,-1.83933573908343) q[6];
cx q[6],q[3];
u1(3.07073061320948) q[3];
u3(-1.39289098928520,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.89226608788633,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.53894053085907,2.44281622100791,-1.86639311495667) q[3];
u3(2.71677923561932,-0.815845101595276,0.982652035608190) q[6];
u3(0.232229032733441,-1.90430195117416,1.74733298350330) q[5];
u3(0.784981089510448,-3.11007748151153,2.72090825177693) q[4];
cx q[4],q[5];
u1(1.57888742235237) q[5];
u3(-2.97463912433184,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.608531927170144,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.57831573269234,2.08916421011145,-0.562789218558983) q[5];
u3(1.05353680070991,2.21557981269252,2.75293060197935) q[4];
u3(2.10530377277986,1.86934505529618,-2.14599158064236) q[0];
u3(2.15942288985947,2.33448570201592,-3.30414309798964) q[7];
cx q[7],q[0];
u1(0.362910166134785) q[0];
u3(-1.47221709244301,0.0,0.0) q[7];
cx q[0],q[7];
u3(-0.135306676978990,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.786516738905788,2.93231060383777,0.758821801351204) q[0];
u3(0.567391936344907,2.80417276584809,-1.81412428766177) q[7];
u3(2.56243215694596,-1.13177574845682,1.61079419473023) q[1];
u3(2.25011378421428,-2.83394533327125,-1.86047013918062) q[9];
cx q[9],q[1];
u1(1.78448938986194) q[1];
u3(-3.28340645748369,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.903007051016475,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.973707306878940,-3.22310593729540,0.790121040211479) q[1];
u3(2.85662570832552,-2.42201731887318,-1.06814469657427) q[9];
u3(0.916886020771390,2.71624922598511,-1.43634733055852) q[8];
u3(1.93329178441236,1.00038064985140,-1.04994388904308) q[2];
cx q[2],q[8];
u1(-0.203441600515716) q[8];
u3(-1.73242442501548,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.963300828077874,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.58492926797191,-0.670378997472413,1.90534186111007) q[8];
u3(0.412181021010253,0.127139625296754,-1.76265870922863) q[2];
u3(1.75834201514620,1.91602249924333,-3.09906044831758) q[8];
u3(1.61537838728783,1.95208366809737,-3.06540099925137) q[5];
cx q[5],q[8];
u1(2.39078639609908) q[8];
u3(-3.14352584677613,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.58031978415840,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.88229263917694,-2.45162322018448,1.53311621218158) q[8];
u3(0.971436504914151,-2.56214446091182,-3.20306983335344) q[5];
u3(2.11624067434936,-0.258981426854764,2.77200286771134) q[9];
u3(2.18066177798608,-2.47138129192071,-1.83839585941963) q[2];
cx q[2],q[9];
u1(1.00670804309753) q[9];
u3(-1.42552536641563,0.0,0.0) q[2];
cx q[9],q[2];
u3(3.00339679764235,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.27355777711596,-0.217687728610824,-3.32599414448968) q[9];
u3(0.719879602703348,3.93770863958078,-1.35643396039405) q[2];
u3(2.45737552738867,-1.66718857303702,3.64907156394541) q[6];
u3(0.390683376359267,1.98638599283071,-1.00450000323595) q[4];
cx q[4],q[6];
u1(2.75364298441229) q[6];
u3(-2.31541259418027,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.804243489760748,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.59645706271090,-0.0415942710356907,-1.90047506359268) q[6];
u3(2.56897221668045,-0.807768458532629,0.783563664966425) q[4];
u3(1.42292294844457,3.18960913928306,-2.48009824466783) q[3];
u3(0.886217505398210,1.69271362243112,-1.81383586410075) q[1];
cx q[1],q[3];
u1(1.98329509216493) q[3];
u3(-2.32266181813104,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.11783463519143,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.16721958684105,-1.69666452824531,-1.08024111292202) q[3];
u3(1.14715517696234,0.643225951287225,1.74256487973144) q[1];
u3(1.47980127541754,1.35660111723190,-3.80189199945556) q[7];
u3(0.727689037698729,2.99204532967044,-2.92445933491986) q[0];
cx q[0],q[7];
u1(-0.182620312259937) q[7];
u3(-2.03528217693965,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.49535883064101,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.25092208045066,-3.56183975923677,1.13503840530697) q[7];
u3(2.32266864834512,-0.304439494626024,4.39962220430048) q[0];
u3(1.88779396058727,0.375826848253417,-1.60352270345136) q[1];
u3(1.23956315477605,0.541229477218799,-3.23114172259139) q[3];
cx q[3],q[1];
u1(-0.433084590937420) q[1];
u3(-1.93833314324449,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.27096490714210,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.99601010322652,-0.459344877496971,-0.616594087516219) q[1];
u3(1.03404660996862,0.429119212619230,-3.53584508951276) q[3];
u3(1.99506427270236,-0.465153720112960,0.813334447274085) q[9];
u3(1.24509304094662,-2.17310042684441,-1.67151545087630) q[2];
cx q[2],q[9];
u1(0.279957016442832) q[9];
u3(-1.60525753746682,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.18989478315788,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.924980102316567,1.40915827379901,-1.82717461321994) q[9];
u3(1.01748236110454,-0.0466294057102943,-0.510104483088728) q[2];
u3(2.10333647348097,0.100822841551445,1.69293542255785) q[8];
u3(2.29919476244909,-1.40297904787294,-2.91865006966563) q[7];
cx q[7],q[8];
u1(1.48573875252330) q[8];
u3(-1.03391964665989,0.0,0.0) q[7];
cx q[8],q[7];
u3(-0.503976728129091,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.29475451312181,-0.384454519404861,2.81288948913918) q[8];
u3(0.699538971825851,3.91295263468304,-0.0326744896168365) q[7];
u3(1.19671738928561,-2.08121670719158,1.60936531660784) q[4];
u3(0.862835984588808,1.79802500279972,-2.87610524320923) q[0];
cx q[0],q[4];
u1(0.697725357701436) q[4];
u3(-3.32747412422847,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.43844728584181,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.997371479252262,3.11264616893631,-1.89755958419124) q[4];
u3(1.19667607450063,-3.11012089474067,-3.10861704658198) q[0];
u3(1.67895924391400,3.26850882314972,-0.723715513330105) q[6];
u3(2.05895392540365,2.99359736100888,-0.961295629969523) q[5];
cx q[5],q[6];
u1(1.72198470527751) q[6];
u3(-3.06639078279227,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.567660836506316,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.16234650274497,0.567562897376748,-2.73935254041152) q[6];
u3(2.52614620429195,0.147181216669159,4.05962795795570) q[5];
u3(2.29153791493789,0.256145663388040,1.69629861197437) q[0];
u3(2.07180729448585,-1.34747050549221,-0.308645396071328) q[9];
cx q[9],q[0];
u1(1.26782152447848) q[0];
u3(-2.81647478247018,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.524686379686491,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.79172583113272,-2.35440955470101,0.439675871779463) q[0];
u3(0.940900759039130,3.62449837750191,-0.114188665064846) q[9];
u3(1.67903478597835,-0.348997584295311,0.830365673691853) q[8];
u3(2.12489329082118,-1.30245488757343,-1.69424099724427) q[1];
cx q[1],q[8];
u1(1.55476999711583) q[8];
u3(-3.74338447873806,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.31974327948496,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.993766686707601,1.18264059918297,-1.34608656909574) q[8];
u3(1.17676142254406,-2.51319356977570,2.22198034897457) q[1];
u3(2.17850397019863,-3.17600393723948,2.75954132683326) q[7];
u3(0.521069149934235,3.12109723576139,-1.34972451253760) q[6];
cx q[6],q[7];
u1(0.507908062566115) q[7];
u3(-0.172021771249016,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.87980547485287,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.739105511341077,2.06059424747320,-1.19220825181680) q[7];
u3(2.23928155279174,-2.26362258037180,1.11266199045244) q[6];
u3(1.29406953805707,-1.35538213648333,1.58375440084630) q[3];
u3(0.492614629919553,1.57071781725788,-2.93667245173928) q[4];
cx q[4],q[3];
u1(0.850408544802932) q[3];
u3(-1.53246792446506,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.200575928216035,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.301862194173020,-0.623387476704480,2.53688734755587) q[3];
u3(1.97774656960364,2.26832934076825,-1.46468626382888) q[4];
u3(1.36865376375114,2.02293408979693,-0.682881701316603) q[2];
u3(2.38999184743994,-0.326199671494237,-4.06808155998920) q[5];
cx q[5],q[2];
u1(1.92152447959957) q[2];
u3(-2.24367826005419,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.35277359446326,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.23692889409973,4.05735516013561,-1.56751448644843) q[2];
u3(0.747640463297582,-4.94684465702920,-0.352595045257656) q[5];
u3(2.76868187896759,2.26387758492289,-3.50905946092334) q[3];
u3(1.11397687286047,-1.79277801155114,2.93995199601775) q[8];
cx q[8],q[3];
u1(1.10102646849594) q[3];
u3(-0.974369496572766,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.93536280381293,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.41355740468013,-1.23729116677654,1.20045267031454) q[3];
u3(0.472557915367224,-4.25152076974881,-1.64751630822192) q[8];
u3(0.657327317202389,2.30452704290227,-2.70797398750623) q[9];
u3(0.628226250673623,0.844226307300361,-2.80556830979965) q[7];
cx q[7],q[9];
u1(1.52152779684239) q[9];
u3(-0.604266811917767,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.90331458844762,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.468268409539307,-2.09657738322368,2.41037498290450) q[9];
u3(0.621848491486776,-1.63504772562028,3.71876832215040) q[7];
u3(0.963189921060959,2.23149026580534,-3.03040918199061) q[4];
u3(1.50847812015840,-2.64166035930615,3.45170281543558) q[5];
cx q[5],q[4];
u1(0.0971560872686017) q[4];
u3(-1.09598483942166,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.37878046832106,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.47574892957445,0.852282414561164,-4.39302390532638) q[4];
u3(2.10172604888887,1.16455107364044,3.43557819106904) q[5];
u3(0.638029722416327,-1.90884619937192,-0.634722468696430) q[0];
u3(1.20430848847427,-2.64452275696880,-0.644816839753476) q[2];
cx q[2],q[0];
u1(1.37704805937437) q[0];
u3(-3.70365963995143,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.07389318953909,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.22287243546136,0.585526564867752,0.583171653742524) q[0];
u3(1.63574762414703,4.23812980970294,1.78518032874372) q[2];
u3(0.970976742233213,-2.29714406512936,2.31956891724505) q[6];
u3(0.438368814578617,1.74233557298833,-3.69051871146238) q[1];
cx q[1],q[6];
u1(2.15278718151420) q[6];
u3(-3.16981060643436,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.64271068784519,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.08833843777262,2.38621771591678,0.180086198580577) q[6];
u3(2.74873172556687,1.21221857796981,1.41998626497924) q[1];
u3(1.15976144520768,-0.787517352247162,-1.22364603076439) q[5];
u3(2.28772455101090,-5.34895383766272,0.716151224178726) q[1];
cx q[1],q[5];
u1(-0.522475538164374) q[5];
u3(1.19091416964682,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.98805344996011,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.77463964303899,-0.401850435333072,2.37385432610472) q[5];
u3(1.77142673074859,-4.40202056895900,-1.73828766416424) q[1];
u3(2.89939550705054,2.24868745403072,0.794354984476707) q[0];
u3(1.78651736993574,0.370973433556220,-3.88620026675779) q[2];
cx q[2],q[0];
u1(2.30049217893637) q[0];
u3(0.345340839778733,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.56418676389068,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.70568973088910,-0.801846055764553,-1.24978573273373) q[0];
u3(0.776868333949803,4.54331761811412,0.368468657665752) q[2];
u3(1.69315773468280,-0.800227634534515,0.185438777652210) q[8];
u3(1.55130829891347,-2.27566644060954,-1.57494800471934) q[9];
cx q[9],q[8];
u1(1.86188646068707) q[8];
u3(-2.17732558330782,0.0,0.0) q[9];
cx q[8],q[9];
u3(-0.282467883287830,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.880641341036877,-0.565731098718259,-2.52386172953455) q[8];
u3(2.23049025816409,0.226948114707912,-4.95224948798708) q[9];
u3(1.87744223339826,2.62958925474559,-2.87138406834475) q[4];
u3(0.925955241855639,3.19660438217771,-2.51449466906685) q[3];
cx q[3],q[4];
u1(2.32578126922768) q[4];
u3(-1.64345462867908,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.10098843088289,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.605219953968900,-2.35945706104125,2.50417942812517) q[4];
u3(1.39284594139168,-0.603536104084840,0.188790992834876) q[3];
u3(1.11077533624980,1.79926655121216,-2.61065439105631) q[6];
u3(2.19581814834708,-2.47542006406049,3.38598534562637) q[7];
cx q[7],q[6];
u1(2.95722891938645) q[6];
u3(-2.24844834876891,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.23916346640971,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.42451414066171,3.69816073804779,-0.587351533489821) q[6];
u3(1.54811865035071,5.17355706070302,0.634457326558254) q[7];
u3(1.68212901995621,0.406920845135864,-3.53373887092396) q[5];
u3(1.97899679796376,4.36804045786878,-1.74585279102719) q[0];
cx q[0],q[5];
u1(1.82204506090130) q[5];
u3(-3.68762215508224,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.72951983611680,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.34981113157941,-3.09080738617468,3.00117121207590) q[5];
u3(2.35195407777207,2.97282547895615,1.64075136072537) q[0];
u3(1.27854086019271,-1.32170258268685,0.911237855012970) q[8];
u3(1.93939556018298,-2.72031907136552,0.311198667531379) q[4];
cx q[4],q[8];
u1(1.67373175060853) q[8];
u3(-2.81806159289150,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.843126112330906,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.16659951827044,0.764441385159968,-1.78806507250832) q[8];
u3(1.14452570763245,1.67606664895951,3.22588779091306) q[4];
u3(2.01116422714303,1.07858876635371,0.619359071932437) q[9];
u3(2.38158359951208,-0.444009223712871,-3.80139358302441) q[1];
cx q[1],q[9];
u1(2.82989986908037) q[9];
u3(-1.88687285889373,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.57841208257451,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.856812301232510,2.13646441200268,-0.422224393129063) q[9];
u3(0.629870015860240,-2.63089526489136,-0.523967177838314) q[1];
u3(1.01283498444236,4.15766879199669,-1.65388007541710) q[6];
u3(1.22643739075815,1.67479157248039,-1.61344421148525) q[2];
cx q[2],q[6];
u1(2.88451315709080) q[6];
u3(-1.77089237554436,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.01647573526862,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.49704455376324,1.68458987005151,-0.192597145242360) q[6];
u3(0.603258969830154,-0.557592315204182,-4.22603497838962) q[2];
u3(1.07320277256778,-0.159632735359300,1.71461103325192) q[7];
u3(1.05516903393796,-1.00468967452622,-1.45860317818671) q[3];
cx q[3],q[7];
u1(1.52839286328452) q[7];
u3(0.959463710350775,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.16311655170069,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.41717407689270,-3.03755895314980,0.575416889109801) q[7];
u3(0.905443290728094,-2.52168015074625,1.85248144041485) q[3];
u3(1.18881020224352,-2.48048598935435,1.54818080060110) q[1];
u3(0.284635486919491,0.309619271083837,-2.06496940945271) q[3];
cx q[3],q[1];
u1(1.69749092353556) q[1];
u3(-0.172045684503377,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.20560494700603,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.21718910094663,-0.206684373785416,0.845776341296002) q[1];
u3(1.00784832524205,-4.33656893968668,-0.645127665420304) q[3];
u3(0.528240717324774,0.764853344348220,-2.93312326684991) q[0];
u3(1.23141519745814,2.75227829739418,-3.34508334416385) q[2];
cx q[2],q[0];
u1(2.12921751298990) q[0];
u3(0.299253041353792,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.33487153787883,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.47498487695563,0.277273878023427,3.19852593020275) q[0];
u3(1.46230577738950,1.12999905473771,4.24897622259377) q[2];
u3(1.46677344695973,-0.572804422745973,1.25423605643008) q[5];
u3(1.40485124958867,-1.77925600640056,-1.69720532735468) q[9];
cx q[9],q[5];
u1(1.69723828976532) q[5];
u3(-2.41882082609254,0.0,0.0) q[9];
cx q[5],q[9];
u3(3.10933543973768,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.20569166498827,-1.28938539018754,1.99165621906373) q[5];
u3(0.745319657629781,0.264188647196580,5.06398492922349) q[9];
u3(2.33482915137685,2.76692236255972,-1.83267200405049) q[7];
u3(1.20448162010011,1.15520884960887,-0.763082038412068) q[8];
cx q[8],q[7];
u1(1.49187557497624) q[7];
u3(0.511403933499573,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.842094709888288,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.47027982250251,1.88701183679888,-3.55087686436755) q[7];
u3(1.14616238600993,-3.10117048882713,-0.560912133576554) q[8];
u3(1.41061664250478,-0.379040446734311,0.781737183990502) q[4];
u3(1.39752546844472,-1.42153886078392,-1.14069972360616) q[6];
cx q[6],q[4];
u1(1.70894489698474) q[4];
u3(-0.601668512521426,0.0,0.0) q[6];
cx q[4],q[6];
u3(3.08592256331571,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.686169481646492,-0.827587634199830,-0.660187166538464) q[4];
u3(1.04228555888042,3.52851851111467,1.96912417717170) q[6];
u3(1.52537433394903,1.93552638123623,-2.58754875198291) q[8];
u3(0.587225094590543,-2.42663274720989,2.61048971498353) q[9];
cx q[9],q[8];
u1(1.64846432763381) q[8];
u3(-2.43910390454214,0.0,0.0) q[9];
cx q[8],q[9];
u3(3.24717295247982,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.65772385095003,-2.78842303682257,1.97150752066430) q[8];
u3(1.48307604600204,-1.86313111929960,2.25379113997367) q[9];
u3(1.33627720475302,0.421049482401388,-1.48106286525691) q[5];
u3(1.76103465343888,-4.35864845183998,1.82411588082222) q[1];
cx q[1],q[5];
u1(1.36622134160281) q[5];
u3(-3.15673036782042,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.62606441965024,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.03523310993394,-0.825240278682528,-1.85008077320898) q[5];
u3(1.39991993415298,1.16751747978100,-1.61787256740630) q[1];
u3(2.15620418101160,-2.19655335823354,-0.376607952894755) q[3];
u3(2.31657537880949,-1.76535447100367,0.289232745047519) q[0];
cx q[0],q[3];
u1(3.26228055166797) q[3];
u3(-1.11414645391493,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.37559296193762,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.453823389682061,-1.49828849807359,-0.927661292735388) q[3];
u3(1.57950217108722,-1.50179017317573,-2.69330321647683) q[0];
u3(2.68092557813377,-1.21601228053802,-1.01417241239926) q[2];
u3(0.761815498935150,-4.68465243431742,0.362145965563343) q[4];
cx q[4],q[2];
u1(2.19005702772828) q[2];
u3(-1.73815528103667,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.48242955032008,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.74979531563376,1.31903249557081,-3.26829107849604) q[2];
u3(1.09849623412603,2.27858544445343,-3.28157281041443) q[4];
u3(1.66060362470718,1.44674034577380,1.53234691763436) q[6];
u3(0.0475085682562467,0.158398110645043,-3.50423862778748) q[7];
cx q[7],q[6];
u1(0.834807160030236) q[6];
u3(-0.299166431324598,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.72661074954907,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.859153426629459,-0.402784871400037,0.786519823110268) q[6];
u3(1.69604629751614,-0.00593305110114439,5.22005993872376) q[7];
u3(2.06218477646449,-2.78413335619866,1.28868272720834) q[2];
u3(1.49785932416780,-1.87122641203476,0.501441706627702) q[3];
cx q[3],q[2];
u1(1.25200134784822) q[2];
u3(-3.77451396441807,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.96234456882821,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.78958026983401,2.84652293537510,0.973761248663695) q[2];
u3(2.15740380008104,2.45214168960385,3.76842506799772) q[3];
u3(0.217778074661254,1.59847374738939,-1.77242483267995) q[4];
u3(0.995652013387837,-0.290956555252983,-0.0956945870113524) q[7];
cx q[7],q[4];
u1(-0.0676289921059499) q[4];
u3(-1.93639346029432,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.55842710266171,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.83964330663871,-3.70779833891213,1.82911473339783) q[4];
u3(1.90667561364010,4.39433114734006,1.24597597131370) q[7];
u3(0.907717416483354,0.0341021501294748,1.84666864081019) q[5];
u3(1.30895790008857,-1.04773828576255,-2.75749438434421) q[9];
cx q[9],q[5];
u1(3.22688426081848) q[5];
u3(-3.75375543651971,0.0,0.0) q[9];
cx q[5],q[9];
u3(-0.894279790992123,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.456158228945544,-1.06721096874090,0.236386038385138) q[5];
u3(2.08831621811774,3.25681591316704,-0.183060899346593) q[9];
u3(1.62089957611508,-0.435420845853471,0.260479521632094) q[6];
u3(1.02341772265079,-3.50771306635591,-0.404864354342131) q[0];
cx q[0],q[6];
u1(-0.371874531364530) q[6];
u3(-2.04485320077819,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.42347059330087,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.31231792335600,2.07629360383576,-1.19146868096388) q[6];
u3(1.30605039470242,-0.393221739074710,4.80754140060942) q[0];
u3(0.719552680368863,-3.41072909396089,2.69972832797289) q[1];
u3(1.15214791421335,-2.81777041516684,2.41872123279774) q[8];
cx q[8],q[1];
u1(0.889598471136980) q[1];
u3(-0.195821491598090,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.62765610778115,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.656027866485738,-2.24217560136408,3.37534462873459) q[1];
u3(0.934827042376074,0.480633809576142,-3.99628055708375) q[8];
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
