OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.07595822479576,0.362446551785178,-2.58344004621579) q[3];
u3(2.75347807436715,3.15784236214321,-1.70767833127177) q[0];
cx q[0],q[3];
u1(4.27129118498092) q[3];
u3(-3.98125387394975,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.885461951306189,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.72754354978375,-0.520488309652289,1.75228336627946) q[3];
u3(0.790063300081208,5.60136329159742,0.200630086476604) q[0];
u3(0.984894955687161,-1.45413832013148,2.07152350903829) q[11];
u3(0.204202070514318,1.81927767845639,-3.90504600324715) q[1];
cx q[1],q[11];
u1(1.65517712347910) q[11];
u3(-2.78625009120866,0.0,0.0) q[1];
cx q[11],q[1];
u3(0.159505409009109,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.66326649619827,0.298173020205345,-3.40380010446468) q[11];
u3(1.49507753893984,0.869057924957161,4.00289960012420) q[1];
u3(1.27217294832648,2.59651561656000,-2.89942761304709) q[5];
u3(1.83451413855982,2.93922956284342,-3.33696238650830) q[10];
cx q[10],q[5];
u1(2.05509030477396) q[5];
u3(0.0402622336016727,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.724975134748676,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.939010020259857,-0.725481525620737,3.69877667471444) q[5];
u3(1.87018701510859,1.19356502167054,5.06212974020732) q[10];
u3(1.80351534535525,0.926091504993999,0.896134033997972) q[2];
u3(0.0635616717277447,-1.17824795700260,-4.12669392833547) q[4];
cx q[4],q[2];
u1(2.00469568941835) q[2];
u3(-2.86683470339579,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.33509434642829,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.45192183891406,-3.83388310924966,0.691961822009828) q[2];
u3(0.510007645983225,-0.704201894692728,-2.41109826239982) q[4];
u3(2.64268332617727,2.24273426894701,-0.935782464648110) q[9];
u3(1.95405448781845,1.51036459915695,-3.75876563489414) q[7];
cx q[7],q[9];
u1(2.03851779060960) q[9];
u3(-2.84643071433928,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.898436512592622,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.40662768889100,3.04162700897828,0.236367467135260) q[9];
u3(1.94777527240083,2.53909694886399,-1.40084620498923) q[7];
u3(0.626946246892311,0.266374027343909,-1.99557648834635) q[6];
u3(0.908272894400596,2.17151129940265,-3.50814777715692) q[8];
cx q[8],q[6];
u1(1.41956008723828) q[6];
u3(0.200241594344691,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.20938059806554,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.570107767105869,3.09461308761203,-0.640604278084619) q[6];
u3(0.639140732925724,-0.863142659693088,-2.16367791847808) q[8];
u3(0.120540449004841,-0.978418885827516,1.89741541111097) q[6];
u3(1.09759999492711,0.581563518107496,-1.71001349961935) q[7];
cx q[7],q[6];
u1(0.885951234894961) q[6];
u3(-3.43148927317836,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.65938749167399,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.731720972487660,1.27140046948553,-0.642375858225138) q[6];
u3(1.11999361095502,-0.431494791302078,3.58305712588751) q[7];
u3(1.71959975200390,2.83335268943565,-1.61626886853301) q[11];
u3(2.10419766521009,0.801708697582529,-2.46689564437010) q[8];
cx q[8],q[11];
u1(1.19043044363414) q[11];
u3(-1.32399008982649,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.74280262401967,0.0,0.0) q[8];
cx q[8],q[11];
u3(0.521721243130464,0.129107412574617,-2.33010391299263) q[11];
u3(1.06768305912420,2.61816191176355,-1.58399911934129) q[8];
u3(1.41671134301703,1.17123422472390,-3.09401707582173) q[3];
u3(2.29460060425377,-3.39240782087350,2.76242948680904) q[2];
cx q[2],q[3];
u1(2.35953819219086) q[3];
u3(-2.88821727683052,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.19550330982905,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.20831224887870,-0.152827455233688,2.36760071917280) q[3];
u3(2.36999454208642,-0.680573247784764,2.30274339441850) q[2];
u3(1.11336532601185,0.352718833299508,1.47552275391954) q[9];
u3(1.57850651345482,-2.37980427823594,-1.31005278562348) q[0];
cx q[0],q[9];
u1(1.32346741642640) q[9];
u3(-0.259577493513695,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.29597023877091,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.64915823175740,-0.430212439785439,2.51418424768715) q[9];
u3(2.08720370673305,-2.65878993423833,-3.16419564817581) q[0];
u3(1.53090392886766,1.48640312127659,-1.18619234495297) q[10];
u3(1.88125452759555,-0.487836189513146,-3.45626986803851) q[4];
cx q[4],q[10];
u1(0.888746136516758) q[10];
u3(-1.15196402611542,0.0,0.0) q[4];
cx q[10],q[4];
u3(-0.306245608620265,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.21251346381664,2.68255957503512,-0.938982853368296) q[10];
u3(1.10183702570324,-3.63539766271654,1.38287246411155) q[4];
u3(2.28715218766654,0.0329900104161978,-0.287325728333718) q[5];
u3(0.581412787444661,-3.60863940570112,-0.772128181984085) q[1];
cx q[1],q[5];
u1(3.51746999920151) q[5];
u3(-4.16774412306032,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.810572844869772,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.11966322692934,2.23050505272494,-1.48276103706256) q[5];
u3(2.62656311437198,0.988345054157623,-3.74469322181735) q[1];
u3(2.74011485703436,1.10780806702184,-1.28599920398278) q[9];
u3(2.30150853887430,1.44904583257177,-2.84135925904681) q[11];
cx q[11],q[9];
u1(0.347401561482497) q[9];
u3(-1.53409014813780,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.23415437046898,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.47907893516209,-2.35149929423680,3.21143093886809) q[9];
u3(2.13981256381062,0.824480485628373,5.04555643557580) q[11];
u3(0.109578435531745,-1.27858675908320,-0.922801159812139) q[4];
u3(1.58371820420053,-5.33686174855250,0.934153079774762) q[8];
cx q[8],q[4];
u1(-0.286972577094434) q[4];
u3(-1.73124120926951,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.22558251627685,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.66113930854471,2.26935151930339,-2.82687046550104) q[4];
u3(0.417703976939094,1.50852766662772,3.78485210348151) q[8];
u3(2.25292879764655,-1.49875818849054,0.227684687861813) q[5];
u3(1.67618724358695,-3.04730750460784,-0.356457183492807) q[6];
cx q[6],q[5];
u1(1.77468108859494) q[5];
u3(-3.22512494884046,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.42348503330476,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.427065660946374,2.66823531276068,-3.40714909941297) q[5];
u3(2.41508716004028,0.0240811796004483,2.14864181314489) q[6];
u3(2.68535066159334,-2.06879614028772,1.14102473307725) q[0];
u3(2.34462513809596,-2.50173001103987,-0.499293160691616) q[2];
cx q[2],q[0];
u1(2.28780288836087) q[0];
u3(0.280653001902792,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.22482640587306,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.764207102792252,2.86907721747950,-0.832426149712933) q[0];
u3(1.42647613196256,0.599826557997011,2.89537717427988) q[2];
u3(0.199277865811645,-0.0383206764695215,-0.274937185465834) q[7];
u3(1.21408920872533,-2.65883066873973,1.96221139664302) q[3];
cx q[3],q[7];
u1(0.703254888622501) q[7];
u3(-0.0779853022406085,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.63792797588228,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.12534007157317,-2.66961834909792,1.07647256940009) q[7];
u3(0.281019751510828,-1.48027783300749,1.94473675982000) q[3];
u3(1.08448103441421,-0.272681029089181,-0.0320867037366586) q[1];
u3(1.74061327630036,-2.94735236990183,0.545896016461063) q[10];
cx q[10],q[1];
u1(3.06968696968141) q[1];
u3(-2.82203657200017,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.03139545156873,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.50232646180770,3.55424542311272,-0.639473696975219) q[1];
u3(1.92528376199766,5.29256554867334,0.669082828164727) q[10];
u3(0.474966860036681,3.47098621735158,-2.64659780227470) q[0];
u3(1.81171904373560,1.36556130014960,-2.07727373426922) q[11];
cx q[11],q[0];
u1(-1.15990650580467) q[0];
u3(0.136369511160248,0.0,0.0) q[11];
cx q[0],q[11];
u3(3.61563158928318,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.26574380831205,1.47542879136538,0.199562467077889) q[0];
u3(1.19282886984514,-1.21924730012875,-1.97306724266068) q[11];
u3(1.75538111447650,2.09292963421566,-2.66080010577833) q[8];
u3(1.56690654026358,2.43801105115784,-3.57041013139829) q[7];
cx q[7],q[8];
u1(4.23355299259826) q[8];
u3(-2.65494570299592,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.0563895692950405,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.43874080480986,-0.950330899757255,-1.48168247437067) q[8];
u3(0.270418556531249,-0.808744336253401,0.203086198985496) q[7];
u3(0.665630177208084,1.03330445491062,-0.00478371294656044) q[6];
u3(1.22566404446719,-0.776288889347728,-1.02710393228747) q[2];
cx q[2],q[6];
u1(1.92147072863628) q[6];
u3(0.408979421335025,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.846217822718107,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.17116818125056,3.30816218557306,-1.42639847812671) q[6];
u3(2.45342549049831,-2.76227244047945,2.63884716581848) q[2];
u3(1.58799671319936,2.95928646653304,-0.218672761520640) q[3];
u3(0.921004434107698,1.02655941128311,-1.15342595125676) q[4];
cx q[4],q[3];
u1(1.08758056047040) q[3];
u3(-0.405737353869708,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.09612625355573,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.80392877158224,-2.45436518601895,1.69375800492888) q[3];
u3(1.96172638340337,1.39870079923893,-0.446449895305576) q[4];
u3(1.26182796000151,0.526410404466022,-2.83850875263489) q[9];
u3(1.65358370666203,-2.85381736756539,3.39460209259694) q[1];
cx q[1],q[9];
u1(0.111879375426860) q[9];
u3(-1.43424088905511,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.78259761917081,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.60391702793020,2.72965015760419,-0.222460727311381) q[9];
u3(1.47067618797839,-0.526836569051370,-0.449389837194912) q[1];
u3(1.00754434540497,0.0383909806043696,-0.560857083791028) q[10];
u3(0.976226500688075,0.0124557162858895,-1.81949917643534) q[5];
cx q[5],q[10];
u1(4.28931778179403) q[10];
u3(-3.47757664159200,0.0,0.0) q[5];
cx q[10],q[5];
u3(-0.421416772905208,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.67368594095794,-0.990868588104315,-1.71101084655047) q[10];
u3(1.96367649931698,-5.09946319872017,-0.990719764065468) q[5];
u3(1.12079836891161,-0.142664171843921,1.42056593406350) q[1];
u3(1.18693646631334,-1.71085099750503,-0.308286033575617) q[9];
cx q[9],q[1];
u1(2.46571274656302) q[1];
u3(-1.79455831643654,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.283974182872154,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.86440195293110,-1.63314330459607,-0.678938766502966) q[1];
u3(0.907237135035368,5.88497414897504,-0.368803580212225) q[9];
u3(1.74694545653470,-0.468848397456499,0.978197224255597) q[3];
u3(1.77021840818959,-2.86787885267460,-0.307910049379587) q[10];
cx q[10],q[3];
u1(1.50071700303500) q[3];
u3(-1.03733984237786,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.72292021823213,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.08295931944951,-1.06872414707559,-2.15280786928128) q[3];
u3(1.06960964703094,0.0982968930640740,1.24606709896428) q[10];
u3(0.965043283464822,-0.123960574166864,-2.21833756556614) q[8];
u3(0.886894701487058,1.76091502839638,-3.99766772169640) q[11];
cx q[11],q[8];
u1(3.53226404513106) q[8];
u3(-1.57146024689939,0.0,0.0) q[11];
cx q[8],q[11];
u3(2.51507247143964,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.49293624369291,2.15770760435602,-2.36963555087103) q[8];
u3(0.883612771967285,1.25166596030632,0.569571919437513) q[11];
u3(1.57288829811405,-1.92387688353116,-0.243037319340725) q[7];
u3(1.60663269836862,-2.23235583875871,0.201999067299749) q[6];
cx q[6],q[7];
u1(2.07647793855661) q[7];
u3(-2.31609780180910,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.73693685009808,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.70798193973115,1.33829142883475,-0.967771937611674) q[7];
u3(2.19648667051227,0.834515589436113,-1.93117631906118) q[6];
u3(1.79604305393189,0.177033977845882,-1.43089545777016) q[2];
u3(2.41485060854290,-3.27008423443322,1.76749841478188) q[5];
cx q[5],q[2];
u1(1.12495968872740) q[2];
u3(-0.506350601509112,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.0614961667552028,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.76936698333507,3.63289857170770,-1.93156173032887) q[2];
u3(2.55917640619684,-3.72829017687471,1.32029847274730) q[5];
u3(2.15188792218698,-3.63662182245705,0.645214704856093) q[4];
u3(1.82213378940174,0.313947898916979,4.00606054670982) q[0];
cx q[0],q[4];
u1(1.63356648958192) q[4];
u3(-0.657242875155612,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.13977000997995,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.02360175589532,-1.00811225761768,1.54629137159710) q[4];
u3(0.377029374316121,-5.32435533518630,-0.770322882336932) q[0];
u3(0.997808438414469,0.718849773053188,-2.96938864540100) q[4];
u3(1.37450055023729,-2.83614667152865,3.35650099822806) q[11];
cx q[11],q[4];
u1(3.11657121331787) q[4];
u3(-0.822376365838680,0.0,0.0) q[11];
cx q[4],q[11];
u3(2.02918018653439,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.61477216824269,-3.96330386248281,0.920391707072519) q[4];
u3(2.15103528333472,0.856661853059416,2.16721447053082) q[11];
u3(1.16417559479415,1.14181508545626,-2.82285949035367) q[7];
u3(1.19082788644077,2.66485203033123,-3.43556834620428) q[10];
cx q[10],q[7];
u1(3.16927606263903) q[7];
u3(-2.01856330173750,0.0,0.0) q[10];
cx q[7],q[10];
u3(1.42430926244422,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.80325160156609,2.73586159906249,-2.21703057453486) q[7];
u3(2.02288847434656,3.43617723240131,0.949668907208912) q[10];
u3(1.14951019305239,-1.11237395637655,0.184085328845886) q[3];
u3(1.40922557665086,-3.58583381721483,-1.22809982044825) q[1];
cx q[1],q[3];
u1(0.224364082537447) q[3];
u3(-1.24617983757167,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.65300751790136,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.71678147403324,0.853843952398519,-0.0840804871854170) q[3];
u3(0.270529959096525,-0.418966525021744,2.41780266164473) q[1];
u3(0.774119748456378,-0.159865489877104,-1.58285866649896) q[5];
u3(1.29072726203422,1.95791190509672,-3.66584573978170) q[8];
cx q[8],q[5];
u1(2.94481417486717) q[5];
u3(-2.09309022652229,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.483986883690146,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.35143738300162,0.723594757605533,2.81074692882139) q[5];
u3(1.00567570071800,3.87790758542580,-0.272036077999133) q[8];
u3(1.59511154999599,0.630734531877862,-2.39202994349806) q[9];
u3(2.52785790686340,-3.53681935754734,2.20599706736253) q[0];
cx q[0],q[9];
u1(2.12509041296891) q[9];
u3(-2.72274890550883,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.32539225620978,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.409964917021315,-0.610008220600030,-2.02812537038774) q[9];
u3(0.948054207391072,-1.70243945459592,0.993943141308875) q[0];
u3(0.637846288701748,1.94036908394864,-2.77717461683232) q[6];
u3(0.923628512178863,-3.53628371024867,2.72483906646598) q[2];
cx q[2],q[6];
u1(0.0916264051774147) q[6];
u3(-2.08697269595719,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.802638092165504,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.70327447855363,1.55060945560623,2.00510041210880) q[6];
u3(1.24886994319916,-3.97629885606444,0.402658567971515) q[2];
u3(2.24092511438121,-1.58502302464641,2.51547195806829) q[10];
u3(2.66488057596152,-2.26157152001217,-0.817105966473399) q[9];
cx q[9],q[10];
u1(1.71039124077153) q[10];
u3(-2.90069737280632,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.855135513249352,0.0,0.0) q[9];
cx q[9],q[10];
u3(0.507091925401420,-0.227283139765260,1.78253795721093) q[10];
u3(2.35998224407466,-2.33869144234496,1.59235796137024) q[9];
u3(0.829775956551375,-1.06117578357430,2.27468378317931) q[4];
u3(0.277742188225751,1.25380859348281,-2.91056140879519) q[8];
cx q[8],q[4];
u1(2.28640249073691) q[4];
u3(-1.50927334686923,0.0,0.0) q[8];
cx q[4],q[8];
u3(3.78686139788750,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.11004265690850,-0.986895619261183,3.12277381707002) q[4];
u3(2.06298189590532,0.329282624492394,4.89505143702242) q[8];
u3(1.19859442519426,-1.27024035417034,-0.828802429825968) q[3];
u3(1.95199152463877,-2.58707568586759,-0.587303719579588) q[0];
cx q[0],q[3];
u1(1.40392507466307) q[3];
u3(-0.999551822391517,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.70654533998493,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.11367218142261,1.66279401868685,-2.09502448355700) q[3];
u3(1.39059983778041,-0.859517467132820,-5.40066097351244) q[0];
u3(1.16741074720885,1.69556089845619,-1.25539934467541) q[1];
u3(0.162037737244661,0.482848928546051,-1.04446177972786) q[11];
cx q[11],q[1];
u1(1.43035682316388) q[1];
u3(-0.837542080741425,0.0,0.0) q[11];
cx q[1],q[11];
u3(2.01217936536458,0.0,0.0) q[11];
cx q[11],q[1];
u3(2.19409830905821,-2.12529312318167,4.03631516469620) q[1];
u3(0.527708840982685,1.86537093441384,-0.0230982719909026) q[11];
u3(0.199999864113970,-0.780844753545010,1.93459118194785) q[5];
u3(0.374497650217778,-2.63921018617231,0.993920416784621) q[2];
cx q[2],q[5];
u1(2.56968772213675) q[5];
u3(-2.88512868998878,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.08645021334141,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.05273752449661,-1.35009409927534,2.82358535467326) q[5];
u3(2.03816175600127,1.86789971917614,0.996398249788632) q[2];
u3(1.35999738449098,0.988219968907246,-2.26834145715001) q[7];
u3(1.70615270307687,2.03344697884250,-3.85611357571170) q[6];
cx q[6],q[7];
u1(0.182822292381968) q[7];
u3(-1.33542896272978,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.45287894215044,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.50796928218280,-0.711045500119459,4.53714191049335) q[7];
u3(1.54903060558761,0.155175460425712,-5.12197184920742) q[6];
u3(1.15904884293269,1.67896681673520,-1.48265726855226) q[8];
u3(0.482321018132229,1.82659172899767,-3.66284847166873) q[6];
cx q[6],q[8];
u1(2.29267577963482) q[8];
u3(-3.00040697094086,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.23776398988232,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.23243394151722,-0.391720808139513,4.80994986036311) q[8];
u3(2.35593381518651,0.00592013929516322,1.07579520900640) q[6];
u3(1.49111955680980,1.75668352507788,-0.898456858134390) q[3];
u3(2.10422449741314,-0.441657085962482,-3.14606164826908) q[2];
cx q[2],q[3];
u1(-0.122833281275034) q[3];
u3(0.326680749590136,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.80299043618774,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.16162468626492,-1.91672679041094,3.38192614727129) q[3];
u3(1.64616247200972,-3.45745243385328,-1.39311248164468) q[2];
u3(2.65138855346234,-2.58076164963624,0.281256693052917) q[5];
u3(2.29164448080827,1.38920784419876,2.46838709399254) q[10];
cx q[10],q[5];
u1(1.34001680999655) q[5];
u3(-0.0262175158873268,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.67438886020647,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.03819692550094,-4.00931884177090,-0.205378336380625) q[5];
u3(2.73434595708544,5.27671825428047,-0.110984213189001) q[10];
u3(0.859070595295487,-2.96325697992439,2.57568875305502) q[11];
u3(1.42023553753980,0.179874259995054,-2.31706486368778) q[0];
cx q[0],q[11];
u1(2.64817964782337) q[11];
u3(-3.07619208844462,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.59906760226930,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.46779579854795,-2.50462360286744,0.974285646864046) q[11];
u3(1.33697423290549,-1.08645415451335,-3.70699425428801) q[0];
u3(1.64328203599361,-0.159528442860429,1.35610186233707) q[4];
u3(1.53766283638821,-1.78060690734055,-2.50015541759127) q[9];
cx q[9],q[4];
u1(2.83659101929693) q[4];
u3(-1.54642328007962,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.978135130892243,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.65418542530984,-0.296124439427755,-1.49583319801727) q[4];
u3(0.431540097135666,1.27047581398142,-2.88837719635993) q[9];
u3(1.85172993040664,-1.08160747770401,3.56466535281413) q[1];
u3(1.51594021439517,2.09279840547537,1.45667358267791) q[7];
cx q[7],q[1];
u1(1.78334827050024) q[1];
u3(-2.40646804205828,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.526277015600991,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.80378878475447,-3.69186178988299,0.379309457500942) q[1];
u3(1.42972718360656,0.129388146833471,-0.599568810608265) q[7];
u3(1.07840265616788,-0.566273579716116,1.28863731817961) q[8];
u3(0.818109346817713,-1.47654443706411,-2.38028521771221) q[2];
cx q[2],q[8];
u1(-0.151491959371030) q[8];
u3(0.693656253472241,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.98364348803546,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.71508881416479,-1.83544029885710,3.69410480057343) q[8];
u3(2.20404368677270,-4.46292142658639,1.55350863858185) q[2];
u3(0.671483834049293,0.856680057387277,-1.36581690814074) q[5];
u3(0.513182129795738,-1.53421490719120,-0.608575491862148) q[1];
cx q[1],q[5];
u1(1.46325190498239) q[5];
u3(-0.646672444999028,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.15147093499312,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.13163812774641,-0.329334798260036,-2.34800496794068) q[5];
u3(1.46375335714250,-0.645675100812473,4.45880043396997) q[1];
u3(0.762722957278021,1.37727082139977,-0.744882665834373) q[0];
u3(0.623807611531091,-0.127152471145647,-1.12939404671804) q[4];
cx q[4],q[0];
u1(2.90624219397588) q[0];
u3(-1.66529147388301,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.627144496119007,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.22476148406648,1.19317282636899,0.0399449147980502) q[0];
u3(0.901333938458846,-1.80879983736057,3.82056994699713) q[4];
u3(1.50855582575244,0.0796125589866989,2.08543015941881) q[9];
u3(1.31299081758263,-1.18456673762034,-1.94455323424130) q[10];
cx q[10],q[9];
u1(3.63797994638201) q[9];
u3(-2.16877327806701,0.0,0.0) q[10];
cx q[9],q[10];
u3(-0.0181718417194985,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.754859121890678,0.990247515573494,-1.69464327947465) q[9];
u3(1.21282757930931,2.41428510329509,3.22388810619046) q[10];
u3(2.32695345810211,-1.47161683710086,-0.913817551827691) q[3];
u3(0.897415444344687,-2.57328975380950,-0.620381519108979) q[7];
cx q[7],q[3];
u1(2.27131901021787) q[3];
u3(-2.93645263139535,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.481308872388793,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.49239550247663,-1.28946182641857,2.78513815027897) q[3];
u3(0.742173861273085,0.813971941750105,-1.68118449769453) q[7];
u3(1.77473711436475,2.19361707355254,-3.31169386356429) q[11];
u3(2.12252367135290,-3.24951484717180,2.92539640244026) q[6];
cx q[6],q[11];
u1(3.34248346125544) q[11];
u3(-1.52022245362094,0.0,0.0) q[6];
cx q[11],q[6];
u3(2.62608195950711,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.58111179459857,-1.75323892889243,0.566504413175245) q[11];
u3(0.813215184743642,1.12735001280222,4.78244228575421) q[6];
u3(0.785567247691280,0.771963957296634,-0.407007299663128) q[9];
u3(1.37589051904438,0.128728845389089,-2.69246908039118) q[5];
cx q[5],q[9];
u1(-1.27605481503596) q[9];
u3(-0.350014207702404,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.97446651584608,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.44015316234933,2.28961190728396,-3.91498274542557) q[9];
u3(1.62123832273604,0.950415294990965,-3.32023961057482) q[5];
u3(1.64464140529309,1.73483611217159,-0.928405055529164) q[8];
u3(1.74277533441334,-1.08891903048761,-3.41862407253529) q[0];
cx q[0],q[8];
u1(2.98834440390973) q[8];
u3(-1.78419286788538,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.431222694125435,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.917266904169453,-1.87943296927749,1.20621203933564) q[8];
u3(0.994514277514789,-0.631825884687082,-1.09684766462289) q[0];
u3(0.743059555156433,0.560000970901284,-1.13544198196206) q[6];
u3(0.938219438869527,-0.379760927945427,-0.208862255937367) q[11];
cx q[11],q[6];
u1(1.91253447570005) q[6];
u3(0.0873717953203672,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.56362798794846,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.727317452548325,-2.93351317752807,1.15981192217033) q[6];
u3(1.26920401946497,3.06103812504309,1.98029039358257) q[11];
u3(1.38561394055385,-2.36311926285013,-0.510890337188103) q[3];
u3(1.48283728229986,-4.30346804455165,-1.64408691317584) q[2];
cx q[2],q[3];
u1(1.63532782956439) q[3];
u3(-3.20668928990931,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.64876814646648,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.575215129463209,-3.71006313178914,0.703875047470052) q[3];
u3(1.68672054199791,3.76572037438238,-1.04442433673084) q[2];
u3(1.02872070061908,0.548372386731224,-0.433224428835326) q[1];
u3(1.14727510390397,-0.939140307471300,-1.52960581533570) q[4];
cx q[4],q[1];
u1(3.32356653983260) q[1];
u3(-0.627514969574082,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.83903900683250,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.77212495710586,2.47622374714245,-0.135967164422265) q[1];
u3(0.935989528923834,0.922440473699619,2.87286813414538) q[4];
u3(1.47920041421594,1.74916711422648,-0.152667214737198) q[10];
u3(0.981075609059426,0.784446415809310,-4.47542651133295) q[7];
cx q[7],q[10];
u1(2.22742046994076) q[10];
u3(-1.84498145382559,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.94770749094306,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.771612319099026,-0.947927298706649,3.34364850167991) q[10];
u3(1.76578319520025,-0.272504671377812,1.79352321154734) q[7];
u3(1.55884643091282,2.93550991815421,-1.57559700761304) q[0];
u3(2.21571766023911,0.893703773151473,-2.83909904565001) q[2];
cx q[2],q[0];
u1(2.64159432769676) q[0];
u3(-2.45584588095257,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.74457329940703,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.36029018298596,-1.86266641285422,2.64669279059539) q[0];
u3(2.37403901153216,-1.24356042850268,-3.73614080342264) q[2];
u3(1.44530131275232,1.50830308944710,0.620909236819071) q[7];
u3(1.32473312274481,0.274644720187751,-3.38428311578632) q[6];
cx q[6],q[7];
u1(1.29795098235954) q[7];
u3(-0.212457244785746,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.63291778140320,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.824359289001706,-1.09000169947953,4.12575847022570) q[7];
u3(1.25262716152255,0.0323818709860340,5.57974149643843) q[6];
u3(1.58387458657931,-1.58254069251930,0.462895005890530) q[1];
u3(0.947407928836398,-2.23808192602720,0.380734409069397) q[8];
cx q[8],q[1];
u1(1.07923767619270) q[1];
u3(-0.134109484274783,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.52910752956849,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.65475396368727,1.18738526374474,0.778569809132964) q[1];
u3(2.12461466570547,4.44312580871343,-0.192216169415220) q[8];
u3(2.15939110543959,-2.31727243939926,0.0792262126316818) q[4];
u3(2.00465647793641,-3.95903855342796,-1.32148095470245) q[3];
cx q[3],q[4];
u1(0.497254254858383) q[4];
u3(-0.290550199363185,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.57430037997687,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.70318591147250,-2.97294220788033,3.20387580820314) q[4];
u3(2.04462327805523,-1.15516839116571,3.63346311544609) q[3];
u3(1.34739273670057,-3.16800001275280,1.20408444605989) q[5];
u3(1.37020806297670,-2.98050093725875,0.191599897324103) q[11];
cx q[11],q[5];
u1(2.85662719593624) q[5];
u3(-1.55205868602121,0.0,0.0) q[11];
cx q[5],q[11];
u3(0.891172341134362,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.11067917468565,2.05822856440282,-1.46060502331720) q[5];
u3(3.10762231291184,1.55243984513044,-3.91955167040724) q[11];
u3(1.40139358372603,0.367985005106925,1.32747900929070) q[9];
u3(1.61319580333770,-1.06657448951498,-2.15551681995631) q[10];
cx q[10],q[9];
u1(1.52710756329275) q[9];
u3(-3.20065195372242,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.260127048111831,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.06455876091297,-3.06837194766372,2.38541014559550) q[9];
u3(2.00277850023555,2.11563091876765,-3.77881860215956) q[10];
u3(1.89579164986558,1.82435869857899,0.369254732117069) q[3];
u3(1.12232179560301,-1.01653344134389,-2.48571680452376) q[7];
cx q[7],q[3];
u1(3.11533451327083) q[3];
u3(-2.16634570696673,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.758478110200850,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.70039217134183,-1.89864630871279,3.99812084699348) q[3];
u3(0.346912182827202,-0.582856424942991,-2.39264227366156) q[7];
u3(2.32437515674697,-0.392687592766005,1.89839344352474) q[8];
u3(2.09638537022153,-1.37690732227615,-1.44446883936582) q[0];
cx q[0],q[8];
u1(0.526029243179842) q[8];
u3(-3.37366022999211,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.46293644406277,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.86332330156525,-0.516414524267211,2.54112671287627) q[8];
u3(0.747416023064932,1.04564186431334,-0.480230710632285) q[0];
u3(1.35239493782051,0.298827112339573,-1.46020139490011) q[4];
u3(2.24656489193909,1.19674614742941,-4.31955483822508) q[2];
cx q[2],q[4];
u1(3.63248672631208) q[4];
u3(-1.35473425300627,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.31812427222559,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.46736014210238,-3.24194268504716,2.49735491803402) q[4];
u3(1.83858153190856,0.473490916151198,1.88478578045265) q[2];
u3(1.83519825211363,4.07378861433519,-1.21300125072954) q[5];
u3(2.38153257030129,2.15644635865531,-0.477020829748298) q[10];
cx q[10],q[5];
u1(3.35403512378255) q[5];
u3(-1.37850655985592,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.19474970747088,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.78061090358968,-0.652570386961990,2.11794216665202) q[5];
u3(1.77841013035241,-4.66472078894618,0.130552408496767) q[10];
u3(2.03828464204765,2.65461133453499,-2.23253152484878) q[6];
u3(1.86619262898775,1.59121053222347,-1.62055539067751) q[9];
cx q[9],q[6];
u1(3.07469181226595) q[6];
u3(-1.31382680552053,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.53552580817599,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.69865523480312,-4.52682087146198,0.151174946034931) q[6];
u3(2.03981474873327,-1.55269081690812,3.08143037224909) q[9];
u3(1.38512120063138,-0.190066107593831,2.11322639563654) q[11];
u3(1.40698993107804,-1.01202889054626,-1.69603184737588) q[1];
cx q[1],q[11];
u1(1.49871223378332) q[11];
u3(-0.0348273297600932,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.94826807591351,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.02552120844678,-2.75304822109847,2.59649983627182) q[11];
u3(1.45819100508288,3.81992517573606,0.254417245791505) q[1];
u3(1.70282829309959,-0.381971768235166,-1.70622658754381) q[2];
u3(1.45170331972159,0.410941809797349,-4.42703297012218) q[11];
cx q[11],q[2];
u1(2.78435637602850) q[2];
u3(-2.37483058796065,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.553156398787018,0.0,0.0) q[11];
cx q[11],q[2];
u3(2.01654911400035,-1.39135704704056,1.75408348908235) q[2];
u3(0.495016878232795,0.636734019688213,-2.38350659972852) q[11];
u3(1.43111371857349,1.00352282629188,1.58706223399852) q[7];
u3(1.37206215138387,-1.76470655193764,-2.52717418385207) q[10];
cx q[10],q[7];
u1(1.73595995026113) q[7];
u3(-3.13394067071924,0.0,0.0) q[10];
cx q[7],q[10];
u3(1.02496489155110,0.0,0.0) q[10];
cx q[10],q[7];
u3(2.64861133014901,-4.30779890310281,0.569335883271925) q[7];
u3(0.952739391516115,-1.50687633156472,4.65134759937532) q[10];
u3(1.63293568580681,0.0980122545687957,1.65414901988991) q[6];
u3(1.78997696818003,-0.591039056931011,-0.875914933699837) q[9];
cx q[9],q[6];
u1(1.79139680600025) q[6];
u3(-3.19749987855554,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.343614071092826,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.50072921687615,3.49985010684817,-1.87439159907778) q[6];
u3(0.452428276103066,0.627752125041421,-5.41865914816521) q[9];
u3(1.81123497075670,0.798347733977552,1.78207966388007) q[3];
u3(1.42512774420471,-2.18265141054580,-1.99965671071650) q[5];
cx q[5],q[3];
u1(1.75883368764576) q[3];
u3(-2.95894777846365,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.55529935800608,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.66552790887804,0.640874263911966,-2.39556854822944) q[3];
u3(2.87239730515848,0.230912371256292,3.12300696469945) q[5];
u3(1.74113240198114,-0.127210636118393,1.82727593220782) q[8];
u3(2.33828184634184,-1.62654176347410,-0.601945505160649) q[0];
cx q[0],q[8];
u1(1.60125223823032) q[8];
u3(0.0550206196211249,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.04081553935998,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.23469508208609,-1.16193343050545,-0.872637830549851) q[8];
u3(1.45680613854510,-0.649081540182863,-2.76254281054479) q[0];
u3(2.12707063951451,1.20640661369917,0.849863544268326) q[1];
u3(0.958489511553133,-4.90254821199338,-0.436671927785483) q[4];
cx q[4],q[1];
u1(3.64227884735406) q[1];
u3(-1.38094998230900,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.12244151534374,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.26720137653282,-1.07095053514253,-0.0344757588259313) q[1];
u3(0.506148869235446,2.13249195567608,-2.02750940534642) q[4];
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
