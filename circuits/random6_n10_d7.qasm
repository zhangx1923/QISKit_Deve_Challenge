OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.50965214724471,1.74243029486228,-0.749355139203119) q[3];
u3(1.80013950974462,5.14250700520822,0.990387735785991) q[6];
cx q[6],q[3];
u1(2.19526971444435) q[3];
u3(-0.154255362623491,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.04247802727913,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.57589616784430,-4.44713312346837,1.35856106011968) q[3];
u3(1.68153699682826,-0.281987685108965,-5.67870843534886) q[6];
u3(2.57230282403157,-1.06043497532805,2.47985660132907) q[4];
u3(2.61409851231419,-1.34450813429761,-0.422119291087462) q[9];
cx q[9],q[4];
u1(3.50367446711493) q[4];
u3(-4.25294295655633,0.0,0.0) q[9];
cx q[4],q[9];
u3(-0.799308627191096,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.544688685700946,0.280539425073043,3.11071465796695) q[4];
u3(1.28793493477237,0.920263264777754,-3.24119440557568) q[9];
u3(1.75712128838008,-0.871446437109072,1.59299634692016) q[7];
u3(2.38163449282143,-2.15437824019666,-0.364234980280401) q[5];
cx q[5],q[7];
u1(1.50027007251078) q[7];
u3(-0.291186414050772,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.23448060237529,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.664617350409867,-0.179068867248088,-2.40963488085580) q[7];
u3(0.639875474361620,-2.09792055004415,2.23991022930048) q[5];
u3(1.26773102924525,-0.158385309263983,1.70644755801633) q[2];
u3(1.52334036246464,-2.56320391411611,-0.613282006214649) q[0];
cx q[0],q[2];
u1(2.01021728931790) q[2];
u3(-2.47995852173061,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.0400701763927027,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.14170163325610,-3.79699458981628,1.53524507927134) q[2];
u3(1.61397424910188,-1.89824834524016,0.760896413343785) q[0];
u3(1.93833490184448,-0.826004775329402,0.0608546356368975) q[1];
u3(0.721595323014284,-2.91607046396594,-1.42764252011860) q[8];
cx q[8],q[1];
u1(0.969206705402980) q[1];
u3(-0.211636068138330,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.43940447852006,0.0,0.0) q[8];
cx q[8],q[1];
u3(3.02779303930885,-0.799371061014585,2.23725399230781) q[1];
u3(1.54920097378714,-0.300622844916666,4.72559392087638) q[8];
u3(1.94354340651216,2.18677232813514,-2.81122663380148) q[9];
u3(1.86972963037431,-3.22734037827474,2.92909764460840) q[7];
cx q[7],q[9];
u1(3.16997094957509) q[9];
u3(-1.81579971040770,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.02303323408932,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.337404595691273,0.381906833649595,0.341679894699623) q[9];
u3(1.83312017349026,2.26276583622706,0.188392779057335) q[7];
u3(1.33349813468123,2.30065466265909,-1.15711253825738) q[5];
u3(1.84249587515569,1.49939014049971,-0.919111052494903) q[1];
cx q[1],q[5];
u1(3.07152281710648) q[5];
u3(-1.87533287487304,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.836944651640010,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.64497273618394,-3.70813048556028,1.50312882642639) q[5];
u3(0.836314700916684,-4.77595473399569,-0.471665009802915) q[1];
u3(2.56750076385810,4.37103855826277,-1.64766800727157) q[3];
u3(1.01602982411464,2.99541538911336,-0.619228653447213) q[4];
cx q[4],q[3];
u1(1.96946173092822) q[3];
u3(-2.76155563715036,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.02609432696871,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.46217337594603,2.12750884581587,1.02158474072449) q[3];
u3(0.495220066150126,0.870228047933074,-1.63889556437970) q[4];
u3(2.24462501304247,-1.38538611638021,-1.04469949442766) q[6];
u3(0.652400230283144,-2.36289401932929,-1.83361144453739) q[0];
cx q[0],q[6];
u1(0.321406678389076) q[6];
u3(-1.16111588870240,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.01895217659093,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.68857343669420,-3.73477267385009,1.55824670724988) q[6];
u3(2.25721396712412,2.01447028710612,1.53186061004617) q[0];
u3(1.66784742107106,0.253264158710128,-2.27600817591921) q[2];
u3(1.24301409125296,-3.64168848849901,1.05655148907836) q[8];
cx q[8],q[2];
u1(0.478662771469030) q[2];
u3(-1.54853984317707,0.0,0.0) q[8];
cx q[2],q[8];
u3(3.03365128852295,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.72539963154538,1.85721708912345,-3.87098958706963) q[2];
u3(0.467315939654475,3.83602155071776,1.64947084947190) q[8];
u3(1.63604559760236,0.673127542413641,-3.63505997523151) q[9];
u3(2.11176374303856,2.51183880703977,-3.27579385561011) q[3];
cx q[3],q[9];
u1(3.59162353270155) q[9];
u3(-1.74956675736257,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.44280113470318,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.85757907496633,1.26579640160895,1.41266989346096) q[9];
u3(1.38525416072216,1.92682367631720,-2.66661761650517) q[3];
u3(0.867627990612070,-1.72448812983314,2.07176438437351) q[7];
u3(0.852085364296794,1.64393026186455,-3.37506008868481) q[5];
cx q[5],q[7];
u1(-0.0544930285932821) q[7];
u3(-2.61729549533255,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.14441099909665,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.04463611673798,-2.31999922929031,2.54161270248595) q[7];
u3(2.10603659451862,-2.51605092495265,-3.42795327870061) q[5];
u3(2.16378133640211,2.02932099380244,-1.96931142105380) q[1];
u3(2.56855252305724,2.28247010434807,-3.03151719965011) q[0];
cx q[0],q[1];
u1(2.82711157585968) q[1];
u3(-1.18381675445044,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.194423378137283,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.28456179493504,0.120500625961023,-1.98646146058058) q[1];
u3(1.04461927491432,4.26059322179281,0.245503722333538) q[0];
u3(1.16999072677633,1.81429425307103,-3.75632473320736) q[4];
u3(0.957268202550630,3.14852151851239,-2.51219485723422) q[2];
cx q[2],q[4];
u1(1.37700923487565) q[4];
u3(-0.957281272811435,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.330117986700653,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.33707729667023,-2.42909242579111,3.24147729368160) q[4];
u3(1.25135800305064,-1.36024582947799,2.73595217079365) q[2];
u3(1.37765200871633,0.221944491715876,0.606517138620066) q[8];
u3(1.22690476672410,-2.28976319334928,-1.02731585994897) q[6];
cx q[6],q[8];
u1(1.49838810415476) q[8];
u3(-3.37405148538969,0.0,0.0) q[6];
cx q[8],q[6];
u3(2.52272477085473,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.08645856241776,1.63248948238627,-0.618143420742213) q[8];
u3(0.531728926098579,2.80482289700836,1.11758514834600) q[6];
u3(1.95418433330909,-1.31706463587913,0.219698082079209) q[1];
u3(1.75007735241458,-2.87608148977901,-1.18604881380507) q[3];
cx q[3],q[1];
u1(3.08036834349314) q[1];
u3(-1.82471081583642,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.518681066858989,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.46166444579712,-0.688452233888920,0.884620722976419) q[1];
u3(1.94478873241536,0.290419805031635,5.23419464734349) q[3];
u3(2.66921659560812,2.21587112589134,-1.40914354723548) q[9];
u3(2.70381923488442,0.431946830613672,-4.92346779337204) q[6];
cx q[6],q[9];
u1(3.38148007653456) q[9];
u3(-3.91068105757386,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.387828720663132,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.78556116605841,0.559330469931668,0.306881617948223) q[9];
u3(1.64297271539215,-1.49555903634233,-2.88383154586518) q[6];
u3(2.07816337038076,0.402633446088976,1.89111808804427) q[5];
u3(1.36694607087502,-2.90906704785700,-2.68092429377129) q[0];
cx q[0],q[5];
u1(2.46557627797585) q[5];
u3(-1.68176379494008,0.0,0.0) q[0];
cx q[5],q[0];
u3(-0.125648615569319,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.19662590361788,0.721434699189548,-3.94732555823263) q[5];
u3(2.31305861997443,-2.55342708255115,0.944106547507913) q[0];
u3(1.35807154291747,4.39529388620327,-1.35669698007952) q[4];
u3(0.970652629664705,2.70583377725818,0.138579668513671) q[7];
cx q[7],q[4];
u1(1.89038108359969) q[4];
u3(-2.99607700122974,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.838582220671621,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.15124738861532,-0.663727920394410,3.28127961096945) q[4];
u3(2.36649479379969,-3.59882445632545,1.84279261645742) q[7];
u3(2.24545413086070,-0.0958628041755450,-0.172831274643977) q[2];
u3(1.40910621493095,0.0768951426823357,-5.27543693360437) q[8];
cx q[8],q[2];
u1(0.690465838825724) q[2];
u3(-0.162218347516167,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.84030683723331,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.76123964460612,-2.28660528059447,1.28284720977198) q[2];
u3(1.62203039864735,1.20620948962422,4.24222479494333) q[8];
u3(1.17905303126667,2.08232266270709,-3.67632242673727) q[8];
u3(1.47317027023990,2.43894367661662,-2.73277924802504) q[0];
cx q[0],q[8];
u1(1.37907731859168) q[8];
u3(-0.912673439548936,0.0,0.0) q[0];
cx q[8],q[0];
u3(3.22068944118922,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.638110109530153,-2.09808038183619,1.51796763607460) q[8];
u3(0.992207068648383,3.52157691166801,-1.68554083941050) q[0];
u3(1.05732097364611,2.49402928836126,-1.72269929884217) q[3];
u3(1.25895500805197,0.624523418615649,-2.63213116035076) q[5];
cx q[5],q[3];
u1(3.48011934304681) q[3];
u3(-0.653929219654484,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.64187955136468,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.67954427137665,3.27095484097622,-0.487189336428972) q[3];
u3(1.30913792709515,1.96684348608013,1.21887305632781) q[5];
u3(2.19551360349338,0.418919864384325,2.40084644804387) q[4];
u3(2.37775344376783,-2.51250382375102,-1.86051159245224) q[2];
cx q[2],q[4];
u1(1.96807974883017) q[4];
u3(-3.12134168863651,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.620137026353637,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.18144725514711,-2.18325980260896,3.34783984921627) q[4];
u3(1.07997559137344,-1.86231079449089,-4.36277757148315) q[2];
u3(0.451614265540686,1.00714853071883,-2.49320687413302) q[9];
u3(1.40888285548424,2.14611587379514,-3.69657686569491) q[6];
cx q[6],q[9];
u1(0.154576534897390) q[9];
u3(-0.514945787132529,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.85237292321679,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.953324442805101,-1.17338386113649,0.0975160620839446) q[9];
u3(0.579923236506337,1.84973321012236,-1.49851412048768) q[6];
u3(1.24389051417611,3.18360633048946,-2.34200219073223) q[7];
u3(0.538744405276150,2.43701753160109,-2.52525083176017) q[1];
cx q[1],q[7];
u1(2.51997596367363) q[7];
u3(-2.87228476689300,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.70456300011262,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.355958004370143,2.01662918312551,-2.40307189965362) q[7];
u3(0.896147826106684,-1.88745870888331,1.10715688831852) q[1];
u3(0.886767802784421,2.30856015427880,-3.01904902694268) q[4];
u3(1.19548445963741,-3.87279029709570,1.89628307564733) q[2];
cx q[2],q[4];
u1(1.02814517734016) q[4];
u3(-3.11489193969796,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.60841535698945,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.18681661754040,0.841466742685820,-2.58807960413215) q[4];
u3(0.806882840558782,0.163320156634629,0.823047649948251) q[2];
u3(2.61313269213664,1.25625585380712,-1.73111718647297) q[9];
u3(2.08086647938624,-4.07623649098228,1.89875341716793) q[7];
cx q[7],q[9];
u1(2.09094774720837) q[9];
u3(-1.77947237155908,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.210245281018607,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.789863660175759,1.00751335058137,-0.391661575521762) q[9];
u3(1.27432971059285,1.01260004110012,5.19198627521614) q[7];
u3(1.96917629095696,-2.94376384509198,3.05665161477853) q[8];
u3(0.734095288513716,2.95131399782834,-2.20765089247113) q[5];
cx q[5],q[8];
u1(3.51342654018460) q[8];
u3(-1.23660725292306,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.17933060017556,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.24502089260664,-1.58693970916688,1.07299690183968) q[8];
u3(0.546059103212493,2.99516454375497,-1.37368221221778) q[5];
u3(2.28366102048859,3.27483313069441,-0.207060101771442) q[6];
u3(1.79549800238921,3.09216102521314,-0.978306896959238) q[0];
cx q[0],q[6];
u1(0.693365712976142) q[6];
u3(-1.56871679179717,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.28408211691771,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.71458955836028,-1.91458851251014,1.77311835121655) q[6];
u3(2.86519011909148,0.814584207168381,-2.66636514570313) q[0];
u3(0.752710407499218,0.984297516864093,-2.53802651436537) q[1];
u3(2.19468946876316,2.16659532118366,-3.19469374570919) q[3];
cx q[3],q[1];
u1(2.24126824191649) q[1];
u3(-3.01340846189885,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.19816056855581,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.937909040684090,2.59899274110112,-1.87133429088841) q[1];
u3(2.22649893179005,-6.01233110271597,0.225077256741036) q[3];
u3(1.58313577413250,-1.20473061464882,-1.01460586965608) q[2];
u3(0.368054337903638,-4.07277925173765,-0.226345775183708) q[0];
cx q[0],q[2];
u1(2.16264806970513) q[2];
u3(-2.99087560866540,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.43107269325390,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.13656515270968,-2.08389760366560,0.727839483822893) q[2];
u3(1.55077063023813,3.66547044393227,1.98313953272429) q[0];
u3(0.631477570225879,-1.66400396947307,-0.846191273864680) q[6];
u3(1.82680839013804,-3.33641118221472,-0.244934648043935) q[4];
cx q[4],q[6];
u1(0.137857148003611) q[6];
u3(-1.02858471313300,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.50016530369068,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.876319457410902,-0.678170846703264,-0.148118660660156) q[6];
u3(0.530045451317212,-1.55390519683017,4.17691279672224) q[4];
u3(1.37507469378113,-0.0553690116140584,2.67668323754204) q[1];
u3(1.62022953403707,-0.633016310085368,-1.61743683066408) q[7];
cx q[7],q[1];
u1(1.63837509902854) q[1];
u3(-2.18677414548719,0.0,0.0) q[7];
cx q[1],q[7];
u3(3.46865066401922,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.40382724490496,-0.994932697054807,0.894916185486791) q[1];
u3(0.984126792123850,3.21050446960047,1.01059589054997) q[7];
u3(0.904251717724113,1.20556493430308,0.0633324851882894) q[3];
u3(1.68414148331121,0.328084686199825,-2.43214709192254) q[8];
cx q[8],q[3];
u1(2.64632674571644) q[3];
u3(-1.99743450016523,0.0,0.0) q[8];
cx q[3],q[8];
u3(0.0980822001269779,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.48531027040319,-1.67763541681770,1.87250159990069) q[3];
u3(1.02338586569329,-3.75211658566533,1.26839723021282) q[8];
u3(0.907809959605024,0.826013839635101,-0.875381769571632) q[9];
u3(0.782000179979363,-3.50616239469863,0.609197388576447) q[5];
cx q[5],q[9];
u1(3.51600582251589) q[9];
u3(-1.53893589515103,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.18601833737992,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.27959735405599,4.00993048865592,-2.22155935389084) q[9];
u3(0.810748373761234,1.85072586414125,0.0846220240447931) q[5];
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
