OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0.952600115759848,0.318758712142087,2.41996479243504) q[2];
u3(0.913732787467228,-1.66628073128335,-1.04474059247331) q[0];
cx q[0],q[2];
u1(1.08716169779791) q[2];
u3(-0.566561694274123,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.132721684635262,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.27043443418575,-1.31238855897045,-2.40844297761311) q[2];
u3(0.595432310751500,-1.41960344245049,-1.13295007953749) q[0];
u3(1.23042693890356,0.169276463827087,0.351646150873296) q[10];
u3(1.35669509681019,-1.77974446898452,-1.64587699440107) q[13];
cx q[13],q[10];
u1(0.844242975252488) q[10];
u3(-1.07443467534163,0.0,0.0) q[13];
cx q[10],q[13];
u3(1.86482603443121,0.0,0.0) q[13];
cx q[13],q[10];
u3(1.55704031178152,-3.14938156422111,1.56029544842794) q[10];
u3(2.40449051518541,-4.49426905882622,0.0665572825598595) q[13];
u3(2.04795952131688,-2.59415973823115,2.91350464476469) q[5];
u3(0.133801272552869,1.89598933501347,-0.231847032120518) q[11];
cx q[11],q[5];
u1(2.64180617268562) q[5];
u3(-1.97276101242578,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.32258369046144,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.76404908324695,1.17385562950812,-3.33697859282919) q[5];
u3(1.96480160209058,-1.76175428124894,-0.288892781711331) q[11];
u3(1.78033305274946,-0.183319940269612,0.358830276888195) q[4];
u3(1.48744200116235,-2.95367241202297,-1.30574939420613) q[6];
cx q[6],q[4];
u1(1.35723204639906) q[4];
u3(-0.692985523079999,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.88970847847053,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.14887231596478,0.919578010569696,-3.13091968426402) q[4];
u3(0.903897218150299,-4.45615882235589,0.393653567381963) q[6];
u3(1.08188826114090,-0.508550869959953,1.52368102239234) q[1];
u3(0.253344419636558,-0.486537545297982,-1.29137750295972) q[7];
cx q[7],q[1];
u1(1.43973616935028) q[1];
u3(-0.361039335129918,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.11611677854784,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.116988435117024,-1.97001467214559,-1.73978635208719) q[1];
u3(0.576392221246751,1.23979312371985,-3.77200969572520) q[7];
u3(1.71847228547936,-1.51955260939886,-1.45321609713252) q[8];
u3(1.30101622062314,-4.12066447001599,0.174288145009313) q[9];
cx q[9],q[8];
u1(1.47246316941475) q[8];
u3(-3.59705166241886,0.0,0.0) q[9];
cx q[8],q[9];
u3(2.34142876496573,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.66371855550491,0.0174133375488036,0.751496669436298) q[8];
u3(2.05446862617677,-1.69936789835861,0.496678680646441) q[9];
u3(0.593809593586495,-0.242216420272996,1.04613366073411) q[15];
u3(0.769341416111544,-1.91309323530766,0.898863972919464) q[12];
cx q[12],q[15];
u1(4.42150224922227) q[15];
u3(-3.92460615393692,0.0,0.0) q[12];
cx q[15],q[12];
u3(-0.527813681732285,0.0,0.0) q[12];
cx q[12],q[15];
u3(0.962745506340889,1.97665126403352,-0.0419637375784530) q[15];
u3(2.73844669301413,-2.75067943681306,0.752681739492714) q[12];
u3(1.99042127308257,3.01827432427128,-2.19075707252183) q[14];
u3(0.788826335646171,1.03475941502072,-1.55645969138227) q[3];
cx q[3],q[14];
u1(1.65543015349934) q[14];
u3(-2.42912760590648,0.0,0.0) q[3];
cx q[14],q[3];
u3(0.863303846718705,0.0,0.0) q[3];
cx q[3],q[14];
u3(2.06594688040598,-1.04562537502725,-1.15315590003848) q[14];
u3(1.25732930549637,-3.08541611621320,1.46305919404775) q[3];
u3(2.73115908286289,1.34029360338871,-1.11841180790121) q[12];
u3(1.64585791109810,5.08803764196711,0.0195349357969641) q[5];
cx q[5],q[12];
u1(3.49123312592319) q[12];
u3(-1.39778530142508,0.0,0.0) q[5];
cx q[12],q[5];
u3(2.11465396038978,0.0,0.0) q[5];
cx q[5],q[12];
u3(1.08898680907057,-4.29251525264169,1.14174774368015) q[12];
u3(2.63182598842770,-0.419472939945073,-3.75853946074822) q[5];
u3(1.17012531140351,1.50566565323665,-0.159134647144300) q[6];
u3(1.08994236371738,0.586224812803951,-4.51600660172865) q[7];
cx q[7],q[6];
u1(1.65489683282714) q[6];
u3(-2.59117801273583,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.46267695028717,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.18729029243492,-1.52668413463213,3.19082235525556) q[6];
u3(0.625995973874529,-1.53042169552166,3.05716256497185) q[7];
u3(0.651145878589656,0.181180251295869,0.918894916332694) q[4];
u3(0.629309796958810,-1.11890913469463,-0.845799441038690) q[15];
cx q[15],q[4];
u1(2.46252820820408) q[4];
u3(-2.05547839879163,0.0,0.0) q[15];
cx q[4],q[15];
u3(0.313261495197385,0.0,0.0) q[15];
cx q[15],q[4];
u3(2.08600589699375,3.48621049099511,0.0733755261401630) q[4];
u3(1.14899143145694,-2.73805764746451,0.807575040116050) q[15];
u3(1.77381732300453,1.82030032845764,-3.94184123048143) q[3];
u3(2.26385834061928,-2.51638265970990,3.62841373376252) q[13];
cx q[13],q[3];
u1(1.54506693168897) q[3];
u3(-0.927106689921759,0.0,0.0) q[13];
cx q[3],q[13];
u3(2.75240955509552,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.32140585557710,-1.69319613284862,2.77750347297968) q[3];
u3(2.29020371811749,0.200619120761221,1.93914987490213) q[13];
u3(1.89756700234807,0.556237432574587,2.36454874991459) q[8];
u3(1.96936750442007,-1.27285995292304,-1.05861547292497) q[0];
cx q[0],q[8];
u1(1.69203256215408) q[8];
u3(0.201363924316359,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.629016431630138,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.920532747146566,0.753802494050862,0.230972945104374) q[8];
u3(2.38823472529532,3.96703642550830,-1.71458612023748) q[0];
u3(1.72905514157752,1.33352926385852,-2.64540793057751) q[9];
u3(1.03297165875496,2.91056719776673,-3.03360825739379) q[10];
cx q[10],q[9];
u1(0.235967940191799) q[9];
u3(-0.975684749460104,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.01079274969025,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.79720957937748,1.72535830102376,-1.91430726012013) q[9];
u3(2.88066079428274,3.11956119376218,1.24570026779082) q[10];
u3(1.95874505889170,2.12791726250733,-3.91263490500634) q[1];
u3(0.185847494434205,1.02233411725954,0.548912126749666) q[11];
cx q[11],q[1];
u1(1.37058646395138) q[1];
u3(-0.252572041093863,0.0,0.0) q[11];
cx q[1],q[11];
u3(2.08283347767777,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.79211805924172,-2.51679169779738,1.16942743521967) q[1];
u3(1.50891840801169,-2.27155193378698,-1.27551242183740) q[11];
u3(1.92035048564325,-2.84460226751039,2.76958498715188) q[14];
u3(0.891653223355963,3.23859146701461,-1.76396486006862) q[2];
cx q[2],q[14];
u1(1.89045783915982) q[14];
u3(-2.75059562427868,0.0,0.0) q[2];
cx q[14],q[2];
u3(1.15572191422554,0.0,0.0) q[2];
cx q[2],q[14];
u3(1.81526901014353,-0.908044539066524,1.97415414136631) q[14];
u3(2.43857395958299,1.17924475431142,-4.86100457920192) q[2];
u3(1.48121918976783,0.0611757523926308,2.06367768605965) q[4];
u3(0.769186490810409,-0.265386270908230,-1.10485076126035) q[0];
cx q[0],q[4];
u1(0.714582628465545) q[4];
u3(-1.14248628003623,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.0465175738444175,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.00567512595880,-0.113143251143404,1.50053454595867) q[4];
u3(2.20914842128975,-5.12641631104323,-0.420625957718383) q[0];
u3(2.01834432898458,0.376457474697125,0.536185734968087) q[10];
u3(1.04459655790777,-2.70993071713027,-0.883058546496036) q[14];
cx q[14],q[10];
u1(1.63193316569728) q[10];
u3(-2.88138700932512,0.0,0.0) q[14];
cx q[10],q[14];
u3(0.551625999521582,0.0,0.0) q[14];
cx q[14],q[10];
u3(0.530992719541874,3.32186668541779,-0.486094454880171) q[10];
u3(1.99858949751899,1.26431905490768,-0.162404521275950) q[14];
u3(2.19805328486193,3.55353002414530,-1.90144300961047) q[6];
u3(0.840125593255450,1.04920105709080,0.214395606742289) q[3];
cx q[3],q[6];
u1(2.89562817857527) q[6];
u3(-1.18491305420776,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.54998248950303,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.14744586860664,3.49604489502891,0.410124293508089) q[6];
u3(1.78123835953304,0.327399367558587,1.09272058376737) q[3];
u3(1.52333478781124,0.352488779966032,2.16520418510275) q[7];
u3(1.15634513775179,-0.841961877494668,-1.62597923296649) q[1];
cx q[1],q[7];
u1(0.742750162380131) q[7];
u3(-1.29798955790731,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.277418105443593,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.989847556257003,1.07616090214669,-2.74780324677361) q[7];
u3(1.37647596294331,2.06639596835608,2.42020299972313) q[1];
u3(2.98650423197403,2.71707804866868,-1.51009801678375) q[13];
u3(2.85958366187556,0.878935220967128,-4.23270728292752) q[2];
cx q[2],q[13];
u1(0.887192951302090) q[13];
u3(0.0279380574666324,0.0,0.0) q[2];
cx q[13],q[2];
u3(2.62110682503908,0.0,0.0) q[2];
cx q[2],q[13];
u3(2.33820057004108,0.821036077042357,3.15846777989534) q[13];
u3(1.11373113754596,-2.01100943950506,-3.40114849547101) q[2];
u3(1.22499228943550,0.462185499994828,0.907331171569201) q[5];
u3(0.777366949586559,-0.523875071170839,-3.42184337634660) q[8];
cx q[8],q[5];
u1(-0.505952470166849) q[5];
u3(-1.74456657427716,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.14529068487235,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.381881206582471,2.47665224122941,-0.483467528792668) q[5];
u3(1.01917782800772,-1.83127570298842,-3.29085396012889) q[8];
u3(1.03543640215187,0.840470917553225,0.825599107295805) q[15];
u3(1.41122366483764,0.295887090561867,-2.37274939190231) q[12];
cx q[12],q[15];
u1(2.69817526761841) q[15];
u3(-1.93069092573357,0.0,0.0) q[12];
cx q[15],q[12];
u3(0.806434366073277,0.0,0.0) q[12];
cx q[12],q[15];
u3(2.10983943690943,0.216553863876585,-1.13922093667074) q[15];
u3(1.08744004945417,2.49179924615367,1.63229150578151) q[12];
u3(2.49744622958943,1.57421976375911,-4.36985200689262) q[11];
u3(0.288666494076798,-1.87209304393866,3.32762781504992) q[9];
cx q[9],q[11];
u1(3.27086635433626) q[11];
u3(-1.89384482871773,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.579983064588701,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.68352374170562,-0.228861727850205,-3.09675872379846) q[11];
u3(2.31889443211333,0.00527135098977471,-1.16516347329747) q[9];
u3(1.63357194363666,0.269071805981965,1.30676082284160) q[9];
u3(1.80327112897008,-2.20223524516959,-0.687558671998874) q[2];
cx q[2],q[9];
u1(2.05581905850661) q[9];
u3(-2.97013468977617,0.0,0.0) q[2];
cx q[9],q[2];
u3(0.855270922773022,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.26661749002055,2.41459459168689,1.06745762093537) q[9];
u3(2.87066552477033,-0.203071725489517,-5.57581668678158) q[2];
u3(0.483844992009839,-3.11253731194531,2.20644210177569) q[3];
u3(0.912478766068660,2.48827967563270,-3.73184536980142) q[1];
cx q[1],q[3];
u1(1.59281697943526) q[3];
u3(-0.319673366046433,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.23559711479253,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.81115563966004,-1.64314421757752,3.68855025508786) q[3];
u3(1.37668988463927,0.124145607136219,-5.58340611954110) q[1];
u3(0.855104843523981,-0.477762797795731,-1.30422190800578) q[8];
u3(2.04574569931645,0.968925126884652,-4.62460856104764) q[5];
cx q[5],q[8];
u1(0.699345702356010) q[8];
u3(-0.991974166463508,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.72690163322275,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.03226545279983,-3.01146646997065,2.01537451431880) q[8];
u3(1.22992544298796,3.92131805059036,-1.96905848299981) q[5];
u3(1.64221583685567,-0.525416728454909,3.51475395460946) q[6];
u3(1.07375313165861,1.87338610763801,1.78236659364032) q[11];
cx q[11],q[6];
u1(1.16324309373823) q[6];
u3(-0.238080280707536,0.0,0.0) q[11];
cx q[6],q[11];
u3(2.19489771482821,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.10255536973042,2.00910806207588,-1.43292419628560) q[6];
u3(1.62558741286091,2.75433428241141,-2.17884508235187) q[11];
u3(1.26166330203264,0.920564778195360,-2.75955037135403) q[4];
u3(2.71261582702885,2.93997661415402,-3.19739484854806) q[7];
cx q[7],q[4];
u1(1.68197316411258) q[4];
u3(-2.04567835271310,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.0323935021551827,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.50073966212820,0.624114822001190,-3.19730535399420) q[4];
u3(1.62356590337460,1.87640076519242,-1.37574292227010) q[7];
u3(1.39701819407564,1.13762103549651,-2.68856435843186) q[10];
u3(0.651020056551575,1.96955559912720,-3.29807813521725) q[13];
cx q[13],q[10];
u1(1.87199448850365) q[10];
u3(-2.54685997047681,0.0,0.0) q[13];
cx q[10],q[13];
u3(-0.194003622212802,0.0,0.0) q[13];
cx q[13],q[10];
u3(1.48435162185000,2.83722779307744,-3.14448556702906) q[10];
u3(1.16073319580941,-2.80960178749090,-0.578393230358841) q[13];
u3(2.19459016616569,3.59796339009235,-1.80050546179260) q[0];
u3(1.46645135828049,2.16870468060768,-1.43981115493372) q[14];
cx q[14],q[0];
u1(2.10709295676564) q[0];
u3(-2.68843997915959,0.0,0.0) q[14];
cx q[0],q[14];
u3(0.955495711712023,0.0,0.0) q[14];
cx q[14],q[0];
u3(2.19550295973210,-0.962053509778821,-2.65389561110101) q[0];
u3(1.60254923120925,0.195987997009117,-3.15603043327770) q[14];
u3(1.62638698000550,0.0536810403100056,2.60949651214952) q[12];
u3(1.80646614509598,-1.71572276765529,-1.45102753112050) q[15];
cx q[15],q[12];
u1(1.99729244951968) q[12];
u3(-2.43437267260791,0.0,0.0) q[15];
cx q[12],q[15];
u3(-0.234492618806965,0.0,0.0) q[15];
cx q[15],q[12];
u3(1.95182026906508,-3.13818680223558,0.826723248952809) q[12];
u3(0.551976799351262,-3.65973245756043,-2.45266195968712) q[15];
u3(1.24034835160980,-0.699084193300202,2.16886687785744) q[2];
u3(2.22464960407988,-2.00236916890537,-1.73602404287475) q[15];
cx q[15],q[2];
u1(2.62913711251526) q[2];
u3(-1.70949596090793,0.0,0.0) q[15];
cx q[2],q[15];
u3(3.46632275374845,0.0,0.0) q[15];
cx q[15],q[2];
u3(0.358290277008894,-0.682063623050985,-0.0506199939655006) q[2];
u3(0.934361942194606,3.58083915743037,1.23835972429612) q[15];
u3(1.62527806024709,0.481199133941953,1.39247870668178) q[4];
u3(1.33995954640972,-2.71215620244124,-0.707836319427790) q[6];
cx q[6],q[4];
u1(0.946857688252911) q[4];
u3(-1.68508670664255,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.493130400398965,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.20867489259751,1.89050288609969,-3.86841677888963) q[4];
u3(2.90824905459417,-0.0408834660317494,3.00060425283812) q[6];
u3(1.64800217825654,-2.02503304985958,0.616889298058610) q[11];
u3(0.847130546208462,-4.47046804936633,0.653860025442921) q[13];
cx q[13],q[11];
u1(1.47622894982159) q[11];
u3(-0.00696604258602451,0.0,0.0) q[13];
cx q[11],q[13];
u3(1.84567374038958,0.0,0.0) q[13];
cx q[13],q[11];
u3(2.24161603372245,1.27185878552922,-1.84386697415898) q[11];
u3(1.04815031912413,1.69782901249005,2.73802613149665) q[13];
u3(0.949950863128520,2.30602791238177,-2.89912851553816) q[8];
u3(1.98724959429292,-2.87594186984013,3.25789028127545) q[10];
cx q[10],q[8];
u1(2.96495780055821) q[8];
u3(-2.46602638748133,0.0,0.0) q[10];
cx q[8],q[10];
u3(1.11996224530464,0.0,0.0) q[10];
cx q[10],q[8];
u3(2.30575769323728,-2.76616524362284,2.25827887054633) q[8];
u3(1.95123544423288,3.39488614200744,2.02215771937279) q[10];
u3(1.56730052480427,0.581950126080219,2.08286893936387) q[9];
u3(1.79938443600027,-2.24782184112402,-1.79156620180431) q[0];
cx q[0],q[9];
u1(4.44384438211587) q[9];
u3(-3.90525354907707,0.0,0.0) q[0];
cx q[9],q[0];
u3(-0.913682328051777,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.63024764506642,1.76830469190041,-1.41202491793165) q[9];
u3(0.691533962586634,2.11799219019108,-0.341767624101743) q[0];
u3(1.85689942924260,0.259830107278133,1.54599677339860) q[12];
u3(1.45804250687805,-1.04763150975642,-1.14695284719742) q[1];
cx q[1],q[12];
u1(2.52743324321808) q[12];
u3(0.0754298632425743,0.0,0.0) q[1];
cx q[12],q[1];
u3(1.29731762849004,0.0,0.0) q[1];
cx q[1],q[12];
u3(0.738711981592318,-2.60270188774094,-0.571596685949341) q[12];
u3(0.821536911334629,2.29794836538821,-0.749695063948244) q[1];
u3(1.74071488918067,0.863459680597493,-1.40757155386038) q[3];
u3(1.30216281847174,0.950246898880296,-4.82100759750148) q[5];
cx q[5],q[3];
u1(2.14770967420390) q[3];
u3(-1.76841850612334,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.79079591324863,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.81810219591944,4.08767972319268,-1.97110133451942) q[3];
u3(0.759679293547979,-1.75027693676041,4.02084823476793) q[5];
u3(1.39037272675166,-0.942467689051604,2.71656876670939) q[14];
u3(1.35886539946370,-1.94445574756196,-1.54828706358871) q[7];
cx q[7],q[14];
u1(1.38958879319097) q[14];
u3(-3.57364285027667,0.0,0.0) q[7];
cx q[14],q[7];
u3(2.27088340446012,0.0,0.0) q[7];
cx q[7],q[14];
u3(1.71079688682097,-1.93102164313247,0.990343996111640) q[14];
u3(2.66672865646088,-1.04940661863704,-2.61397670103969) q[7];
u3(1.91985122517015,-1.87208143651770,0.106389381000368) q[2];
u3(1.97512102327709,-4.05914824381283,0.394376216689718) q[9];
cx q[9],q[2];
u1(2.01589806999751) q[2];
u3(-3.20854784548130,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.57833911338677,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.87914632449648,-1.81146076133307,4.09620906682083) q[2];
u3(2.41712992730270,-1.61537549782810,0.807961853812655) q[9];
u3(0.854909537948167,3.39719625597519,-1.23287588386485) q[1];
u3(0.989040623834299,1.75727313267026,-1.91042280366662) q[8];
cx q[8],q[1];
u1(2.16718029729480) q[1];
u3(-2.41203016875702,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.43921743500858,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.80190073007429,2.03762340556101,-2.10744674376407) q[1];
u3(1.38429568662269,-0.603330097350423,3.18352114383653) q[8];
u3(0.248354295419635,-3.00139505221069,3.03957632441651) q[6];
u3(1.54241551967751,0.818095760574440,-1.55860398354211) q[5];
cx q[5],q[6];
u1(3.07795463339928) q[6];
u3(-1.94679403807767,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.29268667141944,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.794974844627086,0.357085925698839,-1.66501161862130) q[6];
u3(1.25475809383168,-6.01841202039425,0.0364640104182965) q[5];
u3(0.541935156180667,0.0181523845849827,-0.315516531729415) q[14];
u3(0.749339989684846,-2.94953649577804,0.971313488097479) q[12];
cx q[12],q[14];
u1(0.470172959163933) q[14];
u3(-1.13739684729606,0.0,0.0) q[12];
cx q[14],q[12];
u3(3.05602911849354,0.0,0.0) q[12];
cx q[12],q[14];
u3(2.06560878201687,-1.31575159123377,2.61001080886370) q[14];
u3(1.81585338654388,-2.10933393273920,1.77019982077704) q[12];
u3(1.72240877754788,-0.0749315489947331,-1.57209137975950) q[13];
u3(1.30583818925073,0.365211152749046,-2.88300726570452) q[10];
cx q[10],q[13];
u1(2.99783298907931) q[13];
u3(-1.62986464715203,0.0,0.0) q[10];
cx q[13],q[10];
u3(0.626352874332765,0.0,0.0) q[10];
cx q[10],q[13];
u3(1.13057179917029,2.14122605835901,-1.66028941991803) q[13];
u3(2.73221563479751,-1.10649483470046,0.414469669000913) q[10];
u3(2.40755631355354,0.497904335607283,-3.45094487906663) q[4];
u3(2.12382257734422,3.31879163303106,-2.79832515073629) q[7];
cx q[7],q[4];
u1(2.08978936518338) q[4];
u3(-2.98677100462300,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.25126852663125,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.42579355951809,-4.03899576384076,1.30582129413687) q[4];
u3(2.19832447806629,0.764911575385512,-3.53610408024099) q[7];
u3(1.71755279609116,-1.40111493580486,-1.46397772616887) q[3];
u3(1.07235618635517,1.43840363853207,-4.56354596098512) q[11];
cx q[11],q[3];
u1(1.89055577780066) q[3];
u3(-2.77011278599483,0.0,0.0) q[11];
cx q[3],q[11];
u3(0.850726834905266,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.06173168541254,-2.68373563081790,0.358455721239535) q[3];
u3(1.74848430910366,1.72073122927280,0.431215310838234) q[11];
u3(0.990144959521869,1.03639236430275,-1.00664942873775) q[0];
u3(0.550071471044835,-4.01345214775866,1.26824910807370) q[15];
cx q[15],q[0];
u1(2.29354556324152) q[0];
u3(-1.87610541884507,0.0,0.0) q[15];
cx q[0],q[15];
u3(2.99496701586963,0.0,0.0) q[15];
cx q[15],q[0];
u3(1.30715349555312,1.01015482784684,-3.29783606649122) q[0];
u3(1.33991398672467,1.64182198379631,3.03098853586526) q[15];
u3(1.99507159845635,0.186950951867960,0.558949177681827) q[9];
u3(1.59333098594311,-1.99644479571946,-0.992192148441799) q[12];
cx q[12],q[9];
u1(1.69959363924762) q[9];
u3(-0.644418202415845,0.0,0.0) q[12];
cx q[9],q[12];
u3(2.68453170366931,0.0,0.0) q[12];
cx q[12],q[9];
u3(1.66037893261022,0.943599592845225,-0.518248552233253) q[9];
u3(1.91675815332954,-1.62062412827590,0.613431088267841) q[12];
u3(1.83746714730590,-0.808944505781808,-1.22689513776861) q[5];
u3(1.54042091815404,-4.99120091181002,0.732673558797515) q[0];
cx q[0],q[5];
u1(2.70550151867079) q[5];
u3(-3.15382861706021,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.615515182238067,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.20662767413294,-1.38193434869890,2.59518247990934) q[5];
u3(2.46261933842021,2.97668166758614,1.61918952173442) q[0];
u3(1.48973977601039,1.70100142180155,-0.445490476298045) q[2];
u3(1.00656615387397,0.369733821194709,-4.27647797845918) q[8];
cx q[8],q[2];
u1(-0.174921115144747) q[2];
u3(-2.61762043511105,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.23063200272132,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.48742892833747,0.168576322352003,0.877004903183888) q[2];
u3(2.14444017335396,3.57379600846795,-2.39453737567506) q[8];
u3(0.405150600693487,-1.48057399731632,2.10135723375596) q[11];
u3(0.498520158868009,-2.46445129817880,0.574085752391453) q[3];
cx q[3],q[11];
u1(4.23803459696992) q[11];
u3(-3.04172342432760,0.0,0.0) q[3];
cx q[11],q[3];
u3(0.481067053215505,0.0,0.0) q[3];
cx q[3],q[11];
u3(0.948645755606814,3.31172941952940,0.132929791032318) q[11];
u3(1.32035486638873,3.00832814210938,0.142621801214041) q[3];
u3(2.38223919636442,-0.771306942909517,2.56303142667484) q[13];
u3(2.62193615115712,-3.48938507069519,-0.866480850696402) q[6];
cx q[6],q[13];
u1(3.52016563797813) q[13];
u3(-1.17435813490731,0.0,0.0) q[6];
cx q[13],q[6];
u3(2.28425135482167,0.0,0.0) q[6];
cx q[6],q[13];
u3(1.74022906446951,1.45744752962466,1.26692702789080) q[13];
u3(1.35813793912940,3.81691808904531,-1.35456994139631) q[6];
u3(2.23800115857789,1.88513986036967,-0.490227490711788) q[7];
u3(2.14342943976231,0.521548752948878,-3.64343976059407) q[15];
cx q[15],q[7];
u1(2.45906906405943) q[7];
u3(-1.86789727225822,0.0,0.0) q[15];
cx q[7],q[15];
u3(3.23717129514558,0.0,0.0) q[15];
cx q[15],q[7];
u3(2.44062534533060,1.20102169469797,-2.67820034397349) q[7];
u3(0.898527250313457,1.15072137643470,-3.51651935751615) q[15];
u3(0.871517641690863,2.01187922430350,-1.49305011311065) q[10];
u3(0.511718582841238,-0.909936823617740,0.150603504204526) q[4];
cx q[4],q[10];
u1(1.11681334535856) q[10];
u3(-1.37951377548512,0.0,0.0) q[4];
cx q[10],q[4];
u3(-0.0265730879675359,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.44848702166208,-2.14735145513926,2.09413724502161) q[10];
u3(0.957486671023931,-1.20805540790977,-2.06431684860876) q[4];
u3(2.21211919403581,1.52559838412399,0.439819799127567) q[14];
u3(1.94135462160606,1.10380242650637,-3.38268328024862) q[1];
cx q[1],q[14];
u1(-0.953330441329479) q[14];
u3(0.194287096160198,0.0,0.0) q[1];
cx q[14],q[1];
u3(3.29457570444874,0.0,0.0) q[1];
cx q[1],q[14];
u3(0.721237472288068,0.863697428727508,1.16177525375139) q[14];
u3(1.43202055454806,4.39450480190656,-0.562584950067484) q[1];
u3(1.73971085657772,-1.10768674531783,-0.768971637504340) q[2];
u3(1.62091594368235,-2.82818208661141,0.434082680513266) q[9];
cx q[9],q[2];
u1(0.731936447800166) q[2];
u3(-1.22587820340223,0.0,0.0) q[9];
cx q[2],q[9];
u3(-0.267035253680119,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.79300153708135,4.03786876815250,-0.940551195648449) q[2];
u3(1.11488676822660,4.19355937963212,-1.40005775685434) q[9];
u3(2.34127811610722,-1.61259056741166,1.37358532128534) q[0];
u3(2.29065859582485,1.34663687781023,3.54683262950828) q[15];
cx q[15],q[0];
u1(0.542862529593519) q[0];
u3(-1.47415514192434,0.0,0.0) q[15];
cx q[0],q[15];
u3(2.03976269369497,0.0,0.0) q[15];
cx q[15],q[0];
u3(2.17991123037431,-4.13027756272864,2.14305715340677) q[0];
u3(1.36748168079859,2.25006345451515,0.887346079412445) q[15];
u3(2.57533700478657,-1.77077459261729,-1.09537658055192) q[3];
u3(1.10648385747870,-1.18981089176464,-3.71754972135511) q[10];
cx q[10],q[3];
u1(0.595702125309287) q[3];
u3(-1.17797537840164,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.62093313572800,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.87296870617304,2.93369742917956,-1.88147938203371) q[3];
u3(1.90512906125805,0.0811356518205155,-4.60524662449557) q[10];
u3(1.46512355039899,0.862652505132374,-1.51319771162636) q[12];
u3(2.47370495426330,1.30848980398156,-4.54235491287556) q[14];
cx q[14],q[12];
u1(3.33510389619601) q[12];
u3(-1.41186639251629,0.0,0.0) q[14];
cx q[12],q[14];
u3(2.42700447088481,0.0,0.0) q[14];
cx q[14],q[12];
u3(1.00044092949902,-1.73339054426731,-0.329678712217289) q[12];
u3(2.22341583188928,-0.104820498308539,-3.63635561229139) q[14];
u3(1.42313020997332,0.842787974512704,2.21370721172087) q[4];
u3(1.33782841496599,-1.39177346906017,-1.75056914275361) q[11];
cx q[11],q[4];
u1(2.44465518862336) q[4];
u3(-3.11487392835812,0.0,0.0) q[11];
cx q[4],q[11];
u3(0.974583547562529,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.44104118280892,0.755170250671622,-0.246492927314799) q[4];
u3(1.14889580673971,5.01317971998450,0.828026918602148) q[11];
u3(1.87202004583105,1.35717215927315,-1.38028966690479) q[1];
u3(2.22160960518266,-4.70456558326350,1.00948704836454) q[6];
cx q[6],q[1];
u1(-0.455073119824073) q[1];
u3(0.873612555921493,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.25392230950873,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.12024119882507,-1.33369225459767,-0.918288523440375) q[1];
u3(2.22436604923920,-2.34007457610213,-0.553669621732199) q[6];
u3(0.794749954303316,-0.0704695197688390,-1.18844808608704) q[7];
u3(0.276773998531590,-2.82206522944056,1.07544087352482) q[5];
cx q[5],q[7];
u1(1.93797282212957) q[7];
u3(-0.138783499780117,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.677218012161557,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.87367947427839,-1.04502753610892,2.90881347885158) q[7];
u3(1.32189546011455,-3.63078998009587,0.860710952753884) q[5];
u3(2.53339507710775,0.155956216143114,-1.95003184145678) q[13];
u3(1.76652747886652,-3.77445697927224,1.54215443270639) q[8];
cx q[8],q[13];
u1(3.09267900916484) q[13];
u3(-1.69574963075632,0.0,0.0) q[8];
cx q[13],q[8];
u3(0.894602712266336,0.0,0.0) q[8];
cx q[8],q[13];
u3(1.01602727602607,-0.622051541073496,1.18554734215610) q[13];
u3(1.47750031892925,-1.29380872449084,3.52835462565163) q[8];
u3(1.56341234777118,-1.01889513293969,1.32473579352499) q[14];
u3(1.98054093202899,-1.68529629872325,-2.19718026509824) q[1];
cx q[1],q[14];
u1(1.45205529994180) q[14];
u3(-2.22948017202032,0.0,0.0) q[1];
cx q[14],q[1];
u3(0.122979096606763,0.0,0.0) q[1];
cx q[1],q[14];
u3(0.814507263402315,-0.312745093878782,-2.96309366703073) q[14];
u3(1.13824050036044,4.75305164517189,-1.30165557043271) q[1];
u3(1.31757266132381,-2.40063071637617,0.327059749010578) q[12];
u3(1.51224232769881,-3.80652013010451,-0.109555046911306) q[3];
cx q[3],q[12];
u1(1.62328919324136) q[12];
u3(-2.93461199541103,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.845239668147242,0.0,0.0) q[3];
cx q[3],q[12];
u3(0.572412841026281,2.75707457597172,0.248753585803874) q[12];
u3(2.11610669197740,-3.46108996517034,-1.25578274625205) q[3];
u3(1.08902524020195,0.286860620069606,-3.22497609884177) q[15];
u3(1.16410523888238,3.17972560460984,-2.59698105017474) q[8];
cx q[8],q[15];
u1(3.28107721464961) q[15];
u3(-2.12513443678557,0.0,0.0) q[8];
cx q[15],q[8];
u3(1.88766339214138,0.0,0.0) q[8];
cx q[8],q[15];
u3(1.49977833363831,3.32061908654418,-1.49389485857228) q[15];
u3(2.04245607499050,-0.249219534860237,3.32255088064968) q[8];
u3(1.52696734386236,3.00835925001657,-2.54282014019616) q[5];
u3(0.958782391780543,2.75850563582391,-2.65419768136063) q[0];
cx q[0],q[5];
u1(1.48045865919943) q[5];
u3(-0.416062968416978,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.20584973574838,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.77982786189234,-1.32052764153914,-0.457445829642109) q[5];
u3(2.03031849319941,-4.00722838176563,1.70345119237597) q[0];
u3(2.65365560915125,-0.0954616645583442,3.18331610397676) q[13];
u3(1.98506593680334,1.43653611998938,1.88560722831240) q[6];
cx q[6],q[13];
u1(1.01330370170442) q[13];
u3(-3.35993626811299,0.0,0.0) q[6];
cx q[13],q[6];
u3(1.83251229387034,0.0,0.0) q[6];
cx q[6],q[13];
u3(1.12688958354212,1.82291240973420,-2.06725962612151) q[13];
u3(1.57923210525183,0.166617536365972,1.14600361211405) q[6];
u3(1.65774928313819,0.214009544635638,1.75235511088432) q[10];
u3(1.76309569898219,-1.17737842970964,-1.58608144018821) q[9];
cx q[9],q[10];
u1(1.35846103469208) q[10];
u3(-0.664378604311361,0.0,0.0) q[9];
cx q[10],q[9];
u3(2.69163588045420,0.0,0.0) q[9];
cx q[9],q[10];
u3(2.23557413436063,3.60302358981719,-2.33952658632014) q[10];
u3(2.21252656261531,2.08224711199047,-1.71086513506281) q[9];
u3(1.08618236472109,2.71165480609945,-1.64107366527599) q[7];
u3(1.11002489399608,1.74005958671423,-2.56634040868464) q[11];
cx q[11],q[7];
u1(0.0417986483057693) q[7];
u3(-0.713325284810025,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.94273963958177,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.21931208479049,0.892154895186096,-1.12724306109591) q[7];
u3(1.19974089914888,-3.98210998783735,2.22460318199459) q[11];
u3(2.69801230107287,-0.258664848708381,3.06731130578699) q[2];
u3(2.22020729854736,-2.35634340679327,-0.145690401111096) q[4];
cx q[4],q[2];
u1(0.581320259584936) q[2];
u3(-0.327485120945330,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.98299424162960,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.58724766527237,2.99368017232058,-1.23692014113275) q[2];
u3(0.963993434337092,4.83703181985606,-0.746160804962149) q[4];
u3(1.64990352991112,2.80998094218215,-1.75713172890317) q[10];
u3(0.918628282216071,3.15579070294072,-2.48858217445198) q[15];
cx q[15],q[10];
u1(1.10027912435255) q[10];
u3(-3.34189537898486,0.0,0.0) q[15];
cx q[10],q[15];
u3(2.44258621647469,0.0,0.0) q[15];
cx q[15],q[10];
u3(2.77157858037053,3.35550375074916,-1.33387552908457) q[10];
u3(1.76739946327814,-2.33911199461620,-0.305222341493528) q[15];
u3(1.60243903111898,-1.24499838865243,1.26967357369663) q[0];
u3(1.06211445854203,-2.71050023857254,-0.234908659426091) q[5];
cx q[5],q[0];
u1(0.486320524366641) q[0];
u3(-1.47809684784906,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.0599449201469193,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.478735418958019,0.119188695389689,0.900139960194477) q[0];
u3(1.60802173517508,1.17211444759496,-0.464908368275425) q[5];
u3(0.128572477579163,3.08783051314981,-3.11181468680866) q[6];
u3(1.01982673330659,-4.11100801228054,1.35940537418987) q[1];
cx q[1],q[6];
u1(-0.475075009680071) q[6];
u3(-1.59641867755023,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.837364466858510,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.30312073982231,-0.0485159708555493,-1.45424436305019) q[6];
u3(1.53156945518425,-1.64817795880946,-0.541649708008729) q[1];
u3(1.27405778442037,0.215453740112197,1.35193791453589) q[11];
u3(1.19309106324201,-1.47389072042646,-2.84542536997016) q[14];
cx q[14],q[11];
u1(0.344917056261834) q[11];
u3(0.0204967871556523,0.0,0.0) q[14];
cx q[11],q[14];
u3(1.27516134054307,0.0,0.0) q[14];
cx q[14],q[11];
u3(2.07855340203145,-0.698225352052163,1.37418125667857) q[11];
u3(0.978162312258829,-0.728726850938939,-0.483251586763276) q[14];
u3(2.78798773942125,2.40958591882594,-1.28071234321318) q[9];
u3(2.39461032312995,0.0563595504087133,-5.73590350735413) q[12];
cx q[12],q[9];
u1(2.27702481469721) q[9];
u3(-2.60881015543646,0.0,0.0) q[12];
cx q[9],q[12];
u3(1.39657799176911,0.0,0.0) q[12];
cx q[12],q[9];
u3(2.09935786116312,-0.842414035031597,0.729526787361955) q[9];
u3(0.952277968988429,3.29175200107228,2.20532323387914) q[12];
u3(2.54137004644965,-0.273125717197200,-0.565419963673792) q[3];
u3(1.47216303438538,1.03301609142318,-5.20477143450356) q[2];
cx q[2],q[3];
u1(2.91135165333396) q[3];
u3(-1.76151421263863,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.08468014144293,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.26155390973196,1.38593304475098,1.62956208722628) q[3];
u3(0.934880581122844,-2.07893638871629,-0.668945517209811) q[2];
u3(0.669038018537443,0.794217642555063,-0.653532993631931) q[4];
u3(0.698473193100149,-1.09012584441314,-1.45367635553633) q[7];
cx q[7],q[4];
u1(0.763808778080820) q[4];
u3(-1.24441794545391,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.128089437884472,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.62294318132539,1.05993895429682,-1.47690299715737) q[4];
u3(1.59717788715766,0.997652020367702,2.55129403999882) q[7];
u3(1.60090335102872,-0.683720780296624,2.91767859838729) q[8];
u3(1.87254034836701,-1.86578089234723,-1.60693546127753) q[13];
cx q[13],q[8];
u1(1.23538382177110) q[8];
u3(-1.07765104447340,0.0,0.0) q[13];
cx q[8],q[13];
u3(1.75012990412649,0.0,0.0) q[13];
cx q[13],q[8];
u3(2.39193311705734,-2.91918457156153,1.07996156086054) q[8];
u3(1.53541838435975,1.85521056883642,-2.37808488379658) q[13];
u3(1.44903487762204,2.43263145715651,-1.37325789008854) q[14];
u3(1.70164193355842,1.21390961524983,-0.481356793094886) q[5];
cx q[5],q[14];
u1(1.63059712335902) q[14];
u3(-0.437076199839973,0.0,0.0) q[5];
cx q[14],q[5];
u3(-0.347069513896622,0.0,0.0) q[5];
cx q[5],q[14];
u3(2.32503937236107,-4.11004108518742,1.65782202966312) q[14];
u3(1.43131067553993,0.0920765918111566,0.0583015265469681) q[5];
u3(2.75480434546599,1.68636914270830,-0.214311374345756) q[3];
u3(2.93955079501223,5.16745319099634,0.650725972329613) q[8];
cx q[8],q[3];
u1(0.911101519557236) q[3];
u3(-0.0699022018548181,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.83699668729221,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.33430837394519,-1.78409604118787,0.780411204531108) q[3];
u3(2.00336365228741,2.98630746384709,1.63647936395963) q[8];
u3(1.26873016911256,1.75828555084917,-3.53983037863467) q[2];
u3(1.47941638438508,3.46666827983988,-2.63874398319775) q[13];
cx q[13],q[2];
u1(1.47812171994018) q[2];
u3(-2.17903528115416,0.0,0.0) q[13];
cx q[2],q[13];
u3(0.493013972877436,0.0,0.0) q[13];
cx q[13],q[2];
u3(0.598021656271863,-0.117269260749840,-1.87327742230957) q[2];
u3(1.89555114777140,3.28296774699638,2.10844526033807) q[13];
u3(1.90211524927071,1.52367794543601,0.998776916578635) q[11];
u3(0.840551748492657,0.186318874764657,-3.28108676940923) q[15];
cx q[15],q[11];
u1(2.82838378291394) q[11];
u3(-2.01382260704364,0.0,0.0) q[15];
cx q[11],q[15];
u3(1.07083168021438,0.0,0.0) q[15];
cx q[15],q[11];
u3(2.39813305698165,4.39286001435639,-1.35483268829062) q[11];
u3(0.930822188497547,-5.28576090534283,0.403085120950951) q[15];
u3(0.940267544208382,-0.387613099601322,-1.27099599988177) q[0];
u3(1.53694805655208,-2.72962378393563,-0.195584413362883) q[4];
cx q[4],q[0];
u1(-0.354694565397212) q[0];
u3(-1.69725394005745,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.596604005669344,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.99594205596910,-0.525744813794162,0.845358807358764) q[0];
u3(1.22161334418329,1.32750237899300,4.90082437875456) q[4];
u3(1.85047048717112,0.346193452198802,1.38679957826838) q[9];
u3(1.33025560566847,-2.46731691276777,-1.89103197626914) q[7];
cx q[7],q[9];
u1(1.70617119560931) q[9];
u3(0.00575150377638001,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.663578544740868,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.12434438311027,2.07987968330279,-3.39296116817852) q[9];
u3(1.09187265848866,-1.80539016408638,-4.17018420204217) q[7];
u3(2.49879558010230,1.71247127305461,-4.43463964713803) q[1];
u3(1.02217641159260,1.75395263745083,-1.00959960108474) q[10];
cx q[10],q[1];
u1(1.52669746808846) q[1];
u3(0.642835375275064,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.30121781101135,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.384507731425149,-3.32782315943292,1.89933278370764) q[1];
u3(1.80002592184643,-2.92638308891279,1.71783019666327) q[10];
u3(2.20797818797901,1.47374131567012,-4.29830490092662) q[12];
u3(1.10378093889411,-1.59557774859487,2.17227357238492) q[6];
cx q[6],q[12];
u1(4.54175670310544) q[12];
u3(-3.89113834253281,0.0,0.0) q[6];
cx q[12],q[6];
u3(-0.700556738026226,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.39192743850677,-2.50724978694811,1.14996784383657) q[12];
u3(2.69995906884580,-0.289981857732934,4.76179618517812) q[6];
u3(1.37263292282867,1.35914485975476,0.960426917916819) q[13];
u3(0.875146428315904,0.816237348486959,-3.16436451399053) q[11];
cx q[11],q[13];
u1(1.99207510813910) q[13];
u3(-2.95594710595379,0.0,0.0) q[11];
cx q[13],q[11];
u3(0.451496344534292,0.0,0.0) q[11];
cx q[11],q[13];
u3(0.924441032514831,-0.0407191448824961,1.75267626826418) q[13];
u3(1.34149757305825,0.876604719436982,-3.22734879464806) q[11];
u3(2.24823152419452,1.16661400779098,-0.447020801373937) q[9];
u3(2.35815297060231,0.330124243163303,-4.10537819700643) q[0];
cx q[0],q[9];
u1(1.80423955487996) q[9];
u3(-2.68426089037928,0.0,0.0) q[0];
cx q[9],q[0];
u3(3.08547342180309,0.0,0.0) q[0];
cx q[0],q[9];
u3(0.719372966968284,-3.33185732941660,2.24799770078283) q[9];
u3(2.31035286381982,3.45713501207976,2.53087755867789) q[0];
u3(0.668011347401427,-0.992581938278153,1.14774183614747) q[8];
u3(0.886076938821363,-0.532525031563344,-2.18896348247083) q[14];
cx q[14],q[8];
u1(1.50621396386244) q[8];
u3(0.153878540438628,0.0,0.0) q[14];
cx q[8],q[14];
u3(2.01270734818063,0.0,0.0) q[14];
cx q[14],q[8];
u3(1.43724955101211,-0.174790132602539,4.50628678512714) q[8];
u3(1.19089725714933,-1.76346373247382,1.59169167770011) q[14];
u3(1.68257798180440,-0.0919197404994631,1.30151610715006) q[12];
u3(2.24805633115287,-0.839407278870145,-2.25768959556081) q[15];
cx q[15],q[12];
u1(1.99771303048932) q[12];
u3(-1.70818198969452,0.0,0.0) q[15];
cx q[12],q[15];
u3(-0.0650952385954668,0.0,0.0) q[15];
cx q[15],q[12];
u3(1.62629172345187,-0.633556678499330,-1.10429856627862) q[12];
u3(1.40390499963238,2.07561883590082,-1.41297023616687) q[15];
u3(0.782721531970511,0.309562892621418,0.139021228403860) q[7];
u3(0.693924042292426,-2.17319762364899,1.00868943852401) q[2];
cx q[2],q[7];
u1(3.65187121341721) q[7];
u3(-4.33725456093773,0.0,0.0) q[2];
cx q[7],q[2];
u3(-0.226048435294277,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.73007962211562,1.27777346938971,0.743742765773048) q[7];
u3(1.18174665853797,0.414208952999140,-5.21271842992383) q[2];
u3(1.42487813972492,0.379504324524371,-1.34169639880383) q[10];
u3(2.31992968486714,-3.75680743983198,1.99930293317990) q[6];
cx q[6],q[10];
u1(2.28569188816263) q[10];
u3(-2.40802505178497,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.03743339700219,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.610237865503345,-0.794112464527142,-1.15060208510121) q[10];
u3(2.56057776273362,0.248621532960948,5.28804143642789) q[6];
u3(1.86809394743153,-0.0358990343533525,-2.13268428115937) q[5];
u3(1.71163234416884,-3.26955928838757,1.54177978995173) q[4];
cx q[4],q[5];
u1(2.46989773864555) q[5];
u3(-1.57196498046306,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.0622170362595731,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.27877054903996,-1.18887459461594,2.81924949461784) q[5];
u3(1.11617532747765,-0.733550248807935,1.00074843213685) q[4];
u3(1.02931048029032,0.427155431753325,0.824909631345226) q[1];
u3(1.77529207469871,-1.03476697040421,-1.20272248393117) q[3];
cx q[3],q[1];
u1(0.611437932160010) q[1];
u3(-1.64719357650836,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.72635323967820,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.01904215945144,1.47340313665441,0.635041780416553) q[1];
u3(1.44404918956223,2.01431109519716,0.000921264076186867) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
