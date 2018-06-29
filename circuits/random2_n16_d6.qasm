OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(1.71165773960733,-1.76098492761430,-0.873098685010111) q[12];
u3(1.44055294288027,-4.33765801225796,0.0302291275528606) q[2];
cx q[2],q[12];
u1(2.50708957445759) q[12];
u3(-1.92036608702422,0.0,0.0) q[2];
cx q[12],q[2];
u3(-0.120555910575902,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.67250198303288,0.868683131529641,-0.911939091090624) q[12];
u3(0.751106962077019,0.781743886205900,4.90975756346927) q[2];
u3(1.01015488912935,0.0778411900898330,2.00171360394083) q[14];
u3(1.31701133865554,-2.01340748073396,-1.59957818374740) q[1];
cx q[1],q[14];
u1(2.34353432598112) q[14];
u3(-3.28165015077621,0.0,0.0) q[1];
cx q[14],q[1];
u3(1.36132426842965,0.0,0.0) q[1];
cx q[1],q[14];
u3(2.85887072859574,3.82508813554136,-1.98044133884051) q[14];
u3(1.42107091033295,-2.93121993328937,1.52275258739552) q[1];
u3(2.23058897157705,-2.43972987767802,0.0255112855431858) q[5];
u3(1.87006985972735,-3.51264876998083,-1.21185392771391) q[13];
cx q[13],q[5];
u1(-0.917365044084953) q[5];
u3(0.554275405742183,0.0,0.0) q[13];
cx q[5],q[13];
u3(3.56996145090776,0.0,0.0) q[13];
cx q[13],q[5];
u3(1.25148842988260,2.47273488764449,1.36267718615744) q[5];
u3(1.67527272492019,-5.95919525426195,-0.313096436702774) q[13];
u3(0.698398372723359,0.813497809539414,-2.27750428536668) q[0];
u3(1.65834770861463,-5.18006252432142,1.06221758934553) q[10];
cx q[10],q[0];
u1(1.00851923152960) q[0];
u3(-0.623455645064307,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.97285610158686,0.0,0.0) q[10];
cx q[10],q[0];
u3(0.931417486856478,-0.925331946499314,-0.909946211764987) q[0];
u3(1.03475048072454,0.470045281231041,-1.72471633541022) q[10];
u3(0.391302382241817,1.10066194060768,-2.82069199837925) q[9];
u3(1.67159421586000,-2.62217023293563,3.49801419849464) q[4];
cx q[4],q[9];
u1(1.62617957733923) q[9];
u3(0.105492601273430,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.595804728496240,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.625482431830682,-2.19842449055097,3.24504737385215) q[9];
u3(2.60968095563133,-2.38922100752683,3.45555396308643) q[4];
u3(1.82358193081016,2.27453183324654,-2.06101744202619) q[8];
u3(2.27283526053455,1.39355129743192,-1.05393778023110) q[3];
cx q[3],q[8];
u1(3.59905235547439) q[8];
u3(-1.11883298633832,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.68084803793174,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.55546236654460,1.96976453951878,-3.12094389050068) q[8];
u3(2.40149086141816,4.09252309182870,2.17171055105009) q[3];
u3(1.74362986005533,-2.41685794905549,0.541266460757419) q[6];
u3(2.00925024424480,-3.30682401058015,-1.04438076928016) q[11];
cx q[11],q[6];
u1(1.34956999107608) q[6];
u3(-3.54837714760930,0.0,0.0) q[11];
cx q[6],q[11];
u3(2.46529349008097,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.59325370598526,3.60947884321459,-2.29789513261818) q[6];
u3(1.79626055109704,4.91361159023232,-1.09623380301880) q[11];
u3(1.54477340071818,1.84664348044464,-2.99838188284301) q[7];
u3(2.92691166070680,-3.59724062329211,2.04670555393965) q[15];
cx q[15],q[7];
u1(3.25791947641794) q[7];
u3(-1.56235581227109,0.0,0.0) q[15];
cx q[7],q[15];
u3(1.76314313659111,0.0,0.0) q[15];
cx q[15],q[7];
u3(0.813318727352939,2.47440852939287,-0.252230435667791) q[7];
u3(1.61946298682227,-1.39817003908632,-0.384571020688142) q[15];
u3(0.934814896935248,1.45093530616367,-1.26689378324540) q[13];
u3(0.596287795266814,1.04572022389350,-3.13180974413102) q[11];
cx q[11],q[13];
u1(-0.128361000763916) q[13];
u3(-1.53700406881932,0.0,0.0) q[11];
cx q[13],q[11];
u3(0.527398806302239,0.0,0.0) q[11];
cx q[11],q[13];
u3(2.10484831502133,1.32654352301229,-3.18671198791829) q[13];
u3(0.0184163460170579,4.45046179419972,1.24849001597829) q[11];
u3(0.668305354540470,-2.24203372018756,1.97133185373944) q[8];
u3(0.492144245634707,-2.92047975254523,1.07392952759180) q[10];
cx q[10],q[8];
u1(1.70862317223686) q[8];
u3(-0.606716560858133,0.0,0.0) q[10];
cx q[8],q[10];
u3(2.99403027992517,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.61536544497646,0.720751513987140,0.652316656362938) q[8];
u3(2.27653176090471,1.63844283136386,-2.65491367482216) q[10];
u3(2.59649621995916,1.21489235356817,-0.140738594044497) q[5];
u3(1.67939005878452,-0.196609566208364,-1.82495231212784) q[2];
cx q[2],q[5];
u1(0.206504333089694) q[5];
u3(-0.945238541101376,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.76290727438314,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.568841386749096,-2.77669710072069,0.565971648159363) q[5];
u3(1.53813470759455,-0.260831856975168,5.83323732025043) q[2];
u3(0.984425595460247,-3.18985793600029,2.87098205631126) q[14];
u3(0.788797995233657,0.0231404148010645,-1.55265697019751) q[7];
cx q[7],q[14];
u1(0.814544992377447) q[14];
u3(-1.58787460912680,0.0,0.0) q[7];
cx q[14],q[7];
u3(-0.488591333825599,0.0,0.0) q[7];
cx q[7],q[14];
u3(0.891853341934676,3.74357702358754,-1.10762595487690) q[14];
u3(1.05455849722362,2.03568475888120,-4.02909879848973) q[7];
u3(1.93313469633790,3.64462522274274,-2.02347837562755) q[1];
u3(0.312814607111642,1.84665717282155,-1.18569583681478) q[9];
cx q[9],q[1];
u1(4.12943248428553) q[1];
u3(-3.78900577542794,0.0,0.0) q[9];
cx q[1],q[9];
u3(-0.401488631675764,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.23465161883575,0.957070351924521,-0.374520363580847) q[1];
u3(1.31624174172449,4.29160715736401,-1.61112515518032) q[9];
u3(2.34280034278238,0.662086758980458,0.575689347879030) q[12];
u3(0.249243664332019,-3.94320612255415,-0.701023213713744) q[15];
cx q[15],q[12];
u1(2.13651867338510) q[12];
u3(0.763983563639127,0.0,0.0) q[15];
cx q[12],q[15];
u3(1.57889343058944,0.0,0.0) q[15];
cx q[15],q[12];
u3(1.02040703673821,1.34161964909696,-3.17140148924928) q[12];
u3(1.64316836772206,0.798297412780251,-0.944921254047564) q[15];
u3(2.00620646620224,1.19495692794452,-3.21519040515859) q[0];
u3(2.31854610595753,-3.06218607151837,3.08235780575098) q[3];
cx q[3],q[0];
u1(-0.0504097718722307) q[0];
u3(-0.362283962466283,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.07841593519870,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.743876565888754,2.47667631795245,-3.53081129649034) q[0];
u3(0.381930140754037,-1.84700145209087,-4.10558467355881) q[3];
u3(1.37975052970424,-2.32128241768238,1.47944607127388) q[6];
u3(0.891427596294256,-2.21149619468944,0.330005182139657) q[4];
cx q[4],q[6];
u1(3.94925895387035) q[6];
u3(-4.25149164271094,0.0,0.0) q[4];
cx q[6],q[4];
u3(-0.568778453105901,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.05063465949714,0.949435439933522,-1.97115306696118) q[6];
u3(1.09733784286271,1.12379260879056,0.371113804582541) q[4];
u3(2.03368218722030,2.89574136145519,-0.913987871106758) q[4];
u3(2.13267078053423,1.05875059906800,-1.70553514288038) q[8];
cx q[8],q[4];
u1(1.26561932561932) q[4];
u3(-3.34989066928240,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.00425221261289,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.73804459752110,-1.85131480596878,1.49949804421908) q[4];
u3(0.807123286299969,-0.970288961121153,3.05111678306604) q[8];
u3(0.633364989024004,-1.39438385671266,1.01779851665149) q[2];
u3(0.988660125916163,-0.543136784691618,-1.42515993514621) q[7];
cx q[7],q[2];
u1(0.0692263427641613) q[2];
u3(-2.40568684262305,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.18074693553764,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.886930481659938,-1.01023961790634,5.02539229960123) q[2];
u3(0.329228128985608,-3.92015990761470,0.117059154201875) q[7];
u3(1.23684417060444,1.52743600692608,-2.09494550450848) q[11];
u3(0.355745806813701,-2.05862227634274,1.32190491782017) q[3];
cx q[3],q[11];
u1(1.56528942139784) q[11];
u3(-2.39051814074092,0.0,0.0) q[3];
cx q[11],q[3];
u3(0.116210707117035,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.51434790454016,1.19175799413404,-0.891006749970789) q[11];
u3(0.339152081267766,1.45483710741321,2.16276263633036) q[3];
u3(0.432005120733507,2.13092647261236,-0.0120038211339493) q[10];
u3(1.44333511329805,0.415195741585876,-3.66289086337499) q[15];
cx q[15],q[10];
u1(2.97753995865074) q[10];
u3(-2.32621335743450,0.0,0.0) q[15];
cx q[10],q[15];
u3(1.20566364454840,0.0,0.0) q[15];
cx q[15],q[10];
u3(1.72412631565880,2.17386430215656,-0.903261797859973) q[10];
u3(2.67877311936068,5.88485076116999,-0.344880041712390) q[15];
u3(2.87062630261656,-0.515731884689886,2.95636504756934) q[9];
u3(1.85794183355804,-2.52707919733758,-1.98595770957235) q[6];
cx q[6],q[9];
u1(0.575338341852180) q[9];
u3(-1.50896360524710,0.0,0.0) q[6];
cx q[9],q[6];
u3(2.81167777509536,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.64743144993042,-4.34668039309327,1.17427543004302) q[9];
u3(1.91991178922545,0.230927297975141,-4.18856371706603) q[6];
u3(1.55484508002700,3.03080709626470,-1.39148435361108) q[1];
u3(2.71240739127433,1.90315663515983,-0.838789321092003) q[0];
cx q[0],q[1];
u1(2.88881886874043) q[1];
u3(-1.65408356829768,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.985780023873016,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.651438338126822,3.69319769587001,-2.29199869068030) q[1];
u3(0.599665528887635,0.0868263186380549,5.76062676700792) q[0];
u3(1.18073638578915,1.62278003973674,-3.36694970710868) q[12];
u3(0.609913582687328,2.02945175000484,-2.73254975613929) q[13];
cx q[13],q[12];
u1(2.12456170644470) q[12];
u3(-3.16228190351168,0.0,0.0) q[13];
cx q[12],q[13];
u3(1.65908721309701,0.0,0.0) q[13];
cx q[13],q[12];
u3(2.49347323969087,3.13858109854517,-0.849267243020412) q[12];
u3(2.48576107758225,0.0490407265048309,0.469286923156500) q[13];
u3(1.95220241699705,1.55678009309815,-0.858955139088241) q[5];
u3(2.37077418393097,-0.000282974328486452,-3.72700334272996) q[14];
cx q[14],q[5];
u1(2.61503340777789) q[5];
u3(-1.77224671618394,0.0,0.0) q[14];
cx q[5],q[14];
u3(-0.00714966523280136,0.0,0.0) q[14];
cx q[14],q[5];
u3(1.67788997372732,0.855374128983601,-0.955149492088639) q[5];
u3(1.49362864117720,-0.478029086785255,1.03603836973122) q[14];
u3(1.32739049934211,0.787164359786960,-2.82333080822217) q[2];
u3(2.08181804720621,3.33512712771816,-2.83302215938715) q[8];
cx q[8],q[2];
u1(1.71259196413826) q[2];
u3(-0.660141098653301,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.72177958332521,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.17073306842765,-0.226702049556107,-0.167675864428414) q[2];
u3(2.31663490568635,-2.90439251154897,1.93142350697278) q[8];
u3(1.31488433240496,-1.97348465006307,-0.117728607568665) q[3];
u3(1.22120010960012,-4.32477197757370,1.04070104232716) q[14];
cx q[14],q[3];
u1(1.54513603736587) q[3];
u3(0.237530691976484,0.0,0.0) q[14];
cx q[3],q[14];
u3(1.09222427502063,0.0,0.0) q[14];
cx q[14],q[3];
u3(0.584835018449243,1.61623433372183,-1.62508752009961) q[3];
u3(1.40714514406975,1.03597189181026,0.838600809048471) q[14];
u3(0.935090074765386,1.78329830382973,-0.0497247548442858) q[10];
u3(1.24656901582607,-0.0608649054846013,-3.96594074301508) q[4];
cx q[4],q[10];
u1(2.71285772641687) q[10];
u3(-2.56890541299601,0.0,0.0) q[4];
cx q[10],q[4];
u3(-1.41802014013798,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.623290663810191,0.702130309114221,-0.943953008187256) q[10];
u3(2.02313480963645,-1.97978232097065,-2.10736884157865) q[4];
u3(2.40802357675430,-1.18631586743926,0.0531981234525911) q[12];
u3(1.73112159683739,-3.70331881031465,-0.973645993358145) q[15];
cx q[15],q[12];
u1(2.40001533746820) q[12];
u3(-1.81033425761379,0.0,0.0) q[15];
cx q[12],q[15];
u3(0.171148155542915,0.0,0.0) q[15];
cx q[15],q[12];
u3(2.20057945289448,-4.14306107933439,0.780016604655700) q[12];
u3(0.837573156914172,-0.393953449295612,1.93291146817680) q[15];
u3(1.86066329801245,1.99173949061922,0.0906136824099083) q[5];
u3(2.05423084078498,-0.758438415855840,-4.21531224504434) q[7];
cx q[7],q[5];
u1(0.307398383413449) q[5];
u3(-0.844221895823297,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.21535075969649,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.401380901301694,-1.76213501887526,0.665899199195736) q[5];
u3(1.09989363879219,-1.70189466198731,1.59879005510664) q[7];
u3(1.14050686222208,3.31593201400844,-1.50439368026489) q[13];
u3(1.30237068539336,0.857595230839150,-2.30723596665781) q[1];
cx q[1],q[13];
u1(1.72282118836746) q[13];
u3(-2.61569383165723,0.0,0.0) q[1];
cx q[13],q[1];
u3(3.35072154850910,0.0,0.0) q[1];
cx q[1],q[13];
u3(2.27335468033754,2.41556469605505,-1.00584144531917) q[13];
u3(1.50577505008308,2.11777896369970,-3.40323130461972) q[1];
u3(2.34419013200261,0.171865236684418,-1.26935845586057) q[9];
u3(1.56583112454038,0.0402050466286499,-3.86622119575010) q[11];
cx q[11],q[9];
u1(0.156919714096749) q[9];
u3(-1.43342495029944,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.45560069647217,0.0,0.0) q[11];
cx q[11],q[9];
u3(2.66623564603480,4.59462446724560,-1.40465261829020) q[9];
u3(2.02099106816976,1.56578115018937,-0.951999408634145) q[11];
u3(1.87953716804752,-1.58126381288956,1.55662386538307) q[0];
u3(2.47632053645942,-3.47704366035457,-0.104870949615392) q[6];
cx q[6],q[0];
u1(1.62490631096103) q[0];
u3(-3.01622829277234,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.21705140323124,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.08531224499167,-2.49301696247009,1.81268922348975) q[0];
u3(0.811994743124363,1.77410999156071,0.869115131743204) q[6];
u3(1.56483612525758,1.98351497963561,-3.14290725990766) q[10];
u3(2.76562870087795,-3.10747646538496,2.60213374856226) q[0];
cx q[0],q[10];
u1(2.79484163058593) q[10];
u3(-1.99401163863311,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.307502416008423,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.54533187927781,1.46822752192290,-1.44576439560669) q[10];
u3(2.96603312497000,3.92527068894130,-1.07686153483139) q[0];
u3(1.52867399457966,-1.55490531073481,0.832268046823373) q[5];
u3(1.07956397364793,-2.52854903824021,-0.0435554140442589) q[6];
cx q[6],q[5];
u1(0.864526390494404) q[5];
u3(-0.183118510698893,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.95332368412975,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.64227899114391,0.00726012792225206,-1.79036137247212) q[5];
u3(2.76669802450715,1.15206592276772,4.02838388434394) q[6];
u3(0.629575710036214,2.03213386191417,-2.86609894195779) q[1];
u3(0.819986516249650,1.46454938104620,-2.93723497206762) q[3];
cx q[3],q[1];
u1(0.540302451089527) q[1];
u3(-1.63470773292394,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.21308697015946,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.36840722414186,0.0570553391856994,2.67857488533056) q[1];
u3(1.34643443171962,0.689454935795636,4.87938587016957) q[3];
u3(3.07329898919083,0.0382709411793842,-0.0918362928458844) q[9];
u3(1.66504557024544,0.377733304883602,-5.82600148880911) q[12];
cx q[12],q[9];
u1(2.18354943902024) q[9];
u3(-1.97905003638968,0.0,0.0) q[12];
cx q[9],q[12];
u3(0.719694803116433,0.0,0.0) q[12];
cx q[12],q[9];
u3(0.938298702372144,2.91235621041915,0.775923880905636) q[9];
u3(1.44894353098054,1.51383785641539,2.62765009407628) q[12];
u3(0.985081067898289,-0.888787736285855,1.54966025939263) q[15];
u3(0.757204840306849,-1.93017192461550,0.583451379603211) q[13];
cx q[13],q[15];
u1(1.85225128451682) q[15];
u3(-1.60750138279458,0.0,0.0) q[13];
cx q[15],q[13];
u3(0.950958557467046,0.0,0.0) q[13];
cx q[13],q[15];
u3(2.11411206204884,1.00307551582203,-1.43786506272288) q[15];
u3(2.00178610415378,-0.535839667039897,-1.77315877017502) q[13];
u3(1.13678187854029,0.261985868790722,-2.04072698024986) q[4];
u3(1.66422113049208,-4.53067260182423,1.69268862182872) q[2];
cx q[2],q[4];
u1(2.60865194699286) q[4];
u3(-1.69126290404796,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.325111109965041,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.36812842493948,-1.15934667483931,-1.25066153737099) q[4];
u3(0.386464387949054,-2.96536074206051,0.669564637622321) q[2];
u3(1.26986594149927,1.82985135927579,-1.41650592680692) q[8];
u3(0.811211294931582,-1.30497456066302,0.423335672478480) q[11];
cx q[11],q[8];
u1(1.08762299308497) q[8];
u3(-3.49408971966094,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.88760615360303,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.50487701296712,3.07425814647136,-0.826516618916251) q[8];
u3(1.79883991180047,0.789087734805411,-2.14608904282383) q[11];
u3(1.76156791797720,-0.431260777467001,3.28536130227524) q[7];
u3(1.59791435193150,-0.447842349466862,0.464691323924498) q[14];
cx q[14],q[7];
u1(2.94809054422312) q[7];
u3(-2.05321928227887,0.0,0.0) q[14];
cx q[7],q[14];
u3(1.14731754197989,0.0,0.0) q[14];
cx q[14],q[7];
u3(1.83664035745832,-0.540093546829165,-0.596796649846125) q[7];
u3(2.23664864111305,-0.606239747800699,-0.208415748859883) q[14];
u3(0.855982975025065,-0.0909339505261151,0.495882528242186) q[0];
u3(0.732374131766115,-0.317299181780566,-2.68070050998658) q[13];
cx q[13],q[0];
u1(1.54947686665271) q[0];
u3(-1.02630610688018,0.0,0.0) q[13];
cx q[0],q[13];
u3(2.65531043747756,0.0,0.0) q[13];
cx q[13],q[0];
u3(1.73220073002439,-0.441430205452975,1.13007753442342) q[0];
u3(1.26956608383691,4.46698156660223,-1.78731543327409) q[13];
u3(2.30655983609754,-1.72220897410842,-1.15593233549513) q[11];
u3(0.912766914387706,-4.13813524609066,-0.481897763193275) q[5];
cx q[5],q[11];
u1(0.452603395515628) q[11];
u3(-0.834324409241612,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.41079110929869,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.75866360094440,1.34443233708911,0.574272147865046) q[11];
u3(1.82304929093009,0.805124647376536,-2.72984547521524) q[5];
u3(2.37542266878702,0.970435870498070,-1.45016710130681) q[1];
u3(1.81468717066640,2.07797615154085,-3.79422836313373) q[4];
cx q[4],q[1];
u1(1.83962824959289) q[1];
u3(-3.18879568745648,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.587343008569542,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.99653548961710,-0.563925665622659,0.789516874961841) q[1];
u3(0.748912463205430,-4.16423267116790,0.305173076330149) q[4];
u3(1.65877555232813,0.201473683109865,-2.73561122558692) q[12];
u3(2.71046970582817,-0.723486569567429,-4.75478960710259) q[9];
cx q[9],q[12];
u1(3.56275073072653) q[12];
u3(-1.42543310798692,0.0,0.0) q[9];
cx q[12],q[9];
u3(2.42060059818763,0.0,0.0) q[9];
cx q[9],q[12];
u3(2.88049645263889,1.41908552013360,1.13518375829370) q[12];
u3(0.937474465030243,0.547514689700048,1.09405892600529) q[9];
u3(0.301834214000451,2.33201559086489,-1.76115198316358) q[10];
u3(0.910214883870793,0.509105497355690,-2.26910024279639) q[6];
cx q[6],q[10];
u1(-0.467899390782340) q[10];
u3(-1.97357835914186,0.0,0.0) q[6];
cx q[10],q[6];
u3(0.958709437104131,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.95759927332715,1.66846320979369,-4.48892810382744) q[10];
u3(2.60210763876995,0.226486599272818,-4.78377424615782) q[6];
u3(1.56615066983898,-1.08180091609690,0.190460371284451) q[8];
u3(1.65300491959740,-2.65009781240945,0.572156145122441) q[2];
cx q[2],q[8];
u1(2.60197293855178) q[8];
u3(-1.90043372485864,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.109378843034847,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.17296498998892,-4.03774536809942,1.67907245255561) q[8];
u3(2.47186103661238,-0.743620792503176,-4.26114165430317) q[2];
u3(2.39520618635488,2.70733138093750,-1.24948946302891) q[7];
u3(0.965921854374790,1.03194939609776,-1.67089302093797) q[14];
cx q[14],q[7];
u1(1.32713377994306) q[7];
u3(-0.949958165690825,0.0,0.0) q[14];
cx q[7],q[14];
u3(-0.164720558597927,0.0,0.0) q[14];
cx q[14],q[7];
u3(0.845931325609719,1.30483219780692,-3.41614642948141) q[7];
u3(1.41885636054023,0.831259972796579,-4.24332518049675) q[14];
u3(0.476619104811329,1.99230945116247,-3.34039218211351) q[15];
u3(1.20160383925083,-3.00657505172056,3.19745900626286) q[3];
cx q[3],q[15];
u1(2.07372841987141) q[15];
u3(-1.47114443310297,0.0,0.0) q[3];
cx q[15],q[3];
u3(3.70432102157605,0.0,0.0) q[3];
cx q[3],q[15];
u3(1.43933991363839,-2.69837076560574,0.00131325742080746) q[15];
u3(2.60569320010044,1.45812956288494,-2.40522019205231) q[3];
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
