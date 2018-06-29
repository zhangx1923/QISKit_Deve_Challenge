OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(2.04957380876631,-0.383494552819898,-0.0453819329185728) q[4];
u3(1.26015950687530,-3.29366181531635,-1.00342015016806) q[2];
cx q[2],q[4];
u1(1.66928004224845) q[4];
u3(0.118570767755972,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.830600545470845,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.40750067340361,-1.65456642387646,1.47609393833069) q[4];
u3(2.54957875483307,-0.308241809950806,1.09106713552005) q[2];
u3(2.44040245927774,0.589646627518280,-0.957825127705105) q[1];
u3(0.903901448579629,0.305359808558326,-4.43288138314817) q[5];
cx q[5],q[1];
u1(2.64954081666145) q[1];
u3(-2.02575145308849,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.46413204409025,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.907492752365815,3.12611820992845,0.745898544703383) q[1];
u3(0.735934735489032,-2.83114479657415,3.35646464457429) q[5];
u3(2.77560856287843,-1.35361592436578,-1.33564246587913) q[3];
u3(0.692249146613222,-1.13828567557340,-3.95952353622431) q[9];
cx q[9],q[3];
u1(2.36510140010281) q[3];
u3(-1.77047556661152,0.0,0.0) q[9];
cx q[3],q[9];
u3(0.0863497278748735,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.259130038081487,-2.13961387065826,1.30235158486774) q[3];
u3(2.08314328230663,1.20231425059295,-3.63680320358485) q[9];
u3(1.70199111515843,0.467603268143756,-2.21591757337177) q[7];
u3(2.00690790015950,1.81746720154559,-3.84148509400833) q[6];
cx q[6],q[7];
u1(2.86730866461466) q[7];
u3(-1.62467620836013,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.664717991088467,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.16798424965337,-1.42499769103330,-2.13991614807931) q[7];
u3(1.26943475339552,2.71387756619205,2.51967552511212) q[6];
u3(1.44062641563110,-2.28839299885134,0.439954936440349) q[11];
u3(1.99578843196999,-3.22202380422499,0.283266410091743) q[8];
cx q[8],q[11];
u1(0.778084007565278) q[11];
u3(-3.33236444112665,0.0,0.0) q[8];
cx q[11],q[8];
u3(1.97946620335593,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.98679532805537,-0.316813218779370,1.50784234354979) q[11];
u3(2.45685967313299,-0.135386412357275,-3.27337815641387) q[8];
u3(2.46971915730789,-1.44939799956459,1.49137018478003) q[0];
u3(1.69673183391392,1.57153970848029,3.89133120699792) q[12];
cx q[12],q[0];
u1(1.15900117218654) q[0];
u3(-3.35796560597448,0.0,0.0) q[12];
cx q[0],q[12];
u3(1.66353980868958,0.0,0.0) q[12];
cx q[12],q[0];
u3(1.39610849305898,3.07554898966169,-3.15348972417868) q[0];
u3(1.98231902519940,4.13075660886350,0.287112141889911) q[12];
u3(1.84479918236137,-0.931742397641178,-1.36006149563953) q[7];
u3(1.37146939464161,-5.07816201831122,0.670442377677787) q[5];
cx q[5],q[7];
u1(3.48670806590853) q[7];
u3(-1.59913226722158,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.05188897796656,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.12755806386393,2.26115350838035,-0.419982389891959) q[7];
u3(1.67509069728330,-0.539893558835812,-4.63125587218544) q[5];
u3(0.238253802070984,-2.13847453243737,2.67260727858805) q[8];
u3(0.826187281973558,1.54691084913102,-2.80748395608104) q[0];
cx q[0],q[8];
u1(2.45899294032369) q[8];
u3(-3.02964652402115,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.559666167993812,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.71928687461161,-1.26683535722938,2.20061268774719) q[8];
u3(1.49232444073546,-3.74725918123365,-1.95399975887312) q[0];
u3(1.18582772400167,-1.36858696141854,-0.263379405017060) q[6];
u3(2.10277034455004,-2.29419054872852,-0.614914456669465) q[4];
cx q[4],q[6];
u1(0.388688387089755) q[6];
u3(-1.47488256734780,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.24432988384397,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.452191321906555,-2.51707054492223,0.270779239856188) q[6];
u3(1.23644578508506,-0.877640802235283,2.82828687927229) q[4];
u3(2.75643895265860,-0.831498926306990,2.95864073932384) q[2];
u3(2.56400676654913,-0.137206290757163,1.66595396322911) q[1];
cx q[1],q[2];
u1(2.41057182056419) q[2];
u3(-1.50895221636477,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.960532590745167,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.87937635629407,-0.918814868932344,-1.12024572334452) q[2];
u3(1.49002710956795,3.01454546971492,0.762591579937824) q[1];
u3(0.762773814118893,3.37487800901604,-0.650251603608680) q[3];
u3(1.63085230355580,0.593444436264441,-0.695163644476956) q[12];
cx q[12],q[3];
u1(3.77725954259076) q[3];
u3(-3.22938662969658,0.0,0.0) q[12];
cx q[3],q[12];
u3(-0.921372359110347,0.0,0.0) q[12];
cx q[12],q[3];
u3(0.780529942104045,1.41413215368705,0.513727550586222) q[3];
u3(2.59483526662089,0.878193623189627,3.37621327077375) q[12];
u3(1.31673271242951,3.31646057872370,-0.325797620840511) q[10];
u3(2.96848685057463,0.953736876446501,-2.04571526854699) q[11];
cx q[11],q[10];
u1(-0.0147210519288161) q[10];
u3(-0.975733094042209,0.0,0.0) q[11];
cx q[10],q[11];
u3(2.50720370491982,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.44435592644955,3.71524571859231,-1.87326671367510) q[10];
u3(0.220681898045555,-1.56628206387113,-2.10422888083850) q[11];
u3(1.20495158298582,-0.115239456994971,1.92457254788199) q[0];
u3(1.51333217419763,-0.0823953326475555,0.0314011241765251) q[11];
cx q[11],q[0];
u1(-0.347411660575376) q[0];
u3(-1.64151673573032,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.835917897842464,0.0,0.0) q[11];
cx q[11],q[0];
u3(2.00402814201574,3.40366468021135,-0.480782591973442) q[0];
u3(2.15579144643112,1.03383077033008,4.83881756024853) q[11];
u3(0.991106503594148,-3.47778819759130,1.10022302819477) q[10];
u3(1.96886129144426,-5.71618305829869,-0.315888076164133) q[1];
cx q[1],q[10];
u1(0.494940659215640) q[10];
u3(-1.35735705989631,0.0,0.0) q[1];
cx q[10],q[1];
u3(3.00972488722168,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.73209715742972,-1.72888517239274,-0.651496585355979) q[10];
u3(2.63279094212153,-0.905699934248786,2.34172988609306) q[1];
u3(0.210230495693400,-2.70497251917015,3.24332910865410) q[12];
u3(0.419290938120165,0.666783373709621,-1.31285084843350) q[3];
cx q[3],q[12];
u1(1.51544066352788) q[12];
u3(-0.130475856181152,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.640736831151986,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.09506615180810,0.244330566163245,-0.569817903178299) q[12];
u3(2.10939858884412,-1.63113072121189,1.58971392873750) q[3];
u3(1.33747964957795,3.27250877452047,-1.45490496455934) q[4];
u3(1.52250124065794,0.311164807099475,-3.05622302183957) q[7];
cx q[7],q[4];
u1(2.95314591409508) q[4];
u3(-1.39450614498401,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.83870471586835,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.709316247335449,-2.79760901467753,2.05399098900119) q[4];
u3(2.24061964818310,0.631160783450801,2.75230856215707) q[7];
u3(1.41854292937110,2.66497441894601,-1.18985003015507) q[6];
u3(2.60062691856844,1.98285716574093,-0.579788645412400) q[2];
cx q[2],q[6];
u1(-0.848953604521449) q[6];
u3(-1.58745294871885,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.41745099404181,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.17093866215054,1.13950035406746,1.47114478345215) q[6];
u3(1.16408477818527,0.306319132541819,2.42119282301247) q[2];
u3(0.289824628535339,0.0504392398820800,0.303783665856724) q[8];
u3(0.603693634578605,-0.365589272029447,-2.14968086913066) q[9];
cx q[9],q[8];
u1(0.774685286303296) q[8];
u3(-3.46122394056943,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.84810013592816,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.98948032329327,0.863913253057887,1.06507170273292) q[8];
u3(2.26246299323801,-4.24088703508868,1.86441393765919) q[9];
u3(0.718229939193197,3.41412436894692,-1.58161261697961) q[3];
u3(1.39541421247619,1.46985711893878,-2.95600980579474) q[6];
cx q[6],q[3];
u1(0.680306436570856) q[3];
u3(-1.53237683322638,0.0,0.0) q[6];
cx q[3],q[6];
u3(-0.364509972454908,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.453186034138549,-1.12641886019227,0.602231539100260) q[3];
u3(1.08612135187908,2.48717056700255,3.38903768179705) q[6];
u3(2.12805857040867,-0.825037587047663,-0.842184315376402) q[12];
u3(0.462943816354292,-0.715607758533729,-3.77609573398663) q[7];
cx q[7],q[12];
u1(1.43984795096960) q[12];
u3(-1.10688381419376,0.0,0.0) q[7];
cx q[12],q[7];
u3(-0.503869540962008,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.549763835042781,1.65923163232964,-1.70820286537449) q[12];
u3(1.28422088724166,0.494537524446673,-2.50901159688958) q[7];
u3(1.28776931176401,-0.232639957351889,2.70078323619201) q[10];
u3(1.48407084459353,-2.31703844920108,-1.66674320773693) q[2];
cx q[2],q[10];
u1(0.367348486076578) q[10];
u3(-0.562689014415477,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.42646824891181,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.51168973209664,-1.79805102743181,-1.38421607695978) q[10];
u3(1.86745888681350,-0.923727346404918,-4.96048055411586) q[2];
u3(2.69390047854192,-1.00010520099283,2.48729714011543) q[4];
u3(1.80419597644208,-1.58945031077442,0.659188882685548) q[0];
cx q[0],q[4];
u1(-0.164839312960724) q[4];
u3(-2.40834251676961,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.10175006167628,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.45993247247947,-3.40467320255088,2.42708227091744) q[4];
u3(1.11246593126889,3.08765893706634,-1.66217265638772) q[0];
u3(0.966217518679875,0.287368213276214,-1.45585281621591) q[8];
u3(1.94027519306824,-3.64799238971357,1.41541348078902) q[5];
cx q[5],q[8];
u1(1.82549429605070) q[8];
u3(-2.47259968288442,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.0932038250875782,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.24869812425765,-1.31587705019256,-1.14315674500264) q[8];
u3(1.53658154343922,-0.599706555918075,3.21741042966178) q[5];
u3(1.10775097285782,-0.0298604954634288,0.522408048207196) q[1];
u3(1.62440039010947,-0.804981504291553,-1.98887108754870) q[9];
cx q[9],q[1];
u1(0.793622142212086) q[1];
u3(-3.18351399376938,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.76286286118020,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.18182614874992,-0.0449578657105474,2.62873750211862) q[1];
u3(0.710692100208984,-2.45827876761766,-3.52750399666707) q[9];
u3(1.56556611155310,0.829325323514748,-3.26173427265272) q[2];
u3(0.980734526081629,-3.24247625817982,2.16907593308924) q[0];
cx q[0],q[2];
u1(1.69330862760489) q[2];
u3(-2.17794484072863,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.597544500994679,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.575524079514294,0.686611099911267,3.04151698639399) q[2];
u3(2.80222104270590,3.68820157321091,2.03042918260451) q[0];
u3(1.62929552136984,0.261294691273442,-1.30589091990076) q[8];
u3(0.539786649836866,-4.15351790191268,1.40256854259375) q[3];
cx q[3],q[8];
u1(2.20523547588968) q[8];
u3(0.166186255054839,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.49733038108937,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.14598964080778,0.645787938831937,-1.59899330209696) q[8];
u3(1.55605518257763,2.13013511226405,0.0575140128102087) q[3];
u3(1.90473471407717,0.913826544758453,1.19877929781200) q[10];
u3(1.25621040889942,-1.38065782411874,-1.94551884454134) q[6];
cx q[6],q[10];
u1(2.50126781208299) q[10];
u3(-3.19046552003061,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.43472263081933,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.98361810803039,-1.96973761578731,-0.896640905345030) q[10];
u3(1.43932961996787,2.52804399578900,1.42024431461843) q[6];
u3(2.16795012913686,-1.44979745268070,0.177375274323711) q[12];
u3(2.07964847329984,-4.18304991597799,0.691134087817691) q[7];
cx q[7],q[12];
u1(3.12045006046406) q[12];
u3(-2.47866777255856,0.0,0.0) q[7];
cx q[12],q[7];
u3(0.440180146035546,0.0,0.0) q[7];
cx q[7],q[12];
u3(1.86705401015614,3.97078935678936,-1.07281871955171) q[12];
u3(1.46353833690292,-3.41175790813696,1.58041263907030) q[7];
u3(1.00892600655877,0.297914545696884,-1.16778814491862) q[1];
u3(0.428499305302258,-2.70184519749947,0.506204775223774) q[9];
cx q[9],q[1];
u1(0.399210164470226) q[1];
u3(-0.861418889281650,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.32141463219611,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.512837481085758,-1.84398699168257,1.53211306612724) q[1];
u3(1.38444126825164,1.53486173792439,-0.618570640726157) q[9];
u3(2.31145912691219,4.03600155229329,-0.960908064118805) q[11];
u3(1.23736114422933,2.58117724856993,-1.11771754002528) q[4];
cx q[4],q[11];
u1(0.437561567405064) q[11];
u3(-1.30606065080409,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.47779779581200,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.16582842892925,-1.74140426726683,-0.621919221907458) q[11];
u3(1.50428640197053,4.50640569986019,0.618182273470117) q[4];
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
