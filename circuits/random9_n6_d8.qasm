OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.08025275341697,0.00586977101698705,-1.92151026074510) q[3];
u3(1.82039626165179,0.498476980731326,-5.02569203609374) q[4];
cx q[4],q[3];
u1(0.674578853685951) q[3];
u3(-1.24109753192540,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.77476361148824,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.54137034176027,-3.85510280660852,1.83074502873391) q[3];
u3(2.25630767730610,1.54878479157834,2.11650060970317) q[4];
u3(2.03038744499777,3.10624997163469,-0.266901009029060) q[2];
u3(2.18976560866331,2.68615157110175,-1.00296748863038) q[0];
cx q[0],q[2];
u1(2.77364508572650) q[2];
u3(-1.94762467063486,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.0242724085118016,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.85310430237557,-0.566749386191315,-1.45308576008578) q[2];
u3(0.879961368870838,-1.93969886226302,4.03784369962209) q[0];
u3(0.462828732556249,1.76751530815159,0.00573264173470722) q[5];
u3(1.32806012597680,-0.158494108519270,-2.82890013082876) q[1];
cx q[1],q[5];
u1(1.43089777057821) q[5];
u3(-0.164887948255329,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.41386508495940,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.98697530062580,-0.402669538289075,-0.826302409644555) q[5];
u3(2.14347639946180,-2.40549619397061,-2.15244247586498) q[1];
u3(2.59644651349661,-3.60211789801910,1.06951839629457) q[4];
u3(1.90925671267280,-0.356478723659877,3.38674225271110) q[3];
cx q[3],q[4];
u1(1.52278830352298) q[4];
u3(-0.277647884987345,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.64502524058547,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.51874998155477,0.0728269176259954,0.725262793599288) q[4];
u3(2.23630581197853,-2.50260225314044,-0.565870126150836) q[3];
u3(2.25707885213992,-4.11993037188040,1.97113904100626) q[5];
u3(0.567181482320160,-1.21187990476489,1.81439216651394) q[1];
cx q[1],q[5];
u1(-1.08501445126686) q[5];
u3(0.657970630848131,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.84260069012856,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.88724238280003,-2.39331186115069,0.585288536130842) q[5];
u3(0.764539546552742,0.196859362640755,3.25118865370895) q[1];
u3(0.692476352063557,-2.81926399298480,3.03299951675296) q[0];
u3(1.09893267237850,1.00905762454497,-1.82940633119820) q[2];
cx q[2],q[0];
u1(3.27706525314736) q[0];
u3(-0.475321444446070,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.83767935557401,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.39703938667901,3.98875667005137,-1.12338458927342) q[0];
u3(1.08743934636443,-3.13247828790967,1.44868489577128) q[2];
u3(1.32572068068881,-1.18288264275868,0.594895536016163) q[4];
u3(1.78845622040476,-1.39341281288165,-1.61416447911679) q[3];
cx q[3],q[4];
u1(1.74037690020828) q[4];
u3(0.389043026139029,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.00666298894866,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.85400600653878,3.05508040282902,-2.58131876553345) q[4];
u3(2.38076417685414,3.51898799696816,-1.00239572372176) q[3];
u3(1.28281248965935,1.12021325955868,-3.94516632701329) q[0];
u3(0.686581323670133,-2.73878086528599,2.79247862072437) q[2];
cx q[2],q[0];
u1(0.0110745423826346) q[0];
u3(-0.870434769189441,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.66724935333678,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.69456242309211,1.73382798482009,-3.46319994794275) q[0];
u3(1.81076732897974,-0.526770043445009,-4.59430192587475) q[2];
u3(2.14904093675644,-0.278182956306906,1.67318498713191) q[1];
u3(1.89360048108565,-2.34351279118716,-1.55256168158438) q[5];
cx q[5],q[1];
u1(2.84556258224815) q[1];
u3(-1.29075729748644,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.414075088908519,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.10371443431669,4.31362942798966,-0.0446625940089311) q[1];
u3(0.203388019020850,0.524900600312329,-1.30866996013245) q[5];
u3(0.810558620494502,1.32550727835433,-3.13408361989445) q[4];
u3(1.95071689587439,-3.17977595923427,3.06870136579453) q[0];
cx q[0],q[4];
u1(2.09991943946711) q[4];
u3(0.142045129591919,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.41583133104203,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.99779344151113,0.0348258701985817,2.13267423511440) q[4];
u3(0.272775651407694,-0.240849896495845,-1.86675792533230) q[0];
u3(2.08211234347717,-1.59574610004802,-0.302310578787464) q[3];
u3(1.95140037765654,-3.61235467314217,0.611161087352716) q[1];
cx q[1],q[3];
u1(3.31542392318177) q[3];
u3(-1.60833843092051,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.82928383831959,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.45578740125474,-1.81129632779111,1.83617308364459) q[3];
u3(1.02273532127115,-1.75699564479424,-0.589386349800761) q[1];
u3(0.755246720849443,2.33849655258171,-3.76950921776263) q[5];
u3(0.888725667411685,0.0661791134609815,-1.73919856885911) q[2];
cx q[2],q[5];
u1(0.783065188592880) q[5];
u3(-0.223984331469930,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.73197698988107,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.379379665539977,1.66914777967935,-2.88175405518550) q[5];
u3(1.46715994898034,4.66503622390847,1.16486334408695) q[2];
u3(2.85627724734264,-3.42083857010660,2.81745428838995) q[4];
u3(1.43415743664973,3.22455886674066,-2.21564781375807) q[2];
cx q[2],q[4];
u1(2.33967600901655) q[4];
u3(-1.66288796530661,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.0329702224116784,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.802539989632470,3.04246050046981,-0.526020558465693) q[4];
u3(2.07500549291722,1.98134523994035,-4.24677891491483) q[2];
u3(1.02089153912647,0.335548389682576,0.450162635071349) q[0];
u3(0.880124102121061,-1.70727341468207,-1.78844488110414) q[1];
cx q[1],q[0];
u1(2.38889057117692) q[0];
u3(-1.81236815147823,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.128354904241248,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.30463534726101,-1.38855461263255,3.20500558472379) q[0];
u3(1.39334568959418,-3.37035610648822,-1.15631483474499) q[1];
u3(0.780151977597951,0.435916116555071,-3.51943771052670) q[5];
u3(1.70836931843246,-1.25778498961212,4.79262226239491) q[3];
cx q[3],q[5];
u1(2.01665351088660) q[5];
u3(-2.75711816212228,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.765282308497869,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.38375634510258,2.28048804558454,0.0182657211254134) q[5];
u3(2.07915255322815,-5.30354155985128,0.851783589735908) q[3];
u3(2.01354755830796,0.0894358465587291,-0.814272584504998) q[2];
u3(1.11330292899205,0.771544834286734,-5.07565052012883) q[5];
cx q[5],q[2];
u1(1.34854462671180) q[2];
u3(-0.711231913219643,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.87531158353032,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.24802970496977,-1.99647169371805,-0.0134702911662390) q[2];
u3(2.03238220994419,4.76595209819589,-0.776060197037829) q[5];
u3(1.17507016506390,-0.00206676892031132,0.711807039623919) q[1];
u3(1.24085329692917,-2.96180683401470,-0.642598176370215) q[4];
cx q[4],q[1];
u1(0.317594802744989) q[1];
u3(-1.20671955394915,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.34478920576746,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.57880450503404,-0.889347790586518,1.60623534372028) q[1];
u3(1.87875744974426,-0.0149642963349294,1.44697073285294) q[4];
u3(2.33037789665835,3.11560967707881,-2.33986149297771) q[0];
u3(0.520764784272671,-0.717014326987286,1.79710825573456) q[3];
cx q[3],q[0];
u1(2.64959136521857) q[0];
u3(-2.18414878567240,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.0889397877337768,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.34581813790034,2.24270233415948,-2.27331128260529) q[0];
u3(1.59162584521508,3.31527600805795,1.08750049674538) q[3];
u3(0.464613364707881,2.68152302114828,-2.91180015913215) q[4];
u3(0.805265318519210,-0.252993941374768,-1.41063822360999) q[0];
cx q[0],q[4];
u1(-0.0836589504366267) q[4];
u3(-0.955536631494806,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.33855276925373,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.13834826460184,-0.934824068736300,0.294773060509793) q[4];
u3(1.17586586455119,-1.50423767340486,-3.11478039665572) q[0];
u3(2.48618764901005,2.04374038917923,-3.12135076311852) q[2];
u3(1.41999548195232,-2.85226570772939,2.74890941665643) q[3];
cx q[3],q[2];
u1(1.45638710383905) q[2];
u3(-3.12496956074643,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.36069327684379,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.66808451013222,0.114208549867214,-1.36147101617636) q[2];
u3(1.49787207422974,5.15389709856274,-0.0240903558411678) q[3];
u3(1.80241265016088,-0.742001125618165,1.37715857831566) q[1];
u3(2.08543987937113,-2.32346412611947,-2.29717737626166) q[5];
cx q[5],q[1];
u1(2.81811305015731) q[1];
u3(-1.80519836340678,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.976482517253910,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.25512446281306,2.53854032963105,-0.119758388019770) q[1];
u3(2.33740018033323,2.97078191274225,0.986859929989850) q[5];
u3(1.36521333980356,1.75756529829327,-0.219814274638226) q[2];
u3(2.51466940830182,1.24678824589101,-1.86524792193580) q[3];
cx q[3],q[2];
u1(1.44879715333923) q[2];
u3(-2.34600007617998,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.556927176682328,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.58048912469360,2.91394002848202,-2.36004601001867) q[2];
u3(0.948282641521031,4.20034918673704,1.35984732279056) q[3];
u3(1.02653204251542,-1.85493347586886,-0.908453818063755) q[1];
u3(0.797756290126278,-2.98967287686649,-0.140149618691277) q[4];
cx q[4],q[1];
u1(1.82701441807692) q[1];
u3(-2.96268888515689,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.655087325723833,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.24606257786910,-1.87953066816782,-1.09106446161193) q[1];
u3(1.88163728176291,-1.93457714509838,-1.07990943693642) q[4];
u3(1.42609185069535,1.05232183101350,1.02497653399700) q[0];
u3(1.45676846558405,0.0667156052679072,-3.65536152809401) q[5];
cx q[5],q[0];
u1(-0.196227692329884) q[0];
u3(0.655414611796050,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.97144636157032,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.46997252427408,0.909504959964929,-4.53491878083735) q[0];
u3(2.13189734661962,5.11524766398630,-0.0854102901351435) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];