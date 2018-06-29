OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.84556441492058,1.63700938728994,-3.16481816675873) q[0];
u3(1.38443412168883,-2.78469753316823,3.47608004269369) q[2];
cx q[2],q[0];
u1(0.959626306259747) q[0];
u3(0.118469889526408,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.97092695173462,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.67834064202125,-1.54360966709004,0.406581496822358) q[0];
u3(1.87191742371879,-6.21871991906398,0.0532477597196830) q[2];
u3(2.12222058190352,2.82672884000559,-0.760721379653626) q[4];
u3(2.23211465026708,0.654370651116999,-2.78706161103371) q[1];
cx q[1],q[4];
u1(1.87559513323753) q[4];
u3(-2.36883119908353,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.271276143265155,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.36292338537373,-3.65490033934974,1.47938218608106) q[4];
u3(1.36442774829700,0.713392543185946,4.31272549900296) q[1];
u3(0.966719554180849,-0.322667167190336,1.01304471766495) q[1];
u3(0.949832170659581,-2.06187506067214,-1.41770732490910) q[4];
cx q[4],q[1];
u1(2.25312193127797) q[1];
u3(-1.77445267995388,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.03948428411763,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.80563580461637,3.78511010340128,-2.49318660331628) q[1];
u3(2.37146357491020,-3.42472507238006,-2.67709708754444) q[4];
u3(2.30175593753425,-0.849065534638047,2.89805824383288) q[0];
u3(2.98295994861387,1.25248375152512,2.97074118326410) q[2];
cx q[2],q[0];
u1(1.60470577854242) q[0];
u3(-3.14160789985320,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.980305927344998,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.79417528091360,-2.82569462016672,1.28524253967997) q[0];
u3(2.10666272620739,1.21186488630816,-1.57979747476338) q[2];
u3(2.06057055855904,0.909823591161872,-3.21718181039696) q[0];
u3(0.873014231007908,-2.07785192633191,2.49568077778083) q[2];
cx q[2],q[0];
u1(-1.17234905133675) q[0];
u3(0.462463810177311,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.35416947557945,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.90041825527824,-3.67940268191451,0.255919275026179) q[0];
u3(2.39833075330390,-2.48345505842274,-0.896270395106136) q[2];
u3(2.14363405212684,0.0681942238207823,2.02787965992634) q[4];
u3(2.08702455769300,-1.33954716020857,-0.867149297031387) q[1];
cx q[1],q[4];
u1(3.33333613366681) q[4];
u3(-1.40320562153029,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.48162479083440,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.89758905823215,-1.86652216061973,3.84293705378261) q[4];
u3(2.33512991611397,1.65437981603350,-2.56194117570369) q[1];
u3(1.47684625355116,1.03328312544271,-2.14852916695833) q[4];
u3(1.23955201887716,1.93963037379545,-3.75499687390682) q[3];
cx q[3],q[4];
u1(1.53994109036564) q[4];
u3(-0.190543553090072,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.51980413964125,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.90467087935029,0.121209437999958,-0.422480333056470) q[4];
u3(1.87429935932064,-0.0603652974964790,-1.76036131496915) q[3];
u3(2.45312981008941,0.710966453783197,0.583872815399123) q[1];
u3(0.621862469339809,-3.26632571592623,-0.937332990980246) q[0];
cx q[0],q[1];
u1(2.69462775034919) q[1];
u3(-1.73376660116288,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.31310942997695,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.62347167393426,-1.30128851935219,3.02619340380149) q[1];
u3(0.845769408353727,1.68286706119176,4.58407589476423) q[0];
u3(0.575992201376461,-0.872028456782385,1.26078836913498) q[0];
u3(0.932560627919609,-2.72852644553998,1.32049486798124) q[2];
cx q[2],q[0];
u1(1.62499410984212) q[0];
u3(-3.05772976442875,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.937492771528218,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.97483668486397,-1.80098530402997,2.92829234614294) q[0];
u3(1.31751972119564,2.26833599453940,2.72073740083523) q[2];
u3(1.96313857844923,-0.860081210889307,1.73713722019443) q[1];
u3(1.96598457903900,-2.06616056433066,-2.19585582414968) q[4];
cx q[4],q[1];
u1(0.200602148881234) q[1];
u3(-1.31750940835159,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.26359812761958,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.31992277795589,-3.58379573826861,0.426629255716028) q[1];
u3(2.03047516013642,-1.98500834236546,0.199537959564042) q[4];
u3(2.06935875855216,3.20504785706115,-2.99953516655499) q[2];
u3(0.845889062028759,0.457600621567576,1.38733896803272) q[0];
cx q[0],q[2];
u1(1.79320181661818) q[2];
u3(-2.93361419871790,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.616378323366134,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.39931627655319,-0.640697489783171,0.992799325056100) q[2];
u3(0.840246293870303,-0.567541240525698,-1.64054416518543) q[0];
u3(1.60996574157596,0.851267739574406,2.24053662748740) q[4];
u3(2.65365504432426,-1.08900132062710,-1.74425455574405) q[1];
cx q[1],q[4];
u1(0.789566304529619) q[4];
u3(-1.42924053157490,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.88660174854850,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.11623885864406,3.18257526993170,-0.483774451658621) q[4];
u3(1.86440467157442,0.446275707536566,0.0784485969570411) q[1];
u3(0.858185119587551,2.39639379680864,-0.165947947624667) q[1];
u3(0.987647650439488,0.303924619788686,-4.32143550301774) q[0];
cx q[0],q[1];
u1(1.07942246910629) q[1];
u3(-3.21702957727205,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.79583824508330,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.86861336015379,3.39008616014099,-2.71700049804712) q[1];
u3(2.11108289518032,3.38021596700241,-1.27633015433733) q[0];
u3(1.32861708752019,0.0777523674468781,-2.26444340595101) q[4];
u3(1.75457939425685,-0.0824355428994461,-5.07389960118462) q[2];
cx q[2],q[4];
u1(0.905326428583624) q[4];
u3(-3.12487863490085,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.92876371445786,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.74526431830453,1.34855678587411,-0.467575453285526) q[4];
u3(0.985286079736709,1.52852410862134,-3.57988545425295) q[2];
u3(1.62662793913842,0.334329050662220,2.03545010864434) q[2];
u3(2.44335042744639,-0.896551273579989,-2.09073448665975) q[4];
cx q[4],q[2];
u1(0.626742240912131) q[2];
u3(1.35955184483672,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.43699166059023,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.70364797878151,2.08845016196160,-0.918038570338290) q[2];
u3(1.65236547625808,-4.56257628667595,-0.490594438198098) q[4];
u3(1.67004152748783,1.09566965785948,-1.54540044713682) q[3];
u3(1.89369771870254,-5.38235054697111,0.363000666836496) q[1];
cx q[1],q[3];
u1(3.53355125927549) q[3];
u3(-1.29362572559751,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.43348331452884,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.658631545464190,0.665345793653939,-4.26255706340301) q[3];
u3(1.43605389367355,-0.901647392212131,4.64732381097723) q[1];
u3(1.63859242648680,1.98622489276038,-0.187094394815538) q[0];
u3(1.16668832150422,0.153889490728405,-3.63940906385331) q[2];
cx q[2],q[0];
u1(1.54516870832279) q[0];
u3(-0.312583009656484,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.44276371644788,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.16803549590546,3.42230087747501,-0.652134526400848) q[0];
u3(0.626660017437147,1.77405907320880,2.11826745478573) q[2];
u3(0.522843025527049,1.06113589093801,-1.76095839370224) q[3];
u3(0.647506448681120,-3.80797722966022,1.33985770877829) q[1];
cx q[1],q[3];
u1(2.91507649631804) q[3];
u3(-1.29179150659172,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.446863851954998,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.84413041890608,-0.359057138063654,1.63177591589676) q[3];
u3(1.29236990961409,2.26579907240114,3.41153931534827) q[1];
u3(0.987862068017903,0.503826863876788,0.768248232318101) q[3];
u3(0.988057553227239,-0.576320932367154,-1.75736694322535) q[1];
cx q[1],q[3];
u1(3.01301131566870) q[3];
u3(-1.84311722772559,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.675254120565339,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.19227921263466,1.23286161174832,1.82292400408364) q[3];
u3(0.457191785640253,3.44258150526629,2.82934437334048) q[1];
u3(2.57262503267881,1.49012132368719,-2.29185498759835) q[2];
u3(2.49964911841786,0.919583550243755,-4.84519341536866) q[4];
cx q[4],q[2];
u1(0.629264917174604) q[2];
u3(-1.43594397301794,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.157748446972640,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.43504948191942,0.373941863182528,-2.10464768402591) q[2];
u3(1.03733408050653,-1.66985003891343,1.53387432567754) q[4];
u3(1.39480238692512,3.35916303212509,-1.52862164964049) q[4];
u3(2.52884420047429,2.02315276696259,-1.45550419723939) q[2];
cx q[2],q[4];
u1(2.47983676430859) q[4];
u3(-2.04686163148423,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.16442928188494,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.52941201072004,1.59504149642146,-0.326077884085038) q[4];
u3(0.684368088802784,0.231217044210065,5.14444312621675) q[2];
u3(0.924265355105443,0.377097396411472,1.04953917553967) q[1];
u3(0.965621284089501,-1.01873111062765,-3.13211334660602) q[3];
cx q[3],q[1];
u1(3.88455471560049) q[1];
u3(-4.07130849119366,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.347177496354730,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.66935595309963,0.722228215027979,3.26796109587086) q[1];
u3(1.03171599556325,0.141335368942108,2.41068165570457) q[3];
u3(2.88151395896456,1.74586920100510,-2.67734480871805) q[4];
u3(1.93502441855561,-2.50540801107580,2.07333612266037) q[1];
cx q[1],q[4];
u1(2.21776342270438) q[4];
u3(-1.78515028354405,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.228252854229013,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.05939626946822,0.0127060238253968,-1.99883207403317) q[4];
u3(2.78014846951821,3.33043224899535,-2.44256028086598) q[1];
u3(1.20391123801257,-1.02060741291265,2.23609345290232) q[3];
u3(1.23179984974108,-1.61813265268957,-1.98183386736316) q[0];
cx q[0],q[3];
u1(3.10994970334159) q[3];
u3(-1.47192783320670,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.39250020369645,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.47252374504855,-1.12432390273051,0.145091591529498) q[3];
u3(0.814728998451002,2.45476122120899,0.452267042647184) q[0];
u3(1.81237070613576,-1.38381895377723,-0.405775196463231) q[4];
u3(1.28144644295521,-2.27162779056738,-0.107889953649313) q[3];
cx q[3],q[4];
u1(0.311172421456973) q[4];
u3(1.33177415495294,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.02768859496911,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.44381810319219,-0.590439405625515,0.190631591473714) q[4];
u3(0.950114160078458,1.70031045290042,1.16178544695514) q[3];
u3(2.17156747126484,2.51941177772377,-0.0638735577336842) q[2];
u3(2.57501017514887,3.71312501431544,-0.563525118487892) q[0];
cx q[0],q[2];
u1(-0.0604165925496445) q[2];
u3(-1.01134566642693,0.0,0.0) q[0];
cx q[2],q[0];
u3(3.27089312754031,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.75581526538717,-1.99540622531198,3.75637361670014) q[2];
u3(1.56911352574025,4.89662197800573,0.495630159576125) q[0];
u3(0.393033690677449,-0.0733889696712671,-0.646834848894127) q[4];
u3(0.420339935995091,-0.845111085088563,-0.646566521404413) q[1];
cx q[1],q[4];
u1(3.44682015230954) q[4];
u3(-1.55979570769354,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.30084470923951,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.91654834423539,1.55868650918532,0.347149983804266) q[4];
u3(0.738505755164445,0.397880981251592,-1.22357266493049) q[1];
u3(1.43887092795415,2.12166284000189,-0.115130405962349) q[2];
u3(2.36270601774480,-0.931123741371693,-4.03989557665273) q[0];
cx q[0],q[2];
u1(3.33849981251758) q[2];
u3(-1.15569297978634,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.37345499025454,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.19021601459308,-1.83688275271009,-0.0128584430746809) q[2];
u3(1.74341798737126,-5.82018054276630,0.341572466596287) q[0];
u3(2.05623169117169,1.26899762957773,-2.42468520140833) q[2];
u3(1.29748282377526,-1.54215029037889,1.78840753855564) q[0];
cx q[0],q[2];
u1(3.84682593764248) q[2];
u3(-3.38497862197556,0.0,0.0) q[0];
cx q[2],q[0];
u3(-1.10299010061085,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.01428373755180,0.0594410455434470,-1.48835248856211) q[2];
u3(1.38404556215881,1.44423696356255,4.80951537832254) q[0];
u3(1.55090402470979,1.67650777969197,-1.47383255454712) q[1];
u3(1.87262858054281,-0.563070561759414,-3.23149565221146) q[3];
cx q[3],q[1];
u1(1.04637277732552) q[1];
u3(-0.774007378769908,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.98486577884814,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.19453989452134,1.33761299906476,-0.956145235979872) q[1];
u3(2.87025046357122,-5.03646448192888,1.12452003351768) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
