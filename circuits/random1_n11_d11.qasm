OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.27214478980588,-0.488213495275528,1.47176276704169) q[10];
u3(2.16383974254351,-2.49717140736832,-0.435745913606302) q[8];
cx q[8],q[10];
u1(1.80268331609028) q[10];
u3(-2.50443699167937,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.151487020427273,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.67595325217833,0.716297009895002,1.64791422079768) q[10];
u3(1.73032055891896,3.76285339917830,0.562725597382284) q[8];
u3(2.42729680027580,-1.20509276399875,-1.45979611813332) q[5];
u3(0.808742624178128,-5.12054616884789,0.444068320532776) q[6];
cx q[6],q[5];
u1(1.52347311333409) q[5];
u3(-0.949581684993841,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.06862670424355,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.868406924225480,4.64615116475843,-1.46492021350463) q[5];
u3(2.14525612115819,-3.02606291163122,0.676026162546636) q[6];
u3(1.66378684899434,3.07005275566112,-3.02763676209012) q[7];
u3(2.42742933460389,2.10413673905624,-1.24404589702969) q[1];
cx q[1],q[7];
u1(0.00915696622224282) q[7];
u3(-0.938923540053955,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.94343220397808,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.71796317584709,-1.83139639618636,-1.05189649374390) q[7];
u3(1.72030082884326,-0.211193473221935,-0.875319113941161) q[1];
u3(1.77914112569702,0.875282721766262,1.30403196436148) q[4];
u3(1.15764435900553,-0.855750745295262,-1.92210540845014) q[0];
cx q[0],q[4];
u1(3.16624136876938) q[4];
u3(-1.52857355617326,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.560893778769636,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.42215268019395,1.76897197073552,-4.47004149085912) q[4];
u3(1.42930893323978,-1.49012226498997,3.33139502268645) q[0];
u3(2.43705088171686,2.83773908763637,-3.00954213082867) q[9];
u3(1.67394013982971,1.61269300435306,-1.76051789233968) q[2];
cx q[2],q[9];
u1(1.48546713944495) q[9];
u3(-0.382468862469874,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.09100035237209,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.483015940990108,-2.67893783832934,-0.178295819055205) q[9];
u3(1.37589606720667,3.01960221221851,-0.452966626008692) q[2];
u3(0.0930415855594568,0.0941737813717431,1.27747188292348) q[9];
u3(1.27204099700657,-3.12241694132840,1.55922175419028) q[5];
cx q[5],q[9];
u1(3.20391505500770) q[9];
u3(-0.778515766020678,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.61984520834742,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.15952235270195,1.46103967539912,-2.30581218372049) q[9];
u3(0.831736987321519,1.01009734086417,-0.0191092928926786) q[5];
u3(0.948726776543145,0.494424071852726,0.363738268074433) q[2];
u3(1.14206008736828,-1.45498530408155,-1.19617038761291) q[4];
cx q[4],q[2];
u1(1.50299178894626) q[2];
u3(-0.922506815074748,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.53788463288474,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.75882485231169,2.89169098829370,-1.35243766465958) q[2];
u3(1.96521869476323,0.787579856845731,-3.52961588535608) q[4];
u3(1.85970085541201,0.704070111066464,1.12485307565938) q[6];
u3(1.44271969323551,-0.952972565833389,-2.11240083702251) q[3];
cx q[3],q[6];
u1(0.336166866576392) q[6];
u3(-1.60324651204676,0.0,0.0) q[3];
cx q[6],q[3];
u3(3.00849015501151,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.27268536409905,2.15331280386089,1.81571796972134) q[6];
u3(1.71456372898524,-0.700442042442121,-0.210155793646096) q[3];
u3(1.79609969588830,-0.219884732582524,1.47920554128776) q[7];
u3(1.61038837715639,-2.30869323851503,-2.59842520183391) q[8];
cx q[8],q[7];
u1(2.23831947130736) q[7];
u3(-1.98973035366938,0.0,0.0) q[8];
cx q[7],q[8];
u3(3.52584140991934,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.54516040000756,2.17602773484553,0.820417453802850) q[7];
u3(2.18768538753067,0.337451605728604,4.85632505378307) q[8];
u3(1.95886641425023,1.27800569852051,-0.0853866758639020) q[10];
u3(1.69868817763347,-4.65443310966399,1.44608308692187) q[1];
cx q[1],q[10];
u1(3.78693380734200) q[10];
u3(-4.23334383114363,0.0,0.0) q[1];
cx q[10],q[1];
u3(-0.355023458413236,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.25256947896866,-1.21580641704599,2.66270589700181) q[10];
u3(1.80706664138428,-4.32021748456951,-1.19670366242166) q[1];
u3(2.67964833592324,0.967089637593402,1.90769471035306) q[2];
u3(1.44342849973070,-2.01272318140161,-3.05918324736695) q[0];
cx q[0],q[2];
u1(1.84068242687394) q[2];
u3(-2.15505330817405,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.0162922139272321,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.489848844132130,-1.38111790348220,0.644236919177120) q[2];
u3(2.11984473414856,-3.31976479315516,-2.85161860312023) q[0];
u3(1.36949102350449,2.02082699609347,0.102593921060071) q[10];
u3(1.56917792638172,0.517597143151987,-3.88607184825105) q[4];
cx q[4],q[10];
u1(1.44172883803710) q[10];
u3(-3.43775789652582,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.56104890478728,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.440709375408399,3.51712424586646,0.342123146466475) q[10];
u3(2.45527770885189,-3.31481208138506,0.839005002492840) q[4];
u3(1.50272859683332,-1.43228786050796,-0.655627267850714) q[7];
u3(1.05110808967500,-2.83473320818308,-0.222933123822618) q[5];
cx q[5],q[7];
u1(-0.995075306576188) q[7];
u3(0.290899287000314,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.72876186310134,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.10360697051340,-2.66835650877042,2.81790091317724) q[7];
u3(1.51319508200788,-0.0542917432789529,-3.45960903618074) q[5];
u3(0.819598660601930,-0.232465639201389,0.0413509483607731) q[6];
u3(0.482041174265935,-1.11650653048321,-1.09748449814411) q[8];
cx q[8],q[6];
u1(0.572253023838951) q[6];
u3(-1.75785978607493,0.0,0.0) q[8];
cx q[6],q[8];
u3(3.05422275912289,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.35729754065241,-3.04455480962216,1.73366693864220) q[6];
u3(0.734301266271814,0.885895091528553,4.77704457348071) q[8];
u3(1.95653023872577,2.55906540129998,-2.38472020221997) q[9];
u3(1.77784518368745,2.58978524923182,-3.32240909243173) q[3];
cx q[3],q[9];
u1(1.00773411555323) q[9];
u3(-3.16989814448428,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.65716731540968,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.75746843486186,-0.590733315068753,5.27266123881395) q[9];
u3(0.525088884703586,1.99868939949066,0.539175851968725) q[3];
u3(0.817714954164810,-2.92145641110820,0.346188617914797) q[0];
u3(1.15552666876618,-2.78009370303398,-0.953831674431376) q[6];
cx q[6],q[0];
u1(1.51983092557880) q[0];
u3(-0.637798195677606,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.18007144230864,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.90916091519021,2.04247753042017,0.380396170541042) q[0];
u3(2.37668946194147,1.28642033960937,-0.0728878091257370) q[6];
u3(0.752271482412564,0.820926437681924,-2.16914651779004) q[9];
u3(1.32278651099360,-4.19083886647411,1.82315848342863) q[3];
cx q[3],q[9];
u1(2.54074118985621) q[9];
u3(-2.93178943832040,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.55509394411716,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.54106643262168,1.93315421868923,0.131371487215371) q[9];
u3(2.08266070707867,-1.96427717252313,0.505115797752830) q[3];
u3(1.91057265109146,-0.343421779348766,0.584161239161022) q[2];
u3(2.08849444713934,-0.209190550309632,-1.59720716900581) q[10];
cx q[10],q[2];
u1(1.97611600109770) q[2];
u3(0.238787972279260,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.954596973110004,0.0,0.0) q[10];
cx q[10],q[2];
u3(0.789399706279929,3.18963272916562,-1.94444720828227) q[2];
u3(0.844476701765333,-0.365343707994559,0.548335503059559) q[10];
u3(0.421293425341878,2.44450477872503,-1.48266168526831) q[7];
u3(0.894432251953907,-3.35182018844321,1.23822295409443) q[5];
cx q[5],q[7];
u1(1.07542273793708) q[7];
u3(-3.82271135748661,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.67609160361807,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.634323108112657,-3.77820800020688,1.97447899278286) q[7];
u3(1.13078012949865,-1.53843446076317,-0.946142224831949) q[5];
u3(1.19314686449991,2.38809511722997,-2.58632063847834) q[1];
u3(0.279665338696390,1.76938331188070,-2.31786241819824) q[4];
cx q[4],q[1];
u1(-0.0221937611359204) q[1];
u3(-1.44937017262299,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.325323722292535,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.70212775961074,2.09421566383794,-2.05190887678245) q[1];
u3(0.686814950863907,0.873680656680220,3.30757769181354) q[4];
u3(2.48654357233907,0.566279019428696,-2.33972049156680) q[2];
u3(2.26701447727747,2.19336565820266,-3.97712088602903) q[6];
cx q[6],q[2];
u1(2.97271703769752) q[2];
u3(-2.54220235510807,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.40875569828973,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.81834739742130,3.29662794977472,-1.88362412204285) q[2];
u3(1.00300131650850,0.0284508391315152,5.15893397568041) q[6];
u3(2.04421694456051,1.56061992486368,-0.263548436760116) q[0];
u3(2.14023485582478,0.160275222357592,-1.86373285976419) q[7];
cx q[7],q[0];
u1(2.24066281289528) q[0];
u3(-2.55984977532625,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.37359212016631,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.43222249422243,1.77800729924572,-1.05318969831264) q[0];
u3(0.986345348587194,3.21970065952860,1.82209442521408) q[7];
u3(1.27184189614277,4.01327943131977,-1.85121298892246) q[3];
u3(1.94590180577146,1.75716803662435,-2.76084000077352) q[4];
cx q[4],q[3];
u1(-0.250377086360566) q[3];
u3(-1.71815159049651,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.631206789883208,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.95338551656832,-0.996814353737718,-1.12016137140285) q[3];
u3(0.208462872358050,1.79578475168682,-0.562714923165418) q[4];
u3(1.58999220607481,-2.75379999668983,3.49823754886376) q[5];
u3(2.23206910681074,1.80913450584562,-1.36049728447406) q[10];
cx q[10],q[5];
u1(2.31327977087312) q[5];
u3(-1.69572653750984,0.0,0.0) q[10];
cx q[5],q[10];
u3(3.82079969490258,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.58385592676173,2.72694414538784,0.106951382423240) q[5];
u3(0.301545232208501,1.58868821899932,1.31820221153779) q[10];
u3(1.05499603780963,0.509162508185256,-2.06461141712694) q[8];
u3(1.93051715460785,2.09459652114010,-3.59146243106235) q[1];
cx q[1],q[8];
u1(0.771065812707034) q[8];
u3(-0.164550063289751,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.57184150674828,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.601298515982644,3.55805713720407,-1.28632856297749) q[8];
u3(2.90771859211967,2.77615182364818,-1.59070307703729) q[1];
u3(1.67653844010447,1.05382271458544,-3.95723056944405) q[4];
u3(2.25452766720114,3.63386826904083,-1.96297943809173) q[1];
cx q[1],q[4];
u1(2.37082305290516) q[4];
u3(-1.62291363364651,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.44107873803848,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.90343260468129,2.30266012738603,-3.95258460663310) q[4];
u3(1.37283320905673,2.31516531866863,0.849255703218130) q[1];
u3(1.70481349551629,0.821272637956063,0.884019288853493) q[0];
u3(0.579591160374697,-2.12875357125131,-2.55867688201443) q[6];
cx q[6],q[0];
u1(3.96376290998670) q[0];
u3(-3.72687027041550,0.0,0.0) q[6];
cx q[0],q[6];
u3(-0.508380487461962,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.86162390075343,-1.36058833354761,0.140004054793935) q[0];
u3(1.87779627339629,-0.0203319211784723,-3.61252030735353) q[6];
u3(1.79690996117550,0.668830685851115,-3.54241072769689) q[10];
u3(0.876988004341552,-2.28664792334705,2.59734481611048) q[7];
cx q[7],q[10];
u1(3.26807235042007) q[10];
u3(-1.13535664078749,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.17992733992603,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.810198286811296,-0.799384682861780,-0.706343032523314) q[10];
u3(1.15821394716605,-1.67376820752173,-3.50160902144914) q[7];
u3(2.12609476673024,2.08528486219059,0.839414146164655) q[2];
u3(1.43198952345806,0.609379394628673,-2.25470087878829) q[3];
cx q[3],q[2];
u1(-0.105349342638962) q[2];
u3(-1.68034295401397,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.492934552107856,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.39272000571669,-1.31035534580898,0.752016653747223) q[2];
u3(0.953675324252459,4.07247560857521,-0.927810879304631) q[3];
u3(1.72119493694074,1.71005663830974,1.09428211810655) q[5];
u3(2.27838934126673,0.434202873848330,-2.77593521066765) q[9];
cx q[9],q[5];
u1(1.50296634966753) q[5];
u3(0.0249270283732712,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.80408250849823,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.959515892187063,-2.32869646783754,0.280427799875614) q[5];
u3(0.405216241450896,-0.441148050351925,4.44357340618655) q[9];
u3(1.16288075109057,-0.491701813784496,-0.124735543254469) q[6];
u3(1.01045885927110,-3.33767948631256,-0.374420655054426) q[0];
cx q[0],q[6];
u1(2.47008737124889) q[6];
u3(-2.94981275537787,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.61768668085036,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.810612160301298,2.09013090703870,-2.86239244886009) q[6];
u3(1.37260666918128,-3.47002029828819,-2.16184264780156) q[0];
u3(1.79379744407608,2.02730677527190,-3.45039400418958) q[1];
u3(2.70932817791400,3.09570549194203,-2.77378668906278) q[9];
cx q[9],q[1];
u1(2.70681308187547) q[1];
u3(-2.03655943292640,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.562711427407918,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.81991044333604,1.52811622843325,0.119325413675632) q[1];
u3(2.91910501704145,-3.28465379440018,-2.13214346579680) q[9];
u3(1.04762702382481,1.74527928220712,-2.82467426095855) q[10];
u3(0.920554641266439,2.20223195088883,-2.60789026667325) q[8];
cx q[8],q[10];
u1(3.28326864340602) q[10];
u3(-2.35975460267758,0.0,0.0) q[8];
cx q[10],q[8];
u3(1.28616939346356,0.0,0.0) q[8];
cx q[8],q[10];
u3(0.827473178566568,-3.46343212861922,2.01462928272880) q[10];
u3(1.45662714070593,-0.729460016519959,2.30331587692250) q[8];
u3(1.01637017540922,1.99841987065660,-2.97933323378568) q[2];
u3(1.10802551543893,0.899276328265303,-1.70409441771132) q[4];
cx q[4],q[2];
u1(4.01449384512799) q[2];
u3(-4.28732555048418,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.649416938157699,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.60098378713329,-1.66059868826307,0.457873733616706) q[2];
u3(2.18353057567933,0.596157203874505,0.889771364102646) q[4];
u3(2.77010212922830,-3.56698096300707,2.41278458007588) q[5];
u3(1.29200416490509,-0.345606590420448,2.15566613483466) q[3];
cx q[3],q[5];
u1(1.28395378241032) q[5];
u3(-0.747296134618201,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.10563990309700,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.30172372243285,2.71444285624021,-0.602285512307490) q[5];
u3(0.686008233988512,-1.84815413681885,4.33282740512609) q[3];
u3(0.909883222439699,-0.914344898630382,1.67318660400118) q[9];
u3(0.334522451961013,-2.49068910827352,0.667890099331478) q[3];
cx q[3],q[9];
u1(2.97597382555105) q[9];
u3(-1.30809459891579,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.45392617961494,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.810041291857866,2.18257346977828,-2.96784560298221) q[9];
u3(2.84416987993246,-0.833312174235793,-1.08588737312359) q[3];
u3(0.902348600423469,0.924324166349537,-1.16823954043286) q[6];
u3(0.577484017205804,-1.14991860234358,-0.0629487097532432) q[0];
cx q[0],q[6];
u1(2.02685775453227) q[6];
u3(-1.81459779441220,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.92243643495152,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.788446705117047,2.26699888611391,-0.813309888465507) q[6];
u3(0.466553887302992,2.13808960063602,-1.36259016552854) q[0];
u3(1.80659810388314,4.33691495802752,-1.19572456332264) q[7];
u3(1.68998740571526,4.00896521395730,0.0146326459265937) q[10];
cx q[10],q[7];
u1(0.878384244563576) q[7];
u3(-3.37609018663138,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.00659297504111,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.250788094431635,1.80402879571212,-0.733715308134335) q[7];
u3(2.84725789571020,-5.69027454614106,0.277639922374705) q[10];
u3(0.529177011636150,0.593885967029527,2.41329441376808) q[4];
u3(1.46562484989744,-1.46958166009497,-1.67850527213739) q[1];
cx q[1],q[4];
u1(1.13835311448764) q[4];
u3(-0.116628736155495,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.28383180911201,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.44164436868972,0.928084918300598,-4.81507148496645) q[4];
u3(2.74384698684387,2.19660425077749,-0.0551690015512298) q[1];
u3(1.23991060567888,3.41423605869218,-0.524693948564863) q[8];
u3(1.37199769994524,0.575213329981103,-1.58260371839996) q[2];
cx q[2],q[8];
u1(0.0829023743798749) q[8];
u3(-1.12456701099875,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.42507873567538,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.63474591402545,1.19643091187558,0.987973233552440) q[8];
u3(1.42861955740906,2.39145559639798,-3.73774624210449) q[2];
u3(1.01802265624592,-0.191152955195421,1.62309758157172) q[0];
u3(1.54293532404697,-1.93734628873474,-2.60277171939294) q[1];
cx q[1],q[0];
u1(1.46571529816990) q[0];
u3(0.0750394191211814,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.84528094212247,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.04233653277454,2.05450701860702,-0.821981535351814) q[0];
u3(0.788744423117762,4.00226960708195,0.0180314647482513) q[1];
u3(2.12604213576791,-2.91056827011508,3.15796209923784) q[5];
u3(0.526249333661389,2.00786267014914,0.309896217951846) q[9];
cx q[9],q[5];
u1(-0.239440434798061) q[5];
u3(-1.73123482137989,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.09256945562029,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.59203477327176,0.464103231349762,-0.165017282663672) q[5];
u3(2.37778347214576,-0.606484193136291,-5.15916708495696) q[9];
u3(1.57763641981334,0.882895579228729,-1.06442549881279) q[10];
u3(0.348094133218113,1.56062426611958,-3.92991190809056) q[3];
cx q[3],q[10];
u1(0.934130732710973) q[10];
u3(-3.27351313957549,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.81179916290591,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.16418321815744,-0.0332872311605187,-0.157635714043576) q[10];
u3(2.27969767412656,-0.526542148227354,3.08541791596147) q[3];
u3(0.923652555550424,1.61367987151121,0.766293917616997) q[4];
u3(1.61634641503471,0.0560083704589200,-3.65211828232974) q[7];
cx q[7],q[4];
u1(2.58324633527341) q[4];
u3(-1.97870330248772,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.762479775464335,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.86721531196209,-2.65835452565612,-1.40292663191769) q[4];
u3(1.78099253933234,-4.78551117901907,0.563990424933738) q[7];
u3(1.82268609127910,-0.406958319949012,0.908935288144712) q[2];
u3(2.46092700353471,-1.48503067458795,-2.42966317412548) q[8];
cx q[8],q[2];
u1(2.20289353677718) q[2];
u3(-3.16504228336538,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.79674447444765,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.88962368593032,4.71323886272145,-0.393614883164028) q[2];
u3(0.748454718865846,-1.77815637252571,-2.01835975808552) q[8];
u3(1.65393179184408,3.14285933910264,-1.61485860078793) q[5];
u3(2.26659334951325,0.587720698882682,-1.28269374750429) q[0];
cx q[0],q[5];
u1(-0.803673946132013) q[5];
u3(-1.78954233128853,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.07389015453303,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.55308253786839,-1.66978581100260,2.91011495831021) q[5];
u3(2.88979639755663,4.68245763362667,0.662585537067574) q[0];
u3(0.883190581399700,-1.40837674386800,-0.784819338337690) q[2];
u3(0.869089988601695,-3.87738788229050,0.324113894689427) q[8];
cx q[8],q[2];
u1(3.93706240133698) q[2];
u3(-3.51315722283491,0.0,0.0) q[8];
cx q[2],q[8];
u3(-0.305284343280788,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.25910525014710,0.444390624397083,-1.78318484881192) q[2];
u3(2.25117936019620,-0.428189919308530,-0.105865392816049) q[8];
u3(1.65440192611756,1.13376402869682,-2.72829458869347) q[3];
u3(2.73380187809089,4.26967912536594,-0.935178943105412) q[4];
cx q[4],q[3];
u1(-0.465839131174798) q[3];
u3(-2.02904979459095,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.42832485539758,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.68978442584406,0.554927814767037,-4.04245994382677) q[3];
u3(0.828780315333247,-2.20516332978002,-3.01133205527021) q[4];
u3(2.69722208719619,-3.28978764849675,0.217928085235202) q[1];
u3(2.26131788998036,0.801369467341487,2.78433836708047) q[10];
cx q[10],q[1];
u1(0.301763427986276) q[1];
u3(-1.28937141595650,0.0,0.0) q[10];
cx q[1],q[10];
u3(2.53180828736793,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.07520484041646,2.78956045837945,-1.18719553071692) q[1];
u3(1.42867241243761,2.13256715481853,-0.563362919980261) q[10];
u3(0.496652694235000,0.160658294479410,-0.182051993721008) q[7];
u3(0.467554726161151,-1.65181256178268,1.42862323531516) q[9];
cx q[9],q[7];
u1(1.29263473438374) q[7];
u3(0.225537981686657,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.48100426967041,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.826222473438619,0.664525817957926,-2.00359925291266) q[7];
u3(1.11587707781227,3.17772167236965,2.14961156455578) q[9];
u3(2.11909185855360,2.59601060040928,-0.436772455771599) q[5];
u3(1.65081066234830,2.13763502153049,-1.21569463678623) q[4];
cx q[4],q[5];
u1(0.0233494624410777) q[5];
u3(-1.02682432132078,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.85806028056496,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.75633442022368,-2.92876840867340,1.69726814326486) q[5];
u3(0.826655899795468,4.96340597377450,-0.680838377400177) q[4];
u3(1.46945869187940,0.357153580361040,2.24399382778053) q[3];
u3(1.33116177780067,2.71278945795769,3.42708702735726) q[0];
cx q[0],q[3];
u1(1.81127611749143) q[3];
u3(-3.24387470218467,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.795636570587997,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.40986790332799,-2.89858672750404,0.617392534926410) q[3];
u3(2.40668628095786,1.47855362190967,1.26636336420212) q[0];
u3(0.955876732402097,-0.989262024567373,1.65019251316322) q[8];
u3(0.858416306496407,-3.52190209276895,2.69629520712699) q[2];
cx q[2],q[8];
u1(3.65541407136870) q[8];
u3(-3.29417639649582,0.0,0.0) q[2];
cx q[8],q[2];
u3(-0.837733042101958,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.31965045119654,-0.482430711036626,2.07643416784374) q[8];
u3(1.55416412837782,-2.13028025891811,0.0557923753354244) q[2];
u3(1.56183913885734,1.36938400561635,0.240139724782050) q[1];
u3(1.43186297266341,-0.247234862342267,-1.78148395907632) q[7];
cx q[7],q[1];
u1(3.78194438677548) q[1];
u3(-3.66756648956888,0.0,0.0) q[7];
cx q[1],q[7];
u3(-1.03081020559631,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.10997678310248,-1.73753791545408,2.38213098192795) q[1];
u3(1.53195207682057,0.652766992487252,-3.06126548788790) q[7];
u3(2.18095629424120,0.357127898058636,-1.89012736236593) q[6];
u3(2.21326401707256,-3.69278104263394,2.14344523123247) q[9];
cx q[9],q[6];
u1(1.19719976683488) q[6];
u3(-2.94602447577104,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.30719240304056,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.998059401422811,-0.376495249284131,3.27804587507636) q[6];
u3(0.0786568514256029,1.44273273979969,1.26947870261506) q[9];
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