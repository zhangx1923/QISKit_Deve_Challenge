OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(0.879251872857972,-2.00553435120666,2.11504939007455) q[4];
u3(0.845896645504499,1.15174465770283,-2.29518867452324) q[0];
cx q[0],q[4];
u1(-0.196139048130689) q[4];
u3(-1.69479521415988,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.05182131150389,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.84485829011299,-2.18177471184375,2.95866427337923) q[4];
u3(1.16778580774668,1.42313064603638,4.09481785805622) q[0];
u3(2.17297375742779,3.61380559795546,-2.24796227526288) q[2];
u3(2.54031540088252,2.17286411183437,-0.914296414650513) q[7];
cx q[7],q[2];
u1(1.46741090137984) q[2];
u3(-1.12603284178444,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.58349937636195,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.562313476659755,-0.838643684472164,3.56791629251845) q[2];
u3(1.38636393017538,0.758949281403436,-1.78916525222058) q[7];
u3(2.61794259857843,-0.782813653478516,2.76514140493985) q[6];
u3(2.49103353583219,-2.57005757678261,-0.376690328456266) q[5];
cx q[5],q[6];
u1(-0.171006655060468) q[6];
u3(1.11492046043302,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.74362328795623,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.19430308160816,3.04564982005342,-1.45219855798928) q[6];
u3(0.983673610241542,-0.321130549089612,0.407595710889658) q[5];
u3(1.93963291441655,0.472332926680333,-2.59539174494661) q[3];
u3(2.37623991501131,-0.156260681748046,-5.08511324978426) q[1];
cx q[1],q[3];
u1(1.14655363304259) q[3];
u3(-0.636027439562130,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.82115916777733,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.25418071590920,-0.732459181413527,0.993339556860924) q[3];
u3(1.08793562277342,-1.54463716345166,-3.63491500173855) q[1];
u3(1.94987102022222,1.98951991929747,-2.78079571313500) q[5];
u3(1.68642900320879,-2.70225046920825,3.10853085097257) q[2];
cx q[2],q[5];
u1(1.78358577102055) q[5];
u3(0.00419797940003441,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.14654563430611,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.400893255321120,3.27109353571707,-1.33974945548655) q[5];
u3(2.53899767744589,0.347766123144559,5.22386037408590) q[2];
u3(0.684288694318713,1.11842096681188,-0.575973201854305) q[1];
u3(0.825194390096307,-1.35093798362747,-0.106918674931636) q[0];
cx q[0],q[1];
u1(0.614436508531026) q[1];
u3(-0.00404389153494988,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.66138320159106,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.04413148869599,3.67951800001926,-1.21880087058199) q[1];
u3(1.88909370141291,-0.314774796930163,0.518081846745595) q[0];
u3(0.876189583349165,1.33572607168777,-2.92279719770645) q[3];
u3(1.26282659175799,2.49095172806200,-3.70189470269705) q[7];
cx q[7],q[3];
u1(-0.780536296894000) q[3];
u3(1.05766244551927,0.0,0.0) q[7];
cx q[3],q[7];
u3(3.78938324125657,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.47910480413628,-3.58509717174286,2.26060480383099) q[3];
u3(1.47918099632552,2.22998912959161,-3.70218041451272) q[7];
u3(2.00830887483885,0.242746240995647,0.450964804046309) q[6];
u3(2.46174615653398,-1.15561046740971,-1.25027378930116) q[4];
cx q[4],q[6];
u1(3.50172266270857) q[6];
u3(-3.93271084328720,0.0,0.0) q[4];
cx q[6],q[4];
u3(-1.07876500060564,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.67927679799676,-3.30639751252839,2.44898325651910) q[6];
u3(1.94104653904111,-1.07118419850979,2.34448249105996) q[4];
u3(1.33161004235042,0.687034586215642,1.87899739496393) q[0];
u3(2.07736829453677,-2.01001904003777,-0.923003894959595) q[5];
cx q[5],q[0];
u1(3.08110833363413) q[0];
u3(-1.28041194099898,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.79816976268017,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.46331201497597,-0.186347628205214,-2.21068283088747) q[0];
u3(1.56058140214712,-1.60242180071372,-2.41887776701955) q[5];
u3(1.40088519303104,-0.831281386764894,-1.35049081508853) q[4];
u3(1.00674339728597,-4.90452639768811,1.23644303211309) q[6];
cx q[6],q[4];
u1(-0.162192538573140) q[4];
u3(-1.96826625713630,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.833729246638007,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.92867715596458,-1.39904819001721,2.20831875265935) q[4];
u3(2.43083081995645,-0.773715846416695,-0.862759563590108) q[6];
u3(3.06589676594073,-1.47775671316224,4.01931368214692) q[7];
u3(0.642256214370665,-1.12539598323794,3.24865285958670) q[3];
cx q[3],q[7];
u1(1.01001226060327) q[7];
u3(-0.255049012798692,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.11618061822513,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.69502131160651,0.0279968464875691,2.05415979558194) q[7];
u3(0.421284923142699,0.630072647164580,-2.00790613691346) q[3];
u3(2.54652766176315,2.40145697534957,-2.95209923441338) q[2];
u3(1.95955006518557,3.21562310218639,-2.66928400553450) q[1];
cx q[1],q[2];
u1(1.34025106587737) q[2];
u3(0.208439900529380,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.662862187358764,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.25310862427699,0.933783682421464,-3.64551727242085) q[2];
u3(2.16732932072556,1.24015004721044,0.432185966013661) q[1];
u3(1.65598543450925,0.165133506436195,1.43414045965739) q[6];
u3(0.802461478247905,-2.26721483579663,-2.13271251352426) q[7];
cx q[7],q[6];
u1(1.39725582495508) q[6];
u3(-3.26173762431604,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.722062754349630,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.19045825526673,-1.19234664826890,-1.93731044614373) q[6];
u3(0.790825853830974,0.00465329029760980,0.612774806199373) q[7];
u3(2.26994234652459,2.03301017349395,-1.98942740204218) q[1];
u3(2.04299241503535,2.43462458808870,-2.40192633208086) q[4];
cx q[4],q[1];
u1(0.444854929691356) q[1];
u3(-1.51709695104876,0.0,0.0) q[4];
cx q[1],q[4];
u3(-0.160429082517124,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.990495102070270,-1.22782836531070,1.29241579825318) q[1];
u3(1.41555617702968,-3.86364342573980,-0.183972221904876) q[4];
u3(1.85568630063338,2.09346760524360,-3.24123309662215) q[0];
u3(2.17129994074215,-2.71541516977368,3.24847847678902) q[2];
cx q[2],q[0];
u1(0.0747362491641659) q[0];
u3(-0.639065399791694,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.73944409840351,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.19122308576382,1.67750493820416,-2.27985586019670) q[0];
u3(1.91583377568723,-0.724959409343898,3.23012711300420) q[2];
u3(0.924709590773350,-1.27029494073940,1.04697088804407) q[3];
u3(1.60424370060868,-2.66905072758341,0.101823342168567) q[5];
cx q[5],q[3];
u1(0.599866775791146) q[3];
u3(-1.16213918636884,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.07538518380599,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.80455554455409,-0.662946432189143,-0.154349448661828) q[3];
u3(0.410801967951694,5.36316828710570,0.837267001589327) q[5];
u3(1.40065933487780,2.95170748502328,-1.63230538693494) q[6];
u3(0.635186662334255,1.94452610556341,-2.77000321078084) q[5];
cx q[5],q[6];
u1(1.80562138647070) q[6];
u3(-2.07103974449077,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.28900488416022,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.66717325446379,2.85488821310434,-1.20291998543457) q[6];
u3(1.51563273677724,-4.32009586380630,1.14450052943536) q[5];
u3(1.45478836810770,-3.63109081517939,2.53007657488672) q[1];
u3(2.14393354544495,3.52951662428101,-2.56649180902682) q[4];
cx q[4],q[1];
u1(1.68217078398698) q[1];
u3(0.313165850979543,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.03758693056626,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.895197096684247,0.0214539994843939,0.403272337444391) q[1];
u3(1.40264520178071,0.326501840622015,0.776949178300565) q[4];
u3(2.53576120487211,-3.11657427119363,0.569230387819746) q[0];
u3(2.24625057989803,-2.39618598268504,0.0816645631640658) q[2];
cx q[2],q[0];
u1(2.36947393199763) q[0];
u3(-1.58795940883987,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.432848058205888,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.01577660625382,2.42600425290185,-3.55200040752296) q[0];
u3(0.548684641876008,-0.276346521870130,-5.28097145107694) q[2];
u3(1.99592608521532,-0.368373768095095,2.72370023931365) q[3];
u3(2.34403456336118,0.530707853821626,2.24188064290569) q[7];
cx q[7],q[3];
u1(2.65345997685863) q[3];
u3(-1.70233231919230,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.213179658087212,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.82893855092966,-3.55036793682705,1.42322600579630) q[3];
u3(2.79472166992639,2.40543039621813,2.40426536733753) q[7];
u3(1.94111855638895,0.00800743936542953,0.643977562760506) q[0];
u3(0.367005653731102,-4.74688566393720,-0.727606036158885) q[6];
cx q[6],q[0];
u1(0.691986753310061) q[0];
u3(-1.23524119887370,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.18075382191632,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.44005450336486,-3.18455783093497,1.41862073650893) q[0];
u3(2.09155650941606,-2.57816225497093,2.48428873966448) q[6];
u3(0.511196632164479,-0.214158902653726,1.36193398125047) q[2];
u3(1.25669079780479,-0.633002627459584,-1.69566269033540) q[7];
cx q[7],q[2];
u1(2.48686048257401) q[2];
u3(-1.83108531144422,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.359303127514427,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.718693543956431,-0.195941236754833,0.378427603918407) q[2];
u3(1.29575700633661,-3.25931980152683,1.24216269120364) q[7];
u3(2.61244833329364,-0.695789635111332,2.29523475246656) q[3];
u3(1.69379482245744,-1.58789892071867,-1.09818326157595) q[4];
cx q[4],q[3];
u1(0.000885460379190794) q[3];
u3(-1.19754432606123,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.46112509834961,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.18007706521303,-1.19670516980827,3.67100867336726) q[3];
u3(0.847534394081116,0.143667105079354,-4.38480568026158) q[4];
u3(0.665012182795240,2.26644088951782,-2.59893590927502) q[5];
u3(1.09496207714437,1.78421217546085,-4.08056148629631) q[1];
cx q[1],q[5];
u1(1.73399539007460) q[5];
u3(0.0803821658148716,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.659637448365208,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.51437046886230,0.334927456637521,1.94069619710234) q[5];
u3(0.699773727301274,0.260474905817988,-3.84059550204593) q[1];
u3(1.93850423613860,-2.64882561355533,-0.279034002402951) q[5];
u3(2.21819466100487,-3.50880009183767,-0.196466646508970) q[2];
cx q[2],q[5];
u1(-0.642840066370646) q[5];
u3(-2.14457114444813,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.72673511377099,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.30388589359795,0.722729056144685,0.659345886680599) q[5];
u3(0.570056235448083,1.89892571732844,0.506545594725433) q[2];
u3(0.847566800906521,-0.657336775057460,0.168640970710191) q[7];
u3(0.680543270733269,-3.05926844471567,1.08896303969679) q[0];
cx q[0],q[7];
u1(0.147358031767659) q[7];
u3(-0.915209779980169,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.96344604417829,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.89181184514391,2.68965019121107,-2.01839930027325) q[7];
u3(1.93172577126218,0.126273563007707,3.10072693423522) q[0];
u3(0.709478655528180,2.21811025785076,-2.93729983749618) q[4];
u3(0.727894284399852,0.661702880108086,-1.89897831145361) q[3];
cx q[3],q[4];
u1(2.30340572025582) q[4];
u3(-1.62123694735708,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.24737021514331,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.20102591526546,0.766905829330033,0.452802339436519) q[4];
u3(0.227575785533663,-1.80002378493046,-3.91844301889463) q[3];
u3(2.56495063322793,-3.50745913077379,2.24364623417187) q[6];
u3(1.36958601091897,2.73644549172860,-0.967033444039656) q[1];
cx q[1],q[6];
u1(1.55238358782315) q[6];
u3(0.217055540707980,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.743488762172056,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.12852579883202,-3.49303515776785,1.60670964784164) q[6];
u3(2.17384907172615,3.35980037803875,1.61650993232428) q[1];
u3(1.80833192811105,-0.477520570535442,2.37109920259690) q[4];
u3(2.54157748711818,-2.05008632164346,-1.21445722354642) q[2];
cx q[2],q[4];
u1(1.51394979497845) q[4];
u3(-0.0624746536638134,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.51479099672975,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.54629553279024,0.0920815252052275,0.185834561295776) q[4];
u3(1.18674379266351,0.996924100124286,2.40886414692474) q[2];
u3(0.250564046057322,-1.71355461268248,0.766802454236479) q[0];
u3(1.04855807457211,-0.234249401210288,-1.93234102671717) q[7];
cx q[7],q[0];
u1(0.720734864916861) q[0];
u3(-3.18154093927222,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.16562443929862,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.62494403915752,-3.80041156625425,2.02180763760359) q[0];
u3(1.21362944270452,-0.566324493177673,1.67772610997446) q[7];
u3(2.00393442477232,1.63857076230019,1.11425696708665) q[1];
u3(1.24226903069004,-0.422471419169925,-3.25917986906786) q[5];
cx q[5],q[1];
u1(0.777668399405771) q[1];
u3(-1.20500068853080,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.59107771970933,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.81564773185924,1.25927292961899,-0.922892276284372) q[1];
u3(0.837210305048861,-0.0264481980200364,3.43339979731702) q[5];
u3(1.13623988165722,4.22388812372165,-1.79465403908698) q[3];
u3(1.88187538203713,1.46197208405064,-2.81427186574194) q[6];
cx q[6],q[3];
u1(-0.184279112076491) q[3];
u3(1.07786513287132,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.81628173943113,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.594954099625665,2.58545196770566,1.83700509823728) q[3];
u3(1.65614516024478,-2.81491967715920,0.535733188389161) q[6];
u3(3.11151834604086,-1.26199402745825,1.50353850042855) q[3];
u3(2.50161084231112,2.33594076355815,3.01357529134281) q[4];
cx q[4],q[3];
u1(2.51565802675335) q[3];
u3(-2.94270112432704,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.95793527886440,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.74219492016089,0.466036101996559,-0.558864934326520) q[3];
u3(1.11145579596638,-0.914529775100998,-2.32211281975622) q[4];
u3(2.07396065335980,-1.74951880441877,1.05703725263109) q[6];
u3(2.73192489342120,-1.81937254358990,-1.63023525793898) q[7];
cx q[7],q[6];
u1(2.58571485610107) q[6];
u3(0.192040696189006,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.25414634085600,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.35648908470265,2.47669837241578,-2.85820213103680) q[6];
u3(1.53788784803309,2.85754571542263,-2.33118485966903) q[7];
u3(1.11276630445297,1.33184623957782,-2.09932897893961) q[1];
u3(0.552951529669731,1.95027299713506,-3.55446465459390) q[2];
cx q[2],q[1];
u1(-0.0315072336816766) q[1];
u3(-2.54772854384327,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.60904584833442,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.73792578005941,-2.40332985125440,3.01101812873627) q[1];
u3(2.44831141714465,-2.50505748734800,3.21094998724420) q[2];
u3(1.52243884996359,-1.03023954250092,1.59999655908334) q[5];
u3(1.51804363964201,-1.24605199425286,-0.792709199292637) q[0];
cx q[0],q[5];
u1(-0.204768318139242) q[5];
u3(-1.79112178951003,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.628456886804623,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.33689099036418,-3.22837224861834,2.84971838747981) q[5];
u3(0.786016185316015,-0.205630318625603,-3.78520608635173) q[0];
u3(0.471582918929577,1.95137298728388,-3.00874470943998) q[3];
u3(0.389824856538579,-0.365170102066835,-1.59302857534260) q[6];
cx q[6],q[3];
u1(2.48503444486786) q[3];
u3(-1.70302345903240,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.181764981336251,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.64088459904594,-0.0124270995800795,3.31046008727916) q[3];
u3(2.17761290074390,-3.86448188708735,2.13006210094980) q[6];
u3(1.97673247445459,1.37825749990449,-0.490142792429485) q[2];
u3(1.85247316187537,0.200061740952298,-3.87932004878453) q[5];
cx q[5],q[2];
u1(1.96206570554785) q[2];
u3(-2.62882841995968,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.00422370688238738,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.27441042907975,-3.92050545607685,1.35868172695227) q[2];
u3(1.11546600251132,0.462735448616692,-0.325999056457234) q[5];
u3(2.00940586705299,0.0988187487658672,0.528048692393314) q[4];
u3(0.720640593505110,-2.63260587911997,-1.41405840499322) q[1];
cx q[1],q[4];
u1(0.416492821509640) q[4];
u3(-1.39478567254803,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.148372072435027,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.52681824528659,-1.61604587064540,2.94413905468542) q[4];
u3(1.97995367336578,0.230242142168381,3.01046283677894) q[1];
u3(1.39413173021906,-2.03661146531089,3.30194934733223) q[0];
u3(1.63613952346828,1.25245602259920,-1.84245153543516) q[7];
cx q[7],q[0];
u1(1.85470329158432) q[0];
u3(0.685680480904306,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.00599680252907,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.629605258978003,-2.88545996100579,3.09984578522161) q[0];
u3(2.53591463943863,-4.18048057184443,-0.234146209936365) q[7];
u3(1.69254969542857,2.81285329794975,-0.902602368452424) q[6];
u3(2.33503255230049,1.34684767356301,-1.31454645363451) q[4];
cx q[4],q[6];
u1(-0.0595805881270999) q[6];
u3(-1.73581003394841,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.590578379000556,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.40857678262048,0.276555154190963,-4.14270132263010) q[6];
u3(2.38859556533963,4.55521689618964,0.0937205658116933) q[4];
u3(1.44814788948696,-0.919945249567221,0.0757461726694721) q[7];
u3(1.49306375858282,-2.54568660281420,-1.21651839329616) q[2];
cx q[2],q[7];
u1(-0.390759953373492) q[7];
u3(1.08434477817419,0.0,0.0) q[2];
cx q[7],q[2];
u3(4.00722510806368,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.63485539556652,0.514399419261351,1.69223377878582) q[7];
u3(1.73726482421992,-4.20331842832620,-0.447378614511335) q[2];
u3(0.697401692575751,-2.78042154449104,1.79847571863178) q[5];
u3(0.837801450424844,0.868223065598329,-2.84694381129631) q[0];
cx q[0],q[5];
u1(0.212865138077478) q[5];
u3(-0.875347408565077,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.59355801114594,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.95127755257307,0.280043098487220,3.43456191813622) q[5];
u3(2.25354052020759,4.33833911157770,0.0726284946949587) q[0];
u3(1.15749339742189,2.66507318013247,-3.32348423520226) q[3];
u3(1.29290030607804,2.68619032765921,-3.21803247138532) q[1];
cx q[1],q[3];
u1(0.229674943907771) q[3];
u3(-1.44632309473031,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.14809763703284,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.49320879039272,-3.12731046600332,1.14768468524855) q[3];
u3(1.12508381211868,-4.51129127353587,-1.38325328191851) q[1];
u3(2.11633314176039,-1.79815829818566,-1.15883187189270) q[0];
u3(1.85517700402574,-2.39293303321078,-0.00515377481544266) q[5];
cx q[5],q[0];
u1(2.27467682178788) q[0];
u3(-3.16770930797836,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.490945928433375,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.27075049947687,0.310660729008474,1.85688359638572) q[0];
u3(1.71230373112264,1.02260269185341,2.80822882503828) q[5];
u3(2.60863147597926,1.20373306434007,-2.42932344785398) q[7];
u3(2.32313437514678,2.05388252152359,-3.62433570398160) q[6];
cx q[6],q[7];
u1(1.62033948009594) q[7];
u3(0.348679794178719,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.827660212023121,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.69790945810675,-3.43506524467386,2.48474102837592) q[7];
u3(0.458761148349081,-2.79960340694413,3.40866408923452) q[6];
u3(1.67851055682236,1.44335382575500,-2.93094595222820) q[3];
u3(1.16410226599862,2.72366227138403,-3.49756316662755) q[2];
cx q[2],q[3];
u1(1.32532259734810) q[3];
u3(-1.02968857968728,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.581544319034373,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.801281075127104,2.74384443332684,-0.102135738454371) q[3];
u3(2.55177428766233,-1.62246250692196,-0.526754492659618) q[2];
u3(2.30794078180038,-1.89129858054336,1.91947196476118) q[1];
u3(2.73499837693322,-1.68853822271758,0.936097049313980) q[4];
cx q[4],q[1];
u1(1.78767566683136) q[1];
u3(0.391707202658213,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.07074666985913,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.808325910948133,-4.10273264818287,1.94144016566532) q[1];
u3(2.57894065968715,1.97841749642301,2.06526269888198) q[4];
u3(0.356316213097588,-1.54003796310010,0.576174764256325) q[4];
u3(1.00833973639193,-0.954546161277570,-0.689569095964914) q[7];
cx q[7],q[4];
u1(2.74112461339024) q[4];
u3(-1.37102316429440,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.05455014557369,0.0,0.0) q[7];
cx q[7],q[4];
u3(3.07125851739204,1.15219960218872,-4.34054452404967) q[4];
u3(0.643171426481484,0.866530354804456,-1.45383816997155) q[7];
u3(1.28584051402212,1.43592707385878,-3.63894964965499) q[2];
u3(2.55176689387171,2.64810924934400,-2.65740085000343) q[3];
cx q[3],q[2];
u1(1.95212103979135) q[2];
u3(-1.73524927807210,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.83743155639764,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.49452241838813,3.56218681488952,-0.664178153568437) q[2];
u3(0.761324124124615,2.57576436995006,-3.52767441136037) q[3];
u3(1.39922046189816,0.240088768233169,1.22569762138166) q[5];
u3(1.39809895222037,-1.06206410227684,-1.57063192202605) q[6];
cx q[6],q[5];
u1(0.361980934380812) q[5];
u3(-0.812521265287208,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.15621511480785,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.44866861320431,-0.0874792235706390,-0.647674805490405) q[5];
u3(2.32392285710421,3.37848161024044,-0.851212936681425) q[6];
u3(2.05551980018512,0.963886620638293,-0.0895754441012009) q[0];
u3(1.69244895826753,-0.706585725619893,-4.01320536669112) q[1];
cx q[1],q[0];
u1(1.24029363572504) q[0];
u3(-0.682754744807564,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.0753938350251357,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.84477392461512,0.0623823637305725,-0.464332503862588) q[0];
u3(1.27773887883919,-4.45694252370191,0.770100453626382) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];