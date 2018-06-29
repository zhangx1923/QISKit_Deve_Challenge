OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.06287664297086,-1.28614359509035,0.206210670768917) q[10];
u3(1.81245958803158,-2.53138493049470,0.387176939107008) q[7];
cx q[7],q[10];
u1(3.62194752350922) q[10];
u3(-3.18319484465413,0.0,0.0) q[7];
cx q[10],q[7];
u3(-0.821284355592616,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.507953144947083,1.71855769144085,-2.53030956228717) q[10];
u3(0.659229471062082,0.118559006003369,-1.64509363074432) q[7];
u3(1.62069070244823,-2.98926886968242,-0.0381530167413697) q[2];
u3(1.26852077512273,-3.26451983183164,1.23902146712965) q[11];
cx q[11],q[2];
u1(2.08965453339756) q[2];
u3(0.340285356372324,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.989880948730926,0.0,0.0) q[11];
cx q[11],q[2];
u3(0.682039563503984,-0.258968211393971,1.61462437356889) q[2];
u3(1.39355150583751,-4.50272657882103,0.137598217988422) q[11];
u3(0.494621621189021,1.76158562249899,-0.454008783620753) q[0];
u3(0.625407363408443,1.09942008151887,-2.87805134091685) q[5];
cx q[5],q[0];
u1(0.256679437606215) q[0];
u3(-1.60896207553698,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.93985795626953,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.68260279739149,0.284412154587723,-1.65431034926078) q[0];
u3(1.47196123242416,0.858651323692806,-0.916800577567766) q[5];
u3(2.38268068200968,1.79810570473873,-2.44630023592005) q[9];
u3(2.03513976283086,2.53597912547060,-3.47935798073290) q[8];
cx q[8],q[9];
u1(1.53568823684565) q[9];
u3(-1.27131651238101,0.0,0.0) q[8];
cx q[9],q[8];
u3(-0.647125387467177,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.21363559711122,0.258386203814562,-2.46205056185776) q[9];
u3(0.370581839145295,-1.46439590080063,-1.24226758522682) q[8];
u3(2.24088254619244,-2.16608409604692,0.0563120654823825) q[1];
u3(2.23474095524974,2.22677131822807,3.68342595031385) q[3];
cx q[3],q[1];
u1(1.92773790193284) q[1];
u3(0.295901051754483,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.704212343007284,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.860751977352501,-1.68503939864155,2.85972750710992) q[1];
u3(1.66620688713026,2.77889216093737,0.0972270765589225) q[3];
u3(2.10375981065246,2.64292965615526,-1.72607881013844) q[6];
u3(1.17486827438457,0.827335376005495,-0.935732742352232) q[4];
cx q[4],q[6];
u1(0.00237771307860046) q[6];
u3(-2.18256881330721,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.820006486001147,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.41320957403654,-2.09043525863808,2.84010818256445) q[6];
u3(1.62600275628509,-2.64168128053630,-0.755367893440374) q[4];
u3(0.617295538152270,-0.699546860347922,0.879123623260499) q[7];
u3(0.234484121025686,1.66771882196929,-2.36154041841567) q[0];
cx q[0],q[7];
u1(1.58174082865411) q[7];
u3(0.0289140307413767,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.53699895800681,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.28938865298958,-2.27144553933537,1.01145786532860) q[7];
u3(1.07434569431253,1.33076917663754,0.0554883522968670) q[0];
u3(1.67379352670951,1.05238569444842,-1.62376374634304) q[11];
u3(0.464948557160650,0.474835284982163,-3.59644013745290) q[5];
cx q[5],q[11];
u1(-0.335404807412577) q[11];
u3(-1.79386145049072,0.0,0.0) q[5];
cx q[11],q[5];
u3(0.827941645046275,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.17219286063147,0.688069136528005,-2.11161350880298) q[11];
u3(1.95385065041612,-3.72969491541587,2.07718067054820) q[5];
u3(1.01078542738503,-0.292909465404015,2.71916842740665) q[3];
u3(1.43466324136951,-2.64862959335378,-1.42040851762761) q[8];
cx q[8],q[3];
u1(1.41299997976257) q[3];
u3(-2.75273899397674,0.0,0.0) q[8];
cx q[3],q[8];
u3(0.587089191966032,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.348555824494700,1.62332540406852,0.298015599211179) q[3];
u3(1.83366105239836,2.65589716824681,-2.66132372427972) q[8];
u3(1.43838643055297,1.62813222482481,-0.134854416152296) q[10];
u3(2.16368207207855,0.0691168502452302,-3.20519116766108) q[1];
cx q[1],q[10];
u1(0.288864457377238) q[10];
u3(-0.615574756524216,0.0,0.0) q[1];
cx q[10],q[1];
u3(2.26688307560950,0.0,0.0) q[1];
cx q[1],q[10];
u3(0.414623609737310,2.80528489590295,-1.47792525220107) q[10];
u3(1.41759652237821,-1.32214366758083,-4.02711843071158) q[1];
u3(1.25686059091740,0.741680335571516,-0.884131812548164) q[4];
u3(0.891458610628032,0.0912123335664847,-2.42453099410372) q[2];
cx q[2],q[4];
u1(1.77799228336753) q[4];
u3(-2.87078386101008,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.803524277608634,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.675925834055714,-0.893365079379094,1.67753402576272) q[4];
u3(1.73851838341456,-0.579899872842929,1.87636766360321) q[2];
u3(2.23231796244618,1.35399446863118,1.14821394798086) q[9];
u3(0.890014789655658,-1.27738609128064,-2.36808395611331) q[6];
cx q[6],q[9];
u1(0.205256885593258) q[9];
u3(-1.38627198726786,0.0,0.0) q[6];
cx q[9],q[6];
u3(2.45284816417554,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.43595272730712,-2.02681274421888,0.246866592568894) q[9];
u3(2.08453280152492,0.837571507562237,4.42758646287872) q[6];
u3(0.298965890599003,-1.29664224898215,1.19580901424382) q[1];
u3(1.16011553495076,2.78711646521874,-3.14914035002851) q[4];
cx q[4],q[1];
u1(0.0248165502038287) q[1];
u3(0.892067563187585,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.45149712698677,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.57809164380638,2.70138500823754,-0.711876495959885) q[1];
u3(1.50794206272049,1.13615527056823,1.61038453148891) q[4];
u3(0.839420972098550,2.03819851920043,-3.65236779615715) q[10];
u3(1.78361487974919,-2.35220503239381,3.62765345719854) q[11];
cx q[11],q[10];
u1(1.68200459689469) q[10];
u3(-3.05201594176064,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.81866543231316,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.77482714098623,0.766372159326782,-2.61171780057049) q[10];
u3(1.88060418137802,1.48109340265673,-4.50000071443702) q[11];
u3(2.51603378509989,-2.93322915215618,1.32299466877237) q[0];
u3(2.36897317355926,-2.62269273372339,-1.01872144007274) q[6];
cx q[6],q[0];
u1(0.115030961220673) q[0];
u3(1.06235930856951,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.85264798425252,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.96175985242936,-3.14455976095999,1.23675719216043) q[0];
u3(1.66264721811567,5.33169425708837,0.755159477667074) q[6];
u3(2.21044548668482,-0.810933557180724,1.38602065141657) q[8];
u3(2.03846157325557,-1.29435598732371,-1.15049678061858) q[9];
cx q[9],q[8];
u1(1.59671422241355) q[8];
u3(-1.06284248057082,0.0,0.0) q[9];
cx q[8],q[9];
u3(2.68691092961715,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.686620846511931,-1.77252080811857,3.81329775676629) q[8];
u3(2.64471610281931,-3.03091380645665,2.12158275783395) q[9];
u3(1.03719473960385,1.29867008151029,-1.77431796627215) q[5];
u3(1.94369565215060,-2.62955588944917,3.60747422283929) q[2];
cx q[2],q[5];
u1(1.56638715877369) q[5];
u3(-0.654356938521066,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.60553229092985,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.12194969036386,2.64682962925966,-2.60351280959099) q[5];
u3(1.90563003552800,-0.603968304745822,-1.13677444589292) q[2];
u3(2.62607098022478,0.755747221366919,0.990632762113940) q[3];
u3(1.76763026783672,-1.70112637736202,-1.53670677230268) q[7];
cx q[7],q[3];
u1(1.97587470386183) q[3];
u3(-0.0979770485342462,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.834067060356372,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.74720838720590,0.546475706379282,-0.906113760511650) q[3];
u3(1.51031996311552,5.53814656081170,0.178092633421586) q[7];
u3(1.53285883467758,0.319710136429503,1.93663953129825) q[7];
u3(1.34693402882967,-1.36232698544895,-2.14079523221374) q[1];
cx q[1],q[7];
u1(1.75262565746219) q[7];
u3(-0.487321268005186,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.91896765851683,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.52391602553293,2.80902562856088,-0.871410105106981) q[7];
u3(1.44205528509792,-0.769598196334699,2.68610569167578) q[1];
u3(0.806989178760824,0.912411824235210,-2.67677804480034) q[6];
u3(1.58685166123389,-2.35612292768119,3.06656411924947) q[4];
cx q[4],q[6];
u1(-1.21089149235140) q[6];
u3(0.452978712244963,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.79815239628863,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.37616960118962,-2.38312405178914,1.95261814183797) q[6];
u3(0.715780675763194,-3.75564318332700,-2.42486036078250) q[4];
u3(1.34395927717809,1.19139267131689,-1.68827551525582) q[3];
u3(0.674424456367774,-1.27773287889857,0.0335758232750240) q[8];
cx q[8],q[3];
u1(1.96152208324111) q[3];
u3(0.0604570375250839,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.67436256836661,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.60562823709340,-0.0271341917131795,2.09596616967852) q[3];
u3(2.78120715905256,-1.11015524332190,-0.667760443677723) q[8];
u3(1.90671708562577,1.79531312632095,-0.999903888496730) q[10];
u3(2.60743432780247,-0.626674403696023,-4.03246225719016) q[5];
cx q[5],q[10];
u1(0.948926862854574) q[10];
u3(-0.731128094907602,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.50115128142022,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.45345295460117,0.578074415894928,-3.24966991939780) q[10];
u3(1.10391753783880,5.10454661729379,0.537590934861052) q[5];
u3(1.86351339184577,-1.83119867923186,0.426691202659015) q[11];
u3(1.41992098962302,-3.90109257971570,1.17574507201018) q[2];
cx q[2],q[11];
u1(2.14881203139835) q[11];
u3(-3.06028057972464,0.0,0.0) q[2];
cx q[11],q[2];
u3(0.274217207020144,0.0,0.0) q[2];
cx q[2],q[11];
u3(2.09159294840270,2.53189530664122,-1.42531084423661) q[11];
u3(1.67809634817407,-2.81322937020897,-2.70460137676960) q[2];
u3(1.69020319472925,1.56601882844192,-0.545198124937803) q[9];
u3(0.152909498206735,-3.72548895464833,1.61955498295253) q[0];
cx q[0],q[9];
u1(3.00780668454475) q[9];
u3(-1.30954746341932,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.00342866483831,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.14275040609433,3.50697971165747,-1.37295226009903) q[9];
u3(2.24592133680357,-0.426890268407571,4.62464368699411) q[0];
u3(0.644763234157600,2.04098865980801,-2.31797120436651) q[6];
u3(0.640814060382808,0.590294733726508,-3.24301729635643) q[4];
cx q[4],q[6];
u1(3.12509209489606) q[6];
u3(-2.32099129646246,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.77691851843582,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.02311012015942,-0.627976382314842,-2.92785152828596) q[6];
u3(0.603751534941125,3.94013623021764,2.28633072533056) q[4];
u3(1.81443277570416,2.33699287498906,-0.950976752656214) q[1];
u3(1.07574931714799,0.455444023107267,-3.43652569522183) q[11];
cx q[11],q[1];
u1(3.14758814018399) q[1];
u3(-1.92338062732972,0.0,0.0) q[11];
cx q[1],q[11];
u3(0.612077516484122,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.78249181380907,-2.57324885397609,1.82189640033319) q[1];
u3(1.07845162576773,-3.11650742824778,-0.163900531029953) q[11];
u3(1.62893084151682,1.51910603873776,-3.22595592678750) q[5];
u3(1.09226517961598,-2.37920906517843,3.20717634705080) q[3];
cx q[3],q[5];
u1(3.10855933577596) q[5];
u3(-0.700156726852112,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.82099564591674,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.50311024534717,2.73178556379951,-3.38822671242720) q[5];
u3(2.29938707876204,0.577998733941463,0.657548867733400) q[3];
u3(2.65516652265138,0.706898896376997,-0.679777244787871) q[2];
u3(1.71100022573216,-0.0131907530337181,-3.25733187760313) q[9];
cx q[9],q[2];
u1(3.35293747764553) q[2];
u3(-2.04232367007741,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.852488258990452,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.36981678664273,-3.62008834785176,2.35386107593280) q[2];
u3(2.35523457335455,-0.609875591294199,1.28460038801563) q[9];
u3(1.89878211402815,3.46633529048047,-0.617942921562783) q[0];
u3(2.10058685213515,1.51182776619504,-0.959134972911497) q[10];
cx q[10],q[0];
u1(-0.0320019910761693) q[0];
u3(1.16788108604328,0.0,0.0) q[10];
cx q[0],q[10];
u3(3.47017871810926,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.24626661625729,3.41491549146266,-0.680127130252571) q[0];
u3(1.92537303313650,-3.62529953271350,-0.575264534417286) q[10];
u3(1.16874931459810,1.05177344293620,-3.17507583847119) q[8];
u3(2.28134972073478,-3.58922753367608,2.47463319621916) q[7];
cx q[7],q[8];
u1(-1.49526121518078) q[8];
u3(-0.419779873543781,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.75031837337837,0.0,0.0) q[7];
cx q[7],q[8];
u3(3.05298854710278,-2.82262026198461,2.44856073439678) q[8];
u3(1.95618227743440,2.44332319187060,2.79365863410596) q[7];
u3(0.532396205140781,1.22153737594941,-3.91478518069027) q[7];
u3(1.67626317455537,2.89968299046824,-3.26818230386968) q[9];
cx q[9],q[7];
u1(2.02122189792070) q[7];
u3(-2.85455526676941,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.616629660082807,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.90275199702551,-1.76494259516371,1.08684036627894) q[7];
u3(2.37516247911314,-1.55215382430876,-1.46453728402956) q[9];
u3(0.917927953814879,-3.99927365899091,2.12339097800260) q[6];
u3(1.44780638414469,2.85668694610061,-2.47664728800177) q[5];
cx q[5],q[6];
u1(1.54954689324089) q[6];
u3(-2.48124992347606,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.0463698107545105,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.52828030728577,-0.970630094737041,0.443219137154339) q[6];
u3(1.24454101809840,1.91423895220235,0.540171489390431) q[5];
u3(0.439297899399925,0.0608819221440067,-0.383755531588764) q[0];
u3(1.20052343991691,-1.14739672385967,-1.05273836049734) q[11];
cx q[11],q[0];
u1(1.99371974466856) q[0];
u3(-3.10414756705755,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.673111537957729,0.0,0.0) q[11];
cx q[11],q[0];
u3(2.08177113383438,-1.63248297671642,-0.487010952144528) q[0];
u3(3.07674682831303,0.560101084930746,2.49616685802836) q[11];
u3(1.22594462610715,1.34461251640257,-0.576058543889298) q[8];
u3(2.54926492544989,0.318522783914318,-2.88991429480862) q[2];
cx q[2],q[8];
u1(1.63145052621318) q[8];
u3(0.302842843177145,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.538077584128830,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.793572513861060,-2.33875689175532,2.34831253135917) q[8];
u3(2.70379870586103,-1.68162391230016,-0.654723929745170) q[2];
u3(0.879873198065846,-0.419989427740955,0.632822719526496) q[1];
u3(1.19363048585116,-0.474848328804437,-1.33217154071264) q[10];
cx q[10],q[1];
u1(0.741467904203050) q[1];
u3(-3.17749411137294,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.49811712695169,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.15030222343749,-1.27340919400145,-0.669003517009938) q[1];
u3(1.15568839978534,0.987543510067526,4.40417816568329) q[10];
u3(0.120206986253469,-3.79836603668856,2.33708179660108) q[3];
u3(1.09208910886986,-1.76797668208520,0.244483861223376) q[4];
cx q[4],q[3];
u1(0.420591848191788) q[3];
u3(-1.68900238712385,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.205175695731920,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.90820851839050,-1.99766445464961,1.13067477698159) q[3];
u3(2.31943806497830,-3.90722931513372,-0.0948357604325931) q[4];
u3(1.89753216492010,0.988283847197577,-3.31683656039972) q[2];
u3(1.27805578498996,2.69376594921431,-2.78312980759921) q[7];
cx q[7],q[2];
u1(1.45827133831580) q[2];
u3(-0.0624489611310111,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.743143198531045,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.750775792344485,-2.10730845330810,3.68944757008497) q[2];
u3(2.66423685352517,-1.32179792897758,-1.53265226384069) q[7];
u3(1.09144617808186,0.990006136114949,-1.95865239591164) q[11];
u3(0.556412416097329,-1.03945567163237,-0.520890845914269) q[6];
cx q[6],q[11];
u1(0.526480227614764) q[11];
u3(-1.24469612231179,0.0,0.0) q[6];
cx q[11],q[6];
u3(2.97605405996553,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.32347696393041,1.63834460172303,-2.93238998933936) q[11];
u3(2.64745861764855,2.89529088316866,2.99592202144016) q[6];
u3(1.46257554006816,1.46358081390420,-2.54657503245825) q[4];
u3(1.88246831780608,-1.72065063922994,3.28960994089050) q[3];
cx q[3],q[4];
u1(0.943496587696853) q[4];
u3(-1.38586300580883,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.22016315050655,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.951763988269366,0.0289364707525552,-1.81681530804059) q[4];
u3(1.40194666966213,-0.696272850786903,-5.54544642719139) q[3];
u3(1.08061801098759,-1.49479333133644,0.631715087654783) q[0];
u3(0.176152425390580,-0.881026890968979,-0.456537149616518) q[9];
cx q[9],q[0];
u1(1.68189566300488) q[0];
u3(-2.21559751089952,0.0,0.0) q[9];
cx q[0],q[9];
u3(-0.109328881158652,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.62471144275393,-0.874326011290451,3.74229509997224) q[0];
u3(2.27700709500214,0.327760109978108,-1.35716156381520) q[9];
u3(1.50535374579102,0.763394846479936,-3.50806322263167) q[5];
u3(1.32790944666062,4.51844992108273,-1.70266193888685) q[1];
cx q[1],q[5];
u1(0.0484538755063877) q[5];
u3(-0.584832698745579,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.17390507735705,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.926463421830842,3.52967503809254,-1.24270765591631) q[5];
u3(1.70537965641091,-4.86996454296664,-1.05960309177985) q[1];
u3(2.59410206151788,-0.770566211925511,-1.75942144722372) q[8];
u3(1.15416267149512,0.114390793715332,-5.10116433492716) q[10];
cx q[10],q[8];
u1(1.86550278468104) q[8];
u3(-2.38851895885754,0.0,0.0) q[10];
cx q[8],q[10];
u3(3.47805275088708,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.78209918657474,1.82967670961606,-1.64076235886481) q[8];
u3(1.18427214361550,2.08343889739859,-0.180649923096048) q[10];
u3(1.36804354636539,-1.80954252026134,1.83106288550350) q[7];
u3(0.273903895451325,2.33111065463413,-3.60780336235366) q[1];
cx q[1],q[7];
u1(1.16942231643252) q[7];
u3(-0.0344372915918747,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.10768045702515,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.24687307086205,1.27786966450884,-1.20001890262137) q[7];
u3(1.34358783371295,0.459815404164697,-4.46963254090216) q[1];
u3(1.89675216690517,2.97975107015133,-1.55031709096918) q[9];
u3(0.843918075186485,0.775137856455651,-0.325324805434765) q[2];
cx q[2],q[9];
u1(3.71760234857195) q[9];
u3(-1.35506650164976,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.01925554609809,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.12225071389653,2.48952583740427,-1.52940926893856) q[9];
u3(1.57863805870197,4.08332160320261,1.21457194978142) q[2];
u3(1.62132427572067,1.98928379019213,0.868089175289124) q[0];
u3(2.68802801435430,0.547736650052770,-2.60235231951555) q[11];
cx q[11],q[0];
u1(3.36022517709524) q[0];
u3(-1.38035080243446,0.0,0.0) q[11];
cx q[0],q[11];
u3(2.23684425414963,0.0,0.0) q[11];
cx q[11],q[0];
u3(0.909483951295847,-1.21500024619583,-2.92501155638538) q[0];
u3(1.70997774779098,-0.446215115580573,5.65995175129732) q[11];
u3(0.970410426024598,-0.0236945141830152,-1.61449954630191) q[5];
u3(1.52419634383208,0.199753933872447,-5.42765536190884) q[10];
cx q[10],q[5];
u1(2.42287031902133) q[5];
u3(-1.66583503071007,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.979468168322468,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.72791521854724,-0.333033620355525,0.322324295815882) q[5];
u3(2.74800695667057,-0.551774909314585,2.03707792338032) q[10];
u3(1.02986699208075,1.48093948416894,-2.96431585661154) q[3];
u3(1.48242891902860,2.61973827305428,-3.18309503823610) q[6];
cx q[6],q[3];
u1(1.36448784162773) q[3];
u3(0.0420490252823467,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.533646992112845,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.695717494598706,0.190841266522313,-1.59825693000354) q[3];
u3(1.51090768319065,2.00546282652283,1.33872185363679) q[6];
u3(0.939496995120476,-2.04410480688669,2.38441530327758) q[8];
u3(1.02638180885559,-2.30941733502415,1.79201467604160) q[4];
cx q[4],q[8];
u1(3.10276592989983) q[8];
u3(-1.83779895849022,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.683755920935247,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.61661459560051,3.04782539628822,-1.17400234844164) q[8];
u3(0.899135449658777,1.19119107278137,-0.651148049241416) q[4];
u3(1.40742337640294,-0.772247820103646,-1.07421206116463) q[1];
u3(2.36806048883320,-5.28053360362601,0.627102157881249) q[10];
cx q[10],q[1];
u1(2.20381904111441) q[1];
u3(0.349275364240684,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.43801805853420,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.42299082750192,-2.10921399482939,1.45477466499174) q[1];
u3(2.67521777684756,1.55020347736892,1.53220156489320) q[10];
u3(1.95307720287542,2.08072945875422,-1.23013517061580) q[3];
u3(2.57473743916475,1.38339667408038,-1.69990533380269) q[9];
cx q[9],q[3];
u1(1.97102534441362) q[3];
u3(-1.79384856554567,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.02023310529530,0.0,0.0) q[9];
cx q[9],q[3];
u3(2.07745117276902,1.62526720752086,0.112392947626853) q[3];
u3(0.646106148609976,-1.25531405421111,-1.54268299117135) q[9];
u3(1.09800253838221,0.416040201246144,-0.151469854819530) q[8];
u3(0.331734938774017,0.308091966249100,-2.85939618104816) q[5];
cx q[5],q[8];
u1(2.10885275435533) q[8];
u3(-2.86693934955862,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.550049717397104,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.990871478002817,4.19285066305054,-1.14008675187200) q[8];
u3(2.60231370571568,-4.01744620777603,0.324040161936484) q[5];
u3(2.76377656747662,2.14359557074477,-0.509435074736023) q[7];
u3(1.75774602363483,1.82051422477782,-4.24489271961507) q[0];
cx q[0],q[7];
u1(0.0603142859279713) q[7];
u3(-1.02113166203027,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.68978588031244,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.26023320801079,-1.26654076962517,-0.570948554317323) q[7];
u3(2.26302066014755,-0.198329872257085,1.06628451116361) q[0];
u3(2.29260493345894,-2.52651437229703,-0.149943579328365) q[4];
u3(2.74443447338709,2.15130848519458,3.10502284468554) q[11];
cx q[11],q[4];
u1(1.75172409140932) q[4];
u3(-2.38581294447406,0.0,0.0) q[11];
cx q[4],q[11];
u3(3.14931945236108,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.26307750626133,-0.272084215261827,-1.34120014166792) q[4];
u3(1.78607695394193,3.33629375582728,-1.55397933208737) q[11];
u3(2.09028730775568,0.386877878531672,-3.13554295304982) q[6];
u3(1.44802587288835,-3.30219480359560,2.69018290020221) q[2];
cx q[2],q[6];
u1(1.70812856269029) q[6];
u3(0.162323735158397,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.64521274905503,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.67712838967198,3.52468113597392,0.849361671126053) q[6];
u3(2.00871868501643,-2.06305300874541,1.88659108955398) q[2];
u3(2.44487788715748,1.67303179664097,-0.400408566469144) q[7];
u3(1.34849187684828,-0.952275161566300,-1.72015376373791) q[9];
cx q[9],q[7];
u1(2.90641849684325) q[7];
u3(-2.03202774839879,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.604722306180381,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.22425780792179,2.79157530633006,1.52553340144864) q[7];
u3(1.10836942007476,3.14917673371044,-0.508140295413163) q[9];
u3(2.72667185757186,-0.913224493030515,2.55180582063541) q[6];
u3(2.80224049117527,-2.83962539984099,-2.14678902508449) q[10];
cx q[10],q[6];
u1(1.44276069435342) q[6];
u3(0.132155520335391,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.01318784626358,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.74497539348042,-0.0768023663351804,-2.49242423177348) q[6];
u3(1.95464510769187,0.210802351402142,1.70169260842669) q[10];
u3(1.01150492550067,1.36820754349339,-2.44271093039096) q[0];
u3(1.91979899492077,2.22058802794745,-4.01196130358675) q[11];
cx q[11],q[0];
u1(3.18346876542987) q[0];
u3(-1.42098189256672,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.460724071232522,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.25743556768368,1.09528111207907,-2.32729487531246) q[0];
u3(0.700544117805904,0.992843596308982,-4.30888019896453) q[11];
u3(2.65866985804539,0.868618192279897,0.423054672002093) q[3];
u3(0.813993036293095,-3.74018245119097,-1.30903162025706) q[2];
cx q[2],q[3];
u1(1.68921216282001) q[3];
u3(0.241167325901624,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.782251077488025,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.25612406374193,4.20138267797965,-0.0382752519425322) q[3];
u3(2.24988536284604,0.140774158376099,-0.672480184737954) q[2];
u3(0.161554323217086,2.57783226675797,-2.60453602919972) q[4];
u3(0.529881850755755,-0.391981932968932,-0.0294574218256003) q[8];
cx q[8],q[4];
u1(0.425908477681536) q[4];
u3(-1.57007575521260,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.41142513099264,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.03985512916190,-0.585637097505829,2.88596165132622) q[4];
u3(1.40133971473735,-0.640677388158747,0.851033307052076) q[8];
u3(0.818874240725457,1.61715420265721,-2.48405763567716) q[5];
u3(1.63008803928256,2.25831178212708,-3.92840408518966) q[1];
cx q[1],q[5];
u1(2.57444510532289) q[5];
u3(-1.53982021452746,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.177617156196092,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.01186488255156,-0.943693387513048,2.52248554212254) q[5];
u3(2.43061257680589,1.96093580528367,0.126040258603973) q[1];
u3(1.26321662912363,1.76075452097680,-0.264167679445452) q[1];
u3(0.849658383504919,0.0363344413138809,-3.91111242523298) q[11];
cx q[11],q[1];
u1(-0.432302608499425) q[1];
u3(1.22687863871684,0.0,0.0) q[11];
cx q[1],q[11];
u3(3.43344422479266,0.0,0.0) q[11];
cx q[11],q[1];
u3(0.648700450917357,1.37775038090777,1.85527340758675) q[1];
u3(1.79922903790089,-0.416656285900582,2.43503416982670) q[11];
u3(0.751441578481061,-0.199120588695465,-0.588067119711957) q[9];
u3(1.32175524013653,-3.82352282507827,2.20021435775827) q[6];
cx q[6],q[9];
u1(0.819775920322111) q[9];
u3(-1.53994700795045,0.0,0.0) q[6];
cx q[9],q[6];
u3(2.70041713304492,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.15039220981149,2.28967799744928,-0.748746404120834) q[9];
u3(1.55336809827037,-3.60882950697893,0.626774305779339) q[6];
u3(1.92239628897865,0.150486087945510,1.78578360689625) q[5];
u3(1.59283434937618,-3.07510915749106,-2.04825049315850) q[4];
cx q[4],q[5];
u1(0.476514409925033) q[5];
u3(-0.0361052999936562,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.72858690666560,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.135015256371612,-1.99526456285209,3.14905532575065) q[5];
u3(0.0911732578048487,-0.0661900426247359,-2.77739585219037) q[4];
u3(1.57159068559992,-0.0316768340171045,1.53720649412044) q[0];
u3(1.60567338120533,-1.91508247493048,-1.51118777752413) q[7];
cx q[7],q[0];
u1(1.38545723606619) q[0];
u3(-2.42057948625330,0.0,0.0) q[7];
cx q[0],q[7];
u3(-0.0106052055347545,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.668578727519215,-2.45824996552340,-0.0859890789019944) q[0];
u3(2.42494130909823,0.140888375014570,-0.0668417541995168) q[7];
u3(2.04860529346198,0.0945454759899014,2.37090132829745) q[3];
u3(0.782932213879492,-0.734291075472238,-1.72105852698657) q[10];
cx q[10],q[3];
u1(1.76015822453174) q[3];
u3(-2.58089888332155,0.0,0.0) q[10];
cx q[3],q[10];
u3(3.21281666520970,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.48068754874619,3.49619952858910,-2.42687211228936) q[3];
u3(2.42906059236963,0.586521476729413,3.29903136280381) q[10];
u3(2.08878618437999,-0.238654641200449,0.469433499472665) q[2];
u3(2.16579507384301,-1.28335921217650,-1.90573944824825) q[8];
cx q[8],q[2];
u1(1.62577356668319) q[2];
u3(-0.393882928690225,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.24369363981974,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.66024966782310,-1.32587694053307,-0.652636345742399) q[2];
u3(0.347506246070126,1.93819110189424,0.539792130550949) q[8];
u3(0.906340048727773,-0.532971259765525,0.986259063999679) q[9];
u3(0.599794893242082,-1.98746326247900,-1.29871840016698) q[11];
cx q[11],q[9];
u1(1.37521287660480) q[9];
u3(0.147626938617892,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.29661400449456,0.0,0.0) q[11];
cx q[11],q[9];
u3(0.482173880095521,1.60908631012253,0.256527713472804) q[9];
u3(0.797728960476743,-1.75744986510403,-2.70369908191056) q[11];
u3(1.57449458894729,0.669510400428050,-2.66474063807735) q[6];
u3(1.79612423704648,1.73377216360028,-4.14334636416944) q[3];
cx q[3],q[6];
u1(0.198913872282712) q[6];
u3(-1.01409469281069,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.48226589420966,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.67157106980329,2.10772124044688,-0.114886718118000) q[6];
u3(2.61189548302395,0.0290073241050603,0.833903264494565) q[3];
u3(1.10134978426946,2.32664248155294,-3.16198084815359) q[4];
u3(0.612745548934717,1.69716869433856,-2.45331331641982) q[1];
cx q[1],q[4];
u1(-0.275046154731197) q[4];
u3(1.10272633202900,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.01258975039878,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.57675038962557,-0.0818564446352044,2.68521423995188) q[4];
u3(1.21591900996738,-0.882396125713646,-1.17063107323221) q[1];
u3(1.44333437559898,3.49188411856865,-1.33250506087134) q[0];
u3(0.928624747025336,2.59815795338893,-2.92291781002767) q[2];
cx q[2],q[0];
u1(2.12248445728794) q[0];
u3(-2.18966960684936,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.14625305766247,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.147280690052025,-0.976806890661076,4.71994947410100) q[0];
u3(1.42725996129935,-0.673731263041782,3.13488028977151) q[2];
u3(0.900040027904472,-1.21289208270356,2.08090366983960) q[8];
u3(0.598150256716575,1.92897740600286,-3.04697905784719) q[10];
cx q[10],q[8];
u1(1.33768198672546) q[8];
u3(-3.46873736034295,0.0,0.0) q[10];
cx q[8],q[10];
u3(2.29072664244169,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.0787184920841323,0.348340716645219,3.25230497010097) q[8];
u3(1.81920012953763,-3.02809477121641,2.12413336061947) q[10];
u3(2.49205838437751,0.856052960246877,1.72689570258165) q[5];
u3(1.76731781580729,-2.71395792947560,-2.52892644883977) q[7];
cx q[7],q[5];
u1(-0.248464424937263) q[5];
u3(-1.07049817467627,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.80466391867236,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.23421904371033,-1.38964368339637,0.632849621369160) q[5];
u3(1.36853761536540,2.29846134423932,0.899889786543672) q[7];
u3(2.39778631299305,0.107747591067554,-2.02018956953030) q[4];
u3(2.75488280651783,-0.789648420284365,-5.08616256428626) q[5];
cx q[5],q[4];
u1(1.74846190581437) q[4];
u3(-2.00140302378462,0.0,0.0) q[5];
cx q[4],q[5];
u3(3.05935567659138,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.855685269483341,-0.362257656145876,-1.80271664527719) q[4];
u3(1.27613480606809,2.55248433314035,3.57082604392707) q[5];
u3(1.03147025818280,0.0164214001409456,-0.604264815593314) q[9];
u3(2.13156451030996,-3.85167318862203,1.09866314818393) q[6];
cx q[6],q[9];
u1(0.922683388103632) q[9];
u3(-1.23351594607041,0.0,0.0) q[6];
cx q[9],q[6];
u3(3.35397438766498,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.93768533112349,1.61914172168614,-0.476566456987801) q[9];
u3(0.739268714904220,-1.09927520377658,3.88984356556449) q[6];
u3(0.416738164334926,-1.14010112568366,0.567089518928062) q[2];
u3(0.379741535699957,-1.36944772043487,0.0786412383343340) q[3];
cx q[3],q[2];
u1(-0.227223988627625) q[2];
u3(-1.48198447504101,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.983354007164447,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.464009565526585,-0.848970771097389,-1.05249928394657) q[2];
u3(1.31803842607351,4.69229397858335,-0.394003972860426) q[3];
u3(1.02113961563651,-1.65851809131988,1.63419026665551) q[8];
u3(0.807575577549760,-2.38273219429630,0.0837440248877086) q[0];
cx q[0],q[8];
u1(0.608234434142914) q[8];
u3(-3.28649480760620,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.74012490053144,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.23865323096717,1.64102378351738,1.74078214216762) q[8];
u3(1.64811720084520,-0.0324193788871415,-1.68926215577330) q[0];
u3(2.54328213196857,-0.505459071732776,-1.05615083492925) q[11];
u3(0.764769925025683,-0.465030484320424,-3.99523848195897) q[10];
cx q[10],q[11];
u1(0.309292134908509) q[11];
u3(-2.01172132137211,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.28092404329574,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.82276206561142,1.88711799575919,-3.20452648171621) q[11];
u3(2.08327018333165,-0.179474936149387,-4.21560877549907) q[10];
u3(0.983785428203357,1.72315612130108,-0.183881316558001) q[1];
u3(1.32415529253552,-0.441473751527485,-4.18237086127029) q[7];
cx q[7],q[1];
u1(-0.268148165523450) q[1];
u3(-2.48125780323659,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.37975035357625,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.52025243745764,-1.55930010531984,2.90023603058849) q[1];
u3(1.37493296895964,1.12634678709694,-3.03005050006946) q[7];
u3(1.63805949060241,1.09538089951844,1.30065800226989) q[1];
u3(0.696373296585924,-1.33035152495454,-2.74419863796122) q[4];
cx q[4],q[1];
u1(0.265412735336464) q[1];
u3(-0.789417346349644,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.50883763528202,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.05061546054533,0.307888711346570,-2.41517397102050) q[1];
u3(2.50054828475132,1.80202759214308,0.266502288021385) q[4];
u3(2.65909437611437,3.94122113241566,-1.90135600309193) q[8];
u3(0.699334520985169,-0.494309776157437,2.00027997684428) q[2];
cx q[2],q[8];
u1(3.49775411182920) q[8];
u3(-1.38380585527252,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.13624058881797,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.450549146820624,-0.509533198446113,-1.18587960994302) q[8];
u3(1.60805445761663,-2.74086107544703,3.24395998501097) q[2];
u3(1.40448337218088,1.63815305161632,0.814071755468750) q[5];
u3(2.76734599887895,0.133611887990504,-1.97120115346550) q[9];
cx q[9],q[5];
u1(3.70551654264974) q[5];
u3(-1.53785389025073,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.83174949397781,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.37878777805410,-4.04759439238777,1.22700154079850) q[5];
u3(1.97020803037534,-0.580755474590712,-1.97091830893029) q[9];
u3(1.74426189035793,1.19261601084222,0.406951717673573) q[0];
u3(2.55904516154850,-0.391601413780575,-3.15305372740943) q[3];
cx q[3],q[0];
u1(3.13156335293901) q[0];
u3(-1.85433887286213,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.06510092505347,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.90138175159495,3.07717969806958,-0.991375464303518) q[0];
u3(1.99959756102095,2.58169029454386,0.702056625729546) q[3];
u3(2.91646050196291,-2.54625148864532,-0.0414861654142837) q[6];
u3(1.93128089275703,-3.56663203797280,-2.40996076248353) q[11];
cx q[11],q[6];
u1(1.30249331501755) q[6];
u3(-0.566518190471466,0.0,0.0) q[11];
cx q[6],q[11];
u3(3.01837714027680,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.38369517718521,-2.48203250398008,-0.956234436235281) q[6];
u3(0.487797179298681,-1.82919179974582,3.47089869689239) q[11];
u3(0.713491747496595,-0.0227215253969103,0.170226078484800) q[7];
u3(0.867911665459793,-0.666048271067287,-0.646228259538985) q[10];
cx q[10],q[7];
u1(2.94575979234714) q[7];
u3(-1.68925461667518,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.833420277371928,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.29353789319624,-0.765768056026627,3.27344367565887) q[7];
u3(0.723303098595152,0.570998681375636,3.94534310205060) q[10];
u3(1.44018728166872,1.72800680619756,-2.89790392646098) q[0];
u3(0.478018930897319,-1.74417016049659,2.08641065015638) q[11];
cx q[11],q[0];
u1(1.44559920343955) q[0];
u3(-2.67381958342953,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.557266675573487,0.0,0.0) q[11];
cx q[11],q[0];
u3(2.63422695317055,-1.27169314990803,-1.63032737341202) q[0];
u3(2.13054833546975,-1.72742299257224,4.55188356988059) q[11];
u3(2.35558159159983,-3.78747051792475,2.33429956814393) q[7];
u3(0.343287153898850,-2.38310381042393,3.74318530212561) q[1];
cx q[1],q[7];
u1(1.15876667119650) q[7];
u3(-0.611153950159340,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.82762568437218,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.824858333165306,0.113456088881061,1.87133261569889) q[7];
u3(2.84097537043811,-4.38218297791234,1.55338945146553) q[1];
u3(2.12165821425521,1.27696236883505,-1.57252329351453) q[2];
u3(1.96020543941851,-4.47337253100210,1.53771040971253) q[4];
cx q[4],q[2];
u1(1.86172072370331) q[2];
u3(-2.88399751933279,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.793060353934702,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.01487957196002,-1.65063305233794,0.750444865533957) q[2];
u3(1.89342612173017,-1.85787610445388,-2.62063995706781) q[4];
u3(2.09707038473172,-1.52311818435893,0.387525956913385) q[6];
u3(2.38987334704242,-1.64128301335047,0.361827834130394) q[9];
cx q[9],q[6];
u1(0.623609642431834) q[6];
u3(0.0619380818657305,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.15928555053214,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.04876353720127,-1.62139991846211,0.956649977576551) q[6];
u3(1.25460484841858,-4.10426117763993,2.06673823329818) q[9];
u3(2.62682265377400,-0.660296909832585,2.28960189245115) q[8];
u3(2.16067666315240,-1.34448547231206,-1.34785643054938) q[10];
cx q[10],q[8];
u1(1.63630232164906) q[8];
u3(-3.07680222552350,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.574557656835358,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.807158625053898,-1.18593878099646,-0.835113717493980) q[8];
u3(0.690327956717825,-0.285611996071697,5.69746967073904) q[10];
u3(1.37101487951088,1.71766038570570,0.326604885306538) q[5];
u3(2.40861847452463,-0.261640336854827,-4.49397169810359) q[3];
cx q[3],q[5];
u1(0.330173884779323) q[5];
u3(-1.40909759744306,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.99504348644067,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.26112261663674,-0.470386451928992,1.61540445433109) q[5];
u3(1.08757528815646,-1.24326923918191,0.374326191146135) q[3];
u3(1.08390429477841,1.13295913517956,-2.79740769725380) q[3];
u3(1.72469556616032,-3.87456751130570,2.03354115908824) q[8];
cx q[8],q[3];
u1(-0.0300349361210723) q[3];
u3(-1.97770223146755,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.50040002667832,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.51712819917166,3.02105449725622,-0.314692550796476) q[3];
u3(2.03070702164008,-1.00410377738758,3.36457080960595) q[8];
u3(2.90152798980832,-2.76939471591816,1.29663495369781) q[4];
u3(2.41154315842573,-2.14592106736197,-0.800055272647539) q[1];
cx q[1],q[4];
u1(1.85840467143877) q[4];
u3(-3.08752497026130,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.840530499667013,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.24770304204067,3.30862488537232,0.0944325820810494) q[4];
u3(1.61848987912525,5.32601575629652,-0.388199007585128) q[1];
u3(1.88579060219903,1.18852342887430,-2.64443284594130) q[7];
u3(1.69203719322815,2.03752040390325,-3.43978981714424) q[11];
cx q[11],q[7];
u1(1.76583720676356) q[7];
u3(-2.27996234464224,0.0,0.0) q[11];
cx q[7],q[11];
u3(3.70127596208453,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.58718680437298,0.955104052244943,-3.61641351829886) q[7];
u3(2.63805600612772,5.45351677727505,-0.415443101804888) q[11];
u3(1.01411525322338,2.26401619332478,-2.52628208425422) q[0];
u3(0.645517908188058,1.68232963671696,-1.80504828709023) q[2];
cx q[2],q[0];
u1(1.86092702152807) q[0];
u3(-0.0461939774991291,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.883392813422432,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.52064255916500,2.22903106528251,1.49497563963882) q[0];
u3(2.13454169820566,-6.16222729280125,-0.0670717532920810) q[2];
u3(1.35494271953263,-0.127370441704242,1.78587859810203) q[9];
u3(1.24585588263592,-2.43813009851334,-1.90563599224736) q[5];
cx q[5],q[9];
u1(-0.130348030646568) q[9];
u3(-1.70538906100673,0.0,0.0) q[5];
cx q[9],q[5];
u3(0.877292348070476,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.72234586496945,1.14362001767553,-1.26085252912660) q[9];
u3(1.86791292448359,-1.68470858322701,1.44856316278219) q[5];
u3(0.764678738361882,-2.94481532896630,2.87482778742970) q[10];
u3(0.859084083818818,1.05469427528596,-2.84600528183490) q[6];
cx q[6],q[10];
u1(1.20320125386823) q[10];
u3(-0.0101061638633697,0.0,0.0) q[6];
cx q[10],q[6];
u3(2.50557867830008,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.91990805997359,0.533040440017847,-2.34523254312258) q[10];
u3(1.33333585403917,-3.14662830063452,-1.92743632307981) q[6];
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
