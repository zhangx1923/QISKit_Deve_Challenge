OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.17063535861763,1.72309540329942,-0.172526895529872) q[3];
u3(1.25263647772667,-0.253786567469419,-3.95446921269026) q[2];
cx q[2],q[3];
u1(3.96860944420105) q[3];
u3(-4.30943708961146,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.608780436288784,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.59594086004756,-0.623138770487716,4.45378296643608) q[3];
u3(2.67505830703467,-3.95620437450971,0.245617794775590) q[2];
u3(0.993712661743854,2.06313529344837,-3.28913420356428) q[5];
u3(1.15133306287183,2.23496285367209,-3.55620592732640) q[4];
cx q[4],q[5];
u1(0.0277960077174995) q[5];
u3(-2.50237313329996,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.04681955583715,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.93436017755999,-2.17115463056574,-0.857245247606439) q[5];
u3(1.02172243271653,-1.10063296137644,-1.20717688401593) q[4];
u3(2.18938129771272,0.554354551852751,2.36658700114549) q[1];
u3(2.68940199436760,0.416649724906039,2.18716001452621) q[0];
cx q[0],q[1];
u1(1.74397696678226) q[1];
u3(-2.29437375823524,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.46193841863895,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.35664773112842,0.373284414846039,-0.886146634470977) q[1];
u3(2.56792498021814,0.619207284300016,1.10772862567080) q[0];
u3(0.610041261871451,1.09296202441089,-2.57145947171879) q[3];
u3(1.74093922693949,-2.59841739768786,3.25648119809522) q[1];
cx q[1],q[3];
u1(2.22364732905075) q[3];
u3(0.0584888490054425,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.20152079611000,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.46664241544601,-3.85222023188780,1.97540745982903) q[3];
u3(2.26323200251116,1.74275495364323,3.08815700089343) q[1];
u3(2.47763020671774,-0.238721036975298,1.65291511545485) q[5];
u3(1.95576071430983,-3.00104966586622,-1.97860646603262) q[4];
cx q[4],q[5];
u1(2.45955220238170) q[5];
u3(-3.05840647085963,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.21819852003683,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.27871228296693,-1.52843611574268,2.71322102389555) q[5];
u3(2.14226400918231,1.54806945181525,-1.34626177086375) q[4];
u3(1.69758923823346,-0.0869113988721213,1.77730853240763) q[2];
u3(1.83293696158456,-2.80941013592088,-1.53940692539428) q[0];
cx q[0],q[2];
u1(2.36130659430923) q[2];
u3(-1.99488295373826,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.102824963080876,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.72084156081397,0.807747012735710,-3.76330399092407) q[2];
u3(1.33831802174555,0.604389085097426,0.0883515787978009) q[0];
u3(2.96919962887663,-3.67568304082906,2.25137508626469) q[3];
u3(0.397213182097855,1.10973388903781,0.684179388298177) q[5];
cx q[5],q[3];
u1(0.162588301486843) q[3];
u3(-0.781461227867289,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.61733181521340,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.800839176657882,3.71204551696770,-2.38821362989253) q[3];
u3(1.47224394325445,-1.23097264900355,2.18009520922896) q[5];
u3(1.54989016544741,1.58829482304338,-0.361244732619973) q[0];
u3(2.00010935934890,0.403615958175640,-2.12154757681604) q[2];
cx q[2],q[0];
u1(0.585104586407062) q[0];
u3(-1.16286241822804,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.14659874578350,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.76042325140903,-1.47252641597916,1.09535144138044) q[0];
u3(2.33022620583015,-1.97845819418447,1.60050167480060) q[2];
u3(1.18350236846507,1.21515402261567,0.762299203154441) q[4];
u3(0.816373054967539,-0.496504253990343,-2.92224326543869) q[1];
cx q[1],q[4];
u1(0.747085631444343) q[4];
u3(-1.23885118842163,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.93024772657911,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.01649220707536,1.78974459538816,-0.801406487542267) q[4];
u3(2.20120980414096,3.65534402867430,-1.42131120677280) q[1];
u3(1.22871017753931,-0.132078472479101,2.11963717115624) q[0];
u3(2.02413046747185,-2.74568774040928,-1.06433662797644) q[2];
cx q[2],q[0];
u1(1.48003506143010) q[0];
u3(-0.667223542793031,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.00296655891425,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.42156805280166,-4.64963070003972,1.54899499227977) q[0];
u3(0.662126842787468,-1.88086643285673,1.35620840620424) q[2];
u3(1.14187406805521,1.17005623800283,-2.02344893908093) q[1];
u3(1.32041332514166,-4.31280299068870,1.85882719722700) q[4];
cx q[4],q[1];
u1(3.39502220521836) q[1];
u3(-1.44664971698057,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.55697177474287,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.74884866902967,0.270294965743651,-1.86963462691419) q[1];
u3(1.42182914671038,-1.82578993543299,4.25569141807334) q[4];
u3(1.00277433175657,1.40047977902393,1.06274813588677) q[3];
u3(1.21696455943741,-0.949949734412509,-3.03166819439447) q[5];
cx q[5],q[3];
u1(1.50504411532325) q[3];
u3(-0.516931991683764,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.01855120416455,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.57722871803615,-1.28629329514617,-0.544280620829822) q[3];
u3(2.53994093852478,-0.191328738024360,-5.20715166235407) q[5];
u3(1.97743096096746,1.98986674193432,-3.52257041429724) q[4];
u3(0.378764160068255,-0.105199230381882,1.10609462482499) q[1];
cx q[1],q[4];
u1(1.15814579570741) q[4];
u3(-0.881313145289317,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.05470970991544,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.23899602315322,0.335277684006534,-0.312589915715393) q[4];
u3(0.833601932131703,3.52747896025938,0.00673649939178178) q[1];
u3(1.78030270491679,1.06974469869213,0.177163994743888) q[0];
u3(1.50561411741047,-0.587587517867060,-2.71395506063983) q[5];
cx q[5],q[0];
u1(0.996057825957110) q[0];
u3(-1.35607201318340,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.77235202653908,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.06415313413995,-1.87923524252379,1.54739345480214) q[0];
u3(1.26093121540989,-1.18794652536969,3.90563180740666) q[5];
u3(0.470934574472967,-1.55106857779600,1.31755909686561) q[3];
u3(0.508081623867375,-2.95136960954993,0.979192742553532) q[2];
cx q[2],q[3];
u1(0.176124645437037) q[3];
u3(-1.74945808276124,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.19263929488586,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.944325949166124,0.198323724617081,-2.73992082127343) q[3];
u3(1.07055777628594,-0.531128777302777,0.0694895316502959) q[2];
u3(2.26490252828519,0.380567283454265,-2.61587051833239) q[1];
u3(1.75679437763750,-3.72440653321730,1.98695884910013) q[4];
cx q[4],q[1];
u1(1.39815367770040) q[1];
u3(-0.104153971943421,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.29338275007337,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.39345579883430,-3.34002585830203,0.444170342057856) q[1];
u3(1.12060197059266,1.82020375624320,1.31229985299333) q[4];
u3(1.76984829532019,1.47631612158858,-3.08364885430662) q[2];
u3(1.61434233333658,2.29363246263185,-3.21738183195057) q[0];
cx q[0],q[2];
u1(0.606492137756283) q[2];
u3(-1.80033198984183,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.330889846048331,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.04853541862383,1.22786772853203,0.504791211102303) q[2];
u3(2.27871893266333,3.05775007088573,0.0348511515641685) q[0];
u3(2.81629740711162,-1.02732863006649,2.21532294245866) q[3];
u3(2.55741327421199,-2.65931626293396,0.0304759493167732) q[5];
cx q[5],q[3];
u1(2.30480680628118) q[3];
u3(0.127244610923277,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.31982487540423,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.74305925366594,2.16551083090863,-0.906844453675344) q[3];
u3(1.61985075334674,-0.435120625693101,1.76413987005735) q[5];
u3(0.872920199896880,-0.605440356338425,-1.03474689758362) q[3];
u3(2.16658816192778,0.862917117853670,-4.36575825906347) q[0];
cx q[0],q[3];
u1(0.887370737201740) q[3];
u3(-0.0682843909448199,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.61297893048615,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.03206234647157,-1.86424575479239,3.89109855295126) q[3];
u3(0.802339900702435,-0.268413132326391,-4.93704176415958) q[0];
u3(1.47296553107931,1.37377172267948,-4.30955417474015) q[4];
u3(2.32344841169334,2.46947949183727,-2.96935671399223) q[2];
cx q[2],q[4];
u1(1.84949221087108) q[4];
u3(-3.00293257927792,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.462772061336895,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.923913557677493,-0.521366835551268,-1.03172985161469) q[4];
u3(0.555345590000804,2.08369587827146,1.21128124748362) q[2];
u3(0.849928855278469,-1.33928701426443,1.11817231426648) q[5];
u3(0.166731423067423,0.588836450176509,-1.52462371015143) q[1];
cx q[1],q[5];
u1(2.91101024762790) q[5];
u3(-2.28226499501513,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.44834278063975,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.16278497052895,1.82958792660583,0.100262018382652) q[5];
u3(1.62544640495347,1.63231982642208,0.260088023047515) q[1];
u3(0.934216678650263,3.69882380697880,-2.22695816896545) q[2];
u3(0.949092649771882,1.68041079567127,-2.04319887059204) q[3];
cx q[3],q[2];
u1(3.72694236043616) q[2];
u3(-4.43547880028119,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.329169856525205,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.636565626043978,2.08182536849053,-2.52433591611183) q[2];
u3(2.21906426544240,-0.519690787856289,-1.31128416801634) q[3];
u3(2.12992185350886,1.13680954912678,-2.93991287368184) q[4];
u3(2.79487737458723,4.23903376809419,-0.640708576565187) q[0];
cx q[0],q[4];
u1(0.688380513762491) q[4];
u3(-0.166495056838003,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.93007398985974,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.94493015767217,0.471621848080308,3.21002893796210) q[4];
u3(1.92964788574945,-2.22764590719180,1.74169243721101) q[0];
u3(2.15220178520072,-3.22750009586960,2.66241187233578) q[1];
u3(0.899089623963863,3.47076340355644,-2.52993901450737) q[5];
cx q[5],q[1];
u1(0.0250576102121123) q[1];
u3(-1.27705902367600,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.58436354861977,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.73268917587793,2.81600272055097,-0.126017727906930) q[1];
u3(1.32472925933618,2.07351414710968,0.617234344994198) q[5];
u3(0.663881566384891,1.27538412639622,-1.14123370380584) q[2];
u3(0.759346411235749,1.53660770311808,-1.80393497303619) q[1];
cx q[1],q[2];
u1(1.99901127596177) q[2];
u3(-3.15024373783289,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.39921038232066,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.18177200372788,-0.692662980878280,-1.31874513921488) q[2];
u3(0.447866608930489,1.82822723029437,0.696141165631621) q[1];
u3(0.640977066138929,0.996088765223977,0.264061474306918) q[4];
u3(0.950161586849845,-0.208035283250999,-2.36351933098586) q[5];
cx q[5],q[4];
u1(1.54164200282398) q[4];
u3(-2.98882910173923,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.874970964532062,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.48187278293465,0.119891729398120,0.219874708637567) q[4];
u3(2.26609406944107,-2.66545953940319,1.79258645945515) q[5];
u3(1.55245080341060,-0.194198920669551,-2.08721213477520) q[3];
u3(1.72073586538503,-3.88422113124624,2.29291128618002) q[0];
cx q[0],q[3];
u1(3.61917529983074) q[3];
u3(-1.58195598757185,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.19468902757092,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.394572766486326,1.46321107402384,-0.681163354273650) q[3];
u3(0.748360486703488,-2.47677749669926,0.988088297845143) q[0];
u3(2.83281433436880,3.45041622477852,-1.79632228372259) q[1];
u3(0.516574980924098,-0.653207585023908,1.81759568987877) q[5];
cx q[5],q[1];
u1(1.34528083345151) q[1];
u3(-0.334717475081327,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.38622887213511,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.69386853049789,2.84163927486699,-3.36003321664975) q[1];
u3(2.17019197741683,0.997531612980118,3.55391886365449) q[5];
u3(2.39230226204569,-1.68678039611617,0.340751915267251) q[2];
u3(1.68168899522327,-2.45842170522672,-0.270130350774127) q[0];
cx q[0],q[2];
u1(1.76609518992397) q[2];
u3(-3.04586408330107,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.28416129195332,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.68651893259322,1.38435344240395,-2.31438851745040) q[2];
u3(2.36254741696113,-2.97452817373544,0.0667481177072442) q[0];
u3(0.810206647041321,-4.29892748443412,1.21935435091924) q[3];
u3(1.27589796352437,-0.0390456341172229,3.12593427941655) q[4];
cx q[4],q[3];
u1(0.768086993272656) q[3];
u3(-0.330888677385132,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.62767901009560,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.36944662771835,0.692507491461939,1.54709792310616) q[3];
u3(0.516770412892697,-2.20747113284440,0.601915392980738) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
