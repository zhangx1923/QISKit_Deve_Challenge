OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.47033823317139,-0.635959761997133,-0.527634670259708) q[8];
u3(0.687983306572764,-0.155000588074165,-5.96852588265837) q[7];
cx q[7],q[8];
u1(0.806506877780861) q[8];
u3(-0.223739526357586,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.75680603120205,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.74569147598843,-1.36907441554668,1.32995825898957) q[8];
u3(0.353216694910837,5.32506135187238,0.507105620841303) q[7];
u3(2.33729427180304,-1.07317496513026,3.22492695952503) q[10];
u3(2.85134505346622,-3.30667431243475,-1.36834244111554) q[6];
cx q[6],q[10];
u1(0.447464723387449) q[10];
u3(-0.146897955597183,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.62885488270557,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.68637585353538,1.57215222501553,0.654481405873104) q[10];
u3(0.737658836514008,2.78739650178927,1.13670218203607) q[6];
u3(1.22253370842294,-0.736058357976083,0.914459650152153) q[3];
u3(1.56468372297696,-1.60417224747902,-1.77892598848450) q[2];
cx q[2],q[3];
u1(0.699010644805694) q[3];
u3(-0.227173774967096,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.17998104181913,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.58600282941428,1.95073747518133,-3.30320299436698) q[3];
u3(0.666730610956360,-1.03347639801666,0.00987240016444813) q[2];
u3(1.92816823874168,0.714663562109924,-2.58437665141419) q[1];
u3(2.91846721616883,3.44207784495935,-2.46174106782168) q[5];
cx q[5],q[1];
u1(0.767415565682811) q[1];
u3(-1.15515424262287,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.07136961622107,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.831714703895040,-3.55403328568035,1.87105811506556) q[1];
u3(1.54239565953710,-3.61484564475039,1.91712965133960) q[5];
u3(0.855088972662298,1.13650736032826,-2.12182403524810) q[4];
u3(2.02311127183208,2.03541895810519,-3.89464779044957) q[9];
cx q[9],q[4];
u1(-0.192275553466388) q[4];
u3(-1.69059843603262,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.08557035413019,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.49531512518889,0.0470289627663818,1.25751730661478) q[4];
u3(1.86492609866662,1.36338292229651,3.88031482004966) q[9];
u3(1.09229923825991,2.17836396282386,-3.15031412421796) q[4];
u3(1.52538077830430,2.97574188282006,-3.07005954563037) q[7];
cx q[7],q[4];
u1(2.29483684219160) q[4];
u3(-1.72108308268855,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.349251890531345,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.51501468839824,-0.893984783072334,-1.70444018385474) q[4];
u3(1.84250908882365,-1.60399362598696,2.58688543190951) q[7];
u3(1.78445143127372,0.146603236846794,-2.04168071527018) q[8];
u3(1.21820838155646,0.748174758823955,-3.53364350889726) q[5];
cx q[5],q[8];
u1(0.472310107579827) q[8];
u3(-1.32374583222614,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.0963879245449937,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.667268491381732,2.21786803154345,1.68144270422400) q[8];
u3(2.61458881756074,0.683397613307810,4.53106595800666) q[5];
u3(1.71314533678886,1.89207301269584,-2.49938621959118) q[9];
u3(1.05744016418166,2.17762249315355,-3.13461986326429) q[3];
cx q[3],q[9];
u1(1.21735451818102) q[9];
u3(-2.90467643152448,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.35126778683446,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.12370505647377,2.86491774640972,-0.970009246236537) q[9];
u3(1.70969780445407,0.767650522238968,-3.20267773711510) q[3];
u3(1.61271943804929,-1.93504680156303,0.0197386854407805) q[10];
u3(1.57249091472245,-1.96790350362429,0.725393802427227) q[1];
cx q[1],q[10];
u1(1.65619613068329) q[10];
u3(0.0870732612866285,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.12980519945574,0.0,0.0) q[1];
cx q[1],q[10];
u3(2.22602985029971,2.09488883270708,0.0824733484385918) q[10];
u3(1.43440742677039,4.99064448689217,-0.568431632347118) q[1];
u3(1.81318590586840,-1.39329455571005,0.424595160541532) q[2];
u3(0.818619860934810,-2.83084927081958,0.422522982637692) q[6];
cx q[6],q[2];
u1(-0.279922658916310) q[2];
u3(-1.57042668152716,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.05167861201667,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.37966168876975,3.40766779994803,-2.62007862826379) q[2];
u3(1.12358792686053,-2.14726252783101,2.92254927270303) q[6];
u3(2.99918034742785,-1.25178237377158,-1.83703752683839) q[1];
u3(1.01168369214333,-0.237239100380893,-3.67069794785053) q[10];
cx q[10],q[1];
u1(1.76770698086667) q[1];
u3(-3.02700301934056,0.0,0.0) q[10];
cx q[1],q[10];
u3(0.430723961375094,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.831339278977636,1.84053390712844,-4.17079864573058) q[1];
u3(0.341111780746547,1.39113997512597,4.38754602566281) q[10];
u3(1.66389414171315,-2.65767650881867,-0.192162826339587) q[5];
u3(1.73507147317705,-3.34502562342296,-0.795434567676569) q[9];
cx q[9],q[5];
u1(-0.357045921257671) q[5];
u3(-1.77254193302345,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.08158910885340,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.94836167683448,-2.61750710418007,-0.553689196752292) q[5];
u3(1.97657719861237,1.53317436302271,1.06289985659854) q[9];
u3(1.86036056279176,-0.731846317441857,-0.821389504826259) q[4];
u3(1.65452870388109,-3.46546366802843,-0.401410729824652) q[6];
cx q[6],q[4];
u1(2.71848666174931) q[4];
u3(-1.94505022171508,0.0,0.0) q[6];
cx q[4],q[6];
u3(3.18335655158023,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.07144990533651,1.73196886397449,1.26453339332207) q[4];
u3(2.21439758231041,-2.12956258655508,-1.24261786330957) q[6];
u3(2.17568841679962,0.345916533487460,-1.31303930918775) q[8];
u3(2.12600742721184,-4.56823596293264,1.30996930942977) q[3];
cx q[3],q[8];
u1(2.30998700909466) q[8];
u3(-1.88151251528284,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.534667956493570,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.12028934265714,-2.88406751037013,1.67992686818889) q[8];
u3(1.94172914682951,0.990563694014830,-1.46050994862874) q[3];
u3(0.190466235127076,-0.0118722323583339,-1.20766147283577) q[0];
u3(0.669674864744371,-0.215418390220127,-1.07078454816107) q[7];
cx q[7],q[0];
u1(1.69420076500014) q[0];
u3(-2.52515076374158,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.166212629674923,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.630951668790490,1.70870351029404,1.44054326123376) q[0];
u3(2.51757112517276,-1.45537782788869,0.495889610918878) q[7];
u3(2.01744585428116,-0.820125236832403,0.367318830221363) q[9];
u3(1.31666399872084,-3.13201966316562,-0.905377423468499) q[8];
cx q[8],q[9];
u1(0.398954684748361) q[9];
u3(-1.38477744554048,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.25754602827275,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.11572804414452,0.471316637868532,-1.01604356871473) q[9];
u3(2.14967150793012,-2.52082559916945,-1.21897831748504) q[8];
u3(1.25660310654677,3.38998101372612,-0.380737310398376) q[2];
u3(2.04366795240875,2.88675716577549,1.07345264988710) q[0];
cx q[0],q[2];
u1(0.596415663842222) q[2];
u3(-0.271260053685874,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.44368840314883,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.11914457258242,0.620743656312295,-1.81033220856136) q[2];
u3(1.36990066591492,-1.64789569440770,-2.01086813734713) q[0];
u3(1.99024630540794,1.36848620947759,-2.37735037291208) q[5];
u3(2.59904497890319,0.294321504581910,-5.36603394903493) q[4];
cx q[4],q[5];
u1(0.827794447570198) q[5];
u3(-1.51320732582464,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.453407731763241,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.49150782276729,-0.683605657958226,-0.950045830310112) q[5];
u3(1.13888934163034,-5.27223476092131,0.911022424263703) q[4];
u3(0.997181535874500,3.04802282611732,-1.79398545476564) q[7];
u3(1.48213411347809,0.361725749808488,-3.05024374497981) q[6];
cx q[6],q[7];
u1(3.70966971513879) q[7];
u3(-4.25750743976249,0.0,0.0) q[6];
cx q[7],q[6];
u3(-0.172382414546203,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.122920402532496,2.66172696366886,-1.56705978043088) q[7];
u3(1.68752086821305,-0.886633794384079,0.453397750397291) q[6];
u3(2.36990667654152,-1.71450711901020,3.48790524965606) q[10];
u3(0.874187613079524,-0.820328972564512,2.28189540605024) q[3];
cx q[3],q[10];
u1(1.57722891720295) q[10];
u3(0.0396593115638741,0.0,0.0) q[3];
cx q[10],q[3];
u3(0.259527942486499,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.564549697755481,-2.87120277728764,0.972978421234208) q[10];
u3(1.77440157448934,-1.49216533037822,-2.87168643688630) q[3];
u3(0.648409977022571,-0.937065272667254,2.28439133726022) q[0];
u3(1.19514344538598,-2.48412093664905,-1.28120705362561) q[5];
cx q[5],q[0];
u1(-0.544955615074805) q[0];
u3(1.24826451379839,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.86786935555378,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.425454810140031,-1.39859553363171,-0.672296778366807) q[0];
u3(1.18098643977631,-0.555459753209702,4.91485045135449) q[5];
u3(2.15821065174638,1.64509230752346,0.557176088153337) q[10];
u3(2.10454241940033,0.675085928951460,-2.13425617073921) q[6];
cx q[6],q[10];
u1(1.81211558802079) q[10];
u3(-2.29079920940915,0.0,0.0) q[6];
cx q[10],q[6];
u3(3.03621564515814,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.07065908884690,1.07453612613305,-2.89171654566299) q[10];
u3(1.64280628562351,1.51822782786662,-3.61964495071480) q[6];
u3(1.99878475871774,-1.70866129138151,0.387122833477537) q[2];
u3(2.12315749946631,-2.24164785114243,-1.09024966728615) q[9];
cx q[9],q[2];
u1(2.45986250756348) q[2];
u3(-1.65496346808577,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.0718870039680410,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.80354449381356,0.806946389839996,-0.586104658726305) q[2];
u3(1.33062968625748,-1.11011222270444,0.994265640655396) q[9];
u3(1.66482596329310,0.708328427184971,2.35608252900628) q[7];
u3(1.30489904858869,2.30828228794236,3.46610603116262) q[8];
cx q[8],q[7];
u1(2.09904991992453) q[7];
u3(-1.55538618431310,0.0,0.0) q[8];
cx q[7],q[8];
u3(3.64218216802007,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.44657117609226,1.92419628407668,0.461503623955844) q[7];
u3(0.649332027099127,-4.64014154781048,-0.578434251086584) q[8];
u3(1.29898787928299,-1.03156362708240,1.86810778077183) q[1];
u3(0.136291126057921,-3.38159765053794,1.87678385584016) q[3];
cx q[3],q[1];
u1(1.88465021989959) q[1];
u3(-1.77847017241746,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.552459589711155,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.47682188357872,2.31665735500006,-0.421959135664150) q[1];
u3(0.343798036565197,-1.42514118544790,4.65279027812008) q[3];
u3(0.907020367460601,0.970870692255760,-1.27186918444184) q[1];
u3(1.40663571178904,-4.25067628005139,0.761078263402925) q[2];
cx q[2],q[1];
u1(2.14961640896785) q[1];
u3(-2.68124325596867,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.50486410032969,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.69349437452435,0.433928254308715,1.79581106166120) q[1];
u3(1.55161258924011,2.29345616719246,-3.94806469317049) q[2];
u3(2.29101626668770,1.27438074424488,-2.47017335313947) q[5];
u3(2.10795854521849,4.59366081636664,-0.558924719243610) q[3];
cx q[3],q[5];
u1(1.34518393142336) q[5];
u3(-3.66893126986050,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.32978146580358,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.830680270132653,-2.98396219735900,1.47764448825481) q[5];
u3(0.778740732512869,1.26299818097688,3.00721597259832) q[3];
u3(1.52720227713057,-1.27425178045624,0.954914317305542) q[6];
u3(0.901164749006955,-2.35304239018610,0.236168397141725) q[4];
cx q[4],q[6];
u1(1.00352857857061) q[6];
u3(-0.0308564867540506,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.08762918417855,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.709352263314024,-1.14452652954438,4.16464739790829) q[6];
u3(1.26069025162487,-0.784326480104832,-0.429812498275795) q[4];
u3(1.49914472655242,-0.501531759676570,-0.397395007403853) q[8];
u3(1.51845009313799,-3.10171941087221,0.886107846334304) q[7];
cx q[7],q[8];
u1(2.27847760657494) q[8];
u3(-1.75669272492957,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.442486721168913,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.337239684484849,1.92863281026966,-2.21585454454520) q[8];
u3(1.73782371568489,2.00677884715956,-3.58387336285117) q[7];
u3(2.85546261258195,2.97050928027369,-1.70578836844710) q[9];
u3(1.66030116044547,2.55658908785360,-2.89999674039437) q[0];
cx q[0],q[9];
u1(2.54065904652594) q[9];
u3(-1.60558531735382,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.99775377215300,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.26480917762728,2.59057900137748,-0.453263090926274) q[9];
u3(0.148673700880238,2.79320087850256,1.59013974284206) q[0];
u3(2.31304608972649,1.25083623451925,1.57881630990150) q[8];
u3(1.57540910771475,-1.82317986699612,-2.17875640506438) q[9];
cx q[9],q[8];
u1(1.42592925837445) q[8];
u3(-0.952283577668020,0.0,0.0) q[9];
cx q[8],q[9];
u3(3.11475805217423,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.23119531545067,2.16600317281222,-1.52277780364685) q[8];
u3(1.57539267973609,2.08205086201283,2.97862789328057) q[9];
u3(2.58606887358257,2.71595084969677,-3.17462201071302) q[7];
u3(1.33595710098516,1.64207928756695,-1.36469334839398) q[1];
cx q[1],q[7];
u1(1.10728170182702) q[7];
u3(-0.766953661643522,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.200757457825096,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.29949890543295,-3.26637876444078,1.51975025861667) q[7];
u3(2.38353088249102,3.21813608676498,-0.595829162537074) q[1];
u3(1.14478105377190,2.91474012673635,-0.553240190185071) q[5];
u3(0.798458244453610,1.68810541592675,-1.13606948665634) q[3];
cx q[3],q[5];
u1(3.00056753040366) q[5];
u3(-2.43325034407318,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.30992706752112,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.06889829556960,-2.31574625573281,2.39828790775207) q[5];
u3(1.96849421911197,-2.72101494378470,-2.12192459474548) q[3];
u3(0.462221779040213,-0.119437128290452,-1.88957992540163) q[2];
u3(1.97804889689809,0.874571628791093,-4.32881630236766) q[10];
cx q[10],q[2];
u1(-0.0916340058119192) q[2];
u3(-1.68359558319746,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.644381574974226,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.56743943357257,-1.12628728247628,0.641980874865184) q[2];
u3(0.455247225581983,-2.91014286620999,-0.200252426425421) q[10];
u3(2.03259619129584,-0.438252432718767,-0.433788478476232) q[0];
u3(2.04545074073731,-2.88474667171000,0.678460525751534) q[6];
cx q[6],q[0];
u1(1.47618075451547) q[0];
u3(-0.915292184926405,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.21523191098098,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.16186699245988,1.99523538326711,-0.963666528341390) q[0];
u3(2.46752222543852,2.58436187127710,-3.38616324696209) q[6];
u3(1.93836885761347,3.22372238109773,-0.319338048661798) q[2];
u3(2.31810370402128,2.97671473881373,-0.439309462089950) q[7];
cx q[7],q[2];
u1(0.369089046257949) q[2];
u3(-1.03038086782921,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.85156603016034,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.51210440499073,1.53060160946422,0.208585902852660) q[2];
u3(1.53776731649256,-1.64229143663061,-2.67006884782413) q[7];
u3(0.973757577661533,0.628444546465855,-1.09735785833041) q[3];
u3(0.798565253157464,-1.24208054011503,-0.587286963664502) q[0];
cx q[0],q[3];
u1(1.68134888409988) q[3];
u3(0.595356599772709,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.758452859551831,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.15751497901181,2.12195937748704,-1.03007576431578) q[3];
u3(2.41222880559704,2.22198281696521,-1.62593927571517) q[0];
u3(0.966897442773463,-1.01276772681689,1.58218789665625) q[1];
u3(1.08327576612445,-0.842331420900790,-1.70800628311564) q[5];
cx q[5],q[1];
u1(0.848380434942247) q[1];
u3(-1.40888117009464,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.255076559039155,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.83109108059458,-0.354991125923333,2.35427543811212) q[1];
u3(2.94921719676397,0.465642127218793,-4.42976824628898) q[5];
u3(1.32601896282280,-1.09035177217485,2.03411925897024) q[6];
u3(1.68144551115240,-1.58272887567282,-2.73106446500764) q[8];
cx q[8],q[6];
u1(2.43663769832545) q[6];
u3(-2.73626861933548,0.0,0.0) q[8];
cx q[6],q[8];
u3(-1.25321438956021,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.35596080912388,1.14282557079433,-3.43766246457759) q[6];
u3(0.538153825450983,4.22541544298250,0.874033361385455) q[8];
u3(1.44654629965371,-1.35042061931662,1.03243539047768) q[10];
u3(1.74914807014447,-2.24560514263617,0.238689391381509) q[9];
cx q[9],q[10];
u1(1.02903667320943) q[10];
u3(-0.261156563194263,0.0,0.0) q[9];
cx q[10],q[9];
u3(2.33291900752200,0.0,0.0) q[9];
cx q[9],q[10];
u3(2.10836303061828,-0.191472814672206,3.12375746636262) q[10];
u3(2.14903765640291,-0.00344792297418195,5.31914225903168) q[9];
u3(1.89941234467043,0.0996246773793482,1.57824133276592) q[9];
u3(2.01514379253460,-0.755350755646635,-2.51039026388839) q[10];
cx q[10],q[9];
u1(1.34879907363042) q[9];
u3(-0.122813124316542,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.08573463605223,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.81514641138984,3.27767089324099,-1.21071345028537) q[9];
u3(1.82197494408387,-3.52930355907632,-0.242883680420564) q[10];
u3(2.27523544266029,0.847913814141665,1.14912633738094) q[7];
u3(1.56622157162781,-1.69455882102831,-2.01547971621885) q[0];
cx q[0],q[7];
u1(1.28590908453814) q[7];
u3(-0.512089041945275,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.165954057852290,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.90978454461997,-2.01940523667412,2.75783944906341) q[7];
u3(0.786229978750413,3.29287473332236,2.63346162920165) q[0];
u3(2.23077908049243,-1.53204642282388,-0.842346063249857) q[1];
u3(0.364185950080215,-5.20120383182941,-0.0494891698570559) q[3];
cx q[3],q[1];
u1(1.92450800435819) q[1];
u3(-2.25497935625701,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.33722349123759,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.77876848758385,-2.99325996090669,0.0486195795641196) q[1];
u3(1.16970406856205,-0.608373202645285,-1.48981337912361) q[3];
u3(0.385525495540430,-1.52858846359075,-0.761086647769335) q[4];
u3(1.09556564974660,-2.75816486999545,-0.831777782682178) q[2];
cx q[2],q[4];
u1(1.71890846428338) q[4];
u3(0.409032374415819,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.592745199204442,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.01088844090319,-0.147857190307146,-0.949564811892801) q[4];
u3(0.839938710113064,-2.83819549508955,-1.27617374122326) q[2];
u3(2.24827230028186,0.799260572677603,-0.845649933682298) q[5];
u3(1.87457475849236,-5.00718035068820,1.18943111493358) q[8];
cx q[8],q[5];
u1(0.744741697898627) q[5];
u3(-0.555206412208904,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.53106196607544,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.814128825212851,1.98645048218245,1.19815838427512) q[5];
u3(1.64158204486997,-1.99599605706321,1.66454015930878) q[8];
u3(1.56173502609383,0.595949141745052,-3.47825832051785) q[1];
u3(2.74887785644331,3.53130212598417,-1.78608151230471) q[4];
cx q[4],q[1];
u1(0.390083892372170) q[1];
u3(-0.169071684536730,0.0,0.0) q[4];
cx q[1],q[4];
u3(4.29260621841350,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.367640479738960,2.84677383853554,-2.50629778736816) q[1];
u3(1.15745383073281,3.08395322122626,0.749040437560333) q[4];
u3(0.513254120782601,1.31473388643702,-2.64349064465359) q[2];
u3(1.59802601673213,-3.39299314671578,2.76761374970597) q[10];
cx q[10],q[2];
u1(0.478260253008439) q[2];
u3(-1.42366644419088,0.0,0.0) q[10];
cx q[2],q[10];
u3(2.55910928443499,0.0,0.0) q[10];
cx q[10],q[2];
u3(0.890620935858210,-4.44335980277224,0.541504616346628) q[2];
u3(1.16804046359774,0.284048930243201,-3.60161283511140) q[10];
u3(1.29571050003270,0.948434672435227,-2.80408778957357) q[9];
u3(0.913932470334436,-3.04055956525716,2.56444527486885) q[5];
cx q[5],q[9];
u1(0.316489236302470) q[9];
u3(-1.68795690142384,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.32254362435779,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.52304238002948,1.25182929997034,-2.61090896165567) q[9];
u3(2.38877929883782,-1.89572267105544,-1.87408640814196) q[5];
u3(1.50600022736104,-1.52497199899477,-0.514212393478509) q[6];
u3(0.562287585646602,-3.45057168305140,0.236937708350022) q[0];
cx q[0],q[6];
u1(1.46557476511238) q[6];
u3(-0.729162305305857,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.99428221421305,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.81206092678971,0.504155512518868,2.66949863566651) q[6];
u3(0.291263461805180,1.28158001543038,-3.57961241523725) q[0];
u3(2.55619951269404,3.04536222664939,-0.991082113115363) q[7];
u3(1.96416201205815,2.03756352046297,-1.30294887082653) q[3];
cx q[3],q[7];
u1(1.26080941867811) q[7];
u3(-3.39740213712020,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.23454843028850,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.137938351898049,1.84894175960006,-2.07792326275402) q[7];
u3(1.16805415325998,-2.04541433774299,3.54748279652732) q[3];
u3(1.60089202199677,0.626956636293327,0.229911305248185) q[9];
u3(0.873746327189451,-2.24017470928145,-1.68482085272408) q[0];
cx q[0],q[9];
u1(1.39734966967363) q[9];
u3(-0.843018810158679,0.0,0.0) q[0];
cx q[9],q[0];
u3(-0.238339091630925,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.13455229694461,2.48036236583839,-2.75420731749195) q[9];
u3(1.26914562224074,3.42572444242375,-0.706621032666279) q[0];
u3(0.821899714639593,2.73593160019775,-0.609737963820598) q[1];
u3(1.47158270111637,-0.135828080663345,-3.68048371368515) q[5];
cx q[5],q[1];
u1(2.99336063317623) q[1];
u3(-1.82214595145517,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.14586823243865,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.720170380206743,-1.53010403220259,-2.15752362796884) q[1];
u3(2.89564105982197,5.70805158656218,-0.473839117266680) q[5];
u3(0.629215891283279,0.000714815150187082,-0.0349056503917710) q[6];
u3(1.67497502269272,-0.692208701404169,-1.26344802026132) q[7];
cx q[7],q[6];
u1(-0.263778585958732) q[6];
u3(-1.25905080679328,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.27970010543192,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.07020962672840,-4.26881298467839,1.11328386610441) q[6];
u3(0.388815546272175,-2.34435910173466,0.746150440205505) q[7];
u3(1.80263994077148,-1.98228519045468,0.962549300114715) q[10];
u3(2.28451001016172,-2.55740232123533,0.362628729719476) q[3];
cx q[3],q[10];
u1(1.35548951097133) q[10];
u3(-0.0957583983052590,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.60992527928471,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.0920345872659381,2.08188120285979,-0.952098784467339) q[10];
u3(0.861728495019577,-1.20366412127247,-2.91257626532052) q[3];
u3(2.73407099672535,2.71309958888110,-3.42346592875494) q[8];
u3(1.18273121735228,-0.501791568365862,1.63661979385006) q[2];
cx q[2],q[8];
u1(1.48619026834308) q[8];
u3(-0.318365421517491,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.99108953931809,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.67441316452478,-1.66778846583965,0.850245898609355) q[8];
u3(1.65148192612445,0.402507292947764,-3.36037249443996) q[2];
u3(0.854565209819639,1.70851695988948,-4.29679075268559) q[10];
u3(1.04154917661416,2.40872648280415,-2.83614772586349) q[0];
cx q[0],q[10];
u1(1.25151181186675) q[10];
u3(0.0629620995964311,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.17842725640700,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.34943052581125,0.247609037994118,-2.54010767455927) q[10];
u3(2.21519186426990,2.36127489101104,2.59194285645980) q[0];
u3(2.12854565538445,-0.112308152811571,0.0297095563560737) q[9];
u3(0.946622213011287,-2.82122300574410,-1.81967934257227) q[1];
cx q[1],q[9];
u1(4.04990369127791) q[9];
u3(-4.16623835432957,0.0,0.0) q[1];
cx q[9],q[1];
u3(-1.26622412421555,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.09140593004745,-1.35068015812239,3.19118577935104) q[9];
u3(1.91300797420576,-2.24002726284243,-1.89697452321029) q[1];
u3(1.46381724316272,2.17379827791185,0.602241685187781) q[6];
u3(2.08317912314180,-0.280796947413190,-3.17812122006141) q[4];
cx q[4],q[6];
u1(-0.473443416364411) q[6];
u3(-1.82567054877511,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.04112733684346,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.46212654079267,-2.42219003304730,3.80954246632818) q[6];
u3(2.14788461961594,-2.95492192718953,3.02017615438744) q[4];
u3(0.458549643730333,2.67205353259942,-2.96210594596063) q[7];
u3(1.45945748508187,-2.52761167724445,2.08525279759207) q[3];
cx q[3],q[7];
u1(2.22491979536978) q[7];
u3(-2.73162501375725,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.791264098291807,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.06840393913607,2.17457622310088,-0.929599999107300) q[7];
u3(0.740902801951878,-3.10341904436798,1.97469924627095) q[3];
u3(1.63705078906109,0.339575689602359,-2.48430436942147) q[2];
u3(1.89937806269414,-2.91931009429355,2.47547574195750) q[8];
cx q[8],q[2];
u1(1.48370118453617) q[2];
u3(-3.45407995823164,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.59638164399602,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.27208915713822,-0.693554825534108,0.112051585354231) q[2];
u3(1.52621340317700,2.02742563532545,-0.463672349324592) q[8];
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
