OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(0.854464304272954,0.554728319883076,-2.02115923179766) q[3];
u3(1.71736657154265,1.05317611545002,-3.61125928874907) q[5];
cx q[5],q[3];
u1(1.30321926123036) q[3];
u3(-3.67800694062179,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.07121570269573,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.995797247393521,3.58498807870814,-2.39597348465161) q[3];
u3(1.66952136779617,1.91550004669618,-3.10224027790265) q[5];
u3(1.40918236950399,2.16368895761367,-3.51430980075551) q[2];
u3(2.19069867007323,2.85516348562707,-3.02653908759542) q[4];
cx q[4],q[2];
u1(2.01517590309040) q[2];
u3(-3.14957773933510,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.293794102067143,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.97246175805341,-0.927873167425869,2.00280465050960) q[2];
u3(2.32206120408105,-1.43700399401129,4.75814834695414) q[4];
u3(2.34580040867820,0.510679863643496,-3.51451123321757) q[6];
u3(2.85294994554604,2.10618202791821,-2.95263809922819) q[1];
cx q[1],q[6];
u1(3.23245609169185) q[6];
u3(-0.671864366189204,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.83483586266484,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.07888811947606,-2.29628167177107,0.126499707840754) q[6];
u3(1.21610139310280,1.26321123775281,-3.15142808524366) q[1];
u3(2.10366784067170,-0.0926961543277373,-1.74998902483197) q[7];
u3(0.820785686492266,0.678045955751412,-4.67182928584256) q[0];
cx q[0],q[7];
u1(0.292996288071379) q[7];
u3(-0.945185014681441,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.28205959916279,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.484073942127223,-3.23204530654604,2.81067778337284) q[7];
u3(2.59456397774737,1.63248231507317,-4.45607735620547) q[0];
u3(2.13771955968725,3.19642861276408,-0.985811003835603) q[5];
u3(1.00145717480591,0.491108837536709,-0.909225539469801) q[1];
cx q[1],q[5];
u1(-0.369565476660615) q[5];
u3(0.725495610596985,0.0,0.0) q[1];
cx q[5],q[1];
u3(4.17166934848569,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.72375627110950,-2.11475923820330,-0.308441726182503) q[5];
u3(1.31846245559119,0.220963047337407,1.23720699901420) q[1];
u3(1.63670273132667,0.370561645642623,0.988863412671925) q[7];
u3(0.981028808468790,-2.23314895179744,-1.92688252169226) q[3];
cx q[3],q[7];
u1(1.27165862458238) q[7];
u3(-0.870877739627912,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.0699104984813923,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.110145598457917,2.86608190165235,-3.13463173554709) q[7];
u3(0.303607574817482,1.79518335331837,-0.121203882504243) q[3];
u3(1.96228148677217,3.10495697983837,-2.53994225422233) q[6];
u3(0.982475511914167,3.21058473704053,-2.55065875817641) q[2];
cx q[2],q[6];
u1(2.84153223920992) q[6];
u3(-2.24162617996455,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.64750967692049,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.01270845949627,-4.20209907771447,0.783650884278982) q[6];
u3(1.43525893236195,-0.595981158433323,-4.05057132426453) q[2];
u3(1.63173924266479,-1.16594492348375,-0.429798615687504) q[4];
u3(1.63774669150184,-2.71201467505028,-0.981122722990252) q[0];
cx q[0],q[4];
u1(-0.183077417537024) q[4];
u3(0.484721885465786,0.0,0.0) q[0];
cx q[4],q[0];
u3(4.25103085628072,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.58059629280158,-0.941429646652846,2.62558080721552) q[4];
u3(0.0816930686793076,-2.45136070439197,-2.31247391731390) q[0];
u3(2.25184280991946,-0.202466081052776,1.91000162915558) q[5];
u3(1.72211284460522,-2.48571227615218,-2.65446314462957) q[2];
cx q[2],q[5];
u1(3.33190048639630) q[5];
u3(-1.33147816671182,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.487957333515232,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.88167668473355,-0.154644446123712,-0.0264567535498365) q[5];
u3(0.715529279013968,0.841522163954729,-1.47403405037459) q[2];
u3(0.708749891877725,-1.71586587352826,4.47952826094804) q[7];
u3(1.34545284875096,-0.767871425106759,0.548093478932838) q[3];
cx q[3],q[7];
u1(0.633813230859707) q[7];
u3(-0.0868726629997398,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.99696506913175,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.25868178115741,-0.926447806659182,1.74178982790181) q[7];
u3(0.886144605957333,-3.23174268727731,2.22372371541566) q[3];
u3(1.06988675419861,0.0990724968685867,1.32294721694061) q[4];
u3(1.53722315240148,-0.181255128654176,-3.07681259481118) q[0];
cx q[0],q[4];
u1(4.25375545210265) q[4];
u3(-3.44491987982668,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.543294914324811,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.30774381801465,-0.590215431138496,-2.68718349442025) q[4];
u3(1.02234036312115,1.37411591537299,-3.54275402143286) q[0];
u3(1.42676299655644,-1.78491517065305,1.34584887253458) q[6];
u3(0.137654668227735,1.68828467555481,-3.43327162704974) q[1];
cx q[1],q[6];
u1(1.58036316600598) q[6];
u3(-0.543981150768596,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.14783901295853,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.61897258007825,-0.145430244500459,-0.684270309037323) q[6];
u3(1.81152958885339,-1.17539955702969,2.07709949649248) q[1];
u3(2.01608259187985,-0.155896842900858,1.95683446766049) q[6];
u3(2.02331565041503,-0.586023288315679,-1.57519789845959) q[7];
cx q[7],q[6];
u1(-0.437975099931947) q[6];
u3(-1.51922196735064,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.10007968639900,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.76053489641052,-1.25396642129705,1.50129609784729) q[6];
u3(1.22428875564055,-1.02260283739372,-4.26365258503699) q[7];
u3(2.02650798994424,-0.401259310041891,2.27778599234763) q[3];
u3(1.80705054923409,-2.29798766727814,-1.26587014174549) q[5];
cx q[5],q[3];
u1(1.26063710472840) q[3];
u3(-0.635479562911414,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.10769345926364,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.698558179352900,-1.50199475322129,0.909250854779655) q[3];
u3(1.87997521792104,-3.07819851617375,-1.35475029917631) q[5];
u3(2.39407596661800,1.50414743932804,-2.78525635771820) q[1];
u3(2.61432780736108,5.24385463514606,0.458968416832500) q[2];
cx q[2],q[1];
u1(0.160087822884286) q[1];
u3(-1.52495072888132,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.04748606914406,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.790846842832117,0.204464821072612,-0.876730390627333) q[1];
u3(1.86402528229807,-3.14592685288502,1.79511761481310) q[2];
u3(2.32078418632137,-1.38256925035291,2.14052111992258) q[0];
u3(2.46152156532445,-2.62667552301458,-0.565552255005418) q[4];
cx q[4],q[0];
u1(0.902446649823391) q[0];
u3(-1.38087338643450,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.376358478669108,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.508607091040189,-1.35159539068322,2.87394190740399) q[0];
u3(0.723057600779228,-4.12332192171684,-1.40217686450800) q[4];
u3(2.72721233842650,0.229774180818210,-1.43984328706674) q[6];
u3(1.29781900112375,1.07028391236602,-4.59854667110129) q[4];
cx q[4],q[6];
u1(4.32762713328491) q[6];
u3(-3.53834878637670,0.0,0.0) q[4];
cx q[6],q[4];
u3(-0.624028123865806,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.61435108619405,-0.865321936923662,0.805984400856175) q[6];
u3(1.77636896004956,-2.20845450596208,-2.57393222975677) q[4];
u3(0.997300625395125,-0.149255503870938,2.88315132654154) q[7];
u3(1.29753665547595,-1.65434311595164,-1.44280103802907) q[1];
cx q[1],q[7];
u1(0.757787747932027) q[7];
u3(-1.51808970981327,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.0118752301119625,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.683648271177689,-4.07176085637306,0.702901166279085) q[7];
u3(1.84017728102915,5.86576149730952,-0.00198280272176321) q[1];
u3(0.811932498882882,0.758537269421753,-2.07243857049125) q[0];
u3(1.06903318440312,1.33371760388643,-4.57648674116637) q[5];
cx q[5],q[0];
u1(-0.207318309478508) q[0];
u3(-1.92217631439868,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.619935583243174,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.99554899344462,2.64346019289519,-2.71197992807985) q[0];
u3(0.586641910847245,2.84592003782877,2.99100087768582) q[5];
u3(0.478087874439396,-0.309590390961801,1.38714478697919) q[3];
u3(1.05056907977712,-2.05085072128270,0.417511835901456) q[2];
cx q[2],q[3];
u1(4.01021968953527) q[3];
u3(-3.37065536523843,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.692636352999680,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.52022451313922,0.624676475574605,-0.699025282980044) q[3];
u3(1.39975268342712,-5.07093843026110,-1.04010895641442) q[2];
u3(1.44510765298195,0.647125994488176,2.18014745727137) q[4];
u3(2.24011972139941,-3.56699064894699,-2.40517603092569) q[5];
cx q[5],q[4];
u1(2.57064100191671) q[4];
u3(-1.96099129215184,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.790949343025240,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.36990481105982,1.17710248755203,-3.33577495462340) q[4];
u3(0.979464835524744,2.11542776225341,-2.82049185252605) q[5];
u3(1.77801507732548,-2.41620917630042,0.602707660149289) q[7];
u3(1.14389041262547,-4.05212045421939,0.592135052037850) q[6];
cx q[6],q[7];
u1(2.21424404624813) q[7];
u3(-2.41169677124788,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.473928426904284,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.81091406870086,2.19758683072848,0.351572420904230) q[7];
u3(2.32184177987963,1.20656011237796,0.439774157297993) q[6];
u3(2.44653256068766,1.84602904039098,-2.34812202176311) q[0];
u3(1.95485218759922,2.62525216460689,-3.30264108593487) q[1];
cx q[1],q[0];
u1(0.325547921679259) q[0];
u3(-1.00997693715304,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.29442333582126,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.34979627710495,-1.11022254501476,-0.0788854383660026) q[0];
u3(2.04623882823865,2.72794547414233,-3.34459729686612) q[1];
u3(1.96240560941465,-0.690612008024240,0.259650439154110) q[3];
u3(1.22951054588182,-3.17042587977881,-1.31119327456553) q[2];
cx q[2],q[3];
u1(1.45234272142760) q[3];
u3(-1.08296884644665,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.545325352698958,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.54672537357949,-3.36438311039570,0.917841294885484) q[3];
u3(2.16455413906245,-3.32091631734853,-1.14206540022781) q[2];
u3(1.41106900085890,0.847318806608841,-3.25769087228007) q[3];
u3(1.89031161407022,3.60718747583459,-2.44382103480836) q[0];
cx q[0],q[3];
u1(1.63217333009082) q[3];
u3(-1.02007011954854,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.20055649268640,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.742377899147116,-4.56859560836884,0.718546937543772) q[3];
u3(0.936595459194771,-4.06746910421578,-1.47036159058185) q[0];
u3(0.699517575272441,-0.853759177675714,1.00159331454417) q[4];
u3(0.752809360988627,-1.31174810467669,-1.25323036604550) q[5];
cx q[5],q[4];
u1(0.857205266147023) q[4];
u3(-1.71502738351605,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.15821845343798,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.770297928080983,-0.825080151033578,2.85244109276940) q[4];
u3(1.44807607632978,-0.667954748627297,-0.629777041750803) q[5];
u3(1.54737570824028,0.668146789104051,2.17218993097245) q[1];
u3(2.32367913506001,-1.97820003845888,-1.13188379674634) q[7];
cx q[7],q[1];
u1(1.93034839069676) q[1];
u3(-2.56867777925981,0.0,0.0) q[7];
cx q[1],q[7];
u3(3.14761690487979,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.17892463053019,-2.02171335722284,-0.206298325077726) q[1];
u3(0.639026115471110,-5.04677051901916,-0.129887737315022) q[7];
u3(2.63741754087235,-0.977171754840889,2.08293061752154) q[2];
u3(2.01483172856670,-1.45907552046854,-0.792392603661351) q[6];
cx q[6],q[2];
u1(1.17358898008054) q[2];
u3(-0.121718828600082,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.49156269071647,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.567597064882374,1.82909787703771,-3.46205548192971) q[2];
u3(0.842446575876999,2.84110609118644,1.13612470806500) q[6];
u3(0.753246930367553,-1.24303153001940,1.71087181883608) q[5];
u3(1.02540167922448,2.29146697058501,-3.28533267716276) q[4];
cx q[4],q[5];
u1(1.61333605593628) q[5];
u3(0.497201757599740,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.791617529376914,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.675502034862386,1.31891625362089,-0.243194656600755) q[5];
u3(1.30067150383380,-2.39789501531326,-3.17689921438174) q[4];
u3(1.90781555250748,-0.128902339557611,-2.39531729387480) q[3];
u3(2.45400941246273,3.37793836443819,-0.678628508483444) q[0];
cx q[0],q[3];
u1(2.80208061962362) q[3];
u3(-1.88309620071538,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.04330746941096,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.77626855839913,1.80078000713860,-4.30143095992172) q[3];
u3(1.18058680170777,1.57870898685363,-0.00763161721693772) q[0];
u3(0.893198178742416,0.420558831325262,-1.91584166780958) q[2];
u3(1.63689220788358,2.16656103203877,-3.28273994593939) q[1];
cx q[1],q[2];
u1(1.45475893356409) q[2];
u3(-3.36163480432535,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.06681771455285,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.463825865301779,3.52419116413912,-1.77393786363341) q[2];
u3(1.18045104635626,-3.46476940976718,-1.85362995113941) q[1];
u3(1.00176861417465,1.02562433220339,-1.84624505668107) q[7];
u3(1.81384559759473,-4.29636845610679,1.69930680575053) q[6];
cx q[6],q[7];
u1(1.71084465255120) q[7];
u3(-2.97591647814447,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.719751622998264,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.37558057161250,-1.58714830939199,1.91016393492037) q[7];
u3(2.00248626746388,-0.0150615545031156,4.36820620290203) q[6];
u3(1.20988854577803,1.31770440840399,-0.262703052533601) q[7];
u3(2.04413388683052,-0.520857654323115,-4.16114141778765) q[1];
cx q[1],q[7];
u1(2.36718093567840) q[7];
u3(-2.92821587918244,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.04130900048968,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.05047904351015,1.76689669709740,-3.66308701289875) q[7];
u3(0.974869207730757,-6.06125724343462,-0.0676369646563302) q[1];
u3(0.983291564925307,0.644411458685705,-1.63381499067533) q[5];
u3(0.512106600435861,1.95889346232202,-3.93546661171318) q[4];
cx q[4],q[5];
u1(0.926721620683431) q[5];
u3(0.0994972625858557,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.48873634702882,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.18479888510852,-1.37950676646144,-2.04096020850096) q[5];
u3(2.16397327872720,-0.545261902780628,-3.26736187755732) q[4];
u3(3.08690885401593,1.42959702034280,-1.32338675643277) q[0];
u3(2.58511247254069,3.93027964694425,-0.784524769491824) q[6];
cx q[6],q[0];
u1(2.23486678461540) q[0];
u3(-1.95410738955329,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.773077546323693,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.527602286806890,1.11049346247077,2.29257196360488) q[0];
u3(0.874779370210089,3.73975339908034,-0.373495938321222) q[6];
u3(0.656941144210016,-0.408403491345979,0.660240949033877) q[2];
u3(0.864906868960298,-0.738527413486217,-0.832653303455075) q[3];
cx q[3],q[2];
u1(3.29604969358964) q[2];
u3(-1.39838778930986,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.53022079146145,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.29305673405188,2.31118019855599,-1.98238926279568) q[2];
u3(2.46931662588835,0.433549697939163,5.18905473354838) q[3];
u3(0.978085496223349,0.369771615386415,-2.14234358698498) q[1];
u3(1.98011079520908,1.07451935655577,-3.86451497633595) q[2];
cx q[2],q[1];
u1(3.47771546370447) q[1];
u3(-4.04926867936802,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.519772578195894,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.04048141200962,2.52387533651499,-1.05909494456660) q[1];
u3(1.40590169148840,3.24274165778820,-2.70214153559223) q[2];
u3(1.46307165056728,-0.658985078420982,-1.81919965424557) q[4];
u3(0.532342218605291,1.40014429467960,-3.71521366900586) q[6];
cx q[6],q[4];
u1(0.437709216110019) q[4];
u3(-1.20375333448669,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.36291511603428,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.45843351889681,3.87955560807490,-0.852343297118501) q[4];
u3(1.92130684446247,-1.03469936106530,-2.12448638375340) q[6];
u3(2.10768573846227,2.90374970351271,-1.25551060063017) q[3];
u3(1.99514329329478,1.11891976963580,-1.28219994809243) q[0];
cx q[0],q[3];
u1(-0.379467400775609) q[3];
u3(0.195114970110192,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.97920974089482,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.01783899673061,0.328212939571094,-1.35486816833051) q[3];
u3(0.501643691294574,-4.04347351108763,-0.782914421189550) q[0];
u3(2.04207637457161,-0.0330646892251717,1.65504645880741) q[5];
u3(2.12016848336622,-1.51409995915305,-2.07211154647592) q[7];
cx q[7],q[5];
u1(2.92128700802957) q[5];
u3(-2.00375511924362,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.73490615893254,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.832051778957644,3.77561765573762,-0.696340855897902) q[5];
u3(2.90933803133328,-0.196087352685537,-1.43947700808396) q[7];
u3(0.526000413472528,1.97143573742434,-1.18596133244086) q[6];
u3(1.27759961323953,0.108155988248385,-2.42357283029047) q[4];
cx q[4],q[6];
u1(-0.308537444817642) q[6];
u3(-2.24782400399704,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.04400406543518,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.89453177302326,1.73889131418470,1.40636727744181) q[6];
u3(0.891026758600653,0.0418395883399161,-1.47815270252287) q[4];
u3(1.33758680311170,-1.83513239498736,0.744942570002462) q[0];
u3(1.45467552500752,-3.47090372723643,-0.0100122256166197) q[3];
cx q[3],q[0];
u1(1.23628058927021) q[0];
u3(-1.05274280073351,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.08907377106977,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.861346091625610,-1.98929488387230,4.18089154223927) q[0];
u3(2.52432161454044,2.71326103337534,-2.22998250483195) q[3];
u3(0.605223123776622,1.48753155634859,-2.86030698104132) q[1];
u3(1.49051838743371,-2.85698820901250,2.96821471176443) q[7];
cx q[7],q[1];
u1(2.41616680745348) q[1];
u3(-3.19371312444767,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.53056502495426,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.55314168004764,-2.28832371919348,1.62879397711237) q[1];
u3(1.94841473468488,0.858913039070683,1.87588402775677) q[7];
u3(1.91709472115448,-0.502821758188364,1.09660058441901) q[2];
u3(2.01840881486232,-0.388359252604326,-1.08150753066750) q[5];
cx q[5],q[2];
u1(4.24935200381092) q[2];
u3(-3.53139750549089,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.753442146417635,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.497592948669365,-0.0722601733196298,-3.19583436127554) q[2];
u3(2.83470374146743,-1.48612321381995,3.91923770007321) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
