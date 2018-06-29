OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.19331545338740,-1.20262818356576,1.97771567178939) q[0];
u3(1.81526404373082,-1.45163223082155,-2.29180106571759) q[6];
cx q[6],q[0];
u1(2.87018424999333) q[0];
u3(-1.83216904001207,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.38411123115245,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.76945061765917,-0.000997002798129504,0.246443125501940) q[0];
u3(2.10213959532011,-2.32434040452964,3.95055979668587) q[6];
u3(1.01825981629500,0.576548512830961,0.557956326794739) q[5];
u3(2.14143260584997,-0.680263743656154,-2.57547183186889) q[2];
cx q[2],q[5];
u1(1.43142093814379) q[5];
u3(-0.459521768782801,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.311450486515327,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.133737803053155,4.42851747410663,-1.59622313453161) q[5];
u3(2.30626812565195,1.07248698004132,3.84612797276856) q[2];
u3(1.47582879340912,1.54048188245653,-0.604279867617462) q[4];
u3(1.94250045523416,-0.0454217879842180,-3.44832960038131) q[1];
cx q[1],q[4];
u1(0.431232398612186) q[4];
u3(-0.260442254341856,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.11051236140104,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.94894635561677,-1.76049849264162,-0.0820459958488078) q[4];
u3(1.26759762395662,-1.50407771326045,-4.38951931768042) q[1];
u3(1.30380566841943,0.0276942874932233,2.08046684208913) q[6];
u3(1.69046462047223,-0.695481304768155,-0.658218486450999) q[0];
cx q[0],q[6];
u1(2.71544537288657) q[6];
u3(-1.89919702768631,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.906560712848214,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.70607929406137,-1.00219860049871,3.16676140591822) q[6];
u3(1.55404933897350,-2.96882150225235,2.83189820099350) q[0];
u3(1.46675231359594,0.652491198729166,-2.89279735804441) q[3];
u3(1.06009857417591,-3.07791492435953,2.71102266139297) q[1];
cx q[1],q[3];
u1(1.68162823235863) q[3];
u3(-0.498421038929812,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.33655697563558,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.541006377418235,0.510707123380010,0.291311602982709) q[3];
u3(0.690067338984249,-4.47349217135696,0.340608979006994) q[1];
u3(2.00748591681796,1.56614382294820,-4.17526738543873) q[4];
u3(1.97870661216285,1.91846480352591,-3.24024573051411) q[2];
cx q[2],q[4];
u1(2.19922387422981) q[4];
u3(-1.53502275302657,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.561945223926595,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.18310128962134,-0.135973270520945,1.96458497555268) q[4];
u3(2.78291537240351,-2.04099691932409,-0.818487940151076) q[2];
u3(0.626905363937532,-2.05319879127153,1.65851379698522) q[3];
u3(0.915269773844499,-2.79308221639677,0.841389020252254) q[5];
cx q[5],q[3];
u1(2.72375289024254) q[3];
u3(-1.85723480683761,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.100925901996356,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.37555407666348,-2.05785853831814,1.34513558746208) q[3];
u3(2.31878824676294,-1.18534034755906,-3.73534154962099) q[5];
u3(0.992292217639324,0.403934872207393,-2.46879098523057) q[2];
u3(1.33799494591235,3.10291747214711,-2.56931286031713) q[6];
cx q[6],q[2];
u1(-0.0290377620450137) q[2];
u3(-2.40460090038248,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.17528126399888,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.69036846892410,-0.0708381414250822,0.400799072400844) q[2];
u3(1.84718464123198,-2.93996431665003,-2.87376499140191) q[6];
u3(0.480998198917311,0.0645393565857307,-1.91143442896639) q[0];
u3(1.07938868546234,2.10264470661470,-3.98486097084913) q[4];
cx q[4],q[0];
u1(3.60705395627713) q[0];
u3(-0.757468956935838,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.73805881163297,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.09429497729535,-0.453311045328671,2.09045053652410) q[0];
u3(1.25095273801644,1.62845281093939,-4.62025326721224) q[4];
u3(2.59502300433029,1.61129234137486,-2.17337793098114) q[4];
u3(1.56043353924426,2.25975644222930,-2.98608220758274) q[5];
cx q[5],q[4];
u1(-0.110289583732972) q[4];
u3(-2.08237357846928,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.64848440728976,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.71646851855598,-2.98253324835090,-0.398770337886836) q[4];
u3(1.69645445111250,-1.68274432652745,-1.09937104836535) q[5];
u3(1.82298731292849,1.03369865257138,-3.55737661431386) q[6];
u3(1.15017439963099,2.32265711222540,-1.99434526896922) q[3];
cx q[3],q[6];
u1(-0.0189757932648111) q[6];
u3(-1.92251384171836,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.57952754351255,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.909063891132314,-2.05468487034103,-0.528235813784472) q[6];
u3(1.86920982579804,-0.998482183319627,-0.761586922135490) q[3];
u3(1.89574475494844,-2.97865422671276,0.448258027514399) q[0];
u3(2.70278142235255,0.0927727308151862,1.91514561377170) q[1];
cx q[1],q[0];
u1(1.75027970327157) q[0];
u3(-2.53561631025044,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.928806592560611,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.37165205003973,0.179919763909853,-0.389388252380796) q[0];
u3(0.200017019904000,-1.87439771280536,3.98209743606196) q[1];
u3(1.66819074122528,0.631372941350281,-2.44264479410550) q[4];
u3(2.62311313794316,-2.64632415506291,3.61034131064600) q[5];
cx q[5],q[4];
u1(2.38500237996163) q[4];
u3(-1.80298620065418,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.306062851404350,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.20045096364276,-0.258301235638600,2.40791605612142) q[4];
u3(1.03996504947681,4.13104896431691,-2.15174493403525) q[5];
u3(1.25234170965760,1.87892226393234,0.247478659130887) q[6];
u3(1.51320972361562,0.278872381901582,-3.38686472754843) q[2];
cx q[2],q[6];
u1(2.32326108076457) q[6];
u3(-2.09383935900495,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.119907767083580,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.74452015105169,1.66533486482422,-1.70468373187971) q[6];
u3(2.86352461922882,-2.27351995147311,0.200977437523818) q[2];
u3(1.49900771271804,1.18752629128668,-1.49034579455767) q[0];
u3(0.492770061975838,-4.45442895902103,1.17684325989601) q[3];
cx q[3],q[0];
u1(0.981767285621502) q[0];
u3(-0.370598473452227,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.66290933700983,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.75812824320114,0.682589069971041,-2.36506936876453) q[0];
u3(1.67648075768669,3.01175357850318,-2.19118875280526) q[3];
u3(0.792721796873016,0.639595455229572,0.420496164424848) q[4];
u3(2.14392164198009,-1.08467372715149,-3.83634819363690) q[3];
cx q[3],q[4];
u1(3.83816016543874) q[4];
u3(-4.11342711239033,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.111160665974293,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.933618182348122,-2.36630015841110,0.900475432607644) q[4];
u3(2.61874195200573,2.76641765730653,2.27045811986657) q[3];
u3(0.848561742040519,3.12462661134679,-2.78173897064130) q[6];
u3(0.916349528697796,1.28925087235023,-1.59276196361145) q[0];
cx q[0],q[6];
u1(0.888500537019120) q[6];
u3(-0.138767364596789,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.87832279534259,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.46358692011088,-2.28555496699925,3.03505650740917) q[6];
u3(2.32025788891201,-0.467917951116096,-0.618526945554391) q[0];
u3(2.47538190432990,-0.433112920529351,1.88133165986054) q[2];
u3(2.67689735604677,-1.68674179081052,-0.916011761704026) q[1];
cx q[1],q[2];
u1(3.72832346224180) q[2];
u3(-1.40208692917785,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.31414481032495,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.00760201682260,0.495027610770633,-0.546945234903632) q[2];
u3(2.99439963871932,-1.12276362247789,1.90796755967948) q[1];
u3(1.73321617570519,-1.29566979875921,1.31443232769854) q[3];
u3(1.00743357908585,-1.96248774393516,-0.108826224782218) q[0];
cx q[0],q[3];
u1(2.51437412252960) q[3];
u3(-2.90489485282215,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.32595136680957,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.36628054306017,-0.0212879508687160,0.504976958575914) q[3];
u3(1.26793967901000,-0.552486383623683,0.438297384160871) q[0];
u3(0.531246883436610,2.17234379737800,-0.368438669519079) q[2];
u3(1.53379370656985,-0.452618454224497,-4.02201435045608) q[4];
cx q[4],q[2];
u1(0.741365130228083) q[2];
u3(-1.47999820600900,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.100205698538379,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.64068069247848,0.741929627141550,-2.13639567312248) q[2];
u3(1.03080014811069,1.41621722974893,3.76610581409160) q[4];
u3(0.349664222988767,-3.25836636841118,2.85087770035763) q[6];
u3(0.700684468524672,-3.47634368122020,1.53739157282734) q[1];
cx q[1],q[6];
u1(0.0504048304522524) q[6];
u3(-0.852930227094420,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.56054385522288,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.49497116003592,2.28062059232840,-2.35952967118369) q[6];
u3(0.427646089495763,2.48884155357022,2.09890562972417) q[1];
u3(1.75996655506634,0.630583968365534,2.18706135625514) q[3];
u3(0.891642427226421,-2.78211728402853,-3.27705466016535) q[1];
cx q[1],q[3];
u1(-0.137445018959670) q[3];
u3(-1.41288242060921,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.74695801077644,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.77908768342350,3.01220146897899,-2.87742789455407) q[3];
u3(1.26524856493251,2.71288963521443,-3.29864346066002) q[1];
u3(0.514935110774130,2.93083954042036,-0.999080427929067) q[4];
u3(1.96147169548520,3.25270417127484,1.70257063272256) q[0];
cx q[0],q[4];
u1(0.791206606353576) q[4];
u3(-1.43989619582393,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.375336692900898,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.01457984806701,-1.17143933950745,4.05042256154344) q[4];
u3(2.61607894846181,2.64359411554040,0.688383009243973) q[0];
u3(2.09266295687326,1.57038346088616,-3.04640374473262) q[6];
u3(1.85967761665935,-2.26556090129113,3.20554342589279) q[2];
cx q[2],q[6];
u1(0.960911607155556) q[6];
u3(-3.47589729154914,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.01100161221986,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.712096702430572,-1.45621918103495,2.18842376416413) q[6];
u3(1.49248235293130,-2.10119366185435,-2.20060671171877) q[2];
u3(0.428339078247367,-0.888540352569186,1.14861607392671) q[3];
u3(0.721561552579423,-0.819510602204400,-1.44846603987087) q[5];
cx q[5],q[3];
u1(3.24091269357979) q[3];
u3(-1.29081855660746,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.55367320885180,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.19852668760104,0.172253491809362,0.219336086975955) q[3];
u3(2.05549968912619,0.00751485179207068,1.62367776915185) q[5];
u3(0.448858676712036,2.03340177164588,-1.41678413517866) q[6];
u3(1.38753727497249,-2.58997222563421,0.886033168218129) q[2];
cx q[2],q[6];
u1(2.24520367194497) q[6];
u3(-2.60503300731092,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.278493432644803,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.994713572706848,2.71744827990934,-3.01067380015672) q[6];
u3(1.00739260398981,-3.79091002232227,-0.975904427953440) q[2];
u3(2.73604216891689,-1.54599418374579,1.68674800850051) q[1];
u3(1.73990777991928,1.67482323563554,2.47208948972385) q[4];
cx q[4],q[1];
u1(-0.457965790048427) q[1];
u3(1.29854326327735,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.82836754317328,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.19427530933483,-0.162160974294023,-1.82353342130105) q[1];
u3(0.688028800845204,-1.10259520731218,-1.10954720208306) q[4];
u3(2.13949674503655,3.59441611886909,-0.722314472479111) q[2];
u3(1.99351519243243,2.03006201723074,-2.15456229879179) q[3];
cx q[3],q[2];
u1(0.988910748511510) q[2];
u3(-1.09123269239733,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.32714301155249,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.65638738698937,0.487151573540381,-0.151022348418068) q[2];
u3(1.40139131248934,-2.57005421584216,-3.01406071305216) q[3];
u3(1.93636840748037,1.98136948007001,-2.05508644202354) q[0];
u3(2.03608766873535,1.95974425990906,-2.16773116387608) q[1];
cx q[1],q[0];
u1(-1.04131512122011) q[0];
u3(0.159963385304504,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.67532869230420,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.15157995166529,3.56914892252025,0.113625656714930) q[0];
u3(1.52688409771151,-2.05342778912545,3.87417259162590) q[1];
u3(2.15693776097150,0.118067355915257,-2.82475056284144) q[6];
u3(2.78798781416462,4.77826404860207,-0.961481003862042) q[5];
cx q[5],q[6];
u1(-0.291451298044340) q[6];
u3(-1.89422395594944,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.756253802874140,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.41081514891367,0.754307199542622,-2.11068350259444) q[6];
u3(1.31007976688006,-0.109037015191326,-4.00835122177799) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
