OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(0.519995281518056,0.0273346848851256,-0.0311451738904169) q[2];
u3(1.36090149170349,-2.45073047577075,1.06488395898798) q[5];
cx q[5],q[2];
u1(3.42456062906159) q[2];
u3(-1.78719505852905,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.15451624416208,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.92044596893127,1.67163132335225,-3.18469122725381) q[2];
u3(0.868951784603174,0.588934380549165,1.09068519756277) q[5];
u3(0.727247054149703,-1.29669729765165,1.05161442384797) q[3];
u3(0.663496140253228,-2.55779619775088,0.860392040574567) q[0];
cx q[0],q[3];
u1(4.25234537338545) q[3];
u3(-3.55408861899304,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.637993720296903,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.19649249912969,-0.251260524563492,1.83798861317159) q[3];
u3(1.34301567534358,-1.48857380123282,-4.39936425824217) q[0];
u3(2.23346996553635,2.04060009904208,-4.20882800921742) q[1];
u3(0.839998940484308,-0.457483834570160,2.50520483420440) q[4];
cx q[4],q[1];
u1(0.781974616893081) q[1];
u3(-0.579823599799751,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.06012744256741,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.27175752638548,2.68276283142751,-1.12244371465547) q[1];
u3(1.65611146929808,3.16917764957526,0.0523687995877697) q[4];
u3(2.25187694416037,3.31386940315266,-0.552506408318382) q[3];
u3(1.80461964848515,1.46829960815843,-1.17713538371091) q[0];
cx q[0],q[3];
u1(0.432048625852772) q[3];
u3(-1.29788672642196,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.19551690660851,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.95225478029486,-3.93144861424347,2.28779271031570) q[3];
u3(1.93991209986152,4.97248229190773,-0.472136935429382) q[0];
u3(0.656092604496633,0.185365889694643,-2.36834643657400) q[2];
u3(1.80903691096481,2.25093891949346,-3.72915890695452) q[4];
cx q[4],q[2];
u1(0.709042710753956) q[2];
u3(-1.51729105999752,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.0714303947517232,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.40993709127569,4.47563577079618,-1.70529506767279) q[2];
u3(2.39353160644242,1.58969687199821,-3.58537489753656) q[4];
u3(1.01917080072168,0.863562040308321,-2.19678674572076) q[5];
u3(0.826603836244039,-3.67665161626106,2.27458061764667) q[1];
cx q[1],q[5];
u1(2.87960747991206) q[5];
u3(-1.52013259714458,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.594786977895914,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.238987307943187,-2.77956965592460,0.780188317075233) q[5];
u3(0.887687351524102,-1.56346172870803,0.677303447004056) q[1];
u3(1.49320664556329,3.52890448210692,-1.57728389760361) q[5];
u3(1.64489660444657,2.33990376719910,-2.48621210204582) q[4];
cx q[4],q[5];
u1(0.897563933047494) q[5];
u3(-3.41860730234773,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.08847567575309,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.655720762262779,0.643475925917929,-1.66405467219551) q[5];
u3(1.50068128409949,-2.12636792785979,-0.389000734943433) q[4];
u3(1.52163017577852,0.526745316777648,1.94932014264678) q[0];
u3(1.64557133894770,-1.33841950916902,-1.01236041682745) q[3];
cx q[3],q[0];
u1(2.68167951912930) q[0];
u3(-2.97418141916677,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.952671390824512,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.92192799015331,2.37815358412068,0.603800884920434) q[0];
u3(2.27803268731220,3.12236492428199,-0.417235480153180) q[3];
u3(1.26206486565460,-0.0535029373928793,-1.98500645384530) q[2];
u3(1.86864333732425,-3.82559237725839,1.61039687579012) q[1];
cx q[1],q[2];
u1(2.10603054156603) q[2];
u3(0.110065057255454,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.14461000362508,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.03372937385719,-2.68770014114394,-0.202531592235490) q[2];
u3(1.44222216793238,3.91855138957205,1.31970837601204) q[1];
u3(1.68006871037766,0.952756503859912,-3.45688313908682) q[3];
u3(1.29033781490478,2.46090498005107,-2.60186791707037) q[4];
cx q[4],q[3];
u1(0.870888659687371) q[3];
u3(-1.39469224253448,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.72629839372400,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.694395359441573,0.186718210804742,-0.103854609983103) q[3];
u3(1.85441949309888,0.310601237846997,-1.26815647264832) q[4];
u3(2.03331781124086,-0.946285224227096,-0.634571207846777) q[0];
u3(1.06609450719661,0.982095841986196,-5.26360962532938) q[1];
cx q[1],q[0];
u1(1.65008312593647) q[0];
u3(-2.74767287016489,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.0445436755697501,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.47707260633760,-3.08655773462378,-1.11725764465748) q[0];
u3(0.641238060404570,-3.88316116371321,-1.77298715736738) q[1];
u3(1.16169721085964,-1.73012714586709,1.06269930805995) q[2];
u3(0.880919566779217,-1.53485265937079,0.154935697644071) q[5];
cx q[5],q[2];
u1(-0.0434822926942158) q[2];
u3(-1.69191057970427,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.49671388941382,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.92428651777067,2.17026782250637,-1.99821040748436) q[2];
u3(2.03539209250075,-4.63637567528386,-0.0923002185686008) q[5];
u3(2.13333074788795,0.561147166946610,-2.14025182851945) q[0];
u3(2.47196478756790,-4.43334583598565,1.37651095758050) q[3];
cx q[3],q[0];
u1(0.803967448550417) q[0];
u3(-3.13404909772620,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.80131167578502,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.30101069261457,1.78419270967358,1.89161523345501) q[0];
u3(2.32769718617762,-1.69335516449563,3.66084731736328) q[3];
u3(1.71413081149385,2.61812925181022,-1.28661848262645) q[4];
u3(0.805553648822760,1.04079757613726,-0.191670204770087) q[2];
cx q[2],q[4];
u1(1.80985658995573) q[4];
u3(0.181770038506836,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.38612866218316,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.861597248334059,-2.81751950852286,0.799687783148095) q[4];
u3(2.03677748495666,-0.225374976678075,0.369966990390872) q[2];
u3(1.16822319404914,-0.387332460480261,2.13514807629256) q[5];
u3(1.27579727145158,-0.773372766075899,-0.945109796700272) q[1];
cx q[1],q[5];
u1(1.52350761059906) q[5];
u3(-2.78234261493893,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.12013870776784,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.73952121180173,-1.78598112329290,3.54609077677725) q[5];
u3(2.16515981472741,-2.35917075475004,-0.249088042111359) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
