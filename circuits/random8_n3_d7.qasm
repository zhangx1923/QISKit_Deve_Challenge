OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(1.96301575944058,-3.44346844882379,2.71488605689238) q[0];
u3(2.21541239159684,2.84397306331367,-3.18394972544155) q[2];
cx q[2],q[0];
u1(1.87249943398068) q[0];
u3(-2.09493273944997,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.0283168192874879,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.47192323717589,2.53400376696927,-1.62718816995957) q[0];
u3(1.33863706605717,-5.20898675986997,0.357549496508157) q[2];
u3(0.675901923096957,1.45406467931457,-2.35195418694166) q[2];
u3(1.19299492215136,2.04648823421442,-3.51230621070601) q[1];
cx q[1],q[2];
u1(2.17619678388374) q[2];
u3(-3.29166477805298,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.635119375931512,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.19862359477963,-0.714787264688403,-0.801069705753418) q[2];
u3(1.62112669770625,-1.94763589415848,4.32177269578723) q[1];
u3(0.693860904101472,0.0651771409063399,-1.18161598687343) q[2];
u3(0.399225824310258,-1.48378782476828,0.148424368866757) q[0];
cx q[0],q[2];
u1(1.68230532120885) q[2];
u3(-2.11400755911046,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.0220454066457381,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.63817201006017,-0.634345434603960,3.48152705570779) q[2];
u3(2.16914200557587,0.753236344140720,-0.565991917291683) q[0];
u3(0.125124965231160,0.898060699544508,-0.306949466441903) q[2];
u3(0.570001446004770,-2.49107185211797,0.657680633781542) q[1];
cx q[1],q[2];
u1(0.970233095771237) q[2];
u3(-3.07296764485853,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.61804646224465,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.944926914366750,-2.74105637208899,1.31288210389359) q[2];
u3(0.698337486014894,0.680918506568046,0.509873775708278) q[1];
u3(1.69587109764298,1.29788774137120,-0.204842159747904) q[0];
u3(1.87954159465717,-0.593662496570586,-3.88545319893687) q[2];
cx q[2],q[0];
u1(3.30998323177999) q[0];
u3(-0.702682081317110,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.87863316307308,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.95570828353724,-1.32023971000208,-2.41922984599241) q[0];
u3(1.42808288946797,-0.198826952732894,-2.50801648065549) q[2];
u3(1.88862090670197,0.0774838852035491,1.08369620562171) q[1];
u3(1.74570914836876,-1.98688647703206,-1.39931234595594) q[2];
cx q[2],q[1];
u1(-0.0729953277540740) q[1];
u3(-1.52765813366792,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.598050083247094,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.75255757327732,-2.49369649034989,0.452393631941123) q[1];
u3(0.705717485885226,1.51594888311924,4.68693924592265) q[2];
u3(2.93506299812310,-2.74220750194965,0.954903862194645) q[0];
u3(2.43534651033694,-1.05175013150106,0.803021940687052) q[2];
cx q[2],q[0];
u1(0.511707193580959) q[0];
u3(-1.66445396672433,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.42850868248074,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.872437576342543,0.402519327438829,-1.63441313508185) q[0];
u3(1.31749801432504,2.24272049084885,1.95990857047752) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
