OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(1.14586899021739,0.127398451818483,1.91301341617496) q[6];
u3(2.20950175693428,-2.02749548503016,-0.947488153945867) q[10];
cx q[10],q[6];
u1(3.03819483304805) q[6];
u3(-2.42450496216931,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.15139717685092,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.05798423635951,-2.07189623032724,2.67393097205252) q[6];
u3(1.60658559045827,0.148552278215413,0.627470570423662) q[10];
u3(1.84827854170888,3.13407051739360,-3.05483234668303) q[8];
u3(0.552841993773134,2.91000201498009,-2.11675176868149) q[13];
cx q[13],q[8];
u1(0.742815165827642) q[8];
u3(-0.193012249208113,0.0,0.0) q[13];
cx q[8],q[13];
u3(1.60429127126785,0.0,0.0) q[13];
cx q[13],q[8];
u3(1.57711008635055,3.38853787135549,-1.86768184800474) q[8];
u3(2.19434963325996,-2.96704165842992,-3.02091895726964) q[13];
u3(2.14894436930298,0.693093424739702,1.03218731783220) q[9];
u3(0.506014400765562,-4.97119810576810,-0.259580099113004) q[2];
cx q[2],q[9];
u1(2.40645843922898) q[9];
u3(-1.75742690621214,0.0,0.0) q[2];
cx q[9],q[2];
u3(0.642743968401995,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.60588926502942,-1.32812373055317,-1.78214810237851) q[9];
u3(0.809035133894489,4.00121031623018,-0.948462539018538) q[2];
u3(0.902764998666229,-1.33965697452009,1.71752516191072) q[3];
u3(0.404796241875408,-0.501876747416522,-1.56257892860740) q[0];
cx q[0],q[3];
u1(0.476313816002045) q[3];
u3(-1.50894907677037,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.10638604720830,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.83818782045032,-1.96151130600585,4.17834108635424) q[3];
u3(1.68021581856307,-6.19851543105095,-0.0455904560354479) q[0];
u3(0.612810365851840,-0.879001075433275,0.505357649201283) q[1];
u3(0.370095705365426,-3.33294962096011,2.74770550881043) q[11];
cx q[11],q[1];
u1(1.60669868329131) q[1];
u3(0.0782745474500528,0.0,0.0) q[11];
cx q[1],q[11];
u3(2.55351407340880,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.43092591718632,3.61907191412779,-0.315760241990844) q[1];
u3(1.10016122660550,1.19818646516862,-3.46844080993811) q[11];
u3(1.22775058321107,-1.97411582790398,0.670298553916384) q[12];
u3(2.23325042822524,-3.99715626265371,0.467050523955292) q[7];
cx q[7],q[12];
u1(1.26139403733222) q[12];
u3(-0.0528489825127962,0.0,0.0) q[7];
cx q[12],q[7];
u3(2.20138259351225,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.113003445443643,-2.12467617282934,3.65302048560336) q[12];
u3(1.50666535993195,-1.12195838073781,4.17072936620604) q[7];
u3(1.76801558696413,-2.14443102615432,3.58199568405739) q[5];
u3(0.720750025146261,-0.878905401566630,2.34216286284806) q[4];
cx q[4],q[5];
u1(1.62209134661893) q[5];
u3(0.0760716241140917,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.13041164019152,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.54684050162112,3.99579662160835,-0.680847515430760) q[5];
u3(2.19883192185829,3.23083833702079,2.80015854973264) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13];
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
measure q[12] -> c[12];
measure q[13] -> c[13];
