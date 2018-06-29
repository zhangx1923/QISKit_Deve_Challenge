OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.70493324781274,1.29563951381143,0.595644461298640) q[6];
u3(0.270812408459769,-0.198300568111844,-3.86114085763336) q[2];
cx q[2],q[6];
u1(0.000160944640369642) q[6];
u3(-0.549415640928290,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.58115502082450,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.85155425216910,-1.77308725632284,1.38217738033997) q[6];
u3(1.34367993338734,2.64774813740953,2.73555049190341) q[2];
u3(1.43289845746817,-0.892051757979649,0.0656182431472125) q[7];
u3(0.270320467653874,-2.98279846311166,-0.351910344700530) q[8];
cx q[8],q[7];
u1(2.84716225590103) q[7];
u3(-1.88044570381019,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.702720772045127,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.621759005001087,2.39089705228858,-1.60397598910580) q[7];
u3(1.74049091907678,-2.14883372245564,3.85889316793069) q[8];
u3(2.67581774565374,1.40484302792417,-2.76345177050570) q[9];
u3(1.63668392902211,-2.56873915456900,2.42579067173451) q[5];
cx q[5],q[9];
u1(0.146724918737408) q[9];
u3(-1.25967751180480,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.54624883915951,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.35492213754304,0.443113359970896,-3.14620298390151) q[9];
u3(2.28917537704409,0.567715181928003,-3.12094882549567) q[5];
u3(1.77573350351335,3.54888381778530,-1.90971330824749) q[4];
u3(1.16333287130961,2.94747596823489,-2.49714592972001) q[3];
cx q[3],q[4];
u1(0.656781417916701) q[4];
u3(-1.18321129808061,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.57768523410555,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.59150735698600,-1.06244462057638,2.08156979001586) q[4];
u3(0.504783568874457,1.88218160722403,-0.776111734058629) q[3];
u3(2.40463063811761,-3.64090848523131,0.634642854388506) q[1];
u3(2.84641833741291,-1.99765471891501,-0.271956756614605) q[0];
cx q[0],q[1];
u1(3.00156932666872) q[1];
u3(-1.66006657715496,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.769743403549633,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.04058381817206,-0.468463097388442,0.00322802958971732) q[1];
u3(0.852981274286085,0.705553331555713,0.505925383882872) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
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
