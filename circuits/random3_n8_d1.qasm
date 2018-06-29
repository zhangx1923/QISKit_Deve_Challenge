OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(2.35196297375287,-2.38330315303692,2.31789772669229) q[4];
u3(1.84188084479294,0.0500953433572743,0.0910992423053599) q[5];
cx q[5],q[4];
u1(0.324572054620364) q[4];
u3(-1.18033268678029,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.51314247213466,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.09637221377566,2.83286879586072,-2.18409725494163) q[4];
u3(2.25871218586935,-4.65241918153056,0.684134062474325) q[5];
u3(1.76439191676416,-0.171983517739698,1.85614559995745) q[3];
u3(1.51002985259740,-2.12220964571363,-2.77922895636182) q[6];
cx q[6],q[3];
u1(0.843687419434680) q[3];
u3(-0.290921909740110,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.93294518143392,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.65239577957795,-0.466155268325892,-3.71537068233915) q[3];
u3(0.898978262817428,-1.23795502217129,-4.92688015573809) q[6];
u3(0.850508052079795,1.58149616553108,-3.24300017876266) q[7];
u3(1.84058391626695,2.21658834113104,-3.97545822705414) q[0];
cx q[0],q[7];
u1(-0.184811068503369) q[7];
u3(-1.95277127289309,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.666329069012850,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.898441339464624,-1.06453907180611,0.371023876649843) q[7];
u3(1.72908664322732,0.729188618074803,5.24895982293413) q[0];
u3(2.38769962141021,2.36296070625398,-3.54862185517621) q[1];
u3(0.966201594082698,2.42467258422417,-1.53882433824538) q[2];
cx q[2],q[1];
u1(1.46693576025308) q[1];
u3(-0.937402540694290,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.06420578568271,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.74730144734121,-1.22404807791780,1.07961655468533) q[1];
u3(2.44585088359073,0.00344247051832591,0.862690536696566) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
