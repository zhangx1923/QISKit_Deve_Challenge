OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.975080042867909,0.573171168522652,-1.12375699945140) q[3];
u3(0.932112255419601,-3.34147191263796,0.968616111042918) q[4];
cx q[4],q[3];
u1(2.05656505731145) q[3];
u3(-2.56131974293255,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.33164447380694,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.19415288020674,-1.21141326685946,1.30821542044848) q[3];
u3(2.08668684483924,-4.80907395647340,-1.47043457314548) q[4];
u3(0.628948424397409,0.627166098172935,0.0874196479544453) q[0];
u3(1.07237749252246,-0.447917064302853,-1.49169939697926) q[1];
cx q[1],q[0];
u1(0.641616328803957) q[0];
u3(-3.11923722351904,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.59537720332903,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.51777804812157,4.15589798245583,-1.85958768763008) q[0];
u3(1.71418557618384,-3.07526118584761,-1.21856415315403) q[1];
u3(1.14766426035684,1.25993645028100,-1.57616675651047) q[1];
u3(0.236270350516948,0.782317091872796,-2.76899410655368) q[0];
cx q[0],q[1];
u1(1.21219502461644) q[1];
u3(-0.167333612798486,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.21047788783360,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.489483395403522,-1.01435878634157,2.88112117641440) q[1];
u3(1.61497194197631,5.09140888083479,0.688587175241214) q[0];
u3(2.03749557656679,1.21262536925319,-1.08499389131525) q[4];
u3(0.869131748908568,0.638035454466986,-3.00674665981092) q[2];
cx q[2],q[4];
u1(-0.404428828380771) q[4];
u3(1.14959178189036,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.40000889462973,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.66878702816312,1.79269003914497,0.613663230637452) q[4];
u3(1.87450475602642,2.17458192470856,0.375147352212454) q[2];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
