OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0.828025373737698,-0.0361822444089779,1.65337324127222) q[1];
u3(0.757598469953920,-0.890221776947457,-0.756395040224237) q[6];
cx q[6],q[1];
u1(0.630551595309588) q[1];
u3(-0.935009685631675,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.79652428843554,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.34600610190984,1.69877230427432,-3.86735705585482) q[1];
u3(1.18918124497565,-0.125649970503505,0.225974564154014) q[6];
u3(1.68651605564757,1.51806013628549,-3.14974442994510) q[0];
u3(1.89080411342001,-1.71299034879564,3.84307472009509) q[4];
cx q[4],q[0];
u1(0.461545420166565) q[0];
u3(-1.01105314142225,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.06861385052465,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.46508326270679,2.79433105119373,-1.15070749486755) q[0];
u3(2.16277573140643,1.98190094100497,-3.94701206140648) q[4];
u3(1.48440371124693,-0.808905521013425,-1.96314007129054) q[3];
u3(2.26693468865054,0.797117430796298,-4.73476175803719) q[5];
cx q[5],q[3];
u1(3.25683141191960) q[3];
u3(-1.47547135955825,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.507992250950788,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.51934314324547,-2.04945614445163,0.0226210181086415) q[3];
u3(2.25825307267626,2.51983696633592,-3.67732924826438) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
