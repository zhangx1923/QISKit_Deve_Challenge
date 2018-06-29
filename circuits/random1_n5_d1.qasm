OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.28547089794849,0.498635085259115,2.57311786522518) q[1];
u3(1.73115138707912,2.44636290472959,3.05686679251504) q[3];
cx q[3],q[1];
u1(-0.257726716585263) q[1];
u3(-2.17218916732755,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.38000589689513,0.0,0.0) q[3];
cx q[3],q[1];
u3(3.10548940197861,-2.31131161925421,2.73167749355291) q[1];
u3(0.479266146449861,0.473008088126857,-5.28044799577304) q[3];
u3(0.447154133411855,3.33514715635010,-2.13433635038709) q[4];
u3(1.71862929562599,0.836898085882029,-1.93976107544292) q[0];
cx q[0],q[4];
u1(1.66042688909070) q[4];
u3(-2.89031534934732,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.580473855199463,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.24718012066428,-0.149223546580052,0.826817094358188) q[4];
u3(1.87762295144631,-1.92141379880100,4.10044532057567) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
