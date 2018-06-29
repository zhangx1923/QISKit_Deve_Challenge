OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(1.27512492774617,2.43471806146854,-2.19135179366191) q[1];
u3(0.526366112064364,-3.19102377449510,2.84786935525137) q[0];
cx q[0],q[1];
u1(2.36199118075978) q[1];
u3(-3.15418512550790,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.54134552988007,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.23267997890810,2.95902067094742,-2.44606001010232) q[1];
u3(2.23394989049040,2.63395815430357,-1.42089793382282) q[0];
u3(1.36886880301226,3.27897117892735,-0.725843518674730) q[0];
u3(2.44748371162504,2.79330619143244,1.03443578983758) q[2];
cx q[2],q[0];
u1(2.70709897313333) q[0];
u3(-1.95580566090895,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.355034438853412,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.57558368862934,-1.28876661248769,0.540157556743184) q[0];
u3(0.624547425049106,0.459431223407152,4.09736202432452) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
