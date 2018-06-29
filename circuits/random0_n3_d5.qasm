OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(0.383147942490432,0.0913296563478839,-1.76434281511679) q[1];
u3(0.978564320073048,-5.34683445895417,0.719082460486854) q[0];
cx q[0],q[1];
u1(-0.465569366552312) q[1];
u3(-1.65974530806926,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.32172644664938,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.60609150244617,3.13429613379871,-2.34746222566812) q[1];
u3(1.00391825584593,0.501587273810664,-5.09058787435753) q[0];
u3(2.86866189227649,2.41754565962456,-3.68586946269476) q[1];
u3(1.38490898153676,1.32873524742524,1.24631263878006) q[0];
cx q[0],q[1];
u1(-0.393486041436950) q[1];
u3(-1.90748871176963,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.48907339024940,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.04995434196042,-1.11849043669946,-0.875032773707085) q[1];
u3(2.21017207126691,-1.95608151709945,2.16092137213597) q[0];
u3(1.53704829230621,3.56066403152370,-1.57909856303043) q[2];
u3(1.90979874726014,2.07915063399915,-0.0586440895855942) q[0];
cx q[0],q[2];
u1(3.87548370055209) q[2];
u3(-1.42462441957475,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.20807228783200,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.28677485275353,-3.19962194036290,1.60051224917128) q[2];
u3(1.21642224086908,-2.41942398990478,2.82226890227203) q[0];
u3(0.781369521454415,-2.97786863290677,-0.00385022652202505) q[2];
u3(1.96648167200759,-3.63264141442202,1.08192391515890) q[1];
cx q[1],q[2];
u1(2.39334103442754) q[2];
u3(-1.58140907169428,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.54591306271713,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.56338152321454,-3.16309006809538,2.90799719301965) q[2];
u3(1.68833677559215,0.594572799115800,4.63735866837726) q[1];
u3(1.41730217241845,-2.52294083474014,-0.503033559341291) q[2];
u3(1.91931352688627,-3.06056572865320,-1.16117357407454) q[0];
cx q[0],q[2];
u1(-0.813056654033339) q[2];
u3(0.267889068629575,0.0,0.0) q[0];
cx q[2],q[0];
u3(3.53897327189953,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.588469617862164,1.06531198309399,2.13954571346120) q[2];
u3(1.58642441169397,2.16687325474963,-1.62267224467680) q[0];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
