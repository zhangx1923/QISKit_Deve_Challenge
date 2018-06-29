OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.34583382094918,1.55224223037998,-0.582519567913510) q[4];
u3(0.589760709193644,-0.171181237084113,-2.66883983520875) q[5];
cx q[5],q[4];
u1(1.77021486567873) q[4];
u3(-2.20067333776952,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.290448619305378,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.32160672769389,-2.74406001683270,0.926536714964912) q[4];
u3(2.72107521157515,-0.322860895112523,-4.46118058112523) q[5];
u3(1.78279064228493,-0.299955436174209,2.29945002371390) q[0];
u3(1.70533275638426,-1.98279757130086,-1.16604665519265) q[6];
cx q[6],q[0];
u1(1.52396657119578) q[0];
u3(-0.338300906626988,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.14284550544827,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.31986940763687,-2.67206703194895,1.18293306863832) q[0];
u3(0.533352524473744,-2.41544036075461,-2.12220034996936) q[6];
u3(2.47989926069241,1.14036473692763,0.850747100167234) q[3];
u3(1.75606228707433,-0.562448816708368,-3.79645341723301) q[2];
cx q[2],q[3];
u1(3.47449749768341) q[3];
u3(-3.23236306831953,0.0,0.0) q[2];
cx q[3],q[2];
u3(-1.06067111426121,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.21001725652575,-0.411186836228508,-2.45836159438521) q[3];
u3(2.02457741330100,-0.848007026395340,3.40323537255488) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
