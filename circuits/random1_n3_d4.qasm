OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(2.26664243392572,1.93930213595698,-0.737139221571849) q[1];
u3(2.42553934864421,-0.777480810633764,-4.04855045919169) q[2];
cx q[2],q[1];
u1(1.44387289444426) q[1];
u3(-0.169456102905684,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.07837571506625,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.17534823009566,1.04298004522413,-1.92508534442228) q[1];
u3(0.107430041537102,3.60700045818921,-1.12240722416251) q[2];
u3(1.77951764360374,1.74601766540105,-0.316904709195881) q[0];
u3(2.90729709346351,0.584662025898677,-3.53046605423228) q[1];
cx q[1],q[0];
u1(1.56636818924925) q[0];
u3(-0.140216088141050,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.89274304312210,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.36669414810488,-0.875954923332507,-1.47817283246199) q[0];
u3(2.85102272229242,-2.47717141557486,0.355455700167369) q[1];
u3(1.58197711860315,1.22670317996568,1.79991472469521) q[0];
u3(0.977367530758704,-1.50560714537231,-0.511971879610025) q[2];
cx q[2],q[0];
u1(1.00739795987956) q[0];
u3(-3.24362816034736,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.31432635835818,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.68865248343359,-2.64690475907464,3.50305836216586) q[0];
u3(1.61580715773260,-0.846073355057002,-0.466073587481426) q[2];
u3(0.410400322030050,-3.14621356468032,2.83583404348055) q[0];
u3(0.398889518147673,0.317611486724866,-2.53057702782329) q[1];
cx q[1],q[0];
u1(1.46969362962436) q[0];
u3(0.219461975500303,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.486866150263787,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.00984894972939,-0.0248853988326956,2.35502779306159) q[0];
u3(2.73899058853324,-0.905821607571740,1.80305797791997) q[1];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
