OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.89339799193809,0.191197577498529,1.90214109415105) q[0];
u3(0.955614035990764,-1.37398345102104,-2.80684164785093) q[1];
cx q[1],q[0];
u1(-0.180188142451818) q[0];
u3(0.653507585640948,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.94280347299827,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.65590801706525,-1.59208805327277,-1.50747630946300) q[0];
u3(0.318200794818437,-4.13619848677877,0.0172983490770227) q[1];
u3(1.82844245296368,0.482557906044342,0.714475713788872) q[2];
u3(1.29653948596706,-2.32230481645121,-1.60387423253517) q[3];
cx q[3],q[2];
u1(0.428069505026934) q[2];
u3(-1.69900677234583,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.11033653126528,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.521535783707477,-0.502861904827727,0.242736585430664) q[2];
u3(2.64016200313672,-1.67049684493962,1.89206111114581) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
