OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(0.377184757109965,1.68350946848210,-2.47068907501412) q[13];
u3(0.420022245515392,1.24336754109523,-2.62451027914714) q[9];
cx q[9],q[13];
u1(0.711983580564005) q[13];
u3(-3.23678108536677,0.0,0.0) q[9];
cx q[13],q[9];
u3(1.81593510576784,0.0,0.0) q[9];
cx q[9],q[13];
u3(1.90833401543333,0.944761463972583,-4.39561655307634) q[13];
u3(2.23282654894691,-2.91014761782660,-0.0739265264474529) q[9];
u3(1.06917725449022,2.71601816969703,-1.23359125903873) q[6];
u3(0.803137239825629,2.45013699645157,-1.96544283420402) q[3];
cx q[3],q[6];
u1(-0.781339431243775) q[6];
u3(1.19987823975839,0.0,0.0) q[3];
cx q[6],q[3];
u3(4.13093166684740,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.446130791418225,1.72050401223407,-2.00502449367459) q[6];
u3(0.427830907767095,2.01950328118769,-2.42828969552085) q[3];
u3(1.83295516934344,1.55361490449654,0.0117122171344889) q[5];
u3(2.27793992404602,0.641105062105658,-1.50959266345569) q[4];
cx q[4],q[5];
u1(0.604120560111317) q[5];
u3(-1.11824373648876,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.48505686997875,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.956406509896030,-1.12231857195643,-0.242446085515606) q[5];
u3(0.567066890084222,-1.14963988705507,-0.543642797704195) q[4];
u3(2.19686662151541,-0.505605779791812,1.60846320782265) q[12];
u3(2.51771971056783,-1.42573496535136,0.364216378425341) q[10];
cx q[10],q[12];
u1(1.87233350629234) q[12];
u3(-2.85144181806412,0.0,0.0) q[10];
cx q[12],q[10];
u3(0.840097201824279,0.0,0.0) q[10];
cx q[10],q[12];
u3(0.556947113361222,-2.68898233438409,3.30757142242344) q[12];
u3(0.894297072343543,-1.11351137265522,1.01014782665615) q[10];
u3(1.59278398581948,-0.00570445268566455,1.08030479189114) q[0];
u3(1.50524438461910,-2.51123665034162,-1.95854040964241) q[11];
cx q[11],q[0];
u1(0.886764619276286) q[0];
u3(-1.57136560719723,0.0,0.0) q[11];
cx q[0],q[11];
u3(-0.550648016949123,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.34238538707215,3.25421810234378,-2.52549958895940) q[0];
u3(1.10104759066849,-0.674835664386522,1.39918335315398) q[11];
u3(1.79051743662014,2.63243504054394,-2.59609785346220) q[2];
u3(1.66018090833760,-3.02997628617464,2.72813353827514) q[8];
cx q[8],q[2];
u1(0.468341073631396) q[2];
u3(-3.34643695606545,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.75291390526016,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.58972059857611,-0.142922037089006,1.24612998552743) q[2];
u3(2.23870367840403,-2.55208039510021,0.197640483182605) q[8];
u3(1.07562498590612,-3.87619729949214,1.15748852256006) q[1];
u3(1.53893656976815,-1.28468159023910,0.833561694428751) q[7];
cx q[7],q[1];
u1(1.18291360289270) q[1];
u3(-3.46172306906800,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.18937216998258,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.712503090483607,-0.358135460260072,4.30370525091589) q[1];
u3(2.21159688827264,-1.29956506659602,0.658954121478953) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13];
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
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
