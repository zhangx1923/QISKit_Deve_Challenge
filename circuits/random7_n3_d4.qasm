OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(0.917873597574230,-2.41650957600743,3.45677643807559) q[2];
u3(1.10356111412381,1.17779989598869,-2.12647059595192) q[1];
cx q[1],q[2];
u1(1.21507432735955) q[2];
u3(-3.55459772513261,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.50986428968400,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.641456326761334,2.13783823837576,-2.58194264829903) q[2];
u3(1.96997085262262,5.00682713187593,0.989343731948349) q[1];
u3(1.34295709796371,-0.108788084081250,2.74552359118719) q[1];
u3(1.66849418990431,-2.55923127392973,-1.82117789803283) q[0];
cx q[0],q[1];
u1(1.00673080079306) q[1];
u3(-1.49642348800079,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.316720405758504,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.16714906406352,-2.47876000312603,-1.15615166018138) q[1];
u3(2.35773417894346,2.13158224219620,-2.30419917009008) q[0];
u3(1.87363071536090,1.09979437289372,-0.228383028722179) q[1];
u3(1.75673387007083,-0.705694601428989,-4.30489345912135) q[2];
cx q[2],q[1];
u1(1.82862525965724) q[1];
u3(-2.65088276628110,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.245448795054173,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.48848683821152,2.39083726140203,-1.94809512078804) q[1];
u3(0.572946301636230,0.0850222053425467,0.550734186881839) q[2];
u3(2.21972915686405,1.38725478353393,-0.133328802365065) q[0];
u3(2.45821429153552,0.771071732764227,-3.63730757311462) q[1];
cx q[1],q[0];
u1(2.95904565152982) q[0];
u3(-1.60810189818181,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.894351293900812,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.794624247498855,0.519632112249056,2.44206756826146) q[0];
u3(0.562753337220330,-0.411499068778183,0.805695895545268) q[1];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
