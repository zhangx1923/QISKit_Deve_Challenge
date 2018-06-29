OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.687778128965431,0.453650480639642,-3.26279301857920) q[3];
u3(1.48906000688375,3.22798501720115,-2.55695624515345) q[1];
cx q[1],q[3];
u1(3.23378951625018) q[3];
u3(-1.53130117689620,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.36462066056842,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.63717498431192,-1.84747623339412,2.11790387109710) q[3];
u3(1.60565203864941,-3.67007794517872,-1.90157715000933) q[1];
u3(1.43239313990532,0.959147407862324,-3.19133889366722) q[4];
u3(2.03756512480664,2.92270672856294,-2.87214298071654) q[2];
cx q[2],q[4];
u1(-0.0408999083195241) q[4];
u3(-2.27940466616226,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.34471376972441,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.44890692208947,2.10611163705145,1.06917327645364) q[4];
u3(2.17107660843305,3.54893456310601,0.583623047015085) q[2];
u3(1.91102475439785,1.81796639478752,0.401214925996100) q[0];
u3(2.23283154188075,0.288634563123118,-2.02391534611381) q[1];
cx q[1],q[0];
u1(0.480390066056317) q[0];
u3(-1.28881790303558,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.165173490803597,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.66095552893873,1.73180311493585,0.115793862821300) q[0];
u3(2.48572278359685,0.725288174386612,1.76928831574058) q[1];
u3(2.08495985830096,1.87656622112024,-1.14563468687691) q[2];
u3(3.06173545081671,4.77691542696340,1.42641687746788) q[4];
cx q[4],q[2];
u1(1.77647625786722) q[2];
u3(-2.83515952532947,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.0864968002396997,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.94280078517972,1.46317994475029,0.283098363730586) q[2];
u3(0.505341813796992,0.435361292131894,4.30018199994985) q[4];
u3(1.24632567222116,1.03758553358858,-3.15993698348126) q[1];
u3(2.00029386527371,2.76178109463364,-2.98360754713753) q[3];
cx q[3],q[1];
u1(-0.170620538038630) q[1];
u3(-1.87932042233788,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.16094878662352,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.938348553096616,-2.00767809874612,3.39969905890335) q[1];
u3(1.83862988794636,0.745876705268254,5.49756524938378) q[3];
u3(2.02356356655181,-0.259042233053229,1.37590312872109) q[0];
u3(1.45243756908195,-2.11048608094007,-0.414559906451266) q[4];
cx q[4],q[0];
u1(2.88105693213584) q[0];
u3(-1.83004247264359,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.618818844845201,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.97818222499597,-2.77642409034962,1.62025741554718) q[0];
u3(1.92881721072229,0.275417574218318,2.34795175661338) q[4];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
