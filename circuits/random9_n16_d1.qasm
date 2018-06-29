OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(2.53229413142077,1.25763229662930,0.364722988418676) q[14];
u3(1.13592335117419,0.443788823292679,-4.16064051996001) q[13];
cx q[13],q[14];
u1(3.35084669668417) q[14];
u3(-1.25859168709649,0.0,0.0) q[13];
cx q[14],q[13];
u3(2.14730484176900,0.0,0.0) q[13];
cx q[13],q[14];
u3(2.61353789830315,4.12974872426816,-0.384705401408084) q[14];
u3(1.94143468348821,-1.09891130633392,0.0864063573496178) q[13];
u3(0.473080982016240,-0.114691733231767,0.377694370090956) q[6];
u3(0.701503433159529,-2.09821973770388,0.695283163977832) q[7];
cx q[7],q[6];
u1(1.25346717381760) q[6];
u3(-3.46858432342368,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.20817955500006,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.926711158248136,0.425929274669308,1.16985499813764) q[6];
u3(1.29666751056495,-0.757807118115000,0.549890764726760) q[7];
u3(1.30472710531674,-1.97637307319098,0.0567945116201558) q[11];
u3(0.564809695163609,-1.89818826636817,-0.173206865426684) q[1];
cx q[1],q[11];
u1(0.105304071963268) q[11];
u3(-1.41998900260755,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.75189339419576,0.0,0.0) q[1];
cx q[1],q[11];
u3(2.12185645654502,2.85316649695453,-3.13673903078339) q[11];
u3(2.31316132258719,3.36938985196575,-0.746247493337586) q[1];
u3(0.756672493084505,-0.290629069710034,-1.40658914068248) q[10];
u3(1.45649997175339,2.21486111037620,-3.51894202363100) q[5];
cx q[5],q[10];
u1(1.51180461868102) q[10];
u3(-0.205289970674343,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.31717064253862,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.11531259500652,0.793262903310302,0.369761188164844) q[10];
u3(1.26397983416945,-2.80624776009664,0.964374594732803) q[5];
u3(2.05445473049121,2.52858159024264,-0.794210627765389) q[12];
u3(2.75576789701333,-0.436027998138094,-5.02219477402476) q[2];
cx q[2],q[12];
u1(2.48777659358086) q[12];
u3(-1.90943729002682,0.0,0.0) q[2];
cx q[12],q[2];
u3(1.01563677343943,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.76627488856928,-1.05594259888773,-1.32821836757513) q[12];
u3(1.09447578911704,0.970630866087056,-3.60803335918116) q[2];
u3(0.476066625035391,-0.821172522961175,0.577406190217639) q[8];
u3(1.34499845951988,-3.29947362001053,0.328403460547372) q[9];
cx q[9],q[8];
u1(2.01973416610547) q[8];
u3(-2.20174671901409,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.664805450567265,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.864646535617138,-0.615432690230922,2.78405857038860) q[8];
u3(1.94801916997347,-1.93012468158433,-1.57376505209729) q[9];
u3(1.27676987975876,0.955823998409121,-2.25157684153269) q[15];
u3(1.65104039751592,-2.84474301695874,2.69756554529718) q[3];
cx q[3],q[15];
u1(3.71257736777637) q[15];
u3(-4.58718684054619,0.0,0.0) q[3];
cx q[15],q[3];
u3(-0.635972792794034,0.0,0.0) q[3];
cx q[3],q[15];
u3(1.18309591553881,-2.76469574372163,2.47892208358720) q[15];
u3(0.861378646224140,-0.681250783122201,0.971078223727507) q[3];
u3(1.76604204005908,-0.923614748161815,-0.292684890928882) q[4];
u3(2.79089525758346,1.60648342397082,-4.20059995181067) q[0];
cx q[0],q[4];
u1(4.37878485057876) q[4];
u3(-3.36351817970654,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.434718057692804,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.20372317717934,-0.332734574783519,2.92135578826945) q[4];
u3(1.27667959949268,3.52106092055877,0.218772509603800) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[14] -> c[14];
measure q[15] -> c[15];
