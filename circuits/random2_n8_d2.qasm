OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(2.32341794330205,1.81850353056781,-4.38326963759609) q[4];
u3(1.70978366679146,-2.50506455638093,3.52868284639317) q[7];
cx q[7],q[4];
u1(0.936752174431559) q[4];
u3(-3.39550312413839,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.98772953143197,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.89340892503938,3.55280185588090,-2.14852640645179) q[4];
u3(0.768467913937544,-1.33634320018274,1.41257511880096) q[7];
u3(2.53260995615939,-3.24136344505803,2.87956309619043) q[0];
u3(0.907428570388019,2.79059798148949,-1.65967208049657) q[3];
cx q[3],q[0];
u1(1.76110527413839) q[0];
u3(0.641776930943187,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.947097476983555,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.716103754279540,0.554339933175247,0.736440410322932) q[0];
u3(1.09495252493851,3.86477839413636,-1.40590118432481) q[3];
u3(2.23875262990071,2.78406845831951,-0.549026510996905) q[2];
u3(2.58791849527501,1.73716435421854,-1.66020003204531) q[6];
cx q[6],q[2];
u1(3.20464398720990) q[2];
u3(-0.496427702556003,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.93977645625488,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.68638225176393,-0.118961639848861,-1.51147577032012) q[2];
u3(1.03372412477203,-0.511629614935295,-3.42717162910158) q[6];
u3(1.54305450460951,-0.146201921418188,0.261448676653989) q[5];
u3(1.04011457834756,-1.87497723566997,-1.55897765939761) q[1];
cx q[1],q[5];
u1(0.276701990398158) q[5];
u3(-1.43742787877251,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.98375545761970,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.49511803764914,1.58531523707085,-1.26502303221734) q[5];
u3(1.83381431587237,-1.02228030120163,4.85780862071597) q[1];
u3(0.598524221993640,2.13067839153152,-2.23561346718275) q[2];
u3(1.03516535806358,1.42036258119777,-1.95186207721082) q[6];
cx q[6],q[2];
u1(1.54851183074049) q[2];
u3(-0.687076986296273,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.383990019981965,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.18084326684893,2.93578987237580,-3.29500650870949) q[2];
u3(1.62371224426797,1.04554944708031,1.85243806476754) q[6];
u3(1.63403392383366,0.571980793810821,1.20083081105133) q[4];
u3(1.79115106026747,-1.92621877356221,-2.32452054069748) q[7];
cx q[7],q[4];
u1(1.71464339190324) q[4];
u3(-0.817345927794281,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.121794590102287,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.822323760952407,0.00925130467147706,-0.715314191631938) q[4];
u3(1.11376276794168,1.34819559738441,-3.62399346781070) q[7];
u3(1.98024916124640,3.43824138028337,-2.14227693747525) q[1];
u3(2.54157391949486,1.78385575777350,-1.91022439087067) q[5];
cx q[5],q[1];
u1(0.383351500499725) q[1];
u3(-0.858441601809121,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.40971581700453,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.09678548111662,-2.38647758893449,-1.69046008202769) q[1];
u3(1.30974942375404,0.0739504112210603,4.77962394761765) q[5];
u3(2.89376839491470,2.12315535380219,-0.626553213828036) q[0];
u3(1.65949869484875,0.589922661555861,-3.70818900310890) q[3];
cx q[3],q[0];
u1(4.35282329884634) q[0];
u3(-3.80787918007504,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.408991017444326,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.78004895063100,-0.982691343441871,0.875056549832656) q[0];
u3(2.01682274010250,2.92009380124581,2.38981809665607) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
