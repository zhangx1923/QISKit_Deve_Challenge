OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.30564524571414,-1.87273649828459,0.0993143423342970) q[3];
u3(2.14931428865245,-1.69790853356627,0.0851739755623669) q[0];
cx q[0],q[3];
u1(3.52168508325683) q[3];
u3(-0.949453720980630,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.80686978165565,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.37992236338174,-2.15634811596603,3.36411506940878) q[3];
u3(1.85615504727011,-5.17691136366506,-0.496507918584201) q[0];
u3(1.45700456623677,3.47973609901930,-0.394032121738596) q[2];
u3(1.22966371913854,1.46448721563825,-1.38001773037624) q[1];
cx q[1],q[2];
u1(-0.120275333115057) q[2];
u3(-2.28967765272111,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.29161245535723,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.819559205747294,-3.12889599084427,-0.0684969591992111) q[2];
u3(0.911279038367651,1.65545515134450,-3.59569290047742) q[1];
u3(0.712290476876289,-0.834141334410507,1.58537767292277) q[4];
u3(0.349092328916326,1.88922876103682,-2.83090466948866) q[3];
cx q[3],q[4];
u1(2.24624528138919) q[4];
u3(-0.107638685109409,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.37104722414401,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.734583809568608,0.136635706684665,-1.80310563400375) q[4];
u3(1.87719094294013,-3.67910932526308,1.42989122881492) q[3];
u3(2.12651697763215,1.53648483978365,0.698351059070045) q[1];
u3(1.71143974245922,0.821777791077047,-2.22561832608822) q[0];
cx q[0],q[1];
u1(1.48082947842151) q[1];
u3(-3.77342793929330,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.12161869111482,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.75791812329418,-2.11292553548267,3.19787111565907) q[1];
u3(0.722833138483016,-3.63143424408649,1.81433132489103) q[0];
u3(0.993618951729719,0.260538702166426,1.56920986768154) q[0];
u3(1.43528650346123,-1.03259274431179,-2.53570598050189) q[2];
cx q[2],q[0];
u1(2.36460671854537) q[0];
u3(0.0755587116017642,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.36950127968701,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.05389320897189,-1.15617601665344,-0.602960821769166) q[0];
u3(2.40779294076618,-4.40783578960938,-1.74112648105353) q[2];
u3(1.87552326117942,-0.263738446881090,-1.67829281723253) q[3];
u3(0.931104332916917,-4.12373994581521,1.27998259187394) q[1];
cx q[1],q[3];
u1(-0.763430616422817) q[3];
u3(1.13188844050003,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.47956553496086,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.746054763587034,-2.18705712698226,3.34785492145729) q[3];
u3(1.28956501422544,3.09639550577378,-2.90214228709166) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
