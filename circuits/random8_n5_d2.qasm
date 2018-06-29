OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.12365982890964,-0.896718237373297,-0.853087736879083) q[1];
u3(2.29954098430769,1.67514708428685,-4.14180905066578) q[4];
cx q[4],q[1];
u1(2.22694701362374) q[1];
u3(0.213288359162999,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.23818821992944,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.14303485985360,1.25749196644974,-0.0381250134873061) q[1];
u3(2.62892439648157,-3.29707350241975,-1.49068950754497) q[4];
u3(2.49831359074288,2.83950154259890,-0.330728078919143) q[0];
u3(2.85963258117893,1.79851510744754,-1.87606338638628) q[2];
cx q[2],q[0];
u1(1.85338719945504) q[0];
u3(-2.13406176106184,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.97382745763419,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.89457796475562,-4.83091437811374,1.35916496215604) q[0];
u3(1.53769719093844,-4.58298184020946,0.101501458141637) q[2];
u3(0.948504647776825,1.11269557172149,-1.28155135324234) q[2];
u3(0.355654624848898,-4.46477021179998,1.71301738775951) q[1];
cx q[1],q[2];
u1(3.39993812405877) q[2];
u3(-1.64681202904133,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.19655432865347,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.43881191064692,1.39477703368891,-1.03747454642297) q[2];
u3(2.50831066348214,-1.77761830242115,1.19897827157286) q[1];
u3(2.10475186974080,0.218132934121025,0.338337019183138) q[4];
u3(1.66349896717412,-1.67200517129666,-1.74756437284682) q[3];
cx q[3],q[4];
u1(3.15144135348311) q[4];
u3(-2.46323889514335,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.983489824987865,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.20930814827369,1.81536145009421,-1.89674991764702) q[4];
u3(0.968014267472283,-1.62323894727458,3.25219207537717) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
