OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.30171387987975,2.54635505001625,-1.43134249012680) q[6];
u3(2.96874592885837,3.80741316355593,-0.439770263818138) q[4];
cx q[4],q[6];
u1(-0.432105336871379) q[6];
u3(-1.67761209410975,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.652863893267410,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.21985663347541,0.801870924712242,1.83106406054053) q[6];
u3(0.302243413458666,-4.07064481310453,0.688680339843350) q[4];
u3(2.38163868641900,-2.72562743983900,0.949501082361906) q[7];
u3(1.87313287896452,-3.71021400027302,-1.95630459066046) q[9];
cx q[9],q[7];
u1(0.845421328122906) q[7];
u3(-0.246781046112495,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.79657068296615,0.0,0.0) q[9];
cx q[9],q[7];
u3(2.55953741558094,-1.33038036444910,-0.134971708373052) q[7];
u3(2.89991241752065,-1.83622382512232,-0.663821895153502) q[9];
u3(2.91069517245241,1.56948688805667,1.50955193439843) q[10];
u3(1.52622439485364,-1.97620087178180,-2.03001262088313) q[11];
cx q[11],q[10];
u1(3.77836619816502) q[10];
u3(-1.25569126233393,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.54935004436904,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.74220471081931,2.67121538500828,-2.73521014719881) q[10];
u3(1.14617542034480,1.44111507315363,1.50358782580032) q[11];
u3(2.10981853223561,-0.664036788627053,-2.13976532323381) q[2];
u3(1.04741847573434,-3.84348364379306,1.25952750267983) q[0];
cx q[0],q[2];
u1(0.663604015701598) q[2];
u3(-1.22408463329778,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.93627923616107,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.301593291310189,-0.719531857813036,4.06082417338635) q[2];
u3(0.0969958628280253,-2.20314308637890,1.07036112079655) q[0];
u3(2.52759630601989,2.22174005380352,-0.583201050892862) q[5];
u3(1.89522262844624,5.08142336989446,0.673752860681410) q[8];
cx q[8],q[5];
u1(1.92961816117643) q[5];
u3(0.321612581129166,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.11733505940152,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.49266804672315,-1.55239110076755,0.657258682777709) q[5];
u3(1.51009541966712,2.40870683931457,-0.337192350710587) q[8];
u3(0.924872467043372,-0.400812873316924,-1.80172295334365) q[1];
u3(1.81332087923206,1.55471276179658,-3.90506092953746) q[3];
cx q[3],q[1];
u1(-0.440395269362116) q[1];
u3(-1.64314446170335,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.01412138891958,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.745374570040956,1.57915967614341,-1.06201446312179) q[1];
u3(1.87142794316060,-3.71027232892610,-2.01313234008140) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
