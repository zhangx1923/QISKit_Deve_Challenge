OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(0.938020944544951,-0.720460940987001,2.72421649747871) q[11];
u3(0.824057479714065,-1.81281029433956,-1.89973117036750) q[1];
cx q[1],q[11];
u1(-0.742287279459004) q[11];
u3(-1.34981963709683,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.58248917585088,0.0,0.0) q[1];
cx q[1],q[11];
u3(2.22271473242464,3.25495217494620,0.153893298821773) q[11];
u3(2.14579474574358,-3.72492820382473,-2.22731789368079) q[1];
u3(1.57619361939428,0.474074117259611,1.67453683802970) q[6];
u3(1.72334804578619,-2.11827552912878,-1.34085445208794) q[3];
cx q[3],q[6];
u1(-0.127900941492047) q[6];
u3(-0.649245450733083,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.75047521766506,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.17612027784981,-2.11940523593284,0.0874515850808277) q[6];
u3(0.808648197912233,-2.41539558945830,-1.91163037019667) q[3];
u3(2.42753196011957,-2.10520820427709,-0.331350643904580) q[4];
u3(1.77533803916774,-4.95377086293457,-1.13957182940807) q[2];
cx q[2],q[4];
u1(1.61153168193653) q[4];
u3(0.246825734323816,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.821095220745843,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.90845963668159,-3.08492622881262,2.13268409110586) q[4];
u3(2.44722041143051,0.478459761453663,-5.53158115498920) q[2];
u3(1.77484738378681,2.38437102845830,-2.98479627321169) q[5];
u3(1.28195330651952,3.07460920609447,-3.00176936993931) q[8];
cx q[8],q[5];
u1(3.75006849869216) q[5];
u3(-1.58352387679398,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.32192389035441,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.35411194843867,-0.149168515262281,0.685968841638999) q[5];
u3(1.92698194862323,2.57179060352881,1.83726077502194) q[8];
u3(1.29913241179733,1.26018568621237,-0.738612000583007) q[9];
u3(0.216818530094673,-1.77566173936343,-0.799617140330556) q[12];
cx q[12],q[9];
u1(0.181331303710893) q[9];
u3(-1.52553729817291,0.0,0.0) q[12];
cx q[9],q[12];
u3(2.11325352761829,0.0,0.0) q[12];
cx q[12],q[9];
u3(0.625703231504917,-3.01624454686419,2.53066953367727) q[9];
u3(0.749915643133657,5.64305052266505,0.425203408353005) q[12];
u3(0.377824586608184,-0.307967236894987,0.532259028605534) q[0];
u3(1.04362463789928,-3.75529232917057,1.27830634316572) q[7];
cx q[7],q[0];
u1(0.334176520151548) q[0];
u3(-1.42538826786257,0.0,0.0) q[7];
cx q[0],q[7];
u3(3.22044536915301,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.03166028263662,1.21509187239005,-3.91703266984519) q[0];
u3(2.25900618711509,-3.32572102365850,-1.80674192134026) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12];
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
