OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(1.88110388022359,-2.45144110470164,0.602520985392576) q[1];
u3(1.50810594439386,-2.58598756910212,0.448212288634340) q[0];
cx q[0],q[1];
u1(0.760985355296544) q[1];
u3(-3.48042408428743,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.94019547092156,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.00976362409719,1.89065066129232,1.28411879900343) q[1];
u3(2.24742820527797,-0.0600169726127839,3.35483387712811) q[0];
u3(1.33480132021248,0.681424934725610,1.97196488682040) q[2];
u3(0.960784671929735,-0.832026502363888,-1.52805745111252) q[1];
cx q[1],q[2];
u1(-0.256300655145180) q[2];
u3(-2.09978765272207,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.66551764069175,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.52819350460716,-0.599035721281697,0.938660513856894) q[2];
u3(0.179975960185063,2.97361489386617,2.11373506989326) q[1];
u3(1.83731604494566,2.21473541867515,0.105075411869921) q[0];
u3(1.70185086596781,0.564156494926423,-2.32856776178109) q[2];
cx q[2],q[0];
u1(2.53614424655868) q[0];
u3(-1.93362032460075,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.04266187008890,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.28305169577039,-2.23260747374848,1.81832274905665) q[0];
u3(1.90361940574530,-0.907027068175402,-5.31843539504318) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
