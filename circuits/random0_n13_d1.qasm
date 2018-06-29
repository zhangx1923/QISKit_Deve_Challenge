OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.46741584245397,2.75535052478881,-1.52234781552504) q[5];
u3(0.456705853537602,0.707711066662244,-1.14536849618533) q[8];
cx q[8],q[5];
u1(1.88580801763249) q[5];
u3(-2.53545587044942,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.15389238593226,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.33153183374995,0.222204087298083,-0.0140538492738960) q[5];
u3(2.63979367584413,4.73445692819953,-0.703910470612386) q[8];
u3(1.98383211677832,-1.49989928517876,3.78363956913186) q[3];
u3(1.42663704070420,1.97539063319865,1.92402614368689) q[11];
cx q[11],q[3];
u1(0.753196086064123) q[3];
u3(-3.23078652511246,0.0,0.0) q[11];
cx q[3],q[11];
u3(1.78108172537202,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.79961668379733,1.57420903239179,-1.48947580975417) q[3];
u3(0.767475904104038,1.71644504493290,-0.654424883120197) q[11];
u3(2.18719747084908,0.888097059612952,-2.41139552034210) q[0];
u3(2.49965645947825,-0.0949894661766111,-5.17937914835329) q[2];
cx q[2],q[0];
u1(3.73874519323784) q[0];
u3(-1.39185890778044,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.25314760104391,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.471266024859673,-2.36118414938063,0.918684457395112) q[0];
u3(1.76297480921986,1.85669893983428,1.38941087520211) q[2];
u3(1.21645456428753,2.72124941388697,-2.86537959759940) q[12];
u3(1.70937944925581,-3.40581365337752,2.48055578449038) q[1];
cx q[1],q[12];
u1(0.313701438984881) q[12];
u3(-1.69760091756272,0.0,0.0) q[1];
cx q[12],q[1];
u3(2.40291649979221,0.0,0.0) q[1];
cx q[1],q[12];
u3(2.26177806490720,1.01158341672352,1.82837159580528) q[12];
u3(1.94752491960554,0.872010643756365,-3.63132717011266) q[1];
u3(1.57674314628588,-0.0458077109034981,-1.63218417925575) q[9];
u3(0.686930180168901,0.886983908519490,-3.54354973371269) q[4];
cx q[4],q[9];
u1(4.26112745020175) q[9];
u3(-3.47689066947240,0.0,0.0) q[4];
cx q[9],q[4];
u3(-0.488811279361806,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.861141613349208,2.88659480624487,-2.76933646434134) q[9];
u3(1.61298666100602,2.62630782987330,0.107104476893588) q[4];
u3(1.05396918823349,-1.06248103146620,-1.00818718839262) q[6];
u3(0.982312941363160,-4.06610335921112,0.304372444556277) q[10];
cx q[10],q[6];
u1(0.291951502998086) q[6];
u3(-1.46148608457044,0.0,0.0) q[10];
cx q[6],q[10];
u3(-0.0752672827086656,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.13352462526055,2.09048855385147,-1.73533547388399) q[6];
u3(1.32737569381673,-2.23278575735777,3.06960284618716) q[10];
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
