OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(2.45192130614668,3.43959172158169,-0.467270815404653) q[2];
u3(1.69779862869540,3.18414541721938,0.376972722649992) q[3];
cx q[3],q[2];
u1(2.36833190738719) q[2];
u3(-1.70980983290780,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.604682584304978,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.50306226590108,-1.64335390531354,1.40656913842849) q[2];
u3(1.42849108129318,-4.76176282077742,0.363261130368941) q[3];
u3(2.53305679932082,0.699036100256254,1.07105647300075) q[0];
u3(0.581140811984831,-4.34576231438089,-1.38545288893677) q[6];
cx q[6],q[0];
u1(0.653505788459284) q[0];
u3(-0.305726279685848,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.02395426396422,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.18922451221310,0.439636996136545,-2.80615165101330) q[0];
u3(1.54383225143054,0.838942839803274,3.26550083600066) q[6];
u3(1.12103703461405,1.41016220282049,-0.0738707689504825) q[5];
u3(1.32254663129781,1.14776770540337,-2.25253271267807) q[4];
cx q[4],q[5];
u1(3.98131699113058) q[5];
u3(-3.65902030321069,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.707597230936066,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.588802775128922,0.966476059760813,0.497088919170027) q[5];
u3(1.99210110041674,-3.64633368196878,0.673336423668068) q[4];
u3(1.87340781280737,1.48531958503060,-0.0597077688944327) q[6];
u3(1.98441412405333,0.385246516103168,-2.03349958705530) q[1];
cx q[1],q[6];
u1(-0.624183239194374) q[6];
u3(-1.98950676331821,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.42323992972089,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.21756006319148,2.51749741637484,-0.700882629175171) q[6];
u3(1.38573083640013,1.23616474325825,-3.81595159604452) q[1];
u3(1.46204707359420,2.33337608386244,-3.54601006436612) q[4];
u3(2.01638661462106,3.62216000443859,-2.52416351143705) q[2];
cx q[2],q[4];
u1(1.56931805268316) q[4];
u3(0.632703053544458,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.707279611056231,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.52960456783138,-4.14481951074135,1.02668154036328) q[4];
u3(1.49634648270147,1.69856720882274,-2.05063451971499) q[2];
u3(1.75795977684542,-1.29211128409005,-0.338038887237751) q[5];
u3(1.82487905457078,-2.64481731317778,-0.0621913909428389) q[3];
cx q[3],q[5];
u1(3.40471798863485) q[5];
u3(-3.29308650536400,0.0,0.0) q[3];
cx q[5],q[3];
u3(-1.36327319101496,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.819884502744605,2.50018089460747,-3.43881044901172) q[5];
u3(2.77874437920549,0.935118681225876,-2.29989929573725) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];