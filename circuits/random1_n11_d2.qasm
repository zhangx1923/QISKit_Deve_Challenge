OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.27013063334427,-0.0518678276523850,0.466764432629316) q[9];
u3(2.06233625650496,-1.45786539476216,-1.69231615620356) q[3];
cx q[3],q[9];
u1(3.57353166152397) q[9];
u3(-1.40962126670175,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.29246661777140,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.14336747177771,-2.77053233435027,0.508594633906971) q[9];
u3(2.03504089972934,4.31861307112273,0.822494470477959) q[3];
u3(1.51109187793396,1.25231552762133,-2.91113431207844) q[6];
u3(0.879766147770252,-3.03800628262831,2.25660786903066) q[0];
cx q[0],q[6];
u1(3.35224654197542) q[6];
u3(-2.20169837652428,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.50200001429575,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.31204546979111,-0.762834184395847,1.32552196113985) q[6];
u3(1.76390831600991,-4.69559393307572,-0.636626732763352) q[0];
u3(2.27530888936593,0.116493079303111,2.02136246593894) q[8];
u3(1.88075322531497,-1.65646565635958,-2.75829051315240) q[7];
cx q[7],q[8];
u1(2.07480132226003) q[8];
u3(0.0620526441493650,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.954793787218014,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.835369990628771,0.809957640517190,-1.33696372331484) q[8];
u3(1.50036641066148,-2.47049140083478,-3.66039833076204) q[7];
u3(1.28240437881686,0.925906116717408,2.19795118061276) q[10];
u3(0.196296253786865,2.93140805274460,3.32371073808409) q[2];
cx q[2],q[10];
u1(3.90792885543596) q[10];
u3(-3.56772709904216,0.0,0.0) q[2];
cx q[10],q[2];
u3(-0.235578622708384,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.785421665475552,-1.47988385376741,3.80179616380972) q[10];
u3(1.51271146096481,-0.0956865270496288,-0.478318069178529) q[2];
u3(1.41344231318006,-2.83066781925971,-0.204291467939087) q[5];
u3(1.55264490290341,-3.29253696216552,-1.20369637624468) q[1];
cx q[1],q[5];
u1(-0.752624391326115) q[5];
u3(1.22617674982302,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.59023298123917,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.21769802157102,2.19449924498890,-1.37556498689011) q[5];
u3(0.148115555009014,0.263944236997421,1.04054852917231) q[1];
u3(1.61954137142852,-0.212541534773903,0.137905995015391) q[9];
u3(1.27573897668617,-2.75324435150983,-0.886869880335875) q[2];
cx q[2],q[9];
u1(1.26514161383765) q[9];
u3(-2.95596307525392,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.69782811088174,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.42188372051859,0.406210205674996,1.88281635825750) q[9];
u3(0.926195735186122,-4.20088974332312,0.112902398462898) q[2];
u3(0.738540587772388,1.76026892596055,-3.15350950688864) q[6];
u3(0.688044668542599,1.39922517272243,-2.64179514904198) q[4];
cx q[4],q[6];
u1(2.46237518890138) q[6];
u3(-1.72364029529841,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.14304345969839,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.683381042768464,1.16939958304705,-4.24584881219889) q[6];
u3(1.91899638122682,-0.853309839714703,-0.491659968868314) q[4];
u3(2.16427500174793,1.42720749960799,-1.12282798730718) q[1];
u3(1.77078378968859,1.41990497311431,-4.76784427916772) q[3];
cx q[3],q[1];
u1(2.17312174222000) q[1];
u3(-1.94930665039649,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.0294836434948933,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.31706170544093,-0.979414118854892,1.64170484088110) q[1];
u3(0.991024834152530,5.07199882279686,1.09618673153266) q[3];
u3(1.48168774271313,-0.492862046246435,-0.0992884331103528) q[10];
u3(1.48377552295267,-3.56091406709937,0.329148605719864) q[5];
cx q[5],q[10];
u1(0.755011816465853) q[10];
u3(-3.09568561748135,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.90025687636792,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.55567476425631,1.36120745495130,-2.66205610166476) q[10];
u3(2.03816123667959,1.96850941711898,-0.654355497355481) q[5];
u3(0.506666827872286,0.164995562574655,-1.71577515623714) q[7];
u3(1.32203674474730,0.483590854850242,-5.44798678829507) q[8];
cx q[8],q[7];
u1(3.78000381680466) q[7];
u3(-4.00963056836307,0.0,0.0) q[8];
cx q[7],q[8];
u3(-0.269328424273166,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.01301272034541,1.84397969964793,-3.01051298199930) q[7];
u3(0.744576227617177,-0.505916841416322,-5.73079979682806) q[8];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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
