OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.59084476372687,-0.280203731767003,0.834395271192675) q[9];
u3(1.75472423174388,-2.69222972181257,-1.28542953671659) q[4];
cx q[4],q[9];
u1(3.34650213519729) q[9];
u3(-1.70791246548799,0.0,0.0) q[4];
cx q[9],q[4];
u3(2.24200479069178,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.853658381655277,-1.85141928271306,4.30248567159943) q[9];
u3(0.537111174801310,-0.906606989825168,-1.16523210130458) q[4];
u3(1.19425445206086,0.418962727487792,-2.63976343401178) q[3];
u3(1.54657278374796,-2.78338867744576,3.06077481571788) q[6];
cx q[6],q[3];
u1(0.882031760360635) q[3];
u3(-1.39127494418177,0.0,0.0) q[6];
cx q[3],q[6];
u3(-0.118470732699257,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.33746056342012,0.440311723846317,-1.85543350307603) q[3];
u3(0.751939604723901,2.75363755589068,-2.20871178987434) q[6];
u3(2.58673176005909,-1.98283022875207,1.79419409723306) q[2];
u3(2.11361582376258,-3.31729375545326,-2.33956922164962) q[1];
cx q[1],q[2];
u1(3.80657825284793) q[2];
u3(-3.51403918154203,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.776582865569660,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.82893175486770,1.79840644498508,-2.76055277746658) q[2];
u3(1.65843740117150,-3.50450453119147,-0.408013318204311) q[1];
u3(1.60123280801689,1.67345310846894,-2.76645206675796) q[7];
u3(0.180472529390153,1.27975512707986,-2.46938498668532) q[5];
cx q[5],q[7];
u1(1.64856369868618) q[7];
u3(-3.10820758518957,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.40301925981143,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.92741641461830,1.82265841046919,-2.80991257046578) q[7];
u3(1.31590384528130,2.49900986619333,-0.314382521992001) q[5];
u3(0.835265707237791,-0.968673513108360,0.467051410737174) q[10];
u3(1.23861017054018,-2.35643910360004,0.800819694701367) q[8];
cx q[8],q[10];
u1(1.73187502130625) q[10];
u3(0.275738122666796,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.975854444182104,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.07860669758891,-2.07507177696132,0.563795017713161) q[10];
u3(2.99410764755339,-0.813257846069432,0.291964428524230) q[8];
u3(1.51041779544066,-0.776158740367963,-0.642799961622468) q[5];
u3(1.60445741616795,-2.63798336130121,0.293161408738113) q[7];
cx q[7],q[5];
u1(1.29310535888548) q[5];
u3(-0.821136411206452,0.0,0.0) q[7];
cx q[5],q[7];
u3(-0.493606648047993,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.94652045213843,-0.241614736614909,-0.762052925857586) q[5];
u3(1.69538515627632,2.09069840217366,0.771400903664313) q[7];
u3(1.16377855351556,0.300255795804370,-1.45476976794417) q[10];
u3(1.89538689838037,1.95444081410935,-4.00707937976816) q[6];
cx q[6],q[10];
u1(3.53530919689484) q[10];
u3(-1.38438844411426,0.0,0.0) q[6];
cx q[10],q[6];
u3(2.49303542408907,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.557577317815602,3.88269397731713,-2.32743081505358) q[10];
u3(2.01442277357723,-1.02279612342533,-4.28847538507003) q[6];
u3(1.78592928656594,0.108364817520623,-2.05546758653921) q[8];
u3(1.31536561290910,-4.12936963921836,1.53499848243450) q[4];
cx q[4],q[8];
u1(3.30164314411102) q[8];
u3(-1.42417887544687,0.0,0.0) q[4];
cx q[8],q[4];
u3(2.48636694464363,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.62665977744408,-0.232723782638111,-1.64929028429536) q[8];
u3(2.16016495103999,-3.74025058091772,0.546921009271124) q[4];
u3(0.921883718974432,-1.69605091223927,0.207946408895140) q[3];
u3(0.377614730362408,-2.13869127499370,0.603262015641396) q[2];
cx q[2],q[3];
u1(1.13472105362525) q[3];
u3(-3.59867362469134,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.76022354031363,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.78723381329701,1.29247035442785,-1.29350809866722) q[3];
u3(1.89917484316131,3.78254964410919,1.25735288073703) q[2];
u3(0.526036525129153,1.31687531415155,-3.75595132569609) q[1];
u3(1.80952266426704,-1.40783607096068,4.46809022794065) q[0];
cx q[0],q[1];
u1(1.55834914445826) q[1];
u3(0.147613018530045,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.928565841782323,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.79777961168974,-3.56771581773702,0.631752829943451) q[1];
u3(2.15896167120378,3.83544202911685,-1.93367155687233) q[0];
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