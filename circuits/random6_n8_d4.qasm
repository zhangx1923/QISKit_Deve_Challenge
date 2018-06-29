OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.68079709459079,1.31158604832202,-2.88855346016228) q[2];
u3(1.08560316488430,-2.45507168647757,3.00632547826119) q[5];
cx q[5],q[2];
u1(-0.241420273719716) q[2];
u3(-2.33218526592442,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.30262469743746,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.28783763880252,0.218722919818200,1.37389271179467) q[2];
u3(0.598759887765437,-0.455453534366718,2.69192598362501) q[5];
u3(2.59577455103970,0.392366520460433,-1.86337015517946) q[7];
u3(2.68956378158902,4.52112699342817,-0.0490875474169292) q[3];
cx q[3],q[7];
u1(1.77252892156634) q[7];
u3(0.159544881523491,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.672279714067999,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.430345849821622,2.64350961238440,-0.397606049105870) q[7];
u3(0.613856098293385,2.70881947079196,-0.896873427040023) q[3];
u3(0.455092496761643,-1.07653944357868,1.17880472819551) q[4];
u3(0.967543374521160,-0.193930250619400,-1.21706677456819) q[6];
cx q[6],q[4];
u1(1.69182064051855) q[4];
u3(-3.36682271089761,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.32228768965590,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.05879867803947,3.79233714072661,-2.13311251973336) q[4];
u3(2.31141252600167,-2.60838722111538,2.96908260028185) q[6];
u3(1.11338529446326,-1.34546005465576,1.12667950497172) q[1];
u3(0.638397693863644,-1.59038649649799,-0.181669550128022) q[0];
cx q[0],q[1];
u1(1.18951811511755) q[1];
u3(-0.0382731782681296,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.72353046707254,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.12479544673214,4.72160301271998,-1.20304195187058) q[1];
u3(0.579932464532254,-2.56256642837145,-0.0933913917735432) q[0];
u3(1.23079802098890,0.604727630654748,1.28350576770934) q[0];
u3(0.995330347879597,-1.05801826998773,-2.23175254093339) q[7];
cx q[7],q[0];
u1(1.66892523204587) q[0];
u3(-0.877019759721651,0.0,0.0) q[7];
cx q[0],q[7];
u3(3.09549316730089,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.78878447292242,-0.844482257768993,0.348925658619877) q[0];
u3(0.953238100936780,-0.223626458788357,-3.09204589208718) q[7];
u3(0.967836092632520,1.06045758888385,-2.98772227512347) q[5];
u3(1.77231275968189,-3.33258467330563,2.92709708790064) q[4];
cx q[4],q[5];
u1(1.71088896033934) q[5];
u3(-2.91523974989030,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.873423393272098,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.712140472317736,-3.17027087058790,0.322660283782300) q[5];
u3(2.66891494509281,-2.16383435114628,-0.934702193998421) q[4];
u3(1.62384604185506,1.58973362128200,-0.707673896164690) q[2];
u3(0.789425947679228,0.419920424341418,-3.17135479872554) q[3];
cx q[3],q[2];
u1(0.474421590166367) q[2];
u3(-0.142470978926836,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.74406849410159,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.74291120986239,2.29071272431957,-2.08278590011680) q[2];
u3(1.10149181723818,2.46277404329673,2.90464839535662) q[3];
u3(1.28414383564978,2.94701389643672,-1.40047224101377) q[6];
u3(2.46029426036717,2.11949490565252,-1.04628265525112) q[1];
cx q[1],q[6];
u1(1.18417477579604) q[6];
u3(-0.865267162125119,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.301600975042851,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.60335287501741,-1.32032483205490,0.450376246181895) q[6];
u3(2.16042310276751,3.96293422172858,-0.666370082150398) q[1];
u3(1.54512850682276,1.42195876912029,-2.00562693147664) q[2];
u3(0.297445635018052,-1.94765427113405,1.77719550080166) q[6];
cx q[6],q[2];
u1(2.20910085898979) q[2];
u3(-2.87372978252007,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.977942779976174,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.20809655123256,-0.399967306530207,2.44222481003456) q[2];
u3(0.832980243926968,1.79724316829183,-3.44535816809982) q[6];
u3(2.39318719657671,-4.45823978721216,1.51350935305490) q[0];
u3(0.713720752212400,-1.88345463756289,4.23231015926978) q[5];
cx q[5],q[0];
u1(3.40053314880974) q[0];
u3(-4.41607521870918,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.440718367663316,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.973790482814322,-1.03418783908250,3.73877902466310) q[0];
u3(1.21380791444521,1.39283566372437,1.66833533645352) q[5];
u3(1.16815658670145,0.0833830919456235,1.62069777561722) q[7];
u3(1.43565747528952,-2.22474155928669,-1.35983215997017) q[1];
cx q[1],q[7];
u1(1.68081331649063) q[7];
u3(0.194965607244989,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.13378302472118,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.32909484056942,-0.927683471516269,-0.642620426784745) q[7];
u3(1.72235901211486,4.99101411065771,-0.407986983462660) q[1];
u3(2.24715559447712,-0.809935941782867,-0.921616919097996) q[3];
u3(1.23183567371587,0.898828189828990,-4.34036573862921) q[4];
cx q[4],q[3];
u1(0.767664862250329) q[3];
u3(-1.41739954348674,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.73474104695353,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.88756087506943,2.18085096327021,-0.898916060694357) q[3];
u3(2.25044400464551,-0.252231209687938,5.60130463137160) q[4];
u3(1.49355008698772,-0.269729787167411,0.940758362251968) q[3];
u3(1.40177872239859,-1.78516553931620,-2.19493231673962) q[1];
cx q[1],q[3];
u1(2.94849197295250) q[3];
u3(-2.01111315856164,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.58541923335327,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.43274007991043,1.86058681814215,-4.32131381186379) q[3];
u3(1.00591817408157,5.44832277930704,-0.342003524800452) q[1];
u3(1.52226596401128,1.74748252836691,-2.45095333045009) q[6];
u3(0.824051758450150,2.42538178587485,-3.51321392434612) q[0];
cx q[0],q[6];
u1(0.268845727147574) q[6];
u3(-0.752576371798398,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.76038985944333,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.90631481264207,-2.09266273825064,1.60728719855188) q[6];
u3(0.708213002699928,-1.75957420517497,-4.06597201957520) q[0];
u3(1.00250905418854,-0.478380658240940,-0.970509496239132) q[2];
u3(1.15103193842599,-4.08781944856478,1.02305585987232) q[7];
cx q[7],q[2];
u1(0.985862814894839) q[2];
u3(-0.0962609420279177,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.83539844357143,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.382030736338539,0.976194298450749,-0.499989863809298) q[2];
u3(1.98714859234246,2.04021148128826,3.23938002808977) q[7];
u3(1.61346594537041,-1.75239566865731,-0.215829787371997) q[5];
u3(1.38890738671469,-2.00864893710293,0.198552908805067) q[4];
cx q[4],q[5];
u1(-0.977408079014586) q[5];
u3(0.136709308654917,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.67947517868277,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.26796389797667,-1.37345785076276,-0.760295895980506) q[5];
u3(1.96918827958722,-0.146977633598541,1.26373572361731) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];