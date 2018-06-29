OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.87626681056000,3.00761673884805,-0.319957486632723) q[12];
u3(1.80760087975362,2.35159923764836,-1.66164352693113) q[8];
cx q[8],q[12];
u1(3.40627219148793) q[12];
u3(-4.01100522118001,0.0,0.0) q[8];
cx q[12],q[8];
u3(-0.701943468234999,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.06745198812063,-2.89293420175983,-0.382803766277964) q[12];
u3(2.41610746748360,4.30898520873371,1.23475325417945) q[8];
u3(1.73124299367557,-2.18744470086250,-0.0693313818879824) q[3];
u3(1.11385511287869,-3.22491790794233,0.773800805224225) q[6];
cx q[6],q[3];
u1(2.88764573867106) q[3];
u3(-1.87739429282459,0.0,0.0) q[6];
cx q[3],q[6];
u3(-0.00652491467868899,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.45054699220030,-2.94820887990273,2.45153312771252) q[3];
u3(2.60742519672328,0.421423325876013,5.62437319290452) q[6];
u3(2.83286196031190,-2.98908905263899,1.30351090142513) q[10];
u3(2.76355638333343,0.741293546459219,2.23561619258686) q[5];
cx q[5],q[10];
u1(3.25820453849811) q[10];
u3(-0.737716795783929,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.59992041796585,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.96088604183441,-0.797627505487170,3.76314360499231) q[10];
u3(1.66871136345171,1.40456488877366,0.123382531825139) q[5];
u3(2.09820752936817,2.37212518701766,-1.13007456255389) q[9];
u3(1.69044990660253,1.23148495147690,-1.43004183375890) q[1];
cx q[1],q[9];
u1(1.64868722907434) q[9];
u3(-0.758831350943295,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.47480108717200,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.730200708092574,-4.22228878368750,1.52749185787449) q[9];
u3(2.17855010412649,4.02693988614913,0.455250775628674) q[1];
u3(2.77089380147784,2.70026654269422,-0.260683881007926) q[7];
u3(2.70882830045467,4.17815711944319,0.273511191118051) q[2];
cx q[2],q[7];
u1(2.08143865651699) q[7];
u3(-3.42185457104035,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.884669425775072,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.23399340881497,2.14242348307000,-1.81735438533607) q[7];
u3(0.987008192935879,-1.41518971454217,-2.42748676448651) q[2];
u3(1.88188246140053,2.41020717108666,0.686336768928732) q[0];
u3(2.07540536311912,0.375922070318620,-1.84237524035034) q[11];
cx q[11],q[0];
u1(-0.0404044351649606) q[0];
u3(1.17608934478020,0.0,0.0) q[11];
cx q[0],q[11];
u3(3.40982111922347,0.0,0.0) q[11];
cx q[11],q[0];
u3(0.774066373063484,0.0579418660638360,2.16408471442195) q[0];
u3(1.80334005659322,-1.70741941060779,-3.97097664532440) q[11];
u3(0.551229700863175,-0.283671135712446,-1.38238733748619) q[5];
u3(1.54209445293407,1.09527103952202,-4.62480394820396) q[11];
cx q[11],q[5];
u1(3.39293193196176) q[5];
u3(-0.721508848251719,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.51642187204874,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.56343532468978,-0.420536349079787,3.26148553322466) q[5];
u3(1.90812609089844,-3.50134822419002,-1.78932657885413) q[11];
u3(2.91644096332565,-1.74492477705008,1.53852286910721) q[0];
u3(1.61038522597588,1.11172891108012,2.45257469663158) q[9];
cx q[9],q[0];
u1(3.80669691680715) q[0];
u3(-4.28108145348960,0.0,0.0) q[9];
cx q[0],q[9];
u3(-1.00349300250015,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.97581965656532,1.09454098849552,-1.12220783238271) q[0];
u3(0.284275648931984,1.77087385061979,1.81127889150769) q[9];
u3(2.68199201729403,2.01393162276679,0.120392587100086) q[6];
u3(2.46833509725114,1.05193863790637,-3.59140472966787) q[3];
cx q[3],q[6];
u1(1.40647525816941) q[6];
u3(-0.287529209327323,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.49587499805229,0.0,0.0) q[3];
cx q[3],q[6];
u3(3.01622546884266,4.00493843690007,-1.66481265846363) q[6];
u3(2.07979663352230,1.91169422804876,-0.826001274784135) q[3];
u3(1.98514446306224,-0.312375933215710,-1.55341778233438) q[10];
u3(0.969825270141483,-4.18322672986079,0.652324594543594) q[2];
cx q[2],q[10];
u1(3.88247412825518) q[10];
u3(-4.31178426084817,0.0,0.0) q[2];
cx q[10],q[2];
u3(-0.485825493419676,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.40973320049053,-0.626853142250888,2.39716757203724) q[10];
u3(2.24234456004184,-1.64530696560367,4.05324204458206) q[2];
u3(2.57894563372724,0.747395980451765,-2.46695969484531) q[12];
u3(2.01167059240399,5.02346274540386,1.18446229537018) q[7];
cx q[7],q[12];
u1(4.01667717596644) q[12];
u3(-3.85121994647976,0.0,0.0) q[7];
cx q[12],q[7];
u3(-0.292486238600456,0.0,0.0) q[7];
cx q[7],q[12];
u3(2.60303768806576,-3.30849976240652,0.357605526922833) q[12];
u3(1.06752137240358,-1.69435788871034,-4.43515968967472) q[7];
u3(2.04705073992890,-2.85742889867715,0.294270313765889) q[4];
u3(2.16432663421882,-3.16431812941183,-0.654178696714923) q[8];
cx q[8],q[4];
u1(1.13470071959296) q[4];
u3(-3.43161899730776,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.17923698068927,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.58846336525754,-0.771488641186888,1.17173196569920) q[4];
u3(1.97544914274351,-2.79535690715970,3.23791958338004) q[8];
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