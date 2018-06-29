OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.33935387655336,2.23329302625165,-1.34708753995854) q[8];
u3(0.438337069876294,1.14272539196333,-1.76139131911478) q[6];
cx q[6],q[8];
u1(2.48888069205124) q[8];
u3(0.200271113796725,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.33334937483314,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.15337555592487,-1.89192106349766,3.82997111411512) q[8];
u3(1.25280765308376,-0.507177032825745,3.68942152303397) q[6];
u3(1.60175185059849,2.78456527382563,-2.91096856577185) q[0];
u3(2.18637012571693,-3.27478441975908,2.55168599452979) q[9];
cx q[9],q[0];
u1(1.43530909000692) q[0];
u3(-0.286110837189995,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.94043349580732,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.26184340024784,-0.923040103640900,0.125716723129934) q[0];
u3(2.73629121754250,1.25013525238063,-4.81099453609430) q[9];
u3(0.926207247822725,3.38329828864185,-1.41517242597936) q[10];
u3(1.33018676994994,1.23385640772243,-0.384631496670319) q[5];
cx q[5],q[10];
u1(1.77121358728676) q[10];
u3(0.183010798974502,0.0,0.0) q[5];
cx q[10],q[5];
u3(0.489004910238505,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.20714627481743,-0.200514519340473,1.51121188314946) q[10];
u3(1.93689467245460,-0.0403877568410645,-0.927900237558089) q[5];
u3(1.45739462028565,-1.08117531992388,0.864335338243327) q[7];
u3(0.117853800645395,-1.56539945718307,0.618984517887801) q[3];
cx q[3],q[7];
u1(3.09419532420196) q[7];
u3(-1.03608005681729,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.49595228536827,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.887136091051534,0.323845074562873,-1.36792158387243) q[7];
u3(0.722998142109061,0.845740418964613,-5.40929501262417) q[3];
u3(0.245621328159527,-0.476132127465095,-1.34933850166701) q[2];
u3(1.77833873502296,0.493881258829441,-4.78869316147225) q[1];
cx q[1],q[2];
u1(2.73238482719436) q[2];
u3(-1.98203141252282,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.25072275484370,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.873359963931430,-1.66954978402820,2.47500295816488) q[2];
u3(1.97981907936755,0.281081823863332,5.03209416723828) q[1];
u3(2.26938536033037,0.160038412841084,2.11152673351591) q[6];
u3(1.16762350548451,3.39272375064663,2.78282813129682) q[5];
cx q[5],q[6];
u1(2.43472424482332) q[6];
u3(-2.93300949043804,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.54066352405977,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.00683110605434,-4.96705171582474,0.773472100161675) q[6];
u3(1.19481635078769,1.68352177775851,1.45574150469487) q[5];
u3(1.32722039133626,1.12753682939376,0.620691193461973) q[1];
u3(1.22277972798920,-1.54080537124049,-1.15302264576515) q[7];
cx q[7],q[1];
u1(3.30843585854132) q[1];
u3(-0.778646083584576,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.56725964506569,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.73512146412022,-1.23474822837218,-1.10966904041506) q[1];
u3(2.73101468320966,-1.30475185661412,3.17667711418649) q[7];
u3(2.35295756952942,0.0397171033397267,0.779342246393388) q[0];
u3(0.161852414597593,-1.52895909744758,-3.42343645189733) q[4];
cx q[4],q[0];
u1(4.05112813064846) q[0];
u3(-3.57891022903114,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.706787291452792,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.49519817759972,2.77520310277521,-2.37840805667262) q[0];
u3(0.740173596841953,1.12128003613400,3.28423414574800) q[4];
u3(2.29737143560810,0.0876178265314318,1.71053384760948) q[2];
u3(2.31782644296978,-1.29511818335119,-0.614281624668125) q[8];
cx q[8],q[2];
u1(2.55721070756767) q[2];
u3(-1.94758925468094,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.286702634618270,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.47134656085037,1.79537173753224,-0.0609282746692151) q[2];
u3(0.921311891889808,-3.08386626028452,0.434802121304066) q[8];
u3(1.52275632109612,0.415504398700433,-2.05682938483098) q[3];
u3(1.36804217607389,2.80405833025391,-3.43765138443037) q[9];
cx q[9],q[3];
u1(3.43605650430869) q[3];
u3(-4.65006385695905,0.0,0.0) q[9];
cx q[3],q[9];
u3(-0.351918407804622,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.96655316892133,-1.14202283737615,1.25789205110406) q[3];
u3(1.58735256591077,-1.52168093048019,4.13768019873657) q[9];
u3(1.02789867249371,1.14280305322389,-1.94911587526539) q[4];
u3(0.466248616987781,-3.38057918803353,2.36719929424129) q[9];
cx q[9],q[4];
u1(0.305545146193762) q[4];
u3(-1.54901150950192,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.08458799974610,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.854929850745378,0.835535140248376,1.77825114873907) q[4];
u3(1.76092280074044,-1.26318204645030,-1.65865491424837) q[9];
u3(1.13642690634854,0.959251294749365,-1.66624775105014) q[0];
u3(0.413457320679364,-1.59279240268807,-0.215764060018343) q[7];
cx q[7],q[0];
u1(1.08737493800433) q[0];
u3(-0.0365070547102577,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.67639897430918,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.12753804552337,2.65633723134435,0.224705360424395) q[0];
u3(0.788886850174391,-2.30293017777920,-3.97832448016679) q[7];
u3(2.90652872694853,2.64907037084509,0.143070917769903) q[3];
u3(2.17720143874965,0.989866880655640,-3.65179802137540) q[2];
cx q[2],q[3];
u1(1.37552909372895) q[3];
u3(-0.326724046734948,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.35540223794761,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.702858121005632,-2.27520378324145,1.37050375508870) q[3];
u3(2.18158885814687,0.723278141916603,-5.08663419527162) q[2];
u3(2.09615625439932,0.0563797600188065,-0.375842636470558) q[10];
u3(0.487712579387751,-1.17616342957713,-3.72049373700532) q[5];
cx q[5],q[10];
u1(0.991226817485084) q[10];
u3(-0.139009246214215,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.67475729731969,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.437921372503715,-1.03193277296356,2.24045457665499) q[10];
u3(0.994891937100375,-4.40454679548336,-0.0365905237518596) q[5];
u3(1.60714609962293,0.801477343017417,-2.85074775388195) q[1];
u3(1.41565180217542,-2.78761424600611,2.91451452838811) q[8];
cx q[8],q[1];
u1(1.73989398924381) q[1];
u3(-0.599054478547310,0.0,0.0) q[8];
cx q[1],q[8];
u3(0.0888499749180409,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.18119323746724,2.30069386392981,-0.0448199066562263) q[1];
u3(0.906180396965204,-0.594103057311124,3.72966281656741) q[8];
u3(1.72959688112409,1.52301975785346,-3.01197019212099) q[8];
u3(1.18673915704540,-2.75737078174750,3.43317257531034) q[0];
cx q[0],q[8];
u1(2.90380054055148) q[8];
u3(-1.04150488780665,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.59622629713766,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.79503080183390,0.909256637189472,1.92166007433178) q[8];
u3(1.59963493781454,0.535349829158143,1.82426755929994) q[0];
u3(2.25691477584162,2.46440065498760,-1.32472088894779) q[9];
u3(1.64658500433098,0.828656043089382,-1.41375960460125) q[1];
cx q[1],q[9];
u1(2.56290878383169) q[9];
u3(-2.20462566166986,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.37235076069162,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.86220707394450,0.501892665708683,-1.35339922813514) q[9];
u3(1.87913554386602,-2.40095402494277,-2.70556772187047) q[1];
u3(1.90061479588534,1.89101482066357,-0.922039746000962) q[3];
u3(2.72650141204420,2.87945400446194,-1.05539963487219) q[10];
cx q[10],q[3];
u1(0.496814119945661) q[3];
u3(-1.50993320963971,0.0,0.0) q[10];
cx q[3],q[10];
u3(-0.268858976060288,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.408103761457474,1.70523584077947,-1.71338088487317) q[3];
u3(2.21824735306345,0.0881865820892838,1.11353999240259) q[10];
u3(2.00163505001000,-1.71951847671476,-0.240010724148730) q[6];
u3(1.91392046104829,-3.88876609509442,1.02127962366515) q[5];
cx q[5],q[6];
u1(1.52748404341018) q[6];
u3(-0.842181651661086,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.10880922972424,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.965556153321428,-4.33936166442296,1.87254238765981) q[6];
u3(2.40592810053772,2.01592443888321,3.05561933834925) q[5];
u3(2.35949203748165,-0.912049133118548,2.75663275158291) q[4];
u3(2.24220035257588,-2.44815173858388,-1.79782082530030) q[7];
cx q[7],q[4];
u1(-0.515294049085468) q[4];
u3(1.07605587622407,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.42400057869059,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.586953022998031,-2.92953887623247,2.70908661285609) q[4];
u3(2.08482656200791,0.835890309367451,5.18976468837907) q[7];
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
