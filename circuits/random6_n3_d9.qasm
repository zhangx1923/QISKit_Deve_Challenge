OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(1.00652887622803,0.195972122961161,-1.92708468480443) q[2];
u3(1.31429770609559,-3.42218178737084,2.33094924924285) q[1];
cx q[1],q[2];
u1(0.543711025705919) q[2];
u3(-0.120729876259384,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.82083694750251,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.36593172519113,1.31551217160909,0.413613314014002) q[2];
u3(0.329592183734356,-5.98322756205736,-0.0968998326691342) q[1];
u3(1.41038557990990,1.01191513828160,-3.13497468733693) q[1];
u3(0.889457498673405,-3.25358554712402,2.88359882890617) q[2];
cx q[2],q[1];
u1(2.44331022594893) q[1];
u3(-1.85668879627839,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.231425585004168,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.77198105303915,-0.0530045185637184,1.08944849252367) q[1];
u3(2.42967028910443,2.84612642255634,1.18619282111718) q[2];
u3(2.89438947512110,-1.66229094697294,2.45568244808000) q[1];
u3(2.01387285827835,-2.35682397644124,-0.756393797410866) q[0];
cx q[0],q[1];
u1(-0.198816694127555) q[1];
u3(0.724924483636939,0.0,0.0) q[0];
cx q[1],q[0];
u3(4.16686422097353,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.63049220264165,-0.170438273579485,2.85814006099759) q[1];
u3(2.77384547670481,3.58810355009096,1.20507685601851) q[0];
u3(0.815611646731479,0.799050597050143,-3.25881089840168) q[2];
u3(1.81433124212151,-3.37028693647192,2.64928585662135) q[0];
cx q[0],q[2];
u1(1.56148785731048) q[2];
u3(-0.784469376951113,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.341509919977737,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.29710809861342,-3.98315924233323,1.98432103738237) q[2];
u3(2.39410706418224,0.402956572363464,0.763141576863390) q[0];
u3(1.60063116997167,-2.32379120711991,3.92573182938451) q[2];
u3(2.28834777862551,1.39400550994824,-1.24449194992325) q[0];
cx q[0],q[2];
u1(1.06010039024542) q[2];
u3(-0.701671358124787,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.34056353829522,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.619544161487321,-0.0102780129314752,2.86374691308857) q[2];
u3(2.18533216298103,-3.36833221813437,2.69597303936239) q[0];
u3(0.451490379722823,-2.72829400921389,2.86716931184090) q[2];
u3(0.756030241913972,-2.86688089386632,-0.136464711085127) q[0];
cx q[0],q[2];
u1(3.11361602654876) q[2];
u3(-1.86621411721808,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.629135204315276,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.29065687088752,-4.08349938392709,2.01058659265598) q[2];
u3(0.330943472897037,-3.54040435167288,1.87891231651036) q[0];
u3(2.45659308567645,-0.317913692307838,-0.415234255824572) q[0];
u3(0.749572698857424,-0.222814363651198,-4.92186491220337) q[1];
cx q[1],q[0];
u1(1.37590336390020) q[0];
u3(-0.861869180796394,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.00598648323717,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.701502347325651,3.42803965845365,-1.34384865506219) q[0];
u3(0.994859937034623,-0.243446598293577,1.44746828415176) q[1];
u3(1.88929095658466,-1.67664958665151,4.48959204463302) q[0];
u3(0.274709151657154,-2.21426725644910,3.69069471317280) q[1];
cx q[1],q[0];
u1(3.17748924352876) q[0];
u3(-2.45728497575742,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.52779029211526,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.02677758523750,-0.759585060860318,2.40635291873298) q[0];
u3(2.38746571902201,0.703445687633830,1.54387858924118) q[1];
u3(1.37214713253376,-0.0528263757977057,1.67001002726568) q[1];
u3(1.05200852701222,-3.00963398516643,-1.95348798934203) q[2];
cx q[2],q[1];
u1(1.60190721039943) q[1];
u3(-3.45249260581737,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.711822925340502,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.53156513448785,0.180954700792059,0.660821240877076) q[1];
u3(1.01516155983523,-0.181377264268648,4.44906635844338) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
