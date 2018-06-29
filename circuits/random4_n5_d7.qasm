OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.79664197554201,-2.07126147741876,2.62360839424208) q[1];
u3(1.75127931696290,0.0886806527273640,-1.32653083674631) q[0];
cx q[0],q[1];
u1(0.644570257228448) q[1];
u3(-0.173170949516622,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.22548883460929,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.56246570294926,-0.620250608397910,-2.07880456637108) q[1];
u3(0.995673601944774,-4.10811479600884,-0.618394262872261) q[0];
u3(1.53296829618389,-1.95245494085569,-0.418693899183485) q[2];
u3(1.81567799295824,-3.89467306638508,0.856177037180229) q[3];
cx q[3],q[2];
u1(1.66627137025389) q[2];
u3(-2.40934768985589,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.408446595063424,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.04435043573174,-0.198333731705442,-3.81113922700703) q[2];
u3(1.64786548243717,-2.05250243727131,-2.87835160037495) q[3];
u3(1.22510633546240,0.199139550357278,1.15002965907349) q[4];
u3(1.05306135690883,-2.08736648817450,-1.40037186187244) q[3];
cx q[3],q[4];
u1(2.24289904305224) q[4];
u3(-0.249876989404031,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.54576240788498,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.44274045833681,-0.804029258942969,-1.15497883100391) q[4];
u3(0.900945041005622,0.598319441392055,-1.62893207590745) q[3];
u3(2.49021745999771,1.77306961222473,-2.90217494549348) q[2];
u3(2.68300903409418,0.400171552955698,-4.42231350266056) q[1];
cx q[1],q[2];
u1(0.117031683362176) q[2];
u3(-1.23616726037799,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.39784864223179,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.24171425996871,-0.0143967678972660,-1.39797567931699) q[2];
u3(1.95870130974329,-0.682952104999528,-5.32131622828738) q[1];
u3(2.49560317297118,1.96245840340541,0.945112538019327) q[2];
u3(0.626654206676672,-4.90522645281720,-0.105686950793886) q[4];
cx q[4],q[2];
u1(1.64121832215027) q[2];
u3(-3.46879337522817,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.23492316477413,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.17423528405745,2.10249983087436,-3.07298581478784) q[2];
u3(1.47558760272502,1.50147314216887,3.37510974530570) q[4];
u3(0.439566744632850,1.26940329167256,-0.353737306070854) q[1];
u3(0.916463595016964,-2.61969800839504,0.787738973168704) q[3];
cx q[3],q[1];
u1(1.19494857832643) q[1];
u3(-0.211567163730652,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.04921402527180,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.43660765257189,3.32125536248406,-2.57415743212990) q[1];
u3(1.04403380853465,1.15587646048244,0.0994318603343045) q[3];
u3(1.72898602479060,0.787975556550906,-3.60049788274016) q[3];
u3(0.870617525736766,3.18192876065983,-3.00079361926774) q[4];
cx q[4],q[3];
u1(-1.22320317641060) q[3];
u3(0.458041246295561,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.37650492204636,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.42008748461812,1.53073468445288,-2.50073961153479) q[3];
u3(0.588157517057630,-0.203650633540644,0.636054097998208) q[4];
u3(0.0179727928564108,2.97151047210111,-2.83961937776823) q[1];
u3(0.211359325414694,-4.23882223133474,1.58894095147391) q[0];
cx q[0],q[1];
u1(2.66211397435066) q[1];
u3(-2.11702147393755,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.19916545262912,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.41938129410828,0.614286685374681,0.737047494552815) q[1];
u3(2.06181711767154,0.997065841085403,1.66204649008860) q[0];
u3(2.27229669074550,2.10726902450208,-3.62282811204409) q[2];
u3(1.39922383248693,3.09840207271742,-2.26361201609116) q[4];
cx q[4],q[2];
u1(-0.173614096916191) q[2];
u3(-1.81612657044664,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.68601315213667,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.46634700748950,1.38939463389154,-3.80118303237900) q[2];
u3(1.45764668524255,1.40999492384580,3.63669644216664) q[4];
u3(2.77585220273904,1.33045159750425,1.23181458507753) q[0];
u3(0.637357659554423,0.554615170879680,-5.41868562984392) q[3];
cx q[3],q[0];
u1(2.82607022274403) q[0];
u3(-1.89732718474956,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.18361727832174,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.00743124049431,-2.79735215084714,2.08607493080594) q[0];
u3(2.63563680592441,-0.703550334213198,1.77081551331826) q[3];
u3(1.01514013513913,-0.621027689639031,1.16614323956297) q[2];
u3(0.475107629963562,-2.22310081479674,1.45642748719665) q[3];
cx q[3],q[2];
u1(1.82022156710512) q[2];
u3(-2.10281828358126,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.66559073933065,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.515351144053096,-0.424276139346535,-3.06764302692036) q[2];
u3(2.10354936742670,2.90158558383053,1.10863517697404) q[3];
u3(1.13411541460227,-1.29459693205158,0.541606062495560) q[0];
u3(0.915583601768218,-1.62365181469237,0.599706422760737) q[4];
cx q[4],q[0];
u1(1.47761615804570) q[0];
u3(-2.25921712569256,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.29058178039029,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.661582714756526,-0.864260084486185,3.45753905483611) q[0];
u3(2.29560347626537,-0.334461697532582,-5.64725817635895) q[4];
u3(2.59526202106905,0.376449649835097,2.64647224885365) q[3];
u3(2.62286907474407,-0.793149831697014,1.04785216449668) q[4];
cx q[4],q[3];
u1(3.34641011834968) q[3];
u3(-3.59301795526829,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.877708402411083,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.78616751095707,-2.45068445031384,2.42100384198713) q[3];
u3(1.60303484936401,-1.92148216633786,-1.61388150343921) q[4];
u3(1.44946349399731,1.17095217439224,1.22160949308610) q[2];
u3(1.36304300005049,-1.72791431128019,-1.84372475267946) q[0];
cx q[0],q[2];
u1(0.0189675118516701) q[2];
u3(-1.96538479613773,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.808829755710150,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.65721235463839,-1.12881882915120,-0.478473096350964) q[2];
u3(2.87566318918582,-0.338335535698377,-5.40208373273318) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
