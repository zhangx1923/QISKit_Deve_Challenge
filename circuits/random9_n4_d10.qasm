OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(2.41464852828875,2.31600111854618,0.0733991774957787) q[3];
u3(1.90190664221644,0.794754564390297,-3.21201357588748) q[2];
cx q[2],q[3];
u1(0.300343213712928) q[3];
u3(-1.23884765019853,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.58559616985891,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.275330777519273,1.47181863505901,0.0208521173251162) q[3];
u3(1.30167091566776,4.29203478546307,0.502538293492114) q[2];
u3(1.10280981724643,0.140215353842559,0.727779095888235) q[1];
u3(2.16320755920142,-0.994256100881188,-1.90597148312334) q[0];
cx q[0],q[1];
u1(-0.354892133042320) q[1];
u3(-1.85008668518146,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.05053185056596,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.25008505272994,-0.758588890941178,4.72509149163022) q[1];
u3(1.87082868810336,-4.30982986171485,0.0420927891845491) q[0];
u3(2.68125576955400,-4.43565261077991,1.58538467204960) q[1];
u3(0.724003742560829,2.13459000952906,-0.918994694281738) q[0];
cx q[0],q[1];
u1(1.44745625256463) q[1];
u3(-3.06484985827302,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.150696690589730,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.10948010797145,-2.94625108268400,0.509543752155987) q[1];
u3(1.89401585186360,-0.527419518125313,1.66604703637579) q[0];
u3(2.22927297513278,1.45800112343320,-4.32583312851380) q[2];
u3(0.305255698681613,-2.09957201731654,3.92971043999698) q[3];
cx q[3],q[2];
u1(-1.27718895224555) q[2];
u3(0.0674072612436973,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.38625941572425,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.946120062932722,-1.70464321229004,2.71832804868597) q[2];
u3(2.15929699977743,0.174125365102297,5.36636825519660) q[3];
u3(2.89584435220019,2.08785819437863,-1.01636992194609) q[1];
u3(2.09834499853680,0.977762333561940,-4.95559512249476) q[0];
cx q[0],q[1];
u1(2.51246940905941) q[1];
u3(-3.01611210320528,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.781682348979048,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.03571525893276,1.54410652538447,1.64675330739747) q[1];
u3(1.42178668267678,2.51181845973269,1.05574678022873) q[0];
u3(1.55170422906793,0.101268161389799,1.64508711677012) q[3];
u3(1.04016734112639,-1.16076182379161,-0.595039286451578) q[2];
cx q[2],q[3];
u1(1.28924458341322) q[3];
u3(-0.908751094886810,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.90525867219232,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.49419147223687,1.02529262719702,-4.20001633575704) q[3];
u3(2.08312527532716,1.19364280433427,-4.24334626074690) q[2];
u3(1.87084594169810,0.244253957672039,1.58387762572769) q[3];
u3(1.54819881036175,-1.03065873514125,-2.48656072056925) q[0];
cx q[0],q[3];
u1(3.09917608099567) q[3];
u3(-1.57009653924697,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.17700561689879,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.907339935593210,-0.197534231709585,-0.780590587890091) q[3];
u3(1.83183884973093,-2.97478556738495,-1.17293398483607) q[0];
u3(1.63393457927476,-0.712489457425845,-0.00502493649385810) q[2];
u3(0.270048229347405,-2.39653493216013,-1.15741234073454) q[1];
cx q[1],q[2];
u1(1.31844722097708) q[2];
u3(-0.281019936844884,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.45723261042005,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.596557872940653,2.20795372898121,-0.585015658008233) q[2];
u3(0.765311096037219,-3.55546453483703,2.31157544669943) q[1];
u3(2.20892888460819,-0.974023552097972,0.270281985698292) q[1];
u3(0.719449214792267,-2.10297803810743,-2.26107415716519) q[0];
cx q[0],q[1];
u1(1.56666778618155) q[1];
u3(-3.46760392270765,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.36602818258628,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.566163805266295,-2.02182885586977,3.63230431392141) q[1];
u3(1.84433745927350,-3.40975993371438,-2.39419469281047) q[0];
u3(0.679649165281011,-1.65098623681104,1.34392230043117) q[3];
u3(0.130709566796325,-3.15125757968523,2.46254136009033) q[2];
cx q[2],q[3];
u1(1.57626716653899) q[3];
u3(-1.07839121340509,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.434350906136234,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.56700605471531,-1.01263500916787,1.52924143277657) q[3];
u3(0.784293616487850,4.75974000017463,-0.780006020801078) q[2];
u3(2.21005333593808,-0.996054597375306,-1.01034270650488) q[2];
u3(1.80382035649619,-2.54334905273321,0.401911657242544) q[1];
cx q[1],q[2];
u1(2.05089163323418) q[2];
u3(-3.07627388773247,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.22468126437104,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.694065389828617,-2.40399741650809,3.27384584676985) q[2];
u3(1.59610483363511,1.74986846739085,2.66971306624750) q[1];
u3(1.70517638064722,-1.84277842829466,0.540587680027610) q[0];
u3(2.70894315241581,-3.14080485284090,1.00249276806500) q[3];
cx q[3],q[0];
u1(1.75794701740627) q[0];
u3(0.642425257306742,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.00462402708059,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.73089035913394,-2.16198919474061,-1.19255979352208) q[0];
u3(0.538574052825675,-0.377569647117857,-0.898468475818395) q[3];
u3(2.28580938971971,1.28421330664163,-1.73899712714087) q[1];
u3(2.56850011116303,-0.159597784112093,-3.60922304014482) q[0];
cx q[0],q[1];
u1(-0.161109115781318) q[1];
u3(-1.68027484824241,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.00857858443754,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.74572768077743,4.09915804649177,0.482948605063286) q[1];
u3(2.27759804847038,-1.75001764819744,-3.35513665302500) q[0];
u3(1.06843301461282,2.39887529672345,-1.33636786652756) q[3];
u3(1.47842941903259,1.55814021862873,-0.478039860143892) q[2];
cx q[2],q[3];
u1(0.141577996719414) q[3];
u3(-1.14639003689632,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.80460454991704,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.66207843136210,1.68310352961991,-2.40351603428628) q[3];
u3(0.791679043992845,0.292236113264256,-2.03433314229715) q[2];
u3(0.967483718789095,1.91544582747609,-1.36009402933767) q[0];
u3(0.702567179113346,-2.37688838544883,0.995755770415137) q[3];
cx q[3],q[0];
u1(2.61945220742040) q[0];
u3(-3.12257160138751,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.06148077081360,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.887173554942060,3.57352004985679,-1.16178240439360) q[0];
u3(0.402585999766802,1.59046116352694,0.583569292407484) q[3];
u3(1.16363572356785,0.967645726290571,0.370333032586009) q[1];
u3(0.856174278654424,-1.13234567341388,-1.94388703914141) q[2];
cx q[2],q[1];
u1(2.94828210592392) q[1];
u3(-1.33309486092603,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.48784955031410,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.11099706922342,-2.47292037705906,3.39846444899690) q[1];
u3(0.513017562654239,3.75570036377451,0.403550294944963) q[2];
u3(2.57161675913185,-0.313762140817655,-2.26634791150210) q[1];
u3(1.82518354686404,-4.18341253065593,2.02184298138563) q[2];
cx q[2],q[1];
u1(-0.773070751869440) q[1];
u3(1.15604653513979,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.54479924031450,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.64492053568508,-2.37800070688224,0.565515504693607) q[1];
u3(1.98878398612026,4.18831650437795,1.97585866411810) q[2];
u3(2.26442105945789,0.234293604889235,1.89199470373740) q[3];
u3(1.99411915242939,-1.06745470359683,-0.377780016954558) q[0];
cx q[0],q[3];
u1(1.35398099042776) q[3];
u3(-3.48295521209956,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.32366025768062,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.769539964012536,0.287319332709334,-1.19623629423387) q[3];
u3(2.22249757783287,0.962942424003691,-2.44994517096240) q[0];
u3(2.64693546175439,-3.58775997375497,2.48591351947056) q[2];
u3(1.19311097672529,3.54251911044272,-2.07281798751844) q[1];
cx q[1],q[2];
u1(-1.30107142205758) q[2];
u3(0.203395521127263,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.51816117553359,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.85677415093446,-1.74194405562852,-1.04330670612948) q[2];
u3(2.46475030580781,0.317007961715753,-1.42019564337077) q[1];
u3(1.75990296174384,0.410218653241045,-3.54680454905798) q[0];
u3(1.64596916543311,-1.32400385562443,4.83961042280328) q[3];
cx q[3],q[0];
u1(1.60379762895920) q[0];
u3(-2.27432470553424,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.12717373415033,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.93921149250744,-2.87930143175058,0.770572923774979) q[0];
u3(1.30063424927946,-4.87525230302032,-0.481289429417322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
