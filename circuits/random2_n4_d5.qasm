OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.75063563127124,-0.425195168384591,2.51656589813048) q[0];
u3(1.17371833922041,-2.27682983763980,-1.72208698096855) q[3];
cx q[3],q[0];
u1(-0.304059478025590) q[0];
u3(1.00161408188381,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.56827678052684,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.62286802259248,-3.13790583970978,1.77844470201577) q[0];
u3(1.20559444167376,-2.89916214373178,-0.266718200464322) q[3];
u3(1.21915172573341,2.31709334859991,-1.22377055911361) q[1];
u3(1.11130177134650,0.582366924561584,-3.30500570226978) q[2];
cx q[2],q[1];
u1(0.332481877741158) q[1];
u3(-0.552834134896667,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.35517655027614,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.49102737055530,3.36178709616406,-2.89131291073874) q[1];
u3(1.64540472365623,-2.06484057214620,0.943308543644546) q[2];
u3(0.253082377195675,-2.36409914755656,3.00332173197056) q[2];
u3(0.945851852795539,0.459533840055018,-1.86405096537699) q[3];
cx q[3],q[2];
u1(3.22767789673711) q[2];
u3(-1.58297498505925,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.351515969529691,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.190560951045340,2.68113809201335,-2.18119824893708) q[2];
u3(2.53164138976011,2.25009492843477,-4.00627351271277) q[3];
u3(1.39486066999207,-0.244139731830712,-1.51257995680483) q[1];
u3(1.72449860231404,0.398026358033632,-4.74225253512844) q[0];
cx q[0],q[1];
u1(1.65656129436456) q[1];
u3(-2.44112871610775,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.511525098763444,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.15320797168950,-1.82281738421432,0.246208501942187) q[1];
u3(0.606603185870026,0.219910447439867,4.37049456001165) q[0];
u3(2.25733724979245,1.91976694400111,-3.26519260322897) q[0];
u3(0.984289089648438,-1.96883367289434,2.80069353559603) q[2];
cx q[2],q[0];
u1(2.06096242595304) q[0];
u3(-1.60714881909114,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.589988267501134,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.44532742420811,-0.791841815102753,1.37274935522929) q[0];
u3(2.48253039039609,3.12147378623520,-0.349311895820188) q[2];
u3(2.36733758695093,-2.26216431444367,3.40742119065948) q[3];
u3(1.45890778011881,-0.930100201740366,2.21264404263755) q[1];
cx q[1],q[3];
u1(0.0224220292526580) q[3];
u3(-1.28433855633763,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.74670281152797,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.09140513361147,0.602683441532147,-0.0192620260663904) q[3];
u3(0.997148093463938,-5.80132426094851,0.0188830952664447) q[1];
u3(1.71702588983577,3.10151273530713,-1.21810065298638) q[3];
u3(2.23314226198833,1.62500357982162,-2.61134265934691) q[1];
cx q[1],q[3];
u1(2.08433266554054) q[3];
u3(0.595164375394141,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.18502819206799,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.56052076465481,-2.68349749150857,1.90253687782866) q[3];
u3(1.11103861258290,2.60840715583645,-1.03484754430094) q[1];
u3(1.53174907146744,0.914015002600260,1.66063312706383) q[2];
u3(2.03863572412941,-1.75189040016556,-0.599989873916610) q[0];
cx q[0],q[2];
u1(-0.560023395187508) q[2];
u3(0.326827887655442,0.0,0.0) q[0];
cx q[2],q[0];
u3(4.34304457879068,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.24230486568863,3.09962062057529,-1.00215476120168) q[2];
u3(0.866980009404866,5.14180542233778,0.219524611112540) q[0];
u3(1.08510373003941,-1.53771974885400,0.962417317141858) q[1];
u3(0.608939410772803,-3.66263705701345,1.47152629033577) q[3];
cx q[3],q[1];
u1(4.34840859194101) q[1];
u3(-3.53600736852656,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.0361726114867615,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.485490196160807,-2.56473173764503,1.59856787713501) q[1];
u3(2.05504389646326,4.41636009335328,1.82848397564339) q[3];
u3(1.15665816832386,-1.44377875804436,-0.833518914309968) q[2];
u3(1.33506662820249,-3.65428212478375,0.842795818337620) q[0];
cx q[0],q[2];
u1(0.414195397203402) q[2];
u3(-1.99109104582807,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.97896780390971,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.0454521029709102,-1.97138470038684,0.311603007083328) q[2];
u3(2.88457551025029,2.16072220135825,-2.49355145030247) q[0];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
