OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.74016491864921,0.446712953080417,0.565996690587142) q[7];
u3(2.27150115429529,-0.745106708449456,-1.46411172452331) q[0];
cx q[0],q[7];
u1(1.40877419373730) q[7];
u3(-0.0401641228482785,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.81329109839458,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.12840982992941,4.24919665256143,-1.47990587204466) q[7];
u3(1.33377435920897,-0.371417506436046,4.24833763052915) q[0];
u3(2.73261567881830,-1.37022994608218,3.07560060818737) q[3];
u3(1.34040860966394,0.534264034151359,0.812793188982806) q[4];
cx q[4],q[3];
u1(-0.750143144854248) q[3];
u3(0.142675986787792,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.67268716024302,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.23539804175974,-2.32004140903370,2.31256079345777) q[3];
u3(1.22438208067810,4.27896615699194,0.580067975016682) q[4];
u3(2.67129169189159,-2.72567199698471,-0.240676272687299) q[1];
u3(2.52672567933665,1.21987214218140,2.44776355827627) q[5];
cx q[5],q[1];
u1(0.827858175459872) q[1];
u3(-0.529094906637629,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.93317406920688,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.41545233337302,0.653637656511079,2.44573488212418) q[1];
u3(2.59999395563848,-2.63330549848177,1.32910951234721) q[5];
u3(1.86437809198356,0.978893681345957,1.34098002595598) q[2];
u3(0.911567203396530,-0.487766894240314,-3.12734492014527) q[6];
cx q[6],q[2];
u1(3.33794868471446) q[2];
u3(-0.783396793341208,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.76286259259871,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.722863761449949,1.95188522250730,-1.42096778092252) q[2];
u3(1.11293496420962,-1.33776705787539,4.74907021110675) q[6];
u3(1.74697027255530,1.78876850573479,-3.57713305236318) q[5];
u3(1.18186500029885,-2.81383589296464,3.20391046874905) q[2];
cx q[2],q[5];
u1(1.58331583062042) q[5];
u3(-0.295886173023128,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.22601710010555,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.78084529639471,1.39428389735224,-3.26231003635695) q[5];
u3(2.56239809963382,0.217443158151504,2.66782288634352) q[2];
u3(0.711005357592656,-0.527713454269281,0.989004876825145) q[7];
u3(0.828781426345265,-2.08390047952895,-1.16904996843082) q[3];
cx q[3],q[7];
u1(2.21462081112656) q[7];
u3(-1.53767909146273,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.56411664862883,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.14730324939871,-2.36482837638234,-0.971108374524150) q[7];
u3(1.21466679208432,-2.88122798094803,0.871848929801265) q[3];
u3(2.28835424891802,-0.0374231500188358,1.24423260670360) q[6];
u3(2.05564230890694,-1.97952681439983,-2.56744416891680) q[1];
cx q[1],q[6];
u1(0.900152634987114) q[6];
u3(-1.26726072924632,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.105323389740058,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.475163337740429,0.134698686122702,-3.39602415533423) q[6];
u3(1.47581496172444,-0.0806051499826852,0.217361309237355) q[1];
u3(1.22370299633308,-0.989605994989941,-1.54530906791439) q[4];
u3(1.51084990846105,1.27047754797965,-4.56683555502099) q[0];
cx q[0],q[4];
u1(0.109506911229273) q[4];
u3(-0.962218271027944,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.49670935806089,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.03977355807778,2.95811064613992,-1.06681991203456) q[4];
u3(0.446863160467096,-1.01424019956503,0.287019121649927) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
