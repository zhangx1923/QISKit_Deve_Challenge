OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.67258585505495,4.03014252291242,-1.31967037954283) q[1];
u3(2.20494151897613,1.12388688759741,-2.44481780508401) q[8];
cx q[8],q[1];
u1(4.36562401921901) q[1];
u3(-3.40278504596319,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.582510782367347,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.928655798808649,-3.02444668658083,0.904039796575545) q[1];
u3(2.26877022193147,2.95965045622644,-0.934346466361612) q[8];
u3(2.12418021934585,1.21167662566505,-3.64716022872639) q[3];
u3(2.56119235669752,1.56461622719854,-3.17000895491597) q[5];
cx q[5],q[3];
u1(1.48432360756369) q[3];
u3(-0.573161370672168,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.94472214676063,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.17465121846676,0.236072527418321,1.63749129588698) q[3];
u3(2.54004218658194,-0.333324266872466,2.11318541899660) q[5];
u3(1.53743324256938,-0.798827462750023,0.582208726112627) q[0];
u3(2.73570214450105,-0.953243509691547,-1.54786696914436) q[9];
cx q[9],q[0];
u1(-0.365745664728217) q[0];
u3(-1.86216711128289,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.903522264920569,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.60675836560483,0.819772934549019,-0.162530448034955) q[0];
u3(2.70117902609138,-2.66946207683919,-0.962867978023455) q[9];
u3(0.714818786901876,-1.06375255645380,2.52750201875204) q[6];
u3(1.12341059173941,-1.61261573411253,-2.74880325243476) q[4];
cx q[4],q[6];
u1(0.380480825904771) q[6];
u3(-1.54551428511268,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.99298906571314,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.687017740870453,0.882410195506377,0.0961385996796655) q[6];
u3(1.39938344721711,1.82351926185626,-4.29594359702474) q[4];
u3(1.14134940344211,2.46558061119508,-0.603579314336876) q[7];
u3(0.547668527952851,0.751771071168364,-2.30846227674538) q[2];
cx q[2],q[7];
u1(3.55162056494922) q[7];
u3(-3.94726297526225,0.0,0.0) q[2];
cx q[7],q[2];
u3(-1.16562104965029,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.10233035117939,-4.29066409067093,1.61634590039545) q[7];
u3(1.24406268306381,1.42989438069473,-1.47123073700482) q[2];
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
