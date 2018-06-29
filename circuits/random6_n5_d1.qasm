OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.42650774490738,-3.42017812079316,1.01339568037501) q[1];
u3(2.66610659872537,1.20785031037448,2.90864696297077) q[0];
cx q[0],q[1];
u1(1.01466627520271) q[1];
u3(-0.508060879465990,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.98948264049098,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.403049481472487,2.39740183726253,-2.20429274867641) q[1];
u3(1.36993784878966,-0.159962830476568,-4.23350648638130) q[0];
u3(2.54272299121011,0.633824074710255,-3.26086398198558) q[2];
u3(2.12667122106930,3.18442152978144,-2.67450579557959) q[3];
cx q[3],q[2];
u1(0.243958239299737) q[2];
u3(-1.57831066365536,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.42931879952901,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.733612361608273,1.73972975732213,-1.18421696956065) q[2];
u3(2.79516581537544,0.968866746276189,4.29257552071968) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
