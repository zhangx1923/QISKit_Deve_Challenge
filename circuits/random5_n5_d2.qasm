OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.30229048439632,-0.645985327148003,0.834257996487505) q[3];
u3(1.10298724835647,-2.10960592239083,-2.03121779599543) q[1];
cx q[1],q[3];
u1(-1.30192230111776) q[3];
u3(0.460968816132175,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.55214242281011,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.49985401946584,-0.146227997881002,1.51514880731432) q[3];
u3(1.90532135145975,-3.52800094109418,-0.345533850009165) q[1];
u3(0.871752399332938,1.79901425591201,-0.943003217489311) q[4];
u3(0.497485007012352,-1.96370131435848,0.452368257309771) q[0];
cx q[0],q[4];
u1(4.14111300915553) q[4];
u3(-3.31688963144748,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.295798412776409,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.56844706888238,-3.06055032601721,3.00172098733096) q[4];
u3(0.850139680790006,2.30247168121973,-0.522814388612091) q[0];
u3(0.970923209537095,1.07842251374357,-1.54209592757328) q[1];
u3(0.441854396511823,-0.558481227820540,-1.20529281429337) q[2];
cx q[2],q[1];
u1(3.40563683713446) q[1];
u3(-0.394906799871462,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.60106395597448,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.06685312843435,1.89606558638078,1.33364831433836) q[1];
u3(1.19116602537516,4.83409942293696,-0.913957661207834) q[2];
u3(1.15548916966891,1.68161613493347,-0.200965001032224) q[4];
u3(2.27455907458483,-0.104587971570015,-2.46605693651717) q[3];
cx q[3],q[4];
u1(-0.393234248800656) q[4];
u3(-1.75283842452607,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.833411459711384,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.38383788602485,-3.77525837139329,1.41427055944676) q[4];
u3(2.29812675253641,-5.48781300335335,-0.360846794687261) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
