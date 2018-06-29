OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.27471328775107,2.85751408005685,0.0182215866142208) q[2];
u3(1.66271096540061,0.110166624591618,-4.41723993102754) q[5];
cx q[5],q[2];
u1(-0.0320059605848524) q[2];
u3(-1.30467083883685,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.28376194246505,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.59755633192002,-0.301701374686017,3.07484375244203) q[2];
u3(1.42973824381156,-1.74589402924781,-2.39827166860063) q[5];
u3(2.10294740033323,2.44074139164208,-3.08908467388884) q[4];
u3(0.839795009886390,-2.97493210472626,2.79952906969343) q[0];
cx q[0],q[4];
u1(4.32226272862061) q[4];
u3(-3.68936231527986,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.867750335961443,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.67926996530468,-2.40793591284770,-0.194259810040854) q[4];
u3(1.06192341883366,-3.23224534390581,2.82356041877989) q[0];
u3(1.77768565155136,3.50085238208121,-1.22583483836267) q[6];
u3(1.43765628355133,2.42192855251554,-2.97005536484486) q[7];
cx q[7],q[6];
u1(0.839847043993023) q[6];
u3(-3.48634278076340,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.61341614916734,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.94836299320298,-3.60992868358148,2.56440624621868) q[6];
u3(1.87670063135863,-4.67642807848329,0.616688377362674) q[7];
u3(2.78233499344165,3.73645799984662,-2.21359421356780) q[3];
u3(0.737515183621703,0.395557296432897,2.24563845299784) q[1];
cx q[1],q[3];
u1(2.45021961502612) q[3];
u3(-2.81194588529211,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.409517180336804,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.185030754732848,0.262383386189969,1.85572932277416) q[3];
u3(1.19801187128604,-1.34150140998848,1.29314628882871) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
