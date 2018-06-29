OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.73866943787813,-1.78532120618247,4.15695488906556) q[4];
u3(0.515434771713103,-0.955793608179301,1.65991896903458) q[2];
cx q[2],q[4];
u1(1.54878074526536) q[4];
u3(-0.654394758785813,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.32565135335898,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.16347279901337,-0.679268997836289,0.650235480901147) q[4];
u3(1.68128084321768,3.49639709494150,0.306356850705630) q[2];
u3(2.03550247535047,2.21854598009832,-2.34228240293497) q[3];
u3(1.51195650708752,-3.44959596066211,2.63896100867162) q[0];
cx q[0],q[3];
u1(0.551912134897197) q[3];
u3(-0.941565030796930,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.48058383649785,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.11846479759178,3.97555980871524,-1.62528004196117) q[3];
u3(1.52527985229019,1.87881065326811,3.42557872098984) q[0];
u3(1.79780067359907,1.03846057515449,-1.54969989633431) q[2];
u3(0.999353947370932,1.45385080877022,-4.09108066723766) q[4];
cx q[4],q[2];
u1(1.71998532136861) q[2];
u3(0.215450033526110,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.08499044860822,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.05497685546829,1.36853379027272,0.719206459169629) q[2];
u3(2.24956661408314,2.06258165069955,2.88392424367046) q[4];
u3(1.25818911578484,2.10322836069170,-3.15848498306653) q[3];
u3(1.79835952832111,-3.07221874782652,2.83816970978387) q[1];
cx q[1],q[3];
u1(2.92773680382375) q[3];
u3(-1.67494989495118,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.625995946550180,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.624463319160608,0.235550366969648,2.26961354239792) q[3];
u3(1.75382265056451,0.893790411990401,-2.07480284209307) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];