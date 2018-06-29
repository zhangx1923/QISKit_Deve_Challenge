OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.75010050410861,1.62243955447982,-2.46008607655511) q[6];
u3(2.50539208581615,-4.15614560955416,1.98067863733672) q[2];
cx q[2],q[6];
u1(2.27334330176838) q[6];
u3(-1.50730922254985,0.0,0.0) q[2];
cx q[6],q[2];
u3(3.79408465194393,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.46786929353583,-1.74961483822995,-1.26751793404843) q[6];
u3(0.802058652222992,3.01564840512517,-1.24750070985691) q[2];
u3(0.873664167747903,2.01774867970019,-3.57640553761117) q[1];
u3(0.969967424089145,2.80967029846598,-2.55133290928408) q[3];
cx q[3],q[1];
u1(1.20936043645768) q[1];
u3(0.172713345991612,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.26713114676307,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.08976750326441,-0.424780587726308,-3.21895674526332) q[1];
u3(1.14888199845707,-3.15790671095441,-1.47315651313810) q[3];
u3(1.60529812972229,0.215694804059523,0.481858122777024) q[0];
u3(0.614322269694535,-1.74638299624880,-1.88988825493573) q[4];
cx q[4],q[0];
u1(-1.16474059158209) q[0];
u3(0.240911955222038,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.72531654446259,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.903127153942436,5.11743771295837,-0.497998593519887) q[0];
u3(0.845091806984222,2.37068084927667,2.80602045043056) q[4];
u3(2.38472026678427,-0.837916565914234,-0.315258833358667) q[1];
u3(0.450918562568571,0.557605659003996,-4.57715043465569) q[2];
cx q[2],q[1];
u1(3.26605368845857) q[1];
u3(-0.781885638043757,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.99877975620140,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.918658612530876,1.43018816618276,-0.366162641355566) q[1];
u3(1.92268557063637,2.30352626304690,-0.897255450079889) q[2];
u3(0.696074879459524,0.258703648656685,-2.77000020758928) q[6];
u3(1.26503208095419,0.779447420290998,-4.62377704425622) q[0];
cx q[0],q[6];
u1(3.20412545881475) q[6];
u3(-1.56742320767957,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.41524287723129,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.22659995192666,0.708528475436572,-1.08440197893230) q[6];
u3(0.157005491998880,-1.46328301770569,4.69442698539266) q[0];
u3(1.82110901067424,-0.203975331840360,2.26133488586227) q[5];
u3(2.88581033982083,-1.59651167650867,-0.555752097378665) q[4];
cx q[4],q[5];
u1(0.572257361613067) q[5];
u3(-0.952707806972426,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.89998221109042,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.79587444394655,2.71719851332243,-3.43652583133776) q[5];
u3(1.45768391133280,0.639216115949998,-3.60327607348574) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
