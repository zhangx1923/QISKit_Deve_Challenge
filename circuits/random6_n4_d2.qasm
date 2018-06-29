OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.52310750009014,3.69108986996388,-1.61170189737317) q[3];
u3(0.735716000932522,2.16486031358452,-1.10868754564603) q[0];
cx q[0],q[3];
u1(-1.27945773882275) q[3];
u3(0.660013910213916,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.97172376960947,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.816293287063853,0.134991078213370,-4.60990001085621) q[3];
u3(1.80278869315959,1.02922075660062,-3.49764411125070) q[0];
u3(1.85104913798114,-2.54135547988372,-0.124443426413534) q[1];
u3(2.18591461083106,-4.49759801363976,-1.51958809887021) q[2];
cx q[2],q[1];
u1(2.12527863746862) q[1];
u3(-1.97802526782248,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.787562143544699,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.12881440978225,1.91893876551208,-2.64659435060975) q[1];
u3(2.30122978106336,-3.92643315364168,0.820451390653232) q[2];
u3(1.35229637986353,0.692552993918640,-0.123536867242180) q[0];
u3(2.36476749019020,-1.02611681199969,-4.56161475879249) q[1];
cx q[1],q[0];
u1(1.55876606528264) q[0];
u3(-2.23756846043569,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.428935599055507,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.02596277762082,0.639207667006789,-1.66838090639461) q[0];
u3(0.520156906510112,3.25751530877719,-2.31200638418504) q[1];
u3(1.68060968815227,-1.56455171076994,-0.175393182414949) q[3];
u3(0.764028682039435,-4.20677920077372,-0.224466631711191) q[2];
cx q[2],q[3];
u1(1.91237466596710) q[3];
u3(-2.60185486252350,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.681602262766076,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.26594466566602,-1.95728385680781,0.0712067317555340) q[3];
u3(1.14183594829046,-3.13579007011690,-2.60149267836981) q[2];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
