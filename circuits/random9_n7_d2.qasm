OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.57653940620706,-0.894175646016651,0.431422680685710) q[3];
u3(1.89483002812396,-2.80906382982975,0.103094274182761) q[1];
cx q[1],q[3];
u1(2.81358347767748) q[3];
u3(-1.73854867793992,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.22754995319204,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.94630106091668,0.483464331898048,-1.62412145639668) q[3];
u3(1.88453190546072,-0.359026005747095,2.34726243125218) q[1];
u3(1.14460899772838,3.51324221517188,-1.84654232040514) q[4];
u3(1.81131798988323,1.18582020476508,-2.49240143616801) q[6];
cx q[6],q[4];
u1(3.40577795432758) q[4];
u3(-1.04627068687319,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.34162949253573,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.08134967938618,0.951256822527629,0.803414595242472) q[4];
u3(1.85265275343342,3.09505473682565,-0.545422394100585) q[6];
u3(0.272071986943991,-0.697841108528694,0.152698534015811) q[5];
u3(1.20361022323014,-0.493799863434366,-1.02254476800268) q[0];
cx q[0],q[5];
u1(1.47750392570894) q[5];
u3(-0.227213216386334,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.37691430241371,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.982444425968689,-1.39787902327874,2.94141549193736) q[5];
u3(1.24709402329852,1.81938245909476,-1.23372741603773) q[0];
u3(2.43575505678304,-0.589807183176980,1.76902164974055) q[5];
u3(2.12209785011148,-2.62064025610075,-2.13141985402755) q[2];
cx q[2],q[5];
u1(1.37586861145479) q[5];
u3(-2.90858231685607,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.78418653007206,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.60341131067796,-1.24893677067607,0.0275685820009905) q[5];
u3(1.01129880053481,-5.30880899356923,0.600853756444601) q[2];
u3(1.38271840188293,-0.285619234738471,1.04548131578030) q[1];
u3(0.909717196898389,-1.59739369963776,-2.06453982568596) q[3];
cx q[3],q[1];
u1(0.424291920431275) q[1];
u3(-1.27272263767953,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.0241882870652022,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.07779669726285,1.66753609351245,-0.103373497013645) q[1];
u3(0.208352220984354,2.29211531773740,-1.71247535780975) q[3];
u3(2.35831092934099,-1.75839163174615,0.538635785336886) q[6];
u3(2.04328030781904,-3.96126747143457,0.747621107992989) q[0];
cx q[0],q[6];
u1(1.54552604062096) q[6];
u3(-0.780975328453786,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.58795788298263,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.19483164674746,-2.07107960649870,2.91762419333074) q[6];
u3(1.27134685702423,3.60288896275721,-1.41005990506413) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
