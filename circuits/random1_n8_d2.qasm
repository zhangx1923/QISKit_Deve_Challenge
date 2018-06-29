OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.61487738919498,1.18675752064087,-2.67364422166325) q[0];
u3(1.63808250592739,-2.08370463023827,2.87835722762371) q[6];
cx q[6],q[0];
u1(2.13018851720501) q[0];
u3(-2.83502422145020,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.42175314430674,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.77780697199029,1.82861863859308,-4.18266701695739) q[0];
u3(0.951188218139094,-1.14662335407858,0.602847716966160) q[6];
u3(0.329809170905941,-1.60617893400787,-0.931993776062148) q[1];
u3(2.29913828104088,-3.26036287038164,2.35986502523851) q[2];
cx q[2],q[1];
u1(3.41319802491782) q[1];
u3(-1.18162317990443,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.90316692770061,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.40313477261449,-1.82086524679187,3.17159204860778) q[1];
u3(1.18656658518699,0.417458730422722,-3.85381496277409) q[2];
u3(1.36499612901359,3.36862671848611,-2.84009140432805) q[4];
u3(2.61231572011779,0.789781661144339,-2.47099332382250) q[3];
cx q[3],q[4];
u1(0.0137870819594150) q[4];
u3(-0.801457991596513,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.54263518343731,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.24528326333816,-1.16007659595847,-0.106329713794562) q[4];
u3(0.275143546500785,-2.72900333698858,3.33254113186809) q[3];
u3(1.25825416063299,0.648946004200342,-0.835691233215711) q[7];
u3(1.21601710560038,0.179803278810037,-2.45677597346002) q[5];
cx q[5],q[7];
u1(1.56077978178552) q[7];
u3(0.815930503936632,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.06241946067790,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.902885023322838,0.391510220988969,-1.50723054221806) q[7];
u3(1.54485471387144,4.32216348446612,-1.68517001896540) q[5];
u3(0.920220938333771,1.84431866413452,-0.146974301554638) q[0];
u3(1.87883330736454,0.0921831401824356,-3.75555557874896) q[1];
cx q[1],q[0];
u1(1.06044533052526) q[0];
u3(-0.432954284860644,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.16637100729800,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.44576940596134,1.51454054893171,2.63107401028539) q[0];
u3(0.551535598525733,-3.18221675199447,1.24175350344014) q[1];
u3(2.19923259876960,0.309014828872944,0.712104037972544) q[2];
u3(2.15237583432437,-2.10627973571156,-0.941437420313647) q[4];
cx q[4],q[2];
u1(0.457500725600106) q[2];
u3(-1.19518544926261,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.14141602707725,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.44866420429833,-2.93934203997551,0.793778921297555) q[2];
u3(0.625947247846927,1.89123284508830,1.24578883894870) q[4];
u3(2.97479163559805,-0.195609592834983,1.69925480215801) q[6];
u3(2.26237663371386,-1.41105422421996,-0.151320482546012) q[5];
cx q[5],q[6];
u1(0.822453253761153) q[6];
u3(-0.0609803542460055,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.81033476330486,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.709644068348313,-1.87898583833737,1.37939935431469) q[6];
u3(2.43737057164533,-0.338654434256482,-2.43043838995087) q[5];
u3(2.43322674016737,-1.21670382139421,-1.85204580305045) q[7];
u3(0.703415582774196,-0.719967399298953,-3.68090119422603) q[3];
cx q[3],q[7];
u1(0.131425331301668) q[7];
u3(-0.379289997247486,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.48105423525408,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.03657391217254,0.845773732428119,-1.69068559485616) q[7];
u3(2.21485356501240,0.769783857978445,1.99188088127979) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
