OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0.909337059992479,-0.703187370477347,0.0264702854621424) q[0];
u3(1.37658369309170,-3.10764653337249,-1.27818719599519) q[2];
cx q[2],q[0];
u1(1.07691225137450) q[0];
u3(-0.175161576501209,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.36596458412828,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.955366028267195,0.811392348705691,-2.28574311579316) q[0];
u3(2.39262710722981,-1.15709561347456,0.235563048482978) q[2];
u3(1.39768964186083,0.130960394009679,0.762750856863966) q[5];
u3(0.975013760363926,-2.62000138249612,-1.64406140643863) q[1];
cx q[1],q[5];
u1(0.673497097409622) q[5];
u3(-0.198437832514725,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.22472359387825,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.68021182602539,3.54908875359418,-0.715906966806440) q[5];
u3(1.56964217357185,0.964482184074148,-5.15417915061502) q[1];
u3(2.00539098789633,-3.90546281691561,2.29714123211196) q[3];
u3(0.746555381055859,2.65974094945760,-0.679505881717830) q[6];
cx q[6],q[3];
u1(1.78701638297739) q[3];
u3(-0.293345085245924,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.04546086905437,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.72899502908192,-2.66972248845277,2.85778250344359) q[3];
u3(1.56687590625394,2.78538978432622,0.0262911836159425) q[6];
u3(1.50034579732963,3.89688087222806,-0.877692450045382) q[2];
u3(0.563345502611474,1.71771762570165,-1.33565712322369) q[6];
cx q[6],q[2];
u1(0.637444471890433) q[2];
u3(-1.12377760098152,0.0,0.0) q[6];
cx q[2],q[6];
u3(3.12310832568603,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.0857165257372846,0.485662112480141,-0.128070024993018) q[2];
u3(0.674704118430925,-0.391790136776959,4.60300858292882) q[6];
u3(1.63891006953535,-0.473497924367064,1.27109689795768) q[1];
u3(1.52272131001717,-1.77582212578251,-2.63968098184689) q[0];
cx q[0],q[1];
u1(1.01947796366943) q[1];
u3(-0.179628244998193,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.55711017566349,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.62641397872084,0.311815700334594,-0.241649444656492) q[1];
u3(1.87668982707733,3.03254699471202,-2.33383396713366) q[0];
u3(0.988368923663540,1.40101634464045,-1.62043401250900) q[5];
u3(0.321614546884052,-4.49917155113340,1.42182804467439) q[3];
cx q[3],q[5];
u1(1.83256890780705) q[5];
u3(-2.99754343092820,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.757546051146773,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.935456235840116,-1.44936312152755,0.643111824604795) q[5];
u3(1.48964902585655,-4.56142117554715,0.251335630829843) q[3];
u3(2.45897181406276,0.767344979375514,-0.624988266425320) q[0];
u3(1.71824670167629,-5.22241018556707,0.974733070944554) q[5];
cx q[5],q[0];
u1(2.35303554050460) q[0];
u3(-1.90169997599017,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.713082925924857,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.953640742193214,0.0981264086780651,-3.00478991613011) q[0];
u3(1.17159597034886,3.38971690420931,1.95872632004126) q[5];
u3(1.57709680312443,0.394944123232187,-1.60917718655776) q[1];
u3(0.726510507857081,-3.97455722449122,1.62876332355359) q[4];
cx q[4],q[1];
u1(2.06340565796425) q[1];
u3(-0.00360169009620548,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.805589505005168,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.84823191561013,-1.93610584781353,3.46455517508181) q[1];
u3(1.00773602096346,4.50085981650821,0.136412522448240) q[4];
u3(2.34026363588481,-2.30163875250110,-0.275637481139103) q[6];
u3(1.72176678171286,-3.64147690365244,-0.313074381930909) q[2];
cx q[2],q[6];
u1(0.447904569710528) q[6];
u3(-1.43563999204387,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.28934151638749,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.98174001157096,-0.701867983951706,-1.50501160157473) q[6];
u3(1.53976350099919,-3.40523540396703,-1.33655694982706) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
