OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.70413165631329,-1.87976781514861,-0.190865512255909) q[0];
u3(1.55183755781199,-3.84693058195237,-1.41509368835193) q[1];
cx q[1],q[0];
u1(-0.337048208844538) q[0];
u3(-1.81080970407406,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.922479786751642,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.56915708890799,-0.423455305785610,-3.01136782684892) q[0];
u3(1.87791220897265,-1.28154980096636,4.87244119004202) q[1];
u3(2.58687697248269,1.59493793237058,-2.11968719921731) q[2];
u3(2.35994297759450,1.95296399511956,-3.83970013317952) q[4];
cx q[4],q[2];
u1(0.572397381680027) q[2];
u3(-1.65886047187729,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.95801567738332,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.98419469256288,2.08683041690632,-3.39440672524613) q[2];
u3(1.40295121386029,0.0772978554503057,-3.90384987329804) q[4];
u3(1.16418999841603,0.380524508791979,1.85424616887848) q[6];
u3(1.32177219069343,-1.23107299153319,-2.19773046586361) q[3];
cx q[3],q[6];
u1(-0.0988978658491968) q[6];
u3(-1.36365546139899,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.785976213666623,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.59149926256944,1.66586322588280,-3.58383914374526) q[6];
u3(2.11475201150150,4.18633141328613,-0.522252248549439) q[3];
u3(1.11544798479343,2.90098210614071,-1.97650457110358) q[5];
u3(1.75539862358855,1.34384974318804,-1.64511654839001) q[7];
cx q[7],q[5];
u1(-0.0723486785504448) q[5];
u3(-1.91313139087256,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.668310676510574,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.84059848048076,-2.57676186520700,2.24157564109711) q[5];
u3(0.435591852767611,2.99704930683397,2.62032687928264) q[7];
u3(1.73482856908142,0.337224354740542,1.18199560385825) q[7];
u3(2.13798531303564,-0.483259671474237,-1.42119881475118) q[2];
cx q[2],q[7];
u1(3.39645829312687) q[7];
u3(-1.75646568946295,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.51599977575798,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.883850320111216,2.23693388701173,-0.583224464952831) q[7];
u3(1.90461933026255,-2.16325103309686,1.05055998157309) q[2];
u3(1.49340346063128,-0.612354399833837,-2.05414768976179) q[4];
u3(2.22935596142282,-6.11973723529810,0.0256164977871394) q[1];
cx q[1],q[4];
u1(1.08090197918501) q[4];
u3(-0.308934269959957,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.41250869677672,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.69140556220229,2.12334055957983,-1.77049395396236) q[4];
u3(2.12492345091395,-4.32390740559755,-1.36541890040098) q[1];
u3(0.467349689333154,0.797163424627642,-0.164370014148757) q[5];
u3(1.19490666067578,-0.384692183044842,-1.10555700654937) q[3];
cx q[3],q[5];
u1(1.93011896776970) q[5];
u3(0.364656181652814,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.776054556863194,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.66160271504130,-0.821827433565893,1.76929908729203) q[5];
u3(2.49073487518584,-3.64937646671950,1.60561689396217) q[3];
u3(1.99100442655558,-2.56276276162457,-0.304690180627395) q[6];
u3(2.01708679372155,-3.70330805483961,-0.220200841525284) q[0];
cx q[0],q[6];
u1(2.80350037739568) q[6];
u3(-2.26727087613320,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.60718597562107,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.94354994897705,-1.51357296246527,4.70904979924888) q[6];
u3(1.44278280197422,3.80011765906246,-1.62184089653681) q[0];
u3(2.01410785881741,2.10511139800753,0.00636819029096047) q[6];
u3(1.25202571343723,0.782729901091398,-2.83899555242055) q[7];
cx q[7],q[6];
u1(1.59731387832313) q[6];
u3(-0.304601739814455,0.0,0.0) q[7];
cx q[6],q[7];
u3(-0.175805076417590,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.14086810052428,-1.47254273220340,-1.06282460781146) q[6];
u3(2.40216969398008,0.826938666139733,-2.59734208878169) q[7];
u3(1.30321352998512,-0.221718294393601,1.28058156237375) q[1];
u3(1.85400106723833,-1.68314690110269,-1.93052940029815) q[3];
cx q[3],q[1];
u1(1.39511726667266) q[1];
u3(-2.70020276220752,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.307536172365353,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.03942253299293,4.27575010230423,-0.598815541654238) q[1];
u3(1.63405178258829,0.328850372733665,0.824938497927861) q[3];
u3(2.11313100698769,3.25353402080905,-1.88121405139155) q[4];
u3(1.44302564502992,2.75398438642340,-2.86063166678672) q[2];
cx q[2],q[4];
u1(3.80524010927636) q[4];
u3(-1.41174841098684,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.11731286537135,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.50791142304470,-3.26709796049417,0.0164868382966235) q[4];
u3(1.75735367937495,-1.37711773675859,-3.23706569906061) q[2];
u3(1.19004643531082,2.29847080468149,-3.71687858564318) q[5];
u3(1.70291684819103,-2.58142384703644,2.73288595073658) q[0];
cx q[0],q[5];
u1(0.965703476183566) q[5];
u3(-3.27288169996652,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.77775899131093,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.352454513704516,3.73716648472067,-0.0479702199385326) q[5];
u3(0.686487323250422,1.22712007909946,0.824924123134491) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
