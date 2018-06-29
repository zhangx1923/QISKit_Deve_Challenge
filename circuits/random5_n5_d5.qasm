OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.42538077194642,0.149806910012440,-1.26591298086605) q[2];
u3(1.70227990840078,0.925466794284366,-4.58645209135050) q[1];
cx q[1],q[2];
u1(3.43814855178034) q[2];
u3(-0.651105291664774,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.68191650172843,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.473295627222671,-1.65024182331622,1.11172388233936) q[2];
u3(1.99032282379196,-1.72704060057471,-1.74295557294033) q[1];
u3(1.25052875075087,3.56857287919479,-0.478055445186474) q[3];
u3(0.859413518304346,1.79515457714528,-1.74129518007634) q[4];
cx q[4],q[3];
u1(3.38572716185210) q[3];
u3(-1.11568847400622,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.74879100866968,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.42275363475215,-0.844147147276678,-1.77650513690763) q[3];
u3(1.25567423695989,3.48891106892302,0.855403696685345) q[4];
u3(0.478884562065193,0.764949081770533,-0.578790339865441) q[0];
u3(0.665491213838843,-2.24428139670662,1.33232909097877) q[1];
cx q[1],q[0];
u1(2.69278804920576) q[0];
u3(-2.05571175030529,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.43205579980661,0.0,0.0) q[1];
cx q[1],q[0];
u3(3.04343552645911,-1.58930263718757,2.28607363339019) q[0];
u3(1.71238907381504,1.37638365805417,0.424480525164089) q[1];
u3(1.61466629539310,1.13590859233039,-1.16178254135061) q[4];
u3(0.365685962982345,-0.203095161862536,-2.35791874808928) q[2];
cx q[2],q[4];
u1(1.29224590620269) q[4];
u3(-0.815561689277560,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.55427375011600,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.85875922654099,0.00956862673154546,-0.993497143056310) q[4];
u3(1.62659898100296,-0.925668660240118,2.17327776215491) q[2];
u3(0.882268205833222,0.435498618694641,-2.33347141101598) q[3];
u3(1.90909141138885,-2.68890886913241,2.73475488014654) q[2];
cx q[2],q[3];
u1(0.831220388182860) q[3];
u3(-3.22221286129035,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.50855654758223,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.37401508860754,1.61458097991612,-1.00035317331566) q[3];
u3(1.39972316037690,-2.65138650808259,-1.85098072787703) q[2];
u3(1.72174321783524,0.965957832870328,1.85718892047884) q[4];
u3(2.00390628283733,-2.02363577886556,-1.40946210747335) q[1];
cx q[1],q[4];
u1(2.24357538138625) q[4];
u3(-3.21040983668383,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.70484452326362,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.03258918895441,-1.03685386222510,1.11829975806483) q[4];
u3(1.96657142624647,3.92480889355425,0.611734837215654) q[1];
u3(1.43862556037628,-1.39880245493778,-1.00884815354241) q[1];
u3(2.13154095853938,-4.74200793467237,1.22520976422926) q[2];
cx q[2],q[1];
u1(1.25630256232190) q[1];
u3(-3.12546946393737,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.70833101975813,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.881406551719904,-1.66082446742859,3.79797796047873) q[1];
u3(1.90707646739551,-3.10595205134562,-0.422778811051819) q[2];
u3(0.933995964840327,2.21246993987889,-0.503976952493724) q[3];
u3(2.07039377152767,1.30616142924861,-1.43785119891479) q[4];
cx q[4],q[3];
u1(0.484713809017929) q[3];
u3(-1.13031968402968,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.63134889124743,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.16286904249348,1.46112647232502,-0.367864094022520) q[3];
u3(1.95652125541858,2.16820486873801,4.02902642041018) q[4];
u3(1.14971579313096,-1.80107675990589,0.263174521025255) q[4];
u3(1.45048848415536,-4.31221983355152,0.429595953879445) q[1];
cx q[1],q[4];
u1(1.79929232678126) q[4];
u3(-2.98899936696403,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.46015299974710,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.41802029495947,-1.78746982647365,1.03805155770774) q[4];
u3(1.13850172959679,0.753660135955959,0.714454032885850) q[1];
u3(2.75098057484254,-0.554998688419745,-1.21585109653590) q[2];
u3(1.00958675006072,0.473928353752147,-4.70161854595334) q[0];
cx q[0],q[2];
u1(0.837730414653689) q[2];
u3(-0.449710795948805,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.42516580928842,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.51008203174673,4.37533296498178,-0.495175738785704) q[2];
u3(1.37037008711338,-2.33542373632522,-0.219545185399112) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
