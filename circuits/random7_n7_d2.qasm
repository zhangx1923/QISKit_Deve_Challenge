OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.71517193389626,-2.52488766107839,-0.464394313511005) q[5];
u3(2.22606582321899,-3.35521970821416,-1.32277528045748) q[3];
cx q[3],q[5];
u1(1.83575694901989) q[5];
u3(-2.35660494785556,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.571735511620796,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.43077103980226,0.464060533424156,-1.30831544042205) q[5];
u3(2.40434410091162,2.66583307039486,0.135358439773583) q[3];
u3(1.73785445134650,1.40839026263568,-4.26776697652727) q[4];
u3(0.563061088926674,-1.61201626754483,2.98860063107481) q[2];
cx q[2],q[4];
u1(0.911014898668783) q[4];
u3(-0.441097697749417,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.59305638440376,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.804043851444659,-1.89305001199071,2.85026164585817) q[4];
u3(1.24677203635446,2.54928675608603,-0.804814465741375) q[2];
u3(2.70862890459912,0.916559157742631,-2.11904040379141) q[6];
u3(2.02328087621898,4.00210413217471,0.274513929294069) q[0];
cx q[0],q[6];
u1(2.64797893703910) q[6];
u3(-2.96338488042880,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.824520700315745,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.97303722472105,0.944544662596395,-1.59134102394428) q[6];
u3(2.20561902718222,-0.148128032275267,3.01375267698714) q[0];
u3(1.90888687178949,-2.77945377283024,0.135880927948222) q[6];
u3(2.34215617383173,-3.56223720754600,-1.20149760156404) q[1];
cx q[1],q[6];
u1(0.459191414428455) q[6];
u3(-1.51647489821502,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.09001473720302,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.27816964965095,1.94642777981734,-3.34156616358425) q[6];
u3(1.63459707505922,3.55144210827485,-1.67894958208530) q[1];
u3(0.485739252256328,-2.02041959995421,1.63133388100691) q[2];
u3(0.727373159398587,2.58018167716076,-3.61278972941952) q[5];
cx q[5],q[2];
u1(3.35922765673944) q[2];
u3(-1.42170029896774,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.38196728486444,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.329186143816381,-0.926445430274538,3.71751871382261) q[2];
u3(2.02903444378800,0.513990276550635,3.75880381329883) q[5];
u3(1.71876682878241,-1.37273794397201,-1.04749526422942) q[4];
u3(1.24281676408314,-4.44830390953017,0.846928280523936) q[3];
cx q[3],q[4];
u1(1.69348523372917) q[4];
u3(0.570280330813620,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.916198266151387,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.25427289580049,2.02697826006665,-3.66959752538215) q[4];
u3(0.519254408652977,-0.168612642930295,0.383628974028507) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
