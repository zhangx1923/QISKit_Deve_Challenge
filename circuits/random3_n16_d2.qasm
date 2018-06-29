OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(1.69185361558571,2.11650738367645,-2.34563999394612) q[8];
u3(0.617380908105819,1.22077666340403,-2.35926039323758) q[12];
cx q[12],q[8];
u1(1.65190358099775) q[8];
u3(-3.29107591447826,0.0,0.0) q[12];
cx q[8],q[12];
u3(2.65332213835175,0.0,0.0) q[12];
cx q[12],q[8];
u3(0.249742317006215,1.36781213206840,-1.53469968411868) q[8];
u3(1.08869546743599,4.08200162857589,-1.44057524049316) q[12];
u3(2.12741544806097,0.977305179723529,-4.00800611794735) q[5];
u3(1.92678956136489,3.24478248617485,-2.44320663490869) q[15];
cx q[15],q[5];
u1(0.125161812933679) q[5];
u3(-1.67113330523749,0.0,0.0) q[15];
cx q[5],q[15];
u3(2.82443585352108,0.0,0.0) q[15];
cx q[15],q[5];
u3(0.317032233212180,1.73955145756982,-4.49641637473493) q[5];
u3(1.02317592828336,2.08853035203325,-2.93413931202086) q[15];
u3(2.52562394739224,-2.43943950048930,-0.545997339039995) q[11];
u3(1.74586490594316,-3.17482506308392,1.33210145023936) q[9];
cx q[9],q[11];
u1(2.81968001960631) q[11];
u3(-1.60519077021586,0.0,0.0) q[9];
cx q[11],q[9];
u3(-0.0163466413622586,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.29059907474730,-3.57228551028663,1.31725040518014) q[11];
u3(1.21379707926727,2.37518815284621,-3.30344599387359) q[9];
u3(1.32678234961636,3.77316535554836,-1.05922529228234) q[6];
u3(0.437461460014657,1.06586694774152,-1.78369919892210) q[3];
cx q[3],q[6];
u1(0.531318665646554) q[6];
u3(-0.962654164396614,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.05330459024727,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.45085951385784,-1.08102366368627,-0.884968449585634) q[6];
u3(1.93193851538697,0.0795877709092156,4.34016129909900) q[3];
u3(1.94805829649953,-0.778295734409523,1.99790969610740) q[1];
u3(1.78467071829683,-0.600983927990371,-0.556827296459506) q[2];
cx q[2],q[1];
u1(3.89724890239991) q[1];
u3(-4.32627362088365,0.0,0.0) q[2];
cx q[1],q[2];
u3(-1.08597701934119,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.397330470422475,0.388222775921655,1.05823718808493) q[1];
u3(1.23259587862052,0.441452473067762,-2.83483825883340) q[2];
u3(1.07912680026980,2.02836473206551,-2.86111122319995) q[14];
u3(1.41111501531988,-2.22074884346193,3.02500877327166) q[10];
cx q[10],q[14];
u1(2.96695338898325) q[14];
u3(-2.74746948012215,0.0,0.0) q[10];
cx q[14],q[10];
u3(1.05484137972807,0.0,0.0) q[10];
cx q[10],q[14];
u3(1.42825935405809,2.41652445825360,-0.992725148903837) q[14];
u3(2.47501352783555,-3.11850051244772,1.96670054323957) q[10];
u3(1.94684041982979,-0.158153366072815,-2.06706336644357) q[7];
u3(2.28630310244270,0.651689182154008,-4.35058844620255) q[0];
cx q[0],q[7];
u1(1.65875989628236) q[7];
u3(0.693783649956706,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.10597033797602,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.34562606687861,-1.45056101846074,4.29256014416674) q[7];
u3(2.42164577314136,-3.84748363440012,2.24135964058678) q[0];
u3(1.13406232732399,2.28357173137346,-1.87318342093545) q[4];
u3(2.17058895671724,0.491491326793905,-2.84525837510697) q[13];
cx q[13],q[4];
u1(1.42615452893035) q[4];
u3(-0.642795456046270,0.0,0.0) q[13];
cx q[4],q[13];
u3(1.89408479219741,0.0,0.0) q[13];
cx q[13],q[4];
u3(1.23125103030366,-0.281468099788470,0.388966562231926) q[4];
u3(1.04220773627744,-0.952679408088108,-3.85061375628200) q[13];
u3(2.02780811701042,3.03595286638922,-0.469126638284671) q[12];
u3(1.68568754276952,1.51726485840955,-1.34637366684520) q[8];
cx q[8],q[12];
u1(2.36845973863332) q[12];
u3(-1.69860353346085,0.0,0.0) q[8];
cx q[12],q[8];
u3(3.28887732476342,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.13666886485247,-1.72941219529164,3.06680623312731) q[12];
u3(2.14406516589925,-0.0116767759970053,3.94647754837508) q[8];
u3(0.882703489928949,-2.84143344960107,2.46225343524499) q[10];
u3(0.432168513782061,0.504099694278695,-2.96588914782967) q[6];
cx q[6],q[10];
u1(2.29760854231718) q[10];
u3(-1.84978771655300,0.0,0.0) q[6];
cx q[10],q[6];
u3(3.47607844399562,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.31572162806725,-3.55190149199423,0.259350183599850) q[10];
u3(2.51426066329062,-4.60176789149368,-1.07312572642512) q[6];
u3(1.51241638611913,0.507181119667129,-2.47330799127928) q[9];
u3(0.720698754499183,0.883298815258485,-4.01238397929927) q[1];
cx q[1],q[9];
u1(2.86211391263198) q[9];
u3(-1.75476425005724,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.655187481988689,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.42975861929081,0.396697973605080,1.31717995151207) q[9];
u3(2.35517451520028,-1.81747316376418,-3.56275880167654) q[1];
u3(2.19170592510861,-4.20467504503445,1.29270529570363) q[13];
u3(1.73309239914656,0.663834344542848,3.84525930076383) q[7];
cx q[7],q[13];
u1(2.57625731019379) q[13];
u3(-1.87222740780901,0.0,0.0) q[7];
cx q[13],q[7];
u3(3.33167239485673,0.0,0.0) q[7];
cx q[7],q[13];
u3(1.85845626238206,0.560037807021034,1.95476157846641) q[13];
u3(0.513802075786613,-1.45870013474630,-2.55628124487207) q[7];
u3(2.06384224288846,-0.405946090642492,1.34446629439121) q[4];
u3(2.36357564231432,-1.82740591325526,-2.88499584491758) q[3];
cx q[3],q[4];
u1(3.49658092391716) q[4];
u3(-1.10872552597300,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.49374096275169,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.67588995035169,1.65259935101037,-3.84295955163100) q[4];
u3(0.432977566292115,1.26649813048368,-2.56884290192893) q[3];
u3(2.35236377723472,0.717498665767786,-2.73743183738981) q[2];
u3(1.70573088287377,-3.12661075234250,3.10540714766303) q[5];
cx q[5],q[2];
u1(-1.39596977471499) q[2];
u3(1.06220831221621,0.0,0.0) q[5];
cx q[2],q[5];
u3(4.05055837368597,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.38903351386859,1.86606143037688,-0.732340731377681) q[2];
u3(0.524105627543201,-0.0690073207287640,-3.81802646773192) q[5];
u3(1.15773381536024,0.438927312274310,-2.47177738552502) q[11];
u3(1.40766687305424,-3.16880210937282,2.57288226497874) q[0];
cx q[0],q[11];
u1(1.51065656559814) q[11];
u3(-1.15896425516034,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.152362251868274,0.0,0.0) q[0];
cx q[0],q[11];
u3(2.44204771288086,-0.114936535007122,-0.290624432208984) q[11];
u3(1.36455246064205,2.36333647345491,3.80815517028117) q[0];
u3(1.46410250528874,-1.15867480282464,1.12284970706146) q[14];
u3(1.59112409750928,-2.10494642334102,-0.786849542838118) q[15];
cx q[15],q[14];
u1(1.84588571773583) q[14];
u3(-0.0353190553264544,0.0,0.0) q[15];
cx q[14],q[15];
u3(1.00468701813974,0.0,0.0) q[15];
cx q[15],q[14];
u3(0.437446760970227,3.04427129094793,-1.81430038540966) q[14];
u3(2.33764253396071,4.48562547564184,0.343494438129488) q[15];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];