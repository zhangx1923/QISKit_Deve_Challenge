OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.08070466336600,-0.110644784239953,2.03648450943722) q[2];
u3(0.977110388378125,-2.78334036485502,-1.70317443991652) q[0];
cx q[0],q[2];
u1(0.360994986722734) q[2];
u3(-1.03468893320066,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.50324017793384,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.872501017726120,0.542113272156876,2.88205305177650) q[2];
u3(0.903789501675370,2.59899799537241,-1.76803564197138) q[0];
u3(1.56555764354423,2.78224113260527,-0.490292693126464) q[3];
u3(1.46004664479136,1.48949864958968,-1.95042893529971) q[1];
cx q[1],q[3];
u1(2.33019338856850) q[3];
u3(-2.56406133138349,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.656855258729650,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.95167350716914,1.92161993563233,-2.62302874826299) q[3];
u3(0.804226147999524,-5.19775227572786,-0.592969587361683) q[1];
u3(1.14653437583058,-1.34294864613160,1.64504197431184) q[1];
u3(0.326114562831253,2.23796059223754,-3.72419284778773) q[4];
cx q[4],q[1];
u1(1.32485592210163) q[1];
u3(-3.11312395966183,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.364464954742347,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.38915343187116,0.114636225964027,3.55537297014970) q[1];
u3(0.630786212209966,-2.45106906028634,-1.25395105305526) q[4];
u3(1.75900824285810,2.12589260486323,-2.81034224818201) q[2];
u3(0.785348881639804,2.43811610647057,-1.49837622739909) q[3];
cx q[3],q[2];
u1(-0.177700754526265) q[2];
u3(0.815473184603850,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.80135271558480,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.50748198413373,-2.39947669162507,3.17775644435077) q[2];
u3(2.03677229251654,0.531494776966007,-3.00231996241717) q[3];
u3(0.0817472637110536,-1.60630004766286,2.70117187627894) q[4];
u3(1.07197823389210,-2.24374390955071,0.919550265024740) q[3];
cx q[3],q[4];
u1(2.92990383360166) q[4];
u3(-2.54981332053786,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.24490930191513,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.223480801165880,1.15761137414610,0.464767571088269) q[4];
u3(1.29222470717934,1.59109046564669,-3.32010975426636) q[3];
u3(2.99169708633436,2.03598691597162,0.294216106217752) q[1];
u3(2.77766226622154,0.638654660030140,-3.93960215414450) q[0];
cx q[0],q[1];
u1(2.59727014094675) q[1];
u3(-2.93594707686552,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.455734525205431,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.78604231430157,0.562449413531646,-1.49122880589830) q[1];
u3(2.63929208662955,4.57029966019927,1.56608431585443) q[0];
u3(2.34261459923030,-1.08422078629375,2.28517936397820) q[0];
u3(2.36649364993976,-0.642962471413501,0.0489328212449677) q[4];
cx q[4],q[0];
u1(4.31125482030129) q[0];
u3(-3.92309311478475,0.0,0.0) q[4];
cx q[0],q[4];
u3(-1.07716324970875,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.52797308249688,-0.455884634844045,-3.08655243747375) q[0];
u3(1.65345153701594,-0.358095019051946,5.80213101974858) q[4];
u3(0.883846774402787,2.48707727095444,-1.18224497067362) q[3];
u3(1.28025073986526,1.55473227145553,-0.300284911211305) q[1];
cx q[1],q[3];
u1(3.59422061345056) q[3];
u3(-0.807944283287618,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.56608101774590,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.992691927236591,2.78437398171014,0.0171838955714272) q[3];
u3(2.37455669603782,-0.875191649728212,1.38898772685897) q[1];
u3(2.02585546267824,1.44597522270838,-4.43674143492479) q[2];
u3(0.791127975627687,-2.05320722011864,3.23134124184189) q[0];
cx q[0],q[2];
u1(2.15485830905614) q[2];
u3(-3.01601531309753,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.35841722436308,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.78856765295514,-1.67725805975447,-0.894693640276243) q[2];
u3(1.80632901777109,-2.75384401277491,2.14079416743217) q[0];
u3(2.13696025015656,-0.233703188240765,0.972325181956617) q[1];
u3(1.84334246056414,-1.29309313388725,-2.37274991081536) q[3];
cx q[3],q[1];
u1(2.17032381689798) q[1];
u3(-3.12690673775543,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.33910019645070,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.98141742708831,-2.02860885299577,3.52730585169095) q[1];
u3(1.69874313226971,2.22738817380244,1.54401612090723) q[3];
u3(2.31805994586561,1.92023879795861,-2.85389019979156) q[2];
u3(1.73186633348463,2.46444377764814,-2.97318994078511) q[1];
cx q[1],q[2];
u1(1.12540212497774) q[2];
u3(-1.52898096119792,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.119164449787298,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.48207221520726,0.198016389284058,-1.85930211885079) q[2];
u3(0.320359982728345,-2.43946120004760,1.18303795468464) q[1];
u3(2.17970687464256,-0.842337060474209,-0.441560439721462) q[4];
u3(1.24384227451786,-2.61284965465530,-0.810116795581415) q[3];
cx q[3],q[4];
u1(3.55059196906474) q[4];
u3(-4.31985205972065,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.681505511864634,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.17557196489980,-3.92410557662815,-0.653783942084706) q[4];
u3(2.45397590781742,1.88016096012356,-1.23594286259466) q[3];
u3(2.13583487397079,0.699231193668440,-2.05157997909726) q[0];
u3(2.11362809957585,3.99819523389937,-0.156780106435272) q[2];
cx q[2],q[0];
u1(3.23854554703567) q[0];
u3(-1.31004440675044,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.91171187728148,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.52516806534047,-3.90063510741802,1.89309380617361) q[0];
u3(0.717585301561104,-0.439617938937086,-2.93235943722351) q[2];
u3(2.79908235717190,-0.260058879761162,1.86042882492058) q[3];
u3(2.28853171346990,-1.83506773015030,-0.755324241117980) q[4];
cx q[4],q[3];
u1(2.70049501684980) q[3];
u3(-3.18458126455623,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.940636793465024,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.02774736613185,1.84322380033814,0.441236330439319) q[3];
u3(1.91169282601952,0.741980778730331,-1.78260484095863) q[4];
u3(1.02808513472198,1.73900287703837,-2.91466015120529) q[4];
u3(1.49692579906235,-2.09883227978239,3.50689942798261) q[2];
cx q[2],q[4];
u1(0.799514069679963) q[4];
u3(-0.218103905053961,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.24957441850001,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.871728600046254,0.143137562385692,-3.71347995753910) q[4];
u3(1.08323191511627,2.85181786302962,0.0669298693033671) q[2];
u3(1.88917390437245,3.71287419670273,-1.16451840493511) q[3];
u3(2.20097792235188,0.981420801231965,-1.36560232239395) q[1];
cx q[1],q[3];
u1(-1.15993473272535) q[3];
u3(0.484200447274448,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.33363490090313,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.45862484442186,-3.87678124929945,2.24916049288615) q[3];
u3(0.877654526079292,-3.68373415855246,1.84306576416984) q[1];
u3(1.30891158435368,-1.89775244182452,-0.0565574724076239) q[0];
u3(1.05787544978504,-2.48213253505190,0.937448314194650) q[3];
cx q[3],q[0];
u1(2.39611545046397) q[0];
u3(-2.67290669452497,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.46995161269206,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.49059877097735,-0.300603116661393,1.76888662454809) q[0];
u3(1.33069744052460,0.0181140536072111,3.47752545351572) q[3];
u3(0.282925771361282,-2.44146168781214,2.31791304300649) q[1];
u3(0.938572851150805,0.0314254343466809,-2.73833792810734) q[4];
cx q[4],q[1];
u1(1.04689633363314) q[1];
u3(-2.96622261582823,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.83045011043207,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.32276660059624,1.70143949979799,-3.33763571334974) q[1];
u3(2.09676653043394,1.26531372689259,-1.38393551810591) q[4];
u3(1.68791045627023,-0.954678486369356,-0.858927422946434) q[4];
u3(2.31982090932712,-2.75538967664413,-0.0879342177905131) q[1];
cx q[1],q[4];
u1(3.74169250824990) q[4];
u3(-4.33552162873269,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.510665419355981,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.78179990080282,2.02250292672652,0.926296389970377) q[4];
u3(1.46617967984654,1.09795881905740,-0.217311336352223) q[1];
u3(2.27519795812969,1.32775480562187,-2.28743182481458) q[3];
u3(1.62774240074313,-2.10519403241220,2.44729581629326) q[0];
cx q[0],q[3];
u1(0.658608259414272) q[3];
u3(-1.73388295049472,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.303381696978613,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.688708609823125,-0.0351795805860778,0.460692031166611) q[3];
u3(1.29065490100049,-2.44501049153420,-0.112842762280515) q[0];
u3(2.05234393121276,-1.42344153999835,0.454517885020357) q[1];
u3(1.98175589896365,-1.92286497578504,-0.211991204021908) q[2];
cx q[2],q[1];
u1(3.51988065475518) q[1];
u3(-1.11851748896442,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.15251203387529,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.837772453491820,4.16490644700396,-1.91596337340674) q[1];
u3(0.998690123732936,-4.06186898694608,-1.08031910550346) q[2];
u3(0.721052997306657,-2.12147870869881,2.51401564662438) q[4];
u3(0.645356464478012,-2.80181017378807,0.250591798065762) q[0];
cx q[0],q[4];
u1(1.65574315241557) q[4];
u3(-0.562498982898282,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.199746598320100,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.62170076348197,-2.27581030129409,0.332680042644322) q[4];
u3(1.74032632621526,3.39303612729750,1.16623311678092) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
