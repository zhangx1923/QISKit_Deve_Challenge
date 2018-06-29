OPENQASM 2.0;
include "qelib1.inc";
qreg q[15];
creg c[15];
u3(0.189693191782773,-1.77755839315911,1.31723589967900) q[8];
u3(1.00909956044427,-0.212808615426816,-1.65097560041628) q[10];
cx q[10],q[8];
u1(2.64624882678015) q[8];
u3(-2.18543034683303,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.321588987017688,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.71787049031554,1.67519444714786,0.479025203632686) q[8];
u3(0.521506853152290,3.06420807487912,1.59046310317028) q[10];
u3(2.51869287148095,-3.19642202727186,2.93969050611660) q[13];
u3(0.614611491289680,3.20229870781754,-1.99886207351066) q[14];
cx q[14],q[13];
u1(1.58601069382486) q[13];
u3(-2.53943099062529,0.0,0.0) q[14];
cx q[13],q[14];
u3(3.28014400634768,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.41990764017055,0.160852597866427,-0.808199709585344) q[13];
u3(1.47308183397504,-0.0449755841173056,-4.37388568552685) q[14];
u3(1.46637609694495,1.47713530641755,-2.70845137522158) q[3];
u3(1.39600878270597,-4.80491861537134,1.36075603185596) q[11];
cx q[11],q[3];
u1(-0.241488948016009) q[3];
u3(-1.55006344248043,0.0,0.0) q[11];
cx q[3],q[11];
u3(0.472786680249339,0.0,0.0) q[11];
cx q[11],q[3];
u3(0.700749514398964,2.00867279782242,-3.80035824346793) q[3];
u3(1.63075789613611,2.83986117315537,-0.346836729672699) q[11];
u3(2.46134612732177,-1.20837036637861,-1.18951807635781) q[2];
u3(1.30218499205684,-1.59225236854991,-3.55532184749615) q[6];
cx q[6],q[2];
u1(0.772769782746516) q[2];
u3(-1.08858833661993,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.0701070902134597,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.727038265018320,0.457689862575223,-1.21836725009969) q[2];
u3(2.01096083731203,-4.12497753770120,-0.210548089342751) q[6];
u3(2.44403726847349,-0.440131114905870,-1.21616894295387) q[7];
u3(0.613313183987813,-0.490062263878688,-4.31159216780613) q[1];
cx q[1],q[7];
u1(1.44741867677858) q[7];
u3(-2.41853172221731,0.0,0.0) q[1];
cx q[7],q[1];
u3(3.28671512051661,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.12962329547917,2.51764332899219,-3.41103609898363) q[7];
u3(1.70911592232501,1.76813741114930,3.71630561798461) q[1];
u3(2.96618240408909,0.411506611588344,0.147333313791841) q[12];
u3(0.540497250785191,-5.44402268157803,0.705793870964089) q[0];
cx q[0],q[12];
u1(0.654343301849350) q[12];
u3(-0.205105861865363,0.0,0.0) q[0];
cx q[12],q[0];
u3(2.00395920626205,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.30491331673001,0.774214942476277,-3.87342224557057) q[12];
u3(0.923325583685888,0.0186865191238084,-2.65630179322683) q[0];
u3(2.69028791592037,2.69882531268789,-1.80902512276839) q[9];
u3(1.77809999179078,1.81280680950542,-2.26164543674004) q[5];
cx q[5],q[9];
u1(1.44618486210547) q[9];
u3(-3.65156375741670,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.31599662964990,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.92476302410700,2.88009261656469,-1.20225428774773) q[9];
u3(2.16216450359247,-5.33941528723970,0.159926456072006) q[5];
u3(1.66459048227142,3.56486866814920,-1.37688360711363) q[9];
u3(2.15114732544095,2.18432167712876,-2.16988491439407) q[8];
cx q[8],q[9];
u1(0.0918648178305688) q[9];
u3(-0.973171932116662,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.47170579970703,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.98968643919165,-1.97345031831851,3.56395153376810) q[9];
u3(2.73772640753898,1.86498143111365,-3.24643011294718) q[8];
u3(1.39850553866177,-0.575086005965562,-0.779503178137703) q[2];
u3(0.950476669239496,-2.78271356856767,0.254436432520205) q[6];
cx q[6],q[2];
u1(1.88201219196061) q[2];
u3(-3.00822715073283,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.65824063823690,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.90644299183331,1.46616185996376,-1.20587054530712) q[2];
u3(1.28149389493504,-4.72208910499813,-0.310717906960690) q[6];
u3(2.61729833869026,0.593214263064718,0.373817753266068) q[4];
u3(1.06522021253760,-3.64800373920149,-0.200820837521631) q[10];
cx q[10],q[4];
u1(1.96137293244633) q[4];
u3(0.214547974993465,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.837849558881548,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.46466914732323,-1.52992253357298,1.44617250972594) q[4];
u3(1.01479024847879,-4.39140975088818,-1.58404354594740) q[10];
u3(1.55738275710141,0.387251720812559,1.70314688752840) q[13];
u3(1.65734485486501,-2.68798776864513,-1.93359307385734) q[14];
cx q[14],q[13];
u1(0.364874250868593) q[13];
u3(-1.45364589201517,0.0,0.0) q[14];
cx q[13],q[14];
u3(2.62630345606280,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.69483613307654,-2.52564385487666,-0.0164816222691471) q[13];
u3(1.21245117973404,3.70375485434544,-0.358127984615310) q[14];
u3(0.0988123379296390,2.90924715734443,-1.75395332650473) q[11];
u3(0.914054423783070,-2.13233756103816,0.634203844353201) q[3];
cx q[3],q[11];
u1(1.57098437366629) q[11];
u3(-2.49878222151366,0.0,0.0) q[3];
cx q[11],q[3];
u3(3.30389014125122,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.78648091868817,-1.64424073043624,-1.26591338308395) q[11];
u3(1.00171008578608,4.83482168276364,-1.41055729440471) q[3];
u3(1.58349512908447,-0.0380182720021027,2.75690867853278) q[0];
u3(1.71123842273599,-1.92224732178741,-1.83451694489090) q[5];
cx q[5],q[0];
u1(1.42648757937691) q[0];
u3(-0.727608445443774,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.16193307154142,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.25353633625404,-2.24454670401671,-2.04109158510191) q[0];
u3(0.962320768814825,2.98829338707379,1.43796738384694) q[5];
u3(1.73456314030184,3.69851195078800,-2.24743992517467) q[1];
u3(2.51426277875479,1.45988083913102,-1.98337170558199) q[12];
cx q[12],q[1];
u1(2.02207507436234) q[1];
u3(0.507957190219453,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.80914206003999,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.98075870102789,-0.0422663610285079,-1.25845931454974) q[1];
u3(1.34961586563101,-4.50687518304113,-0.363970957677004) q[12];
u3(0.600676039547410,-1.82251829979907,-0.656545636028684) q[3];
u3(0.984071417106516,-3.59380772538456,-0.0185687857201660) q[11];
cx q[11],q[3];
u1(0.553573767405825) q[3];
u3(-1.62601332434339,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.94921343984071,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.85533809596339,-2.12010680161971,1.59229601156687) q[3];
u3(2.48480575221635,2.06971488910257,-3.56122774911476) q[11];
u3(2.35237823491289,1.46938255369459,-2.60286393099495) q[4];
u3(1.20429319212206,2.54325578157222,-3.71066264239402) q[9];
cx q[9],q[4];
u1(0.0365085475125122) q[4];
u3(-1.67013122856688,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.556471979100774,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.14457164790652,0.570292640278878,-0.0141662410010752) q[4];
u3(1.73657542588726,2.29441189189391,2.45150646159521) q[9];
u3(1.49614345528096,1.68933776768964,-2.81463295548290) q[10];
u3(1.99699364686165,-1.98309654464516,3.35566769161870) q[1];
cx q[1],q[10];
u1(1.80871247206439) q[10];
u3(-2.48608542014282,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.107954443587267,0.0,0.0) q[1];
cx q[1],q[10];
u3(2.65158423481060,-0.569940290063230,1.08513044436242) q[10];
u3(0.833536922988262,-0.625824833379641,3.77831566024622) q[1];
u3(1.11861541116236,1.48134823831809,-2.77302303157756) q[5];
u3(1.71602497704860,-3.62614794778499,2.43269951873355) q[6];
cx q[6],q[5];
u1(3.08567552170880) q[5];
u3(-1.30861984502640,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.21788947565973,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.28634414736203,-1.83423739310654,2.42008479093283) q[5];
u3(1.58662617823113,0.104432600671829,4.69360328415126) q[6];
u3(1.42839249623114,0.520547363936262,-2.02758277234979) q[2];
u3(2.52049082673882,-3.60744023586480,2.60390030043186) q[7];
cx q[7],q[2];
u1(1.48187459915288) q[2];
u3(-0.217172438507297,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.39275060946802,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.803901667820742,-1.00426532037426,1.53660316857989) q[2];
u3(2.08262656851897,-1.84618442304931,-3.18227000635274) q[7];
u3(1.55644016358425,-2.03979204119719,-0.354555082318340) q[13];
u3(1.79902093342846,-4.45999869490190,-0.943647772671844) q[0];
cx q[0],q[13];
u1(0.369396393439950) q[13];
u3(-0.956542889952109,0.0,0.0) q[0];
cx q[13],q[0];
u3(2.58133038199197,0.0,0.0) q[0];
cx q[0],q[13];
u3(2.90898230386230,2.83399519088180,-0.820510442993537) q[13];
u3(0.681930794280604,3.74997534848623,1.61035749035357) q[0];
u3(2.60868559973777,2.32700639648296,-3.50447076023561) q[8];
u3(0.992100819645034,3.06819859822207,-2.86326064558728) q[12];
cx q[12],q[8];
u1(1.11131111187958) q[8];
u3(-0.631222569829921,0.0,0.0) q[12];
cx q[8],q[12];
u3(2.80350387293646,0.0,0.0) q[12];
cx q[12],q[8];
u3(0.962059439322576,-2.45937539526766,0.912475105267540) q[8];
u3(1.76255116330006,-3.48170442072452,0.954131219354912) q[12];
u3(1.37968569305932,-0.0665158716216303,2.60624156359156) q[5];
u3(1.29173892307233,-1.05946241011200,-1.63548086945812) q[12];
cx q[12],q[5];
u1(0.156249740308526) q[5];
u3(-1.41217879618698,0.0,0.0) q[12];
cx q[5],q[12];
u3(1.23006068906915,0.0,0.0) q[12];
cx q[12],q[5];
u3(2.27911977736723,-0.328575516229829,1.78117698402428) q[5];
u3(1.89046144863914,-1.32789265097937,1.79203758350543) q[12];
u3(0.884677335277089,0.462045588667196,0.828939498043250) q[3];
u3(1.23598006416573,-0.488840063086320,-2.00829128692207) q[6];
cx q[6],q[3];
u1(2.14881800402383) q[3];
u3(-2.56420484396472,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.36093425293355,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.53647975021128,1.00194664295160,0.340916399373300) q[3];
u3(2.17502037040834,-2.46610248096413,-2.03051827293403) q[6];
u3(1.73137549623369,-1.82731328181934,-0.402228652112974) q[8];
u3(2.56055535968870,-2.79153596810002,-0.733816484380809) q[14];
cx q[14],q[8];
u1(0.662671943945430) q[8];
u3(-0.452254331431678,0.0,0.0) q[14];
cx q[8],q[14];
u3(2.20681910398470,0.0,0.0) q[14];
cx q[14],q[8];
u3(0.179516301083130,-1.98250519736094,-1.09040265883535) q[8];
u3(2.33891858294801,1.23484888055529,0.373279786642819) q[14];
u3(0.982758793791667,2.65976328639953,-3.16351585874026) q[4];
u3(0.780524955647461,2.88514389997731,-2.92916469285582) q[13];
cx q[13],q[4];
u1(0.516477166950317) q[4];
u3(-3.09517653299352,0.0,0.0) q[13];
cx q[4],q[13];
u3(1.53689007829412,0.0,0.0) q[13];
cx q[13],q[4];
u3(0.711233252448271,1.36064096597512,-0.484653418763263) q[4];
u3(1.58884775882760,-2.95111099779214,-0.444580425343891) q[13];
u3(2.49043273394602,-1.10756011859922,-1.62092287238192) q[9];
u3(0.732537100160075,-1.57522795356716,-2.86025223920355) q[2];
cx q[2],q[9];
u1(1.34027426967857) q[9];
u3(-0.726999207209376,0.0,0.0) q[2];
cx q[9],q[2];
u3(-0.297280025664235,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.39406751560793,-1.96516284793705,-0.355962708095427) q[9];
u3(1.20202539828872,-1.70198128296077,2.97021823129111) q[2];
u3(1.89101934002235,1.17200456783556,1.48455142899950) q[10];
u3(1.50982248465693,-1.57705213354312,-1.59742548607438) q[0];
cx q[0],q[10];
u1(-0.190636788778359) q[10];
u3(-1.50741935689113,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.749312844583650,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.78933449439365,2.37042885371734,-2.57590859375069) q[10];
u3(1.48756555588936,2.31821076044151,-2.75778656632223) q[0];
u3(0.879419096513511,2.67862470907051,-0.124504727878627) q[11];
u3(1.62594502298687,1.56175527829228,-1.69636220198536) q[7];
cx q[7],q[11];
u1(1.00209113459302) q[11];
u3(-1.36954517087167,0.0,0.0) q[7];
cx q[11],q[7];
u3(-0.111261254193610,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.50774039408886,0.629972425745717,-3.76476218985948) q[11];
u3(1.69734632815881,1.37441727334743,-0.675692444314142) q[7];
u3(1.96647207629976,1.47520688591423,-2.77240802455649) q[10];
u3(1.58951757788655,-1.88402343670733,2.03303653133895) q[2];
cx q[2],q[10];
u1(-1.36531753367545) q[10];
u3(0.241087099043299,0.0,0.0) q[2];
cx q[10],q[2];
u3(3.53575393700740,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.12466882046969,-2.08666551228026,1.39383126240351) q[10];
u3(0.771845383920386,2.01101723172107,-4.05313126828862) q[2];
u3(3.01394069155958,0.767245960871779,0.0243529968382606) q[6];
u3(1.48796789091130,-0.748590404228011,-3.39529314099553) q[9];
cx q[9],q[6];
u1(1.68938915010027) q[6];
u3(-2.39700353453598,0.0,0.0) q[9];
cx q[6],q[9];
u3(3.17522416031372,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.663621826885270,-1.01443599031117,4.86530407421374) q[6];
u3(1.74808944946989,-0.587294185804780,0.162737562676502) q[9];
u3(1.40036625508558,1.51198514244016,-2.71475856430141) q[1];
u3(2.06640094619493,-3.35501153430106,2.60643630687927) q[12];
cx q[12],q[1];
u1(1.41705668655723) q[1];
u3(-3.28603713085536,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.54413500852019,0.0,0.0) q[12];
cx q[12],q[1];
u3(0.661800971688349,-1.27600459663725,0.885030584757966) q[1];
u3(2.49302327794598,-3.43757171977510,0.181494369301188) q[12];
u3(1.34268307441639,-0.172662257284668,0.577798610560219) q[13];
u3(0.711844085143567,-2.46923158668350,-1.13997565261353) q[3];
cx q[3],q[13];
u1(3.35098093615318) q[13];
u3(-1.39687462207860,0.0,0.0) q[3];
cx q[13],q[3];
u3(2.48231752072843,0.0,0.0) q[3];
cx q[3],q[13];
u3(0.485434126266176,-0.862080754974081,-0.295332931854171) q[13];
u3(1.82207446284066,2.34475444978697,1.97703763401169) q[3];
u3(1.33380045389864,3.53006123463318,-1.70579970296540) q[7];
u3(0.389504691035544,1.73765752464504,-2.71115158846820) q[4];
cx q[4],q[7];
u1(0.115182507574746) q[7];
u3(-0.912707905674791,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.66407175526816,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.38890200245395,-2.00741296957888,2.75317631992407) q[7];
u3(2.19604747876799,-2.29785395194476,1.83223481466960) q[4];
u3(2.74523835313225,1.86141783779722,-2.68007460076711) q[14];
u3(1.44923893027277,1.93976863000915,-3.33250421539897) q[11];
cx q[11],q[14];
u1(1.90901112542356) q[14];
u3(-2.30134293080202,0.0,0.0) q[11];
cx q[14],q[11];
u3(3.52202756620166,0.0,0.0) q[11];
cx q[11],q[14];
u3(0.359602759417134,0.237488242756600,-3.51698114596948) q[14];
u3(1.45177179691048,-2.74612796174647,1.83105050128641) q[11];
u3(1.76908665090573,1.88071234919413,-0.216966794175300) q[8];
u3(0.764252691225924,-0.393022867155286,-3.60198699762078) q[0];
cx q[0],q[8];
u1(1.53270880674875) q[8];
u3(-2.52554388433781,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.06785633885323,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.349242273303719,3.00393701833286,-1.03326005245296) q[8];
u3(0.608534824397881,1.37417081571735,-4.86281820067034) q[0];
u3(1.50886836571260,4.24846548223004,-1.85153453252040) q[6];
u3(1.38612269631366,2.05283727028595,-2.31103318633180) q[3];
cx q[3],q[6];
u1(0.701700019512911) q[6];
u3(-1.48351839230033,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.79322628087105,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.19513384026926,-2.39322088463409,0.0477555526518536) q[6];
u3(1.61660869628694,-3.56987611529646,-2.61799830447498) q[3];
u3(0.482548447126318,2.66127191647026,-2.97690806397871) q[13];
u3(0.644930969299063,0.267447145310157,-0.812611943529653) q[11];
cx q[11],q[13];
u1(1.78129282022417) q[13];
u3(-3.13770773600779,0.0,0.0) q[11];
cx q[13],q[11];
u3(0.381104262837483,0.0,0.0) q[11];
cx q[11],q[13];
u3(2.82843180285552,-0.360774544861017,0.555495717429095) q[13];
u3(2.66184099786742,2.23836555796305,-0.841749944148845) q[11];
u3(1.87401498458564,-1.32700352562592,4.23865542858565) q[7];
u3(0.575216609612956,-1.70927948636778,3.48362908201108) q[12];
cx q[12],q[7];
u1(1.55038418609090) q[7];
u3(0.190567430879855,0.0,0.0) q[12];
cx q[7],q[12];
u3(0.571381607817611,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.05291536867770,-1.54048559342610,4.07861867128737) q[7];
u3(1.00598598066779,-1.74566947062619,-3.42230050254529) q[12];
u3(0.187277085845031,1.99208907385053,-1.68833064525501) q[10];
u3(1.07299230577560,-3.16110252439427,0.369824546512231) q[1];
cx q[1],q[10];
u1(1.61801328509474) q[10];
u3(0.285203173863602,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.689874732444222,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.22207480643889,2.07254315230596,-2.00874325680730) q[10];
u3(2.68791936087423,-0.000147100071011597,3.28609740251273) q[1];
u3(1.84842214945591,0.714979349341262,-3.39797277069585) q[2];
u3(2.03314373488739,3.13404199439326,-2.20893399807045) q[4];
cx q[4],q[2];
u1(0.585164399532272) q[2];
u3(-1.31268977735912,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.180680843018004,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.07039472936841,-2.86504067239320,2.61973286186805) q[2];
u3(1.57896906994134,-3.82710581739982,1.36759495558505) q[4];
u3(1.60276103546590,1.85649844917746,-2.85636989291563) q[14];
u3(0.721609752298637,2.23299104178377,-2.91231033946200) q[8];
cx q[8],q[14];
u1(3.14992843818524) q[14];
u3(-1.27148363644089,0.0,0.0) q[8];
cx q[14],q[8];
u3(2.08194333988852,0.0,0.0) q[8];
cx q[8],q[14];
u3(2.15972911290964,-1.74932888284360,0.0957552060297332) q[14];
u3(1.29938613117145,-0.131593885749578,-1.29444791225020) q[8];
u3(2.09112097245681,3.88433893330430,-1.49804084113088) q[9];
u3(2.37497928003829,2.56412908351506,-0.0666044769575345) q[0];
cx q[0],q[9];
u1(0.981279827582509) q[9];
u3(-0.725083665658992,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.86763081309024,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.49201840305499,2.94719002055752,-3.22568459368679) q[9];
u3(1.75802651668774,-3.76375868704928,0.0199618389379219) q[0];
u3(1.37860715177774,0.540831312642419,2.37460019655909) q[1];
u3(1.78998363809943,-3.61843901096184,-2.55854411460098) q[5];
cx q[5],q[1];
u1(2.31638910741203) q[1];
u3(-2.87950863608734,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.981179226647277,0.0,0.0) q[5];
cx q[5],q[1];
u3(3.05353037953235,2.16081773605711,-1.58664105039583) q[1];
u3(1.70229684529656,-0.359275553713835,5.56863164426036) q[5];
u3(1.13648759493369,-0.553493877311965,1.59728887385583) q[3];
u3(0.847932143366858,-1.51848078082512,-0.301574577666078) q[11];
cx q[11],q[3];
u1(0.0713811955589221) q[3];
u3(-0.802577450552867,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.49591936826338,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.20581134903636,3.27239827491092,-2.63016854183603) q[3];
u3(0.972107730813300,1.14596562753319,-2.06974878094550) q[11];
u3(1.79000200434649,-1.44348220382124,0.709715277242205) q[12];
u3(2.11330973516443,-3.99525354298501,0.277109618802682) q[10];
cx q[10],q[12];
u1(1.83486896675192) q[12];
u3(-2.40389459427356,0.0,0.0) q[10];
cx q[12],q[10];
u3(3.39708615694412,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.08463629206617,0.314590350096179,1.00276350351648) q[12];
u3(1.54727992560682,0.827385262701135,0.438539551742229) q[10];
u3(2.12218942525841,-0.982427757352764,2.21997955492451) q[14];
u3(2.76064883263742,0.0792409819815862,2.45938164862319) q[0];
cx q[0],q[14];
u1(0.135095580630262) q[14];
u3(-1.64961599211994,0.0,0.0) q[0];
cx q[14],q[0];
u3(2.13668720055335,0.0,0.0) q[0];
cx q[0],q[14];
u3(1.06819236443566,0.945445977314832,-2.49623724861942) q[14];
u3(0.929339894811829,1.95618529982319,0.123242124446279) q[0];
u3(1.48338102378437,0.698150334849822,-3.38712215863932) q[2];
u3(0.906828052316502,3.15516257312249,-2.86886568518910) q[4];
cx q[4],q[2];
u1(-0.0517546944627207) q[2];
u3(-1.31872346271411,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.22376232335860,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.34641448562400,-2.95716349975076,2.30482519831555) q[2];
u3(2.02089880725022,0.733164355644651,-2.79452974558598) q[4];
u3(0.655072423781820,-1.33572684542418,1.23619934823709) q[9];
u3(0.693867272674596,-2.93830840677706,1.39480592470662) q[8];
cx q[8],q[9];
u1(3.28351933368412) q[9];
u3(-2.01515199954083,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.48996494183634,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.90139911604806,-3.91958604122346,0.750340063471147) q[9];
u3(1.06511450290783,-2.43524513043029,2.10287824515099) q[8];
u3(2.33970086973162,1.72664700261565,-2.58337667882297) q[7];
u3(1.85679992595782,1.88683914183229,-1.87937831107348) q[13];
cx q[13],q[7];
u1(1.24330942438451) q[7];
u3(-3.40284980310114,0.0,0.0) q[13];
cx q[7],q[13];
u3(2.10770892683448,0.0,0.0) q[13];
cx q[13],q[7];
u3(0.425962041591350,-2.68663915376437,-0.396136344283990) q[7];
u3(1.97067212728460,0.238037539125362,-5.69659400508549) q[13];
u3(1.30810010211994,2.90843005750857,-2.32204325472663) q[11];
u3(1.63018367291141,0.657012041018110,-2.21948637232198) q[10];
cx q[10],q[11];
u1(0.849684227106801) q[11];
u3(-3.44918636215690,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.03610303698886,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.61577913760663,2.55818992159862,-2.09050248057183) q[11];
u3(0.920241980883877,-0.281671427564193,-1.01076391230553) q[10];
u3(2.49237194788286,2.33808974397887,-3.50592810918373) q[7];
u3(1.45354831834939,2.09438487954394,-1.79959425775953) q[0];
cx q[0],q[7];
u1(1.40516663794851) q[7];
u3(-0.0370456865829953,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.82420375201841,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.752293558425553,0.609091743731443,1.85466507048627) q[7];
u3(1.57421099961642,1.60016888924720,0.398549945157074) q[0];
u3(2.44192100325531,1.86701639722260,-4.39675745437891) q[14];
u3(0.482652228101032,-1.02359579366777,3.03794912736016) q[8];
cx q[8],q[14];
u1(1.66196637578319) q[14];
u3(-2.34277510101084,0.0,0.0) q[8];
cx q[14],q[8];
u3(3.67164828536508,0.0,0.0) q[8];
cx q[8],q[14];
u3(1.33473678371325,-1.74829019868359,3.33217529551314) q[14];
u3(0.795175427305146,4.48554681127494,0.839204572439838) q[8];
u3(1.04379473022619,1.78244998511504,-3.41759620679413) q[4];
u3(2.30506532447415,2.70454205781857,-3.38423148635221) q[5];
cx q[5],q[4];
u1(4.26835690479845) q[4];
u3(-2.94534988289814,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.555448780775669,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.43814046195970,2.82545366709753,-1.04212022096764) q[4];
u3(2.40238635271093,-0.814567688337696,-0.980842570816775) q[5];
u3(1.70639245629121,0.189269984395745,2.49542797465638) q[3];
u3(1.52447155063584,-0.285146663261837,-1.46382927305600) q[1];
cx q[1],q[3];
u1(3.44849975530094) q[3];
u3(-4.57971617501396,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.338639682149631,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.66648349796288,-3.73295068246217,1.54874888635660) q[3];
u3(1.51437396794431,0.474296826037784,-1.75793739237210) q[1];
u3(0.705059339568248,-0.895094215974160,-1.83201308522099) q[13];
u3(2.05869649320605,1.27113890732360,-4.52277650807401) q[2];
cx q[2],q[13];
u1(4.26279632472245) q[13];
u3(-3.15127496420799,0.0,0.0) q[2];
cx q[13],q[2];
u3(-0.409716449083805,0.0,0.0) q[2];
cx q[2],q[13];
u3(2.86489510293838,-3.78376072521163,0.523303903843707) q[13];
u3(1.07251953679175,-3.44679782690186,2.27343096121077) q[2];
u3(0.905352375804355,1.72193421621963,-2.02445998502671) q[12];
u3(0.463449355736168,0.599024539123360,-1.30670256681557) q[6];
cx q[6],q[12];
u1(1.59059139375453) q[12];
u3(-0.986450261749662,0.0,0.0) q[6];
cx q[12],q[6];
u3(2.70822678512832,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.43194221851544,1.57812284850382,0.976400389298054) q[12];
u3(1.19455982290378,3.96373972131217,1.85345985306865) q[6];
u3(2.93930848355788,2.63060609673439,-1.67818394223913) q[9];
u3(2.21652102173799,5.12311863445134,0.700487605378255) q[4];
cx q[4],q[9];
u1(2.59306885078014) q[9];
u3(-1.76980784310457,0.0,0.0) q[4];
cx q[9],q[4];
u3(-0.0181725992499255,0.0,0.0) q[4];
cx q[4],q[9];
u3(2.83133483384796,-0.469798543478384,-0.631167623826124) q[9];
u3(1.71607632392847,-2.68616430317185,2.72166176985194) q[4];
u3(1.04221091277670,1.20006940190011,-1.25646257545157) q[13];
u3(0.910393594178476,1.28682957485085,-4.70502072247365) q[0];
cx q[0],q[13];
u1(2.95406780682272) q[13];
u3(-1.64160758424894,0.0,0.0) q[0];
cx q[13],q[0];
u3(0.669003927466739,0.0,0.0) q[0];
cx q[0],q[13];
u3(0.640998354696723,-0.516150490045730,-0.434794442080119) q[13];
u3(1.62781386968088,2.79312470237408,2.36876178331773) q[0];
u3(2.05756575889167,-1.92579113596930,-0.449678123554436) q[5];
u3(1.83355215647490,-1.91650439673285,0.646321011848800) q[10];
cx q[10],q[5];
u1(-0.484480937301662) q[5];
u3(-2.05751958048787,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.74797520412583,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.21070219874782,-1.41095094511221,-0.0583177826651190) q[5];
u3(2.37379371842075,0.778594934687590,0.296061916414852) q[10];
u3(0.254505438999196,1.68992595248049,0.0876151522617241) q[12];
u3(1.48850487169398,1.50148101974726,-0.391307403270322) q[2];
cx q[2],q[12];
u1(0.573966600286635) q[12];
u3(-1.34100092705392,0.0,0.0) q[2];
cx q[12],q[2];
u3(-0.118869178389278,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.69528930473332,1.81600967184326,-0.521937275697206) q[12];
u3(1.27940502521101,-0.460936756773245,-1.35358871196588) q[2];
u3(2.78458272370356,-0.967762727835862,2.20424428510222) q[8];
u3(2.94989526118718,-3.00345054407688,0.0468257836752255) q[11];
cx q[11],q[8];
u1(2.51637888538954) q[8];
u3(-1.65105322665356,0.0,0.0) q[11];
cx q[8],q[11];
u3(3.40927243035758,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.63685141196662,-1.39115515725057,3.52473872327660) q[8];
u3(1.27449650455282,-3.14360924237736,-0.975238878866980) q[11];
u3(0.572493798108264,-2.83266026196080,2.14267017718749) q[14];
u3(0.656668554041030,-2.12293460502409,-0.251229967585972) q[7];
cx q[7],q[14];
u1(1.28185770125452) q[14];
u3(-0.0620065392313864,0.0,0.0) q[7];
cx q[14],q[7];
u3(2.20376431986970,0.0,0.0) q[7];
cx q[7],q[14];
u3(1.80485246975608,0.465281528135540,4.05081266660445) q[14];
u3(1.74786961332677,1.00174651153531,-2.46836077596479) q[7];
u3(2.41714378848955,2.29498375041410,0.555042648947964) q[1];
u3(1.57074467101884,0.568995487236587,-2.48425729698610) q[3];
cx q[3],q[1];
u1(1.34329771843807) q[1];
u3(-2.65252490147361,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.16067072312889,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.16793831055159,-3.96323112271888,2.09553477168525) q[1];
u3(1.65198765708769,1.98270257893887,-3.08463267499131) q[3];
u3(0.833053550694581,1.21659003887080,-1.48269979255248) q[10];
u3(0.0355993298320251,-2.81250036488903,0.426436396999691) q[2];
cx q[2],q[10];
u1(2.60770841103753) q[10];
u3(0.0881978305536002,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.39525219012467,0.0,0.0) q[2];
cx q[2],q[10];
u3(2.38054061041526,-0.482580456958369,-2.06802541393579) q[10];
u3(1.81187840596746,-3.06724376372049,-1.20720460327462) q[2];
u3(1.37609288012294,-0.973798607431722,0.380004924662203) q[14];
u3(1.28751602425081,-1.41465925140513,-1.27226130698614) q[5];
cx q[5],q[14];
u1(1.00383663577413) q[14];
u3(-0.141940470892503,0.0,0.0) q[5];
cx q[14],q[5];
u3(1.91329927787216,0.0,0.0) q[5];
cx q[5],q[14];
u3(2.12343275423024,3.17503169492832,-2.18414113997342) q[14];
u3(2.00197340979036,-3.35230041593883,1.81921375312097) q[5];
u3(1.35531399194736,-0.545635052246061,-1.20027725080645) q[4];
u3(1.11780386385597,-3.67176673091649,1.25994809883051) q[7];
cx q[7],q[4];
u1(3.04751757608583) q[4];
u3(-1.69708300372920,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.331700716364977,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.74330453306938,-2.57602310707580,-0.278484644972300) q[4];
u3(1.51977847723100,2.07785890992837,0.818534765162642) q[7];
u3(1.49836124301721,-0.823038469569821,-0.291656978850433) q[9];
u3(2.43531618407261,-2.82888400690988,-0.352056252551104) q[3];
cx q[3],q[9];
u1(0.359736541313239) q[9];
u3(-1.49818464367947,0.0,0.0) q[3];
cx q[9],q[3];
u3(-0.163189663954372,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.49402163573199,-1.28988284771473,0.977751807549732) q[9];
u3(1.38792284782223,-3.10197817106726,1.41333608817245) q[3];
u3(1.64905673568504,-0.969218719757964,2.41549501237496) q[8];
u3(1.11750688678759,-1.55129218716292,-2.12883395908344) q[1];
cx q[1],q[8];
u1(2.21738904557926) q[8];
u3(-3.21909010516631,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.474431545534428,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.36589389500143,-0.887796972539314,-0.775378374914604) q[8];
u3(2.31081230820535,-2.65557052629725,2.37664787949614) q[1];
u3(1.83474573517747,-0.345598482004455,1.23179904103028) q[6];
u3(1.91178420743596,-1.60171711559257,-1.96866762817877) q[11];
cx q[11],q[6];
u1(-0.428824686579083) q[6];
u3(-1.48138604019075,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.99796151753744,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.39270705583388,1.56117617540030,-1.65620046784832) q[6];
u3(2.06964724979982,0.411192432928621,1.88989619579746) q[11];
u3(1.89824125266706,-1.38928329797361,4.43879553822894) q[12];
u3(1.42715823014795,1.33084308853641,0.233779106297451) q[0];
cx q[0],q[12];
u1(-0.0251786083516325) q[12];
u3(0.907931120536137,0.0,0.0) q[0];
cx q[12],q[0];
u3(3.82182849626861,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.49410556875479,-3.43366993562881,2.07570550677153) q[12];
u3(1.04157495531250,-2.18405465336755,-1.26046119589410) q[0];
u3(0.498249458426636,3.27633444306597,-2.08101012570132) q[1];
u3(0.767870370291835,-3.22721263511501,1.24525567285677) q[14];
cx q[14],q[1];
u1(-0.376101564986660) q[1];
u3(-1.46735175461666,0.0,0.0) q[14];
cx q[1],q[14];
u3(0.714906672710632,0.0,0.0) q[14];
cx q[14],q[1];
u3(1.90207506799819,1.13679787122315,0.275476221052947) q[1];
u3(1.38980997091656,-2.20442243112592,-2.87077337606787) q[14];
u3(2.17618663210547,2.53769299606290,-1.33797833109594) q[6];
u3(2.36476580261244,0.321843339415521,-1.91799317856901) q[9];
cx q[9],q[6];
u1(1.59869802374079) q[6];
u3(-2.81367945131502,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.355976314229035,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.75673041844795,-4.07842873511054,0.965604998797623) q[6];
u3(1.69773739725755,5.64742007772058,-0.297235152269389) q[9];
u3(2.26407149903098,0.575300672457636,-3.03999867413677) q[13];
u3(1.43042530651396,-3.06074081178886,3.16205341576643) q[0];
cx q[0],q[13];
u1(3.30468302360306) q[13];
u3(-1.46397902239823,0.0,0.0) q[0];
cx q[13],q[0];
u3(2.78285611661565,0.0,0.0) q[0];
cx q[0],q[13];
u3(0.303019952162489,0.321504855548721,-1.01139713669281) q[13];
u3(0.564291740759678,-3.84059811201780,0.907407869240630) q[0];
u3(1.36039209197946,0.183653701163601,0.897329983613951) q[10];
u3(1.91442535089966,-0.962894425252379,-1.82690365395949) q[7];
cx q[7],q[10];
u1(0.282375383909354) q[10];
u3(-0.826020618280987,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.10910494985511,0.0,0.0) q[7];
cx q[7],q[10];
u3(2.03050066945890,1.43503452503771,-1.00973291985875) q[10];
u3(1.03054732740008,1.35229656942859,-1.14388001507141) q[7];
u3(1.11635020590935,-0.562654988261754,0.794213966547001) q[4];
u3(0.839695715716267,-2.20029249975177,0.340503985907028) q[11];
cx q[11],q[4];
u1(1.15369124869781) q[4];
u3(-3.33482068325270,0.0,0.0) q[11];
cx q[4],q[11];
u3(2.43513422979313,0.0,0.0) q[11];
cx q[11],q[4];
u3(2.28368274549431,2.20043830946636,-2.96827595176752) q[4];
u3(1.08798045052969,2.04823424073682,4.16751258385945) q[11];
u3(0.487575234090848,-0.612692448862342,0.840403580533690) q[8];
u3(0.557642076041875,-1.41797694377245,-1.09361879894234) q[3];
cx q[3],q[8];
u1(-0.443142560033433) q[8];
u3(-2.02551053352670,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.39655497541456,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.65338742275220,2.06024857548725,-1.09595025021802) q[8];
u3(2.71448155574111,-2.44862667069912,-1.60438044216882) q[3];
u3(2.86975004609519,1.73193732377016,-3.78456222866445) q[2];
u3(1.01277135209376,0.741392835601557,0.383202154016104) q[5];
cx q[5],q[2];
u1(3.38426096011878) q[2];
u3(-1.47073070710945,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.83683707247068,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.834211598122114,-1.76653088373645,-0.201952321239286) q[2];
u3(0.599352737230176,-0.197454159955335,-2.18933066157726) q[5];
u3(0.392421045055318,2.28013880565537,-2.86662277593669) q[13];
u3(0.918420991740656,-3.53484587110783,1.88270405264954) q[9];
cx q[9],q[13];
u1(0.0851271063255277) q[13];
u3(-1.26205725067537,0.0,0.0) q[9];
cx q[13],q[9];
u3(2.33808167469797,0.0,0.0) q[9];
cx q[9],q[13];
u3(2.10414575530768,1.48977832960946,-0.536200930403553) q[13];
u3(0.836749720125114,0.836316181705542,-0.159588300479965) q[9];
u3(0.442931958192711,0.515528457726670,-0.690005825076874) q[4];
u3(0.970847311980667,1.89210025429981,-4.22348610126803) q[0];
cx q[0],q[4];
u1(-0.573098865990134) q[4];
u3(1.45779231654215,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.79450500013428,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.88144825118201,0.997100595742002,-2.11352383808272) q[4];
u3(0.937943732913999,-1.29377715356442,-3.37913993713891) q[0];
u3(1.16925214737678,2.52688916373747,-0.272838004664757) q[2];
u3(1.02397800871117,0.612829292332651,-2.32800150257627) q[12];
cx q[12],q[2];
u1(0.173072688223434) q[2];
u3(-0.875580765297841,0.0,0.0) q[12];
cx q[2],q[12];
u3(2.69355149804653,0.0,0.0) q[12];
cx q[12],q[2];
u3(1.74680712720962,-0.192276603657432,0.791916806935562) q[2];
u3(0.686574718324287,-3.08560116042348,-1.11217129100327) q[12];
u3(2.44258566512502,3.21078063026358,-0.260084518052500) q[8];
u3(2.32051673921494,2.83168271364665,-2.11472920536704) q[3];
cx q[3],q[8];
u1(0.968596084366270) q[8];
u3(-0.699245580691270,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.00604539771190,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.570543921981182,3.62978531168016,-1.84070168561051) q[8];
u3(2.23544386629015,-3.46213837052780,0.287842277231430) q[3];
u3(2.55242135489567,2.99618119736788,-2.35014499721919) q[10];
u3(1.44181366537312,2.00526759825579,-1.80272716104259) q[1];
cx q[1],q[10];
u1(2.02121590864758) q[10];
u3(-2.77098414120683,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.47506684907404,0.0,0.0) q[1];
cx q[1],q[10];
u3(2.69221069279752,1.22735084266649,-4.33808510398425) q[10];
u3(1.11179698468941,4.79282962339637,-0.901717748303637) q[1];
u3(1.58837957944251,-0.457240003557932,1.57037459031506) q[6];
u3(2.02310648804115,-2.08612537409269,-2.76657668076927) q[7];
cx q[7],q[6];
u1(0.266992051509862) q[6];
u3(-1.47215871535421,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.40584712916229,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.77584264299336,0.833918540199526,-1.50165412604701) q[6];
u3(1.95002748865418,-2.43271312172303,-3.42914770891131) q[7];
u3(1.72665804994758,-1.93208662823254,-0.311223591943424) q[11];
u3(1.53794229182968,-3.85226095870667,-0.306503372454861) q[5];
cx q[5],q[11];
u1(2.11955335369251) q[11];
u3(0.0604402156660462,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.56892147096928,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.81807744282830,-4.28459756721122,1.83639953512210) q[11];
u3(2.57147080559383,-1.79875602557160,3.83040966613314) q[5];
u3(0.728140119220777,2.10560757326228,-1.72809350153816) q[6];
u3(1.00102713045891,2.27927388492808,-1.66265228093339) q[8];
cx q[8],q[6];
u1(0.296811793563372) q[6];
u3(-1.17924781446179,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.60002272063432,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.614484737610148,3.45970376822799,-1.96264661725751) q[6];
u3(1.63378823699270,-1.74759553800169,-0.353295788378492) q[8];
u3(1.16886540488346,-4.12974902197599,1.08886973716619) q[7];
u3(2.09125023793352,-5.25835580659689,0.0529607487105319) q[5];
cx q[5],q[7];
u1(0.896867232890609) q[7];
u3(-3.21887676077325,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.07417308337217,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.47698636402310,-4.61543695418878,0.998024574244869) q[7];
u3(1.08375093899366,-1.76577227803989,1.11579801192066) q[5];
u3(2.47950776609432,-0.385721345766985,-1.67644815263107) q[12];
u3(1.97221150416503,1.27799591711798,-4.61025605540494) q[9];
cx q[9],q[12];
u1(0.618981395161054) q[12];
u3(-1.52008820425584,0.0,0.0) q[9];
cx q[12],q[9];
u3(2.74320406514588,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.85355227464073,0.984241991021107,-3.96921568817894) q[12];
u3(0.841573945205331,-2.44041658477218,0.406404377961374) q[9];
u3(1.24239561673485,1.74886414048364,0.971001747439391) q[11];
u3(2.74768433527211,-0.456816771814528,-3.03521919911473) q[1];
cx q[1],q[11];
u1(4.36433323235918) q[11];
u3(-3.86021148163013,0.0,0.0) q[1];
cx q[11],q[1];
u3(-0.710242038753528,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.13236000560474,-1.30049813147260,3.31179757717260) q[11];
u3(1.36464428044927,-0.281561479068831,-2.51041693800550) q[1];
u3(0.579382317460228,0.392161874098468,-2.60628673677072) q[10];
u3(1.59116347373655,2.47882308935114,-3.17746421592243) q[4];
cx q[4],q[10];
u1(2.09425482358325) q[10];
u3(-1.88797220129879,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.305197978992426,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.62076739289181,1.02394666441746,2.68614541576143) q[10];
u3(1.66775865234493,3.29017967687872,1.14891272703443) q[4];
u3(0.796423225040962,-1.71206414067250,2.22424757400438) q[0];
u3(0.465781744454866,-2.04018458794094,-0.168302643998254) q[14];
cx q[14],q[0];
u1(1.54263479353350) q[0];
u3(-2.39121070884538,0.0,0.0) q[14];
cx q[0],q[14];
u3(3.40710437085978,0.0,0.0) q[14];
cx q[14],q[0];
u3(1.78281580456117,-0.0153113235166122,2.73304068310684) q[0];
u3(1.18127162278399,1.37434115183848,-2.06898530524323) q[14];
u3(1.15457275328500,2.78941474621821,-1.70650006834320) q[13];
u3(1.67026035978224,0.992292456095971,-2.77341402313718) q[2];
cx q[2],q[13];
u1(3.34046050148205) q[13];
u3(-0.847239970385513,0.0,0.0) q[2];
cx q[13],q[2];
u3(1.76469807539161,0.0,0.0) q[2];
cx q[2],q[13];
u3(2.24914596175567,2.84462570307644,-1.44018490150670) q[13];
u3(1.12994845564161,0.797964634709649,-3.37727447890926) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14];
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