OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(0.383824731999968,1.85848743097718,-0.923193685096286) q[5];
u3(1.27545445411932,0.246881452234029,-2.96123820645508) q[3];
cx q[3],q[5];
u1(1.83221378444308) q[5];
u3(-2.61284073096125,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.936965922418506,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.776319427169728,-1.23293599792503,-1.96277779191324) q[5];
u3(1.69087190820660,3.52009985862103,-1.43568079237321) q[3];
u3(1.58321242459898,3.32988188766697,-2.02986358965018) q[4];
u3(0.995776972007453,1.06169358931712,-0.721749639626722) q[1];
cx q[1],q[4];
u1(-0.0547909345917583) q[4];
u3(-2.10698855382812,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.849049917058339,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.29007823940406,0.635365284127704,0.697745782755670) q[4];
u3(0.707433936035816,0.275196462411947,-0.393363931609194) q[1];
u3(0.689151051764536,2.36293159292421,-3.01142278702107) q[2];
u3(0.291429481327411,-0.743243806099744,-0.787341826001340) q[0];
cx q[0],q[2];
u1(0.632520938286125) q[2];
u3(-1.67104772916662,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.219262654076592,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.43276086546646,2.57085194785027,-0.798374183787509) q[2];
u3(1.99585997443318,1.52960887786100,-1.54053009209085) q[0];
u3(0.842265281615917,0.820179636887532,-2.12245965465026) q[2];
u3(0.424226261453600,-1.70520172574211,0.130684886051523) q[4];
cx q[4],q[2];
u1(1.42765528770378) q[2];
u3(-0.930512288388228,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.0238920143991608,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.560173631791268,1.66423144061010,-0.926030806363942) q[2];
u3(0.929989694873875,5.42356054992320,0.809345994810085) q[4];
u3(1.87751429095357,1.66681117572721,-3.01667305884975) q[0];
u3(0.907050510936757,-2.86276693822534,3.31360477695410) q[3];
cx q[3],q[0];
u1(0.281669367065745) q[0];
u3(-1.79461693294442,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.982771879462051,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.42357319744002,-1.85448867355199,4.33478081259563) q[0];
u3(1.98565254310010,3.11126878953311,1.10411489735171) q[3];
u3(1.27351364276001,1.53288620652976,-2.59673014344655) q[1];
u3(1.36473635609074,-3.85484513972344,2.15596097102823) q[5];
cx q[5],q[1];
u1(0.973303875741875) q[1];
u3(-3.32197976773712,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.07641053456435,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.19466639748686,0.793993774723354,-1.34821492727101) q[1];
u3(1.45359368976888,1.91716317036500,2.69009080457001) q[5];
u3(2.44105913779209,-1.69090070709913,0.612797286528868) q[0];
u3(2.31557944430035,-3.14585478780228,-1.17105637287377) q[5];
cx q[5],q[0];
u1(1.66829996241053) q[0];
u3(0.707353361262258,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.27173021180939,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.90957093102910,-2.54007222877059,-1.16446129751931) q[0];
u3(2.19394691701234,-0.436941348448315,-3.23547581441985) q[5];
u3(2.00609588628877,3.30955346172354,-2.07667362078287) q[1];
u3(1.58789214997175,1.85462989902791,-2.09581686719514) q[3];
cx q[3],q[1];
u1(0.703955171279338) q[1];
u3(-1.35508953122091,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.76334966888768,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.23969431723788,-0.768810486744889,0.640240238476777) q[1];
u3(1.96696980813893,-0.218448969764328,-5.05225235283075) q[3];
u3(1.36839219722757,-0.702153471339964,0.197573383343396) q[4];
u3(0.498496618598099,-2.67387350249044,-0.567978351695543) q[2];
cx q[2],q[4];
u1(3.41608242486300) q[4];
u3(-1.55766123530634,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.55285514081029,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.26527613377739,3.94645298804065,-1.63725744684847) q[4];
u3(1.70609099525561,4.82289391610341,-1.03225010999058) q[2];
u3(3.00437115653739,-2.39544538169490,3.82857208600381) q[5];
u3(1.10177954668123,-0.661765780913166,2.89413521442565) q[2];
cx q[2],q[5];
u1(3.79102436800288) q[5];
u3(-4.34962636627458,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.359435533275606,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.11835835853448,-1.35755743885092,0.586551503747049) q[5];
u3(2.17070521734752,0.873333956226582,-0.476576977271353) q[2];
u3(2.65123191079145,-0.958960593082272,0.381877340749941) q[4];
u3(1.98942147672166,-2.90304073925058,-1.13528703579528) q[3];
cx q[3],q[4];
u1(1.51568737487150) q[4];
u3(-3.52445626331403,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.15480542060276,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.20913756503547,-2.70376656335883,3.34739477274897) q[4];
u3(1.76948603619978,0.492348079487182,-5.63273160909096) q[3];
u3(2.57911014960012,-3.65790490968423,1.82911270757939) q[1];
u3(0.858428213678491,1.67601666991474,-0.195132540756579) q[0];
cx q[0],q[1];
u1(3.47535096416748) q[1];
u3(-1.46578190757312,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.86073296197251,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.87593321678431,-2.47600336505328,1.18952322158539) q[1];
u3(2.50493220341757,4.67318682854506,-0.0213599910091977) q[0];
u3(1.67087790063469,-0.589596682461751,-0.265917853720976) q[0];
u3(0.708688126693857,-2.32281413783732,-1.18204408805136) q[2];
cx q[2],q[0];
u1(4.04568323315006) q[0];
u3(-1.35166928791485,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.58282797411651,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.07360010434733,3.67868056910102,-1.91919779037552) q[0];
u3(1.50264709269433,-1.11695704079643,2.65449071779074) q[2];
u3(1.70121239852643,-0.492160787074758,0.354923035566312) q[4];
u3(1.55476290633374,-3.15805502619932,-0.597715405841909) q[3];
cx q[3],q[4];
u1(2.19049839913023) q[4];
u3(0.265074909161245,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.39322706454978,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.97180498235888,0.728592191294519,-2.77920395460681) q[4];
u3(2.83362749851578,4.68391380631129,0.0668805231819278) q[3];
u3(2.68624831480438,3.44406467045351,-0.686217834387517) q[1];
u3(2.28654499447329,2.65937895807429,-2.25373492847583) q[5];
cx q[5],q[1];
u1(2.09071219473348) q[1];
u3(-2.56805094532546,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.22564530035902,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.391026004555584,1.53583830882186,-4.40432788075507) q[1];
u3(2.41776649942909,-2.37977272261287,3.29306204354595) q[5];
u3(1.65561577512260,-0.581434694814940,2.38481679357750) q[4];
u3(0.796620411054835,-0.718031909285691,-1.57092131476580) q[2];
cx q[2],q[4];
u1(1.97561791389809) q[4];
u3(-3.13348801556303,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.79244533428263,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.17786457718549,-2.04165846406785,3.73357184952280) q[4];
u3(0.837185395097511,-1.45568513577936,1.34367022615734) q[2];
u3(2.39134596144722,2.67951484707837,-0.474892967054501) q[5];
u3(2.52916672180310,1.51374630793657,-2.45900198402994) q[0];
cx q[0],q[5];
u1(2.89606728694928) q[5];
u3(-1.75620600584356,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.931298405476309,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.21432407971383,0.643684278423894,-0.747059616308071) q[5];
u3(0.691524433761461,-2.80055694971532,-2.03759450784596) q[0];
u3(2.23547736335405,0.929642827555718,-3.45264217691445) q[3];
u3(2.95886797830257,1.56027584812501,-3.38169953367211) q[1];
cx q[1],q[3];
u1(2.30384350721837) q[3];
u3(0.621823993155638,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.66216843586695,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.16586009047662,2.64519496334600,0.187474873012217) q[3];
u3(2.16890211988257,3.86623211277568,-0.343855051862270) q[1];
u3(0.932134848377261,1.94458635591571,-3.91759892003361) q[4];
u3(1.82972285895003,2.36308795963407,-3.23584873976300) q[2];
cx q[2],q[4];
u1(1.77063062678159) q[4];
u3(-2.13090221668657,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.192772832370326,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.30315076138037,-1.57015071666872,3.21951159815291) q[4];
u3(0.841381466968607,-2.47782894594107,3.66075055689429) q[2];
u3(1.57439065494551,-2.15765248117598,0.343375473051380) q[5];
u3(0.957536559495765,-4.28937529719360,0.704356495184570) q[1];
cx q[1],q[5];
u1(-0.298259373400143) q[5];
u3(0.999958739325005,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.55067185826812,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.520287273155836,2.71964510617080,-1.32345297817976) q[5];
u3(1.35476985107209,-3.16740116795862,2.50705089804593) q[1];
u3(1.29959206569065,1.09229825724017,1.08214584575314) q[0];
u3(0.459117359024423,-0.442251601263720,-2.89563279038234) q[3];
cx q[3],q[0];
u1(1.42575975161119) q[0];
u3(-3.51436524849920,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.12353795183964,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.32271582013781,-0.769894793216452,4.60800308039340) q[0];
u3(0.914222546240896,-1.12044209668566,1.33327620862347) q[3];
u3(2.24546682590860,3.23555589963794,-0.462331375948911) q[5];
u3(2.58889505386854,3.27869822204110,-2.04125395478807) q[2];
cx q[2],q[5];
u1(4.53553328288318) q[5];
u3(-3.69876640009773,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.600222218481427,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.78079821864620,1.25057826206681,-2.56136031788229) q[5];
u3(0.500537517805709,0.0862002733435576,-2.55054079703476) q[2];
u3(1.91649891016214,-0.780348001717189,3.29099169454321) q[1];
u3(2.88994418120637,-3.29016362514161,-1.23579236075600) q[0];
cx q[0],q[1];
u1(3.00546336948225) q[1];
u3(-2.24154286336288,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.12758887033751,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.12611318396533,-0.633393617568103,0.403245676781527) q[1];
u3(1.83862467251199,-4.42538195117690,1.67699981512629) q[0];
u3(1.12090198627273,2.04286221355318,-2.01215797142866) q[4];
u3(0.406668112938987,-2.79893460062846,1.71915497821382) q[3];
cx q[3],q[4];
u1(2.39729501508939) q[4];
u3(-1.73908979711976,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.329681007466804,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.96281054282752,-2.24272714448646,4.03713816886244) q[4];
u3(1.89633479605589,-3.62393788885563,2.24971135673709) q[3];
u3(1.38331279169285,0.907800490297888,-3.83506226173821) q[1];
u3(1.05169021636780,2.26327457964132,-2.21130221604684) q[3];
cx q[3],q[1];
u1(0.829425994894040) q[1];
u3(-1.33327586195042,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.00796715685770,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.17777877693810,1.39578232314354,0.764969762818059) q[1];
u3(1.46986933516214,5.13740934198150,0.243547027410567) q[3];
u3(2.05694635314485,3.03170360317371,-1.74802853285168) q[2];
u3(1.69790352979649,1.54036708169846,-0.448776831595734) q[0];
cx q[0],q[2];
u1(1.53219754633197) q[2];
u3(-3.43097112672982,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.36971789951198,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.25096063395480,0.208694182388882,2.16634292519379) q[2];
u3(1.13829051411383,-0.443142329771741,-5.53537498837201) q[0];
u3(0.138853699383896,-2.11206766383173,3.41696041780864) q[5];
u3(1.26079576489524,0.666163790598730,-2.37855964451707) q[4];
cx q[4],q[5];
u1(1.02955512212035) q[5];
u3(-1.14396397155686,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.860366772282345,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.67105140883288,3.58056566975990,-2.56973199737434) q[5];
u3(1.40380516361705,1.33538945918720,-3.79078749815364) q[4];
u3(1.19598719896345,-0.854130945809976,2.36789635065368) q[3];
u3(1.30011915812933,-2.05713177625036,-2.06512396413358) q[2];
cx q[2],q[3];
u1(2.17919718230689) q[3];
u3(-1.59327534159529,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.441130009753586,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.71315744059528,-2.53922484486155,-0.374617848760340) q[3];
u3(1.11350078988787,2.60341162607315,0.963714627636746) q[2];
u3(2.65931585715257,-1.48301126291087,4.48760813373090) q[0];
u3(1.58459054130743,0.947208710075791,0.862155223405295) q[4];
cx q[4],q[0];
u1(1.64939913235892) q[0];
u3(-2.47651981465031,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.108716538510806,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.66571703921006,-0.622000904602349,-2.19522799209114) q[0];
u3(1.43992529289098,-0.111746424249906,4.26440744790751) q[4];
u3(0.844450282308965,2.79906204908281,-0.957079300054263) q[1];
u3(1.72992541649532,0.249044665978621,-3.81309595939752) q[5];
cx q[5],q[1];
u1(-0.828232866743178) q[1];
u3(0.117730017091242,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.80580799755361,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.863277145224468,0.0552723640240638,0.986512278651087) q[1];
u3(1.49542381658562,-0.339244007166644,1.66650697667751) q[5];
u3(1.38312067588574,1.41663585099500,-2.86845478801068) q[0];
u3(0.750181470382166,-2.35453503185355,1.88401905248987) q[3];
cx q[3],q[0];
u1(1.04419503150205) q[0];
u3(-2.69273228430048,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.283649612985428,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.08845913614932,0.721222409648424,-3.68959675791932) q[0];
u3(1.90942684894728,-2.53269351276980,-3.36386524390562) q[3];
u3(1.56021370967791,-1.61761646482404,-0.0489351884483719) q[1];
u3(2.27265804318974,-3.12971317663976,-1.07676294208311) q[2];
cx q[2],q[1];
u1(3.32660901731321) q[1];
u3(-1.21309529307718,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.25648463930596,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.04910723248371,-1.33120967394535,2.93026804412672) q[1];
u3(2.95036245354344,1.57057901626372,3.62875644494706) q[2];
u3(1.76283820876013,-1.05482725056046,0.209377310388000) q[5];
u3(1.60456485485294,-2.42572237469694,0.813085239042513) q[4];
cx q[4],q[5];
u1(0.970844605541341) q[5];
u3(-0.437998317254760,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.38756577455009,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.33532293837025,0.834130137004826,-0.786308322472743) q[5];
u3(2.03430354408551,-0.869394421361398,2.11641463253364) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
