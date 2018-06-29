OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0.926612486661964,0.334535365547991,2.18822995339968) q[13];
u3(0.790243615601293,-0.439067136317212,-2.77387067777406) q[4];
cx q[4],q[13];
u1(1.58537581113023) q[13];
u3(-2.20911722814672,0.0,0.0) q[4];
cx q[13],q[4];
u3(0.326919843440078,0.0,0.0) q[4];
cx q[4],q[13];
u3(0.684716104260289,1.18361451376944,0.336235167977638) q[13];
u3(2.05435412828562,-1.44335655963375,-4.30971402927632) q[4];
u3(1.47374580231717,1.11947795016417,-0.0992008602610986) q[15];
u3(1.92645314847875,0.405553387002306,-4.35332167068890) q[0];
cx q[0],q[15];
u1(-0.0843924417624573) q[15];
u3(0.816795193737886,0.0,0.0) q[0];
cx q[15],q[0];
u3(3.73344676939604,0.0,0.0) q[0];
cx q[0],q[15];
u3(2.05105273337391,-1.44070903387294,3.30812425683184) q[15];
u3(2.79494567263955,1.42661576113961,0.180368119592334) q[0];
u3(1.25765438284152,1.14022871553658,-3.53289684455499) q[6];
u3(0.516741522046894,1.74790715078884,-2.22782069639276) q[11];
cx q[11],q[6];
u1(1.43941591228781) q[6];
u3(-0.0112724334738727,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.26021860119654,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.296445704199084,-0.477716766417625,-2.85167518255505) q[6];
u3(1.54180183589136,-0.727456711418955,3.56727801314341) q[11];
u3(0.513242716592139,1.84629250309736,-2.52958608448031) q[2];
u3(1.55193411182018,3.83408952219594,-2.41809891886493) q[7];
cx q[7],q[2];
u1(0.459852438744524) q[2];
u3(-1.40880903336820,0.0,0.0) q[7];
cx q[2],q[7];
u3(3.03933280290410,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.472896725809127,0.797081043905856,-0.945953074147258) q[2];
u3(1.32446316277423,1.39942828798606,1.65030682020948) q[7];
u3(2.29956370211881,1.44993581505559,-0.937796385818231) q[3];
u3(1.45973196060168,-4.48363069325872,1.59401928837569) q[9];
cx q[9],q[3];
u1(0.313394563877017) q[3];
u3(-1.46585541745296,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.53969648263638,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.60074959195304,-0.889083441499602,-0.778535657483170) q[3];
u3(1.56500817994673,1.22917556786356,-4.88265668856868) q[9];
u3(0.435614101812310,-3.47116467133345,2.69253512728449) q[1];
u3(1.18879618679561,-0.0167088440219814,-2.35402590014308) q[14];
cx q[14],q[1];
u1(2.17284591318098) q[1];
u3(-2.85749150389280,0.0,0.0) q[14];
cx q[1],q[14];
u3(1.39199052276588,0.0,0.0) q[14];
cx q[14],q[1];
u3(2.46105359390513,-1.09223345483974,3.61064227325473) q[1];
u3(1.11651219917316,-2.27600202856713,-0.106644252479767) q[14];
u3(1.58293029785017,2.09991058577181,-2.69487975698424) q[12];
u3(1.00664903498077,-2.28196229284744,2.53104347272944) q[8];
cx q[8],q[12];
u1(3.14570686983616) q[12];
u3(-1.42297897508283,0.0,0.0) q[8];
cx q[12],q[8];
u3(2.63456056352596,0.0,0.0) q[8];
cx q[8],q[12];
u3(2.01153985204762,0.588484328153723,-0.698179650470414) q[12];
u3(1.64692315020026,0.0883414868722806,2.08239902664084) q[8];
u3(0.847080156489057,-0.637804662095740,-2.23624649307059) q[10];
u3(2.20894095339901,1.83313802778196,-4.41351126698182) q[5];
cx q[5],q[10];
u1(1.47177008673499) q[10];
u3(-3.70113031904954,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.03317278259618,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.863180735941614,-2.60343591464365,1.64998203771628) q[10];
u3(0.224337650102777,0.0408400457237534,5.76839224595487) q[5];
u3(1.74116361941023,1.92165174027992,-0.596712628458340) q[7];
u3(2.38080696788652,5.25368732103930,0.685286798889622) q[2];
cx q[2],q[7];
u1(1.11530450573778) q[7];
u3(-3.33644500173379,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.07531184581825,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.721030260760104,2.66210287051748,-1.34433065720516) q[7];
u3(2.98184851507564,-3.38235494978346,-0.746152639280928) q[2];
u3(0.410802007369724,-2.37946939369771,3.46553509154976) q[14];
u3(0.340914089497046,-3.87071292416668,2.27300586878866) q[1];
cx q[1],q[14];
u1(1.31911925963125) q[14];
u3(0.203416620823137,0.0,0.0) q[1];
cx q[14],q[1];
u3(2.58593446871005,0.0,0.0) q[1];
cx q[1],q[14];
u3(0.836428818348620,-0.671971830034125,0.641631391536425) q[14];
u3(1.19991782438439,-1.18666064512141,2.84933529719089) q[1];
u3(1.57370136686121,1.36770603261326,-3.78432056837763) q[4];
u3(1.26743519666428,2.92779498176749,-2.60146143909628) q[0];
cx q[0],q[4];
u1(3.28410810313236) q[4];
u3(-1.56078394534726,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.734287605110107,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.19789811655093,0.918809938866682,-2.38901135022255) q[4];
u3(2.22630287310153,4.25381830848836,0.0131311611717004) q[0];
u3(0.702079960237297,0.480592679065253,-0.920094230938103) q[13];
u3(0.923127664086138,-3.30319613794939,1.29533435326417) q[6];
cx q[6],q[13];
u1(2.40083736250768) q[13];
u3(-1.46998732098316,0.0,0.0) q[6];
cx q[13],q[6];
u3(3.37234227094943,0.0,0.0) q[6];
cx q[6],q[13];
u3(1.77798108559155,1.43142218397449,2.52621875402862) q[13];
u3(1.35741696934532,0.0195789247485239,-2.64323621219100) q[6];
u3(1.28545040882205,-0.652672672344946,1.35362518675849) q[3];
u3(0.931185561646461,-1.63277222698955,-1.31318560243450) q[11];
cx q[11],q[3];
u1(2.06705768780428) q[3];
u3(0.147156211440483,0.0,0.0) q[11];
cx q[3],q[11];
u3(1.15905194165921,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.59490097789868,1.09885418228575,-3.47567640195203) q[3];
u3(1.63405091279929,-3.62434314129094,2.11402772361462) q[11];
u3(1.76177239618319,-0.901939204366000,2.34538729607642) q[9];
u3(1.07576869335000,-1.33291597524784,-1.48412832934194) q[12];
cx q[12],q[9];
u1(1.50886078031871) q[9];
u3(-3.06460013504857,0.0,0.0) q[12];
cx q[9],q[12];
u3(0.785099634562065,0.0,0.0) q[12];
cx q[12],q[9];
u3(2.30874155263359,-2.90356557617848,1.40704748380377) q[9];
u3(2.75046481582968,-2.03639217284113,-4.19225584417618) q[12];
u3(0.719957912200558,-0.513438737194430,-2.01818309277149) q[10];
u3(0.987997436792502,1.69892717367435,-4.26526953925844) q[5];
cx q[5],q[10];
u1(1.28351770650182) q[10];
u3(-3.65560951874237,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.10898656066071,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.59739156475968,-0.0188290854753179,-3.15168850285208) q[10];
u3(1.58736323223800,2.77220298270251,-1.64354482976360) q[5];
u3(1.74002736126560,3.87245552904838,-1.11113126929035) q[15];
u3(2.02124224114725,3.53455212003298,0.478094698867159) q[8];
cx q[8],q[15];
u1(4.15652167094208) q[15];
u3(-3.25839598533970,0.0,0.0) q[8];
cx q[15],q[8];
u3(-0.666753720780228,0.0,0.0) q[8];
cx q[8],q[15];
u3(1.23938565028430,0.0428080991698686,1.19713051208350) q[15];
u3(1.02649154552393,-1.99998809523666,-3.65640366093708) q[8];
u3(1.29833851042391,3.44237174787001,-2.12543318359421) q[13];
u3(2.17060388794672,2.13997283270019,-1.09552343123775) q[14];
cx q[14],q[13];
u1(-1.42200717623848) q[13];
u3(0.274343976912759,0.0,0.0) q[14];
cx q[13],q[14];
u3(3.26949278123920,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.62357509774858,2.83990725877541,-2.95284128551051) q[13];
u3(1.17517569176201,-4.91394625405886,1.30860536488676) q[14];
u3(2.60868954685169,2.93677763942229,-0.137479647801377) q[15];
u3(1.87763344603978,1.15461855751757,-4.12641190402321) q[11];
cx q[11],q[15];
u1(0.164285962048153) q[15];
u3(-1.59953410982165,0.0,0.0) q[11];
cx q[15],q[11];
u3(2.11077324435709,0.0,0.0) q[11];
cx q[11],q[15];
u3(2.78707672093369,0.648029926157679,0.880357640333773) q[15];
u3(1.44691625881801,4.64501173059631,-0.266664488182029) q[11];
u3(1.29252893010070,1.65790909058801,0.977100784913419) q[1];
u3(0.834950436690865,-0.0584037520694267,-3.13301334305564) q[5];
cx q[5],q[1];
u1(2.34788591530748) q[1];
u3(-2.98452411366002,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.53177086074444,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.06321962376949,-1.66420643180622,0.0772698972838840) q[1];
u3(1.43710953141631,-3.93914245889114,-0.750223031449087) q[5];
u3(1.48450152126594,-2.30281840116940,-0.348361022495402) q[9];
u3(1.34231872208419,-4.35144377025089,-1.64427436542594) q[12];
cx q[12],q[9];
u1(1.37516071460252) q[9];
u3(-3.37044943452111,0.0,0.0) q[12];
cx q[9],q[12];
u3(2.37156556128165,0.0,0.0) q[12];
cx q[12],q[9];
u3(1.11477683505507,2.55286604309560,-3.19073759614167) q[9];
u3(1.01182073928967,1.01772777188266,1.45924471253315) q[12];
u3(1.82073966097737,0.469164719035575,1.09721413067347) q[6];
u3(1.44213181900268,-0.618676423299432,-2.16197041521389) q[10];
cx q[10],q[6];
u1(1.48521983114405) q[6];
u3(-0.629157758653331,0.0,0.0) q[10];
cx q[6],q[10];
u3(-0.465300204495927,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.64636281006604,-3.53101877220634,0.113624725478770) q[6];
u3(2.39940689379025,3.83862802610010,-0.600105519087495) q[10];
u3(0.539862919151970,-2.39579347101177,1.81996189109310) q[8];
u3(1.02096683236807,-3.91977806361495,1.98213745373712) q[0];
cx q[0],q[8];
u1(0.953966721365830) q[8];
u3(-3.63559419085405,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.77729532740569,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.79510933475709,1.37576155055209,0.970778723925616) q[8];
u3(0.573129727388788,1.26084193071603,-0.151127716253841) q[0];
u3(1.37248199800735,0.815737043773334,1.53065034138529) q[4];
u3(1.75125041418529,-1.10313720267351,-0.394935035585971) q[7];
cx q[7],q[4];
u1(3.15266652526979) q[4];
u3(-0.744011871313877,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.01048560904574,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.40207566571970,3.70048269455448,-2.17567564577368) q[4];
u3(2.23252912145351,2.13001442816654,2.53626913931597) q[7];
u3(2.21143204046302,3.73855888832742,-1.58468625722196) q[2];
u3(2.10548581408415,1.56890676021001,-2.67088677712459) q[3];
cx q[3],q[2];
u1(0.366879669119202) q[2];
u3(-1.28785202717118,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.84675129513782,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.62679843792713,0.408482051741396,2.84260999362147) q[2];
u3(1.85911200868554,-2.85680449163391,0.154243478967809) q[3];
u3(1.01035345280299,0.165769987382684,2.16475975707169) q[2];
u3(0.816263464212845,-0.799577782590861,-2.59754489124988) q[11];
cx q[11],q[2];
u1(2.19861174953091) q[2];
u3(-0.234040747482894,0.0,0.0) q[11];
cx q[2],q[11];
u3(1.51918541659615,0.0,0.0) q[11];
cx q[11],q[2];
u3(0.724945037832755,-3.59933838825048,2.55183195457948) q[2];
u3(2.97117525399314,-0.370428480760825,2.50255716780721) q[11];
u3(1.79809993423367,-0.179406804652588,-0.868996516856425) q[3];
u3(1.15606351621790,-4.34249751651053,0.988617317746629) q[1];
cx q[1],q[3];
u1(2.26156957640115) q[3];
u3(-1.89552337213579,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.40173526379898,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.05859245485654,-0.524004184140046,-3.44606664740587) q[3];
u3(1.09958357543022,4.93268330569393,-0.363268024720373) q[1];
u3(2.80090840100965,0.842704971506356,-2.28572654970770) q[4];
u3(2.65928408165709,5.88364062324891,-0.293984868937056) q[8];
cx q[8],q[4];
u1(2.88259148748486) q[4];
u3(-1.67820458099804,0.0,0.0) q[8];
cx q[4],q[8];
u3(3.20703508784886,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.14658854635964,-2.02625649086564,2.01730922243307) q[4];
u3(1.03437060211118,0.661331226996970,-4.21546853337167) q[8];
u3(0.830798864080012,-0.784368973369230,1.14192020544190) q[14];
u3(0.785385955128842,-2.70066204597485,0.409547620179870) q[6];
cx q[6],q[14];
u1(0.0542952999740944) q[14];
u3(-1.02638135398860,0.0,0.0) q[6];
cx q[14],q[6];
u3(1.61502958754216,0.0,0.0) q[6];
cx q[6],q[14];
u3(0.510693009232525,3.56936937341456,-2.21342983824698) q[14];
u3(2.16728598562797,2.36985009254488,-0.569923579041209) q[6];
u3(1.70599417392261,0.616576981409986,-1.61270345064447) q[10];
u3(0.567655530040565,0.406010730495276,-2.99804132534289) q[13];
cx q[13],q[10];
u1(-0.453419768725846) q[10];
u3(0.838637619058213,0.0,0.0) q[13];
cx q[10],q[13];
u3(4.08602513194895,0.0,0.0) q[13];
cx q[13],q[10];
u3(2.48225551427039,-1.01689219089910,2.19627239281382) q[10];
u3(1.32910279468909,1.53766241926350,-2.17896934808221) q[13];
u3(2.71571478978616,1.12178068653255,-0.435276343859132) q[5];
u3(1.76783878567841,-1.04417806902591,-2.54997259443778) q[15];
cx q[15],q[5];
u1(1.15784440236244) q[5];
u3(-3.65469150600077,0.0,0.0) q[15];
cx q[5],q[15];
u3(1.58348196143482,0.0,0.0) q[15];
cx q[15],q[5];
u3(1.11377856382645,-2.39987964588206,-1.39377465723697) q[5];
u3(2.58708675464058,2.85393264395235,1.70025850541179) q[15];
u3(1.30825395657615,-1.19815243366965,1.06680067156640) q[9];
u3(1.11513914079766,-2.95558976960256,-0.555955375840677) q[12];
cx q[12],q[9];
u1(2.34255210101516) q[9];
u3(-1.65969065823003,0.0,0.0) q[12];
cx q[9],q[12];
u3(3.38273444254001,0.0,0.0) q[12];
cx q[12],q[9];
u3(1.95580895897167,-1.44921463774424,-0.304776318830294) q[9];
u3(1.43575830601153,2.18181969239475,0.947274539823610) q[12];
u3(2.54273781472425,1.53619188671009,1.28713661230262) q[7];
u3(0.613838653307615,-4.59270044082719,0.231845935278624) q[0];
cx q[0],q[7];
u1(3.51848042657554) q[7];
u3(-3.22825726616868,0.0,0.0) q[0];
cx q[7],q[0];
u3(-0.831757663828052,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.16081120501014,-1.44283469991507,-1.66931871300157) q[7];
u3(1.27241102856116,-0.151713147823791,3.36130713460828) q[0];
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
