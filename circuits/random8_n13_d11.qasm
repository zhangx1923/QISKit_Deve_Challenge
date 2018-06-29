OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.70922482607690,2.99461700661026,-1.51060986322316) q[3];
u3(1.88266621828592,2.12343158891731,-0.0119231031407512) q[7];
cx q[7],q[3];
u1(2.62772623125352) q[3];
u3(-1.69967456573768,0.0,0.0) q[7];
cx q[3],q[7];
u3(-0.0127101855634326,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.79243610339652,1.32764132632192,1.68103492083134) q[3];
u3(1.50189974936126,-2.35282641738683,-3.63701946837202) q[7];
u3(1.87665076824775,-1.09670970097156,1.22663156892767) q[2];
u3(2.35830266670440,-1.41886742384876,-1.36951753904848) q[11];
cx q[11],q[2];
u1(0.585905113194614) q[2];
u3(-1.35583076900405,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.93557416797213,0.0,0.0) q[11];
cx q[11],q[2];
u3(2.93335236696712,0.736155226623074,-0.336493460425627) q[2];
u3(0.831996420215777,-2.64245202009313,1.01583758242669) q[11];
u3(2.10678979272341,-2.56245104431925,0.304697258031914) q[5];
u3(2.05515908427047,-3.22972354305926,0.670100163199854) q[0];
cx q[0],q[5];
u1(1.93438027977890) q[5];
u3(-3.08858319719122,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.58664947195745,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.30513555141915,0.834207839782218,-0.671134111068650) q[5];
u3(2.12543702913879,-3.16533572089160,1.98131319853831) q[0];
u3(0.587836297582155,-2.11731881282766,2.32581641301361) q[9];
u3(1.07232635746792,0.927520817651997,-2.04916013992910) q[10];
cx q[10],q[9];
u1(-0.172891757611749) q[9];
u3(-2.06657743900324,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.57107499817707,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.681095025781294,1.17468427817714,-4.34282847596179) q[9];
u3(2.76734150546949,-1.18510863093916,-0.229704379359660) q[10];
u3(2.78434062307563,-3.11371430742291,0.163121396420386) q[12];
u3(2.04619846240681,-1.40715532668282,0.933099863243745) q[6];
cx q[6],q[12];
u1(0.994268191648513) q[12];
u3(-0.367180123984385,0.0,0.0) q[6];
cx q[12],q[6];
u3(1.82415788284815,0.0,0.0) q[6];
cx q[6],q[12];
u3(2.40504913965663,1.34063589691539,-3.52446975067222) q[12];
u3(0.953059794118576,2.88770765842949,-0.401171458789098) q[6];
u3(0.442926165548928,-0.153782842719748,1.00149600414746) q[4];
u3(1.15764095228834,-0.426176513175858,-1.67569713170834) q[8];
cx q[8],q[4];
u1(2.54797610223206) q[4];
u3(0.208772115826314,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.24738643947123,0.0,0.0) q[8];
cx q[8],q[4];
u3(0.790908349196343,4.03471707332438,-0.167473858961255) q[4];
u3(0.437810319439774,-2.47328423750986,-1.96299196734976) q[8];
u3(0.398170342363440,-2.58151514579579,3.56412895892695) q[10];
u3(0.309602110688995,-3.83781615925404,2.25862668488734) q[1];
cx q[1],q[10];
u1(3.43434828280410) q[10];
u3(-1.74980650987084,0.0,0.0) q[1];
cx q[10],q[1];
u3(2.45362555770854,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.66959078343399,0.925535336471272,1.85563273278954) q[10];
u3(1.47299185735557,-3.39474537327879,-1.43186851764815) q[1];
u3(1.67648807599338,-0.666261983854102,3.51879276007649) q[2];
u3(1.00459050704661,1.66358297369174,1.53336537800649) q[6];
cx q[6],q[2];
u1(2.62902134392907) q[2];
u3(-1.60906287698936,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.902858987817902,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.41243856515551,-2.17017941795045,2.44665538572471) q[2];
u3(1.40514280212440,-3.33754041303605,-2.03917008183065) q[6];
u3(1.60191149426509,-2.26875283350029,-0.763926085784712) q[12];
u3(1.89888659825456,-5.04677142302051,-0.886389728996721) q[3];
cx q[3],q[12];
u1(4.38726174350449) q[12];
u3(-3.46455468948484,0.0,0.0) q[3];
cx q[12],q[3];
u3(-0.373338909566020,0.0,0.0) q[3];
cx q[3],q[12];
u3(2.05655689947250,-0.493271420872126,2.48545863019841) q[12];
u3(1.35874942801437,1.88264370367572,4.18376009233263) q[3];
u3(1.17887816846784,0.149058005934573,-1.45206828516225) q[5];
u3(2.46231105664895,-3.42168334735658,1.83694889955835) q[11];
cx q[11],q[5];
u1(2.43865688092332) q[5];
u3(-2.11886681163990,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.26750535349688,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.95918393554839,0.629138205565648,1.67382199390914) q[5];
u3(1.31229953398658,-0.273965276003370,-5.65077550262784) q[11];
u3(1.89688623380133,2.85452346353005,-2.07688309685845) q[4];
u3(0.262954672101245,1.37511031297879,-0.124028195032529) q[9];
cx q[9],q[4];
u1(3.37698499831047) q[4];
u3(-1.42275541470286,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.26711524050292,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.978024352473354,-3.00978253241912,-1.14102003467453) q[4];
u3(1.50169840400949,-2.05424309454803,3.87916645307133) q[9];
u3(1.07704283395596,2.68334361213169,-1.11156702313336) q[8];
u3(0.829647742725597,1.68965970279484,-2.38438843050757) q[7];
cx q[7],q[8];
u1(1.47712030994268) q[8];
u3(-3.33433179010119,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.10698496526026,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.10606490765567,0.746809870124137,0.178854851755355) q[8];
u3(1.04030030264436,1.96864171414109,-3.60870097607889) q[7];
u3(0.944004128464313,0.437948606115302,2.67386985237007) q[0];
u3(1.61775365235975,-1.58387264618875,-1.60479717392424) q[10];
cx q[10],q[0];
u1(2.82161914510071) q[0];
u3(-1.59317109403852,0.0,0.0) q[10];
cx q[0],q[10];
u3(0.126726052509570,0.0,0.0) q[10];
cx q[10],q[0];
u3(0.314278178468378,-0.773745823504685,1.50022750155957) q[0];
u3(1.95540163018651,-0.955005703200723,2.56715532618407) q[10];
u3(1.59687487561870,0.895347418801590,0.100360832729846) q[3];
u3(0.481210944103371,-0.549865875009611,-1.74412601356812) q[9];
cx q[9],q[3];
u1(1.16044144727192) q[3];
u3(-0.354112248112681,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.52736704553460,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.534841648471819,3.41251486708998,-1.93540695856088) q[3];
u3(1.94721252775555,-5.06655928021562,-0.369554021019216) q[9];
u3(1.99798028078229,0.436760355727928,2.10494663479407) q[6];
u3(1.74089764164869,-1.97667992897733,-2.57617117645932) q[1];
cx q[1],q[6];
u1(1.56066806630478) q[6];
u3(-2.14418847997900,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.49027855595286,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.47342684176866,0.855689033532773,1.84045729608846) q[6];
u3(1.28775892815773,-4.86506234711900,-1.15504908201289) q[1];
u3(1.79899524457776,-0.994018021001506,-0.0179201495152029) q[8];
u3(2.30642128549352,-4.07409167219314,-1.10878865451457) q[7];
cx q[7],q[8];
u1(1.89872667474827) q[8];
u3(0.485885105045120,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.950699232477278,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.13417681202267,-3.19630394374264,2.23550392446924) q[8];
u3(1.77839434134838,-0.221024770606634,-4.41651854180884) q[7];
u3(1.54125808511696,0.997656409201704,0.511050531409934) q[12];
u3(1.02520218989314,-1.13060296318509,-2.52792112797579) q[4];
cx q[4],q[12];
u1(0.888714761771909) q[12];
u3(-1.43673997707840,0.0,0.0) q[4];
cx q[12],q[4];
u3(3.42408795883688,0.0,0.0) q[4];
cx q[4],q[12];
u3(1.18050442112704,-0.398539146885103,-0.306339084351726) q[12];
u3(0.336386026292907,3.50141574136440,0.309793869754404) q[4];
u3(1.48063448052559,-2.07661529629818,0.841281155374702) q[11];
u3(1.47046845311890,-2.43557817698959,0.706743270408386) q[5];
cx q[5],q[11];
u1(-0.250825730072091) q[11];
u3(0.976804081379498,0.0,0.0) q[5];
cx q[11],q[5];
u3(3.73719750819401,0.0,0.0) q[5];
cx q[5],q[11];
u3(2.74157544279069,1.67441809512225,-4.20373332743462) q[11];
u3(1.58245039677905,-1.19939733125747,0.102085921612608) q[5];
u3(0.690740901170918,0.947234457220221,-2.42702406864719) q[1];
u3(1.66943333990739,-2.32450325159221,3.21799412442753) q[7];
cx q[7],q[1];
u1(3.52189099270420) q[1];
u3(-1.57941303015520,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.30284016502365,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.63172405640075,1.19799887466619,-2.52544776433173) q[1];
u3(1.86828062170191,1.72264577741211,-0.336968219596771) q[7];
u3(1.78366390545945,2.51258588970891,-1.21652612122481) q[5];
u3(1.22542846378369,1.17079163289296,-2.45619486163814) q[4];
cx q[4],q[5];
u1(-0.0482115964342444) q[5];
u3(-0.767051377297074,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.05060235616664,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.04896625125680,-2.41922284789756,1.75819318480280) q[5];
u3(1.57219123484194,-2.23497045901349,-3.31140976173209) q[4];
u3(0.802867086112134,2.30630376589555,-3.46189894127180) q[12];
u3(1.16233999058314,2.33638574091382,-3.77823927530949) q[3];
cx q[3],q[12];
u1(1.69325820244059) q[12];
u3(-0.547610572651895,0.0,0.0) q[3];
cx q[12],q[3];
u3(-0.151136344141954,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.88323903022969,-1.60259623966230,-0.671812653982729) q[12];
u3(2.65472586429230,-2.59140470951763,1.75584892720343) q[3];
u3(0.987877879695404,0.762034882206673,-1.58180260639090) q[9];
u3(1.04618545659181,0.855135213534607,-3.69986323341469) q[10];
cx q[10],q[9];
u1(2.98668091666425) q[9];
u3(-2.11005196922700,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.209838111772321,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.638968054159744,2.54803520581067,-1.96047774642641) q[9];
u3(0.488498433577176,1.10859217510217,-0.770860186686120) q[10];
u3(1.01460559719246,1.01125232853354,1.48480936876365) q[2];
u3(0.648390013888624,-1.15935829148572,-2.10923867730260) q[0];
cx q[0],q[2];
u1(1.36478288100972) q[2];
u3(-0.464435999875108,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.0576418972307047,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.14784676679911,0.615981509359639,-0.519242386637372) q[2];
u3(1.60891979729463,2.34805413508089,2.66530880754573) q[0];
u3(1.83376792896307,0.394292875304174,2.69352259477289) q[6];
u3(0.786738886466978,3.19925622319511,2.88188115137973) q[8];
cx q[8],q[6];
u1(1.25273972819916) q[6];
u3(-0.747015007329270,0.0,0.0) q[8];
cx q[6],q[8];
u3(3.09476434778809,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.52306188144512,1.21325639381323,-0.395528482813886) q[6];
u3(1.76210799067227,1.05759631310670,4.23640777448883) q[8];
u3(1.97816745104717,-0.999924255518177,0.722177121863623) q[4];
u3(1.71258710421575,-1.54091237257214,-1.14893056705545) q[3];
cx q[3],q[4];
u1(1.73272765039343) q[4];
u3(-3.11842182280594,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.03419918678656,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.53002988537642,2.01744778477072,0.907342574746518) q[4];
u3(1.14070043200231,-1.83943580236484,3.55888753315836) q[3];
u3(0.914668124987688,-1.24230667480377,1.90804809883389) q[6];
u3(0.741793577051431,-1.65263852377442,-1.23386584935074) q[8];
cx q[8],q[6];
u1(1.48264456647364) q[6];
u3(-0.360090701440064,0.0,0.0) q[8];
cx q[6],q[8];
u3(1.18984264603288,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.89413617363301,-0.475909959066930,1.43014923835791) q[6];
u3(1.88443951839336,0.113602764991964,3.59860626366173) q[8];
u3(1.27644591013559,0.793014105184220,-0.282930445549283) q[9];
u3(0.921196607536481,0.320924666887379,-4.44439608238755) q[1];
cx q[1],q[9];
u1(2.08577548142517) q[9];
u3(0.174731874256029,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.750863447977325,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.00455379423439,-1.47294769751690,-0.281361230297785) q[9];
u3(1.77415917880988,1.96853295084035,-3.72632766456252) q[1];
u3(1.59676564667844,3.45594085339840,-0.424957174529394) q[7];
u3(1.43215872500943,1.59011493078047,-1.35284514617093) q[5];
cx q[5],q[7];
u1(3.23473445601364) q[7];
u3(-2.43162470896159,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.30496718412459,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.76451631867147,-2.49438115680108,1.75330712890545) q[7];
u3(1.24568100758324,0.737498104747388,-0.632648736140165) q[5];
u3(1.30689760502872,-2.14452231770453,-0.199697330982318) q[11];
u3(0.564061782037503,-3.22729257897426,0.773128845993350) q[10];
cx q[10],q[11];
u1(0.984988181918646) q[11];
u3(-3.55832094729842,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.94114442151206,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.26603307442873,0.843271329281480,-1.39201590418691) q[11];
u3(1.90688303493012,-0.327831579444463,1.59876987294221) q[10];
u3(0.870213423843009,-2.86500191437281,2.91082857762910) q[2];
u3(1.02927578821657,-3.01538879274893,1.39992068804894) q[12];
cx q[12],q[2];
u1(3.29125700517220) q[2];
u3(-1.22202971953989,0.0,0.0) q[12];
cx q[2],q[12];
u3(2.31354208259176,0.0,0.0) q[12];
cx q[12],q[2];
u3(2.24920770550674,3.17941872578490,-1.28304093717114) q[2];
u3(1.48666528124965,0.365324812200004,0.650114865179527) q[12];
u3(1.20065081729416,1.26420513234133,-3.81207176649803) q[2];
u3(1.82153513635740,2.25576092635887,-3.18464340895081) q[4];
cx q[4],q[2];
u1(0.235204294160871) q[2];
u3(-1.06150982063549,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.84760309914164,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.91968341961570,-1.73655373173019,1.41826050857761) q[2];
u3(2.30487570124496,-0.249473604198429,5.75125153116093) q[4];
u3(1.16759404498997,-0.520972792757839,2.36677435670144) q[1];
u3(1.85663603898834,-1.99911112257326,-0.971168942016742) q[9];
cx q[9],q[1];
u1(3.49589693789316) q[1];
u3(-4.21133278816799,0.0,0.0) q[9];
cx q[1],q[9];
u3(-0.638886510945853,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.828083388369675,-1.04074313433513,-0.882457297022030) q[1];
u3(1.59446707411165,-4.35071513198451,1.49049263641132) q[9];
u3(1.32673243704843,2.93005120899928,-2.28244103308725) q[0];
u3(0.980978303001466,2.96645440301223,-3.27094648513055) q[3];
cx q[3],q[0];
u1(3.08910583127952) q[0];
u3(-2.38253109708138,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.52110140723368,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.83967391303596,0.466535132147348,-1.08253803107746) q[0];
u3(1.43071985816782,-0.842695311499783,-4.48481495443562) q[3];
u3(1.39802543659240,2.70330634686197,-1.81048758594371) q[5];
u3(0.918495227893987,-3.25909273430632,2.72403726385225) q[12];
cx q[12],q[5];
u1(0.373603440446656) q[5];
u3(-1.45332517854724,0.0,0.0) q[12];
cx q[5],q[12];
u3(1.82477413053913,0.0,0.0) q[12];
cx q[12],q[5];
u3(1.43855350710295,-2.29532637112755,1.77064095026347) q[5];
u3(1.03326902500158,-2.28147277887354,-0.797421962232703) q[12];
u3(1.16518603181027,-2.10978866956733,3.19583924998046) q[7];
u3(1.96090563142034,2.48830458640410,-1.60096640485085) q[11];
cx q[11],q[7];
u1(2.07444963629208) q[7];
u3(-2.48027724819242,0.0,0.0) q[11];
cx q[7],q[11];
u3(3.18638418745755,0.0,0.0) q[11];
cx q[11],q[7];
u3(0.947097849810472,-0.325714236656240,-2.15591098614463) q[7];
u3(1.48405597219659,5.47360204681787,-0.561784993260795) q[11];
u3(2.37621797191613,2.33745503976061,-0.569542171602109) q[10];
u3(2.02442090858151,0.754843255497896,-2.31265215049799) q[8];
cx q[8],q[10];
u1(1.52882419467246) q[10];
u3(0.411515075310471,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.807189550877769,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.00415556733032,2.00744810970248,-3.44496047268617) q[10];
u3(1.76688735365718,-0.965372647540095,-0.287003037815165) q[8];
u3(1.75783745024200,1.38750000347012,-4.49429760136905) q[2];
u3(1.15251337462758,-1.73805400296831,4.17242458350146) q[12];
cx q[12],q[2];
u1(0.402142988953970) q[2];
u3(-1.14288152115808,0.0,0.0) q[12];
cx q[2],q[12];
u3(1.63867340972571,0.0,0.0) q[12];
cx q[12],q[2];
u3(2.77762345452509,0.897149190057200,-3.87807495651044) q[2];
u3(1.58714875403063,0.703851314965320,-5.33137995927828) q[12];
u3(2.65517288323493,2.92186262193621,-2.39260440792415) q[8];
u3(0.846805922171812,0.209969323476953,1.32614472231156) q[11];
cx q[11],q[8];
u1(1.64496151896463) q[8];
u3(0.0665791968226781,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.634869625675998,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.77458651727619,2.62360928130120,0.165741275341773) q[8];
u3(1.21653777301554,1.86768543987788,1.82830601121706) q[11];
u3(1.57111542637549,0.552333368000649,-1.08295820837439) q[1];
u3(1.58839781712033,1.41186355403438,-4.48594300181109) q[5];
cx q[5],q[1];
u1(0.735166292450216) q[1];
u3(-3.16070569470440,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.12886790051766,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.58156041745794,-3.31022811302524,2.12032153825165) q[1];
u3(0.396595757803086,3.25038180416355,0.444936048255244) q[5];
u3(1.93880024392851,1.00144807423900,-3.38642324292681) q[6];
u3(1.04082964105988,-2.27940434826003,2.91132301306543) q[9];
cx q[9],q[6];
u1(2.90591472697366) q[6];
u3(-1.82212634821435,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.911031818486223,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.20705641308079,-0.430388909897489,-2.82336097869298) q[6];
u3(2.03301120420975,-4.36822603772580,-1.34117302904896) q[9];
u3(0.656636205687342,1.93323212779582,-2.70415057952463) q[0];
u3(0.423707866380698,0.582596006458676,-1.81316277199659) q[10];
cx q[10],q[0];
u1(3.81363560818406) q[0];
u3(-3.33718007659279,0.0,0.0) q[10];
cx q[0],q[10];
u3(-1.09160352665000,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.70538137999789,-2.32281548445285,3.35765662973056) q[0];
u3(1.58626312856902,-5.86563454793169,-0.336171746652712) q[10];
u3(1.65893520600284,2.71746262339904,-2.27803118617503) q[7];
u3(2.97767025502610,0.411572387518816,-1.97260240841088) q[3];
cx q[3],q[7];
u1(1.32322013221384) q[7];
u3(-3.55038223723438,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.12419146918683,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.16908510781681,0.223572199721127,-1.56653152458562) q[7];
u3(1.56073192467941,-1.14698996671834,4.64411339948770) q[3];
u3(2.21166356722865,0.973206570581518,-3.15576467550221) q[11];
u3(0.966177863204746,2.55462990597202,-2.79471583866467) q[0];
cx q[0],q[11];
u1(2.18930599403387) q[11];
u3(0.485952430727689,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.45834900564098,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.35426189371197,3.57646818702862,-1.52157991749034) q[11];
u3(1.63909308674014,2.86159062703831,1.56849669627535) q[0];
u3(1.25226247180440,-0.737961598842262,2.05470036380538) q[3];
u3(1.16919586449359,-1.43018712931790,-2.40675077028380) q[5];
cx q[5],q[3];
u1(-0.0305678712986475) q[3];
u3(-2.55748342923620,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.66194422367172,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.42833236435245,2.28316397459871,-1.71056618253689) q[3];
u3(1.00503522746278,-3.55669147396988,2.33646542205971) q[5];
u3(2.05273961418179,-0.675063532594679,0.152912316680320) q[1];
u3(2.15880668633721,-2.52347853373953,0.0375141977703366) q[8];
cx q[8],q[1];
u1(1.51253627103506) q[1];
u3(-3.15960115184644,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.11829245142571,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.34707376371390,1.01045264415726,-3.23384040693188) q[1];
u3(1.93400573798065,2.34925002697608,-0.476506903874236) q[8];
u3(0.817823102202497,1.54062177142505,-3.81013420981085) q[9];
u3(1.14126141391497,2.39426510797679,-3.00892662066719) q[7];
cx q[7],q[9];
u1(2.55971864295937) q[9];
u3(-1.97927638680374,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.702445806598719,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.348914380715244,2.54049617566017,-3.30519996471446) q[9];
u3(0.837297480726754,-0.509343834707342,-3.25350112872684) q[7];
u3(1.50167573890654,1.68139247385148,-3.55803340264324) q[12];
u3(2.44303504433852,3.96380232062514,-2.25894944037250) q[6];
cx q[6],q[12];
u1(2.83495825946325) q[12];
u3(-2.53792306099087,0.0,0.0) q[6];
cx q[12],q[6];
u3(1.66688930730201,0.0,0.0) q[6];
cx q[6],q[12];
u3(2.00253901486732,1.27139831391410,-1.79783218497833) q[12];
u3(2.68891844076418,-5.29527627361874,0.903083307915614) q[6];
u3(1.90552580971776,-1.76661329492828,0.379559403150334) q[10];
u3(2.88463058645758,-1.10956382079960,-0.0821353808034114) q[2];
cx q[2],q[10];
u1(2.26090740929638) q[10];
u3(-1.61992665044006,0.0,0.0) q[2];
cx q[10],q[2];
u3(3.61305725089814,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.21630518467781,2.42735368381905,1.22117408082373) q[10];
u3(1.90827808568875,2.69788085520877,3.29438626476181) q[2];
u3(2.43365273833453,3.63807507566805,-0.865918249508386) q[10];
u3(1.43989278205306,1.91057963284924,-2.59510121359228) q[5];
cx q[5],q[10];
u1(1.11196910258740) q[10];
u3(-0.302549708239193,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.80324964718010,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.672810574378514,-1.10626977221175,4.48617010830048) q[10];
u3(2.48722771610984,0.749499813116114,-4.96901918243876) q[5];
u3(1.60997387396961,0.195268549251299,2.18220323260566) q[12];
u3(0.851193353465989,-0.676795752883007,-1.50105600406266) q[11];
cx q[11],q[12];
u1(1.16315052706223) q[12];
u3(0.0410477047377122,0.0,0.0) q[11];
cx q[12],q[11];
u3(2.64063854511811,0.0,0.0) q[11];
cx q[11],q[12];
u3(0.810088816421477,-2.09842075258322,3.44305759519243) q[12];
u3(0.905488165176204,-2.17019190971969,1.64438506187946) q[11];
u3(2.46084877179670,1.85130148818900,-0.871546121369902) q[8];
u3(2.63971925494098,1.03918071992457,-3.07588237347002) q[4];
cx q[4],q[8];
u1(3.23812330273145) q[8];
u3(-1.51911918876382,0.0,0.0) q[4];
cx q[8],q[4];
u3(2.66291444292699,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.248164182464446,4.03862820277643,0.0331055035370325) q[8];
u3(1.65024024293592,-0.157469621397671,-5.54672182735572) q[4];
u3(2.47919840011649,3.03714754816802,-2.99004134049836) q[3];
u3(2.08433754040950,2.87079463421409,-3.26224275642234) q[1];
cx q[1],q[3];
u1(2.02011217190565) q[3];
u3(-1.56525004999326,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.140632933458710,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.16282234298865,1.63479655486599,-2.98627648412001) q[3];
u3(0.453139776422224,1.67729532508834,-1.75326656068702) q[1];
u3(1.01109843859681,-0.418296549933871,1.02485855283876) q[9];
u3(1.05502163152590,-0.569890796818285,-1.29710611311740) q[7];
cx q[7],q[9];
u1(0.217080864882673) q[9];
u3(-1.04227063326590,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.55538867740323,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.60867170354342,-0.317912389972031,-1.53705568451148) q[9];
u3(1.55492896816716,-2.21685291769139,-0.706695041983369) q[7];
u3(1.53078524757176,2.14658660985956,-2.90918358159808) q[0];
u3(2.42645052764160,2.11885779737040,-3.69123338009652) q[6];
cx q[6],q[0];
u1(0.0642283077788506) q[0];
u3(-0.382762990480790,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.35892745259250,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.76689864171966,-0.662770910850198,1.56053590128268) q[0];
u3(0.245825680950258,-1.97022162215118,-1.45448102645402) q[6];
u3(1.69194962793700,1.81113916896831,0.00863384465111450) q[0];
u3(0.522610215719371,-0.723237641676983,-1.59193149584539) q[4];
cx q[4],q[0];
u1(-0.138501917772377) q[0];
u3(-2.53110790959414,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.23719396674426,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.14609364800492,0.753228992800711,-0.317509050962759) q[0];
u3(2.18214565012161,1.79536138218948,0.250487986360808) q[4];
u3(0.909532793660469,2.95449973732074,-0.953594516983251) q[12];
u3(1.13741108198606,1.67484213289375,-1.51505173647779) q[2];
cx q[2],q[12];
u1(1.23420826804063) q[12];
u3(-3.39986850580781,0.0,0.0) q[2];
cx q[12],q[2];
u3(2.18707497953166,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.54608703802504,0.729891837845439,-3.07792561088837) q[12];
u3(0.425279052610229,-2.36290549000449,-1.32803695726690) q[2];
u3(1.17033728621751,-0.216978436682862,-2.21224557474976) q[8];
u3(1.56847015214807,0.357137457464916,-4.89044235584578) q[11];
cx q[11],q[8];
u1(2.49339899918848) q[8];
u3(-1.75891743708558,0.0,0.0) q[11];
cx q[8],q[11];
u3(3.31929603439586,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.18195193006431,-3.28140541150337,-0.113172872792964) q[8];
u3(1.68407199266688,4.96526816644327,0.679794874128645) q[11];
u3(2.53935322749084,-2.38640471761948,3.22420121707944) q[7];
u3(0.916674943493324,2.11985319021063,0.0155884004308577) q[9];
cx q[9],q[7];
u1(-0.151584694193652) q[7];
u3(-2.01487245742797,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.10128195240176,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.968022483382539,2.46079984287541,-3.69824434599761) q[7];
u3(1.99102817256612,-2.18030497811999,-2.55036088712681) q[9];
u3(0.918936833343693,-1.10266219052877,3.48852463629191) q[6];
u3(1.55549547168686,-1.42014390897478,1.50233529198968) q[1];
cx q[1],q[6];
u1(1.61623044894115) q[6];
u3(0.194769082223729,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.867934797832891,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.77400476664787,1.67142740558178,-4.08243378657929) q[6];
u3(0.376701644467903,1.07393646392237,-4.97088905730233) q[1];
u3(1.73521244812170,1.76589152593172,-0.346337182566973) q[5];
u3(0.618928678911519,0.840801523736621,-4.05739150947696) q[10];
cx q[10],q[5];
u1(2.29532013728531) q[5];
u3(0.115570070341555,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.46831056453930,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.46840520442810,1.73273802550074,-4.07659266750243) q[5];
u3(0.421042408612123,3.78016685902837,0.354631852780290) q[10];
u3(0.443414480967963,-0.341637345238582,0.177772473716776) q[6];
u3(0.730876480416085,-1.27757723968969,-1.25986859618290) q[4];
cx q[4],q[6];
u1(0.787467067784335) q[6];
u3(-1.40616893755175,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.92994125380893,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.00964419510272,-1.32907434172827,-1.78065943304345) q[6];
u3(0.240219223179904,-0.308322604081067,1.73058539589358) q[4];
u3(2.27115135363562,-1.42749131218834,3.41497148470647) q[0];
u3(2.25118337374371,-1.84237330700061,-0.275177835584407) q[3];
cx q[3],q[0];
u1(1.09822145850792) q[0];
u3(-0.227746586962018,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.68515860585142,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.58186645424787,-0.871790675526733,2.08103895550724) q[0];
u3(1.50897428172007,-1.43085701305929,-1.86357194123503) q[3];
u3(1.67663268405562,-0.503203053237850,2.03646464896362) q[1];
u3(1.35729340563079,-0.700562389320325,-1.01300748605539) q[2];
cx q[2],q[1];
u1(1.45098633083311) q[1];
u3(-3.03620858629264,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.280009152027248,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.869621714685654,2.04549452797444,-0.0259612163329155) q[1];
u3(1.63320687316023,3.66207716110435,2.60441329434486) q[2];
u3(2.26153089228714,-1.62843481656312,-1.47350444416753) q[10];
u3(1.89383813767662,-2.57992421564508,0.345971861852708) q[5];
cx q[5],q[10];
u1(0.00291443346702591) q[10];
u3(-1.23459330442194,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.66142260905188,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.74660083105284,1.34017576843927,0.986793845867673) q[10];
u3(1.24758275772558,1.93210265382247,3.73589261449529) q[5];
u3(1.61622784801460,3.22075658184088,-2.35113276530523) q[8];
u3(0.995306739892649,2.32444468645340,-2.08580102163406) q[12];
cx q[12],q[8];
u1(2.89220997329220) q[8];
u3(-1.49272932089680,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.23407123182923,0.0,0.0) q[12];
cx q[12],q[8];
u3(1.37279634664749,1.75978478461831,0.534725777974633) q[8];
u3(0.985741323695631,1.36030438532034,-2.48777527313582) q[12];
u3(0.836623806757436,2.77015403877407,-1.65229476853420) q[9];
u3(0.520150131840500,-0.163997205253843,-1.24534995565888) q[11];
cx q[11],q[9];
u1(1.41350128346748) q[9];
u3(-0.0231077738275549,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.31288374419411,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.11282498862305,1.30160708518070,-1.46234926126584) q[9];
u3(1.90607497539152,2.38111850517375,1.79671217360434) q[11];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12];
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
