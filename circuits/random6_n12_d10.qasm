OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.10072175165044,0.642970878740077,0.621159897986718) q[3];
u3(0.422327842131152,-2.54474305509230,-1.83954246666368) q[10];
cx q[10],q[3];
u1(1.49702782045527) q[3];
u3(-0.428819222891600,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.25920786971338,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.09911440127969,0.0850413904722387,3.45155873339953) q[3];
u3(1.49052977026468,4.18338338119087,-0.633396393621459) q[10];
u3(1.99475958951502,1.14866361333262,-3.78548713572523) q[4];
u3(1.68909584471634,-1.63638850178002,3.70406507083171) q[2];
cx q[2],q[4];
u1(3.16736090010766) q[4];
u3(-4.34097931894000,0.0,0.0) q[2];
cx q[4],q[2];
u3(-0.344512271765638,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.62116385139376,3.09339832648471,0.684425339571745) q[4];
u3(2.21486229734150,0.477698656732461,1.15564071107549) q[2];
u3(0.881870554523622,3.31925685957354,-1.02634356168745) q[7];
u3(1.38916074328633,2.17046330960893,-1.68314685829085) q[6];
cx q[6],q[7];
u1(1.45539212314940) q[7];
u3(-0.411442799235693,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.95108447657422,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.76309073353326,1.96776952483183,-3.59128492366113) q[7];
u3(1.62838808899136,-4.67833973925343,0.00700746368318628) q[6];
u3(1.53487860070800,1.47656424842735,-3.09422283661699) q[0];
u3(0.219122945690991,1.43032726493291,-2.19505528905612) q[11];
cx q[11],q[0];
u1(2.17583790462469) q[0];
u3(0.0397026901682183,0.0,0.0) q[11];
cx q[0],q[11];
u3(1.67461116465292,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.06131072748935,0.131385797135990,0.610541257563452) q[0];
u3(2.65305807082771,-1.30393567681759,-1.55380363781557) q[11];
u3(0.204269571017474,0.529886290624810,-1.48225478894739) q[8];
u3(1.55990086711913,-3.69147810022053,1.82939674837617) q[9];
cx q[9],q[8];
u1(1.46051752306661) q[8];
u3(-0.595598772122559,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.86441857478521,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.48717687589983,-0.327624559242455,1.86125396725526) q[8];
u3(1.74512885245875,1.17441625113538,-1.97524120924240) q[9];
u3(0.661338782879606,-2.02269383726462,1.97108633496220) q[1];
u3(0.240067237664944,-4.03682545678014,2.12774277066352) q[5];
cx q[5],q[1];
u1(1.70766752822182) q[1];
u3(0.134071378024890,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.753524778252649,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.42069456038298,-0.484609801205385,2.03586478982281) q[1];
u3(1.49694677517780,-2.69249470379195,-0.462380464090673) q[5];
u3(1.97495069371361,2.02702152383329,-2.94730284071713) q[11];
u3(0.837848067332031,1.77226394529923,-2.16319848070823) q[6];
cx q[6],q[11];
u1(1.82578480094532) q[11];
u3(0.249468668987473,0.0,0.0) q[6];
cx q[11],q[6];
u3(0.862375808916918,0.0,0.0) q[6];
cx q[6],q[11];
u3(2.09587424219276,-2.32083683831571,-0.674943060900565) q[11];
u3(2.55181379099195,-2.60834521721371,3.23146022763797) q[6];
u3(0.126025176017947,0.338812089935660,-2.63958097307218) q[10];
u3(1.47769263553557,2.54837693436118,-2.97894313225943) q[9];
cx q[9],q[10];
u1(1.69809229482560) q[10];
u3(-2.57144505136493,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.150660713111457,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.26646826587055,2.26618957024571,0.403781833747498) q[10];
u3(2.45397196596911,2.56782362479397,-2.26083013700990) q[9];
u3(0.850249982271015,1.31257886389532,-1.57508277771809) q[5];
u3(0.534482826834888,0.0288875725693637,-1.83550673948225) q[4];
cx q[4],q[5];
u1(2.33178032158604) q[5];
u3(-1.63793844297995,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.338868841565608,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.269538639513242,-1.68462902641740,1.76744111412168) q[5];
u3(2.23894523656135,0.975387922890040,-0.690543656567122) q[4];
u3(1.98120988400741,0.376471992111821,-0.984010439876485) q[3];
u3(1.11375520385537,-3.67164336565974,0.780410545255151) q[8];
cx q[8],q[3];
u1(0.707617719409040) q[3];
u3(-1.27569060154971,0.0,0.0) q[8];
cx q[3],q[8];
u3(-0.242448613072036,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.01564427062443,-1.97582576630923,2.71139495876344) q[3];
u3(1.40928726791471,1.62923185304925,-1.60421399628245) q[8];
u3(2.29775962921142,-1.00611094082799,1.88285991187252) q[1];
u3(1.70879318771242,-1.35167766913131,-1.00193660622127) q[2];
cx q[2],q[1];
u1(4.02788274714408) q[1];
u3(-1.54121609183898,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.99285424354087,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.75222569259218,3.43372946504843,0.529125941071832) q[1];
u3(0.582685355195862,2.53726380988505,-2.83362189134967) q[2];
u3(1.13523196182012,-0.486443016962259,1.46543408213853) q[7];
u3(1.63684822816785,-1.70763457028633,-2.03986693264218) q[0];
cx q[0],q[7];
u1(1.79061996886441) q[7];
u3(-2.75621830335908,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.413683042363526,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.32765976402713,1.65784260749402,0.462145867817050) q[7];
u3(0.626869140512995,4.20295535251893,0.378523703811736) q[0];
u3(1.11024850198854,-0.266487893235997,-1.02451291270909) q[10];
u3(0.818382790527626,-3.22121039086642,0.990718885721823) q[3];
cx q[3],q[10];
u1(1.45413013253583) q[10];
u3(0.155979898682878,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.18525593647901,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.28610694484175,1.79664328153083,-2.89589628534285) q[10];
u3(1.03285157809010,0.215403224125128,-3.04727551748807) q[3];
u3(2.99823629290600,-3.11123653582710,0.138172443547479) q[0];
u3(2.85261238063229,2.21842746569080,2.64625230580761) q[4];
cx q[4],q[0];
u1(1.75528113506371) q[0];
u3(-2.93113420873799,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.61490461561687,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.79537244503969,-2.39397926558451,2.72034969922681) q[0];
u3(2.09729101862363,2.85313843483566,-2.83311444243160) q[4];
u3(1.28958841055704,-0.684139975518097,2.04507809724107) q[7];
u3(1.25590600243901,-1.64968362541470,-1.06842669572362) q[2];
cx q[2],q[7];
u1(1.50715430251466) q[7];
u3(-3.15272319288915,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.32229983637093,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.43263532392143,0.406415381913024,0.789168036060866) q[7];
u3(2.51749621829814,0.807012351006065,-0.505303032464955) q[2];
u3(2.15336392474953,-1.33787589370099,-1.03155859584140) q[1];
u3(0.601813222764662,-4.25110337619390,-0.277637112897713) q[9];
cx q[9],q[1];
u1(2.07900097689455) q[1];
u3(0.270558124483745,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.26353313658871,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.20294284892132,-2.10362524775867,0.289124064932740) q[1];
u3(1.35650141882507,1.43803772115593,2.72582051688516) q[9];
u3(1.69999543075402,0.715230901523664,1.47085114917896) q[6];
u3(1.72714571105390,-1.37457060451483,-2.40164833683001) q[5];
cx q[5],q[6];
u1(-0.0830312046495028) q[6];
u3(-1.71721900484598,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.721327316743440,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.21759382142359,-1.24937208299739,3.21028166356160) q[6];
u3(1.62679587649528,1.19761067665047,4.89251095717389) q[5];
u3(0.973790084056395,-2.92205092892908,2.63519609360858) q[8];
u3(0.779007577280411,-0.397205680167905,-2.50668490307578) q[11];
cx q[11],q[8];
u1(2.92742681660906) q[8];
u3(-1.17048697658401,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.462305111632238,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.41657693924179,0.0275202075540691,-0.414345381330232) q[8];
u3(2.21688525568808,4.83704058048917,1.13653028944737) q[11];
u3(1.98974332362602,2.14494167404652,-2.61187881249941) q[1];
u3(1.85420991234769,-2.93924990136677,2.76986211781772) q[0];
cx q[0],q[1];
u1(2.63072567322529) q[1];
u3(-1.42735547555837,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.58863397174589,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.63927115472387,2.49764822352585,-0.554424901680325) q[1];
u3(0.757705764045656,0.521901375828038,-0.210305993585986) q[0];
u3(1.85292061488226,1.89057907507673,-2.18644316094374) q[6];
u3(1.21532749342290,-2.56931824646819,2.25384720707957) q[10];
cx q[10],q[6];
u1(2.50134271069807) q[6];
u3(-3.00635191174872,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.12828445844781,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.94472669940663,-1.76465487216638,3.86227355664382) q[6];
u3(1.11038795337887,-1.25059428374507,4.47907388872712) q[10];
u3(1.10198962868990,-1.07091628270206,0.883706255738973) q[9];
u3(1.77566334545908,-1.39130205071841,-2.27611368046405) q[2];
cx q[2],q[9];
u1(-0.0515387722340683) q[9];
u3(-1.97997320254580,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.69906613806281,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.02301981962102,-0.831456989696199,2.27001908167491) q[9];
u3(2.27301149923255,-0.320835329954217,-4.97618468338416) q[2];
u3(0.968256556165262,-1.05986666313815,1.74580955559397) q[7];
u3(1.31647120826899,-1.37515465530714,-2.09228891749208) q[5];
cx q[5],q[7];
u1(2.24325450112657) q[7];
u3(0.153448539255423,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.35129436013235,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.99408786215014,-0.715134046411000,1.64247483937951) q[7];
u3(0.763394362834077,-1.00940670710827,-4.18362031211304) q[5];
u3(1.28905512825837,1.84946658692325,-3.29290382617476) q[3];
u3(1.54564488564934,-3.11307712496724,3.15281481066806) q[4];
cx q[4],q[3];
u1(0.634621774913181) q[3];
u3(-1.47163739420396,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.211902349544095,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.23089425132101,-2.35629137818123,0.155335037647253) q[3];
u3(2.01489895301875,0.137112920175171,4.76530664961418) q[4];
u3(2.45375608573538,-1.93652920540675,3.79857770956477) q[8];
u3(1.01315791475223,1.32188769536489,0.708211510443498) q[11];
cx q[11],q[8];
u1(2.85771506025186) q[8];
u3(-1.75497821488078,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.752204243727441,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.93444930051081,1.07840989886879,-0.231107393003001) q[8];
u3(0.964127042089046,4.10752226875855,-0.321528498058635) q[11];
u3(1.82653920118794,2.20690881374495,-2.21745992362109) q[3];
u3(1.05008596379790,1.88960457778339,-2.86052914298415) q[8];
cx q[8],q[3];
u1(0.262099044373505) q[3];
u3(-1.94770795005886,0.0,0.0) q[8];
cx q[3],q[8];
u3(0.359782095201713,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.452309037313985,1.68614107026668,-1.88245015830739) q[3];
u3(1.46802085717107,-3.62491003392336,-0.386472046508278) q[8];
u3(2.37909637332922,1.35731920058901,-3.35050583757279) q[4];
u3(2.05986047249229,2.80802977118581,-2.85068375776546) q[9];
cx q[9],q[4];
u1(3.17919488343074) q[4];
u3(-1.63027261760493,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.56198100163289,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.37166656712769,3.27543365510369,-1.02967232375998) q[4];
u3(2.02095165570959,-2.09183376883286,-0.723026598499036) q[9];
u3(2.55521222986801,1.59997538074434,-0.451103805978540) q[0];
u3(2.39339696683306,2.12256875822321,-2.42051042722223) q[10];
cx q[10],q[0];
u1(0.733879348859448) q[0];
u3(-3.30108316868989,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.06706193513609,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.93944572034815,-0.290071526237378,2.85651454808737) q[0];
u3(2.63457060977740,-1.14227492905467,0.257334344489096) q[10];
u3(2.73866956817829,3.04292200815213,-1.27513135566219) q[2];
u3(1.24931313149205,2.45446823229497,-0.857664799931843) q[11];
cx q[11],q[2];
u1(0.928149683926646) q[2];
u3(-0.350019011715319,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.39086170816108,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.80368378125920,-2.93727065931424,2.86623222839414) q[2];
u3(0.652996338656682,2.96566052778691,3.08189502901915) q[11];
u3(0.942559734875640,2.45646081657679,-0.336656260829292) q[5];
u3(0.935513557675105,-0.00486665660722152,-2.56289972115855) q[1];
cx q[1],q[5];
u1(3.32052153871861) q[5];
u3(-0.775778692712314,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.86138632148586,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.23541480877876,-3.75146603727566,2.44965850266816) q[5];
u3(1.74847641053021,-0.899577952172276,-3.26535505841878) q[1];
u3(1.48954346638749,-1.94341217275731,0.586001201204864) q[7];
u3(1.92210307057069,-4.07285260343419,-0.194847733961186) q[6];
cx q[6],q[7];
u1(0.0576941173106593) q[7];
u3(-1.41649699062538,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.66569425552113,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.694079801169140,-1.50367488251941,2.31592741067444) q[7];
u3(0.432323074716047,0.347760807998731,3.12095761901247) q[6];
u3(2.03696659329572,2.42375160505014,-1.82694942265748) q[7];
u3(2.60843331036563,1.16407657356323,-0.715464936617018) q[4];
cx q[4],q[7];
u1(2.76464470079381) q[7];
u3(-2.25417865182766,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.22351749745442,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.69618551307165,2.15491644802065,-3.67243519271664) q[7];
u3(0.238457774039239,4.78121834137460,-1.43982777739162) q[4];
u3(1.39688483251071,2.33411367701184,-1.46823145263741) q[8];
u3(2.23754320838455,0.289202419843356,-3.33409488410661) q[1];
cx q[1],q[8];
u1(1.33277670569406) q[8];
u3(-1.09576466654134,0.0,0.0) q[1];
cx q[8],q[1];
u3(-0.511109447903054,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.50486039079041,-1.79802706975823,4.22814854791510) q[8];
u3(1.13585643103162,-3.81290377521500,-1.06296177324172) q[1];
u3(1.79376044212250,1.61309055917299,-2.65668780925292) q[3];
u3(1.84839947221701,2.56273345135111,-2.97884029435234) q[2];
cx q[2],q[3];
u1(1.78819979839944) q[3];
u3(0.105532853263700,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.19065515338069,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.98214410902456,-1.34848174438228,-0.827813473267463) q[3];
u3(1.66499668687754,2.19100147840623,2.94468754979843) q[2];
u3(1.60054679906742,1.39520953086499,-2.53380787243591) q[0];
u3(1.23710626124049,2.15229284168845,-3.66943567498450) q[10];
cx q[10],q[0];
u1(1.63092608072932) q[0];
u3(-3.09414807922716,0.0,0.0) q[10];
cx q[0],q[10];
u3(0.912853065669996,0.0,0.0) q[10];
cx q[10],q[0];
u3(3.07276577535488,1.90728689813493,-3.04136546830154) q[0];
u3(1.58344937500442,-1.18941880058837,-2.82313561866157) q[10];
u3(2.12013949240878,-3.33414227978983,0.271940778618952) q[5];
u3(2.96540546649877,2.23207480641323,3.72914901033231) q[9];
cx q[9],q[5];
u1(0.277288618322945) q[5];
u3(-1.28405331393320,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.48802484272361,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.986347847513850,2.83786839457276,-0.995163424936136) q[5];
u3(1.70018326527230,-3.55039121677102,0.462392586421882) q[9];
u3(1.79390075916089,-1.22560654569127,0.912266036078284) q[11];
u3(1.52327158732152,-2.25043028392611,0.391678941382600) q[6];
cx q[6],q[11];
u1(0.755386953759919) q[11];
u3(-3.52874315353726,0.0,0.0) q[6];
cx q[11],q[6];
u3(1.50633439424359,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.90914023015251,1.98769938886577,-1.92685025546663) q[11];
u3(2.26217860947261,-4.23951468302997,-0.536293539812910) q[6];
u3(1.26700875045635,1.04942726157018,-2.70626705475135) q[8];
u3(2.51549008021585,3.61148830842162,-2.41965142680398) q[2];
cx q[2],q[8];
u1(0.0991815521554429) q[8];
u3(-2.05236371336570,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.935132763242468,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.785802699413647,3.45286807536671,-1.82109960879574) q[8];
u3(1.09010901484693,1.09302531053229,1.91096943536297) q[2];
u3(1.89049106339769,-3.70868672460415,1.50246544083915) q[6];
u3(1.35650019853178,2.37192841113515,-0.254803536125469) q[5];
cx q[5],q[6];
u1(3.21300456801394) q[6];
u3(-2.04661707392557,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.20503744764844,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.22536619438393,3.80388143565534,-1.24271589349043) q[6];
u3(2.65564739092859,0.952873658341575,-3.16808468938167) q[5];
u3(1.87068803997497,0.533239893255010,-3.05901766827915) q[9];
u3(1.54282677423276,-3.52844055470187,2.14469252235956) q[7];
cx q[7],q[9];
u1(1.24303747107314) q[9];
u3(-0.202861929353681,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.21979055169450,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.57373467865989,3.32216250239208,-1.30457561966746) q[9];
u3(1.97902486552470,-4.75645796764058,-1.48613579436983) q[7];
u3(2.28272043523243,-0.594634708561782,2.00511592825043) q[11];
u3(2.81833644450932,-1.35326262286429,-0.775843106995844) q[4];
cx q[4],q[11];
u1(2.07600854272121) q[11];
u3(-0.0745388699894753,0.0,0.0) q[4];
cx q[11],q[4];
u3(0.780144376836947,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.83782885848849,-2.98475712548316,2.21510655701994) q[11];
u3(1.84773658014488,-1.17378292810940,-2.46598859459152) q[4];
u3(0.363431950548011,-3.28794583037252,2.14543657597423) q[1];
u3(0.461982050689995,0.142877014953082,-2.17414968130073) q[10];
cx q[10],q[1];
u1(1.24578411456898) q[1];
u3(-3.21886010022138,0.0,0.0) q[10];
cx q[1],q[10];
u3(2.34341839066344,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.21705861711039,-1.73825871769412,3.49656153764539) q[1];
u3(1.20204587832157,2.97516115300285,0.434581838646827) q[10];
u3(1.92375519901367,-1.62571698532095,3.85657661367282) q[3];
u3(2.73625814388926,-0.719673364681966,-1.19182447038087) q[0];
cx q[0],q[3];
u1(-0.853586282898790) q[3];
u3(0.162933926138205,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.84130385119209,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.87681044613770,-3.21097258859099,2.72349128735533) q[3];
u3(0.679687959323405,-0.307759830156340,5.64181311817558) q[0];
u3(1.44653394895220,-1.78576466045839,4.43725899769656) q[0];
u3(0.175773572502291,0.246196353116664,0.181346796341476) q[5];
cx q[5],q[0];
u1(3.17879303304422) q[0];
u3(-1.09027198435137,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.54956028044968,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.16707973349254,-2.18899461244955,1.77942443791742) q[0];
u3(2.32907956313823,-0.971975988247516,3.96766685095040) q[5];
u3(0.801895267267544,-0.846755161513537,1.03540607661588) q[7];
u3(0.723022918421188,2.63663388669559,-3.25835528576056) q[8];
cx q[8],q[7];
u1(1.41911135028220) q[7];
u3(-0.516757298180479,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.0290177642384577,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.47430570976924,-1.34356826407269,3.25527655183595) q[7];
u3(0.0168943055952913,-2.63536186199962,-2.90489554482436) q[8];
u3(1.32919858414538,1.86030120006720,0.172089826703074) q[9];
u3(2.31048934031169,1.00789453937254,-1.54160334184670) q[3];
cx q[3],q[9];
u1(-0.00946453391546065) q[9];
u3(-1.47919092944845,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.53312066978964,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.08642750745393,-1.12793961215593,-2.27190964305916) q[9];
u3(0.937973256000019,-1.66185407739112,-1.64222060996663) q[3];
u3(1.66144946344273,-0.530534971167497,0.284498195878907) q[10];
u3(1.87344985282409,-2.82038889804337,0.154529785901303) q[11];
cx q[11],q[10];
u1(2.27226689954557) q[10];
u3(-1.51142449837423,0.0,0.0) q[11];
cx q[10],q[11];
u3(3.37857037246934,0.0,0.0) q[11];
cx q[11],q[10];
u3(0.702986718499120,2.41667434080967,-3.61291532645736) q[10];
u3(1.62053322639440,0.0444122216232865,1.04962822946008) q[11];
u3(1.22507852473707,1.30196869744180,-0.681852954671393) q[2];
u3(0.951659586027554,0.839680329464897,-3.75710421203123) q[1];
cx q[1],q[2];
u1(-1.24051490345302) q[2];
u3(0.390190513351502,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.66593027256442,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.07999063441420,0.791192672642200,-2.69763092253161) q[2];
u3(2.04891334058482,4.61800700847772,-1.04478109063313) q[1];
u3(1.40837890830101,2.37094648077330,-2.51386523645689) q[4];
u3(1.93312645989982,-2.89756539361122,2.78610420222805) q[6];
cx q[6],q[4];
u1(2.46652542602891) q[4];
u3(-3.11082225159689,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.47729617453472,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.52475128730938,-3.68556344478923,-0.365864143480049) q[4];
u3(1.43648831871646,-2.22443322537658,-0.894869712634432) q[6];
u3(1.95914846688287,-0.482080595723158,0.356397529059632) q[10];
u3(1.50608977230142,-3.76096145960862,-0.409463049546385) q[0];
cx q[0],q[10];
u1(1.31740810364122) q[10];
u3(-3.35887881627984,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.19414234490936,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.04226175564860,2.68577324132212,1.32881983354305) q[10];
u3(1.96810403940103,0.922942183525932,-0.179521433405924) q[0];
u3(1.91291142572023,1.81397681538127,-3.88832088325159) q[1];
u3(0.957085598602423,1.67054444958467,-1.37754159828113) q[7];
cx q[7],q[1];
u1(0.351844174173680) q[1];
u3(-1.12100101426399,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.27652565040554,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.924924462877639,-0.734328676519289,1.21794801510854) q[1];
u3(0.971837205484744,-1.42485303462624,4.20511157847200) q[7];
u3(1.91778556682154,-0.393528830089365,1.77923273617389) q[6];
u3(2.06165060563063,-1.86266810347185,-2.71048291683619) q[4];
cx q[4],q[6];
u1(2.11507935339944) q[6];
u3(-1.73961113450594,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.22926195701829,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.22723435237109,1.63836121154500,0.365830632201207) q[6];
u3(2.23992003314988,5.67334334936699,0.148993721129801) q[4];
u3(1.87099032660426,1.44886960659575,-2.51883559962189) q[9];
u3(2.21302215953232,2.23767699771682,-3.84450328682849) q[8];
cx q[8],q[9];
u1(-0.576982986678239) q[9];
u3(0.243769888950665,0.0,0.0) q[8];
cx q[9],q[8];
u3(4.15756920021335,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.40855081135004,1.67202722190065,-2.51298150805232) q[9];
u3(1.38787473106176,5.55006670134277,-0.643836574060084) q[8];
u3(0.641458979501496,0.974249795597295,0.0397020319492353) q[2];
u3(0.656197833363524,0.0820676089743984,-1.66250061043925) q[5];
cx q[5],q[2];
u1(3.36058557844988) q[2];
u3(-1.34281580715835,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.40945101030376,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.28373509866034,-4.46916699308722,1.78072583955624) q[2];
u3(1.38055838375281,0.529527383681243,1.29510955217399) q[5];
u3(1.92505814095207,-1.59269938143459,-0.950488938305546) q[11];
u3(1.95582284209776,-2.17211102524973,-0.215865336620219) q[3];
cx q[3],q[11];
u1(1.10054696540737) q[11];
u3(-1.41869223220771,0.0,0.0) q[3];
cx q[11],q[3];
u3(-0.345092195983273,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.22463415894391,1.29785065991742,-2.01067555207981) q[11];
u3(0.695182496026861,4.76595715945133,-0.654378803348538) q[3];
u3(1.78334109657188,0.736098103072527,2.03993617677239) q[2];
u3(1.63762079263663,-1.45535542373001,-1.57136907363888) q[4];
cx q[4],q[2];
u1(3.24337408453711) q[2];
u3(-2.14121533640040,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.46609818064292,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.80097992128810,2.56571981177008,-1.23778693862363) q[2];
u3(1.14234506357419,-1.02221794520098,-1.89595082433847) q[4];
u3(2.17964706938659,0.787031766244636,1.14382615371797) q[7];
u3(0.851363439847054,-2.85325117955494,-2.06462661577653) q[10];
cx q[10],q[7];
u1(1.47882397199407) q[7];
u3(-0.223118405622985,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.21345393513909,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.03744739875604,-3.90584780830102,2.02092093190899) q[7];
u3(1.65451592355545,2.11509601156030,3.47790312788553) q[10];
u3(1.49014798382592,2.95511171314287,-1.26063930372501) q[11];
u3(1.78665439107159,1.76997256647637,-0.633502859700559) q[5];
cx q[5],q[11];
u1(3.10428451275974) q[11];
u3(-1.97395970004195,0.0,0.0) q[5];
cx q[11],q[5];
u3(0.791963562360117,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.16195524064818,2.23288508716217,-2.17654788346974) q[11];
u3(0.456019763428385,-3.65500943614204,0.697509636915526) q[5];
u3(1.71228573119761,1.11632859565365,1.15893301346924) q[8];
u3(0.773478241784808,-2.09129499512333,-1.36058256584162) q[1];
cx q[1],q[8];
u1(3.04962209987717) q[8];
u3(-1.31484724794572,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.15709018774491,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.44076729796417,0.291611091994348,3.23526120288179) q[8];
u3(1.19911106951256,-1.05122422433139,3.96544282339129) q[1];
u3(1.97861532413098,-0.451143371692971,-2.32316733005044) q[9];
u3(2.35399143929518,5.58707072157865,0.234499029279975) q[3];
cx q[3],q[9];
u1(-0.403583016674564) q[9];
u3(1.39247994444609,0.0,0.0) q[3];
cx q[9],q[3];
u3(3.37280592694056,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.68240628969067,0.960817609222983,-1.01786356330336) q[9];
u3(1.80296953631411,3.56178002669814,-2.48806930647639) q[3];
u3(1.75703517912130,-0.425626325025430,1.00627723631933) q[0];
u3(1.76471630394629,-0.867987621557171,-1.34405577987679) q[6];
cx q[6],q[0];
u1(1.83757829934976) q[0];
u3(-3.14658198664468,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.816092136846917,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.21858338027948,3.03573123837556,-0.175433073921360) q[0];
u3(0.270929409376451,0.680502567800521,1.27496440328801) q[6];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
