OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(0.901415636392812,1.12055117267830,-0.158903403744502) q[3];
u3(0.294969408362141,-0.141710534014527,-1.82544012859363) q[5];
cx q[5],q[3];
u1(1.42272207917865) q[3];
u3(-2.98743502354972,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.913715009311873,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.33384085109438,0.522948946294170,0.0256269107723112) q[3];
u3(1.73347700804325,0.0509757487483649,4.91568572767777) q[5];
u3(1.30707208039066,-1.45751904847353,0.482649573082133) q[1];
u3(0.993285842706849,-2.73944098604215,0.940936440450390) q[4];
cx q[4],q[1];
u1(0.662992551026638) q[1];
u3(-1.31350747645429,0.0,0.0) q[4];
cx q[1],q[4];
u3(-0.168388952493518,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.57698724008621,1.57444235534605,-0.656167320450391) q[1];
u3(1.48906521278119,-2.81723979730290,-2.02879079660957) q[4];
u3(1.53463326205091,3.28864775465776,-2.16219191403079) q[0];
u3(0.399711982669026,2.59124009667896,-1.97143757848753) q[2];
cx q[2],q[0];
u1(2.29761878862086) q[0];
u3(-2.01887016914164,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.173082880805981,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.20742230097281,-2.26568831566812,-1.86031481918661) q[0];
u3(2.57067294444774,1.48480613273413,-2.34308198404872) q[2];
u3(2.01700112557825,1.99990858646435,0.383938266511987) q[5];
u3(2.18728537636605,0.249295732111345,-2.33866645353985) q[3];
cx q[3],q[5];
u1(2.72793488906897) q[5];
u3(-2.05264187130507,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.690471185288657,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.08970122862190,0.373254926618365,0.842049566420064) q[5];
u3(2.41890049620886,1.79032164509511,0.162127794174958) q[3];
u3(0.913234411481686,0.321022620561494,-0.495547831844371) q[2];
u3(0.816525602063009,-3.37199220524674,0.411372806146367) q[0];
cx q[0],q[2];
u1(-0.343631771736728) q[2];
u3(-2.35160170016043,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.59451901889354,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.02118799276221,2.08471447168570,-3.08199474418188) q[2];
u3(2.35715074520972,-1.03509525671171,-3.81314442345799) q[0];
u3(2.19321264268706,-0.678139218954325,0.215345135951419) q[4];
u3(1.29599884306125,-2.27826423535878,-1.15577402171906) q[1];
cx q[1],q[4];
u1(1.87501958230982) q[4];
u3(-2.41979509503884,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.322075982963715,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.28343903461748,-0.743481572471275,-0.787323967848615) q[4];
u3(2.46937317266405,4.77617318497998,-0.185724203529517) q[1];
u3(0.352502176411873,0.556199222131969,-1.29820021835314) q[1];
u3(1.69673338538564,-4.46098564566856,1.63003018696402) q[2];
cx q[2],q[1];
u1(1.72511140259418) q[1];
u3(-2.25764980123497,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.28334698736423,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.45427610751199,-1.97230044970905,2.52460666163013) q[1];
u3(1.35963856811287,3.47679679751070,-0.177037535697790) q[2];
u3(2.55038048422649,2.24832136913034,0.556116286437514) q[3];
u3(1.80160477959726,-0.0566161019021051,-2.01290294950065) q[5];
cx q[5],q[3];
u1(1.34461371863948) q[3];
u3(-0.968688675041822,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.139238928634067,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.11678161733208,1.61912044219457,-4.27137086843504) q[3];
u3(2.03773738326368,4.69274410264853,1.40745932358725) q[5];
u3(2.26756604190162,1.84072271573106,-0.419715433402625) q[0];
u3(2.57629170729138,-0.137494468544379,-4.23870271668198) q[4];
cx q[4],q[0];
u1(2.68275649464751) q[0];
u3(-1.85615647855507,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.836546821493917,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.67329759804964,1.42083844411480,-2.33492678231010) q[0];
u3(2.03223957881539,-1.45180256390949,1.09172547928241) q[4];
u3(2.60897202607679,0.995864379906027,1.35322983554820) q[0];
u3(0.931019289547601,-3.14434571296962,-1.37203635239463) q[1];
cx q[1],q[0];
u1(4.43096344834337) q[0];
u3(-3.00478290041194,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.401027405211100,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.95448459431000,-3.97832748677800,0.712437099954304) q[0];
u3(1.69517226564186,1.43986625248820,-2.15090115822470) q[1];
u3(1.65722726044401,2.35672442741149,-2.11027593904369) q[4];
u3(0.758269318453775,1.71431763915033,-1.54192001457999) q[5];
cx q[5],q[4];
u1(-1.16908040520073) q[4];
u3(0.454595987744503,0.0,0.0) q[5];
cx q[4],q[5];
u3(3.70976881190508,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.23493723942613,3.64027751721565,-1.45301146702187) q[4];
u3(1.26699941211419,-4.50484036069813,1.65039934195681) q[5];
u3(1.67461001901520,-3.93623084007599,0.868417465997815) q[3];
u3(1.64731079180981,0.101934748958347,3.71249354915915) q[2];
cx q[2],q[3];
u1(2.45462220159715) q[3];
u3(-2.86944417020167,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.03651364752361,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.89820011145842,1.60624076553225,-3.45407846998495) q[3];
u3(1.77475347980225,2.01835859212370,-2.44212530127653) q[2];
u3(2.07231346766391,-1.03098567204411,-0.809204254049321) q[1];
u3(0.885262583950905,1.26175335935236,-4.74578730474081) q[2];
cx q[2],q[1];
u1(1.74863829169545) q[1];
u3(0.172531698204515,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.12016687753786,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.05518576875807,1.14007254886188,-0.233852721134961) q[1];
u3(0.729354293144833,0.511616200992822,4.83667236357463) q[2];
u3(0.751066116898310,0.779710574893294,2.21346533364678) q[5];
u3(1.78979369787813,-1.26941217357063,-1.85152700740407) q[3];
cx q[3],q[5];
u1(1.48087258907311) q[5];
u3(-0.165131361272329,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.78882108019996,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.60751579226151,0.956124349555149,-1.62566761786486) q[5];
u3(2.27922549073732,-0.415920576775661,-3.44476281698177) q[3];
u3(2.51675012891228,1.47767675461728,1.24704894696910) q[4];
u3(1.89310528490275,-0.208660308891052,-3.47671201125153) q[0];
cx q[0],q[4];
u1(0.700692339309138) q[4];
u3(-3.43798226487123,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.72089717014558,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.71398263566363,-0.630377395123101,0.369788134012147) q[4];
u3(2.37308660786772,0.933588838335710,-1.15738333080320) q[0];
u3(1.91348376520615,0.418360559134948,2.13374284032827) q[5];
u3(0.715453683811127,-2.64327275717560,-2.24847270234990) q[2];
cx q[2],q[5];
u1(-1.11992117490548) q[5];
u3(0.173418630614084,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.41489768921809,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.85056421165157,-2.07506292895144,1.67860161495039) q[5];
u3(0.775856418460797,-0.886421874371378,-4.38948598363481) q[2];
u3(0.668290036987845,-1.41883208731906,0.350833755205659) q[1];
u3(0.768554251315631,-2.11717048907044,1.14428295211671) q[0];
cx q[0],q[1];
u1(-0.0450898342287194) q[1];
u3(-1.39474451948220,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.28763067043917,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.76052398763627,2.61078083629265,-2.55747830303007) q[1];
u3(1.02179653038291,3.08736409770099,1.59496944334905) q[0];
u3(1.61178754317958,2.50566476765979,-3.44393043370336) q[4];
u3(2.54393273971451,-2.96102160853633,3.09869054427948) q[3];
cx q[3],q[4];
u1(2.46383958715997) q[4];
u3(-1.76456925971809,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.0152484684571901,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.925124275922081,1.31875930095395,-2.92443013807657) q[4];
u3(2.88316571753182,-5.01815207453487,-0.199803219611230) q[3];
u3(1.47374998813494,-0.353388620958206,0.710047014283772) q[4];
u3(1.42605121361176,-1.71480536111371,-1.24832000040894) q[1];
cx q[1],q[4];
u1(-0.539913035230233) q[4];
u3(-1.66564339825229,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.19334018159951,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.07463595120280,3.73026215614824,-0.790730936211318) q[4];
u3(1.02894473210824,0.717111399216273,3.00790730722114) q[1];
u3(1.51043326868267,0.895903474920070,-3.07497012222880) q[0];
u3(0.995130206320897,-2.99669762800707,1.99770966442343) q[2];
cx q[2],q[0];
u1(1.53119346454718) q[0];
u3(-2.38939800124215,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.37302850395556,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.29231594969655,-0.846858706405442,4.10107805485028) q[0];
u3(1.33548761408415,0.795501994830984,0.599765325713415) q[2];
u3(1.60715730487242,1.14321288932145,-3.84109635526181) q[5];
u3(1.26488647384131,2.20088945759925,-2.35663650872152) q[3];
cx q[3],q[5];
u1(3.19084852668199) q[5];
u3(-2.73551290587664,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.15843049168343,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.69281046595104,-1.99442405702063,0.820721430632692) q[5];
u3(1.15988857779832,-3.89197309879939,1.55758320543199) q[3];
u3(2.53163139310229,0.121959525459499,-0.464634427319670) q[1];
u3(1.28075524102367,-0.345605745569232,-4.25288539608361) q[5];
cx q[5],q[1];
u1(0.00242776982953608) q[1];
u3(-1.87301416420460,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.600673400478958,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.58334237457737,-1.25278272563973,-1.71106387495325) q[1];
u3(0.928626434083619,-1.58847686373633,2.59954470653170) q[5];
u3(1.82799250081441,0.148236295051120,1.36623937784459) q[4];
u3(1.54270380301109,-0.459478631332779,-0.412661962436863) q[0];
cx q[0],q[4];
u1(3.09521029025014) q[4];
u3(-1.92458043798695,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.857806486222595,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.41237289820416,0.959409566749345,2.36266475905562) q[4];
u3(0.908683691477934,3.38118359775559,0.229796168531590) q[0];
u3(2.41877901731119,0.148471162713933,-2.60399907452365) q[3];
u3(2.53152637967640,4.47485538162491,0.808335237467508) q[2];
cx q[2],q[3];
u1(2.52691709488781) q[3];
u3(-1.66528428138548,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.60156345448330,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.331595336016710,-0.585287907102498,1.94392920374493) q[3];
u3(2.23348676855325,3.85178843008694,-1.13303563439872) q[2];
u3(1.57922255050645,-1.49146632729531,-0.625939835255984) q[4];
u3(1.23847933715255,-1.98564256333868,-0.952725727563349) q[3];
cx q[3],q[4];
u1(-1.05427917506739) q[4];
u3(0.657526776147455,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.97635423513788,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.87375672647375,2.96744638208197,-1.90182207462910) q[4];
u3(1.95270926686329,-1.92068732447248,3.86934739900347) q[3];
u3(0.870268949744039,0.609025921024491,0.320831240589713) q[0];
u3(1.94659772831508,-0.520715350117135,-2.92079324045904) q[2];
cx q[2],q[0];
u1(3.03861568124862) q[0];
u3(-2.28133774927706,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.707460109809297,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.906988849892965,-0.672355903083219,0.305324923501893) q[0];
u3(1.88265129481251,-0.979522082758014,1.70701192539091) q[2];
u3(2.88247266084877,0.627422044857671,-0.0467678385325505) q[1];
u3(1.40631982009269,-3.75667527823222,-1.06856582555384) q[5];
cx q[5],q[1];
u1(1.17287050677634) q[1];
u3(-0.671354125087371,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.65425949424381,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.23156823736277,-0.478119983340053,0.165701174147465) q[1];
u3(1.94413937286708,2.07340885184908,1.75628008942360) q[5];
u3(2.36463711875741,-1.12567995873727,-1.72220623252724) q[3];
u3(0.854075693638302,-1.25992645007748,-3.68168837681798) q[4];
cx q[4],q[3];
u1(-0.0119903388925151) q[3];
u3(-1.26999142445869,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.18771657168523,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.77278047051445,-1.95916278359188,-0.256299468645506) q[3];
u3(0.705752681038878,-0.364888285114730,2.31916594245850) q[4];
u3(2.41952039513514,-1.68344407968515,0.430648273713820) q[5];
u3(2.82278805752561,1.88686370242131,3.46625083596525) q[0];
cx q[0],q[5];
u1(1.98970989380740) q[5];
u3(0.110329599233896,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.33045143323853,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.02366933054657,-1.33642150995680,1.87254940382724) q[5];
u3(1.42836697908526,-2.65291017872992,-0.723299779414441) q[0];
u3(1.11831938530487,1.13997879450566,1.01587955978512) q[2];
u3(1.78201565015077,-1.51904541743691,-0.278296761889585) q[1];
cx q[1],q[2];
u1(-0.174797919437024) q[2];
u3(-2.47209449453377,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.69492625903932,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.35693729175908,2.85750107676078,-1.56336138499846) q[2];
u3(1.81194593989387,-1.44024990153262,1.04612762875436) q[1];
u3(2.09444774553933,1.61035063719569,-0.743415242608827) q[3];
u3(1.86750475125971,-4.50534632174548,1.35476259001433) q[0];
cx q[0],q[3];
u1(1.59797026323588) q[3];
u3(-3.17400839505360,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.78922376920901,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.77363102189701,-3.45401032843888,2.43634952087486) q[3];
u3(2.64093150846420,2.79795281936329,0.784422361635313) q[0];
u3(0.342299541351091,-1.53016890814095,2.47926811821601) q[1];
u3(0.778131505267908,-2.85157115033271,1.01914085500517) q[4];
cx q[4],q[1];
u1(1.54237899196108) q[1];
u3(-0.888695059230537,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.82490030034170,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.61625879569842,-1.11928208762382,2.79331600097041) q[1];
u3(2.39304603140929,-0.322521089183376,-1.95135348454612) q[4];
u3(2.21965930542037,2.32781127961763,-0.343989270459931) q[2];
u3(2.97445743415578,4.14875688140075,0.103139194857010) q[5];
cx q[5],q[2];
u1(2.22161111942957) q[2];
u3(-1.86580210289549,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.27318861082052,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.06508067164413,0.0747423421568115,0.564465893838335) q[2];
u3(1.14412352365917,-2.71211593754490,-1.60144197139962) q[5];
u3(0.231838046784911,-2.04240427670183,0.101819836605436) q[3];
u3(2.02797078054159,-5.07262135331797,0.145415099749166) q[4];
cx q[4],q[3];
u1(-0.0398958105223166) q[3];
u3(-1.78545954583383,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.07943057697102,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.475143015594951,-1.65809847105704,3.67280123465195) q[3];
u3(0.149556693802635,1.05871988751780,-2.54927895348412) q[4];
u3(1.40583103700399,0.813589870181532,-3.29851580159231) q[2];
u3(1.57941378555734,3.42981905171224,-2.35758004922837) q[1];
cx q[1],q[2];
u1(1.75863329137319) q[2];
u3(-0.712321847917196,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.124866850422489,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.724834049847318,-0.208512136677144,-1.59597974594933) q[2];
u3(2.47655134611870,2.02358420257200,-1.33795553072224) q[1];
u3(1.54572701992754,-1.63477190952152,-0.116530685989287) q[5];
u3(1.82027252003671,-2.15118534810992,-0.0635851306963879) q[0];
cx q[0],q[5];
u1(0.814396612365918) q[5];
u3(-0.230747978901800,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.42527400208530,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.51400657178717,-0.361464867645227,0.768551927596336) q[5];
u3(0.938092626688565,-5.22187570515153,-0.393970134056624) q[0];
u3(1.81392651903525,0.745532681068154,-3.84910448105468) q[2];
u3(0.929489299885559,-2.34054976324454,3.21236080538639) q[4];
cx q[4],q[2];
u1(1.35296908562192) q[2];
u3(-0.582442609278761,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.314034228680861,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.56806633259732,-1.66995122304689,1.38480670166057) q[2];
u3(1.49932791486754,-3.44913796802201,1.54164913137991) q[4];
u3(1.76478024131869,1.59881577518939,-2.64936452293327) q[0];
u3(1.86317181249719,2.01393762499062,-3.27371387208659) q[5];
cx q[5],q[0];
u1(1.24251275376242) q[0];
u3(-3.32031072242131,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.78367289091054,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.07277689531660,0.588403938695056,-0.380887460999528) q[0];
u3(0.620153549466140,1.32763352809548,-2.61237786317294) q[5];
u3(1.72212668202402,2.65138472678334,-2.45514129431871) q[1];
u3(1.76867417765860,-3.11137489462284,2.45350762809369) q[3];
cx q[3],q[1];
u1(3.01184061903806) q[1];
u3(-1.78272391661217,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.45458427997560,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.72258153503435,-1.78138051007293,0.248602513403601) q[1];
u3(1.02724717114578,-3.09937162558817,-0.659923677853288) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
