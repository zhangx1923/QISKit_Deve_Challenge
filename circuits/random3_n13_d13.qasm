OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.31678588855707,2.60966340793212,-2.06026698048557) q[5];
u3(0.960024047268578,2.00434254456352,-1.98730909893446) q[7];
cx q[7],q[5];
u1(3.07008507172560) q[5];
u3(-2.16168411393125,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.48809422951790,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.964436393635551,1.77738111204176,1.96612189401100) q[5];
u3(0.218156173912150,1.78771405837993,3.08979502474053) q[7];
u3(1.72082958950686,-1.28527921769562,-1.61739731263740) q[8];
u3(1.79029178278520,1.90596517019109,-3.80601555699718) q[1];
cx q[1],q[8];
u1(1.77571787310139) q[8];
u3(0.133570140589876,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.951972750413856,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.564371591658056,-1.37726836634897,1.79923684479105) q[8];
u3(2.91624580073889,-4.84255887614556,0.554305081854496) q[1];
u3(1.31736933188887,0.0483330102891536,-1.78318755989918) q[10];
u3(0.524495839065813,-3.53304136066506,1.41633100110707) q[11];
cx q[11],q[10];
u1(0.660945230487182) q[10];
u3(-1.49478908332998,0.0,0.0) q[11];
cx q[10],q[11];
u3(-0.301330096077768,0.0,0.0) q[11];
cx q[11],q[10];
u3(2.90607259336173,1.61420482444793,0.351332419737800) q[10];
u3(2.07432521112373,1.57399733968757,-4.62867624027431) q[11];
u3(0.204288691722919,3.08264695565545,-2.88834782796440) q[0];
u3(0.660663349614256,0.606914066115770,-1.72326172478861) q[4];
cx q[4],q[0];
u1(0.0651267197315049) q[0];
u3(-0.895162379320519,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.58756641989506,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.48279552009831,2.60859842628916,-0.590440983832702) q[0];
u3(2.68130992319307,5.15480150511299,0.993129845006000) q[4];
u3(2.93455938653105,-0.937660460245566,-0.574044248090056) q[6];
u3(0.968725482371825,-4.88241455970841,0.173534968587092) q[12];
cx q[12],q[6];
u1(1.09811575053766) q[6];
u3(-0.0340767738929064,0.0,0.0) q[12];
cx q[6],q[12];
u3(2.73268024907817,0.0,0.0) q[12];
cx q[12],q[6];
u3(2.98949690591316,1.66938588966465,-2.77559659724913) q[6];
u3(1.09069857210314,-1.22468277609806,1.44778620470926) q[12];
u3(1.70590836134224,-1.85646077791761,2.59404067736352) q[9];
u3(0.779516630218408,-0.918829278440171,0.740544202214164) q[2];
cx q[2],q[9];
u1(1.77228768906631) q[9];
u3(0.810412658777456,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.10812097031550,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.849172494166034,3.50577820301365,-2.75453725990143) q[9];
u3(1.42018054413768,1.78714270322649,1.69095415768579) q[2];
u3(1.00628798159553,2.72847418327722,-2.10660767359832) q[12];
u3(0.678046064340774,1.18523176411791,-1.57752467990699) q[6];
cx q[6],q[12];
u1(-0.444401704381354) q[12];
u3(1.19393797225272,0.0,0.0) q[6];
cx q[12],q[6];
u3(3.85086304800741,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.75026495169063,0.0170165363483270,2.88325363348645) q[12];
u3(1.46994524548815,-1.52858903094007,3.82784483840032) q[6];
u3(1.69923502604966,1.35830429405176,-2.69828530850186) q[5];
u3(1.75265991521222,-2.06571787074391,2.76056076457312) q[3];
cx q[3],q[5];
u1(0.409116341502282) q[5];
u3(-0.232901428029461,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.66711417593971,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.04864683418453,2.48875864268854,-2.99902970025400) q[5];
u3(2.64167156110090,0.920129364937576,3.11352809946855) q[3];
u3(1.10670952298988,0.0929411049629211,-0.877110170885924) q[9];
u3(1.52769238352596,0.956828627610884,-5.15374150047662) q[2];
cx q[2],q[9];
u1(1.97941512279847) q[9];
u3(-2.99139977891695,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.45632754181428,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.54842808310178,-1.00660807073392,-2.11590422738099) q[9];
u3(2.34550777086238,-0.567781401933555,5.47766507779474) q[2];
u3(2.86279119651992,2.74211458406586,-2.85025629940548) q[1];
u3(1.20228979710106,0.939696187947038,-0.0481493843183189) q[7];
cx q[7],q[1];
u1(2.28164370933548) q[1];
u3(0.441474319028762,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.75131424686977,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.15434612480683,-0.0553767858122529,0.305707351753280) q[1];
u3(0.957371766799914,-1.42336989388836,-2.38299120848963) q[7];
u3(1.34073558595633,0.580507385265984,1.79906235300994) q[8];
u3(1.66497813665260,-1.14502274863085,-0.500499259561741) q[11];
cx q[11],q[8];
u1(2.88437911286858) q[8];
u3(-1.57684095277196,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.11473950386703,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.69924220136466,-1.38602034860377,0.751144020498965) q[8];
u3(1.69535399475702,-0.674414226215656,1.30829858559570) q[11];
u3(0.780969924893787,1.86385052528291,-1.89964537861455) q[10];
u3(0.497832622330860,-2.55251782170232,2.07942600953576) q[4];
cx q[4],q[10];
u1(1.57427597444313) q[10];
u3(-2.40880063292564,0.0,0.0) q[4];
cx q[10],q[4];
u3(3.30782356437633,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.379067016432860,3.27868200461032,0.262951929757218) q[10];
u3(1.45257249484120,-0.0772121988027688,-5.70639025773475) q[4];
u3(1.56717481769477,3.39035118994764,-1.99934503562625) q[5];
u3(0.579245658899307,-0.0574412743680952,1.11171390615076) q[9];
cx q[9],q[5];
u1(1.48981866501660) q[5];
u3(-0.696779815952874,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.84759163179934,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.279771198106139,-0.888429555019975,-0.290053661448350) q[5];
u3(2.55128841683757,-4.12901159182829,0.342850204998103) q[9];
u3(2.60805624820198,-0.935093020215255,1.97617728272816) q[12];
u3(1.91295119514864,-2.15149384544321,-0.710397799383592) q[10];
cx q[10],q[12];
u1(1.20724409891308) q[12];
u3(-0.703174529451349,0.0,0.0) q[10];
cx q[12],q[10];
u3(-0.103996881400843,0.0,0.0) q[10];
cx q[10],q[12];
u3(1.20146384762737,-1.69411551761548,1.12812491079837) q[12];
u3(0.789883415330012,1.56431382189225,1.98711626359277) q[10];
u3(1.45082939624847,2.28774134819739,-2.59022291634373) q[8];
u3(2.02671435215586,-3.43832392923250,2.69216972868684) q[1];
cx q[1],q[8];
u1(0.122835903756751) q[8];
u3(-2.45529776852577,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.59994472395856,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.39288211622820,3.10324396168488,-1.63308345005511) q[8];
u3(1.73578031035660,0.0541178096553582,-0.529469046262352) q[1];
u3(2.45024782110917,2.55142124659372,-3.40822379073661) q[6];
u3(0.837089607080788,3.11755545059248,-1.58148513721334) q[2];
cx q[2],q[6];
u1(0.395418076445679) q[6];
u3(-0.872275308415955,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.92782464281055,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.73776231922498,-1.78547490684726,2.31918606351303) q[6];
u3(2.63057020302154,-1.39691617560783,3.60364829611624) q[2];
u3(1.57166629008264,0.282817965501775,0.899975306883188) q[11];
u3(0.506914129684670,-2.65869102136605,-1.87773120920724) q[7];
cx q[7],q[11];
u1(0.518055409203342) q[11];
u3(-1.23589406805736,0.0,0.0) q[7];
cx q[11],q[7];
u3(2.93768247033547,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.52284781251379,-3.43462143378458,-0.541501287689700) q[11];
u3(1.73388869638070,-3.16300636557457,1.37532234155008) q[7];
u3(1.08467115568674,2.41439882552131,0.353498932859033) q[0];
u3(1.55621074851196,0.664917912204167,-3.60364089011017) q[3];
cx q[3],q[0];
u1(1.38738663251830) q[0];
u3(-0.273580768638456,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.72594141428500,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.31744247633893,-2.21658463554439,2.64499070357986) q[0];
u3(1.22325373971035,-0.936258809277682,4.75774489959739) q[3];
u3(1.87046094465137,-2.53429842830186,-0.226760922740733) q[9];
u3(1.74196347493847,-3.52851252419402,-0.515210701188680) q[0];
cx q[0],q[9];
u1(1.05801704343338) q[9];
u3(-3.10977896356931,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.60639270912048,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.42874898948769,3.81774762101714,-1.12973736219812) q[9];
u3(1.14593328068647,-1.58741698820105,2.14824240094566) q[0];
u3(2.62453975878943,2.04664184687217,-2.99713566365133) q[12];
u3(1.56617717942067,2.49711969358064,-3.13514039563628) q[7];
cx q[7],q[12];
u1(0.0342748558440260) q[12];
u3(-1.35279255441439,0.0,0.0) q[7];
cx q[12],q[7];
u3(2.01948493329280,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.532645769475730,-2.34854673652196,1.87391396438962) q[12];
u3(1.48381920593806,-0.501745324880908,2.93294195391352) q[7];
u3(1.07009469834989,0.952508801904625,-1.55850787886121) q[1];
u3(1.11162116389362,-0.735520983671397,-0.0574094776147352) q[2];
cx q[2],q[1];
u1(0.0235088792709495) q[1];
u3(-1.08528234416448,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.66515064967947,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.39118260001426,-4.08206425660439,1.07414653114286) q[1];
u3(1.53325571228986,0.0626998163041224,4.49997186327277) q[2];
u3(1.36803400828518,1.73535861538814,-1.46061891014763) q[3];
u3(0.877941900323334,0.424206923356697,-3.29240029665354) q[4];
cx q[4],q[3];
u1(2.18229411761162) q[3];
u3(-1.66511896232071,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.31172009762146,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.97299180366763,-2.27381096899969,2.78418503478489) q[3];
u3(2.65569215604501,-4.78689115657475,0.453969619219324) q[4];
u3(1.56152494774169,1.66197470833833,-0.203615569693269) q[6];
u3(1.08387636392056,1.21841756715037,-4.46381057869502) q[11];
cx q[11],q[6];
u1(2.57433894060658) q[6];
u3(-2.22082921520398,0.0,0.0) q[11];
cx q[6],q[11];
u3(0.364332794320774,0.0,0.0) q[11];
cx q[11],q[6];
u3(2.44846741167324,1.57810981263745,-3.46221948711036) q[6];
u3(2.37911575493103,1.34526723654305,-2.07551285328647) q[11];
u3(2.36056159713508,1.28050121527231,-4.01336196919956) q[5];
u3(0.838365992283840,-2.22560422382816,3.56961708378860) q[8];
cx q[8],q[5];
u1(1.63503175995423) q[5];
u3(-2.25781616531608,0.0,0.0) q[8];
cx q[5],q[8];
u3(3.24097428824998,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.995811054936512,0.160352765966428,-0.617336127379124) q[5];
u3(2.48879671110906,-0.401727348392938,-1.59362453026788) q[8];
u3(1.47392341103277,1.39041212147099,1.43411170027507) q[10];
u3(0.426270487405593,-0.412543312027724,-2.92088977124270) q[3];
cx q[3],q[10];
u1(3.23443907155663) q[10];
u3(-1.34735921620546,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.77694988430061,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.78297751752438,-0.153022834027846,-0.928574681954516) q[10];
u3(2.45890017251534,0.912778256885270,-3.73260847136630) q[3];
u3(1.32071800161213,0.953706708554113,-2.43359633738381) q[1];
u3(0.151779791702034,2.58760595867779,-3.61817231538085) q[6];
cx q[6],q[1];
u1(4.40181426219683) q[1];
u3(-3.24764184741113,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.395980030490025,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.78417363378388,-0.143478767732809,2.51623182999750) q[1];
u3(0.635817220111726,-2.87086259352484,-1.51083046871182) q[6];
u3(1.41405352170961,-0.148016574507499,1.36628544952484) q[12];
u3(1.80315196361953,-0.242822047453213,-1.58305380487074) q[2];
cx q[2],q[12];
u1(3.26692637788743) q[12];
u3(-1.24219497328997,0.0,0.0) q[2];
cx q[12],q[2];
u3(2.59623888723574,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.60269083283545,1.20686280538361,-1.76490908553843) q[12];
u3(2.04258903154611,-0.382394564104369,-3.60918987389983) q[2];
u3(1.12927057692195,-2.21024910740024,0.914449531533676) q[5];
u3(1.22167958919661,-2.38839152375964,0.704140759140868) q[8];
cx q[8],q[5];
u1(1.67729905217857) q[5];
u3(-2.68680235426540,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.694978683203428,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.662664299411430,-2.04849752717237,4.08605382026116) q[5];
u3(2.58690147998746,4.85843636134305,-0.0882307315455853) q[8];
u3(0.971190179714532,-0.673911276064920,1.04886646708679) q[4];
u3(0.505118282609572,-1.55735251911462,0.677881322781940) q[7];
cx q[7],q[4];
u1(0.0949365732785448) q[4];
u3(-0.523356920411143,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.70154129965340,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.08403927580296,-3.36916106444386,1.95158692014923) q[4];
u3(1.85518575414937,-6.11796997902400,-0.134670541760335) q[7];
u3(0.660821686035048,1.37434085033581,-0.748000458166012) q[11];
u3(0.702626930263699,-1.98066943850636,1.07674107256940) q[9];
cx q[9],q[11];
u1(1.88452918537937) q[11];
u3(-0.269713321484368,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.524114639513725,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.17139335994932,-3.44475606483548,1.40517241696028) q[11];
u3(0.237541303631066,-2.11277162353465,0.131538827368672) q[9];
u3(1.88484737297289,-0.504529234455005,2.49594090294957) q[11];
u3(2.29833134998452,-1.85007280150628,-1.74446337713258) q[5];
cx q[5],q[11];
u1(1.35599282397870) q[11];
u3(-3.04331287056856,0.0,0.0) q[5];
cx q[11],q[5];
u3(2.43531230882623,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.32592826435559,-0.554722271620500,-0.613251278972272) q[11];
u3(1.98011553348762,3.93974481433666,-2.11443957034760) q[5];
u3(1.38727646518509,1.10689097123980,0.746564003690012) q[4];
u3(1.31048927309365,-0.475955661110160,-2.60134040813153) q[7];
cx q[7],q[4];
u1(2.33485193221083) q[4];
u3(-2.70280645403009,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.96264039112964,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.857572558193542,-2.84209803040723,0.718721820186481) q[4];
u3(0.437371678582277,-4.39653907483189,-1.02181655343417) q[7];
u3(1.84590331308034,1.53448970000004,-2.39787141652516) q[12];
u3(0.825423697986460,-2.29963447199968,1.98382753057952) q[1];
cx q[1],q[12];
u1(3.33859970552334) q[12];
u3(-1.29313341022065,0.0,0.0) q[1];
cx q[12],q[1];
u3(2.22069961682320,0.0,0.0) q[1];
cx q[1],q[12];
u3(1.94212702462861,1.25692076236946,1.86877903628845) q[12];
u3(1.78694788916478,-1.80670555452461,0.364305774356162) q[1];
u3(0.105114335644915,1.95977538398989,-0.765935792111830) q[2];
u3(1.14245659628269,0.567510621478322,-2.11136447533499) q[6];
cx q[6],q[2];
u1(3.17649278002662) q[2];
u3(-1.71209776691174,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.667213086047658,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.163046049654846,3.04522209882160,0.382319202834094) q[2];
u3(2.01236863941626,1.28396887649823,3.93857451127761) q[6];
u3(2.59248250604998,0.704748507249285,-3.50843296277039) q[10];
u3(2.55738286280325,4.71813685775597,0.246964306681664) q[3];
cx q[3],q[10];
u1(2.08626963256713) q[10];
u3(-1.82459615170991,0.0,0.0) q[3];
cx q[10],q[3];
u3(3.68720905747233,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.14881111582322,-3.15786580628531,0.162718581898468) q[10];
u3(0.643288090406107,2.96997083539972,-0.265262720513988) q[3];
u3(0.460819228961298,-2.73853159248439,3.16018298339616) q[8];
u3(0.320699137936393,-3.37765139576564,0.629377664240387) q[0];
cx q[0],q[8];
u1(3.19138736760127) q[8];
u3(-1.85980011964724,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.06718271541526,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.57980201633100,-3.46598401374347,-0.889277113881753) q[8];
u3(2.45700038813657,2.52469744322795,-1.22060952705932) q[0];
u3(1.86768214689663,1.62970446457152,-2.32786167437308) q[3];
u3(1.44020677612410,2.47421515415314,-3.23616021607021) q[11];
cx q[11],q[3];
u1(0.792807365130151) q[3];
u3(-1.63978469040310,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.91377708117840,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.60242992541158,-2.46859363085954,1.17230058972020) q[3];
u3(1.28215376450890,3.28074966807415,0.343437598556324) q[11];
u3(0.539812411635017,1.84318614108834,-2.23523216815779) q[4];
u3(1.58807900775859,1.87370979615568,-4.30799977720317) q[0];
cx q[0],q[4];
u1(0.232607339261392) q[4];
u3(-1.69861027474247,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.10899907127061,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.23300329431091,-0.406006834984411,0.566719439510832) q[4];
u3(2.24591889369327,-2.71082991244255,0.174407319057749) q[0];
u3(2.38861816822604,1.13246302602094,0.732222712896124) q[12];
u3(0.614460834251022,-3.90474243252937,-0.987041001441368) q[1];
cx q[1],q[12];
u1(1.73786714804908) q[12];
u3(-2.67409085694959,0.0,0.0) q[1];
cx q[12],q[1];
u3(0.0506634973984661,0.0,0.0) q[1];
cx q[1],q[12];
u3(1.51088258654733,-2.35602602400337,0.527613591265502) q[12];
u3(1.45799975518302,-1.51552245336250,3.96373942389015) q[1];
u3(1.19367608655449,1.02352438317493,-3.14336511529449) q[6];
u3(0.508536118749013,-2.98887778733987,2.57807995950297) q[10];
cx q[10],q[6];
u1(0.387743164731614) q[6];
u3(-0.892847413021556,0.0,0.0) q[10];
cx q[6],q[10];
u3(3.19530657649209,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.38888713935894,-0.398803184938955,0.180433614690484) q[6];
u3(0.548710514766169,0.143719624099835,4.24151638675234) q[10];
u3(0.860903001420966,1.86263292745619,-2.01980462992658) q[8];
u3(0.0926382296034513,-1.44767846288580,-1.28138804168069) q[5];
cx q[5],q[8];
u1(2.63000097263955) q[8];
u3(-1.38544822274703,0.0,0.0) q[5];
cx q[8],q[5];
u3(3.16718357676481,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.40227187199836,-2.04916456170246,2.22964245926018) q[8];
u3(2.24217649805425,-4.20381745655257,0.305482689291531) q[5];
u3(2.63069625266800,1.92973011445066,-0.779016356063329) q[2];
u3(1.88948913790908,4.47597924878747,0.329169286595902) q[7];
cx q[7],q[2];
u1(-0.164624932863088) q[2];
u3(-2.43215597782098,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.44357282154554,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.65022579438657,1.35224120643100,-0.460304438601530) q[2];
u3(2.83757951080214,1.61074257132157,3.91777231363404) q[7];
u3(1.92234366213677,-1.96949516027113,-0.545786839997897) q[9];
u3(1.93298167165779,-3.90928048259758,0.827346290473147) q[1];
cx q[1],q[9];
u1(1.74149187469048) q[9];
u3(-3.23308290926910,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.776902223388614,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.70462346217199,2.15118885226587,-3.45758206606409) q[9];
u3(1.66204306582415,-3.14392798961925,1.06418101161120) q[1];
u3(0.735256499199692,-0.174987552117568,-1.92948992220538) q[2];
u3(2.00430163533654,0.176497819172982,-4.93269227139136) q[11];
cx q[11],q[2];
u1(0.0684888398777987) q[2];
u3(-0.933799815254157,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.40127900842079,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.16765871253454,-0.177022100324365,-3.67988769799139) q[2];
u3(0.563462211119690,0.239179791669179,0.322472192035553) q[11];
u3(2.19358780394699,-1.67051034476200,4.49952998287908) q[5];
u3(0.976389702700752,-0.969572717152261,2.56232474986483) q[7];
cx q[7],q[5];
u1(-0.153694218396381) q[5];
u3(-2.22915593397029,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.67549816199077,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.540207406188798,2.55386914249457,0.210169448664310) q[5];
u3(1.98575016507298,-3.98068566763973,-0.285996664498690) q[7];
u3(0.626662984061578,1.95296197682958,-1.03133866784053) q[12];
u3(1.93944102483781,0.786780577392717,-1.83779891383319) q[10];
cx q[10],q[12];
u1(1.32080591957617) q[12];
u3(-3.52794391163462,0.0,0.0) q[10];
cx q[12],q[10];
u3(2.34055675182257,0.0,0.0) q[10];
cx q[10],q[12];
u3(1.06928140063702,0.833265988322791,1.35144764963554) q[12];
u3(1.45207487220878,-4.65472883488701,0.976534873559270) q[10];
u3(2.83783615236595,-2.35274879668212,3.20438033406289) q[6];
u3(0.873040972766478,-0.966685122382346,2.43216922280909) q[4];
cx q[4],q[6];
u1(1.51722645416195) q[6];
u3(-3.22424237177440,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.58807579302786,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.92005153443948,1.36667421093467,1.83803160612047) q[6];
u3(0.908856248777663,-2.06143636431535,1.77541901898395) q[4];
u3(1.91041146786621,-0.671073759892475,-1.89870134144237) q[0];
u3(1.79389227628225,1.56088093364492,-4.26497569248215) q[8];
cx q[8],q[0];
u1(0.912849912872202) q[0];
u3(0.210798624378478,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.97963783316220,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.07642197387718,-3.91810201842311,1.57979674963723) q[0];
u3(2.57798365526401,2.21460008229544,-3.15953991374919) q[8];
u3(1.44237683243707,0.664376319204519,1.80987191295088) q[4];
u3(1.48611633647536,-2.05275532409007,-1.64936559524016) q[5];
cx q[5],q[4];
u1(1.39067191090696) q[4];
u3(-2.37445836079954,0.0,0.0) q[5];
cx q[4],q[5];
u3(3.14536460493687,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.785624882761532,0.299933778795317,3.06324677362468) q[4];
u3(1.50310744306300,-1.05675516002059,2.58478407774453) q[5];
u3(2.50964065737807,3.69205041411730,-1.06972703427940) q[7];
u3(1.60479067836598,1.76028414850148,-0.407662281222382) q[0];
cx q[0],q[7];
u1(0.210756007102706) q[7];
u3(-1.13916453091401,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.17171440050097,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.34521944168289,-0.481605499958555,2.27961873340534) q[7];
u3(2.02024099318289,-3.80374888193136,1.36944658536177) q[0];
u3(1.44924993608305,-0.0141489097906071,1.82907733899773) q[9];
u3(1.42999806176203,-1.45661914117810,-2.00271713361288) q[2];
cx q[2],q[9];
u1(1.53553787565996) q[9];
u3(-2.16307325186444,0.0,0.0) q[2];
cx q[9],q[2];
u3(3.08895753786717,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.08598715771152,-2.15706319418869,-0.837799340072501) q[9];
u3(3.06368624641902,2.49669519487318,0.672575282949064) q[2];
u3(0.686739562531833,2.99341754925211,-1.37614765753636) q[10];
u3(1.43864263739752,1.68144132445283,-1.43437826979686) q[11];
cx q[11],q[10];
u1(2.15421675613464) q[10];
u3(-2.68177101777460,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.38513067536446,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.29651919039911,-0.906799024505174,-1.42212503989330) q[10];
u3(2.61587711055010,1.54387272164598,-3.52786850057617) q[11];
u3(1.77428602902682,1.24483552387223,-3.09602329475647) q[12];
u3(1.42448847814402,2.85634951010099,-3.13409675811163) q[8];
cx q[8],q[12];
u1(2.75555523016130) q[12];
u3(-1.50343184868316,0.0,0.0) q[8];
cx q[12],q[8];
u3(0.561991218045549,0.0,0.0) q[8];
cx q[8],q[12];
u3(1.01875690753063,1.40096061742069,0.218070684722011) q[12];
u3(1.14947602386041,-0.416543248112649,2.30132019970646) q[8];
u3(0.733873039247207,1.40500698764010,-1.77815122820715) q[3];
u3(0.262648896540288,0.698552967582633,-3.09105293663192) q[6];
cx q[6],q[3];
u1(-0.467547762240462) q[3];
u3(-1.95077252676550,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.73618057176466,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.43591956731935,-3.94039823883638,0.420361894017821) q[3];
u3(1.34692261540890,-0.749405595404766,3.90953939772153) q[6];
u3(2.41107387530587,-0.487568498686119,2.56881245504379) q[2];
u3(1.96868667626371,-1.61219317436342,-1.17415964199289) q[3];
cx q[3],q[2];
u1(1.87812523992074) q[2];
u3(-2.38873655767221,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.06397018846374,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.06022747770827,1.22288965533053,-0.257263023299503) q[2];
u3(2.79147020065458,-2.48089753808443,0.986718162262634) q[3];
u3(1.36360209513447,-1.01621978175579,2.71502196974952) q[6];
u3(1.05889835472672,-1.81317794574256,-1.61073891085339) q[0];
cx q[0],q[6];
u1(2.32842597665391) q[6];
u3(0.453533306543452,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.85043736157271,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.83675344239366,2.87118015097052,-2.65861099752361) q[6];
u3(2.07008660395696,-1.27881868889673,-4.37878470714206) q[0];
u3(1.61371332051372,-1.06185643147366,0.137157730204675) q[7];
u3(1.55328525291845,-2.58081042859212,1.37563503311854) q[5];
cx q[5],q[7];
u1(2.40471798016529) q[7];
u3(-3.08700367326882,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.999669419038016,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.62896631574163,1.78256606303183,-3.73602258284801) q[7];
u3(2.31233302326868,-3.46859294052457,1.32093780488704) q[5];
u3(1.85304995999013,1.87029055027396,0.131742716430254) q[9];
u3(0.781623509896677,0.361895976933764,-4.06562669480825) q[10];
cx q[10],q[9];
u1(-0.482041509530746) q[9];
u3(-1.77704943245199,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.796254611757364,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.94980641135737,1.08527705590347,-2.06889700755949) q[9];
u3(0.624240682191067,-0.816031350970872,0.833559472322270) q[10];
u3(2.22626484497455,-0.191676409291019,2.05507248219692) q[1];
u3(2.38637389558233,-1.81648731797363,-1.81788648911237) q[4];
cx q[4],q[1];
u1(2.72171495782186) q[1];
u3(-1.65178888644099,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.526422586535981,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.38399058735212,3.42351670600433,-0.456649858081437) q[1];
u3(0.849075529376812,-0.494203588247646,-4.35695568000398) q[4];
u3(1.24265005808380,1.82523034842996,-3.65984019660936) q[12];
u3(0.977984919472547,-2.53475004784336,3.66618657749540) q[11];
cx q[11],q[12];
u1(-0.777744524413917) q[12];
u3(-1.61751693864645,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.00841609894505,0.0,0.0) q[11];
cx q[11],q[12];
u3(0.633551392615936,-3.07827786617038,3.04871009779080) q[12];
u3(1.05510642755317,1.55563935484286,1.20753422064567) q[11];
u3(1.46975828022973,-4.00948747600435,2.26468033540025) q[11];
u3(2.29364086933898,2.92341817024873,-3.10003511771818) q[10];
cx q[10],q[11];
u1(-0.239052632313434) q[11];
u3(-2.41631712212220,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.46196033907218,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.33666519083493,-2.52942059232992,0.559027141091910) q[11];
u3(1.14585218424717,3.54966924996488,1.44468440316579) q[10];
u3(2.08916111401780,-1.78189600810515,-0.518526253075150) q[5];
u3(1.59261056572462,-3.19726778540541,0.155339450939373) q[3];
cx q[3],q[5];
u1(4.08133662333436) q[5];
u3(-3.53008215472182,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.764810634619826,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.46481305929094,-2.07362280137677,-0.225616202433140) q[5];
u3(1.52627272161742,-1.61482106766740,4.29522457689859) q[3];
u3(0.275564659062665,-2.11762593960191,2.67862008079628) q[8];
u3(1.10091524155361,-2.71076309235025,0.649298275988804) q[1];
cx q[1],q[8];
u1(0.112768337427699) q[8];
u3(-1.19473128327141,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.03535378471111,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.822741577105143,-2.59507287717064,1.11664671184814) q[8];
u3(3.06766419900653,0.853057593967595,1.08689022665247) q[1];
u3(2.96593699544749,-3.00013264787215,0.0241401944799162) q[7];
u3(2.88796705865028,-2.42023354200832,0.402878420408513) q[4];
cx q[4],q[7];
u1(1.07590587858716) q[7];
u3(-0.416117192966427,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.79914216697748,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.12283568516183,1.80887296191514,-0.169391687605278) q[7];
u3(1.25581134403107,-3.94457651361880,-0.992355656452268) q[4];
u3(1.83677533814087,-0.394002009729299,0.413404288889654) q[2];
u3(2.27380514832335,-0.932106692916025,-1.27190908177490) q[9];
cx q[9],q[2];
u1(1.58296523024067) q[2];
u3(-2.18273314573968,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.0329422391986218,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.24450593763146,-0.798479237144183,-0.461008609209870) q[2];
u3(2.30181564291456,-3.59066649374826,1.82657343273287) q[9];
u3(2.68536046922759,0.421989945116679,-2.73129108494666) q[0];
u3(1.96703314547090,2.26316420389665,-2.55577624796771) q[6];
cx q[6],q[0];
u1(1.82584431489411) q[0];
u3(-3.05381432010112,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.00939614382665,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.839291311925443,1.88096109950922,-2.49827043070467) q[0];
u3(2.22498816152517,-0.812275480744719,-2.83480657458540) q[6];
u3(1.88229043340167,0.0509146904899924,1.03046727952102) q[5];
u3(1.85084588774779,-1.51508063764695,-1.22156907048330) q[10];
cx q[10],q[5];
u1(0.572647987881277) q[5];
u3(-1.53284887551139,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.49368646601839,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.38506038999512,-1.99431385454270,0.530729295475911) q[5];
u3(0.735068941602556,-4.73540809797688,-1.13070278413208) q[10];
u3(2.15411757315827,0.419159605291064,1.11195326575679) q[3];
u3(1.96099047630062,-1.71644010682943,-2.42088879220893) q[2];
cx q[2],q[3];
u1(2.23672555922905) q[3];
u3(0.336220732066233,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.52762236211911,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.14552080149415,0.310379107234355,2.68044238191172) q[3];
u3(0.864500321966995,0.682193522095604,-2.96419069727793) q[2];
u3(1.88209182422022,3.39603800250625,-1.66545953947074) q[4];
u3(0.633255917635508,1.35963359853578,-0.294711409704550) q[6];
cx q[6],q[4];
u1(0.131552575401040) q[4];
u3(-1.84623721351979,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.617967782244950,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.94451524790476,2.37580684971275,-3.23049037011445) q[4];
u3(0.865960716591412,3.90317329910211,-0.799478776278184) q[6];
u3(1.08255427430376,-0.805766525253410,-0.929815257609004) q[12];
u3(1.93207757468591,-5.00718418124053,0.640697177129860) q[11];
cx q[11],q[12];
u1(-0.309871634275233) q[12];
u3(1.04520328385188,0.0,0.0) q[11];
cx q[12],q[11];
u3(3.19822649798879,0.0,0.0) q[11];
cx q[11],q[12];
u3(0.304275835568273,1.80933362068873,1.21151556264748) q[12];
u3(2.67997384375252,2.40454922904430,3.10844228672198) q[11];
u3(2.04134488685718,-2.23209416699158,3.33568496949974) q[8];
u3(0.297004724810047,0.810500967334020,1.02423181847754) q[0];
cx q[0],q[8];
u1(-1.26485336405992) q[8];
u3(0.712036215172143,0.0,0.0) q[0];
cx q[8],q[0];
u3(3.56497220268087,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.84068525581853,-3.26816132810956,1.06108860240669) q[8];
u3(1.17951689506070,-0.0617290280850357,-1.46478626049503) q[0];
u3(2.45368760739367,-3.03416576815949,1.09290833918971) q[7];
u3(2.59364918796807,-0.0760018431745547,1.92255897340249) q[1];
cx q[1],q[7];
u1(1.88684049524752) q[7];
u3(-2.99073172633385,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.612711091449661,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.17100952618608,0.551588563219575,2.08081522109697) q[7];
u3(2.23103160079568,4.08032185857708,1.20493576207711) q[1];
u3(1.40403312125958,3.38212434030193,-2.32123100006582) q[12];
u3(1.34266431801487,2.76000470160947,-2.36600263815409) q[11];
cx q[11],q[12];
u1(-0.248956442173424) q[12];
u3(-2.02261457736414,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.58867859186775,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.02776518856014,-1.47930316150428,3.26754650651589) q[12];
u3(1.24781138114021,-2.06029618275722,2.04029118198266) q[11];
u3(2.14557802717663,-1.29798332665953,-1.22614654726002) q[10];
u3(0.481746862708598,-4.77123562568140,0.371733013264468) q[2];
cx q[2],q[10];
u1(2.59564895622348) q[10];
u3(-2.98195001120254,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.12135809192812,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.606379691572881,-0.0478941263935049,-2.47021326805172) q[10];
u3(0.987244689594927,-0.135398744333237,-1.90203019554512) q[2];
u3(1.23631320490740,1.10818026621037,0.227443134006013) q[9];
u3(0.909323553997732,-0.749470855152235,-1.15338993339792) q[0];
cx q[0],q[9];
u1(2.90768974844891) q[9];
u3(-1.60974276615991,0.0,0.0) q[0];
cx q[9],q[0];
u3(0.813847064169042,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.36351145228078,-0.535686485091671,2.28537112897266) q[9];
u3(1.20355841084139,3.29606693067748,2.39809595270571) q[0];
u3(1.47471338713781,0.536550084031528,-3.65456419929564) q[8];
u3(1.77360719002082,-0.843929408607695,4.50246865729298) q[3];
cx q[3],q[8];
u1(0.717488676647478) q[8];
u3(-1.27921817828529,0.0,0.0) q[3];
cx q[8],q[3];
u3(3.02418979093287,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.20693634778819,0.216868752333389,-1.99903333844026) q[8];
u3(1.71358220625785,1.09307607446465,0.792554482581920) q[3];
u3(0.948052875591371,-1.29312303292073,0.767893292203083) q[6];
u3(1.11955838362467,-1.89005068042943,0.371319554429032) q[7];
cx q[7],q[6];
u1(3.06063318105455) q[6];
u3(-1.72355311060276,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.928181790915931,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.33371657764816,0.327921755248957,0.808879450941001) q[6];
u3(2.62025948012992,-2.37042760107545,-2.73508004537918) q[7];
u3(1.08661803478840,0.0180279184616342,0.782062539054494) q[4];
u3(1.43351416894835,-1.02319648082861,-1.78963919960434) q[5];
cx q[5],q[4];
u1(1.09763443987736) q[4];
u3(-1.60918092552251,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.222328864257892,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.28321622200109,-2.28392076937209,-0.776229272290060) q[4];
u3(1.06059335270762,1.22157365623654,-2.51879549426672) q[5];
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
