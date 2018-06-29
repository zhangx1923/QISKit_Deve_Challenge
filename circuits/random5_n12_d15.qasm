OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.70362332402437,2.13252431339230,-2.24861100323434) q[7];
u3(1.21598455047687,2.23604671362552,-3.25343025830404) q[8];
cx q[8],q[7];
u1(2.27234447662100) q[7];
u3(-2.77953706211455,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.00831558874768,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.240445468350136,3.76188521593251,-1.72662408701368) q[7];
u3(2.46732375258038,-2.05223300237694,-0.220485639965455) q[8];
u3(2.35071294806214,3.12315501079768,-2.35904378808282) q[2];
u3(1.92196431625779,2.33385297994713,-2.02158078383230) q[5];
cx q[5],q[2];
u1(2.25478815721704) q[2];
u3(-1.69834440803278,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.55106592430097,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.28470289357124,-0.892377197415539,3.73781078366109) q[2];
u3(1.26038784870216,0.441346398658728,-3.25454777163468) q[5];
u3(0.577162446583234,-0.219076622281943,-0.688180087817559) q[11];
u3(1.65155126865234,1.85363589889723,-4.13933526777936) q[0];
cx q[0],q[11];
u1(1.81841963338049) q[11];
u3(-2.71725068509721,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.877181617627095,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.23014098459268,1.10100295114763,-1.38632544670626) q[11];
u3(1.44730960347144,-0.460250613407596,1.91897253520268) q[0];
u3(2.98753710307420,-0.656966331000237,-0.399967921906418) q[1];
u3(0.357387775470213,-3.36942456933952,-1.58262404044544) q[10];
cx q[10],q[1];
u1(0.824740493790075) q[1];
u3(-1.21236364289245,0.0,0.0) q[10];
cx q[1],q[10];
u3(2.75437409388743,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.777037848883424,0.111645244417057,-0.862966269789445) q[1];
u3(1.08172731999412,-1.81119580028185,4.12191187332802) q[10];
u3(1.38170276444385,1.12329968376114,-1.78820214209268) q[6];
u3(1.47672156082739,-4.55246645672515,1.63120571071316) q[4];
cx q[4],q[6];
u1(1.18003760447719) q[6];
u3(0.0704843943367361,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.46942740212878,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.65394991969190,0.763109618183936,-1.39633473412937) q[6];
u3(0.989954634860140,1.56613489331648,-1.95851619609118) q[4];
u3(1.06938977561677,-1.30493038507009,0.278621596972800) q[9];
u3(1.08998344874754,-2.09558511588578,1.00372264722178) q[3];
cx q[3],q[9];
u1(3.10846003046633) q[9];
u3(-1.54990869176745,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.338614933264481,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.97814924713190,-1.93093194275678,0.268394241140327) q[9];
u3(2.33599110822882,0.526317709871735,0.719440452899283) q[3];
u3(0.453662546061733,0.970440025207030,0.435880831650387) q[4];
u3(1.08871658832255,-0.690490809249906,-1.39600946045593) q[9];
cx q[9],q[4];
u1(1.47816635107900) q[4];
u3(-0.655256873871791,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.86369450017349,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.80051273548372,-2.76040760232810,2.99398405057102) q[4];
u3(1.20999356402334,-0.292223455746127,1.79998268193201) q[9];
u3(1.52367186625190,1.49043475593352,-0.000470934879136187) q[3];
u3(0.427334217132810,-0.293911070837867,-2.74197443090338) q[1];
cx q[1],q[3];
u1(0.739523968298718) q[3];
u3(-1.31347179369772,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.0221692955785544,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.96709946576238,-0.710197572965925,-2.88838056861343) q[3];
u3(2.21965367949538,2.75318963512534,3.52237810660762) q[1];
u3(2.45838247064821,0.179305671897476,2.78659679603684) q[11];
u3(2.71959643555391,-1.76141336593173,0.592970622557159) q[7];
cx q[7],q[11];
u1(2.14983728446544) q[11];
u3(0.355806772326084,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.35874831527698,0.0,0.0) q[7];
cx q[7],q[11];
u3(2.57869859264498,-1.73906622789295,-0.837142813702151) q[11];
u3(1.35360823781516,-2.57821666009054,1.86836105520967) q[7];
u3(2.75052325320336,-2.97424752496838,0.228166369443007) q[5];
u3(2.65986701422812,-3.37263557782043,-2.30274806316018) q[8];
cx q[8],q[5];
u1(0.699361002724838) q[5];
u3(-1.37264105140500,0.0,0.0) q[8];
cx q[5],q[8];
u3(-0.216376258854961,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.70282739066503,0.871767459097008,-2.69769474901453) q[5];
u3(1.40342459106741,0.792492366005543,-2.78947575565642) q[8];
u3(1.03370445707144,0.413217083088397,-2.89396786891867) q[0];
u3(1.27071791857693,3.00096190525772,-2.96061686538278) q[2];
cx q[2],q[0];
u1(-1.03529255038428) q[0];
u3(0.0628948248963612,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.21126259438426,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.02319266456014,0.0128110414080804,-2.08384036567758) q[0];
u3(0.990594238078927,-5.49611863929958,-0.101849668463143) q[2];
u3(2.49328441331255,2.13772485906508,-2.55150518368720) q[6];
u3(2.15650753582793,1.93807854854523,-3.18881350340762) q[10];
cx q[10],q[6];
u1(3.65037669992009) q[6];
u3(-1.31911068104093,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.09432175896404,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.64051036457458,3.66635919303227,0.227112589843064) q[6];
u3(1.87238811961664,-0.500960052013741,5.30060349343436) q[10];
u3(0.925767954268323,0.918651251151478,0.888139424632275) q[6];
u3(1.29681455144660,-0.410492654514490,-2.60722618500847) q[7];
cx q[7],q[6];
u1(2.50672034844844) q[6];
u3(-3.05147477118383,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.57509935821578,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.887729832785256,0.928013068570800,0.228444446160144) q[6];
u3(1.23251637554195,-3.31139357716990,1.94540288210494) q[7];
u3(2.10095608944580,1.14915022082214,-2.85629365092774) q[1];
u3(1.35528187348139,-3.04585019152140,2.77212954486374) q[8];
cx q[8],q[1];
u1(3.90646437862020) q[1];
u3(-4.57316840484815,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.898129948405678,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.33192267995005,1.34427634094527,-0.942435472136557) q[1];
u3(2.78730039504758,-0.537163293474074,1.05489682919918) q[8];
u3(1.46116810161953,1.06968588617097,-3.50471734370012) q[5];
u3(2.63710582762164,-1.89160201146712,4.22523834135972) q[10];
cx q[10],q[5];
u1(1.62089774177689) q[5];
u3(-2.69737687617030,0.0,0.0) q[10];
cx q[5],q[10];
u3(3.31259856828352,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.65291661848365,-3.59565142817170,1.98095222531704) q[5];
u3(2.42143308653816,2.47637605612038,1.36996945239852) q[10];
u3(0.940863518886608,3.64600896604856,-1.35726674398973) q[3];
u3(0.889484721122513,1.62382879251854,-0.742195176226693) q[2];
cx q[2],q[3];
u1(2.57643307749155) q[3];
u3(-2.98396097857673,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.50920364629347,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.43893171190625,-1.41902623581151,3.03180845686567) q[3];
u3(2.66164446378629,-3.90998647684010,1.84069205346162) q[2];
u3(1.79436477380563,1.18704980765601,-0.539496582541894) q[0];
u3(2.60202643674883,-0.654599233203027,-3.95085803935265) q[9];
cx q[9],q[0];
u1(1.40277810862119) q[0];
u3(-0.631083074100119,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.94727634106688,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.19761893048231,-0.274298736403182,1.29964825919322) q[0];
u3(1.22790294929409,-2.11297513680988,2.03386229917637) q[9];
u3(0.734257019126947,-0.498890743826540,-0.275464327286260) q[4];
u3(1.08037706649771,-2.01004562009482,0.638583981236335) q[11];
cx q[11],q[4];
u1(-0.641112446253428) q[4];
u3(1.26961932101767,0.0,0.0) q[11];
cx q[4],q[11];
u3(4.05749027238666,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.75104131181962,2.13281516287043,-1.59965632494666) q[4];
u3(0.835176842762193,-1.15599955231853,1.34611118283814) q[11];
u3(2.38518071392386,2.39890278546794,-2.69069851052747) q[6];
u3(1.78797643941673,2.49947642585149,-2.89323907252376) q[9];
cx q[9],q[6];
u1(1.96720746470854) q[6];
u3(-3.11250511338916,0.0,0.0) q[9];
cx q[6],q[9];
u3(1.84344069602178,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.84646135679517,-2.34701337934381,-1.63077564854524) q[6];
u3(1.04421114184924,-1.97346035330326,-0.895075212395355) q[9];
u3(1.84439005300223,0.386862995171818,-3.10243191576608) q[1];
u3(2.80463149697812,2.88320445387440,-1.93226452662734) q[5];
cx q[5],q[1];
u1(3.33590030461906) q[1];
u3(-1.33619097936793,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.47361324867457,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.818982920607086,-2.28624841565182,0.248682897215540) q[1];
u3(1.20168382352837,-2.65928327269005,2.79966301211755) q[5];
u3(1.87184569604606,3.46471748310836,-2.33242258852286) q[4];
u3(2.48440781911290,1.36785799421263,-2.45526163849302) q[3];
cx q[3],q[4];
u1(-0.661766356732823) q[4];
u3(0.346090583326249,0.0,0.0) q[3];
cx q[4],q[3];
u3(4.16376902188906,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.88460233525638,3.32898899333948,-2.26266366531529) q[4];
u3(0.424828914430415,0.133444688417910,-0.233162522144793) q[3];
u3(2.87566702494325,-2.44941122411828,1.85388313514722) q[0];
u3(2.22843554646797,0.648766923433673,2.48171412225426) q[8];
cx q[8],q[0];
u1(0.859089874548153) q[0];
u3(-0.209127747704924,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.88904051187120,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.921481531130540,-1.97726056952348,-0.536587490378950) q[0];
u3(1.46199326600735,2.26065061244259,2.37611911070145) q[8];
u3(1.93583998754900,2.05020940660095,-1.88086590398191) q[10];
u3(2.84963340700486,0.371708513820884,-4.76260219176841) q[7];
cx q[7],q[10];
u1(1.91996917712951) q[10];
u3(-3.44784483246217,0.0,0.0) q[7];
cx q[10],q[7];
u3(1.12870267218906,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.62228454237049,0.573529655788002,-3.54182532362659) q[10];
u3(1.22406456292754,-0.289011964485368,2.26166952910353) q[7];
u3(0.970502524555416,3.05373442443348,-1.70479747770515) q[2];
u3(1.84478269806576,2.12628123872460,-2.02246571065668) q[11];
cx q[11],q[2];
u1(1.76386989422634) q[2];
u3(-3.00505881286086,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.325157743955648,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.56909414960659,-0.479738472303301,2.34656536484302) q[2];
u3(1.84158516031653,-0.267086569815434,-2.18242147921420) q[11];
u3(0.211778754125989,1.60127635221768,-1.15171304852901) q[6];
u3(1.06543178275588,-3.87950084884339,1.19876605132371) q[10];
cx q[10],q[6];
u1(0.597101910705824) q[6];
u3(-3.15838287781362,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.66657542586190,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.19999583490396,3.80911842059646,-1.59491015729857) q[6];
u3(0.816093230599771,-0.492130988409835,-4.38140765972898) q[10];
u3(1.81599619923935,1.47374280057890,-0.0562556095333746) q[1];
u3(2.19559845630531,1.00958678348729,-1.75467517284102) q[5];
cx q[5],q[1];
u1(-0.471087110271953) q[1];
u3(-2.27692441711715,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.23367805653570,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.69386067707958,-0.697373244203456,2.52609446725450) q[1];
u3(1.05051709481715,-1.82203639875628,2.44146608568443) q[5];
u3(1.99414983962073,1.12283673489778,-0.461272770506628) q[11];
u3(1.35026040645158,-4.35131400389467,1.06770447822120) q[0];
cx q[0],q[11];
u1(1.29328302239085) q[11];
u3(-0.460837618635256,0.0,0.0) q[0];
cx q[11],q[0];
u3(3.03908595260099,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.964505038102163,-1.79159401218855,2.02591065412606) q[11];
u3(0.631062222935705,-1.21238689766462,4.53590862027015) q[0];
u3(1.44974336794503,3.17987378684130,-0.354961647196665) q[9];
u3(0.643324414942779,1.43266865699315,-1.39466208957069) q[8];
cx q[8],q[9];
u1(1.71118809727398) q[9];
u3(-2.34045401131125,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.250298400093432,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.25761843255734,-4.13958164334391,2.03357296785713) q[9];
u3(1.09259693894183,-0.427209291681059,-4.11959914935065) q[8];
u3(1.59338657513683,-1.00309433536027,1.59669497875012) q[3];
u3(1.89024745806069,-1.15430655023224,-2.40981139095168) q[7];
cx q[7],q[3];
u1(1.20898909701223) q[3];
u3(-3.23482364840821,0.0,0.0) q[7];
cx q[3],q[7];
u3(2.48223708043286,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.20054984318730,2.17447149256722,1.76262401477899) q[3];
u3(1.31796168107031,-1.34604994636405,3.13208735583201) q[7];
u3(1.61183102241405,-0.371929046631088,0.794890875677707) q[2];
u3(1.49581652593619,-2.71592227026746,-1.16968989409877) q[4];
cx q[4],q[2];
u1(1.59276683939081) q[2];
u3(-1.93103768300377,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.291969345427993,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.20274086187548,-4.78176297995091,0.898035811602127) q[2];
u3(2.13779844545323,-2.70410880718138,-0.537796707307327) q[4];
u3(1.90112523671696,2.79107078739563,-0.459818527796118) q[9];
u3(2.79735273282784,4.54944219966263,-0.674527553341992) q[8];
cx q[8],q[9];
u1(1.26909175988690) q[9];
u3(-0.165188322793323,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.22820312808953,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.28076440437565,-2.32871345509514,0.865014271895128) q[9];
u3(1.06229118193564,-0.579722927834460,-4.43359824271215) q[8];
u3(1.88251199736804,0.928246673204712,-0.0428314630246173) q[7];
u3(2.14637869313100,0.0980563059579387,-4.34971471354471) q[0];
cx q[0],q[7];
u1(1.83479406461187) q[7];
u3(-2.68421970362135,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.585931884482769,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.87520399291238,2.06113463550642,-1.20036124645893) q[7];
u3(2.03672619343305,0.930370738036596,-3.29802535326177) q[0];
u3(1.24026274126013,-0.324564210588256,0.212720205257500) q[2];
u3(0.813028148187882,-2.16039085516991,-1.25852330348737) q[5];
cx q[5],q[2];
u1(1.42064668609110) q[2];
u3(-0.422498218808261,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.230217835421820,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.87947272517798,-2.36183057843950,0.164462486574787) q[2];
u3(1.15475472992311,-0.428814333480984,3.33843293784295) q[5];
u3(2.16126395538811,-0.0693263279453061,-2.05431755991996) q[11];
u3(1.85394994172291,0.668830657616365,-4.15823567468333) q[10];
cx q[10],q[11];
u1(1.26851485361857) q[11];
u3(-0.966725445104750,0.0,0.0) q[10];
cx q[11],q[10];
u3(-0.109682453105307,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.67102593025909,2.80741060177061,-2.82943101526522) q[11];
u3(2.07245541527358,1.31407909834476,-3.25073489141646) q[10];
u3(3.02834752107610,-1.59562949972654,-0.577598953873241) q[3];
u3(1.34162029337192,-0.483358995427945,-4.12469981740859) q[4];
cx q[4],q[3];
u1(1.67661034092241) q[3];
u3(0.643105719481772,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.834833473574459,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.65417656427060,3.07329299387961,-2.39084840902500) q[3];
u3(0.370030368092338,-1.49552501286069,4.21974262363006) q[4];
u3(1.52255467255771,-1.47648329683711,-0.418723266461381) q[1];
u3(1.28764978740930,-4.18575769978536,0.281196099395888) q[6];
cx q[6],q[1];
u1(1.12380825514633) q[1];
u3(-3.30603351567866,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.44578284790978,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.547935541350847,1.86603266512231,-3.76701856540204) q[1];
u3(0.860220986157920,3.49031176197219,-0.0868606152291860) q[6];
u3(1.73713409886405,-0.360691914680641,-2.64000605757049) q[1];
u3(2.30626352658118,0.531545885823005,-4.74938899026714) q[9];
cx q[9],q[1];
u1(2.01125304565717) q[1];
u3(-2.69712180436555,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.663935715422027,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.33368978694539,3.06138312873472,-2.71106479747597) q[1];
u3(0.985575112737528,5.12508799506679,-0.641024803102533) q[9];
u3(1.11181588163259,0.414082982533603,1.33590484857008) q[5];
u3(1.38832683292516,-1.25348182744509,-0.339633939082854) q[0];
cx q[0],q[5];
u1(-0.151800542348345) q[5];
u3(-1.14437992991941,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.52616402199735,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.42672141896105,3.00748277408505,-2.44795137214709) q[5];
u3(1.83144676019190,1.83276063950802,0.207552648577989) q[0];
u3(2.63277660213948,3.30013692585684,-0.844838765987814) q[8];
u3(2.45424114990283,0.859075809577730,-3.53317763903566) q[3];
cx q[3],q[8];
u1(1.28213393456366) q[8];
u3(-0.755063366593005,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.60038826242965,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.01616324965130,-1.96135208944849,1.67551365743865) q[8];
u3(1.41365463211746,-1.55905030470860,-2.32182132364649) q[3];
u3(0.451230576364369,-0.170297138220299,-0.0602498754675558) q[11];
u3(0.953288963279939,-2.20275558819783,1.75221362664715) q[2];
cx q[2],q[11];
u1(3.04446421094839) q[11];
u3(-1.48216395286863,0.0,0.0) q[2];
cx q[11],q[2];
u3(2.15387330801973,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.45620372011074,-1.48861830801154,-0.0481124143983611) q[11];
u3(2.16814250233019,-3.63986084869599,-1.35603481838843) q[2];
u3(2.19597590634154,0.188024868418208,-2.49855622053529) q[4];
u3(0.452303381006681,1.05146680323887,-3.97023308541525) q[6];
cx q[6],q[4];
u1(0.0943790069454544) q[4];
u3(-1.01609228317676,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.52564755724824,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.677716314287932,-2.94513260590343,-0.122825843240409) q[4];
u3(1.40876776485864,-3.04943400274909,2.60623925342337) q[6];
u3(0.400791953406733,-0.377851386550273,0.565335873751108) q[7];
u3(0.563460179543520,-3.63511862665590,0.960329564679002) q[10];
cx q[10],q[7];
u1(2.05088998535950) q[7];
u3(-2.99641101482238,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.614935878832908,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.90430183454792,-1.11117798205795,0.549594311772390) q[7];
u3(2.24298498555804,-0.518514437997993,-0.653537542039609) q[10];
u3(0.957327259390202,1.00183005507531,-2.60140014158479) q[3];
u3(1.45992779394432,2.54977123728936,-3.00434255119025) q[11];
cx q[11],q[3];
u1(0.240042888615821) q[3];
u3(1.11509386074178,0.0,0.0) q[11];
cx q[3],q[11];
u3(3.06220254670253,0.0,0.0) q[11];
cx q[11],q[3];
u3(0.540176860826539,-4.45499150583619,0.580975584564425) q[3];
u3(1.71880753086150,-2.47035236879172,0.257927824627625) q[11];
u3(0.263396592068981,1.73184177274007,-1.26642109399382) q[6];
u3(1.19069051820175,-2.60332933517197,1.25912213692534) q[1];
cx q[1],q[6];
u1(0.294252668146034) q[6];
u3(-0.555434248529890,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.98804583065521,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.94252653785314,-0.168755356184464,-3.82049909611794) q[6];
u3(0.866353287429091,-0.559320498304597,-0.769505400718914) q[1];
u3(0.609311673056537,-1.99440404337960,2.58587238484940) q[0];
u3(0.196321486235913,-3.32827064031388,2.11332854210553) q[2];
cx q[2],q[0];
u1(-0.0920777198585894) q[0];
u3(-2.23256081215070,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.52070300198331,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.33680244715080,3.63900200772463,-0.548151765890213) q[0];
u3(2.48694648827197,4.13619440350273,1.28158575024328) q[2];
u3(1.12897970265726,-2.77062787620755,0.449779970481185) q[10];
u3(0.695830766091448,-3.02922094513484,-0.245970944904187) q[7];
cx q[7],q[10];
u1(1.47959121967605) q[10];
u3(-3.35254004829609,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.25753152151925,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.47411837136116,-2.78049323512796,1.86024804814388) q[10];
u3(1.98357060342700,0.385323905199485,4.02740651105367) q[7];
u3(2.20802195295463,0.854643472220930,-3.33389537376365) q[5];
u3(2.44049170381670,3.52567844060730,-2.26029848925349) q[9];
cx q[9],q[5];
u1(3.02750873971787) q[5];
u3(-2.06389483447307,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.986752027670757,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.519268200118218,-0.119654609691193,-1.87988941505975) q[5];
u3(2.53842037588518,2.08409592323901,-3.07657104501011) q[9];
u3(1.83311157814736,0.469066635171756,2.28672933792234) q[4];
u3(1.39285752986347,-2.99981287364693,-3.27243627237519) q[8];
cx q[8],q[4];
u1(1.58581819704362) q[4];
u3(-0.0398362221821804,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.00282173192212,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.36035930496009,2.90189790036903,-3.01611316016558) q[4];
u3(0.397045295213092,-2.32601078809551,3.86207763836954) q[8];
u3(2.05845193917816,-0.412359876777434,1.75621796452945) q[3];
u3(1.43703157121056,-2.58156137466336,-0.274182215584586) q[4];
cx q[4],q[3];
u1(3.45057910984437) q[3];
u3(-1.10117524158298,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.87983099869647,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.30870327753604,0.117701099643122,1.60788917979269) q[3];
u3(2.72890943681061,0.405094775166550,-3.49792316760065) q[4];
u3(0.874241211754577,-0.958625960669079,1.64846249446380) q[2];
u3(0.577420994158089,-0.608226542362527,-1.20291820162638) q[7];
cx q[7],q[2];
u1(2.61098983863418) q[2];
u3(-1.38453435039962,0.0,0.0) q[7];
cx q[2],q[7];
u3(-0.0715925274792029,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.89016418385853,-4.16114671775452,-0.161875076280809) q[2];
u3(2.26838688484369,2.57947851985390,3.26992895983583) q[7];
u3(2.90983559388871,-1.05938959822905,-0.201788366287515) q[5];
u3(1.40620058590148,-0.313418112710832,-4.49199886232065) q[6];
cx q[6],q[5];
u1(2.92854423466558) q[5];
u3(-2.17399594921925,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.621551600997438,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.56340828519854,2.01613376799232,-1.71120527366177) q[5];
u3(1.02517870125941,-1.06769824136457,5.01933656471616) q[6];
u3(1.83291796328261,1.22664462491068,-3.93372699714654) q[9];
u3(1.55937419901192,-2.16536622027860,3.39363143794338) q[0];
cx q[0],q[9];
u1(0.0621902425556571) q[9];
u3(-1.15676489054407,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.51055140696025,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.76768790354064,0.544017248634166,0.134927499901171) q[9];
u3(0.810052182551031,0.536966768557339,0.486584403794925) q[0];
u3(2.21483737097280,2.53142011239463,-0.133234486486614) q[10];
u3(2.26655773559786,-0.393618466674606,-3.42182543093426) q[11];
cx q[11],q[10];
u1(4.26903426936465) q[10];
u3(-3.89483750019832,0.0,0.0) q[11];
cx q[10],q[11];
u3(-0.515479369391189,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.81109824424980,-0.525841867452669,3.50278369411442) q[10];
u3(0.465894719591923,3.58057743245660,-0.340231887276200) q[11];
u3(1.93384843591517,1.21691279814417,-2.64585571756483) q[1];
u3(1.89101399321769,-2.71264045813416,2.41829469907039) q[8];
cx q[8],q[1];
u1(0.924132182284716) q[1];
u3(-1.69963002132396,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.91095973417930,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.96219139664848,0.247160986816012,-2.34429802021073) q[1];
u3(1.82115635690041,-1.69222774285033,-2.36323398313384) q[8];
u3(1.76936967493555,2.29231600865614,-0.609644599199975) q[4];
u3(2.01088355029535,0.149562477934668,-2.90945878103616) q[0];
cx q[0],q[4];
u1(0.965815362866227) q[4];
u3(-0.621924605469206,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.74622670988175,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.37659831532844,0.712363106863913,-0.262221089964895) q[4];
u3(1.03038757783321,0.702396783702790,-3.72470738851812) q[0];
u3(2.24268419842177,0.425730761424592,-1.54300391849080) q[2];
u3(2.51909108789351,3.91114448860670,-0.815883119082617) q[10];
cx q[10],q[2];
u1(3.28980020316344) q[2];
u3(-1.51600522979669,0.0,0.0) q[10];
cx q[2],q[10];
u3(2.05580079084734,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.25638493068440,0.0154880845577325,0.357168320671386) q[2];
u3(2.16796147032574,1.58160658349838,4.68108003825539) q[10];
u3(0.824393669987784,1.63859475257811,-0.285092355670521) q[6];
u3(1.05081412906292,0.366127855259611,-2.50335113455660) q[1];
cx q[1],q[6];
u1(1.59878859490704) q[6];
u3(-2.34539122733048,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.0863793579156111,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.04512803857112,-1.39308433313531,0.511192863790636) q[6];
u3(0.719693449237037,-1.27816701812467,3.00498288286453) q[1];
u3(0.885466258727636,-0.297583309532820,2.61187841233444) q[5];
u3(1.94924001227483,-2.43729005187227,-1.80056172014748) q[7];
cx q[7],q[5];
u1(3.26272587325597) q[5];
u3(-1.43980824446230,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.63735307582984,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.885102696849237,1.14200822290476,-0.593989870073506) q[5];
u3(0.928236709775145,1.40925219000992,3.40131967337221) q[7];
u3(0.729979782301758,1.74984846130410,-0.869245893184901) q[9];
u3(0.517902181086374,-2.12383186701753,0.947061585538294) q[3];
cx q[3],q[9];
u1(0.434507579046550) q[9];
u3(-0.152941569039109,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.97133756957393,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.16451181173017,4.08455117081699,-0.837450033007834) q[9];
u3(2.86776262548854,1.74334764456114,2.09444454872893) q[3];
u3(0.484190244817751,0.0232649445795906,-1.71223087238034) q[8];
u3(0.905068770004298,0.725234737187336,-5.07109355777467) q[11];
cx q[11],q[8];
u1(2.62833788826565) q[8];
u3(-2.13274121939105,0.0,0.0) q[11];
cx q[8],q[11];
u3(-0.123384176532503,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.91551230056277,-0.613241896367416,-2.47328005971020) q[8];
u3(1.07861208542418,2.37593784112905,-1.43692116723032) q[11];
u3(1.41829269242617,0.138258320530241,2.25866096521916) q[11];
u3(1.70699715462680,-1.64923264540390,-0.855438448154063) q[0];
cx q[0],q[11];
u1(1.90670416516433) q[11];
u3(0.244585116991309,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.613492222493291,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.15436502632051,1.23711116077251,-2.10606177708143) q[11];
u3(0.998835199350965,1.74079945604509,-3.14925859028570) q[0];
u3(2.66946322602063,0.0162497478369809,-0.770246660690272) q[7];
u3(0.998504189826902,0.0979907800505724,-5.31307817975078) q[6];
cx q[6],q[7];
u1(-0.381088420480574) q[7];
u3(-1.77377869081024,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.22342446022049,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.68729353388431,2.90972850597746,-1.98380866044866) q[7];
u3(1.74023804141556,4.62980125557481,0.895524138447176) q[6];
u3(2.71940739917292,1.43303313536369,-2.46384855479070) q[9];
u3(1.94010657880341,2.14410479584164,-3.23861088051953) q[1];
cx q[1],q[9];
u1(2.10911502231259) q[9];
u3(-1.60517195552000,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.40289243229386,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.99900788620247,-1.61770165800430,0.585877386131582) q[9];
u3(0.698231641282164,-4.45392222220806,0.340624741718813) q[1];
u3(1.97498585978274,2.39977645586681,-3.10983108214174) q[2];
u3(1.79807619971262,1.46573484586679,-1.61103871637379) q[10];
cx q[10],q[2];
u1(2.44578830008644) q[2];
u3(-1.79234490129234,0.0,0.0) q[10];
cx q[2],q[10];
u3(3.27898983902057,0.0,0.0) q[10];
cx q[10],q[2];
u3(0.848527657418328,3.36264496927856,-0.887539800250504) q[2];
u3(1.44010117376489,-0.417067814112605,3.65666109892695) q[10];
u3(1.14470333155334,0.355495842562916,1.47118024310888) q[3];
u3(2.24887332019637,-2.71709527186289,-1.22305501232972) q[8];
cx q[8],q[3];
u1(0.971976708596430) q[3];
u3(-3.08080157151823,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.57058183432398,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.35787200840784,1.27818869997688,-3.60279834088447) q[3];
u3(2.52413625959963,-4.95237302593422,0.0257365074696163) q[8];
u3(0.482482335821546,1.43473431797307,-0.584077165278645) q[5];
u3(1.01232838084597,0.177992094517451,-1.54488745869483) q[4];
cx q[4],q[5];
u1(0.129985483517948) q[5];
u3(-1.54038516564827,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.62857510741070,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.93562990909692,-1.93483171104541,-1.29763098419704) q[5];
u3(0.983505836797524,3.21508991389535,1.16082751493386) q[4];
u3(2.99171192329806,-0.696225573354491,1.44921008289542) q[5];
u3(1.83167800928226,-1.39154680030416,0.239456905082525) q[6];
cx q[6],q[5];
u1(0.730934805115847) q[5];
u3(-1.02648068270899,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.0418416570660105,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.594542218418383,0.213734279514165,-0.934867671932262) q[5];
u3(1.25296621324021,-0.0497123114529564,-0.239293278295307) q[6];
u3(2.98591700665428,0.595252724999076,-2.42857116214777) q[7];
u3(2.24507147995156,0.865878243828879,-2.61252141339635) q[11];
cx q[11],q[7];
u1(1.79891919148089) q[7];
u3(-2.60493279347964,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.211713242564481,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.23693943212259,-1.48231337029532,4.32681076235104) q[7];
u3(1.89501033301613,1.58049777436923,-2.08448532181656) q[11];
u3(1.77984466860609,-0.578089281024972,-1.33463309555311) q[10];
u3(0.862349986650461,-4.94724702547320,1.15313948077646) q[9];
cx q[9],q[10];
u1(0.661844689598852) q[10];
u3(-0.313007437436466,0.0,0.0) q[9];
cx q[10],q[9];
u3(1.37205788937469,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.71214726445173,-1.20633389380063,3.31084269351327) q[10];
u3(1.75404884174373,1.02268804011161,1.86839335119165) q[9];
u3(2.73270566449742,1.84606943991130,-3.19168804374033) q[3];
u3(1.82561161997137,2.37158647255617,-2.50859870789057) q[0];
cx q[0],q[3];
u1(2.64993991455171) q[3];
u3(-1.58305117557341,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.0867163705809926,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.95035756380736,-3.82887211590304,2.39212807366865) q[3];
u3(1.35790203091699,4.45850384058635,-0.00416548030530706) q[0];
u3(2.48872578818422,1.51327647422874,-4.25125173745731) q[4];
u3(2.00625356698372,-2.39199596982902,3.62261164864536) q[2];
cx q[2],q[4];
u1(3.65090378814798) q[4];
u3(-1.46836196191512,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.50115332697125,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.17567739258656,-1.06765392208542,-1.50696678564713) q[4];
u3(1.27442097725293,1.37172300220115,-0.983900144439307) q[2];
u3(1.53788639035512,0.839507597496567,-3.06570770799049) q[8];
u3(2.22451579755565,3.02420244397638,-2.76857375960729) q[1];
cx q[1],q[8];
u1(3.56309875438720) q[8];
u3(-4.55402407295115,0.0,0.0) q[1];
cx q[8],q[1];
u3(-0.538154391088022,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.21747799744720,-1.18860829493179,3.37138559734225) q[8];
u3(1.64452142991674,-2.10638660485771,-1.47703595161423) q[1];
u3(2.27671635406813,-0.561751969244924,1.91392891778685) q[10];
u3(2.16823872480627,-1.39477520172895,-0.00983686823427110) q[9];
cx q[9],q[10];
u1(0.333933060928523) q[10];
u3(-1.27931196351341,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.0433916007064095,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.73935062898298,0.582611797355631,-0.820727198046292) q[10];
u3(0.825859584675321,-1.75232870703899,-0.303558488775160) q[9];
u3(1.15443503646213,-1.12520036105235,0.330538863921787) q[0];
u3(1.47211006018695,-1.96590091418917,-0.773382813002052) q[6];
cx q[6],q[0];
u1(0.689744993178434) q[0];
u3(-0.103471920302319,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.41116519334249,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.78129626479159,-2.88498436683873,1.25188650978873) q[0];
u3(2.03599690919663,-0.115600073782958,4.79176474991779) q[6];
u3(0.578394925411000,-1.63382400392227,1.76111021166567) q[5];
u3(0.0835618525374520,-3.98865652390540,1.23824954826543) q[4];
cx q[4],q[5];
u1(4.14589605158062) q[5];
u3(-3.26131089819201,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.547797535622459,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.995541667315637,-2.77841568843786,-0.257007844115251) q[5];
u3(0.788047858936579,-4.97098797100566,-0.417226341074396) q[4];
u3(1.56753378754711,2.27170685278379,-0.255863418674026) q[1];
u3(0.883931046944518,0.747445362765575,-4.40907753083883) q[8];
cx q[8],q[1];
u1(0.0438928532996066) q[1];
u3(-0.645923669327401,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.06402041739904,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.74488103483605,-2.67028889066212,0.304020939876707) q[1];
u3(0.521125727517423,-5.96219816924050,0.211203766241300) q[8];
u3(1.78601009795912,1.66136791218523,-3.54023812768443) q[7];
u3(2.32630595208336,3.74077159771699,-2.34104860399654) q[11];
cx q[11],q[7];
u1(3.31284454802133) q[7];
u3(-3.54991163329872,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.979066386488746,0.0,0.0) q[11];
cx q[11],q[7];
u3(0.402229077081877,-1.23563644199266,2.78918982140940) q[7];
u3(1.00293024165505,-1.36955084563587,-1.50343814717542) q[11];
u3(0.889499201966700,1.35444658422128,-1.98559139525412) q[2];
u3(1.54491162743943,-4.91339526167355,0.743707245209802) q[3];
cx q[3],q[2];
u1(-0.166534398605140) q[2];
u3(-1.81471452365030,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.734930575969817,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.22592841950417,-3.63543937615338,1.98499575539685) q[2];
u3(2.29624000606251,4.31907961820564,-1.21671762482463) q[3];
u3(1.71548063435244,-0.634185690771895,-1.55861660058447) q[2];
u3(0.606914041687123,-4.17626502574252,0.788727250456135) q[10];
cx q[10],q[2];
u1(0.537129028419782) q[2];
u3(-1.07202729675129,0.0,0.0) q[10];
cx q[2],q[10];
u3(2.35420343374785,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.32376544496752,1.23050360742999,-2.77372785226244) q[2];
u3(2.21697300094337,4.28686621975094,-0.602167935671052) q[10];
u3(1.50255018594433,0.603091138149320,-1.17681009769795) q[4];
u3(0.984883514555039,-3.69585050842743,1.32293885299696) q[3];
cx q[3],q[4];
u1(3.12632774535941) q[4];
u3(-2.08488595278305,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.22042145714244,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.40214142312245,2.54854100748962,-2.18401577463541) q[4];
u3(1.02210012844917,1.70613131540672,2.47677626832013) q[3];
u3(2.32259944280599,1.00252777421524,-0.416368971790443) q[9];
u3(1.12764447902450,-0.277144015923546,-4.39341755240665) q[7];
cx q[7],q[9];
u1(1.42629523189057) q[9];
u3(-3.45911747219076,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.07386986360648,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.97914343525018,0.429827919618645,0.147327817162936) q[9];
u3(2.36232422727513,4.00441434929571,0.501100875560925) q[7];
u3(1.12654859161461,2.25756590935698,-3.58695036126400) q[11];
u3(1.78446172241983,-2.60011231350153,3.39138221540138) q[8];
cx q[8],q[11];
u1(0.101609543065427) q[11];
u3(-0.808414143103326,0.0,0.0) q[8];
cx q[11],q[8];
u3(1.46155364433816,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.15057959023534,4.18724536173581,-2.06072484018856) q[11];
u3(2.16431138234222,-2.96010783879599,3.31322047472284) q[8];
u3(0.372339356884663,-0.570122212609378,0.918063719241312) q[1];
u3(0.245023320230005,-1.68159369596216,0.327942049288999) q[0];
cx q[0],q[1];
u1(0.0779703494350161) q[1];
u3(-1.83330284489140,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.73471330065615,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.63131035490129,-3.75924589413757,1.88495607838439) q[1];
u3(2.07694588140479,-3.50178480435228,0.709114988377068) q[0];
u3(1.22181943211810,2.37717489409506,-1.67266820530224) q[5];
u3(0.514105563487960,1.44875702232701,-1.13809311355075) q[6];
cx q[6],q[5];
u1(-0.890180813379766) q[5];
u3(0.266365336716874,0.0,0.0) q[6];
cx q[5],q[6];
u3(3.75943540732677,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.76480914164426,-2.54279957386334,-0.0750801698854302) q[5];
u3(2.03795810604669,2.64854052856921,0.559416013787107) q[6];
u3(0.290164335264874,2.95490845296393,-1.68835485551544) q[2];
u3(0.542602982729772,0.661309838679858,-2.38465739504023) q[1];
cx q[1],q[2];
u1(1.33186760719182) q[2];
u3(-1.53629796158382,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.770659467296973,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.564357467560321,-4.54503719187640,1.23040660841236) q[2];
u3(1.96713609523494,-1.59636529643597,-0.589295536300153) q[1];
u3(1.36934439368160,-0.846227202534541,0.360943466030544) q[5];
u3(1.52328441335900,-3.37442318562148,0.0739638081170064) q[3];
cx q[3],q[5];
u1(2.04359959572402) q[5];
u3(0.509492560511043,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.44468640290938,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.817112085925335,2.65012651015969,-1.98425657853362) q[5];
u3(2.31857547531567,5.45164323185555,-0.762867076167529) q[3];
u3(2.70566662631753,1.16978079055294,-3.35619732290821) q[8];
u3(2.44893971520560,1.77618837438287,-2.84666234405147) q[10];
cx q[10],q[8];
u1(0.664484925361845) q[8];
u3(-0.0416404903739214,0.0,0.0) q[10];
cx q[8],q[10];
u3(1.75759629603880,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.77841424971407,-0.00653617251817873,1.50086156885317) q[8];
u3(0.899527704810489,-1.02983664541317,0.0790421603563644) q[10];
u3(2.11592980641423,-0.701761426012178,2.12407839965888) q[0];
u3(1.73993391036467,-1.21207646342630,-1.15745380184696) q[6];
cx q[6],q[0];
u1(3.15916774246206) q[0];
u3(-2.33666420417388,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.860621346311648,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.380685126985417,1.65409101473799,-4.40760514671086) q[0];
u3(1.30696879755084,-4.37723857493318,-0.895197594910361) q[6];
u3(1.53612996318782,0.335149026808677,1.77990391184446) q[7];
u3(1.83253523279535,-1.33668947117265,-0.826155963138124) q[9];
cx q[9],q[7];
u1(1.14779378023762) q[7];
u3(-0.325670230465613,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.47926992873508,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.40103784283014,-3.67795783337495,1.20320639806539) q[7];
u3(1.94675987406127,1.95896816540312,2.22881603982139) q[9];
u3(2.50442775195316,0.947848244415285,1.37998181301655) q[4];
u3(0.962030328696702,-5.62696718246869,-0.104789140268837) q[11];
cx q[11],q[4];
u1(-0.890008840639206) q[4];
u3(0.180315177478378,0.0,0.0) q[11];
cx q[4],q[11];
u3(3.55267308617120,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.32638620429666,-2.16610576968939,0.896382418880605) q[4];
u3(2.04546312136644,3.51001616221282,0.863679364212706) q[11];
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
