OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.61223643056255,-0.197827718252671,-0.0587302262129262) q[4];
u3(0.872742571452392,-0.949744768813129,-4.17833823355318) q[5];
cx q[5],q[4];
u1(0.867105773814080) q[4];
u3(-1.29611001969056,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.226868696321367,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.54929808948573,2.01328159924362,-3.38215465217637) q[4];
u3(1.97450684864285,3.81553260466986,-2.35263393037920) q[5];
u3(1.84663307845467,3.21693398107064,-1.64256215725154) q[9];
u3(1.69777723171707,1.67928019128921,-1.93383595318337) q[11];
cx q[11],q[9];
u1(2.67840872585265) q[9];
u3(-2.44787742962383,0.0,0.0) q[11];
cx q[9],q[11];
u3(1.36727066917314,0.0,0.0) q[11];
cx q[11],q[9];
u3(0.455647048867434,-2.50357679889786,3.19473636964024) q[9];
u3(2.79921929665085,-0.445978807610744,1.77556812109893) q[11];
u3(0.429158545612923,-1.51874476875173,-0.178068856210937) q[3];
u3(1.60855146084189,0.970294516733916,-5.11986338051555) q[0];
cx q[0],q[3];
u1(0.723692253950993) q[3];
u3(-3.43630280644846,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.52423785038013,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.984091515478094,1.31646255045398,-0.912260369407078) q[3];
u3(1.74290513127974,-0.877192583854736,-2.14410418272044) q[0];
u3(2.25551702310734,-2.28459662626259,3.37769707613343) q[7];
u3(0.672847251192629,-0.375292972644565,2.55175143740965) q[8];
cx q[8],q[7];
u1(2.46465634524997) q[7];
u3(-2.13933180352235,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.09865915732972,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.47239108965204,-1.35924510596101,2.41785669250462) q[7];
u3(2.49446247685790,-3.23605993629310,-0.839494204516494) q[8];
u3(1.40980436196641,0.284819800513868,-2.38816071082010) q[10];
u3(0.347955781677389,-3.59151921094744,-0.206977939855055) q[2];
cx q[2],q[10];
u1(1.34915649211949) q[10];
u3(-0.460081554535380,0.0,0.0) q[2];
cx q[10],q[2];
u3(2.25905635201438,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.25189805386253,0.981537391376256,1.03350571772814) q[10];
u3(2.56911608318111,0.229621512125429,4.62894880299266) q[2];
u3(1.89776905191098,-0.717173437297709,2.77337708769872) q[6];
u3(2.00587614900498,-2.51681264566182,-1.76622726001159) q[1];
cx q[1],q[6];
u1(2.24406301406187) q[6];
u3(0.0186035827021831,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.952479961869840,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.81076738223887,-0.776743834515366,-0.638469359265664) q[6];
u3(0.580414345819705,-1.44959860151669,0.611456577552232) q[1];
u3(1.86426596161045,-0.0790949567229928,1.49896028278171) q[8];
u3(1.44546318028867,-0.529716563620215,-1.10607508237048) q[6];
cx q[6],q[8];
u1(1.52118144403917) q[8];
u3(-0.169350039770448,0.0,0.0) q[6];
cx q[8],q[6];
u3(2.18956484971045,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.920109069660406,2.63462569179395,-2.48121012946821) q[8];
u3(1.21901816644148,1.05641886466489,-4.18718755969669) q[6];
u3(2.77425983243433,2.11677104451147,-2.31088901448810) q[2];
u3(2.39181150007051,-0.243373980845282,-4.82708646485884) q[3];
cx q[3],q[2];
u1(0.306089049450314) q[2];
u3(-1.12660221301487,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.18497389515073,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.313806566280615,0.347389604942840,-0.0806441222374430) q[2];
u3(2.41744129317920,3.51220653607844,0.119569989611769) q[3];
u3(1.30059534072311,3.09671405878785,-0.522460492707301) q[11];
u3(0.919421102207722,0.776216936937927,-1.00608858646540) q[4];
cx q[4],q[11];
u1(2.49930466212893) q[11];
u3(-1.95605714904202,0.0,0.0) q[4];
cx q[11],q[4];
u3(1.56804827441798,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.45712196527786,1.63714207503733,-3.13356929636434) q[11];
u3(1.84732882570169,1.10740132668013,3.30995871219856) q[4];
u3(0.845444691091577,-2.26838348909939,2.01184500666224) q[0];
u3(0.505480804786158,1.99931856178067,-3.91033960680031) q[7];
cx q[7],q[0];
u1(1.39954401767447) q[0];
u3(-0.592169301804333,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.17338452702373,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.62467328033704,0.810791752918714,1.61912179173892) q[0];
u3(2.18762200055170,1.83754089359167,-1.05873994906778) q[7];
u3(0.704681649764220,-0.783884034938211,0.964645802902841) q[5];
u3(0.896109710969874,-0.00226661831261249,-0.595610376087897) q[9];
cx q[9],q[5];
u1(2.07762810302832) q[5];
u3(-2.94322813027922,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.20670750860742,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.89523617406842,-2.90019636113570,2.51933035350254) q[5];
u3(0.849123188345313,2.12633501140935,-3.52521874109669) q[9];
u3(1.06847565635009,1.25588333135053,-3.11881365694459) q[1];
u3(0.654368229782084,-2.80327203454753,2.85815876098666) q[10];
cx q[10],q[1];
u1(-1.13333493390038) q[1];
u3(0.489969574348314,0.0,0.0) q[10];
cx q[1],q[10];
u3(3.57421885994278,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.929924812862993,-1.12492592283011,0.245045360341050) q[1];
u3(0.770957095303350,-1.91789036292861,-1.49640664640606) q[10];
u3(1.28644665891084,-0.963837685678299,2.36972483829095) q[11];
u3(0.294229174693653,-1.16812411044568,0.198698286054133) q[2];
cx q[2],q[11];
u1(2.14207277942922) q[11];
u3(-2.57242238468762,0.0,0.0) q[2];
cx q[11],q[2];
u3(0.367110004279197,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.07591011708783,-1.51963176336619,1.62462502914486) q[11];
u3(2.88766929272580,0.987534534421389,2.15521284204949) q[2];
u3(0.930571908397128,1.53869821691123,0.382632357405014) q[8];
u3(0.999521910388648,-0.434407753677507,-2.23618031891962) q[5];
cx q[5],q[8];
u1(3.57637085575237) q[8];
u3(-1.35578988315686,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.95503721090945,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.54817488528594,1.29096320546631,0.870824316888878) q[8];
u3(0.678763538437332,5.35901418089207,-0.695781920777367) q[5];
u3(2.00662551802088,1.46581461941078,-3.78932014911884) q[1];
u3(0.374389070230358,-2.10979781540122,3.83960946602486) q[4];
cx q[4],q[1];
u1(1.47450753697259) q[1];
u3(-0.170852176850061,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.43270595974946,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.83534962259837,1.03301991269425,3.24881132437503) q[1];
u3(1.04807517265032,-0.693655713382916,-3.45682211256692) q[4];
u3(0.457378533237319,-1.97878160351902,1.66480163809623) q[0];
u3(0.601722438741646,-0.0241737799059630,-1.52269089274163) q[3];
cx q[3],q[0];
u1(3.19741218198204) q[0];
u3(-2.24103418764746,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.74152373099624,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.52196109375436,3.36528509151983,-1.22271629935690) q[0];
u3(2.98274012580524,-1.34110615498190,0.246096133752356) q[3];
u3(2.86158750913342,1.53719651362106,-2.90058558685286) q[10];
u3(1.98773391798496,3.20271367003049,-1.60904631800949) q[7];
cx q[7],q[10];
u1(3.45200433058107) q[10];
u3(-4.48162053767802,0.0,0.0) q[7];
cx q[10],q[7];
u3(-0.143191244157398,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.33406099123579,-2.02142116936894,2.68755163349874) q[10];
u3(0.761692559355017,-5.09050017102712,-0.620209314931691) q[7];
u3(2.96683492028598,1.79674801278046,1.09346828964225) q[9];
u3(1.41376156848081,-4.27550650342788,-0.431751949493691) q[6];
cx q[6],q[9];
u1(-0.0745084367623456) q[9];
u3(-1.66438923361508,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.322928656577423,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.99308213432711,-2.26884300692977,3.02705637610012) q[9];
u3(1.38768669213390,-4.97305835936510,0.223130380194741) q[6];
u3(1.15487896426774,-1.33960128154742,0.541628196126872) q[1];
u3(0.957164856171682,-1.66798107177665,0.907585707021276) q[8];
cx q[8],q[1];
u1(1.79772235954440) q[1];
u3(-0.217171182279611,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.76336425192900,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.54454289212511,-0.491004457333120,2.64467577547091) q[1];
u3(1.81256345379241,4.64679182455782,-0.381823765909407) q[8];
u3(1.82790513946450,4.11182697000225,-1.22442615273935) q[2];
u3(0.754359642647761,2.15382133423675,-1.41451700704675) q[7];
cx q[7],q[2];
u1(3.03274126825876) q[2];
u3(-1.70361649716447,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.973546071670076,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.91826929007796,0.0511807572444658,-1.90896044940704) q[2];
u3(0.818531006614829,-0.788244726882387,-1.94144448281657) q[7];
u3(1.34205444947200,-1.49185971781963,1.65693511463668) q[11];
u3(0.560022218831858,-1.91960466907461,-1.08294206546548) q[0];
cx q[0],q[11];
u1(1.67026859814249) q[11];
u3(-2.88709645818389,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.766602385392677,0.0,0.0) q[0];
cx q[0],q[11];
u3(2.23701886119481,-3.13364733728846,1.62793674512247) q[11];
u3(1.76479873761557,-1.47844162502688,4.67799085410926) q[0];
u3(1.11987612317376,-0.476536999637936,2.12431435167216) q[9];
u3(1.90045130970161,-2.37119265566270,-1.45984629433101) q[3];
cx q[3],q[9];
u1(2.56205769698275) q[9];
u3(-2.01841356080236,0.0,0.0) q[3];
cx q[9],q[3];
u3(3.24950526602663,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.329323247168112,0.221979339941738,-2.35357845288063) q[9];
u3(1.35266189445470,4.20689698863247,1.16038676507207) q[3];
u3(2.32602656774611,2.18390609197898,-3.00631141392297) q[5];
u3(1.91858189257004,-2.25330991974997,3.32377545853258) q[4];
cx q[4],q[5];
u1(3.64721981375660) q[5];
u3(-4.37196667073473,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.831034144990585,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.07403070742501,0.964802387755573,-4.36434897242959) q[5];
u3(2.18119455332677,4.58804793642239,-1.67479912374956) q[4];
u3(1.73031517036924,0.992495991296445,-1.19349574157793) q[10];
u3(2.37084648323301,-4.67063450059033,0.195989299661712) q[6];
cx q[6],q[10];
u1(2.57492291495462) q[10];
u3(-1.69903148556809,0.0,0.0) q[6];
cx q[10],q[6];
u3(-0.0745683862904261,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.89261739403980,-1.89002385702519,-0.0897123261777361) q[10];
u3(1.48160757568983,-4.09001036378077,-1.31205714312632) q[6];
u3(1.61699842714364,0.980326003045780,-2.50494240927843) q[4];
u3(1.61977338071374,-4.98403207910735,1.07443484442777) q[10];
cx q[10],q[4];
u1(2.98845104999186) q[4];
u3(-1.90119150424298,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.773462603634419,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.44302393566943,-1.89273957660298,-0.577240003904679) q[4];
u3(1.99801339005437,-2.60641937053732,-1.45270303810373) q[10];
u3(1.93492580870783,-2.54218585877899,0.561021591187451) q[9];
u3(2.68774225952651,2.01270659602719,2.90618755370096) q[5];
cx q[5],q[9];
u1(0.629259369548265) q[9];
u3(-0.484756086322095,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.95499379012010,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.592010529802661,-1.48392193041284,1.06447863814775) q[9];
u3(0.682214651168471,3.24195794469929,1.29105761443878) q[5];
u3(2.63453843326948,-1.60181341596800,1.50331091216734) q[6];
u3(2.82254719161297,-1.60401193943928,-1.55346010758264) q[0];
cx q[0],q[6];
u1(-0.503655517275861) q[6];
u3(-1.78371877411520,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.14049560957267,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.84864111871642,-0.0272052958283191,1.72922280635913) q[6];
u3(1.13556058374169,2.80580201561221,-0.700933891701754) q[0];
u3(1.61041628614746,1.82368706388943,0.911096523738585) q[2];
u3(0.475666360531503,-0.739204053421405,-2.46996911055084) q[8];
cx q[8],q[2];
u1(0.814280105535610) q[2];
u3(-0.541902373806494,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.36386833765827,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.41459889238929,2.25963620681041,-0.721306194721842) q[2];
u3(2.24990162993194,2.65528374758638,-1.47453609518558) q[8];
u3(2.77855109457717,1.38392003128906,-0.628909161980977) q[7];
u3(1.92583074689107,4.85658949612753,0.108902522482220) q[3];
cx q[3],q[7];
u1(2.05210335967642) q[7];
u3(-2.62035012878111,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.17367314092080,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.20998650052203,-2.34182573515376,-0.135865168777415) q[7];
u3(1.23714922462234,-0.758328315308189,-0.174645846270021) q[3];
u3(1.85600065134045,-3.48419403091938,2.44929619209913) q[11];
u3(0.239947031059815,-1.85527674515013,3.69876903908756) q[1];
cx q[1],q[11];
u1(0.769045446007219) q[11];
u3(-1.41129119091234,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.32729938365407,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.34495361647364,-1.27695601901482,0.974260811488931) q[11];
u3(1.29258434267069,-1.87301354155405,-1.58802317878372) q[1];
u3(0.557199221988308,0.979231061418773,-0.450835503656874) q[9];
u3(0.210033362855472,0.139641672254308,-2.78755628370049) q[7];
cx q[7],q[9];
u1(-0.0238322377076237) q[9];
u3(-2.34619788668691,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.19661173429718,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.59030974820855,-1.11930636869867,3.28124878592183) q[9];
u3(1.21680527480964,-3.09464280243515,0.474392625200596) q[7];
u3(2.21585280434262,0.366011165471845,-2.83783262317978) q[2];
u3(2.69076993100816,0.730400403010112,-4.38245024888116) q[11];
cx q[11],q[2];
u1(2.58345041992014) q[2];
u3(-1.82052330471333,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.851506875560445,0.0,0.0) q[11];
cx q[11],q[2];
u3(0.885767699269539,-2.28037175756561,2.81259572161102) q[2];
u3(1.63990953689372,0.658102427081139,4.14924745221725) q[11];
u3(2.01710975853747,-2.13580528107183,-0.860666794931872) q[6];
u3(1.27121280432310,-3.56054582436480,0.102192958556625) q[8];
cx q[8],q[6];
u1(0.642017007207326) q[6];
u3(-1.54131518472249,0.0,0.0) q[8];
cx q[6],q[8];
u3(-0.280626992417230,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.68220242601363,-2.81986449726217,3.41565211723424) q[6];
u3(0.653346183161492,-0.354067021437611,-2.09635594429462) q[8];
u3(1.28333466957662,1.17398278176684,-3.47907438854336) q[4];
u3(0.780635651929042,-2.65439243928435,3.23345957136917) q[0];
cx q[0],q[4];
u1(0.426042970966710) q[4];
u3(-3.10857896107421,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.27536160698152,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.07054913226051,-3.55048668113559,2.14411254579414) q[4];
u3(0.762121098118783,-2.48706473732709,2.04084401305489) q[0];
u3(1.48957689056230,-1.56502993382188,-0.877695491176252) q[10];
u3(1.06581108217719,-3.38304864633169,0.275137145536453) q[1];
cx q[1],q[10];
u1(1.05099186982309) q[10];
u3(-0.0594092956648462,0.0,0.0) q[1];
cx q[10],q[1];
u3(2.67402933267830,0.0,0.0) q[1];
cx q[1],q[10];
u3(0.699207987328256,-4.65609057353958,0.916950283269116) q[10];
u3(1.72496686893825,0.291432631718115,4.73988679144085) q[1];
u3(0.959557618924080,2.36556885113488,-0.801792120792238) q[3];
u3(2.13792762916018,0.571722366229583,-2.71055004189300) q[5];
cx q[5],q[3];
u1(0.934228496056431) q[3];
u3(-0.794976774841288,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.71421685531359,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.91768417636547,-3.38101176343090,1.61933145673946) q[3];
u3(1.93548216532920,-1.29217831044462,-4.82149958215630) q[5];
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
