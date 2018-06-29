OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.27959994546419,3.19408573367485,-0.856162493249351) q[5];
u3(1.93026161186592,0.238723557489088,-5.48968142873601) q[10];
cx q[10],q[5];
u1(3.01296989444879) q[5];
u3(-1.55704352707807,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.732896571798783,0.0,0.0) q[10];
cx q[10],q[5];
u3(3.05882556669788,-0.690162466066637,4.60916045158875) q[5];
u3(2.31192082095394,3.83193297640488,-1.78444264374661) q[10];
u3(1.79636875385145,0.603698893411974,1.74461502037459) q[1];
u3(1.98021675982626,-1.73386856300552,-1.52549096476061) q[8];
cx q[8],q[1];
u1(2.70368439417413) q[1];
u3(-1.79331570972586,0.0,0.0) q[8];
cx q[1],q[8];
u3(0.108006508932831,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.224613471222055,-2.11952493517079,0.171987184936274) q[1];
u3(0.704471930460305,1.36473095262731,2.93444802876530) q[8];
u3(1.25500049798387,3.30678300578291,-0.589802327243664) q[6];
u3(0.588590959099399,0.679259803168625,-1.09124931488931) q[4];
cx q[4],q[6];
u1(1.61398078003573) q[6];
u3(0.756762717180540,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.17975840411084,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.93789660732523,-2.05314301284997,0.751780916270162) q[6];
u3(1.02065491559320,-2.45946208413834,2.46796803189787) q[4];
u3(1.26507919291245,-0.390601986158353,-1.41715990787298) q[2];
u3(2.46973546067506,-0.0414259453495607,-5.59600959813544) q[0];
cx q[0],q[2];
u1(1.52806773328174) q[2];
u3(-0.0346271104857774,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.13059828135487,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.79313134312002,2.51476207649008,0.697126762709061) q[2];
u3(2.10165227416111,0.356762439919356,0.197084834305946) q[0];
u3(1.77368168513502,1.46182667879389,-0.691534700492426) q[9];
u3(2.56261488334575,-0.479864675107569,-3.86062936240184) q[7];
cx q[7],q[9];
u1(0.602859231492573) q[9];
u3(-1.21315744774421,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.234551161090540,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.45631180138396,-1.77818821138615,3.07776100992251) q[9];
u3(1.70781033154012,0.303838847659621,-5.59190641487948) q[7];
u3(2.35560653286892,-2.00891028035568,-0.194467546739156) q[10];
u3(1.60746498443059,-4.14326285690942,-0.927142067582784) q[2];
cx q[2],q[10];
u1(0.637670930333746) q[10];
u3(-1.36166073905788,0.0,0.0) q[2];
cx q[10],q[2];
u3(2.13231019387539,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.18646906488451,1.75748738852411,-3.60452424280876) q[10];
u3(0.603305486220102,-3.93901934988936,0.634751556123092) q[2];
u3(3.04191910728916,0.330981082111116,-2.82312453226186) q[7];
u3(2.26137417962933,0.980792792219713,-2.48294340280786) q[3];
cx q[3],q[7];
u1(0.0532432558422711) q[7];
u3(-1.01288556129484,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.70077605970884,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.97349705078236,-2.81190163444665,2.06602035548364) q[7];
u3(0.812529024364156,0.111331538165417,-1.30678686486034) q[3];
u3(2.19485611220377,-0.189936091128485,1.30861077297860) q[9];
u3(2.21621177833916,-1.56406247542968,-1.49462188014767) q[8];
cx q[8],q[9];
u1(1.56821830803749) q[9];
u3(-0.705741069573439,0.0,0.0) q[8];
cx q[9],q[8];
u3(3.22660738517475,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.10056674955601,3.43537147564084,-1.86283939538911) q[9];
u3(1.73429988027582,2.32168149822167,-1.86243473472518) q[8];
u3(1.81067586482682,1.66388591464105,-3.15398068930326) q[5];
u3(2.82053430682939,2.59627360578973,-3.54642007772708) q[1];
cx q[1],q[5];
u1(0.325182292639742) q[5];
u3(-0.187646683663186,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.00906675483991,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.17340019679707,3.20379703516997,0.800776436577028) q[5];
u3(1.83480048558720,0.402037497670885,4.96107410736005) q[1];
u3(1.04276497507689,1.99870216891008,-3.36916969892490) q[4];
u3(0.755241890331236,-2.82711769805218,3.40554429120514) q[6];
cx q[6],q[4];
u1(1.22385128600831) q[4];
u3(-1.65171890669292,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.24939394665808,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.81851244351946,4.57245471126384,-1.24776179740369) q[4];
u3(1.43759983573210,0.633910342445922,-1.79627851903904) q[6];
u3(2.58014273664349,-0.181499164277708,2.31140001351202) q[7];
u3(2.01861802454243,-1.39263032852636,-0.558031103482685) q[9];
cx q[9],q[7];
u1(0.566926265876386) q[7];
u3(-1.41491489769514,0.0,0.0) q[9];
cx q[7],q[9];
u3(3.13932681387471,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.61520118518958,0.921615139181186,-2.77390098843654) q[7];
u3(1.28716909152722,1.45714094983046,-3.54064899533190) q[9];
u3(1.50363592485519,3.64486446038963,-2.26869290399585) q[5];
u3(1.56279133906665,2.55985234882646,-2.02933926189330) q[8];
cx q[8],q[5];
u1(3.67613686574589) q[5];
u3(-1.46098246064327,0.0,0.0) q[8];
cx q[5],q[8];
u3(2.09677879337900,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.931273832535940,3.46584519496869,-2.55276800816624) q[5];
u3(0.747554702769916,-2.34885968799938,2.02433068330748) q[8];
u3(0.708953264840999,1.41563004053442,0.717476745531442) q[10];
u3(0.954913186164748,0.254369788986990,-3.64399623696567) q[0];
cx q[0],q[10];
u1(3.21607187118626) q[10];
u3(-1.19389701802523,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.61330905249061,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.71699265086572,-2.72620937969332,2.11031441459591) q[10];
u3(1.47531635111456,-0.573336824040128,2.23127449776550) q[0];
u3(1.69120328493820,1.59761237140661,-2.72918465335966) q[4];
u3(1.28067587375300,1.68122546104321,-2.46316335223291) q[2];
cx q[2],q[4];
u1(-0.675983967803746) q[4];
u3(1.02654775236647,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.96652700769570,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.23311568032948,2.82453781705866,-3.11355847592755) q[4];
u3(2.25007139868531,2.84290231979775,-1.23454695509039) q[2];
u3(2.21815558130620,0.126457963024025,1.21449873150083) q[1];
u3(0.936147073678271,-2.09086397868470,-2.35137549922931) q[3];
cx q[3],q[1];
u1(0.431029459316204) q[1];
u3(-1.71725816736573,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.14365735748339,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.78268144586774,-4.25618721098562,1.45157461538225) q[1];
u3(0.573846270491060,2.01668772483228,0.328803226033541) q[3];
u3(0.847577199172950,-0.597508183580815,-1.18301108988407) q[4];
u3(2.34483159664690,-3.60509884854033,1.54035628568047) q[3];
cx q[3],q[4];
u1(-0.290758215021453) q[4];
u3(-2.25490913102339,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.07881207131252,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.67562978491494,-2.20139694834151,2.51112873322471) q[4];
u3(1.15637386699972,0.782171213235610,-4.46878148284964) q[3];
u3(2.33652869260465,-0.917645104948749,-1.22105440454305) q[2];
u3(1.34111992345247,-5.15374140597037,0.779766088705530) q[5];
cx q[5],q[2];
u1(4.01994854754440) q[2];
u3(-3.52459965311718,0.0,0.0) q[5];
cx q[2],q[5];
u3(-1.05568367663424,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.86165555791401,-0.607069419309139,3.81108597303849) q[2];
u3(2.06978648225366,1.32899748982485,0.214667844722039) q[5];
u3(1.28182063881241,2.38026880674593,-0.372678503710528) q[8];
u3(1.27907156965253,0.701867133454764,-3.91600584134337) q[0];
cx q[0],q[8];
u1(1.64003323334849) q[8];
u3(-0.572519866016258,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.75935640352611,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.21949742094951,4.05462588538438,-2.20252875481121) q[8];
u3(0.411899901971756,0.459545416704613,-4.00514302616064) q[0];
u3(2.15927505511087,0.0425737569297817,2.49986130068962) q[7];
u3(2.24989774836987,-0.858972154238689,-1.49259951081993) q[9];
cx q[9],q[7];
u1(0.126515481782146) q[7];
u3(-0.462021225683654,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.61092522963180,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.69681136407686,4.03043836254037,-1.44124953769185) q[7];
u3(0.705635659071491,-3.18997906605459,0.668817127399712) q[9];
u3(1.53421006664753,-0.109382300273391,0.683725550167556) q[10];
u3(1.51420386562888,-2.79409915295244,-1.24706920517856) q[6];
cx q[6],q[10];
u1(0.951683409291852) q[10];
u3(-1.39323779825309,0.0,0.0) q[6];
cx q[10],q[6];
u3(-0.194988967185977,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.58306266856712,0.197775437224389,1.64890981283609) q[10];
u3(1.31139695369964,1.01180094885572,-2.53162369513634) q[6];
u3(1.43837393771519,0.619703580370827,0.540104797129180) q[9];
u3(1.40546897855474,-1.58120495461363,-1.11790108836089) q[1];
cx q[1],q[9];
u1(1.15087007977935) q[9];
u3(-0.547809013639775,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.48364765306165,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.08196953318521,0.296418249932040,0.454914618094666) q[9];
u3(0.916033586156399,0.208277755803877,1.91005848579198) q[1];
u3(1.96326460010614,-1.75621294220428,1.76753820373916) q[2];
u3(2.63074308818132,-1.50068290358275,1.11346745070233) q[4];
cx q[4],q[2];
u1(2.98337753224806) q[2];
u3(-2.12962424835054,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.569099971069448,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.96090845840059,3.63381932522668,-0.0878634331689405) q[2];
u3(1.11951711601514,2.27656929027242,-2.86933353471320) q[4];
u3(0.446719526655008,0.439904235977224,1.42582885628848) q[10];
u3(1.71960373964023,-1.55360299332229,-0.303589505003820) q[5];
cx q[5],q[10];
u1(-0.199435672953762) q[10];
u3(-1.65370630688758,0.0,0.0) q[5];
cx q[10],q[5];
u3(0.931860359768356,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.65036024301523,-0.321817990960391,-1.13832976544843) q[10];
u3(2.21501102514870,1.47147198578932,-1.60363153074488) q[5];
u3(2.80566954268542,1.98304454317817,-4.10813098823600) q[8];
u3(1.61989310224901,3.49799069839717,-2.57846897210341) q[7];
cx q[7],q[8];
u1(0.962130068886668) q[8];
u3(-1.51515715695460,0.0,0.0) q[7];
cx q[8],q[7];
u3(-0.470851336318077,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.04042680125264,-2.61604738322600,2.65268366362083) q[8];
u3(2.29687342701382,2.71376630863753,1.06314989309057) q[7];
u3(1.03349220671162,1.99921102248844,-2.84034846902811) q[6];
u3(1.12332333507118,2.38318814838570,-3.56487641938849) q[0];
cx q[0],q[6];
u1(1.67669301429221) q[6];
u3(0.274273150268754,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.661665496540955,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.127281536748764,-1.93420109587014,1.77433069144680) q[6];
u3(0.505673114162230,-0.0169638168256456,-3.32973065767018) q[0];
u3(1.28337297153344,1.35692815405253,1.19687005364523) q[9];
u3(1.45525531136511,-1.36716100467945,-2.52388308771850) q[10];
cx q[10],q[9];
u1(0.0572724707962313) q[9];
u3(-1.72160083362040,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.08657517707149,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.80557693049639,-1.79633882036610,3.91239500732300) q[9];
u3(1.49860663897474,5.83028618265471,-0.195761718815580) q[10];
u3(1.49383837372751,-0.238308676102540,0.889946003725877) q[2];
u3(2.68896567757781,-0.728588969152562,-1.93586421146123) q[0];
cx q[0],q[2];
u1(1.79750391677009) q[2];
u3(-3.09034929268917,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.311073137696293,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.03924794873878,-0.508884132597007,2.72075435636615) q[2];
u3(1.29400288918000,-0.632436937582744,-4.12933486165801) q[0];
u3(0.850743985097741,0.193034107608799,0.483209411742747) q[1];
u3(1.51827100263935,-0.479714389763007,-0.919132283958314) q[3];
cx q[3],q[1];
u1(0.795968506673507) q[1];
u3(0.0153878599960890,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.96849616912088,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.661259050539735,1.86220945574238,-0.914664622391279) q[1];
u3(1.27171550652287,4.03086124378607,-0.502756227911989) q[3];
u3(0.896636688559680,1.80142319718007,-2.03310604175669) q[4];
u3(0.464793462893951,1.16082626099561,-2.89202625144124) q[7];
cx q[7],q[4];
u1(0.523272846076203) q[4];
u3(-0.804186015109993,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.23967951841854,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.90328077531097,-1.95567231294408,-0.0184739415957880) q[4];
u3(1.56320973048908,2.24155207357876,3.45944674817682) q[7];
u3(2.50988268597451,2.32593083237405,-1.56481105876306) q[8];
u3(3.04343417325188,5.34780168969993,0.899793600994764) q[6];
cx q[6],q[8];
u1(2.64327232736657) q[8];
u3(-1.92108698628599,0.0,0.0) q[6];
cx q[8],q[6];
u3(0.688499533103328,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.75126686225918,1.92620959745979,-2.40531239348500) q[8];
u3(1.26198726769631,-1.52081619458984,-0.785982216223962) q[6];
u3(1.95023112325020,3.17009659525606,-0.304810897587285) q[0];
u3(1.58932357310587,2.29145074188075,-1.42681050115682) q[8];
cx q[8],q[0];
u1(3.82909922016046) q[0];
u3(-3.39084297497273,0.0,0.0) q[8];
cx q[0],q[8];
u3(-0.528368443384568,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.53388985739321,1.04715425343122,-1.36726419409703) q[0];
u3(0.258564063442721,3.89796162667204,-2.11175431749960) q[8];
u3(1.43737546003177,0.734467131708875,-1.24484633666141) q[7];
u3(0.786712894939523,-4.15082215182995,1.79854560791680) q[10];
cx q[10],q[7];
u1(-0.556805883346708) q[7];
u3(-2.04426589375545,0.0,0.0) q[10];
cx q[7],q[10];
u3(1.29059144026280,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.67285696605645,0.935596144183989,-2.23874087233009) q[7];
u3(2.41246713320457,-5.55698352842718,-0.608379293235233) q[10];
u3(0.357874202822039,-1.62719621934770,2.42046459710862) q[6];
u3(0.763197280324096,2.39900870934779,-3.30833066478473) q[9];
cx q[9],q[6];
u1(0.915658166619626) q[6];
u3(-1.65200238538968,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.84877691702647,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.99958138712158,-0.902149006026606,1.57227547963862) q[6];
u3(1.26304530373854,1.13764832421865,-3.30905377872633) q[9];
u3(0.979291402306870,-3.41148652382014,2.09449543164496) q[4];
u3(1.68771349082475,3.14440048707549,-2.58123646267845) q[3];
cx q[3],q[4];
u1(1.68253793254634) q[4];
u3(-2.18993082045670,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.88740093103961,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.95140204436153,2.57623535331409,-3.50339009994612) q[4];
u3(0.336840857585726,2.04047607906615,1.38630211474698) q[3];
u3(0.504254787571708,1.33302537869395,-3.96877886762971) q[2];
u3(2.58797551159367,-1.47433478193740,4.28817961298039) q[1];
cx q[1],q[2];
u1(-0.335167834599486) q[2];
u3(-1.82523206035763,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.19157712974946,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.538633441784355,2.94987673789457,0.435435739430182) q[2];
u3(1.81640299088837,1.02907069146953,-1.27808731484476) q[1];
u3(2.24193653188719,0.916370030661628,1.57704112584276) q[9];
u3(1.85133622713055,-2.29021593523953,-2.84764607643640) q[1];
cx q[1],q[9];
u1(1.87506993607512) q[9];
u3(-3.07735216070562,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.448297117643279,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.44237350248574,-3.89962147359550,1.93024367168073) q[9];
u3(2.90600426837052,3.84130710708769,-0.0222263247202243) q[1];
u3(2.22743795104017,0.000840628778041586,2.72636123141057) q[0];
u3(2.27239021251217,-1.87464299957358,-2.03195392153324) q[8];
cx q[8],q[0];
u1(1.64752682637306) q[0];
u3(-2.66745464075628,0.0,0.0) q[8];
cx q[0],q[8];
u3(3.22046989910539,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.30849907909999,0.979735520864154,0.283107814514242) q[0];
u3(0.787586544440852,-3.14399642017162,-0.389740838258754) q[8];
u3(1.48460682108178,1.15354351476937,-3.63116355608233) q[6];
u3(1.56278001126712,2.95956128326144,-2.63135467800729) q[4];
cx q[4],q[6];
u1(-0.162122432889249) q[6];
u3(-1.77156480751075,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.32213329474909,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.516877863719843,0.238659169546868,-1.62360591734682) q[6];
u3(0.999215589015521,-4.30423349572979,-0.587095823930442) q[4];
u3(0.271491745459192,0.476153113819620,-0.408301565277138) q[7];
u3(0.834766699097991,0.245687041379119,-2.65916994615888) q[5];
cx q[5],q[7];
u1(1.54152984367685) q[7];
u3(-2.32543884170592,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.33389685162140,0.0,0.0) q[5];
cx q[5],q[7];
u3(3.00302488158568,2.32456827141141,-2.76569756140911) q[7];
u3(1.57241171313393,-1.64815966838392,3.04321850348971) q[5];
u3(1.72268372189833,0.773580838700070,-2.16525082688720) q[3];
u3(1.96261581352621,-4.13384712922118,1.57496178719571) q[10];
cx q[10],q[3];
u1(2.43486498555771) q[3];
u3(-1.66810840145949,0.0,0.0) q[10];
cx q[3],q[10];
u3(3.45979089951295,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.70231541905837,4.35309817269876,-0.679764040502090) q[3];
u3(0.786445004498686,-1.74492146284047,2.45982224044005) q[10];
u3(2.39133758707209,2.60421648482692,-2.12701576993945) q[1];
u3(1.85407156295203,2.15710828223092,-3.68288703313579) q[6];
cx q[6],q[1];
u1(-0.250082825871534) q[1];
u3(0.794001431826405,0.0,0.0) q[6];
cx q[1],q[6];
u3(4.10236366470584,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.54484241226696,1.57554606015722,-0.602720160961087) q[1];
u3(0.712236863302795,-2.24864138961296,2.21330412597232) q[6];
u3(0.731544474262626,-3.07122232798607,2.75555668085552) q[3];
u3(1.02902041884687,1.84063888164022,-4.09414269540978) q[4];
cx q[4],q[3];
u1(2.90231088740103) q[3];
u3(-1.82546665453501,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.840238452680204,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.13463800751231,0.585007261985890,1.30656455783376) q[3];
u3(1.75228738862664,-0.742588587307748,-3.80583463981315) q[4];
u3(1.54924709834096,2.12155054156885,-2.49524004534075) q[0];
u3(0.330166083237442,2.68987987641312,-3.02538158008256) q[8];
cx q[8],q[0];
u1(0.523733448052420) q[0];
u3(-1.41623433460416,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.91098932949247,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.38554907827326,3.11038978152156,-1.88819082427352) q[0];
u3(1.31052397035109,-5.31518702972600,-0.487665851952419) q[8];
u3(1.18421299696605,0.831176568995717,-1.86207459620908) q[9];
u3(1.67282421667166,-4.10386890430239,1.74897890126226) q[7];
cx q[7],q[9];
u1(3.11822032075754) q[9];
u3(-1.86704520591678,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.70642701721331,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.54249348339102,0.420373072077710,-0.420892409873727) q[9];
u3(0.998308490082716,0.269296109984008,5.65216680876265) q[7];
u3(1.20561405651335,-0.617526448253875,1.87109345021281) q[5];
u3(1.46484277795715,-0.683607843226880,-2.10606614325009) q[2];
cx q[2],q[5];
u1(1.48961850195285) q[5];
u3(-2.35757705349862,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.11345993532756,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.75405393169991,3.19306535807013,-1.45446567869213) q[5];
u3(2.60820524431541,4.69439216538343,-0.244745707409871) q[2];
u3(1.61338471038491,-3.51446578595951,2.50026128519854) q[2];
u3(0.186095353135626,-0.371706327291846,2.10681009294958) q[8];
cx q[8],q[2];
u1(0.646196389272276) q[2];
u3(-1.26280215975827,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.0471241960545761,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.60805361173337,2.23726711630008,-2.77440982943322) q[2];
u3(2.56409952884596,-1.79896668585029,-4.02436810672677) q[8];
u3(0.813899861202219,1.50111919269198,-2.87689663394222) q[10];
u3(1.21139869545673,-2.09385753031978,3.22646321748242) q[9];
cx q[9],q[10];
u1(1.65644869049531) q[10];
u3(-2.48715231443522,0.0,0.0) q[9];
cx q[10],q[9];
u3(3.18659702500037,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.46446711052722,0.886718321755770,-2.02840247931378) q[10];
u3(2.29938419326787,-4.36443797163465,0.192403215630604) q[9];
u3(0.254212484804081,1.14347163882869,-1.86513793223086) q[4];
u3(0.497081231148721,-0.397531378154548,-1.31282872962642) q[1];
cx q[1],q[4];
u1(2.94199091412267) q[4];
u3(-1.32453340716174,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.00768620379062,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.63880174023777,-3.85044996560757,1.82979508235112) q[4];
u3(2.37490793833015,3.71853717304947,-1.54200112828809) q[1];
u3(2.13454146424852,2.83722768808400,-2.45172771008013) q[7];
u3(1.35649956372359,-2.83055215519641,2.57867366248669) q[0];
cx q[0],q[7];
u1(2.27253414959683) q[7];
u3(-3.19494083324135,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.54123806684983,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.47608208784533,-0.0863984925293444,-0.774836802014803) q[7];
u3(1.35498172319777,2.03862761956475,1.61878387228794) q[0];
u3(1.14142738463817,-0.511502070364685,-0.111453630723172) q[5];
u3(1.09398427367662,-2.53259357806582,-0.263666317527557) q[6];
cx q[6],q[5];
u1(0.396401929131049) q[5];
u3(-1.12900575505670,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.59229268021238,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.44541157176925,-2.38400064323790,-0.100108417121950) q[5];
u3(2.83613272575304,1.39374942868247,2.91580648325095) q[6];
u3(1.76995022811229,-3.31793497973255,0.328979380005707) q[4];
u3(1.55182675016790,0.0370948087794003,3.56509675937802) q[6];
cx q[6],q[4];
u1(3.11343639428008) q[4];
u3(-2.46837333948182,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.868976234203099,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.436194746199743,2.36535525016058,-2.93069124974534) q[4];
u3(2.53234265007837,0.763550361887775,-0.458295061861997) q[6];
u3(1.68998014016141,0.542336336597469,-0.786970686334937) q[7];
u3(2.25008914882489,0.721553329729983,-5.09788668476504) q[9];
cx q[9],q[7];
u1(1.24456077379243) q[7];
u3(-2.96835660821763,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.01914058780296,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.838207114230113,1.89388565855562,0.0633696178458645) q[7];
u3(1.45481305422078,4.74375375431975,-0.210091542045515) q[9];
u3(1.18238391856381,0.723682619994939,-2.32079469352995) q[0];
u3(1.59552049978697,-3.78971265821557,1.84633353143929) q[1];
cx q[1],q[0];
u1(0.715089453577028) q[0];
u3(-1.66340499680908,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.99531551670538,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.574309531039367,-2.06849198288646,3.23890737967232) q[0];
u3(1.85257705972034,2.99737225764253,-0.355750504884520) q[1];
u3(2.06515573158313,-4.23106431697111,1.95289637064523) q[8];
u3(0.479819970102022,1.33837052316866,-0.0257768380104542) q[10];
cx q[10],q[8];
u1(0.278947221756997) q[8];
u3(-1.22527793351552,0.0,0.0) q[10];
cx q[8],q[10];
u3(1.65728183555760,0.0,0.0) q[10];
cx q[10],q[8];
u3(2.36094105482785,1.53395635011447,-2.23630504374701) q[8];
u3(2.20348374709843,-1.46008414378789,2.45734084502838) q[10];
u3(2.13061253518047,0.472702402255791,-0.115190221985659) q[5];
u3(0.886624689044174,-3.00276753747969,-1.15577097169022) q[3];
cx q[3],q[5];
u1(0.736329388725762) q[5];
u3(-3.39105472576843,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.78837341726063,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.658878209663977,-1.23550909022439,2.20784647678899) q[5];
u3(0.869420880359771,-2.59573198145330,-2.10630322583036) q[3];
u3(1.55013125623554,-0.908200202912459,-1.27603431160952) q[9];
u3(1.60983072578621,-3.27605463332844,0.165261457350915) q[7];
cx q[7],q[9];
u1(2.89977125971484) q[9];
u3(-1.89626285643452,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.00759140554755,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.56753846278167,-2.51840936919720,2.02813635815307) q[9];
u3(2.74866386895233,0.547900924628432,-3.53734394526768) q[7];
u3(2.44520662344211,3.14739144587157,-2.43064888394854) q[10];
u3(2.18767402007182,1.77064430804421,-1.35816860031151) q[3];
cx q[3],q[10];
u1(0.238196762664806) q[10];
u3(-0.954319448870063,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.11054621658715,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.37658391920988,0.443069294025343,2.09355975175562) q[10];
u3(0.411047278553614,0.131222783341737,1.63736005256594) q[3];
u3(2.03826907654958,-0.282548096722668,1.70983106665356) q[4];
u3(1.96612548912272,-1.05749821913046,-1.87509927983004) q[2];
cx q[2],q[4];
u1(-0.491496922241883) q[4];
u3(1.25068813492764,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.54493041682268,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.47300051566690,-0.0963698386771197,0.274708836314677) q[4];
u3(2.46069666177779,-0.994500267692038,2.63106739145328) q[2];
u3(1.75538638033486,1.15358256373557,-1.55421469524885) q[5];
u3(1.18798594694752,1.10114328633748,-4.00909527643801) q[1];
cx q[1],q[5];
u1(0.285898754979104) q[5];
u3(-0.972126064649005,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.48115507553683,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.42901848728161,-0.0690913954431751,2.21371739588585) q[5];
u3(1.83712885332180,1.02613284463173,2.97156972054058) q[1];
u3(2.62977416712781,0.0315042856226561,-1.23313168944805) q[6];
u3(1.60443581580564,-4.16918623533326,0.942460807179459) q[8];
cx q[8],q[6];
u1(-0.502243032321484) q[6];
u3(-1.72877661173302,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.924975810839569,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.27758142480196,3.19593534268578,-1.64137947035106) q[6];
u3(2.06187756758842,0.788515113626127,0.909379164200391) q[8];
u3(1.34136519906689,1.28023338053445,1.02526159172709) q[8];
u3(2.40233382799573,-1.77635042717690,-0.538684260091219) q[6];
cx q[6],q[8];
u1(0.0200518316286302) q[8];
u3(-1.75045024638777,0.0,0.0) q[6];
cx q[8],q[6];
u3(0.364558331430245,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.662465217673950,-1.29975776357645,2.06795128151951) q[8];
u3(0.540510429915447,0.0304571901436761,3.47617013478956) q[6];
u3(2.28790484981881,2.48628691841757,0.293035221643458) q[10];
u3(1.94812415335719,-0.992470573593168,-5.22806054314674) q[4];
cx q[4],q[10];
u1(1.95990815575353) q[10];
u3(0.383749205046709,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.956548244245860,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.87001658824416,-1.09103103616167,3.74758526317867) q[10];
u3(1.44339972134292,3.39635467994533,-1.34019494404942) q[4];
u3(1.67370838350205,1.36032303205558,1.13158590583138) q[1];
u3(0.836738321029219,0.307065984535007,-2.77878106975543) q[3];
cx q[3],q[1];
u1(1.27913286606566) q[1];
u3(-0.486267903015111,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.02772233617955,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.69904995934488,1.44565207790279,-0.546743632389734) q[1];
u3(2.25969653163135,2.48211078599500,2.83131932866658) q[3];
u3(1.86425349456830,1.06648904826846,-0.149929903855853) q[2];
u3(1.41781432852492,-0.288351128443622,-4.27680605325069) q[7];
cx q[7],q[2];
u1(-0.671541826197361) q[2];
u3(0.344524416996875,0.0,0.0) q[7];
cx q[2],q[7];
u3(4.37328404448146,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.76519726584113,-3.83132276514175,1.98262945257515) q[2];
u3(1.86148712528528,4.11558082050043,-1.83553271810898) q[7];
u3(1.61013424862750,3.31759315772347,-2.58135677633156) q[0];
u3(0.866035361364608,2.25497739264721,-1.55381321948118) q[9];
cx q[9],q[0];
u1(0.360877773620552) q[0];
u3(-1.28298031153101,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.18127929377335,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.79007196925294,0.909450295591818,-4.09868009627343) q[0];
u3(2.69617495761861,-0.770435063752371,-3.73324114660618) q[9];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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
