OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.72694468788570,0.358187951663903,-2.06934130487191) q[4];
u3(2.18977001828971,-3.19199221457417,2.41392168728541) q[2];
cx q[2],q[4];
u1(3.12035069328176) q[4];
u3(-1.81458914854734,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.56410965307457,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.614175335697059,1.56820652986475,1.06262010583081) q[4];
u3(1.94966509827053,-0.552992271795295,-0.863616355821593) q[2];
u3(2.23678318674234,1.69992755463276,0.410193162506458) q[0];
u3(2.22091578350812,-0.980665348686634,-2.80185343645902) q[1];
cx q[1],q[0];
u1(1.87698760070983) q[0];
u3(-2.98777040988437,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.472194449036864,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.754723286884067,-0.390138825814446,1.10941920094514) q[0];
u3(0.532341841099936,-2.57592591509796,-0.493243554992639) q[1];
u3(1.36392752112471,-0.849659705433466,1.80119508560175) q[4];
u3(1.26514030106054,-1.27815599574179,-1.40503694701452) q[3];
cx q[3],q[4];
u1(3.16616388928758) q[4];
u3(-0.761407946123211,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.43072312918792,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.30988280809543,-1.70867583038405,-2.27619646636503) q[4];
u3(2.42819991204856,0.437677923461461,-4.24977083910473) q[3];
u3(2.03304128655675,-1.10366451846366,1.58943607735082) q[2];
u3(1.76871465471836,-1.27899360059084,-1.29800486081515) q[0];
cx q[0],q[2];
u1(1.17536909321492) q[2];
u3(-0.495458667085994,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.91863989661495,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.36353125467483,1.22121918040048,-3.18042895776979) q[2];
u3(1.58434375642610,1.16896842159061,-0.247012073442545) q[0];
u3(2.45210618290452,2.33305778140442,-2.11215918255935) q[2];
u3(1.95071973315631,2.11857055515082,-2.21324169305443) q[4];
cx q[4],q[2];
u1(1.95813285067632) q[2];
u3(-2.91117406271754,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.569201430066193,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.21590068461959,1.40219901518973,0.254507408704496) q[2];
u3(1.80189817039144,0.447004738137838,-4.46840095284046) q[4];
u3(2.23759266376268,-2.14525129438432,0.799435579082915) q[0];
u3(1.83693121675981,-3.11219185471766,0.590808041841355) q[3];
cx q[3],q[0];
u1(1.22578296124293) q[0];
u3(-1.06666891792432,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.00293623604885,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.71463461725198,-1.76867300479754,3.42017340582760) q[0];
u3(2.44454995869098,-1.55239479672443,2.35206327734305) q[3];
u3(3.01620236612434,-2.88663896248968,0.126419288489105) q[4];
u3(2.06486364120004,2.17821479420068,3.45696697353282) q[2];
cx q[2],q[4];
u1(2.31284486427021) q[4];
u3(-1.80628611145814,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.37147043231138,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.62934233437702,1.65069634282109,-2.80608434216602) q[4];
u3(1.70976633029130,1.54372687051899,-2.50950389537537) q[2];
u3(0.435321295846955,2.07817699319501,-2.19577608149686) q[3];
u3(0.146462435818645,1.90598573574174,-3.19062727697926) q[0];
cx q[0],q[3];
u1(2.19186654410481) q[3];
u3(-1.62589736535560,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.71319614779958,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.77541407887716,-1.11575357317968,0.0310011085635694) q[3];
u3(1.52692745908111,1.38359952123428,-2.09309813605062) q[0];
u3(1.35534482799198,1.43017330194209,-0.186215852621172) q[4];
u3(2.38670934701600,-1.14878799122825,-4.07166805982128) q[2];
cx q[2],q[4];
u1(0.953500734875152) q[4];
u3(-3.42308789366492,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.81374203329584,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.18992864320198,-2.09792379839285,0.727376725411788) q[4];
u3(1.67066962138173,-2.18864966834757,-4.07597362588749) q[2];
u3(0.900774820200814,2.45503688449634,0.136759002326484) q[1];
u3(1.59973277138699,-0.0503346668424218,-3.09238951819393) q[3];
cx q[3],q[1];
u1(3.59224759547063) q[1];
u3(-1.17924006390498,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.89009931635011,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.31543001631623,-1.43273107567981,4.60569334602862) q[1];
u3(1.32016637448699,0.312687578980402,-2.79294855847013) q[3];
u3(1.77989305303024,-1.99986520999714,0.629679856949867) q[4];
u3(2.01353369801812,-3.88143304513530,0.148144679099420) q[1];
cx q[1],q[4];
u1(1.03928260055990) q[4];
u3(-0.510261906435253,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.143729300904373,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.22523150955222,-0.863381059143878,4.30207170150191) q[4];
u3(2.39054232917219,1.58122512241314,-0.798954738053114) q[1];
u3(2.16801441511380,0.934402976054122,-2.25262444484999) q[2];
u3(2.43143383703322,1.78583900717275,-3.92552094981920) q[3];
cx q[3],q[2];
u1(1.00375321926358) q[2];
u3(-1.40861808447942,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.526344269627828,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.764047646086926,0.348905102496411,4.19574395544908) q[2];
u3(0.323773196674633,-0.730921639791767,2.80441818413815) q[3];
u3(0.860470619456475,0.334662912485046,-2.71075121683793) q[4];
u3(2.01185164335189,-3.91081057368528,1.75338455497362) q[2];
cx q[2],q[4];
u1(1.91122813334881) q[4];
u3(-3.05979867639579,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.706576278711741,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.655952387730208,2.17802024679732,1.59961605876994) q[4];
u3(1.86776148905860,-0.440908547874343,5.05149538755562) q[2];
u3(2.37003448392626,0.965344926953490,-0.0619484883384890) q[1];
u3(1.91702749578286,-0.469416658647380,-3.48808381840028) q[0];
cx q[0],q[1];
u1(1.48585790209784) q[1];
u3(-0.966949951267680,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.0299413348297850,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.68407643966502,1.77425765631861,-1.51782058660219) q[1];
u3(2.14819419523647,3.43864543403576,0.784711180751890) q[0];
u3(1.76192280190055,-1.18959303303607,1.50048594605212) q[4];
u3(1.93049490577285,-1.34114037706416,-2.26766311032564) q[3];
cx q[3],q[4];
u1(1.50221923076975) q[4];
u3(-2.11781965962951,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.526946363222094,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.23915416809368,2.67858815123219,0.0392991433787557) q[4];
u3(1.36990522285746,-2.39843729908847,3.40989103666732) q[3];
u3(0.903095239950673,1.70213469330057,-0.772517999799510) q[0];
u3(0.611076328649728,-1.22408682268840,-0.454963751272507) q[2];
cx q[2],q[0];
u1(3.12198303647409) q[0];
u3(-1.77455595525873,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.877977539939160,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.419946700633742,0.767797905020747,-2.31520834335700) q[0];
u3(1.81816755324503,4.62728727055675,1.19587556876138) q[2];
u3(0.915228544954019,1.58092964532926,-3.59528074606856) q[4];
u3(1.46645779643679,2.20868152595379,-3.00829573754248) q[2];
cx q[2],q[4];
u1(2.33097422114088) q[4];
u3(-0.0927012919010439,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.54478934322289,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.55671055085450,1.23740081466314,-0.760745018638528) q[4];
u3(1.81573842415688,3.38721420777412,-1.06854376927086) q[2];
u3(0.793721866745921,1.31127920114422,-0.816852106414436) q[3];
u3(1.86691467918824,-0.997928399614342,-4.22392272007671) q[0];
cx q[0],q[3];
u1(1.00108293975687) q[3];
u3(-1.28263133131265,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.202253366575549,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.919647684376001,0.134429397592981,0.752584322988576) q[3];
u3(0.843145500142086,-0.0111032069385398,0.269034144138086) q[0];
u3(2.34816377255574,3.76213539767909,-0.972940756062672) q[1];
u3(2.35795772265126,2.48306269854249,-1.47329619414718) q[2];
cx q[2],q[1];
u1(0.401201786715713) q[1];
u3(-0.849495532047471,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.52984611682607,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.97882763563814,-1.23010337192613,2.68701124110602) q[1];
u3(0.613584716790454,-4.37323162239209,-0.285338654041137) q[2];
u3(0.948071253781475,-0.812450748090017,0.625104187583705) q[4];
u3(0.471375175196360,-2.64808118656872,0.165657699538491) q[0];
cx q[0],q[4];
u1(2.64844710257705) q[4];
u3(-1.39260191511580,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.299769290897989,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.80865958200583,0.148396583947146,-0.720636428092939) q[4];
u3(0.185763417891056,4.16256735453744,-1.88225513760478) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];