OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.08879677056825,-1.67894159057380,0.431944009912680) q[2];
u3(1.28153024495171,-2.30689047589873,-0.618100103902147) q[3];
cx q[3],q[2];
u1(4.45090386946469) q[2];
u3(-3.55984450889220,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.665318541209946,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.13258052836657,1.51000795097673,-3.70562082127348) q[2];
u3(2.26865744582546,0.370559956981862,2.35426744958733) q[3];
u3(2.29054122899552,1.56362116368973,1.13458398446013) q[0];
u3(0.665282968016713,0.164631404686079,-3.95132074445093) q[1];
cx q[1],q[0];
u1(4.49763100064130) q[0];
u3(-3.66368939239899,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.420120028421296,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.313956778340352,-0.0735116215439746,-0.719069704918734) q[0];
u3(1.69335071569321,2.25151163489378,-0.0223167268725604) q[1];
u3(1.98598840679440,0.694630878438073,-2.30551593322831) q[1];
u3(2.42879086700918,2.47959960260358,-3.17339049002768) q[0];
cx q[0],q[1];
u1(3.24931344220932) q[1];
u3(-1.43653550657710,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.54703696706392,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.71693625598893,-1.25279538472538,0.809122733655682) q[1];
u3(0.802125319623357,-3.30963728435565,0.784554602392026) q[0];
u3(1.75768148294220,1.25209153910603,-3.38173078858461) q[3];
u3(1.76396802056069,-2.74737121334550,3.22409045589626) q[2];
cx q[2],q[3];
u1(0.636805472224401) q[3];
u3(-1.92527528402563,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.78262465361696,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.91388869711592,0.830489597081048,-3.00641563402899) q[3];
u3(2.03734389348677,-2.64194092175202,-0.371328478428709) q[2];
u3(1.78569797141116,-2.58044600861950,3.29469318762910) q[3];
u3(0.221897277503301,0.770132593934062,1.03824516462641) q[2];
cx q[2],q[3];
u1(4.22450545034269) q[3];
u3(-3.56365724341958,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.792619791993970,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.93022845130012,3.47019735734294,-0.0790603165386712) q[3];
u3(2.84181285460818,-1.10717743551391,-2.85430968629629) q[2];
u3(2.04195802830309,2.84685860323092,-3.07269579900557) q[1];
u3(1.44685326662423,2.67839848710374,-1.99056956034655) q[0];
cx q[0],q[1];
u1(1.63495474357062) q[1];
u3(-2.48308755994327,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.35405601004859,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.781637939714191,-1.05901153651323,-1.25101446625385) q[1];
u3(2.01946254328668,-4.14600178998696,0.643013967225811) q[0];
u3(1.62822996876164,3.08288851312193,-1.42029873291657) q[1];
u3(1.59742482787320,0.0872594937835987,-2.94952796489201) q[3];
cx q[3],q[1];
u1(1.42712207144874) q[1];
u3(-0.0811193981561544,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.31204726588201,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.61508293806984,-0.485374077683451,-1.68826284860627) q[1];
u3(0.433398353045437,1.77927685184579,3.48489582008769) q[3];
u3(1.48224524622875,-0.433270870558503,0.735499262328348) q[2];
u3(1.79643607276870,-0.821408598007739,-2.01855994272189) q[0];
cx q[0],q[2];
u1(4.18789505845458) q[2];
u3(-2.89517261630945,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.0428224665408186,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.12952665306440,0.501111102940958,0.127448194856131) q[2];
u3(1.99254242648583,1.25937865294895,3.19899707774384) q[0];
u3(1.73747248408818,0.0223207071124862,0.639978575800065) q[1];
u3(1.65666533484725,-0.355941979323976,-1.92785957644873) q[0];
cx q[0],q[1];
u1(1.13230970912442) q[1];
u3(-0.673678381408524,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.78848074609648,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.349391995326837,-2.37506502982667,2.08195454204593) q[1];
u3(2.43745781685567,-1.02555137955505,-3.32561426541559) q[0];
u3(0.246398195111556,-2.58371953593518,3.29079476135408) q[2];
u3(0.995669270772660,-2.87317364311150,1.47138073060164) q[3];
cx q[3],q[2];
u1(1.60005459910146) q[2];
u3(-3.14061767137365,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.18822304270357,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.25341152111559,-0.535343558484873,-1.49387195311895) q[2];
u3(0.990740317116169,2.62814013286792,2.26136178611572) q[3];
u3(2.59132831739939,0.574002266307539,-2.12272384978494) q[0];
u3(2.47480062195076,1.81538995073120,-3.85838966145629) q[3];
cx q[3],q[0];
u1(0.887902697812093) q[0];
u3(-3.08669353339736,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.88986689703907,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.34231784381814,-0.0527668249335557,0.0960673177013246) q[0];
u3(2.93471528401792,-5.36167244321439,0.258284382663159) q[3];
u3(0.585235123487515,1.12975881207100,-3.65473042818748) q[2];
u3(1.44497208911388,-3.17758047988848,2.50751094712269) q[1];
cx q[1],q[2];
u1(2.53443668053750) q[2];
u3(-1.63280379500501,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.38465086696293,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.59795842024590,-1.95584159266683,-1.21155549776208) q[2];
u3(2.67834838915476,1.08942794363058,-1.31794352249062) q[1];
u3(1.06125163730402,2.27015174859622,-2.01395718923338) q[0];
u3(0.125991265137919,1.52664983074429,-2.00991976997385) q[3];
cx q[3],q[0];
u1(3.70155486650637) q[0];
u3(-3.36315896270602,0.0,0.0) q[3];
cx q[0],q[3];
u3(-1.02950251433373,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.336213109689448,-0.658256508781157,0.132908150319246) q[0];
u3(1.64264116153019,0.594602829740499,-1.68519164998359) q[3];
u3(2.67435897479079,0.961461710317557,-2.98896323359528) q[2];
u3(1.88410991788348,-2.67601973854007,2.59177728215277) q[1];
cx q[1],q[2];
u1(0.566943093643980) q[2];
u3(-1.68873722964672,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.361147885515884,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.13339301362959,-2.48282251693093,2.16371638770356) q[2];
u3(2.21321844876988,0.990772980704299,2.19828452171996) q[1];
u3(1.52684566448683,-2.12281083728901,-0.616274366766248) q[2];
u3(2.33510166550881,-2.31875190754023,-0.535563886181791) q[3];
cx q[3],q[2];
u1(2.28709694667929) q[2];
u3(-1.98600415370751,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.09884152521872,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.34926842652915,-0.0221546850687462,-3.21044096951970) q[2];
u3(2.10663455818080,0.556535777748120,-2.79124214341358) q[3];
u3(1.58520066390696,-1.42433318704883,-0.308454731287373) q[0];
u3(1.63167369939328,-4.21041008402609,-1.19471876400383) q[1];
cx q[1],q[0];
u1(1.07537979020302) q[0];
u3(-0.502696481842113,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.17792087622911,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.646685833660388,-0.122343798789013,-0.751532827538013) q[0];
u3(2.79038822195622,-0.302629353724567,-1.01188352965231) q[1];
u3(1.08959017214425,1.68650943889736,1.11818301345579) q[2];
u3(1.77357740665623,-1.60444838789365,-1.39817845331935) q[1];
cx q[1],q[2];
u1(2.63578207711442) q[2];
u3(-1.73633751167931,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.233419196893207,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.84407874498176,-3.16063644381956,2.06728207377102) q[2];
u3(1.56465589137056,-4.06181004509413,1.33781278605934) q[1];
u3(0.941577610323314,1.47503069090618,-0.508842310466331) q[0];
u3(1.35955888091966,0.851161230764900,-3.77794835265547) q[3];
cx q[3],q[0];
u1(2.64864400768705) q[0];
u3(-2.03074962124385,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.38999736337553,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.25723805192851,-4.27145273847700,1.64311341306989) q[0];
u3(0.852018930909798,0.103093253637721,-2.09577725018768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
