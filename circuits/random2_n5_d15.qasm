OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.09162069120059,1.58206854522029,-3.10416141606649) q[4];
u3(1.44365145727943,2.37794850627079,-2.47308880140517) q[3];
cx q[3],q[4];
u1(0.310647003342824) q[4];
u3(-1.48411568123117,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.02644648192610,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.50396464273996,-0.572409156869180,0.436404908198563) q[4];
u3(2.20722829898556,1.47735963928339,-4.36413590626195) q[3];
u3(2.58197518607449,-2.42070039521149,2.03325090852648) q[1];
u3(2.58600122725964,-1.50970718328834,0.0813575060963023) q[2];
cx q[2],q[1];
u1(0.304832403971412) q[1];
u3(-1.56970835546530,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.90890955589399,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.47366660989626,-1.99049798769277,3.79785698837161) q[1];
u3(1.78179850571271,-3.98509573611259,-1.31586456447099) q[2];
u3(0.793387064248973,-0.677050714650782,0.629274478630220) q[4];
u3(0.763001927466926,-3.05586816577500,2.30096038181276) q[3];
cx q[3],q[4];
u1(3.12137862586155) q[4];
u3(-1.10197012159706,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.58741321183699,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.82941256480254,-1.32831045822355,-1.88682765650367) q[4];
u3(0.762729931042676,-0.738431695682222,2.58642330160432) q[3];
u3(2.76016413164475,0.384346848816294,-1.91818643858959) q[0];
u3(1.97874149719848,3.58140765683382,-0.347374297982514) q[2];
cx q[2],q[0];
u1(1.77205234598441) q[0];
u3(-2.32149556946506,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.385668281819762,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.29594340055203,0.593488898703083,-0.0746700355455128) q[0];
u3(2.49789815976643,0.902718355497138,3.69024640241280) q[2];
u3(1.42914214958014,2.28735292740478,-2.84831796532113) q[3];
u3(1.01984999476870,0.976809535513148,-1.70043955346690) q[4];
cx q[4],q[3];
u1(1.35861921402474) q[3];
u3(-0.119382784517658,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.58870832209627,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.812903104301281,-3.02429751500234,0.293542923588740) q[3];
u3(1.11739779719386,-0.199063107224061,-4.72356035046796) q[4];
u3(1.54499747178047,1.52724847681112,-3.05687330888831) q[0];
u3(2.93985264880467,2.22704523190394,-3.87451479849979) q[2];
cx q[2],q[0];
u1(0.576288445505085) q[0];
u3(-0.0507566934380768,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.22906144790868,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.95104842362870,0.508576945303662,0.104878280332539) q[0];
u3(1.82270829210905,3.05368048600697,-0.486556687703123) q[2];
u3(2.87659000886041,1.35691954295504,1.68251913369503) q[3];
u3(1.29300367977510,-0.0812787381971880,-5.36992302172499) q[2];
cx q[2],q[3];
u1(1.68007531311244) q[3];
u3(-0.875659896067558,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.688728400810393,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.32681664557586,-3.74473155532452,2.06054767070482) q[3];
u3(0.881689244557280,-2.80216225088081,-1.49615298869577) q[2];
u3(1.73759644770259,0.767495269731478,1.65101793695422) q[1];
u3(1.72282351375514,-1.39036816470214,-1.86777024966470) q[0];
cx q[0],q[1];
u1(1.67926538423248) q[1];
u3(-2.95085878959807,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.919633443121968,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.922152328917081,-2.11401172560533,1.27722432521579) q[1];
u3(2.32160386517169,1.76874506045201,-4.37451038761029) q[0];
u3(2.25614967185293,0.217452806212996,-1.17559532650454) q[1];
u3(1.58793369648505,0.794282108902339,-4.41112766205235) q[0];
cx q[0],q[1];
u1(1.17952807096583) q[1];
u3(-0.660778602913229,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.91632126118340,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.11260345973541,-2.67448815534030,3.55764767194036) q[1];
u3(1.37394756723997,1.00279810081206,3.95490178981010) q[0];
u3(1.22930737876985,-1.02699087095292,1.21723323281822) q[3];
u3(1.15169829641415,-1.90505862111838,-0.864007343463319) q[2];
cx q[2],q[3];
u1(0.00611195885032245) q[3];
u3(-0.660873613208804,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.02229892830801,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.110197925925321,-1.88632563512809,3.58244818497398) q[3];
u3(2.62250653640369,-4.63959323165705,-0.480099881590637) q[2];
u3(0.675552165714814,-1.93796058741904,2.32397726315969) q[2];
u3(0.483984591686241,-3.01676663638164,1.16546740037263) q[3];
cx q[3],q[2];
u1(1.46884635439122) q[2];
u3(-2.64565858080099,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.738662963423098,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.38978166488535,0.662619395801987,-4.80352078264041) q[2];
u3(2.03759856786940,-3.12029063464672,-2.40358064002391) q[3];
u3(1.50529205414543,2.20512730464783,-3.25374517409567) q[0];
u3(1.53287557542419,2.71763541203941,-3.13646836760276) q[4];
cx q[4],q[0];
u1(1.65479250335391) q[0];
u3(-1.09920393964380,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.81791746889790,0.0,0.0) q[4];
cx q[4],q[0];
u3(3.05410285053078,-1.57212091011780,-2.81424163056721) q[0];
u3(1.95897215232921,-3.97519674395115,-0.320046676335343) q[4];
u3(1.55110470196977,1.75671232072367,-2.86669866094088) q[1];
u3(2.31863606000297,1.82365014360223,-3.84661585414071) q[3];
cx q[3],q[1];
u1(1.15713054293417) q[1];
u3(-0.861681601251826,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.96139849357279,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.35793342031571,-1.90228454644100,3.39015200689580) q[1];
u3(1.83535963398274,0.521375894739217,2.49978221604605) q[3];
u3(0.425070424671566,-1.49827057484518,2.43703952085539) q[4];
u3(0.897056980249959,-1.74938964232985,0.593841075957959) q[0];
cx q[0],q[4];
u1(1.73217035742061) q[4];
u3(-2.60580019358177,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.247599718067754,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.41080403239667,-2.40576432091201,2.45441342180903) q[4];
u3(1.22687052378174,3.92375449119904,-1.74411109079701) q[0];
u3(2.37740967401234,0.615250174643487,0.671262114123139) q[1];
u3(1.34247499579281,-0.736062781213629,-3.95241995372535) q[2];
cx q[2],q[1];
u1(1.67543214782655) q[1];
u3(-0.0141911242482202,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.79434608067580,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.21985395554582,0.596395897225705,2.76509708191967) q[1];
u3(0.577456076939009,-1.43029568190396,2.64632048669859) q[2];
u3(2.74659011300293,-1.91625761970430,3.95868466710985) q[4];
u3(1.33836696428357,0.974309504175693,0.762435914001990) q[3];
cx q[3],q[4];
u1(0.597674890665125) q[4];
u3(-1.39856974352790,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.04804138439634,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.75402870323376,3.84055335636048,-0.669230644872371) q[4];
u3(0.994545125134237,-5.59006822670715,-0.649874773773325) q[3];
u3(1.62559083368489,-0.475031213827868,1.21226347816853) q[0];
u3(1.43220306273613,-0.720191698496994,-1.12902045969922) q[2];
cx q[2],q[0];
u1(1.90802060931102) q[0];
u3(-1.60878314399476,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.86885862920064,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.18533155956682,-0.295484369770980,0.219747839188304) q[0];
u3(1.46883644103334,1.64460552519062,2.54260420017008) q[2];
u3(1.54983954805858,0.824649680055755,-1.04003140784798) q[4];
u3(0.755368804022661,0.305614880030239,-3.30533696103965) q[3];
cx q[3],q[4];
u1(1.46599177426955) q[4];
u3(-1.01475804556764,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.64606233708717,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.41203125542207,0.404953536351989,-1.10804575553676) q[4];
u3(2.30122800256821,-3.84660042232513,0.0478095184308929) q[3];
u3(0.321585664611919,-1.50723789543308,1.03613502471585) q[3];
u3(0.829414129322464,-2.81095652957917,2.08482478557037) q[2];
cx q[2],q[3];
u1(1.27616688800650) q[3];
u3(-0.841407059638679,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.88840306476245,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.790222247404757,2.75352563911686,-0.482003217981956) q[3];
u3(2.17016730867938,-0.644478789873339,-3.99547975605481) q[2];
u3(2.07985907371485,2.31408983717317,-2.72926086823012) q[4];
u3(1.45158789009712,2.06812252267997,-1.61411988156260) q[0];
cx q[0],q[4];
u1(1.31772421459881) q[4];
u3(-0.364252272486082,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.40102926533099,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.75697771222810,-0.162398965204421,1.00088070251585) q[4];
u3(1.89923591256339,0.0586180314788316,3.08562310487535) q[0];
u3(2.62312715569562,-3.81723271089727,2.45743329634277) q[1];
u3(0.858638944693629,-0.264947124492772,3.15833062617324) q[0];
cx q[0],q[1];
u1(1.06146848164640) q[1];
u3(-0.238923703851718,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.22555391328040,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.650660152377408,-1.61252889224568,1.97281771890305) q[1];
u3(1.07747036773830,-1.03922433649993,-0.235628989108539) q[0];
u3(0.834530969315382,-1.11024332855203,1.55819431573886) q[2];
u3(0.343030648460753,-3.97468967345828,2.27731953397256) q[3];
cx q[3],q[2];
u1(1.70181309928031) q[2];
u3(-0.380272557549477,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.40618392386318,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.20129983160509,-1.19573023473952,1.96534292650675) q[2];
u3(2.67177951483187,-2.81076603883787,-0.632391151097498) q[3];
u3(1.35297528822742,0.00763439463184995,0.344417407156091) q[3];
u3(1.30242143771048,-1.89387050305736,-1.81412866134809) q[1];
cx q[1],q[3];
u1(-0.0320598394976972) q[3];
u3(1.34459215355858,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.35086716001383,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.917223766763089,1.00938426778207,0.427717221565873) q[3];
u3(1.17434007494018,3.75706775896337,-2.51122227338218) q[1];
u3(2.56748347857504,2.58741634362497,-1.03139225109898) q[4];
u3(1.95121866765816,1.90786544723167,-2.44079242083467) q[0];
cx q[0],q[4];
u1(-0.324118424948318) q[4];
u3(1.19240220836218,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.70967030451515,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.865811375225027,-1.19725710139940,1.30887673424991) q[4];
u3(2.14656813154530,4.10683807647506,1.68024497710807) q[0];
u3(2.08794487431910,3.76228537098069,-1.50067153127412) q[1];
u3(1.72746455520635,1.80615441172613,-2.13180346526158) q[4];
cx q[4],q[1];
u1(1.56517011951200) q[1];
u3(-3.03669527301182,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.409663607376963,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.83511127570863,-1.17897937414980,3.99921344764982) q[1];
u3(2.16467317314367,0.816952256866668,1.28109007753880) q[4];
u3(0.682073257715289,0.348539728525273,-2.17323471313678) q[2];
u3(1.15852162354965,1.78566425909054,-4.36449277115357) q[0];
cx q[0],q[2];
u1(2.91888651432437) q[2];
u3(-1.68724792635934,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.524669385218808,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.868945910711650,2.33761385296287,0.976550229497022) q[2];
u3(1.68660587316501,0.420742831973471,3.22143636178216) q[0];
u3(1.85444649002701,0.381537458510054,0.616649269110801) q[4];
u3(1.11212423749869,-1.63740324744540,-1.81656945065109) q[1];
cx q[1],q[4];
u1(0.352247792291535) q[4];
u3(-1.80406071041064,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.21324580690403,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.04185078566295,1.91208575673700,-3.49989258637951) q[4];
u3(2.26428852852427,-1.89052023817343,0.212737621680049) q[1];
u3(1.37390719295741,0.792940678881509,1.95931397519800) q[2];
u3(0.459898424986997,-0.881121746204054,-2.55211494943958) q[3];
cx q[3],q[2];
u1(1.16803333228902) q[2];
u3(-0.275710726812213,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.90191375437270,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.26673266253831,1.12092493692492,1.88157193707761) q[2];
u3(2.94889892709146,-2.83394060106765,1.18243307512054) q[3];
u3(2.07356734334656,0.913212767097329,2.22701230734820) q[3];
u3(1.75204340646224,-2.52459205516963,-2.53916861233759) q[4];
cx q[4],q[3];
u1(-0.830334172331328) q[3];
u3(0.282835686419601,0.0,0.0) q[4];
cx q[3],q[4];
u3(4.11668531095060,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.67666503480850,-3.07192596066781,1.38280959494424) q[3];
u3(0.874061227306398,2.20339776796867,-1.24389283627046) q[4];
u3(0.605851074837271,0.607569717573710,-1.38548740380939) q[0];
u3(0.820628991638835,-0.0758366423739180,-2.01596223480537) q[1];
cx q[1],q[0];
u1(1.46207851775989) q[0];
u3(-3.28408043473139,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.19824430584376,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.794630494014872,-0.542600256943646,-0.127912756045313) q[0];
u3(1.68912983853518,1.32496098971009,-3.66239815457049) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
