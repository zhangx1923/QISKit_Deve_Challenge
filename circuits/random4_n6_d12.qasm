OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.98667348442900,0.558498160028733,-2.86933708710482) q[1];
u3(0.661215669175861,-3.16390855085156,2.61015756230353) q[0];
cx q[0],q[1];
u1(2.98787271050485) q[1];
u3(-1.89081086524981,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.35559617532581,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.75709743185114,-1.39602193229422,2.93585078757566) q[1];
u3(1.58917924114317,-0.635996524940733,-0.281144857603916) q[0];
u3(1.49957696776331,-1.94368248654356,-0.153406348206434) q[3];
u3(0.584818274296636,-3.60234838649245,-0.901825858912975) q[2];
cx q[2],q[3];
u1(1.67521227235588) q[3];
u3(-2.98737851868540,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.537892018289109,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.606631880383162,3.31977631436474,-0.537823589969083) q[3];
u3(1.33634044266248,-0.417963767317812,-1.65914704242739) q[2];
u3(0.928093004614450,0.380001233785322,1.43835564651259) q[4];
u3(1.26263225858193,-0.922200270178353,-2.21237768811676) q[5];
cx q[5],q[4];
u1(1.42744486724113) q[4];
u3(-0.812632196785659,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.99938778325344,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.30666936364052,-1.08667726218267,0.175407774313925) q[4];
u3(2.17770948628538,2.01393366078281,-0.671417439982654) q[5];
u3(2.13854598212581,-1.11005080040051,-0.376497438147023) q[1];
u3(0.864854071713477,-4.30152390916653,-0.390986461284158) q[2];
cx q[2],q[1];
u1(2.48842873370327) q[1];
u3(-1.81044330738531,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.31997478813840,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.85869935674249,-1.63156362171880,2.86769245488318) q[1];
u3(0.669930521737558,-1.04841808015751,4.41480057443256) q[2];
u3(1.68123023374700,-0.174540466111029,1.56005714738778) q[0];
u3(2.04533717621394,-1.50597153732947,-0.213953483351452) q[4];
cx q[4],q[0];
u1(3.20015253880059) q[0];
u3(-1.79602637824626,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.952043637795506,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.03050126395438,-0.717876515001910,2.25011729999963) q[0];
u3(2.40389526556410,1.20880693057249,-4.45359931326784) q[4];
u3(1.23883223616221,-0.838484859588885,-1.87825930724978) q[3];
u3(0.856294973706950,-4.93128765926732,0.944373832882795) q[5];
cx q[5],q[3];
u1(1.84459439306599) q[3];
u3(-3.15306167507368,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.55429567973093,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.41116056991091,-3.06055636090528,2.57039647818799) q[3];
u3(1.41163284604806,-3.88969066244627,-1.52303419408114) q[5];
u3(1.27063871424522,1.13288396611506,-2.80734488092255) q[3];
u3(1.51883943775964,-2.82678150051961,2.90847749666472) q[4];
cx q[4],q[3];
u1(1.54089649085505) q[3];
u3(-3.23318758220770,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.88401849586962,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.80542385758926,-1.27172151787036,-1.91729146139603) q[3];
u3(0.882985270946342,2.30230458728287,0.979216988159652) q[4];
u3(2.78237619809674,-1.79752285164374,1.37663495212810) q[5];
u3(2.66540739950813,1.67755341137111,3.75299458632169) q[2];
cx q[2],q[5];
u1(2.13501578243681) q[5];
u3(0.0149918043881045,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.30326128921270,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.71760747903059,3.51971104810848,0.121474170377454) q[5];
u3(1.28306913064857,1.03186745335516,4.92264772617534) q[2];
u3(2.31568513273900,-0.601868137911826,2.40670897027323) q[0];
u3(2.60823003944749,-2.96259282693302,-2.10790694872950) q[1];
cx q[1],q[0];
u1(1.17290792906366) q[0];
u3(-3.50014000021758,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.62035046682043,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.80041150648219,-2.49477953090611,-0.307751469659643) q[0];
u3(2.97442153857793,-0.434293982211142,-3.15430025941873) q[1];
u3(0.986454990682616,2.92462874636410,-1.92803968307548) q[4];
u3(0.638205660544429,1.70224183128752,-3.38618667577696) q[1];
cx q[1],q[4];
u1(0.285292667463923) q[4];
u3(-0.606808851192530,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.66035201228058,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.15865316763478,-0.0942585602963418,0.279192420605493) q[4];
u3(3.00893725712486,-1.78757334152387,-0.188410855042393) q[1];
u3(2.14521160787430,-0.721854789350366,2.36864134857116) q[3];
u3(2.76029827957141,-1.06127577005679,1.17385745694223) q[2];
cx q[2],q[3];
u1(1.63555793203115) q[3];
u3(-2.74864983021692,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.944418446542684,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.45898305376081,-1.30566374575961,0.735789637591018) q[3];
u3(0.429577553681186,-4.15911539990235,0.763037422988702) q[2];
u3(2.54362979967311,1.27546250712672,-1.23121103707702) q[5];
u3(1.94688354096252,0.807358770462207,-3.00450422078323) q[0];
cx q[0],q[5];
u1(3.04724847813130) q[5];
u3(-2.23581774230128,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.08121945333239,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.90158811875090,-1.31627464502830,-1.13901965015792) q[5];
u3(1.33847558341907,1.09531379664024,-3.80949489352817) q[0];
u3(1.92218922912120,0.319809255576742,2.01604815356712) q[3];
u3(1.29921564326135,2.94539759241042,2.92381535197429) q[0];
cx q[0],q[3];
u1(0.669320607161161) q[3];
u3(-1.28571021739147,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.85313822355765,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.83855906468877,-0.799232819438644,0.454625138861792) q[3];
u3(1.10614639341248,-2.96519454624035,-2.44573755473110) q[0];
u3(1.78406068310569,0.343922299069154,2.65676904648179) q[2];
u3(1.78993963708158,-2.82000323971953,-1.68141708922574) q[1];
cx q[1],q[2];
u1(0.281456492928772) q[2];
u3(-1.37209208919030,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.04190426528161,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.60638066245330,1.35749048345429,-4.27823291282614) q[2];
u3(1.41487949158174,0.288556661518344,-1.38764357316384) q[1];
u3(1.77054034435050,2.46990756406103,-0.469758569185558) q[4];
u3(2.66669538741750,0.786590400686193,-2.31515745144081) q[5];
cx q[5],q[4];
u1(2.86808304953644) q[4];
u3(-1.94314755946350,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.43991139276593,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.742048160805452,0.209672890730818,0.736763359277484) q[4];
u3(1.77492156808624,-4.92839407164907,0.707863712850670) q[5];
u3(1.81518112098397,2.68171094555869,-3.38261613525718) q[0];
u3(0.825613754402018,3.40710474092526,-1.89419693859100) q[3];
cx q[3],q[0];
u1(1.13439455846646) q[0];
u3(-0.630278216245581,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.11076944178020,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.05210516464382,2.20870185613144,0.788080449955217) q[0];
u3(1.45392969409939,1.26238000594497,-2.04720877407898) q[3];
u3(1.79505690610208,3.63502962906713,-1.16700589305978) q[2];
u3(1.67467432533223,1.80935907374272,-0.699227649336256) q[4];
cx q[4],q[2];
u1(2.12173167170936) q[2];
u3(-1.65313445922599,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.294073141082243,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.06239346692279,-1.12950652573602,-1.01838012882957) q[2];
u3(2.18958545637654,0.792398560771933,3.50765672950386) q[4];
u3(2.33594847208816,2.76477013696230,0.285126147022078) q[1];
u3(1.11858419573215,-0.261314863403523,-2.32253969619689) q[5];
cx q[5],q[1];
u1(0.198329605727959) q[1];
u3(-1.21860344872552,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.30112748673798,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.23864460847695,-3.59405151394159,1.90200877952381) q[1];
u3(1.97849489663738,1.06260865726328,2.13062157824011) q[5];
u3(1.77331051581648,0.679429676591144,2.42554976936411) q[3];
u3(2.31466255968424,-2.46099192793046,-2.44788542445089) q[4];
cx q[4],q[3];
u1(-0.383286120762504) q[3];
u3(-1.06166190597249,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.60933107780134,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.73311111285777,1.96936966986130,-3.18189516637492) q[3];
u3(1.48658339009308,-1.09325005130754,4.22040146822865) q[4];
u3(0.103865685724235,1.64045031940204,-0.593809047362114) q[2];
u3(0.443482629555309,-2.28047811588745,0.793343659060499) q[5];
cx q[5],q[2];
u1(0.536097335743802) q[2];
u3(-0.917883667594334,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.52334534174215,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.77062634378088,-4.13106432312929,0.155478811631634) q[2];
u3(0.647708626534881,-0.860365700684140,-0.822287844505998) q[5];
u3(1.72355335152926,-1.71114697923109,0.781760153558235) q[0];
u3(2.00653386674110,-4.07981240566680,-0.330398165275523) q[1];
cx q[1],q[0];
u1(0.613957088356249) q[0];
u3(-0.900833472587904,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.46401550931696,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.41246806226867,0.627818859759810,1.71954184644660) q[0];
u3(0.800654140221138,0.734749583997712,3.56180282851923) q[1];
u3(1.96055460683832,0.862659532535469,2.09573168295832) q[0];
u3(1.79885639751028,-0.706147093245452,-1.19693378045569) q[2];
cx q[2],q[0];
u1(2.25775664792064) q[0];
u3(-2.47410831961363,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.28138873568884,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.11606780106276,-3.73538068214270,1.29614438242562) q[0];
u3(2.12329404559430,-0.480383444853613,-4.22724959749289) q[2];
u3(1.17137178469680,1.52067453645909,0.431754525569446) q[1];
u3(2.37475131235251,-0.114522187968148,-4.07781930224827) q[5];
cx q[5],q[1];
u1(-0.563778480998584) q[1];
u3(-1.96278086120438,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.54533337566628,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.69290446061946,3.26181891894656,-2.31052640161426) q[1];
u3(1.96733964229914,-0.602043575349058,-1.37692173864485) q[5];
u3(1.87124124898607,-0.333960473848505,1.57411473978347) q[4];
u3(1.37404225734430,-2.40241543427691,-2.36550287069335) q[3];
cx q[3],q[4];
u1(3.33399327865618) q[4];
u3(-0.489788335749547,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.38410384684802,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.14566134201955,-3.10402913244926,2.09770221241262) q[4];
u3(1.19153331956931,0.624405009352512,-4.43525836159584) q[3];
u3(1.14330402568606,2.78923073601757,-1.77580976081764) q[3];
u3(0.413760913427346,-2.12249177272949,0.434039428763190) q[0];
cx q[0],q[3];
u1(2.68930745490141) q[3];
u3(-1.50842918399491,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.38725973362959,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.943223226354310,-2.02142930790287,2.14168487495512) q[3];
u3(0.911715224310936,-3.13292084584173,1.57535773984455) q[0];
u3(1.68799162834127,-2.14994684483446,-0.199265614026020) q[4];
u3(1.79295825636929,-4.13069752895609,-1.03601380500741) q[1];
cx q[1],q[4];
u1(1.40981088862520) q[4];
u3(-0.545623913177741,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.91257800813557,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.46830144265746,-2.33747301847034,-1.48907671935883) q[4];
u3(0.236425472676046,1.85213995857649,1.82184686536646) q[1];
u3(2.62765081178115,-1.29904174159258,2.23664533989045) q[2];
u3(2.18690721400514,1.84957576944259,3.64599604656220) q[5];
cx q[5],q[2];
u1(1.85213362954553) q[2];
u3(0.0164832312664684,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.20029944547530,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.64869930125556,2.11759538281050,1.14141204196679) q[2];
u3(2.51070889371913,0.893551987109906,-1.21931107301068) q[5];
u3(1.43381704599919,-0.200722858473570,1.37395599570202) q[0];
u3(1.05963917495155,-2.79615381715117,-2.02074115195518) q[1];
cx q[1],q[0];
u1(0.528362621027389) q[0];
u3(-0.783479828801241,0.0,0.0) q[1];
cx q[0],q[1];
u3(4.28862273562983,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.68858306659635,0.322368841673019,-2.42259760472805) q[0];
u3(1.00417904733287,2.60822697308049,-3.04361764098370) q[1];
u3(1.74251775690714,-0.784405706504798,1.26462778767753) q[5];
u3(2.57226860023910,-1.02193578221139,-2.01855918234720) q[3];
cx q[3],q[5];
u1(1.55331219876518) q[5];
u3(-0.683068762545757,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.11465321825356,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.62513213097665,-1.98968196544525,1.01270901832440) q[5];
u3(2.36612172297656,3.93136294109072,-0.179786278170259) q[3];
u3(1.37574655728187,0.743893425147937,-3.65616567899124) q[2];
u3(0.984255400268601,5.00596313467342,-1.17484461902079) q[4];
cx q[4],q[2];
u1(2.68157338565422) q[2];
u3(-1.78886987005780,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.21044388764730,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.68056196100536,-0.0802142185253648,-1.33525054079326) q[2];
u3(2.59840478431715,1.66912482574295,-1.12589449159814) q[4];
u3(2.27028462830128,0.147731357381733,1.17798207548066) q[5];
u3(0.775060084498907,-4.47206310239560,-0.303680929263350) q[4];
cx q[4],q[5];
u1(0.163826778543686) q[5];
u3(-1.23236363408016,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.32584998079741,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.638713186980355,-1.51072923986665,3.34074460922510) q[5];
u3(0.522438130006668,1.80928636511633,2.88744824757216) q[4];
u3(1.80156138767236,-0.788099446425447,-0.803447545077577) q[0];
u3(0.722733400556511,-5.29384859281914,0.783899484249714) q[1];
cx q[1],q[0];
u1(2.90929438960257) q[0];
u3(-2.85108171640277,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.06740696257295,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.64179963140565,1.86157109812208,-1.60924251664394) q[0];
u3(2.28262128609064,-0.558886478135555,5.20943235523649) q[1];
u3(1.83047941379129,0.932364435955006,0.830997927845969) q[2];
u3(0.288397193127221,-4.56579804591771,-0.745209546805349) q[3];
cx q[3],q[2];
u1(2.13144492604271) q[2];
u3(-2.80251978334437,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.34622430236425,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.43179789886917,2.38971406031980,0.124143462391291) q[2];
u3(1.62521155535426,0.586794509443082,4.94095715758472) q[3];
u3(0.875193992674954,-0.930915545488513,2.01850429576955) q[4];
u3(0.512021121248045,-1.71813641895542,-0.287654880796886) q[0];
cx q[0],q[4];
u1(3.07272307603582) q[4];
u3(-2.68046717876197,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.754697562360705,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.31273582958184,2.62852516161538,-1.00378335321790) q[4];
u3(0.850969548513002,2.92993124716882,-2.61457483317392) q[0];
u3(0.835276318862238,0.171522661922470,-0.532768086833264) q[5];
u3(1.13019430396970,-3.91443469363291,1.11456915271894) q[3];
cx q[3],q[5];
u1(3.14100615883306) q[5];
u3(-1.67091101434953,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.715153662496449,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.37398488233931,-0.0751100856173518,-1.23493270795516) q[5];
u3(1.52669322431690,0.376954665064380,-1.85004542824691) q[3];
u3(1.70414129888549,2.81252503975558,-1.47724517530830) q[1];
u3(1.84544327902255,0.669252311951275,-2.93029220019187) q[2];
cx q[2],q[1];
u1(1.06358383201502) q[1];
u3(-1.35982552805771,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.223461378624833,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.361038852218177,0.425594460010759,1.27153410048812) q[1];
u3(2.71787231369341,1.60466198739903,2.31042295027676) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];