OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(0.362235102394486,-0.727892232719277,0.437844204132748) q[0];
u3(0.622636106951338,-0.571094457773547,-1.90372228775227) q[7];
cx q[7],q[0];
u1(0.524660411460160) q[0];
u3(-1.60319058415006,0.0,0.0) q[7];
cx q[0],q[7];
u3(-0.155904478429141,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.13790418081902,0.891411452001624,3.31207616941233) q[0];
u3(1.96128732105566,2.62312861902505,-3.12911088640268) q[7];
u3(1.87187810284829,1.26834361082571,-3.21211768179215) q[1];
u3(3.03046210971676,3.72741387039885,-0.899969725104091) q[3];
cx q[3],q[1];
u1(0.894123275454403) q[1];
u3(-0.177364572528668,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.69540427386205,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.30090013423047,-0.453730137423709,1.73962581967968) q[1];
u3(2.95793277619174,4.73506050426357,0.485785403882230) q[3];
u3(1.68732489259493,2.07420576381010,0.0376395167332035) q[6];
u3(2.85521692788326,0.759285486316156,-3.96625298856554) q[11];
cx q[11],q[6];
u1(3.50555216175875) q[6];
u3(-0.896766727674819,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.60641109413537,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.31788180187956,-0.238418309894324,2.15807383721992) q[6];
u3(0.986558124222800,0.924845659893299,-1.77193855017566) q[11];
u3(0.591932822839610,-0.462185586847639,-0.574026834498068) q[4];
u3(0.602700556190742,-0.621621164549753,-0.872597858437236) q[12];
cx q[12],q[4];
u1(1.89638191028897) q[4];
u3(-2.67467952838683,0.0,0.0) q[12];
cx q[4],q[12];
u3(-0.0338800072902441,0.0,0.0) q[12];
cx q[12],q[4];
u3(2.22150489669231,1.63980707598314,-2.19934475460762) q[4];
u3(0.453247139468112,3.25198531430451,0.927376013616680) q[12];
u3(2.65075820729741,0.413117328275545,-2.40284111954865) q[10];
u3(2.11980928243034,5.74055704927713,0.276561856885859) q[8];
cx q[8],q[10];
u1(0.567573522538166) q[10];
u3(-1.05995601618933,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.75411826675993,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.32812546597696,0.657647160931091,-1.21690344908655) q[10];
u3(1.93536933638251,1.07619700554396,1.27288321651737) q[8];
u3(0.821323247452032,-0.607185012518161,1.39247840823151) q[9];
u3(0.775245169529578,-1.32984192958000,-1.30816745082525) q[13];
cx q[13],q[9];
u1(0.713173209795297) q[9];
u3(-3.26287697278689,0.0,0.0) q[13];
cx q[9],q[13];
u3(1.65075871052470,0.0,0.0) q[13];
cx q[13],q[9];
u3(1.66430146679369,-1.04640469991998,-0.334351602752654) q[9];
u3(2.20418866200167,1.11596673982250,1.81771522814741) q[13];
u3(1.84922636777034,-2.19346486439018,1.08702825622155) q[5];
u3(1.62649789976806,-3.46402663928055,0.235558693272842) q[2];
cx q[2],q[5];
u1(2.11377419901134) q[5];
u3(-2.51878057430789,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.02293960362446,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.21473488436160,1.24260620559115,1.04008598450750) q[5];
u3(0.929187776715455,0.411780323077879,-5.15042471635733) q[2];
u3(2.42393240856941,-1.08712214815998,0.509381015071950) q[7];
u3(1.43532114488152,-4.28339209370542,-0.582601477177125) q[10];
cx q[10],q[7];
u1(1.17246784404532) q[7];
u3(0.0357449785270727,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.58992929676527,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.842317957237689,1.51673374876769,1.55290487507976) q[7];
u3(1.83112662601951,-3.54091720440243,1.72275416394445) q[10];
u3(0.516142908116782,-1.76262495032302,1.26108896166031) q[3];
u3(1.04070357659926,-0.786181504774366,-0.825826312529686) q[1];
cx q[1],q[3];
u1(1.71817835989619) q[3];
u3(-2.38618457759620,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.520495431343514,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.36438322080907,1.86091844501475,0.592802983580509) q[3];
u3(2.18500718982086,-0.573456538523438,0.947171339157404) q[1];
u3(0.396307680925216,0.113960727244583,-1.41813484507176) q[13];
u3(1.50710372636885,-3.71461908393965,2.13006400870442) q[6];
cx q[6],q[13];
u1(1.81711181705686) q[13];
u3(-2.56986593102738,0.0,0.0) q[6];
cx q[13],q[6];
u3(0.870084661973668,0.0,0.0) q[6];
cx q[6],q[13];
u3(0.702791904802999,2.45400530713641,-1.74528437938145) q[13];
u3(2.16011355572257,-1.35438712857764,-2.07140800572975) q[6];
u3(2.20937160385558,-0.0881537804360231,-2.10457615342092) q[9];
u3(1.47713812190636,1.07806549654232,-4.55420537422695) q[8];
cx q[8],q[9];
u1(2.50085140666103) q[9];
u3(-1.74997590274950,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.334581919384963,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.20243350017313,0.993857331024457,-3.94023680894047) q[9];
u3(1.69222405266807,-0.0620869029994333,-5.79331263743969) q[8];
u3(1.15640402845069,3.55680185727800,-0.901816020757610) q[12];
u3(1.16228760414194,1.48969547006851,-1.30612809066804) q[2];
cx q[2],q[12];
u1(0.495135240128241) q[12];
u3(-1.43059884729701,0.0,0.0) q[2];
cx q[12],q[2];
u3(-0.216138927338456,0.0,0.0) q[2];
cx q[2],q[12];
u3(0.965438359292337,0.184891104584098,-2.78123992457025) q[12];
u3(1.26897695233142,-0.939111515824693,2.96212019578499) q[2];
u3(1.21255940614263,3.10257128640351,-1.88782813936795) q[5];
u3(1.50035305474500,0.657675266070991,-2.82485033583410) q[0];
cx q[0],q[5];
u1(1.55793723576456) q[5];
u3(-2.24390301756227,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.60440485423643,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.10399577306025,-2.55179264961957,2.87000654928185) q[5];
u3(1.04934367280034,-0.813313020193318,-0.944874384495248) q[0];
u3(1.90840366736880,-2.00793247798322,-0.911265218573982) q[11];
u3(2.31554431501078,-3.27931035007526,-0.363148027130897) q[4];
cx q[4],q[11];
u1(2.36742316633302) q[11];
u3(-2.05139200714557,0.0,0.0) q[4];
cx q[11],q[4];
u3(3.34048203355100,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.12599181369726,1.09029989995718,-3.25599310221084) q[11];
u3(1.88702270814668,-0.890700989378640,4.96018874079759) q[4];
u3(1.66798870188815,2.00899601602240,-0.386077537009523) q[7];
u3(1.32801243058180,0.799846743512288,-3.35233507953838) q[10];
cx q[10],q[7];
u1(1.75939969612732) q[7];
u3(-2.37864849393380,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.268786038731289,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.92633664163478,-3.16531987662245,1.82419537742196) q[7];
u3(2.01319505077332,3.55390347133784,-2.71030682488519) q[10];
u3(1.68403739654724,1.60615702170655,-1.16643549434117) q[8];
u3(1.50239583392029,-1.15686662357131,-3.26879659820631) q[2];
cx q[2],q[8];
u1(1.74310293844116) q[8];
u3(0.0543931419516090,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.621413512624507,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.51905679185081,-1.09590802317683,-2.34461068143806) q[8];
u3(1.89401882851933,1.99595568606302,-1.26854176127612) q[2];
u3(1.99568374177811,-0.369913330831145,-0.448346684456248) q[9];
u3(1.57173135832357,-3.02651608303010,-0.141563746408434) q[6];
cx q[6],q[9];
u1(3.94976632553612) q[9];
u3(-3.73775626749045,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.215594635528777,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.01529613788765,-3.27521212074202,2.92919394632835) q[9];
u3(1.93789870905535,-1.31110300411537,4.45821184396200) q[6];
u3(1.51108365727567,1.77859300867148,-3.09613161085862) q[13];
u3(2.32061798449934,-2.58772891576448,3.04347209699580) q[3];
cx q[3],q[13];
u1(2.79298677961965) q[13];
u3(-2.25166868232100,0.0,0.0) q[3];
cx q[13],q[3];
u3(1.39420099151661,0.0,0.0) q[3];
cx q[3],q[13];
u3(1.12863604279208,-2.76977926307672,1.17986051505122) q[13];
u3(1.02804530808469,4.85820586770775,-0.489169777455356) q[3];
u3(1.48888968374108,-0.452079480440776,1.73090506586120) q[11];
u3(1.80339256814198,-1.34154243867059,-2.34814824650303) q[12];
cx q[12],q[11];
u1(3.42133037742518) q[11];
u3(-2.12964563790175,0.0,0.0) q[12];
cx q[11],q[12];
u3(1.49658687109216,0.0,0.0) q[12];
cx q[12],q[11];
u3(0.0864309492203234,-2.91630084253810,1.67233816322727) q[11];
u3(1.07422338430065,-1.08861917348977,2.02571176099567) q[12];
u3(2.17560772277898,0.858051482294211,1.11981202316906) q[1];
u3(0.602312189523793,-4.72859571712512,-0.401520791822761) q[4];
cx q[4],q[1];
u1(2.15993532747766) q[1];
u3(-2.73894905265006,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.17788772763284,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.64278413856550,-3.67829037645858,2.19481937791371) q[1];
u3(2.80656043976951,2.26248888817376,0.0952466173440472) q[4];
u3(0.959488107568707,1.34723253688130,-2.20667167782126) q[5];
u3(2.23219296794526,-4.08643522296668,2.15499530246257) q[0];
cx q[0],q[5];
u1(0.591091904308388) q[5];
u3(-1.27284993547552,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.10283023604103,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.81868822204663,-3.06362385911885,2.83250000534333) q[5];
u3(2.12303879000881,-0.416801148866850,1.74859354355355) q[0];
u3(0.251433606618930,2.69095712029450,-3.24811596441683) q[10];
u3(1.14810323944632,-2.70166046334296,2.06570982432125) q[8];
cx q[8],q[10];
u1(-0.0783235566381086) q[10];
u3(-1.68088323205181,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.857722790918431,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.85948542138345,-0.809581687228784,2.37079454287037) q[10];
u3(0.686324828330826,-1.27809057373158,1.08483100732302) q[8];
u3(1.43434383487899,0.533782154377932,-2.61519126481543) q[12];
u3(1.80594053047690,-4.15696978588540,1.73876440455181) q[9];
cx q[9],q[12];
u1(1.02810785052211) q[12];
u3(-0.0399226311715994,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.89193883325566,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.03848142540242,-1.75440998327491,2.08793990346378) q[12];
u3(2.60763940326982,3.25984387366480,1.86486960194357) q[9];
u3(0.824128545694576,2.58134740291720,-0.0948444397781281) q[7];
u3(1.26873074351437,1.26282471661179,-1.26420944954491) q[6];
cx q[6],q[7];
u1(-0.682321253882580) q[7];
u3(1.39416246643441,0.0,0.0) q[6];
cx q[7],q[6];
u3(3.75963279401495,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.53716698318551,-2.72804064418111,2.37553180263798) q[7];
u3(0.906761374387449,4.62413925973326,0.241506074781031) q[6];
u3(0.689256051482983,-0.747545669341009,-0.538079743488238) q[1];
u3(1.63630487396284,-2.92817196803290,1.14546656065155) q[13];
cx q[13],q[1];
u1(2.26479842763967) q[1];
u3(-3.04429380508947,0.0,0.0) q[13];
cx q[1],q[13];
u3(1.49422226179048,0.0,0.0) q[13];
cx q[13],q[1];
u3(0.698980233111809,-0.398620861387648,0.151743593435067) q[1];
u3(2.31791431193939,0.306308043686463,-0.743008021069527) q[13];
u3(0.985995037330092,-1.90771049238109,0.351280194684781) q[5];
u3(1.28580924232675,-1.75041843969755,0.841357118930385) q[4];
cx q[4],q[5];
u1(0.931109772485291) q[5];
u3(-0.236759915851401,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.95287660001152,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.08972417984016,3.07809446329642,-1.76398740831741) q[5];
u3(0.699136476000280,2.58489231667410,3.12915147164730) q[4];
u3(2.47743644186713,-1.02480283096370,-0.651896516705059) q[0];
u3(0.686126938288737,-1.12160915208100,-4.59986867378272) q[3];
cx q[3],q[0];
u1(3.03930010903459) q[0];
u3(-1.95858186899020,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.974987339507918,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.824615641093163,0.755962687260453,-2.84933788713730) q[0];
u3(2.97624339192919,-2.85607596684660,-1.87358507607672) q[3];
u3(1.55141104488358,2.84802102590439,-2.45964113886765) q[2];
u3(1.39726648823579,2.92988604066700,-2.69192022363676) q[11];
cx q[11],q[2];
u1(1.61954090468226) q[2];
u3(-2.96779556744900,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.925832309583176,0.0,0.0) q[11];
cx q[11],q[2];
u3(2.16276308613176,-1.36788685440778,1.68163737270783) q[2];
u3(0.914503057846917,2.03145738464574,2.19672222644580) q[11];
u3(2.78252166034250,-0.420271634341645,1.81310497369356) q[8];
u3(2.76954352389032,-0.864621934486104,0.798070295411409) q[1];
cx q[1],q[8];
u1(1.39126542559708) q[8];
u3(-0.838757714070206,0.0,0.0) q[1];
cx q[8],q[1];
u3(-0.495149137434608,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.00866142366385,-0.470378817571573,4.55844309341140) q[8];
u3(1.82242168407847,-1.57846311548915,1.89509016469310) q[1];
u3(1.54255077881007,0.721194749242579,1.51075217598954) q[6];
u3(1.92252322203303,-1.35657671171129,-1.01377647591422) q[13];
cx q[13],q[6];
u1(2.39245079828689) q[6];
u3(-1.44092515825638,0.0,0.0) q[13];
cx q[6],q[13];
u3(0.0870708699488365,0.0,0.0) q[13];
cx q[13],q[6];
u3(2.05702313629550,-4.78576876082752,0.902927833439418) q[6];
u3(0.489794778945523,-0.0934870779297818,4.41191523520154) q[13];
u3(0.348091565252664,-1.56856249700831,0.832489365268433) q[9];
u3(1.06315618210723,-0.242307716185224,-0.937660111089170) q[11];
cx q[11],q[9];
u1(-0.0373635797860110) q[9];
u3(-1.70050699512432,0.0,0.0) q[11];
cx q[9],q[11];
u3(0.530356912024628,0.0,0.0) q[11];
cx q[11],q[9];
u3(0.593402972820856,-1.28665364184230,0.336015601419218) q[9];
u3(1.14443377306479,-1.17825863924304,-0.748818404696342) q[11];
u3(1.54828396137159,0.746704107837531,-0.685338243389305) q[4];
u3(0.971651274609086,-4.69332158339188,1.40136645198327) q[10];
cx q[10],q[4];
u1(3.61281487905893) q[4];
u3(-4.44631273721528,0.0,0.0) q[10];
cx q[4],q[10];
u3(-0.699785271987589,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.51856833794556,0.598724883068843,2.66000679478594) q[4];
u3(2.60113622486032,0.326425854122189,2.04857559760175) q[10];
u3(1.49108125191670,1.40195825577281,-1.44336740177839) q[7];
u3(0.269536429690040,1.49826093248533,-4.50528218915780) q[3];
cx q[3],q[7];
u1(2.32077294867552) q[7];
u3(-1.65312424005382,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.388363561680745,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.32234824805403,0.477442273872791,-0.486317631920145) q[7];
u3(1.07787843683749,-2.98903130307778,-1.35898993987034) q[3];
u3(1.43656456624237,1.41206337129465,-2.58978889764626) q[2];
u3(1.25665276638412,-2.17610990418772,2.24845474956854) q[12];
cx q[12],q[2];
u1(-0.170966879125304) q[2];
u3(-2.31031869892911,0.0,0.0) q[12];
cx q[2],q[12];
u3(1.40526897158416,0.0,0.0) q[12];
cx q[12],q[2];
u3(0.821861009757983,-2.25873720368811,0.931073991650256) q[2];
u3(1.62750542415115,0.836887334790309,0.241552900893785) q[12];
u3(1.23949924464346,2.09352100347366,-0.619324106890433) q[5];
u3(1.40076158797421,0.514066201430960,-2.62366219918664) q[0];
cx q[0],q[5];
u1(0.889967385528138) q[5];
u3(-3.68614043586168,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.80871647392247,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.71253524859779,1.54721927504503,-0.862671397951050) q[5];
u3(0.568928058304992,1.19462471564902,-3.57213258021516) q[0];
u3(1.60747086421068,-0.261397920154582,1.76652989204142) q[6];
u3(1.12240553770859,-1.05780622964353,-0.747743538646833) q[1];
cx q[1],q[6];
u1(0.323621255243787) q[6];
u3(-1.73460292954692,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.10715463898386,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.740966292062598,0.907365711162786,-3.91234077288092) q[6];
u3(2.35973185031232,0.422060808479533,-1.00080178899879) q[1];
u3(0.666152599070632,-0.413157253740059,0.338562316229555) q[12];
u3(1.03026001216172,-3.14125900677503,1.39236902896015) q[0];
cx q[0],q[12];
u1(1.47781864035508) q[12];
u3(-0.729766343012134,0.0,0.0) q[0];
cx q[12],q[0];
u3(3.06156783722572,0.0,0.0) q[0];
cx q[0],q[12];
u3(2.18692310133608,-0.638985736748350,1.94369888551926) q[12];
u3(1.96589826106175,-1.67963412535014,4.50528974740808) q[0];
u3(2.63693381688123,2.64465670411691,-2.92235197021818) q[4];
u3(0.588259657943047,-0.378792333698743,2.17244914423540) q[5];
cx q[5],q[4];
u1(2.06004335277257) q[4];
u3(-2.87214422116536,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.03367338275223,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.51522549782493,0.633456030153310,1.50468731049102) q[4];
u3(1.91232042167458,-5.64704192926214,0.538461939212560) q[5];
u3(1.95528255435704,1.58089633020700,1.01261542213659) q[9];
u3(1.97599450206819,-0.254234014396089,-3.56163288724108) q[2];
cx q[2],q[9];
u1(1.72297848725455) q[9];
u3(-2.78559210118242,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.01734973817645,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.14275079266045,-0.262687520312935,1.65636247795930) q[9];
u3(1.99907924455628,-1.74084669555442,3.40760246232319) q[2];
u3(0.726885006467087,-0.975482349331651,0.901202506371287) q[7];
u3(0.225873257656067,-0.869471281661985,0.614627391636861) q[8];
cx q[8],q[7];
u1(0.246975230574000) q[7];
u3(-1.45432596112254,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.22219997285725,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.96284034197925,1.18370106916072,-2.07452360712319) q[7];
u3(1.42450742169631,2.68628976377696,3.28256381825378) q[8];
u3(1.67183310084052,3.34748628495721,-0.398326711452807) q[13];
u3(1.38717539581822,1.72384446379613,-1.28017364294874) q[10];
cx q[10],q[13];
u1(0.118414242989653) q[13];
u3(-1.35106021698126,0.0,0.0) q[10];
cx q[13],q[10];
u3(2.32552705345590,0.0,0.0) q[10];
cx q[10],q[13];
u3(1.56833810288131,0.974822837550143,-3.57495636111592) q[13];
u3(2.09390581285838,-1.50143037761169,-1.38445500312031) q[10];
u3(2.13997143705854,-3.97849085772217,2.17175632306765) q[11];
u3(0.411768715962982,-1.28720019216353,3.82838326000473) q[3];
cx q[3],q[11];
u1(3.25129872277852) q[11];
u3(-4.27134245909267,0.0,0.0) q[3];
cx q[11],q[3];
u3(-0.502367215359306,0.0,0.0) q[3];
cx q[3],q[11];
u3(0.789227887364967,-0.567518862393408,-1.30778203388241) q[11];
u3(2.31145708829703,3.50800151974550,-0.467222455739386) q[3];
u3(2.41297345727376,-2.44785496459004,0.617396978389833) q[5];
u3(2.33731837458077,-2.12192450946489,1.17306375886831) q[8];
cx q[8],q[5];
u1(0.215730785360009) q[5];
u3(-1.21528821709956,0.0,0.0) q[8];
cx q[5],q[8];
u3(2.48306031615532,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.07860041166904,-0.991842606284634,2.78436778290182) q[5];
u3(1.53036931791928,3.45768564017770,-0.559278891930933) q[8];
u3(0.309676938141947,3.42624680177176,-2.42943226838793) q[4];
u3(0.762777111540875,-3.04123441703277,0.972365031773332) q[3];
cx q[3],q[4];
u1(2.57991306919058) q[4];
u3(-1.48953030348606,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.94788127969418,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.29221460133872,-0.111151992633225,4.23330506114639) q[4];
u3(0.852796758103564,3.44976155602817,2.35643873219170) q[3];
u3(1.93342656937248,0.646420848873653,2.06959995920518) q[2];
u3(1.10891535334389,3.13653049723014,3.04572424419885) q[13];
cx q[13],q[2];
u1(3.23804334678853) q[2];
u3(-1.15766653588758,0.0,0.0) q[13];
cx q[2],q[13];
u3(2.57851517100026,0.0,0.0) q[13];
cx q[13],q[2];
u3(2.23350471056070,-0.274939923191450,-2.77181795061782) q[2];
u3(0.808154150886191,0.219329800077924,-1.61290058503392) q[13];
u3(1.15339116916577,-0.728828546084969,2.86576836947025) q[6];
u3(1.73577358924257,-2.09807541156792,-1.59250167538141) q[11];
cx q[11],q[6];
u1(1.52353506441062) q[6];
u3(-0.637324401499030,0.0,0.0) q[11];
cx q[6],q[11];
u3(-0.0953453674875335,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.50482113564883,3.56397085709217,-1.63876560831120) q[6];
u3(0.516850230528343,0.149664214572151,4.10111098258315) q[11];
u3(2.28762446392585,-0.363731072892842,-0.454818441942119) q[0];
u3(0.773573661992376,-1.04064416135242,-3.80698421084330) q[1];
cx q[1],q[0];
u1(2.00829025350733) q[0];
u3(0.399311816791007,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.10094039110237,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.76077364309440,-0.297781950135788,1.36277268841873) q[0];
u3(0.726308768073180,1.21764002538320,-2.15755402606628) q[1];
u3(2.12905406890082,-0.141392263976179,0.252895850555024) q[10];
u3(2.36048382612733,-1.53521746690146,-1.44604639852505) q[9];
cx q[9],q[10];
u1(0.449402772042934) q[10];
u3(-0.974609563520319,0.0,0.0) q[9];
cx q[10],q[9];
u3(1.39390600989968,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.49321104510411,0.534568317535800,-0.609976987739393) q[10];
u3(0.342509503148962,1.13750999030619,-4.75057108953565) q[9];
u3(1.91970455439432,-0.820271538191631,0.0269559885171432) q[12];
u3(1.54785694138661,-2.56762107275097,0.880686708367162) q[7];
cx q[7],q[12];
u1(3.31276274913898) q[12];
u3(-1.22613911306672,0.0,0.0) q[7];
cx q[12],q[7];
u3(2.42322594443491,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.604888584743693,1.09747291089140,-1.44195720599609) q[12];
u3(2.53883466211872,3.21488712220765,2.85758974595658) q[7];
u3(2.08914197593574,1.16534650905139,-3.39755043281504) q[1];
u3(2.21123850007875,3.51083897704772,-2.67670197592410) q[5];
cx q[5],q[1];
u1(2.35900865599397) q[1];
u3(-1.51385015420753,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.37414322236739,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.99411755910280,1.93673643137038,-1.32904629385673) q[1];
u3(1.37980699304062,-0.699662920003069,3.59668780187720) q[5];
u3(2.25829922466365,0.725723801706236,-3.14102313494098) q[4];
u3(2.77150301881607,2.96010703605668,-2.50584615501827) q[3];
cx q[3],q[4];
u1(-0.179511508066098) q[4];
u3(-2.65726839453742,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.55934447384846,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.79752150086353,-1.81957029105819,-1.34060352307240) q[4];
u3(2.00761037992587,0.0476652999847689,-6.03612082461969) q[3];
u3(0.970929328534777,-2.21271673338682,-0.652347042463558) q[7];
u3(1.73223790244467,-2.36760129169666,-0.513166233030064) q[2];
cx q[2],q[7];
u1(2.04607458177651) q[7];
u3(-2.54354880068263,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.208357956628621,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.282655494530682,0.490977409272020,0.771334156313986) q[7];
u3(2.41379340719984,0.471437189426643,4.80043443687139) q[2];
u3(2.05341055446664,1.09759938888960,-3.63136589495820) q[12];
u3(1.79770494758870,-2.18073581588356,3.56077964521372) q[9];
cx q[9],q[12];
u1(0.555656270764138) q[12];
u3(0.0290139101511973,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.24518562255896,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.03995063034908,-0.206796250648018,-2.11396609468684) q[12];
u3(1.52273102492642,4.60298254634034,-1.42530374245722) q[9];
u3(2.58008914684416,2.43388741398725,0.327491642933501) q[8];
u3(2.48861189050017,0.568803291711659,-3.94901705821755) q[0];
cx q[0],q[8];
u1(1.49549024192045) q[8];
u3(-0.941924978645016,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.94015646661995,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.51948634861139,2.34005761821135,-2.53897194760351) q[8];
u3(1.00615756301499,5.05449744032901,-0.979858144120655) q[0];
u3(2.46798387814438,1.76181720911486,1.21194550375720) q[11];
u3(0.544903892292844,-2.63776443045091,-1.73253238443498) q[10];
cx q[10],q[11];
u1(1.50103715519553) q[11];
u3(-2.36234703832339,0.0,0.0) q[10];
cx q[11],q[10];
u3(3.41901378206491,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.23380208300379,-0.300769580934585,2.61209914178178) q[11];
u3(1.67670724894182,0.395147242896174,-5.70146329308787) q[10];
u3(1.95216132714526,0.466694905001111,-2.54146046842287) q[13];
u3(1.10453594109463,2.94094876455295,-2.90138411146527) q[6];
cx q[6],q[13];
u1(0.723042655401850) q[13];
u3(-1.55148942534616,0.0,0.0) q[6];
cx q[13],q[6];
u3(2.66714482741532,0.0,0.0) q[6];
cx q[6],q[13];
u3(2.16291086096816,1.55434692818950,-3.32109097147612) q[13];
u3(1.77806405023113,-0.777952300869880,0.276295468575673) q[6];
u3(1.70964316612447,1.10788595997961,1.67062920617871) q[11];
u3(2.10464649848504,-0.961442908175245,-0.964609449445425) q[9];
cx q[9],q[11];
u1(0.182597943063293) q[11];
u3(-0.368312880747514,0.0,0.0) q[9];
cx q[11],q[9];
u3(2.32519481782231,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.21350324838826,0.381579718521991,-3.30773456725453) q[11];
u3(0.846972703625309,-0.105552017820823,-1.44445655983482) q[9];
u3(0.547289378358049,0.820491515136475,-0.225060541519395) q[10];
u3(1.12255999834379,-2.37653188110011,0.875730663524236) q[3];
cx q[3],q[10];
u1(-0.356520855680370) q[10];
u3(-2.36094678900404,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.32212883144273,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.587195691505650,1.36208627317444,-0.428493436500567) q[10];
u3(0.900030080742858,-1.63279505930095,4.39722373253423) q[3];
u3(1.43675991258140,-2.61537236141169,0.818467052083293) q[8];
u3(1.88630312824682,-3.26946174535836,-1.04281076174118) q[13];
cx q[13],q[8];
u1(1.65228937008029) q[8];
u3(0.381750975089374,0.0,0.0) q[13];
cx q[8],q[13];
u3(1.03077276976439,0.0,0.0) q[13];
cx q[13],q[8];
u3(1.87706327262675,-2.30032647436140,2.33899699950736) q[8];
u3(1.23562878922065,-1.14645284336106,-4.18794045929156) q[13];
u3(0.918374509949195,2.05856367545575,-3.12723703597679) q[12];
u3(1.62103483011895,2.79146148687578,-3.35614075241727) q[4];
cx q[4],q[12];
u1(3.00375875461116) q[12];
u3(-1.04773199141960,0.0,0.0) q[4];
cx q[12],q[4];
u3(1.80110676183729,0.0,0.0) q[4];
cx q[4],q[12];
u3(0.862455034972939,1.84128185857399,-0.662384741148645) q[12];
u3(2.98455226952743,-4.06759483892980,-0.242402127242361) q[4];
u3(1.27985654461199,3.56914636795252,-2.49168103422042) q[5];
u3(1.87646114815011,1.47889892177220,-2.07101130609446) q[7];
cx q[7],q[5];
u1(1.39858662310159) q[5];
u3(-0.658917763066385,0.0,0.0) q[7];
cx q[5],q[7];
u3(-0.345636884953498,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.14855490146517,0.474387536286696,-1.40432711342529) q[5];
u3(1.72217434791334,-5.62456631830621,0.377670477154008) q[7];
u3(2.43677863351207,0.791171362639146,-3.36148696480672) q[2];
u3(2.06011227885556,2.56366182527434,-1.92745508354823) q[0];
cx q[0],q[2];
u1(2.22483399886945) q[2];
u3(-3.17629670465997,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.20545903838750,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.92080901756505,0.358521104869689,-3.50145330685851) q[2];
u3(1.35780720789401,-0.364591289656211,2.82102222753820) q[0];
u3(2.43898850591260,-2.87643860557711,-0.136437538377699) q[6];
u3(1.42812472318397,-4.30215567867975,-1.61170241755507) q[1];
cx q[1],q[6];
u1(4.09145164750294) q[6];
u3(-3.54215672973527,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.235057897775987,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.874474392119445,-0.261474017519064,4.47847109211305) q[6];
u3(0.956943629147991,-2.18641973770035,-3.30582504016009) q[1];
u3(0.791083195860115,-0.738668316473659,-1.49297509121008) q[1];
u3(1.04068898008857,-4.29127607361380,0.468722894985430) q[7];
cx q[7],q[1];
u1(0.573512162610358) q[1];
u3(-1.64706781041248,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.46370504110096,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.95353654828764,-0.411466282121691,1.58170500001281) q[1];
u3(2.21965289454414,-4.57335042310530,-0.561052117359029) q[7];
u3(1.22486638014207,-2.60450097503365,-0.395109219730683) q[13];
u3(0.847488265512352,-2.80558954264182,-0.452368072449598) q[5];
cx q[5],q[13];
u1(1.28296367368029) q[13];
u3(-3.10043587758544,0.0,0.0) q[5];
cx q[13],q[5];
u3(2.17166576191416,0.0,0.0) q[5];
cx q[5],q[13];
u3(1.33151098847510,-0.747421543073587,1.64300797851815) q[13];
u3(1.80107026249723,0.541400178773811,-1.11997153650512) q[5];
u3(1.67824634880457,1.04228371321994,-2.06261046142071) q[4];
u3(1.96546344797589,-4.19834256100765,1.63053921446137) q[6];
cx q[6],q[4];
u1(2.31474715569398) q[4];
u3(-0.0743582375332819,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.12597690320950,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.37714760136407,-0.764591686155759,-2.56143507596783) q[4];
u3(0.758248478956407,3.30479313875112,-1.43109580088083) q[6];
u3(2.25526344951272,0.687346286839359,-0.628456429470963) q[3];
u3(2.07299841613608,-0.0841929089986102,-3.94651514659383) q[10];
cx q[10],q[3];
u1(0.0914582949972176) q[3];
u3(-1.32613476049614,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.39459193748224,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.07648747210345,-0.101488471479352,-1.85114950775479) q[3];
u3(1.45530580628408,-0.677205358155256,0.632819370595619) q[10];
u3(0.896477509192508,3.12707330995959,-2.16159920902932) q[11];
u3(1.18840266814704,1.98319172566756,-1.94219074066048) q[9];
cx q[9],q[11];
u1(2.38792288889861) q[11];
u3(-1.87070615292784,0.0,0.0) q[9];
cx q[11],q[9];
u3(-0.0953301904551394,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.70595667258448,-0.574240603146142,-1.34279756054775) q[11];
u3(2.22037649274025,-0.621865623343328,3.06161541259281) q[9];
u3(2.63848049005755,-2.34967488591405,0.746822572415196) q[0];
u3(2.05615783188702,1.63730166649315,1.95754636388204) q[8];
cx q[8],q[0];
u1(-0.204410576177044) q[0];
u3(0.932422497454993,0.0,0.0) q[8];
cx q[0],q[8];
u3(3.72567206580775,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.54229854758069,-1.22387141875318,2.40989037410968) q[0];
u3(1.87564588292457,-3.08505164144206,-1.11975473226890) q[8];
u3(1.69231095162196,-4.52830118468832,1.50217720022417) q[12];
u3(2.48721083658308,-0.372693462316114,2.89634541044626) q[2];
cx q[2],q[12];
u1(3.70439467589132) q[12];
u3(-1.50056531655495,0.0,0.0) q[2];
cx q[12],q[2];
u3(1.71451395607281,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.58555597767035,0.843141557840613,-1.49789773761187) q[12];
u3(1.31412111938771,-1.33948789482429,-0.388562198313884) q[2];
u3(0.685290135308610,0.620390116897418,0.759435092940882) q[11];
u3(1.72720994985642,-0.194110025368641,-3.14907376479840) q[1];
cx q[1],q[11];
u1(1.49150509973799) q[11];
u3(-3.39010642506609,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.46739091350711,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.727134052384596,1.20849551811602,1.63864691354456) q[11];
u3(1.73069970599785,-2.14637616803295,1.85230908977469) q[1];
u3(2.90486422484497,0.339453392018029,0.411561815980916) q[10];
u3(1.39713811314452,-0.133829967491871,-3.19239054767878) q[13];
cx q[13],q[10];
u1(1.30427277850034) q[10];
u3(-0.550996645157019,0.0,0.0) q[13];
cx q[10],q[13];
u3(2.84520500368745,0.0,0.0) q[13];
cx q[13],q[10];
u3(1.71572078896902,3.40495572300207,-0.231879255141112) q[10];
u3(2.49966065324646,-0.584015417434298,-5.34421192233900) q[13];
u3(2.39188519817337,-0.0154293622455961,2.31677322583338) q[9];
u3(2.81255685757168,-2.99126381294832,-1.72639483846577) q[2];
cx q[2],q[9];
u1(2.06185787631574) q[9];
u3(-0.0781621981694183,0.0,0.0) q[2];
cx q[9],q[2];
u3(0.729462620375959,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.00660080422723,0.906702364878681,1.97202932385879) q[9];
u3(1.91131385502720,-0.283940346046442,-4.67423246146457) q[2];
u3(3.10648346862482,-2.56241156950962,3.71384267093184) q[4];
u3(1.12685089207730,-0.580319992050334,2.07701994275734) q[6];
cx q[6],q[4];
u1(1.54800384787562) q[4];
u3(-0.679154789782060,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.87957309238492,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.852343061368906,-3.36640399403771,0.192047683070033) q[4];
u3(1.08046047436034,-2.65935839812277,-2.87340028364251) q[6];
u3(2.74755964833450,1.59155093801059,-1.31126005689546) q[7];
u3(2.20572502527139,0.964032516885741,-4.04207190865311) q[8];
cx q[8],q[7];
u1(1.75587216542848) q[7];
u3(0.633077662578715,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.05775336989703,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.77013304276334,-3.72700187235365,0.735031163004020) q[7];
u3(2.27982096938265,-2.43001551996634,-0.574974577731653) q[8];
u3(1.53161842623349,-0.197226081932307,2.02697560939420) q[0];
u3(1.55024475294257,-2.27455467872746,-2.66030646473333) q[12];
cx q[12],q[0];
u1(2.27027434126538) q[0];
u3(-1.90915308626398,0.0,0.0) q[12];
cx q[0],q[12];
u3(0.342003581911587,0.0,0.0) q[12];
cx q[12],q[0];
u3(1.33927817371368,-1.74712828098253,1.42360734612337) q[0];
u3(1.68748329148550,-1.51333333193331,3.65229085644531) q[12];
u3(1.52941348949419,1.73844572924836,-3.04734677840752) q[5];
u3(2.33795112472524,-1.35784171279221,4.13644762735384) q[3];
cx q[3],q[5];
u1(0.439494606706359) q[5];
u3(-1.50793047418493,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.23518832360188,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.04071964254963,-1.71332063588137,-0.649129333766139) q[5];
u3(2.47275481911708,-2.14325952234094,0.764067651547293) q[3];
u3(2.15580810444672,-1.79604121584316,4.14193877708407) q[1];
u3(0.908947506724682,1.85197632372183,0.285083931012808) q[0];
cx q[0],q[1];
u1(1.26483132186918) q[1];
u3(0.00500212131546296,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.94322488800913,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.43572321469864,0.214823483683971,0.979185630500171) q[1];
u3(1.30779172653368,3.22715901549757,3.02563903841408) q[0];
u3(2.96013552609400,-3.01120693854285,2.67008278541788) q[5];
u3(1.19353001819293,0.705441901580043,0.986655489441037) q[11];
cx q[11],q[5];
u1(1.65507167624733) q[5];
u3(-2.54108424782397,0.0,0.0) q[11];
cx q[5],q[11];
u3(3.17234074098390,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.67299497668769,-1.14448412464840,3.47444231567661) q[5];
u3(2.97707776537684,2.92322072286363,-2.06885873347847) q[11];
u3(0.292267770215541,0.734444396894679,0.188522390360553) q[2];
u3(0.458774357169212,0.354892188428313,-1.99710321870956) q[13];
cx q[13],q[2];
u1(3.04026194081773) q[2];
u3(-2.55761416433335,0.0,0.0) q[13];
cx q[2],q[13];
u3(1.65559408497072,0.0,0.0) q[13];
cx q[13],q[2];
u3(0.189417627874159,1.50261021672186,1.05971432615432) q[2];
u3(1.81393945841967,4.45404719312893,-0.323056885351208) q[13];
u3(1.61746025064138,3.88624193001713,-1.44292086618208) q[4];
u3(0.338167217889837,0.515921777735556,-0.511477551301031) q[6];
cx q[6],q[4];
u1(4.26275352071990) q[4];
u3(-3.44836114039083,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.191589846995907,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.04367136455393,0.157895918956189,-1.12595633421433) q[4];
u3(1.05965994512596,-0.812320034955236,4.39025441736651) q[6];
u3(2.09604967905180,1.42082909689685,-3.55803830827381) q[7];
u3(0.682784569433071,-0.921454000311585,1.96359059423006) q[12];
cx q[12],q[7];
u1(-0.474616839296927) q[7];
u3(-2.07985118330162,0.0,0.0) q[12];
cx q[7],q[12];
u3(1.47946128721703,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.81730475879349,0.700681178627734,-1.34673442210199) q[7];
u3(2.15889600268539,-1.83415124091919,-3.89634781427820) q[12];
u3(0.503062345283355,-2.49137447701695,2.32519123112991) q[3];
u3(0.747985869976534,-3.78254394632270,1.99145194465480) q[10];
cx q[10],q[3];
u1(1.70540813278442) q[3];
u3(-2.18583270872715,0.0,0.0) q[10];
cx q[3],q[10];
u3(3.35742644359697,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.26887165347722,-3.95487998666839,0.685744055941778) q[3];
u3(1.66420930296844,-3.35278068089935,-1.70699476898846) q[10];
u3(1.13755007528205,0.178836517192332,1.13142655426070) q[9];
u3(1.33463325800019,-1.34833200119989,-1.54185118954735) q[8];
cx q[8],q[9];
u1(1.64667444744670) q[9];
u3(0.182751834323771,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.324027879986089,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.81567047097105,-2.68835755275334,0.547696129347743) q[9];
u3(2.49538191215244,1.95264436731208,0.755693390890098) q[8];
u3(2.07885811880028,2.76119644057094,-0.541690375015289) q[6];
u3(1.18251748460503,0.974822227756615,-1.52565124142706) q[8];
cx q[8],q[6];
u1(-0.102046442473273) q[6];
u3(-1.21026937522471,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.59325710098355,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.02460743421035,0.521195135450983,-3.70944346905324) q[6];
u3(1.91620647576461,-2.09539860167152,-2.18749334262026) q[8];
u3(2.81962251379952,1.48962349047148,-1.89184103549055) q[11];
u3(2.68916900508615,5.19631886378885,-0.173746661970717) q[0];
cx q[0],q[11];
u1(3.66433283206298) q[11];
u3(-1.16435764860980,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.45513395008645,0.0,0.0) q[0];
cx q[0],q[11];
u3(2.99670893988718,-1.54387691105746,2.02298140744999) q[11];
u3(0.870375426812469,2.60483108569755,-2.82817050984182) q[0];
u3(1.51827151913771,-1.51711518367554,-0.925091361428472) q[4];
u3(1.99380690056865,1.52930177782488,-4.69692133015750) q[10];
cx q[10],q[4];
u1(0.459547186425127) q[4];
u3(-0.826571527855113,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.12597506445008,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.161818413322129,-1.73068632372887,1.71560825979469) q[4];
u3(2.78126384348769,-1.20988894216342,-4.54195818506057) q[10];
u3(0.816517085623006,2.10114262876384,-3.28383215674258) q[7];
u3(1.43115949403063,3.90377819977021,-2.36317675529088) q[1];
cx q[1],q[7];
u1(0.530531733203402) q[7];
u3(-1.53067846192864,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.265674101433892,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.20630028147457,3.39841652835708,-2.16047607156605) q[7];
u3(2.71629852052347,-0.00915583128247155,-4.41362946285268) q[1];
u3(0.874207669140776,1.91584857464793,-2.78208860256409) q[3];
u3(0.554150420452581,-0.0551112600378498,-1.29451729228742) q[13];
cx q[13],q[3];
u1(2.93737712680495) q[3];
u3(-2.03511628494048,0.0,0.0) q[13];
cx q[3],q[13];
u3(1.25874467051243,0.0,0.0) q[13];
cx q[13],q[3];
u3(0.608160099831623,-1.32366823231286,1.54542736549556) q[3];
u3(2.82768426495907,-4.04017453865803,-1.74747739170917) q[13];
u3(0.472961780635183,-2.18989907949329,2.68823351680735) q[9];
u3(0.776413241042936,0.290214274805394,-2.65048870113849) q[2];
cx q[2],q[9];
u1(1.07919376587104) q[9];
u3(-3.44885922314625,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.45293030973434,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.563599220807325,-0.934490817229900,4.51359348480040) q[9];
u3(2.09849193765266,-1.35150559469808,1.95835601038831) q[2];
u3(1.77657557373934,1.61116848773607,-0.291445782396025) q[5];
u3(1.19658810026675,0.0756326689525393,-3.45685823023254) q[12];
cx q[12],q[5];
u1(-1.06603274895540) q[5];
u3(0.498374573626087,0.0,0.0) q[12];
cx q[5],q[12];
u3(3.37144608609410,0.0,0.0) q[12];
cx q[12],q[5];
u3(2.24415139230984,0.197943907934381,-0.165815762270564) q[5];
u3(0.457509358934574,0.691023387281182,-1.43811677427844) q[12];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13];
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
measure q[13] -> c[13];
