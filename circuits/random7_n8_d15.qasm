OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(2.28880497065046,1.56357662726013,-3.04054423239718) q[0];
u3(1.49417900461491,2.25319479137275,-3.39194986939084) q[6];
cx q[6],q[0];
u1(2.21931911514384) q[0];
u3(-1.65217429142463,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.27460507120046,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.29010454467842,-3.12931261759905,2.70140593821840) q[0];
u3(2.38057447301086,-2.49431224058552,1.70190690550105) q[6];
u3(1.53479412146423,-2.06595930397740,0.171949600762927) q[5];
u3(1.56615722901114,-3.82086097741282,0.760409620779908) q[3];
cx q[3],q[5];
u1(2.53899743477211) q[5];
u3(-2.23336755653940,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.60812873431457,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.60261244266468,-3.32757611520057,1.80045554489730) q[5];
u3(2.13250561174320,-0.429255703029276,5.69089706955634) q[3];
u3(1.46347912884184,-0.322250429221243,-1.27305965114421) q[2];
u3(2.00670805335253,-3.47051304777614,2.14636982586531) q[1];
cx q[1],q[2];
u1(1.89925259779900) q[2];
u3(0.703536353049629,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.46198499945687,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.19312570225607,-0.351103867955963,-1.27780492283112) q[2];
u3(0.511355906365657,-2.58879317955798,-3.01081633411812) q[1];
u3(1.69841641424671,2.32953387963638,-2.25926747649051) q[4];
u3(0.489142682648034,-2.60507014053811,1.56270014354873) q[7];
cx q[7],q[4];
u1(2.36709824458087) q[4];
u3(-1.70695230875909,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.33380672968328,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.466448575733698,-0.00176354309279136,2.20915762563464) q[4];
u3(1.09687588781632,0.0689120811407271,-2.83185393836082) q[7];
u3(1.30239548887154,-1.37969245020197,1.31296487312524) q[5];
u3(0.968825150425938,-1.68307600072499,0.107057285851508) q[0];
cx q[0],q[5];
u1(1.59957962824819) q[5];
u3(-0.539353091676339,0.0,0.0) q[0];
cx q[5],q[0];
u3(-0.116153917691116,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.84222194355592,1.26922652517609,-4.00465525383592) q[5];
u3(1.88525321289156,-2.94600777110429,-3.29393405060431) q[0];
u3(1.89083494287155,1.86756530818567,-2.73471935353815) q[4];
u3(2.24916319391720,1.61279381376376,-4.18803450804266) q[3];
cx q[3],q[4];
u1(3.99095279539018) q[4];
u3(-4.26163067622337,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.424487206279900,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.96919981081053,-3.25262165364688,2.35114538716866) q[4];
u3(0.449507164423762,-1.72807496327773,-0.594608392026726) q[3];
u3(0.833135984829437,-0.373260478184666,-0.721690743570171) q[1];
u3(1.42820551116445,-3.01295098881641,1.01683156783907) q[7];
cx q[7],q[1];
u1(2.10986716915162) q[1];
u3(-2.80637390997449,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.49083747450980,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.58111203368829,0.305631208671598,1.13146214073262) q[1];
u3(2.40031299032246,3.09145746532694,1.87256951098279) q[7];
u3(0.541382440525510,-1.78822543793722,0.615604081621328) q[6];
u3(0.813431898511342,-3.38034411419007,1.89301787778620) q[2];
cx q[2],q[6];
u1(1.87370882389626) q[6];
u3(-2.35707002606226,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.0347174661105376,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.75971150604195,1.75970143680213,-1.27403594052383) q[6];
u3(1.94427956334181,-0.596486851424255,1.47623315383218) q[2];
u3(2.23055503414961,1.46612417940932,-3.30950287949995) q[5];
u3(1.35857230529495,1.88588979317035,-2.32868464840990) q[2];
cx q[2],q[5];
u1(2.83189490904674) q[5];
u3(-2.63078650337626,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.20005071161181,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.17529074408314,1.11588082995011,-2.48609017660125) q[5];
u3(0.718752267039846,4.22501906209432,1.00815100526200) q[2];
u3(2.08757134494579,0.840596623073223,2.00382781575481) q[6];
u3(1.98244044783785,-1.20545069060826,-1.16357678623415) q[1];
cx q[1],q[6];
u1(1.43346963771197) q[6];
u3(-0.702845981522543,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.11119837820387,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.617240511228073,-1.28419194572554,0.317278199960542) q[6];
u3(2.16598728757817,-4.48861881586507,1.14053972425317) q[1];
u3(2.51823846980551,3.72786338331913,-1.56743656197656) q[7];
u3(2.14722150044842,1.31445469824107,-0.940793520788400) q[4];
cx q[4],q[7];
u1(1.75761041806009) q[7];
u3(-2.15901780954818,0.0,0.0) q[4];
cx q[7],q[4];
u3(0.489543105052167,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.91059251701843,1.87368833642396,0.818253741860241) q[7];
u3(0.170900053059850,-2.81175571069561,-1.34762607298993) q[4];
u3(1.04406017586055,-0.761596005062097,3.81123633884962) q[3];
u3(0.967141785756131,-1.13601529634434,0.793980909778454) q[0];
cx q[0],q[3];
u1(-0.00792601552269989) q[3];
u3(-1.58768325713020,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.614270486148209,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.924955267261313,1.97861568092537,-3.32333353646339) q[3];
u3(1.08509981125328,0.921545357609002,1.89231166254944) q[0];
u3(2.58678044294119,1.64126127185059,-4.10195758117164) q[2];
u3(0.697606007486820,3.58129292561863,-2.46297318060397) q[5];
cx q[5],q[2];
u1(0.988653564275678) q[2];
u3(-3.29636616225875,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.88828067227788,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.21823469232739,0.306305490528838,-2.84310664916241) q[2];
u3(2.76932755384721,-1.36489410006756,1.91611718135199) q[5];
u3(1.84697234808762,1.78300860687797,-0.243590742389932) q[1];
u3(2.58490204971357,1.78058344416925,-2.36965755247811) q[0];
cx q[0],q[1];
u1(1.51075243326083) q[1];
u3(-0.177293574828698,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.38956207519450,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.42027815568730,3.47736690318530,-1.23887400878131) q[1];
u3(1.48059996749616,0.344206685725335,1.77776863273516) q[0];
u3(1.50672951389321,0.371184014307487,-1.26277485769260) q[6];
u3(2.30155883623646,-4.12331242078554,0.779433918816065) q[7];
cx q[7],q[6];
u1(-0.253844011559107) q[6];
u3(-1.53750737658885,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.25700169340607,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.42743487729541,-4.89126410495649,1.25469269314794) q[6];
u3(1.50137357629786,1.11034789299832,-1.85270481039226) q[7];
u3(1.63229319509701,2.35868820849117,-2.70993874414672) q[4];
u3(1.78103382145672,-3.69791736828590,2.46661127033885) q[3];
cx q[3],q[4];
u1(3.44077344835711) q[4];
u3(-1.90493945344366,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.45873952257241,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.316482710666583,-1.17771618607024,0.904670168700415) q[4];
u3(1.97449858181753,-2.78059951707810,-0.447900338669082) q[3];
u3(2.52376584444220,0.589593211957395,1.60332215335413) q[3];
u3(1.03263287921780,-5.33496592281978,0.0644775925573140) q[6];
cx q[6],q[3];
u1(3.96564844860560) q[3];
u3(-3.52275882424285,0.0,0.0) q[6];
cx q[3],q[6];
u3(-0.878483847187972,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.69527764973716,1.24149451026271,-3.85057804680246) q[3];
u3(0.528561002116453,-0.182866628979793,5.41984743574072) q[6];
u3(1.26603585501535,-0.446177726014755,1.43723753784807) q[2];
u3(1.20243365048058,-0.530857836871690,-1.33811545398310) q[1];
cx q[1],q[2];
u1(1.01393041263751) q[2];
u3(-1.52851957122823,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.00854215752176102,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.33105761498236,-2.68210264239092,2.86838903444272) q[2];
u3(1.58468178153172,-2.06943581392645,2.22997525997970) q[1];
u3(0.160226427554378,-1.95420021785322,0.0951392517748336) q[5];
u3(1.98011355732031,-4.93380685783930,-0.130393459015925) q[7];
cx q[7],q[5];
u1(1.49458055426707) q[5];
u3(-2.55006667457670,0.0,0.0) q[7];
cx q[5],q[7];
u3(3.33721967597189,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.06142445092254,-0.361535452934017,0.155888834860847) q[5];
u3(2.95876838641592,-1.07264549049071,0.237501968168801) q[7];
u3(1.13534961263664,3.82533689217039,-1.31077320205826) q[4];
u3(1.92014484727226,2.40529028216032,-0.311331536630069) q[0];
cx q[0],q[4];
u1(2.41233632482981) q[4];
u3(0.104736135829800,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.05274272597157,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.02384567648592,0.140429128908962,1.10122275664599) q[4];
u3(1.60780860790946,2.63634568519777,-1.05960234719264) q[0];
u3(1.13410062532872,4.04371982624069,-1.67967883469274) q[5];
u3(0.784846568738177,0.400511920690357,-0.215669017287884) q[3];
cx q[3],q[5];
u1(0.863731592691275) q[5];
u3(-0.414929260829659,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.00725661466077,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.943116531018446,-1.28350309640257,-1.93229827097464) q[5];
u3(0.481250415615465,1.38194043497386,2.78344174738532) q[3];
u3(2.53291157350119,1.10755727409951,1.71755531102364) q[6];
u3(1.41549650020748,-0.909995544045598,-2.33119395352702) q[1];
cx q[1],q[6];
u1(-0.214371072003984) q[6];
u3(-1.91188257988085,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.580233416601931,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.607622871257142,-2.13159675433034,1.12525669189457) q[6];
u3(2.36800124359610,5.26028310014462,-0.0969161146563859) q[1];
u3(0.399221259135586,2.51731648076082,-2.92168507334626) q[7];
u3(1.00505884505004,0.284515660947027,-2.36967399622105) q[4];
cx q[4],q[7];
u1(0.810795334736806) q[7];
u3(-0.195947579257205,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.78181898757350,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.25520805859910,3.11751390891412,-1.09985261757458) q[7];
u3(2.52693324118924,-2.84429179563660,-0.561529655673250) q[4];
u3(0.552992864511027,2.22715147725343,0.665462933188184) q[0];
u3(1.50537274648071,0.582030341875428,-3.89562600789345) q[2];
cx q[2],q[0];
u1(1.17067904779090) q[0];
u3(-1.40584836870920,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.53382636242684,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.385433315330181,-2.16614203894667,2.57485333520980) q[0];
u3(1.81102643557292,1.46736193671011,-0.979476465085598) q[2];
u3(1.60888435426012,1.85195023531792,-2.87208613101062) q[1];
u3(2.02281544841469,1.80128743825078,-4.47173031946954) q[5];
cx q[5],q[1];
u1(3.25541572024972) q[1];
u3(-1.89068488969216,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.673339134903724,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.892105229138294,-4.02446432800251,1.22569167007856) q[1];
u3(1.40936841542133,0.156437904183697,4.08715678894700) q[5];
u3(1.36041397044688,1.67198174504440,-2.45770115965883) q[0];
u3(2.33369383823191,1.82723597774759,-4.11483783171153) q[4];
cx q[4],q[0];
u1(1.08949742924417) q[0];
u3(-3.17635963160293,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.63372180083802,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.33703333525576,-1.57037275869056,0.553400824841293) q[0];
u3(2.11089720539344,0.710726088841655,4.68731357406393) q[4];
u3(1.22559472725764,1.21848998355441,-2.55189154465764) q[2];
u3(1.94833042184128,-2.05467476160602,3.74998124867687) q[3];
cx q[3],q[2];
u1(2.50246195892744) q[2];
u3(-1.84570052376347,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.40717699539687,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.40835720908777,-1.19605691753681,0.0248521063085264) q[2];
u3(0.926554830393586,1.96079064565214,3.97977071122374) q[3];
u3(0.740909632190976,-1.76090014086203,2.52570619768619) q[7];
u3(0.818482919363237,-2.34978779169166,1.14498378462087) q[6];
cx q[6],q[7];
u1(1.64059475963946) q[7];
u3(-0.0404010632572880,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.424188981895131,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.84609840882769,1.59402974940764,-4.11396493794803) q[7];
u3(0.804443515064910,3.64410869464956,2.33320495032820) q[6];
u3(1.17993089484349,-0.235794182513115,1.37347985005690) q[0];
u3(1.84486517020252,-0.817186401864843,-2.48715918867520) q[3];
cx q[3],q[0];
u1(-0.375072063868019) q[0];
u3(1.02257694585276,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.07447238172251,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.53589158202994,-2.28038447493904,3.86452015766236) q[0];
u3(2.10749448431211,-5.25977753711450,0.932852550850988) q[3];
u3(0.380253546999570,1.76871190444574,0.258227527762120) q[2];
u3(1.20931103681420,0.117106289042256,-3.32463041744083) q[6];
cx q[6],q[2];
u1(1.68237086792589) q[2];
u3(-3.23968168530910,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.27481806363857,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.507930779237924,-3.13830685331679,3.11273354018260) q[2];
u3(1.02770202413414,2.59672283271331,1.74209858607626) q[6];
u3(1.97032522414320,-2.31743920702295,-0.450283734229863) q[7];
u3(2.34555267570905,-4.58224268937040,-0.984532316750568) q[5];
cx q[5],q[7];
u1(0.405445689943930) q[7];
u3(-0.996906511898005,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.72124755531072,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.57418527571760,-2.15956433043075,3.68732383574352) q[7];
u3(0.389018997274496,0.497375615816816,3.50573033034819) q[5];
u3(0.428000590638512,1.28136483704953,-1.14463437075044) q[1];
u3(0.980979447035544,-3.07155394009259,1.69263635696244) q[4];
cx q[4],q[1];
u1(-0.143846490414950) q[1];
u3(-2.42652138909405,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.27551068614845,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.698984503826593,0.526882723763965,-1.75510139692287) q[1];
u3(2.69820911923027,-4.06069018802373,2.09691638860096) q[4];
u3(1.66690768192796,0.495264247283500,-3.53239453662799) q[6];
u3(1.37199982521379,3.26114938019719,-2.64942249354821) q[7];
cx q[7],q[6];
u1(2.21165759343340) q[6];
u3(-0.192609262185236,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.28634609253363,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.809583609332592,2.77776022462013,-2.60050053149768) q[6];
u3(0.659445584917283,-2.88917016967753,0.970244287371058) q[7];
u3(1.21159118550103,2.76328910174607,-1.51547101262095) q[0];
u3(0.829576664940021,0.995074793815250,-1.90395159869713) q[2];
cx q[2],q[0];
u1(0.942659381225710) q[0];
u3(-1.23628634431172,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.86225234919698,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.347663051344324,-3.10573394683506,1.88078449437939) q[0];
u3(2.40554983147443,1.41262481178589,0.203730364711251) q[2];
u3(2.24698551999016,1.58665092569754,-3.36464130217823) q[5];
u3(1.22863700608241,1.91863488847394,-2.36098501568832) q[1];
cx q[1],q[5];
u1(0.875740424374602) q[5];
u3(-0.358887904279988,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.46954439453973,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.55222230525407,-1.36348059353306,3.17340919141902) q[5];
u3(0.757992632061217,-0.739941588841210,-2.16950523950063) q[1];
u3(0.461063035697648,-1.77809903584635,1.63489810113076) q[3];
u3(0.377737373030605,1.39581372808382,-1.90215940015616) q[4];
cx q[4],q[3];
u1(2.56907475745400) q[3];
u3(0.0962963476583245,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.25902209124942,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.47066196573622,1.06600340154633,-1.00644912838355) q[3];
u3(0.608895761345916,1.64503581607475,-2.64686197897935) q[4];
u3(2.07124704140047,1.45103105052674,0.598584423980991) q[7];
u3(2.43458063874700,-0.123409368539560,-3.34172702774801) q[5];
cx q[5],q[7];
u1(0.813908892675587) q[7];
u3(-0.206451925146432,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.67189717246096,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.63308694056555,4.00903968761900,-2.09224748658085) q[7];
u3(2.06837490707594,-3.49561129970160,-1.79326638325925) q[5];
u3(1.52414269564930,0.628371239340360,1.18175790536542) q[1];
u3(1.40346330689438,0.545492450771858,-3.25897761085298) q[3];
cx q[3],q[1];
u1(3.63270975126802) q[1];
u3(-1.33385400243615,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.36845855271825,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.842352171257083,3.33960203660069,-2.59767152399502) q[1];
u3(0.996266642023756,-3.32837889600694,-1.62311804508894) q[3];
u3(1.20870793120384,0.676614704962203,1.79252715227881) q[4];
u3(1.20746586364308,-1.71357527289761,-2.07072913135682) q[2];
cx q[2],q[4];
u1(1.54531029672930) q[4];
u3(-3.09737058210969,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.546200921028830,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.50714987285061,-2.64758292532042,0.622760704871446) q[4];
u3(0.153588057935226,-2.30648112624532,-1.06545988752708) q[2];
u3(1.07165859010598,1.92417786563359,-2.69244549933252) q[6];
u3(2.00510813947883,3.95402413476788,-2.24940717020735) q[0];
cx q[0],q[6];
u1(-0.644149674255509) q[6];
u3(-1.84317752517984,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.01656192820992,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.05064240588974,-0.984246451556399,1.88784342543723) q[6];
u3(1.12619918635208,2.17365444947476,1.63602136140836) q[0];
u3(1.69564592904823,0.195193998162566,-1.60437651545531) q[6];
u3(0.658833829271386,1.36689669401219,-3.38429333819062) q[0];
cx q[0],q[6];
u1(-0.344980531715424) q[6];
u3(-1.82063590991062,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.819889028942926,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.73284557357759,2.88070019547088,-2.38043759969461) q[6];
u3(2.05870997411808,1.12581578744710,3.86040134775525) q[0];
u3(0.940916112752881,-0.444309314416603,1.48697002968774) q[5];
u3(0.609464208359724,-0.0955845692169041,-1.51945376579467) q[4];
cx q[4],q[5];
u1(1.61466893288898) q[5];
u3(-3.57088896233757,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.07013100237848,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.91293309821667,0.280942968137548,0.262824735419249) q[5];
u3(2.07312149641873,0.273823077985758,5.58173948835736) q[4];
u3(1.69435469991005,0.0952713385640325,-2.30278465125461) q[1];
u3(2.64831505842351,0.293735172121079,-5.04922488299403) q[2];
cx q[2],q[1];
u1(2.12600216351558) q[1];
u3(-2.82968782871353,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.15317337152164,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.763938865880888,3.53606962120955,-1.92098945781186) q[1];
u3(1.57264647408418,2.42113342011571,-1.35671128453625) q[2];
u3(0.885487581056274,1.62919421152846,-1.95066511992492) q[3];
u3(0.208364015480097,-2.52970209910166,0.715870567949384) q[7];
cx q[7],q[3];
u1(2.29821716844829) q[3];
u3(-3.02630266764005,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.50626291987088,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.63421533979120,-2.97019519805275,2.52652264045678) q[3];
u3(1.31934421914679,4.86148552725464,1.32984029419098) q[7];
u3(1.17051961549463,0.574891530271319,1.52369091370238) q[3];
u3(1.40090794808378,-1.46927258714467,-1.17945281269949) q[5];
cx q[5],q[3];
u1(3.41982039348662) q[3];
u3(-1.60298638771996,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.71884802843857,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.43288475058199,4.62786069115548,-1.40400966878889) q[3];
u3(1.60844785155460,3.01521447070766,2.91920949661441) q[5];
u3(1.98697268966772,-1.59747240378228,-0.453056973653604) q[2];
u3(0.194214546904513,-0.0311343192774585,-4.78144980601377) q[4];
cx q[4],q[2];
u1(1.56077999605377) q[2];
u3(-2.77988679962304,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.325917272259518,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.64572153281400,-0.286042708990499,1.34043056520429) q[2];
u3(2.34739794601673,1.42476942310247,3.26438896266502) q[4];
u3(1.99797514865527,-0.848110165070838,0.698829744398602) q[0];
u3(1.99616146514590,-2.12814295708292,-2.32236721766217) q[1];
cx q[1],q[0];
u1(0.0287025145310675) q[0];
u3(-1.41725770401144,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.95776438159494,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.11948442945405,0.876771381058671,2.02274997074585) q[0];
u3(2.61840291973248,-0.139319779702607,1.26929512266226) q[1];
u3(0.978824905473376,-2.06501942942076,0.929470345507503) q[7];
u3(0.393295508338651,-3.30839818729239,1.32086031595290) q[6];
cx q[6],q[7];
u1(-0.0448441388806775) q[7];
u3(-1.52166087002045,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.74638225788817,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.646079583279208,3.38906685337669,-1.61560504839099) q[7];
u3(1.58877263401615,4.66898952050877,-0.680666968830682) q[6];
u3(1.69753753792272,1.96529993959270,-0.938181846237191) q[4];
u3(1.52769049723837,0.847215919906601,-3.24874530195325) q[6];
cx q[6],q[4];
u1(2.18679405807491) q[4];
u3(-2.79322521247722,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.259958138944178,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.61418378460664,-2.13971924547550,2.39993409882900) q[4];
u3(1.27190265946318,-2.15489531361220,0.764127951466100) q[6];
u3(2.23685044639982,-2.25901069731121,-0.0977229708267848) q[2];
u3(1.12301638719673,1.62194633661843,4.38876243351846) q[1];
cx q[1],q[2];
u1(1.59273189909508) q[2];
u3(-2.60102507011026,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.34609015483777,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.779977761421350,-1.40494024011108,0.463608994587384) q[2];
u3(2.18470120160021,-0.640986651907924,0.0435422804525817) q[1];
u3(1.65426338368647,2.52007611446201,-2.53869002119293) q[0];
u3(0.783093094721264,-3.10961233687766,2.47022126136254) q[5];
cx q[5],q[0];
u1(2.24948424896013) q[0];
u3(0.248421818871046,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.24408042452435,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.57573671554459,2.25338640228036,0.645418092708645) q[0];
u3(2.72170200916525,0.200748856738411,-1.60378996776064) q[5];
u3(0.660441545034159,1.45386618433172,0.631409509875814) q[7];
u3(1.32196153771862,0.677607122031927,-3.63633263630062) q[3];
cx q[3],q[7];
u1(-1.48636683517612) q[7];
u3(-0.617601197993497,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.60228726547046,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.85493719026655,3.38934796243944,0.616458436234121) q[7];
u3(2.34090005622421,-2.24641574171186,-1.17967369222998) q[3];
u3(1.42994724817197,-0.128447623456979,1.98210294972457) q[5];
u3(1.12067520136312,-1.06316231115678,-1.41759067336644) q[1];
cx q[1],q[5];
u1(3.30169627894806) q[5];
u3(-4.30250389939019,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.467340756374973,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.477945982817217,2.19233214006837,-2.48017635833295) q[5];
u3(2.19431588424987,1.45288691738506,1.28637752929921) q[1];
u3(0.317419235434735,2.95935894649575,-2.83792391620586) q[0];
u3(1.14944399692721,0.438333199793353,-1.03908985578335) q[4];
cx q[4],q[0];
u1(1.72024901887744) q[0];
u3(0.0160997113219630,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.71817765460944,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.17405483095219,1.89890372176087,0.596655956084568) q[0];
u3(2.22099251337212,5.47298212574616,-0.795804904531512) q[4];
u3(1.46734867462534,-0.0105816889869815,1.05441422195658) q[2];
u3(0.938739452215183,-3.21337747484219,-1.34899477891508) q[3];
cx q[3],q[2];
u1(1.87599952071434) q[2];
u3(0.445993909874400,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.63296434345097,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.24009949038140,-3.07748907940697,2.63889085975589) q[2];
u3(2.19461648949165,0.701156219106532,-1.56442579551738) q[3];
u3(0.502646300874799,2.21454791910017,-3.86317731978501) q[6];
u3(1.34135905570919,3.50577078602720,-2.69410386027624) q[7];
cx q[7],q[6];
u1(1.02111154180450) q[6];
u3(-3.25093667115239,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.65858439732423,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.80556993062022,-2.36669338640250,2.93300343161666) q[6];
u3(1.87597271332282,4.87270219919512,-1.01137719479878) q[7];
u3(2.21772217821919,-1.05284493503624,3.85443802273643) q[4];
u3(1.25735944129149,1.11981387612517,1.33646474237354) q[5];
cx q[5],q[4];
u1(1.61262749737661) q[4];
u3(-0.00371766516957850,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.797592940851276,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.58367817227408,-0.190556458336459,2.48447029450803) q[4];
u3(1.77810143851005,0.530359846459285,-2.56472448572841) q[5];
u3(1.02486690311944,-1.35951056855216,-0.697283161586547) q[3];
u3(0.933939434881074,-3.64753216785978,-0.566367849854191) q[1];
cx q[1],q[3];
u1(-0.353512358695720) q[3];
u3(-2.37571473506162,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.77889822528813,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.12418991947061,1.17744200941713,-4.73882284678450) q[3];
u3(1.88880950598173,-3.57319564253129,-1.82719674038777) q[1];
u3(1.08778582859136,0.178627141472937,2.10316398519204) q[7];
u3(1.66864653150245,-0.242737913730988,-2.70379463879790) q[6];
cx q[6],q[7];
u1(1.74507959732544) q[7];
u3(0.605820798084632,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.05768932289042,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.816289414362393,-0.589882206563336,-0.853627970814995) q[7];
u3(2.28166697942843,1.41017134234023,2.03498838849812) q[6];
u3(1.38172134503100,0.995916567303518,-2.61993728016918) q[0];
u3(2.02134636533031,2.00821493570444,-3.93283069807355) q[2];
cx q[2],q[0];
u1(-0.219044621882128) q[0];
u3(-2.12949549475159,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.903025837407414,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.109888028019760,2.14112693644321,-3.24474409118335) q[0];
u3(1.60072229135466,-3.65354597867875,0.335749377648606) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
