OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.19787875511803,0.498729195855994,0.953112813116890) q[7];
u3(1.85942600886196,-0.828265566265124,-2.78816146078824) q[1];
cx q[1],q[7];
u1(3.21145464488733) q[7];
u3(-3.98240526995797,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.368327085115756,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.49584612067443,2.85986616571655,-1.78455354399922) q[7];
u3(2.02949639039013,-0.378194429245470,5.79290251234062) q[1];
u3(0.810980588709180,2.41791275171570,-2.09361659884253) q[4];
u3(0.0808939507238932,0.720858714875526,-1.81628308216266) q[9];
cx q[9],q[4];
u1(1.56482113047480) q[4];
u3(-3.14017928119023,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.711340568407868,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.08617719981491,3.83604868146343,-2.37332254062609) q[4];
u3(2.53805078624872,-1.72450081482170,3.09333560788224) q[9];
u3(2.12618379778383,1.71769235396681,0.636488953157415) q[0];
u3(1.82286556887070,0.0225635160148070,-3.12591946624862) q[5];
cx q[5],q[0];
u1(1.73565959679436) q[0];
u3(-2.23309614615408,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.37381602862219,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.477638203369074,-0.133701418457872,-0.561545161650621) q[0];
u3(2.25724864987791,-4.10346899589382,1.84039525596941) q[5];
u3(2.95794710515285,-3.40762746012916,2.78075923144652) q[3];
u3(1.43312686092889,3.29270274875112,-2.78635509600129) q[8];
cx q[8],q[3];
u1(1.28982321700499) q[3];
u3(-0.0553083996194506,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.23845338944251,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.72898258576454,-1.22882023100434,0.861653949582788) q[3];
u3(1.20963994163294,-2.04904814744451,3.04587947729776) q[8];
u3(1.45564600115485,-2.31000098526268,-0.368837293947564) q[6];
u3(2.45529908402828,-3.63410339118079,-0.155606265621434) q[2];
cx q[2],q[6];
u1(1.37398109221625) q[6];
u3(-3.44481070889075,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.21779441858899,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.61737950187735,0.718581953308506,0.939665604017616) q[6];
u3(0.900760852839306,-1.52814135622528,2.71167358031131) q[2];
u3(2.09103432871476,-0.110647479847088,2.24724705871090) q[4];
u3(1.89946215378669,-3.07252265731332,-2.71816265898199) q[7];
cx q[7],q[4];
u1(1.38514337392597) q[4];
u3(-3.49778443378248,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.19477810507667,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.40032009750856,-0.0136615788820142,-2.10055706064238) q[4];
u3(1.27754697447448,-1.10685786328766,-2.17057954962707) q[7];
u3(2.10947569802048,-0.169857995638459,-1.18088004001878) q[9];
u3(0.699933691537436,0.150066867088107,-4.76948834536403) q[8];
cx q[8],q[9];
u1(1.55089935581995) q[9];
u3(-0.447853938407192,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.56183809983068,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.49407085823025,-0.573505313902908,4.95010112183440) q[9];
u3(1.19469491423272,4.75142072640788,-0.0543869374931671) q[8];
u3(1.95660523878953,-0.328555571854800,1.28133147838777) q[3];
u3(1.56760499614201,-0.814903569679584,-1.04009501018017) q[5];
cx q[5],q[3];
u1(4.21926083549365) q[3];
u3(-3.71170351898347,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.164750467762134,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.42262261195352,1.00045833518285,-0.0550030862696236) q[3];
u3(1.79238998376030,3.74765359702904,1.66238035708693) q[5];
u3(2.57067873014855,0.787917144079604,-0.693379132601611) q[2];
u3(1.77408123777927,1.13110889489604,-4.52266120326483) q[0];
cx q[0],q[2];
u1(2.19999619305795) q[2];
u3(0.296521493505516,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.16402725394998,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.286553158209258,3.91013554578675,-0.0110426317074990) q[2];
u3(1.21117240460623,-0.695626730535523,1.10393147048750) q[0];
u3(0.999304378311846,-0.962401558694434,1.18483734145438) q[1];
u3(0.872587187831340,-2.03988117360353,-0.0623931320588650) q[6];
cx q[6],q[1];
u1(2.04474451826680) q[1];
u3(-1.76094329294194,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.17568593271003,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.35920207525062,-2.27480875430957,0.942509749004720) q[1];
u3(0.433116557557253,3.86351006833380,0.456780660988133) q[6];
u3(1.08584100793829,0.238000209321881,1.22146765736303) q[9];
u3(1.57287115423705,-2.82120855473883,-1.14058857493033) q[1];
cx q[1],q[9];
u1(0.500014332113435) q[9];
u3(1.28728686366633,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.89017357214111,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.08577093887496,3.37691351895953,-1.54937986089373) q[9];
u3(1.22195248352719,0.0735781601833099,-0.714660390529978) q[1];
u3(1.79029478926609,0.339683749673614,1.56399531380067) q[6];
u3(1.18055043032693,-2.51940525076371,-2.11352651153076) q[0];
cx q[0],q[6];
u1(1.45517368832359) q[6];
u3(-1.01513874518641,0.0,0.0) q[0];
cx q[6],q[0];
u3(-0.389633305757314,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.26814774760473,0.365356380428376,-0.513626067537162) q[6];
u3(1.21402261790641,-1.53063293537860,0.553746011379051) q[0];
u3(1.37087140131027,-1.09741413642557,1.71802745509888) q[2];
u3(0.167995956446570,-1.18399389995273,0.144793808502360) q[3];
cx q[3],q[2];
u1(-0.391377129251844) q[2];
u3(-1.92856376474930,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.49294240857188,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.50951603468969,-2.76983902385198,0.397484343266572) q[2];
u3(2.86291315615663,-3.61691453494229,-0.912947266024661) q[3];
u3(0.446451122614485,-3.09155246033480,2.78647392826151) q[7];
u3(0.616439987223009,-0.560111923857775,-2.24743353430782) q[5];
cx q[5],q[7];
u1(1.55389533844075) q[7];
u3(-1.00626424185313,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.59531045134323,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.62762292423038,-2.40477860337501,1.87121085438308) q[7];
u3(2.06044529734255,3.27055283075057,-0.163474658013370) q[5];
u3(2.08016497609821,3.68484155082454,-1.60371040015171) q[8];
u3(1.21615932329678,-0.542680849827508,2.71796465612761) q[4];
cx q[4],q[8];
u1(2.45504148994071) q[8];
u3(-1.61865022887464,0.0,0.0) q[4];
cx q[8],q[4];
u3(3.57681488054584,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.83494607798058,-0.930385707929597,-3.13437233793909) q[8];
u3(1.71431880536988,-1.11514590954099,3.94931710375360) q[4];
u3(0.844931413254166,-2.20316356846733,2.15251997583786) q[9];
u3(0.385502837238785,0.713904503434731,-3.63726375673836) q[5];
cx q[5],q[9];
u1(0.470370727477659) q[9];
u3(-1.54584269406895,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.20839114392362,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.84828278230889,-0.708197287819902,-1.85828646055319) q[9];
u3(1.57847496555487,-4.46326272552017,1.61988199264926) q[5];
u3(0.548473637334048,-2.36615386098318,3.42729987766591) q[0];
u3(0.947476974679615,-2.93311842279467,1.28515998220798) q[6];
cx q[6],q[0];
u1(-0.129612406211073) q[0];
u3(-2.02830952791189,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.63075106265071,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.60894070269868,2.15979907344896,-3.00824639920440) q[0];
u3(1.59652729164279,-1.33947542099481,-1.51492420266453) q[6];
u3(1.88570875456399,1.54613985470068,-2.64728249565015) q[2];
u3(1.94282771320746,2.07447192569321,-2.63285902071464) q[1];
cx q[1],q[2];
u1(0.808987079337166) q[2];
u3(-1.44581305548891,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.67898669204038,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.78350745254862,-3.58319098639197,1.79146479159617) q[2];
u3(2.55920689195137,3.35316156007359,-0.820855967822981) q[1];
u3(1.94548363054065,1.37377867315663,-2.91957058978482) q[4];
u3(2.32361208416136,-1.81583658793986,3.81011060592133) q[7];
cx q[7],q[4];
u1(-0.183015967004470) q[4];
u3(-1.42687269388247,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.618662636187793,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.55876121385399,2.24008687682026,0.341181441679173) q[4];
u3(1.94611349029042,0.217637920005187,4.58751575560237) q[7];
u3(1.90234964306779,-0.708670132084288,-0.795015882119176) q[8];
u3(1.97396885766344,-2.56687280027072,0.486679921011203) q[3];
cx q[3],q[8];
u1(1.70459048458456) q[8];
u3(-3.01561777560558,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.907606354731120,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.26743467398946,-2.89762997905563,-0.296382771148436) q[8];
u3(1.65175280116445,0.293617256323861,4.37258705127313) q[3];
u3(1.04035659765371,2.45047970921981,-0.973366582258053) q[2];
u3(1.30333224189298,1.19876142335865,-1.02711946368946) q[9];
cx q[9],q[2];
u1(1.38111265216837) q[2];
u3(-0.783662741140334,0.0,0.0) q[9];
cx q[2],q[9];
u3(-0.542750801276306,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.86859167370427,5.15058851709483,-0.746574382804607) q[2];
u3(2.69400210525148,-0.00460137044985198,3.41113705059067) q[9];
u3(2.06945427802031,-0.677723602702220,1.03578916953784) q[4];
u3(2.26336231137640,-1.58883468795405,-1.66988583981543) q[3];
cx q[3],q[4];
u1(1.68623301178162) q[4];
u3(-2.96239327085601,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.505034180439916,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.47691667064444,-2.70189795270997,2.47980129444382) q[4];
u3(2.56480904372067,-2.85575912478354,1.28832477663006) q[3];
u3(2.27663368293530,-1.89743756234232,3.82135090593148) q[7];
u3(0.686705347698293,2.16167940243113,0.556926193331735) q[8];
cx q[8],q[7];
u1(2.49268032221008) q[7];
u3(0.236259094608195,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.53061277557801,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.33344735998443,0.287409234723637,-2.68019172584403) q[7];
u3(0.254989483605652,5.77485688954416,-0.205409032205920) q[8];
u3(1.85766732561292,2.33179834017069,-1.68337747199646) q[5];
u3(2.43400838696712,1.07201381554594,-0.895501420135358) q[0];
cx q[0],q[5];
u1(1.79200917780231) q[5];
u3(-2.37341880308536,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.06641069169531,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.77534445401210,1.42335242891323,-2.36771347826101) q[5];
u3(1.47742945369897,1.69589295308469,-1.35819548067671) q[0];
u3(0.962422667314029,3.27333460663039,-2.95074723147044) q[1];
u3(0.663984641618679,1.36773715125563,-1.81813491592196) q[6];
cx q[6],q[1];
u1(2.23778949258120) q[1];
u3(-2.92840128646308,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.34560504321688,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.996873398132949,-0.163215556083470,1.11617794708522) q[1];
u3(2.47423521271589,-4.09268867172215,-0.878229993365665) q[6];
u3(1.84210266257617,-0.00290941826952040,1.06218327149622) q[6];
u3(2.16599500336101,-1.06911807779127,-1.79514981464932) q[3];
cx q[3],q[6];
u1(2.74210933981560) q[6];
u3(-1.40109845230117,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.0748691911537589,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.13829921261733,-2.67325625877937,3.60537113450781) q[6];
u3(0.839137001646865,0.973334852801160,4.15991548794901) q[3];
u3(2.50673225602590,1.40063013047745,-0.655191573210084) q[8];
u3(1.92692116159494,0.941788845629697,-4.00221458283070) q[9];
cx q[9],q[8];
u1(2.81945771406622) q[8];
u3(-2.28837428556205,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.212360374592182,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.14092313955640,0.468774887172784,-2.15060315031429) q[8];
u3(2.48651315322462,1.94262098884487,0.176684969682868) q[9];
u3(0.752430364565973,-1.33400857318725,-0.771972097523825) q[0];
u3(1.39451376544178,-3.71973502271516,-0.586796712316151) q[5];
cx q[5],q[0];
u1(1.86437276073225) q[0];
u3(-3.13564564484583,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.814944660385515,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.789020882935989,2.71727480507178,-1.15369632950662) q[0];
u3(0.698660223409116,-2.75001662938803,-0.567118777096608) q[5];
u3(1.95675549388956,0.185390678029150,2.49828767676499) q[2];
u3(1.88368582047897,-0.725186297134078,-0.954598048278723) q[1];
cx q[1],q[2];
u1(3.26437522554511) q[2];
u3(-1.62787046691246,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.21795482264746,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.96789548095607,1.12243911950080,-1.97169828032096) q[2];
u3(1.89911556515158,-3.30566581639949,-1.77213240083880) q[1];
u3(0.728637217136571,-3.90632034890319,1.10685339551685) q[7];
u3(2.30520855373847,0.337515629819015,4.67255865605961) q[4];
cx q[4],q[7];
u1(2.15768588160545) q[7];
u3(-0.100602148701210,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.14946471706815,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.960496240426710,-0.535170480623402,-0.0798432297332653) q[7];
u3(1.57682467750625,-0.266755540726372,-2.70110149406271) q[4];
u3(1.15347269659611,1.70769030133157,-2.85421708847514) q[9];
u3(1.29978318440288,2.39251933547857,-3.66886916422160) q[8];
cx q[8],q[9];
u1(2.86831609272771) q[9];
u3(-2.09798838821500,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.19504147298089,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.20218345244524,0.426249667907806,-2.40192817825702) q[9];
u3(1.31622958066337,-2.23516271100751,-4.03411452931389) q[8];
u3(2.93947291841114,0.0518641710660672,2.68754548072449) q[1];
u3(2.44907248177762,-3.65392818001352,-2.26141494130255) q[2];
cx q[2],q[1];
u1(1.26853726452593) q[1];
u3(-0.224415691627960,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.15473476745801,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.12368729635769,1.61558525373327,-3.86947068051927) q[1];
u3(2.53957144797578,-0.0925397007324045,1.15539244144411) q[2];
u3(0.734387044797357,0.223921203420820,1.16489939486156) q[6];
u3(1.52916605354138,0.148246963538696,-2.78348176776895) q[3];
cx q[3],q[6];
u1(1.62918858223807) q[6];
u3(-3.06039357679284,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.772471766200457,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.841771718320665,-0.0689087744410993,-1.13341944563496) q[6];
u3(2.13723769570757,0.214223085846893,5.19886653186681) q[3];
u3(0.447133361842789,-1.12927137965471,1.33655999379572) q[0];
u3(1.34044702510269,-1.22123661641330,-1.70926288448774) q[4];
cx q[4],q[0];
u1(1.61091291979603) q[0];
u3(0.499829075593848,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.566674797519689,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.38782381852731,0.240124258211322,2.48907856633281) q[0];
u3(1.46726871490178,1.06911145996464,1.46922999187373) q[4];
u3(2.09391454578585,-2.07212623744612,-0.903931786658325) q[5];
u3(2.02245545777573,-3.39711547001314,-0.800975041120755) q[7];
cx q[7],q[5];
u1(-0.244885434248036) q[5];
u3(-1.88210762835861,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.731499217406300,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.60497616485717,0.170251888884039,2.75027414593486) q[5];
u3(2.08622160140556,-2.49590641184111,3.13965839628805) q[7];
u3(2.39679263793940,-1.12540234858762,-1.93695381349438) q[0];
u3(0.798569631397591,-2.14449234782442,-2.97743448091037) q[1];
cx q[1],q[0];
u1(-1.33930285233912) q[0];
u3(0.236200073287032,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.28775510474175,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.47652553421915,1.44133054424407,-1.86206294945577) q[0];
u3(1.64928904077711,0.920010112159387,-4.18232737212025) q[1];
u3(0.638897188482287,-1.24556710737484,1.31848053569390) q[2];
u3(0.737305908627403,-0.790747833230589,-1.50708357493907) q[7];
cx q[7],q[2];
u1(3.17216693351679) q[2];
u3(-1.89739318004726,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.544497541866051,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.51837321989999,-1.28983202014336,4.09144066388756) q[2];
u3(1.99926589906163,-1.82540662175523,-3.24580232020823) q[7];
u3(1.92197475452879,1.67998410361033,0.286680970323381) q[3];
u3(1.51189010009016,0.500570263084771,-2.90974885388785) q[6];
cx q[6],q[3];
u1(1.34463648826166) q[3];
u3(-3.38477483153219,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.32091345431270,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.51092162298739,2.31797293929413,-2.45259091974607) q[3];
u3(2.01664235095639,0.965733959552938,4.31734888540114) q[6];
u3(0.972776960882146,-1.42739012381863,0.241035509515947) q[9];
u3(1.98992994614972,-4.36783264674267,0.561055744155208) q[8];
cx q[8],q[9];
u1(1.13042466461647) q[9];
u3(-0.616647836623594,0.0,0.0) q[8];
cx q[9],q[8];
u3(-0.0625735109098586,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.08379728827237,2.05265389735011,1.03632818620046) q[9];
u3(0.974956772661827,-2.22904518206316,-2.94841069569625) q[8];
u3(0.998423464553055,0.155428683183817,1.50071640409880) q[5];
u3(1.52499100947271,-0.505286854455179,-0.148962088124015) q[4];
cx q[4],q[5];
u1(1.32976808880236) q[5];
u3(-3.34086706446640,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.22764775052912,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.94109681012066,-0.323739348540349,-1.04785245204838) q[5];
u3(0.988021616629667,-0.578622521910583,5.31429170028148) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
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
