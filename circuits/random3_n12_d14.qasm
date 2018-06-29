OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.31623087811472,0.721143185573704,1.82046787313236) q[0];
u3(1.56045704977838,-1.57631923436488,-2.66825631242661) q[1];
cx q[1],q[0];
u1(0.438117902066069) q[0];
u3(-3.24703817028459,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.72721166323136,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.97434089111628,-0.500578256751928,1.03005433986559) q[0];
u3(0.144446879491535,4.45691597366079,-0.888803239630183) q[1];
u3(1.89279923132352,0.582211888965595,-2.51826244872652) q[2];
u3(1.82396152221182,-3.51326886162697,2.51292667024534) q[6];
cx q[6],q[2];
u1(3.21887882825848) q[2];
u3(-3.42808683952836,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.15822893185111,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.49339940018435,-1.54577074591727,-1.50760657057907) q[2];
u3(2.10830700596049,-3.98342424776562,-2.10773413637904) q[6];
u3(1.63792681032503,-1.19201224173238,-1.51286416804894) q[4];
u3(0.267032522988728,0.988548676567178,-4.21610931284402) q[9];
cx q[9],q[4];
u1(2.04295784695111) q[4];
u3(-2.43036753023251,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.885108938603835,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.89018575559632,-2.65336344728916,2.62840719617014) q[4];
u3(1.00755663964295,-1.11218925587458,4.61997664511252) q[9];
u3(1.09173580481285,3.54438988679594,-1.79623679238742) q[7];
u3(1.20735260568631,1.92974860791965,-2.22932449746223) q[11];
cx q[11],q[7];
u1(3.82090773130158) q[7];
u3(-4.15948574198624,0.0,0.0) q[11];
cx q[7],q[11];
u3(-1.20506649780750,0.0,0.0) q[11];
cx q[11],q[7];
u3(0.857443986967126,-0.952119263476349,0.847482558646401) q[7];
u3(0.824344880296457,-1.69235559407257,-3.18371225232151) q[11];
u3(2.40453184069727,1.69780762645595,-1.05642510374528) q[10];
u3(2.32502357737710,1.53682139799045,-4.43180647949517) q[3];
cx q[3],q[10];
u1(0.428445951434664) q[10];
u3(-1.37759415353095,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.03116542581002,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.562072437136987,0.706776765551343,-2.22320132751050) q[10];
u3(1.87233620655000,0.446234220465240,3.92814660499049) q[3];
u3(1.38658186981375,-1.85733819730295,-0.789043142091290) q[5];
u3(0.859653043224674,-2.93253036282829,-0.377300911087842) q[8];
cx q[8],q[5];
u1(1.74109439241518) q[5];
u3(-2.86145349502840,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.03192745917708,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.16304280805165,-1.38693993489219,0.573216763737008) q[5];
u3(1.64728430051286,1.42580101421495,4.81256366005296) q[8];
u3(0.878363423364701,0.954110634281364,-1.13596221998323) q[4];
u3(1.75462560328797,-4.64381376795978,1.15736706969024) q[5];
cx q[5],q[4];
u1(2.91307824712800) q[4];
u3(-1.63038888102654,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.02697351007172,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.05482081779102,-3.69574218642478,1.23485255884866) q[4];
u3(1.70309806188118,3.01110829650224,-2.73930501286771) q[5];
u3(2.48970548735786,1.63584549094373,-1.54744718676009) q[8];
u3(2.09820984424785,1.09750226808581,-3.86782481222861) q[7];
cx q[7],q[8];
u1(0.648238850735021) q[8];
u3(-1.13164710004061,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.19914426701180,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.33974319093857,-0.624860952137672,-1.04702325981374) q[8];
u3(2.50769221467721,-1.79817947586741,-3.60006742989721) q[7];
u3(1.34645349276937,0.381297250888309,-1.00130595455581) q[3];
u3(0.976799200485485,-0.471593470603050,-0.343532689104916) q[6];
cx q[6],q[3];
u1(1.69293404545045) q[3];
u3(-2.70457204826075,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.762067879774571,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.824152227496291,0.605713667998185,-0.910443313936522) q[3];
u3(0.732072768711484,0.743035307111045,4.11905844247241) q[6];
u3(2.54493842416383,-0.871546324742374,-2.17937722916939) q[11];
u3(1.81817023132948,0.811140047227657,-4.46917024575334) q[2];
cx q[2],q[11];
u1(0.144929742345377) q[11];
u3(-1.38626992248542,0.0,0.0) q[2];
cx q[11],q[2];
u3(2.26280095482170,0.0,0.0) q[2];
cx q[2],q[11];
u3(2.86645038664307,-0.387949627851711,-0.221552566117179) q[11];
u3(1.99597691173915,1.02102914719327,1.08546980182201) q[2];
u3(1.91004830474584,0.0172561876278806,0.709953014885079) q[9];
u3(0.315935931090106,-2.75385501602917,-1.11380734933576) q[1];
cx q[1],q[9];
u1(3.30345635198103) q[9];
u3(-0.891391508559092,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.81031137884128,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.21972591507447,0.468240279176515,0.265596408813046) q[9];
u3(0.432628689794709,-5.64644076639681,-0.413288864720858) q[1];
u3(0.469943142523327,-2.69626948995665,2.56940935395794) q[10];
u3(0.948700643831624,1.65331160542281,-1.92215973074340) q[0];
cx q[0],q[10];
u1(1.54105906574179) q[10];
u3(-0.834723469844153,0.0,0.0) q[0];
cx q[10],q[0];
u3(3.07712175517744,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.74981137657333,0.196149614952694,2.36099483551920) q[10];
u3(1.10352266369263,0.124299581347122,-5.09506131825655) q[0];
u3(2.47824886115384,0.839748199056409,-2.81275537708307) q[1];
u3(2.63417585192072,3.58228509895026,-1.76361262398696) q[9];
cx q[9],q[1];
u1(1.89460383273548) q[1];
u3(-3.13371798895242,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.07373844199991,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.61576575938593,1.46079661469289,-2.13275865281534) q[1];
u3(2.01988246198140,2.73628709886875,-1.59988446639361) q[9];
u3(1.48101266359514,-0.982235350603142,1.33889328597559) q[10];
u3(1.52600764954972,-1.01875670638537,-1.37196412070781) q[3];
cx q[3],q[10];
u1(1.29269874790400) q[10];
u3(0.0627756865134359,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.62812204848733,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.727238742617817,-1.96024990800240,3.05303025754198) q[10];
u3(2.50509324480407,1.51458192471511,3.31697964137157) q[3];
u3(1.65574267220463,1.51842166195375,-2.84120948475568) q[0];
u3(0.956570804567809,-2.04249035690061,2.65093917451190) q[11];
cx q[11],q[0];
u1(1.19797535533792) q[0];
u3(-0.192710348592734,0.0,0.0) q[11];
cx q[0],q[11];
u3(2.98674189123184,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.96366208314446,1.01954208822465,-0.197708409649295) q[0];
u3(0.951477762434353,4.98493990325656,0.970075192826185) q[11];
u3(2.33641713058671,-0.206646357497017,-0.315244444528514) q[8];
u3(0.842580727329404,0.720835294637015,-4.69398512023170) q[6];
cx q[6],q[8];
u1(-0.349873882135486) q[8];
u3(-2.35837817791149,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.25336777046210,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.79270014953377,3.51566634503681,-2.51033043962333) q[8];
u3(2.38583704791364,1.99084092626078,1.93248174745951) q[6];
u3(1.57345717210530,-0.761037923278353,3.63547617523515) q[7];
u3(0.781906088633342,1.81388972273148,2.23104680567155) q[4];
cx q[4],q[7];
u1(1.05308348935932) q[7];
u3(-1.34377252020811,0.0,0.0) q[4];
cx q[7],q[4];
u3(-0.258544513696427,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.01589849921906,0.231702053236436,-1.91607743705181) q[7];
u3(0.736300268744495,-1.91310764073229,2.74418728677845) q[4];
u3(1.68483707459677,-0.190103455685190,2.65771318288197) q[5];
u3(2.04381653864218,-2.09759143352770,-1.73542568742475) q[2];
cx q[2],q[5];
u1(1.98384064555241) q[5];
u3(-1.67191133945341,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.21675297134125,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.80573082618984,1.80681464815607,-3.79772985532862) q[5];
u3(2.38943077349268,3.17939833297928,-1.48666183865688) q[2];
u3(2.21594182124451,-0.767417500229570,0.173720953188368) q[0];
u3(1.21920709031945,-2.61511496477180,-0.763306557156101) q[3];
cx q[3],q[0];
u1(3.25155985628111) q[0];
u3(-1.58919805671830,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.77682800999839,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.88575493072985,1.03598542843419,-2.45374626006057) q[0];
u3(1.53005402898536,0.206261160131196,4.84176778583777) q[3];
u3(1.62200348149998,0.781565577671509,-3.53804178417102) q[5];
u3(2.04356431926963,2.68746694052604,-2.70835257551890) q[9];
cx q[9],q[5];
u1(1.98964659975143) q[5];
u3(-1.53832549397884,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.273133467976605,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.660377680356138,2.67228758356484,-2.50435213277168) q[5];
u3(0.823704584189852,3.48236951377615,2.18420985288257) q[9];
u3(1.87149425967208,3.07215725251387,-0.540841373755472) q[2];
u3(2.52868891536627,1.29010844093133,-1.47284473463677) q[4];
cx q[4],q[2];
u1(0.672725933560061) q[2];
u3(-1.27049363119136,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.0575721448069169,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.46916811453475,-0.187207383658860,-1.96643881903346) q[2];
u3(2.60931910166088,3.29636180248739,0.105575576069447) q[4];
u3(1.77010996031380,1.84257692738425,0.763505896353933) q[1];
u3(0.807991715324878,0.920408649860420,-3.36363439433597) q[10];
cx q[10],q[1];
u1(2.49780900633848) q[1];
u3(-1.77834747245247,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.19808043815109,0.0,0.0) q[10];
cx q[10],q[1];
u3(2.83420479890942,1.13268613991303,-1.05286267434563) q[1];
u3(0.993846036586617,-2.92182127865895,-2.42076110220733) q[10];
u3(2.88660461095939,3.78657287658677,-1.73112205060222) q[7];
u3(0.964577278076409,1.20375678562824,-0.218532254947972) q[6];
cx q[6],q[7];
u1(3.43570685879016) q[7];
u3(-1.10744215721504,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.85823379670755,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.867047263203495,-2.98005931529693,2.11404416684389) q[7];
u3(1.60206543203200,2.48092259380041,0.0807077506523615) q[6];
u3(0.336984667981806,0.519501457639368,-0.489670839346098) q[8];
u3(0.241779893069807,-0.373482958442497,-0.893716131245272) q[11];
cx q[11],q[8];
u1(0.471732629243222) q[8];
u3(-1.08133987198014,0.0,0.0) q[11];
cx q[8],q[11];
u3(2.88559741457774,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.31849987483811,3.96922458049442,-1.86826106830571) q[8];
u3(1.50720795117646,2.79898533565152,-0.152922893438411) q[11];
u3(0.557814618215629,0.257161566152057,-1.78464202688017) q[11];
u3(1.52753417369316,-3.88689981029888,1.03997514686464) q[6];
cx q[6],q[11];
u1(3.78672767713357) q[11];
u3(-1.45011258939131,0.0,0.0) q[6];
cx q[11],q[6];
u3(2.11552953842079,0.0,0.0) q[6];
cx q[6],q[11];
u3(2.28069650043244,-4.39218167670042,0.0352243035411757) q[11];
u3(1.54670208910778,2.69486220425691,2.87020850824137) q[6];
u3(1.69937224100634,-0.187405431162921,1.58966132779519) q[9];
u3(1.62057693256008,-2.00146726121109,-1.51543753185023) q[0];
cx q[0],q[9];
u1(2.82887743118693) q[9];
u3(-1.06696420275787,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.56611875326930,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.92591516417399,1.22827772349383,1.79582942989821) q[9];
u3(2.59738418409802,-4.17335224401570,0.147679131657048) q[0];
u3(2.06761447934008,-2.61197310036149,0.0816406283238547) q[5];
u3(1.90743925029466,-4.43315771410128,-1.31284077422621) q[7];
cx q[7],q[5];
u1(3.56131131225070) q[5];
u3(-1.12333679364149,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.63405638260287,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.50783473900038,1.95792622531944,-0.578668513308828) q[5];
u3(1.89483342950092,-4.22255689593610,1.91747946595187) q[7];
u3(2.22703734075552,-0.334349982228950,-0.773241814318966) q[10];
u3(0.906052668924361,-5.23496595513088,0.898729938177044) q[4];
cx q[4],q[10];
u1(1.78557863522058) q[10];
u3(-2.60518874294084,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.0909845523022668,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.77143349432284,-2.29485188073440,0.357855283086029) q[10];
u3(2.64355521400441,-2.33591650879761,1.66400475699001) q[4];
u3(1.82321400019689,-1.85140412378913,4.36136905719626) q[3];
u3(0.827120770793130,-1.50718157965541,2.67583471707494) q[2];
cx q[2],q[3];
u1(1.64796871716920) q[3];
u3(-0.102219082015343,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.06092681581461,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.25754031665529,-1.79676477037278,-2.27759008512984) q[3];
u3(1.11064275322116,-0.580330383447068,-4.81305384667342) q[2];
u3(0.377174505232807,1.97402135592764,-1.66185864973035) q[1];
u3(0.995250292312358,0.233376453113485,-1.54941411557679) q[8];
cx q[8],q[1];
u1(0.812111438281280) q[1];
u3(-1.15213654540832,0.0,0.0) q[8];
cx q[1],q[8];
u3(3.50197926849501,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.914713581706520,0.104924607113637,-0.174107413773268) q[1];
u3(0.567499947180664,3.94830491112701,0.469658344034465) q[8];
u3(2.84327849038983,-0.864012525080615,0.178689316872333) q[8];
u3(1.42955379526046,-3.55223625013552,-1.49753115313200) q[9];
cx q[9],q[8];
u1(-1.43108100362093) q[8];
u3(0.786409550673924,0.0,0.0) q[9];
cx q[8],q[9];
u3(3.84480314858818,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.48903779324543,-0.982271531591081,-0.113118759879348) q[8];
u3(2.05754293269840,-0.213289815966245,-3.85659111340080) q[9];
u3(2.14377230674409,-0.113195326142332,2.93148671691421) q[3];
u3(1.97906831288546,-1.62543314849992,-1.59098266092915) q[5];
cx q[5],q[3];
u1(0.871277842475854) q[3];
u3(-1.51289920509695,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.93804525591716,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.96665437229160,0.629618856625690,-4.57311168051985) q[3];
u3(1.48272034509565,3.78772212941393,-1.76853618663179) q[5];
u3(3.06661123465558,-0.492302074077690,1.54537142320340) q[2];
u3(2.29962342219292,-3.51213698392001,-2.39326442975509) q[10];
cx q[10],q[2];
u1(2.47646625199707) q[2];
u3(-2.99798521734037,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.842763700502514,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.32933629196536,0.0742869959503060,-1.92457408486541) q[2];
u3(0.979469172883205,3.29721263295599,0.613104420818197) q[10];
u3(1.64696960592515,-2.69989532580924,0.0418300423134250) q[0];
u3(2.23783741654263,-3.91583710357505,-1.18757998088846) q[6];
cx q[6],q[0];
u1(0.786248642608757) q[0];
u3(-0.273010103250567,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.48505817225118,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.47143918874130,-1.70772319665729,0.890244207768582) q[0];
u3(1.20738638216162,4.17355576254433,-0.888497733440414) q[6];
u3(2.19559209065914,-1.47865550341956,1.62215329463951) q[11];
u3(1.74955900948503,-1.62265267419331,-0.852989465145790) q[4];
cx q[4],q[11];
u1(1.47488084395374) q[11];
u3(-3.41827321661308,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.06953126584960,0.0,0.0) q[4];
cx q[4],q[11];
u3(0.377735583678389,0.197206162517802,-0.969763329985057) q[11];
u3(1.38321431599489,-1.23788034862011,-1.26879770147841) q[4];
u3(1.99256053290034,0.409493381036029,-1.85640203759040) q[7];
u3(2.31363146341072,-3.74474051003078,1.15514827161036) q[1];
cx q[1],q[7];
u1(0.855121293315130) q[7];
u3(-0.0183526157837499,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.38962627746509,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.09636987642608,-2.97349943443015,3.22728988557205) q[7];
u3(2.53038703776733,2.32820450529158,-3.23570748326255) q[1];
u3(0.319560883009227,-0.0633887293461539,-0.984137822800023) q[2];
u3(1.11124872303006,-1.22577646372268,0.129786649024904) q[11];
cx q[11],q[2];
u1(1.60254758446576) q[2];
u3(-0.774176601458892,0.0,0.0) q[11];
cx q[2],q[11];
u3(-0.442916509749888,0.0,0.0) q[11];
cx q[11],q[2];
u3(0.747980455864696,-0.950144520059083,1.34378021158508) q[2];
u3(1.13074944425997,-0.416100503812434,-1.43895472962843) q[11];
u3(1.85976675840750,-0.0290232451703103,2.18454941551606) q[9];
u3(1.78579320167644,-2.12530082628153,-0.861910999364739) q[8];
cx q[8],q[9];
u1(0.481519359019640) q[9];
u3(-1.29303514431897,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.38176020716652,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.16846007694364,0.108626718346289,1.49448078516241) q[9];
u3(0.785211295482545,-3.16125570211758,-2.38613676303170) q[8];
u3(1.77531812284448,0.981434274722619,0.292505257776701) q[4];
u3(0.990867772417825,-0.350089594134607,-1.76540310026753) q[6];
cx q[6],q[4];
u1(1.87583969772464) q[4];
u3(0.658827157284607,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.08602128713642,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.31064193238040,-0.342413002162476,-2.28873388435662) q[4];
u3(0.322037929331436,4.07029929299621,-0.179106130702256) q[6];
u3(1.23860839129620,-3.74458368588404,1.59689824990822) q[1];
u3(2.76814457601813,-4.78855784522517,0.262895665917668) q[0];
cx q[0],q[1];
u1(-0.168018833518625) q[1];
u3(-1.76897422313576,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.12370285128792,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.63915559999943,0.0279113715977536,1.14580465628481) q[1];
u3(1.62347812353075,4.72598495357257,0.424385704655072) q[0];
u3(0.249659222312325,0.233134136928821,0.182280757206122) q[10];
u3(0.729336594618301,-2.91595719698667,0.795761622535064) q[7];
cx q[7],q[10];
u1(1.63937263839315) q[10];
u3(-3.46884029977166,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.39602943231133,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.47864493838479,-0.559525804460770,4.32250532712387) q[10];
u3(0.944844032472805,5.23467260484757,-0.792238110390040) q[7];
u3(0.495747522637424,0.470838813701002,0.313168322879984) q[5];
u3(1.25748656338194,-2.62297353266114,1.46051208834005) q[3];
cx q[3],q[5];
u1(-0.396064648283046) q[5];
u3(-1.51156733354586,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.782304782285654,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.98992905419160,-1.35128828534647,3.36973518968165) q[5];
u3(2.11116267531035,0.0163273260620347,-3.53373133195700) q[3];
u3(0.995519848220365,0.535923981762607,-3.57841545741352) q[4];
u3(1.16667583315684,2.53200543157912,-2.78829603544291) q[6];
cx q[6],q[4];
u1(2.37446523031935) q[4];
u3(-3.13303212280343,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.54286084430565,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.73093574532333,1.95075967406053,-1.79189535512353) q[4];
u3(2.20846375730113,-0.599445292850807,-1.67853924809326) q[6];
u3(1.68666432234318,1.89833381965357,-0.936142425024222) q[3];
u3(1.23955257219548,0.189273971533520,-3.59579611683041) q[7];
cx q[7],q[3];
u1(3.28907857412007) q[3];
u3(-1.20555200004188,0.0,0.0) q[7];
cx q[3],q[7];
u3(2.50326744100312,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.674396876493371,2.72568492883382,-1.54347433524993) q[3];
u3(0.822433780898879,-1.65761247722289,-0.798789301649709) q[7];
u3(0.863124373808918,-0.668595331931145,1.82231400043384) q[8];
u3(1.09862402789375,2.57094622428418,-3.69667434344867) q[1];
cx q[1],q[8];
u1(0.215367354651902) q[8];
u3(-1.34724776068689,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.39359829292023,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.52954572892323,1.90130418390195,-4.05229386374845) q[8];
u3(1.33230811876700,3.15542771813254,0.616513665032937) q[1];
u3(0.442838400960376,1.06398838673855,-0.466934225496947) q[0];
u3(1.71870002636121,-0.0402995994097861,-2.73735293096430) q[2];
cx q[2],q[0];
u1(3.08170850734604) q[0];
u3(-0.744747775189457,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.88873023374224,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.05810299995660,0.781132305666223,0.991367534009723) q[0];
u3(0.875784938333013,-2.75270989146008,2.73265332871854) q[2];
u3(1.16455650499047,-0.0442826324895747,0.903773821220599) q[5];
u3(0.849406714738697,-2.27253887912435,-1.64274041741430) q[9];
cx q[9],q[5];
u1(1.84946713850623) q[5];
u3(-2.22503694753264,0.0,0.0) q[9];
cx q[5],q[9];
u3(-0.281851640514158,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.511768480586662,-3.56769892366534,1.69525144091609) q[5];
u3(0.402003009497381,1.19471654890362,0.541408781153977) q[9];
u3(1.45607520414548,-0.520726848843212,-1.60916249392855) q[10];
u3(1.32995037218646,-4.57626772527351,1.47660897785375) q[11];
cx q[11],q[10];
u1(0.575593165946971) q[10];
u3(-1.65887508206743,0.0,0.0) q[11];
cx q[10],q[11];
u3(-0.223975401206998,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.67704026553343,2.88952698684492,-1.88699633969041) q[10];
u3(0.656254539086133,2.73573553966071,3.44086574743259) q[11];
u3(1.01389161029057,2.44120011438295,-3.36452436704888) q[1];
u3(2.05846537663218,-2.22517924342388,3.41435731285641) q[5];
cx q[5],q[1];
u1(1.31232827384126) q[1];
u3(-3.19914633720003,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.49548768408547,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.09126622556111,0.183749976763041,-1.39676480986342) q[1];
u3(0.534541314236381,1.02218402390770,-1.99797553899210) q[5];
u3(1.88745765840173,0.160019984721036,-2.35022545887069) q[6];
u3(0.750814653933511,-3.83166152470777,0.937636249289419) q[2];
cx q[2],q[6];
u1(3.59808657166207) q[6];
u3(-1.37937562386849,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.44101509531533,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.73354589359957,3.07295373505598,0.0600886967195009) q[6];
u3(0.832213553795608,0.603855615999832,-2.58232095394437) q[2];
u3(1.79285556154807,2.75986379287557,-1.66665279880391) q[4];
u3(1.58050909360874,1.70053748176202,-2.29990227205083) q[0];
cx q[0],q[4];
u1(1.44053649954311) q[4];
u3(-0.314353012521384,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.25242234770913,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.495227476692952,2.94064104430378,0.453071375082458) q[4];
u3(1.21640369096097,2.86131015291633,-0.170521176780845) q[0];
u3(2.12228569943347,1.14539398224613,-0.238184656915175) q[9];
u3(2.25213192256150,0.0905659183680148,-2.09946838492217) q[3];
cx q[3],q[9];
u1(1.75376430520562) q[9];
u3(-1.94252295858589,0.0,0.0) q[3];
cx q[9],q[3];
u3(3.86403884620815,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.09862565679451,-3.14113565659189,2.07002623560744) q[9];
u3(1.23920433705810,-1.09060509134837,-2.81674879932586) q[3];
u3(1.06973629737283,-0.438778248285661,0.223684620722726) q[11];
u3(1.41799163774419,-2.56632505778009,-0.919084352579130) q[10];
cx q[10],q[11];
u1(1.52133964136885) q[11];
u3(-3.47376631272952,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.15217188494009,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.67153353225170,0.964286271521186,-2.64884007681868) q[11];
u3(2.04796454392176,0.471364015289232,5.34237188445233) q[10];
u3(3.01087953522074,-3.21664878988693,0.755060052561144) q[7];
u3(1.94573704638037,-2.89521030384885,-0.781388665973443) q[8];
cx q[8],q[7];
u1(3.04371318288056) q[7];
u3(-1.84370319784116,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.04420321464637,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.650959604146667,1.76973193663949,-1.32089541540227) q[7];
u3(1.13963948071596,-1.28970490404783,-3.77805791202017) q[8];
u3(2.87069038134183,2.03535977951383,-3.76285544999923) q[9];
u3(0.845040191751223,-1.25628245928825,1.78988523212839) q[1];
cx q[1],q[9];
u1(3.61956645051091) q[9];
u3(-4.22013827689809,0.0,0.0) q[1];
cx q[9],q[1];
u3(-0.910735052083696,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.01014640821102,-3.89465777810301,2.04775423275960) q[9];
u3(1.07236581667331,-0.195496429546247,-0.369535792632883) q[1];
u3(1.03459295224458,-2.05337671149937,1.46126080127226) q[8];
u3(0.377928300759241,-0.797168300131819,-0.739913390868292) q[7];
cx q[7],q[8];
u1(3.13673673475268) q[8];
u3(-1.77952316575265,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.53579615591985,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.94643866159563,3.23335964354148,-2.62409934361137) q[8];
u3(0.977271847659755,2.14179675604415,-1.63608930833889) q[7];
u3(1.95332147944977,-0.351649485688895,0.277153127963131) q[5];
u3(2.33577984767598,-3.24379642720181,-0.152630568910500) q[4];
cx q[4],q[5];
u1(-0.936017308866025) q[5];
u3(0.256386366908309,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.80669903006147,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.85760353545630,2.79471231699040,-1.47201414756015) q[5];
u3(1.71545816375037,1.23236275190364,1.24568925276954) q[4];
u3(1.83033958132453,-0.517391132647509,-2.06456786497274) q[3];
u3(2.07138317829173,1.34087219651424,-4.12191823181741) q[2];
cx q[2],q[3];
u1(2.65281343123286) q[3];
u3(-4.63883810121450,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.529314575925739,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.28108168192614,-2.91406396843006,0.519832503203814) q[3];
u3(0.974424702272946,2.51439992913270,-1.54395620917550) q[2];
u3(2.29078132817956,2.45354716077963,0.502988066707487) q[6];
u3(2.15959539109522,0.320588382733353,-4.30071966921309) q[0];
cx q[0],q[6];
u1(1.21283793247253) q[6];
u3(-0.347430515519203,0.0,0.0) q[0];
cx q[6],q[0];
u3(3.06056503381152,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.50180893977360,0.469788265122381,1.82063309890057) q[6];
u3(0.647425369942975,1.16884463466607,1.55460101740592) q[0];
u3(1.30192071372738,1.82474324423979,-0.317493334699384) q[11];
u3(0.566605282525483,-1.07267499134104,-3.25602358210525) q[10];
cx q[10],q[11];
u1(3.00722331339708) q[11];
u3(-2.01937230195131,0.0,0.0) q[10];
cx q[11],q[10];
u3(0.910450217052976,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.41906886918550,-0.523881288859966,0.0976831659264606) q[11];
u3(0.991876371629111,-2.39525546651767,-3.86534244419385) q[10];
u3(1.94685453227031,0.996488604666422,2.00329612362132) q[3];
u3(0.267730280776466,-2.83021067249643,-2.31049569596220) q[7];
cx q[7],q[3];
u1(-0.308287598297450) q[3];
u3(-1.41216893054345,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.86150637810821,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.49784243946798,-0.719838410589985,0.325312946516695) q[3];
u3(1.70344506013727,0.391670361110797,-3.40756669810064) q[7];
u3(1.25145432849331,1.63180217526152,-0.0718036759800064) q[10];
u3(1.34904077487721,0.892437275287920,-4.00630730311871) q[5];
cx q[5],q[10];
u1(1.15496982423801) q[10];
u3(-0.567909581360090,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.27922977620632,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.51462066787801,0.865381861595866,-2.34984794170087) q[10];
u3(2.43195769216422,0.0785370828108465,-3.98101928209225) q[5];
u3(2.01042431489772,-1.05213078921672,2.21324040217249) q[0];
u3(2.84932605169281,1.36315591920492,3.64579693594570) q[8];
cx q[8],q[0];
u1(1.60671283873301) q[0];
u3(-3.21123210830678,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.941875686729126,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.73728986229256,-1.96491406875265,3.53648474895358) q[0];
u3(1.80945276057909,1.61491945796340,2.09536688632713) q[8];
u3(1.91960200942978,1.24565240851147,-0.331845650774936) q[11];
u3(1.43894960134690,-4.75408829941487,1.22791059949393) q[4];
cx q[4],q[11];
u1(3.12538087633321) q[11];
u3(-1.62186220721132,0.0,0.0) q[4];
cx q[11],q[4];
u3(0.575061212082843,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.82470380169474,-2.48731893942957,2.89792169058652) q[11];
u3(3.06889543428877,2.15193972335404,3.71107911346468) q[4];
u3(0.853299150567614,0.893999749718321,1.94943012832118) q[6];
u3(1.84532870379611,-0.834324580675972,-0.560474514571923) q[2];
cx q[2],q[6];
u1(1.57565953706633) q[6];
u3(-0.169678769019915,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.10628224448839,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.67644320110986,0.310885983454893,3.00121516838471) q[6];
u3(1.47264187229004,-2.60658802684144,-2.82055079861669) q[2];
u3(2.35581202045413,1.37693192935633,-3.40439182696763) q[9];
u3(1.69955469974859,2.50895849481765,-2.82489342742673) q[1];
cx q[1],q[9];
u1(-0.206613341422855) q[9];
u3(-1.24930971412817,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.40666899539828,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.50010191507155,2.78797534340413,1.40788430443424) q[9];
u3(2.29800263932813,-1.85860212019003,1.31605482142564) q[1];
u3(2.71569812335495,-3.43830009839065,2.82475829854569) q[2];
u3(0.712149342107413,3.75402216539915,-2.39964501713298) q[5];
cx q[5],q[2];
u1(-1.00957455880537) q[2];
u3(0.0729544571465768,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.66026009535877,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.37668223213707,2.21606477154943,-2.86646809325821) q[2];
u3(1.99510412739520,3.20721088226608,0.271649829846598) q[5];
u3(1.00189866963437,0.0816652448702240,2.12338792767598) q[6];
u3(1.06923839044387,-0.712642151819022,-2.42030243015830) q[10];
cx q[10],q[6];
u1(1.15692495840561) q[6];
u3(-1.40583549229447,0.0,0.0) q[10];
cx q[6],q[10];
u3(3.93751799316183,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.22185586536103,-2.41426593486365,3.09379319241318) q[6];
u3(1.55593289046649,3.14342699986229,-0.427549312066376) q[10];
u3(1.82053200863931,3.65318654024246,-0.634634742592230) q[0];
u3(1.57523771472762,3.25952700194712,-0.466631823126699) q[9];
cx q[9],q[0];
u1(1.84062174286000) q[0];
u3(-2.85301449923695,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.892950941466419,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.31414985946228,-3.29843404005776,0.960866186691748) q[0];
u3(0.995243938112315,-1.02222173184595,1.57427719265099) q[9];
u3(2.08848779481528,0.535783326730087,0.259963168978301) q[3];
u3(0.579788268467520,-0.887040521099117,-3.63689099555843) q[8];
cx q[8],q[3];
u1(1.45061702790661) q[3];
u3(-0.636159434058924,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.04222868394974,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.67024882568508,1.23580173966436,-2.11984820367592) q[3];
u3(2.27330429411970,5.98639276375151,0.0823979725084456) q[8];
u3(1.92619588070240,1.28081834283189,-3.62741857013496) q[7];
u3(1.09633728965966,2.05533490738939,-2.67280032909636) q[11];
cx q[11],q[7];
u1(0.674924059584335) q[7];
u3(-0.156855161124152,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.71549578044641,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.36336954671012,-1.64649585900403,2.49589126277476) q[7];
u3(2.37226573265614,4.27849284012430,-0.579565615473771) q[11];
u3(0.457539103614898,-1.63195483127219,-0.634456299333777) q[1];
u3(1.97646874937258,-4.20912277590408,1.42453196682985) q[4];
cx q[4],q[1];
u1(1.87691679259740) q[1];
u3(-3.05361141869094,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.608743269600445,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.07342195113663,1.67194135080046,-1.72862326659306) q[1];
u3(2.31010394346166,0.938101350920864,-0.473184673802543) q[4];
u3(1.02895340699385,3.30806150718954,-1.79724490547546) q[7];
u3(0.841066589253915,1.20120162399520,-1.55470284490133) q[1];
cx q[1],q[7];
u1(0.649505402276855) q[7];
u3(-1.29950548379849,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.154769614477492,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.82267274743419,0.660257518585295,0.275327633870295) q[7];
u3(2.29091355674883,3.00652615161468,-0.145660410578595) q[1];
u3(1.64556251133503,1.16273923454404,-2.12904291982107) q[9];
u3(2.79485984687502,-3.60823494814369,2.52745762063229) q[3];
cx q[3],q[9];
u1(-1.29085912227416) q[9];
u3(0.269220576158558,0.0,0.0) q[3];
cx q[9],q[3];
u3(3.26071032926808,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.62780172046384,-3.79347587015463,-0.748321904078884) q[9];
u3(1.58464207551697,1.48465187278444,-4.39202292226327) q[3];
u3(2.51418472891631,3.10112688258424,-2.95391621916653) q[5];
u3(1.06978394796875,1.94472406486299,-0.960615676540749) q[6];
cx q[6],q[5];
u1(0.365141712465939) q[5];
u3(-3.29835492647710,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.64069666221297,0.0,0.0) q[6];
cx q[6],q[5];
u3(3.01402987554633,0.0604281271731358,0.315920880090341) q[5];
u3(1.83295709842924,5.05701775459466,-0.277013069778638) q[6];
u3(1.61959110815112,0.303276289441585,-2.33167861685728) q[4];
u3(2.20932624025300,2.11568097266541,-3.92538593612407) q[8];
cx q[8],q[4];
u1(1.93071530838212) q[4];
u3(-2.26282762925540,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.484639317700753,0.0,0.0) q[8];
cx q[8],q[4];
u3(0.379261074682776,-2.92405796086864,2.68878473478606) q[4];
u3(0.501040853869679,2.78197837304833,-2.83384932358214) q[8];
u3(1.46512439623849,-0.487959437365656,0.761485171926943) q[2];
u3(1.01113611241029,-1.47324180905893,-1.58709715429419) q[11];
cx q[11],q[2];
u1(3.28065457708249) q[2];
u3(-0.762698017108533,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.00103475966229,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.28888355416515,2.61081094028964,-3.61720259784028) q[2];
u3(1.09426180919078,1.15678574513379,2.71970736316807) q[11];
u3(0.884181992677464,1.95013021809825,-3.54468526278121) q[10];
u3(0.935509994672099,-2.92703832332954,3.29634419985106) q[0];
cx q[0],q[10];
u1(1.96917226579706) q[10];
u3(-3.24179550719877,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.870253010819254,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.63931387114215,-2.36231454295144,2.21529695019004) q[10];
u3(1.16821820656275,-0.925739851678548,-1.02383413035010) q[0];
u3(1.37362665846322,0.856131344508585,1.14455415160067) q[9];
u3(0.667262503417366,-0.666809379470970,-2.37429010708387) q[4];
cx q[4],q[9];
u1(-0.00394518080150452) q[9];
u3(-1.79132715619975,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.696178894925165,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.713127291598888,1.37248843416837,-1.76412649794871) q[9];
u3(2.24802628246549,4.18793074019160,1.84069096568827) q[4];
u3(2.07481016059164,1.70136107203799,-2.83250475624359) q[11];
u3(1.00071318865044,-2.21081324645379,2.83358240896502) q[1];
cx q[1],q[11];
u1(1.68312671865013) q[11];
u3(-2.99792986879225,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.44115053518103,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.31650394896051,1.38926139044677,-4.87391745443993) q[11];
u3(0.695345143414870,-0.0433185729849224,-6.13051811861102) q[1];
u3(1.47901056420592,1.43099117930825,-0.705555224979191) q[8];
u3(0.572656508590844,1.09275779597482,-4.43359557572100) q[0];
cx q[0],q[8];
u1(2.93224874748440) q[8];
u3(-1.79091819810970,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.07827688489535,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.19063996151933,-2.64190309252821,3.20727159398466) q[8];
u3(2.03923089002224,1.56095947392600,3.85539145358239) q[0];
u3(0.935167584217461,-1.88970545750452,-0.581067972458931) q[7];
u3(2.00043334166098,-3.61475577928480,-0.633370056057729) q[5];
cx q[5],q[7];
u1(3.55194482261554) q[7];
u3(-1.21202514431375,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.19808550190265,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.624185501366411,0.232624943762500,-1.88672803209704) q[7];
u3(2.21561412167829,1.17752304123697,0.545863644116540) q[5];
u3(1.77685481581814,1.04177608989675,1.08179965564252) q[3];
u3(0.630271775410684,-5.31786702011777,-0.0100063359922293) q[2];
cx q[2],q[3];
u1(1.23540593370527) q[3];
u3(-0.165059652472382,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.62834234290545,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.09719493331003,-0.233255917384367,-0.999890662375274) q[3];
u3(2.37945749496035,-1.47310235742269,4.03570056587464) q[2];
u3(1.18226766764036,-0.143686288057633,-0.268515148263799) q[10];
u3(1.66715451760867,-3.79876565806058,1.62235384448232) q[6];
cx q[6],q[10];
u1(0.356884918311091) q[10];
u3(-1.15683417244841,0.0,0.0) q[6];
cx q[10],q[6];
u3(2.34533395737538,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.25623204625914,-1.63295397802585,-1.30414481422647) q[10];
u3(0.881607519764494,-2.01434094423477,-2.51365135277285) q[6];
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
