OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(2.22179282783637,0.119925650452045,-3.09504868878674) q[1];
u3(2.51594671134987,-0.966738322223421,-4.69008389853248) q[9];
cx q[9],q[1];
u1(1.59730258991676) q[1];
u3(-0.0577635635315485,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.86799170640980,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.976933635160737,4.44965555066193,-0.991954063708528) q[1];
u3(1.12621415708945,1.04917009898109,1.57981153858308) q[9];
u3(2.55526844757421,-1.64342705997051,0.572524358293254) q[2];
u3(2.68370083353496,-2.37709679654910,-0.0995437742694989) q[11];
cx q[11],q[2];
u1(2.27924070049910) q[2];
u3(-3.00595596326180,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.544423668538616,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.60086943283385,1.19810006408162,-1.95875849303222) q[2];
u3(1.37309840884658,-3.52920265285003,-2.43807224766726) q[11];
u3(2.24894485910959,1.94517758902929,-3.62586597474811) q[6];
u3(1.71986173540152,-2.24928719402604,3.36927006888406) q[14];
cx q[14],q[6];
u1(-1.31925464271219) q[6];
u3(0.282184269772454,0.0,0.0) q[14];
cx q[6],q[14];
u3(3.33851323366521,0.0,0.0) q[14];
cx q[14],q[6];
u3(1.86113672883991,-1.78816978073739,-1.43721160475830) q[6];
u3(0.562053464792198,-3.23264232701645,-0.152872540216066) q[14];
u3(1.52150096560600,-2.30060768300264,-0.108556558760575) q[3];
u3(0.879983043240260,-3.33966908343079,0.465714752115796) q[4];
cx q[4],q[3];
u1(0.440473634675450) q[3];
u3(-1.18412648095743,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.07380384523773,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.66009051306557,-0.736118082455835,3.67715390687094) q[3];
u3(1.75414306379571,1.76984939374881,-1.79316762166700) q[4];
u3(1.97863929670796,0.0161367466300456,1.45328844882512) q[13];
u3(1.85777169416602,-2.45213681212186,-2.00149998814117) q[12];
cx q[12],q[13];
u1(0.376998586355594) q[13];
u3(-1.64100954529001,0.0,0.0) q[12];
cx q[13],q[12];
u3(2.10230492563406,0.0,0.0) q[12];
cx q[12],q[13];
u3(1.85983255358882,-1.18789732724586,3.06675453873159) q[13];
u3(0.816376975090152,-2.97871983400631,2.45016736316671) q[12];
u3(1.34869049004386,-0.660333748466211,-1.50369370850315) q[8];
u3(2.43665762416699,1.66830225388968,-3.71825604853259) q[5];
cx q[5],q[8];
u1(-0.973628572945376) q[8];
u3(0.265795003149819,0.0,0.0) q[5];
cx q[8],q[5];
u3(3.83911048048208,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.65601585971208,0.165423681379099,2.54290125198581) q[8];
u3(0.721814556986914,0.991167122104676,-3.73540435291726) q[5];
u3(2.08231353982013,4.49879416184213,-1.57782472661785) q[0];
u3(0.742791039680305,1.20131335155398,0.641761423574118) q[7];
cx q[7],q[0];
u1(0.660552252279611) q[0];
u3(-1.36385297480642,0.0,0.0) q[7];
cx q[0],q[7];
u3(3.03195227970804,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.783369641940820,-0.535330001408262,3.86868360378245) q[0];
u3(1.43283141433652,-0.684679769007412,-3.91036239513457) q[7];
u3(0.876273846159507,0.985021929139898,-1.45401852533501) q[10];
u3(1.12265200492679,-0.628276036327167,-0.0264399255204247) q[15];
cx q[15],q[10];
u1(1.57498879530081) q[10];
u3(-2.20839968884175,0.0,0.0) q[15];
cx q[10],q[15];
u3(3.59720290764286,0.0,0.0) q[15];
cx q[15],q[10];
u3(1.55614222497559,-1.45719066105570,-2.05968450281788) q[10];
u3(0.933698462898735,-2.17631249708328,-3.45399245516465) q[15];
u3(2.70181526769946,1.52037510788045,0.0853918728523044) q[15];
u3(1.58186423665056,-0.605549651051101,-2.92517984702184) q[5];
cx q[5],q[15];
u1(-0.00100129422123296) q[15];
u3(-1.05120619409181,0.0,0.0) q[5];
cx q[15],q[5];
u3(1.78147883135452,0.0,0.0) q[5];
cx q[5],q[15];
u3(0.115253132120303,-3.47911058904989,0.761010061148352) q[15];
u3(0.677014853704618,-0.731716902551871,5.48029854194148) q[5];
u3(2.69241408892446,0.0790022569761055,-3.08920373066214) q[3];
u3(2.29703864206863,5.19944798825203,1.06254927044251) q[13];
cx q[13],q[3];
u1(2.05417257287387) q[3];
u3(-2.57221033416666,0.0,0.0) q[13];
cx q[3],q[13];
u3(1.31094089384354,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.39776699155167,-2.59681069165374,0.980176533130173) q[3];
u3(1.89424890911151,-1.76895695803829,-1.78320045789037) q[13];
u3(0.379943276465442,-1.46526727865622,2.13590423220769) q[14];
u3(0.342120346181319,-2.84260011674160,0.498256979431734) q[7];
cx q[7],q[14];
u1(2.31194140121299) q[14];
u3(-3.03038112679467,0.0,0.0) q[7];
cx q[14],q[7];
u3(1.41074272370284,0.0,0.0) q[7];
cx q[7],q[14];
u3(2.13174384992349,0.195428463967767,-2.38117087228563) q[14];
u3(2.36646060331308,-0.923989317691656,-4.89768101539981) q[7];
u3(1.82737232512603,-0.565690315444470,-0.889554664121535) q[2];
u3(1.32224388291018,-3.29686322697309,-0.00510226940596792) q[8];
cx q[8],q[2];
u1(-0.388668702249260) q[2];
u3(-1.91486527772214,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.50283596875814,0.0,0.0) q[8];
cx q[8],q[2];
u3(0.357859963473671,-2.97560078919365,2.57897720006785) q[2];
u3(1.06424737490346,-0.622480183474074,-4.16187219948821) q[8];
u3(1.41813694807698,-2.63112990583237,0.0558747533314405) q[0];
u3(1.81061285331857,-3.60202755788098,-0.899164786985247) q[6];
cx q[6],q[0];
u1(2.68881068014425) q[0];
u3(-1.63268882562262,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.771485966018625,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.562919835961741,0.361147970189378,0.420127563238091) q[0];
u3(1.70531423532485,-0.372836323613244,1.92495916088063) q[6];
u3(1.11113741582514,2.41478127726344,-2.18531396769960) q[9];
u3(0.528982586340207,1.62726527955371,-2.07975200455779) q[12];
cx q[12],q[9];
u1(0.747909272969596) q[9];
u3(-0.333682951963137,0.0,0.0) q[12];
cx q[9],q[12];
u3(1.46595782396171,0.0,0.0) q[12];
cx q[12],q[9];
u3(0.496809007483388,3.42383837097293,-2.21004286299515) q[9];
u3(1.19165431214847,-0.0522170757270648,5.00022518059407) q[12];
u3(1.21101128782966,0.271147232172939,2.78090757013878) q[4];
u3(1.77143310224781,-2.66718347082859,-2.21367883148699) q[10];
cx q[10],q[4];
u1(0.984603935031132) q[4];
u3(-1.35465902880531,0.0,0.0) q[10];
cx q[4],q[10];
u3(-0.0691135671041712,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.32304909878093,-1.62402806809334,3.54073919543709) q[4];
u3(0.813385893599823,-0.253794673812369,3.78598750776148) q[10];
u3(1.95039994476609,-0.662181720437659,1.26525031662828) q[1];
u3(2.12812604998236,-1.91150377295506,-0.586666522477663) q[11];
cx q[11],q[1];
u1(3.39420293983382) q[1];
u3(-3.70037109575348,0.0,0.0) q[11];
cx q[1],q[11];
u3(-1.17702688397999,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.42346915694380,0.283230667203306,-3.43723780094360) q[1];
u3(1.83069791343174,0.994654653808208,-4.68684327787666) q[11];
u3(0.212284413990357,-2.85487798256004,2.76793106988138) q[14];
u3(0.375773204385269,1.27499938226978,-2.29742433995973) q[11];
cx q[11],q[14];
u1(1.47948400370987) q[14];
u3(-3.04615120911717,0.0,0.0) q[11];
cx q[14],q[11];
u3(2.78076786670135,0.0,0.0) q[11];
cx q[11],q[14];
u3(2.37170460015776,1.56517203951256,-3.59342097908237) q[14];
u3(1.70387127342175,-0.176911555177913,5.00813830952379) q[11];
u3(0.963560449063253,-2.46201470356284,2.71235852206590) q[6];
u3(0.559672146184915,-3.53311957873865,2.15032091655655) q[9];
cx q[9],q[6];
u1(3.24408377477957) q[6];
u3(-2.02317008802914,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.859239514885510,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.04447014367966,4.13582411653113,0.122118424666052) q[6];
u3(2.31608450192946,1.23772621625821,-0.822825769176381) q[9];
u3(2.81565230732861,-0.345989097447534,-0.0191424440402669) q[1];
u3(0.886426061846273,-3.69546326782375,-1.44461830530168) q[2];
cx q[2],q[1];
u1(2.24535838222917) q[1];
u3(-2.63002357461249,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.06271788662182,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.91932733742590,-1.76859040014838,3.23050828164556) q[1];
u3(0.319655463679333,-0.535182072566873,0.846334031885145) q[2];
u3(1.25839750752704,2.62231175113071,-3.16872745046126) q[15];
u3(1.92787043949738,3.33219357294701,-2.39083692605971) q[13];
cx q[13],q[15];
u1(0.258559319512409) q[15];
u3(-1.28501997240864,0.0,0.0) q[13];
cx q[15],q[13];
u3(2.76279530704106,0.0,0.0) q[13];
cx q[13],q[15];
u3(0.928211621304555,0.952590003452534,0.455196898095151) q[15];
u3(0.786380079331384,2.36235521209323,-2.76805544713504) q[13];
u3(2.89937445494852,-2.20216746986661,3.43177710910652) q[7];
u3(0.991953114270134,1.97786577063369,-0.175863472572526) q[12];
cx q[12],q[7];
u1(3.53807535008470) q[7];
u3(-3.85766949174472,0.0,0.0) q[12];
cx q[7],q[12];
u3(-1.10499963143379,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.22901507006629,1.59708203539516,-1.70022869580176) q[7];
u3(2.52811425794930,1.74415224430505,4.08232063456595) q[12];
u3(1.46789176310363,1.93881447784336,-3.44493965795687) q[0];
u3(0.922742074169939,2.06391516963181,-1.98456673546721) q[10];
cx q[10],q[0];
u1(3.36745525791628) q[0];
u3(-1.18541110143795,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.59675617622568,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.37653112664449,0.658623329338186,-2.92978698965894) q[0];
u3(1.32313975965416,3.22254934223224,0.637424250381763) q[10];
u3(2.17210617362684,1.11363122378060,1.18449188782210) q[4];
u3(0.942811820154528,-4.99475261976272,-0.655015499143286) q[5];
cx q[5],q[4];
u1(4.27427742182633) q[4];
u3(-3.17722323512011,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.333056146337973,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.90533605670316,3.26489699695899,-2.09857703078803) q[4];
u3(1.97536829189591,2.91034562361372,-1.85156601291065) q[5];
u3(2.54231526279822,3.41036757767122,-2.41922585037022) q[3];
u3(1.10366320087275,2.70381347289451,-2.16782707246234) q[8];
cx q[8],q[3];
u1(-0.0370846000653706) q[3];
u3(-1.31922478070669,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.36493917479206,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.20047008131333,-0.858781533025612,1.52699915742731) q[3];
u3(1.68346066984771,1.08303224507983,0.694610886901807) q[8];
u3(0.318287141733368,2.16317406076246,-0.704075081463503) q[13];
u3(1.52405548803727,1.80598776264607,-0.959677239836030) q[14];
cx q[14],q[13];
u1(1.85830267508033) q[13];
u3(-2.72929448495081,0.0,0.0) q[14];
cx q[13],q[14];
u3(1.59211436893181,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.48545328617925,-2.03157457544222,-1.78821645903059) q[13];
u3(1.92116078504317,0.0449825991658002,3.28021941708180) q[14];
u3(2.40863338847824,-0.0788681912419559,0.0153789906873164) q[1];
u3(1.13497980544330,-0.211518245183105,-5.67382048355539) q[9];
cx q[9],q[1];
u1(0.310577885250591) q[1];
u3(-0.925001573397611,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.58170020842117,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.78661880788480,0.761609217750732,-0.0590109420479371) q[1];
u3(1.13040263684247,2.37061017432877,-1.18154841807494) q[9];
u3(1.27538741475092,-1.20856006017638,1.72440763933775) q[3];
u3(0.488191566492109,-0.990883976633113,-1.76825993556816) q[2];
cx q[2],q[3];
u1(0.221189654762387) q[3];
u3(-1.32772403443940,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.48546143571664,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.41511873599595,-1.63087007410256,3.36357541411203) q[3];
u3(2.25723728612017,-1.25091076305947,-1.44497724993215) q[2];
u3(1.96207446393147,0.700018962064119,-0.996700873335359) q[4];
u3(2.49547267407248,-3.48959169367531,2.21293247055369) q[0];
cx q[0],q[4];
u1(-0.615356014766390) q[4];
u3(-1.74786207094193,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.920084606844989,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.45840314475052,-0.787834527107320,-0.736104598930003) q[4];
u3(1.73730788319949,-0.977586341374246,1.72911594674951) q[0];
u3(2.21501644017971,-1.33819420589477,0.747562115951137) q[12];
u3(2.21676754502381,-3.44816271366474,-0.401321676830313) q[10];
cx q[10],q[12];
u1(1.16468909588911) q[12];
u3(-3.55605091272687,0.0,0.0) q[10];
cx q[12],q[10];
u3(2.22667985057530,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.28798668960970,-3.45209655304115,1.28520936904543) q[12];
u3(1.73941241820991,-0.157794121547687,3.37346762784876) q[10];
u3(2.32661103265230,-0.900296097552325,-0.692987508316500) q[5];
u3(0.966233073552614,-3.46094449973909,-0.687707173482257) q[7];
cx q[7],q[5];
u1(2.16224147985887) q[5];
u3(-1.85770406215222,0.0,0.0) q[7];
cx q[5],q[7];
u3(3.40664077076574,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.45402740329820,0.412019403270722,1.23090252736052) q[5];
u3(0.605653033346171,1.52260534504428,0.530187605406261) q[7];
u3(0.510577494787937,-0.188923255931231,0.0836864024347266) q[8];
u3(0.769080247776710,-1.15443908904192,-0.857159289592011) q[15];
cx q[15],q[8];
u1(0.000508359148702953) q[8];
u3(-1.59603097066754,0.0,0.0) q[15];
cx q[8],q[15];
u3(2.10679143151675,0.0,0.0) q[15];
cx q[15],q[8];
u3(1.72718031832996,-2.19067722798283,1.14252042455439) q[8];
u3(1.17084928579733,-0.437399400321002,-5.40335916732811) q[15];
u3(0.821633670356253,1.89308715503175,-3.20970654962650) q[11];
u3(1.44927791912206,-3.29920895737703,2.62120296153651) q[6];
cx q[6],q[11];
u1(-0.202233119373884) q[11];
u3(-0.558372605763879,0.0,0.0) q[6];
cx q[11],q[6];
u3(2.09952007170182,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.05683784218669,-4.13310276141707,2.10271050139350) q[11];
u3(2.43099442123458,1.12410019589667,-0.0304103599018812) q[6];
u3(1.47588427758435,0.560437047277836,0.253998880578512) q[3];
u3(1.31737060042774,-0.891129299120484,-1.66548500594186) q[9];
cx q[9],q[3];
u1(0.801057313785369) q[3];
u3(-1.39227084811831,0.0,0.0) q[9];
cx q[3],q[9];
u3(3.15420400691258,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.67410247802071,3.14188082992026,-2.26912785780712) q[3];
u3(1.39003691584289,-0.813618587056564,-1.75124076425044) q[9];
u3(1.91725344036488,-0.809039021732158,1.70326248192815) q[14];
u3(1.31794822232638,-2.07382013855674,-2.85030068061873) q[13];
cx q[13],q[14];
u1(0.428178336789542) q[14];
u3(-1.45954926985870,0.0,0.0) q[13];
cx q[14],q[13];
u3(2.41290663807282,0.0,0.0) q[13];
cx q[13],q[14];
u3(0.918333930928037,0.117159872306181,-3.40915271823283) q[14];
u3(2.27929593226460,-1.89855734696775,-3.19074900395064) q[13];
u3(1.70507066275870,2.08405708048530,0.394239251943245) q[11];
u3(2.61343280897085,0.798225150483942,-3.42992593359331) q[12];
cx q[12],q[11];
u1(2.40971370910449) q[11];
u3(-1.93429260956530,0.0,0.0) q[12];
cx q[11],q[12];
u3(3.23203180096823,0.0,0.0) q[12];
cx q[12],q[11];
u3(1.43757684046014,-0.261286811430557,3.60306340222244) q[11];
u3(2.22442578221123,3.00888545013401,-2.61398533852647) q[12];
u3(2.53536108051178,3.15868921968080,-0.357932510306330) q[10];
u3(2.37072213255875,1.73675817152199,-3.58316050369725) q[1];
cx q[1],q[10];
u1(-1.39679376975316) q[10];
u3(0.549832240059452,0.0,0.0) q[1];
cx q[10],q[1];
u3(3.86122100972660,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.36805699694948,0.295959690858639,3.63697731449370) q[10];
u3(1.03457570242215,-0.355714168755148,-3.09222718822187) q[1];
u3(2.48371447350669,2.50785870972158,-1.74372648913568) q[6];
u3(1.86931935530976,0.475280577952331,-2.10337616941135) q[15];
cx q[15],q[6];
u1(2.05012071149846) q[6];
u3(-2.44171787273413,0.0,0.0) q[15];
cx q[6],q[15];
u3(0.938359017458903,0.0,0.0) q[15];
cx q[15],q[6];
u3(0.627487528793132,-1.55878623590753,-0.653779942343273) q[6];
u3(1.01096061789069,0.951683391897521,-0.390366604395879) q[15];
u3(1.54114418662805,1.45662045488435,0.762489545338576) q[2];
u3(1.19556826711254,-0.0809920720309856,-2.78042993573962) q[0];
cx q[0],q[2];
u1(3.27870404415301) q[2];
u3(-0.732755611988235,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.97233896759744,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.47861704868790,-3.13880466077439,1.42790550539905) q[2];
u3(3.00687892753364,-1.35700147359059,-1.43556450547746) q[0];
u3(0.717441488736736,-2.54062952310574,2.42726715696062) q[7];
u3(0.375010117770778,-3.42468996133192,0.654579342474893) q[5];
cx q[5],q[7];
u1(0.964331079848349) q[7];
u3(-0.239548386945772,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.94765338255801,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.920280012276496,-2.01368612376968,2.44207583170148) q[7];
u3(1.76837022966593,0.407605058245169,1.35510206667776) q[5];
u3(1.34319645875098,0.252511595827121,1.09455542038128) q[4];
u3(0.982019452213460,-1.48349470372691,-2.57858871767492) q[8];
cx q[8],q[4];
u1(1.79609178105965) q[4];
u3(-0.631750609724241,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.97789165310154,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.67444710909907,-0.309709660710307,-1.86903647266318) q[4];
u3(1.28149538220436,-0.342081733159653,-2.52993018548162) q[8];
u3(1.81139207639019,-0.0679247002176127,1.20321998836548) q[0];
u3(1.79314568638640,-0.949798181856448,-1.93503267346799) q[12];
cx q[12],q[0];
u1(0.0100222418325890) q[0];
u3(-1.34191291277975,0.0,0.0) q[12];
cx q[0],q[12];
u3(2.20876738475357,0.0,0.0) q[12];
cx q[12],q[0];
u3(0.703861443342817,-0.735922501327391,2.80936780380730) q[0];
u3(0.886398319665794,-1.86664866973658,-2.48735220856402) q[12];
u3(1.91270353204450,1.40907422255430,-1.67198202851475) q[7];
u3(2.27433722752778,-4.99072429220868,1.17883230263936) q[11];
cx q[11],q[7];
u1(2.93366148025389) q[7];
u3(-2.16109800886931,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.11059194983653,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.64632870827106,1.51258471408834,0.846380074963251) q[7];
u3(1.94642073173229,4.48218822463876,0.445802461289817) q[11];
u3(0.745850752471278,2.33381801882594,-1.87681069815809) q[14];
u3(0.620672499206329,0.159526317307795,-1.49227280576433) q[13];
cx q[13],q[14];
u1(1.47744354759640) q[14];
u3(-0.831566136311193,0.0,0.0) q[13];
cx q[14],q[13];
u3(2.97408621395308,0.0,0.0) q[13];
cx q[13],q[14];
u3(2.10198406157836,2.05232425661975,0.521856539416639) q[14];
u3(0.306901719487976,1.48542238965274,3.77841659807481) q[13];
u3(1.41824436019053,3.43957917869936,-0.907119045563707) q[1];
u3(1.67102881169076,2.63504267798060,-0.995353362937238) q[15];
cx q[15],q[1];
u1(3.28179523854413) q[1];
u3(-3.69157457796588,0.0,0.0) q[15];
cx q[1],q[15];
u3(-1.08937200590512,0.0,0.0) q[15];
cx q[15],q[1];
u3(2.09451367786388,2.65789202066927,-2.11665330788747) q[1];
u3(3.01082600769520,3.43976736573999,2.81157752257100) q[15];
u3(2.38063955567486,1.19110867400083,1.84254077211631) q[3];
u3(1.80622539077065,-2.56253525872381,-2.48681258630464) q[5];
cx q[5],q[3];
u1(2.55314486307648) q[3];
u3(-1.57489972775419,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.282762395790047,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.56213457892728,-0.0758586814071844,0.417046271342062) q[3];
u3(2.64866128019286,-0.355865329707284,-0.223477645770333) q[5];
u3(0.418976547926038,1.40378111833043,-1.61446302871757) q[10];
u3(0.279617271472477,0.673294086830021,-1.24638577133934) q[4];
cx q[4],q[10];
u1(1.68659857135127) q[10];
u3(-2.59762976486311,0.0,0.0) q[4];
cx q[10],q[4];
u3(3.55233328020372,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.32735064482069,1.93291409102714,-3.79382360944103) q[10];
u3(2.29187895075807,-2.74723749845862,2.43974953179597) q[4];
u3(1.47584343418581,0.395966051395278,1.07334497662029) q[2];
u3(1.38166778871812,0.480731257622032,-2.74103498623518) q[6];
cx q[6],q[2];
u1(-0.241259069283588) q[2];
u3(-1.72411074844639,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.808662178901218,0.0,0.0) q[6];
cx q[6],q[2];
u3(3.01200675314718,-2.81017612649676,2.83524540658747) q[2];
u3(1.34210288808633,0.836051120345800,0.342656877630100) q[6];
u3(1.63084600140389,1.22807459840642,-3.22757844878864) q[8];
u3(0.526071414782193,-2.16401793112173,2.81370283100120) q[9];
cx q[9],q[8];
u1(1.85030741121687) q[8];
u3(0.258392031767631,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.784547336588944,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.89440191219914,-1.08591108286844,0.503868775168598) q[8];
u3(1.88329291046881,0.511706864628496,3.49657854473166) q[9];
u3(1.51866270990801,1.38964523135449,1.70116861273902) q[3];
u3(1.18895955769538,-1.39573878407115,-2.67977808890146) q[13];
cx q[13],q[3];
u1(2.44223218809554) q[3];
u3(-1.59152992046677,0.0,0.0) q[13];
cx q[3],q[13];
u3(3.50703951052542,0.0,0.0) q[13];
cx q[13],q[3];
u3(1.82560332615974,2.98684751865756,-1.57121248454233) q[3];
u3(2.21727877283410,-3.16147018699924,-0.615903578641181) q[13];
u3(2.41293686879170,2.16174812195093,-3.48281148215242) q[5];
u3(1.13556019212997,-2.54273809211665,2.79415395632382) q[1];
cx q[1],q[5];
u1(0.109148509158378) q[5];
u3(-1.12458960929982,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.59162219789620,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.02705136209282,4.66332260011233,-0.902825811077249) q[5];
u3(1.43497350445794,0.859512326233377,4.96962545876821) q[1];
u3(2.44909166615557,1.43499966791958,-3.80421491445590) q[4];
u3(1.77250687601800,2.66358816785738,-2.86040324844396) q[14];
cx q[14],q[4];
u1(4.49878713867415) q[4];
u3(-2.77811147875317,0.0,0.0) q[14];
cx q[4],q[14];
u3(0.566570191523831,0.0,0.0) q[14];
cx q[14],q[4];
u3(1.02547404714784,0.736916958143293,0.942580957520937) q[4];
u3(1.08981390523371,-5.82018551848930,0.127026155318889) q[14];
u3(1.65222993846180,2.36714187686276,-1.86315504019041) q[9];
u3(1.11716612375319,-2.65327584387176,2.64856460365495) q[7];
cx q[7],q[9];
u1(1.82516091424024) q[9];
u3(-0.0312976639938194,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.935490341557811,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.04980238094657,4.03729797434404,-1.82183948246079) q[9];
u3(1.94002181843123,3.87273186208464,2.02691100881338) q[7];
u3(2.47272288171593,0.285174160595449,2.18765150786484) q[12];
u3(2.06769139203501,-1.74654355449649,-0.553075323792595) q[6];
cx q[6],q[12];
u1(1.44014823357032) q[12];
u3(-0.409276819429217,0.0,0.0) q[6];
cx q[12],q[6];
u3(1.93110978239892,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.00875817553943,2.12312715340651,-3.47264595792267) q[12];
u3(1.77073008006930,0.733635037827774,1.95663200035994) q[6];
u3(0.860572149230284,2.73813607499790,-1.84085806339640) q[0];
u3(1.13599806633251,1.61750250542887,-2.32723153764319) q[10];
cx q[10],q[0];
u1(0.280611660899258) q[0];
u3(-1.50855924189962,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.27150981602819,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.06005854788343,-3.45992744600683,1.36058644442729) q[0];
u3(2.09704466096853,2.75727053339510,1.01317703125500) q[10];
u3(0.644112592369830,-2.49351169309395,3.33177075817480) q[15];
u3(0.556311338039620,1.83753562101389,-3.70228873203331) q[8];
cx q[8],q[15];
u1(1.46681561627639) q[15];
u3(-2.42933388257993,0.0,0.0) q[8];
cx q[15],q[8];
u3(3.15114636117457,0.0,0.0) q[8];
cx q[8],q[15];
u3(1.71387948198937,-1.18070588154411,2.15780409029351) q[15];
u3(1.05167190539353,3.03929861161670,-1.87870149520666) q[8];
u3(1.65211770136267,1.31332829525833,-1.04246343572969) q[11];
u3(1.26199444891649,-0.606998714584368,-3.12384027294221) q[2];
cx q[2],q[11];
u1(1.10445963246289) q[11];
u3(-0.535843388198851,0.0,0.0) q[2];
cx q[11],q[2];
u3(2.36915654729110,0.0,0.0) q[2];
cx q[2],q[11];
u3(2.15619616333949,0.919202441534205,-3.01861674263409) q[11];
u3(2.49972681346713,-0.838958503117101,-4.28004867435728) q[2];
u3(1.52108333593866,-0.579322698496056,-0.784193696699185) q[6];
u3(0.698308045334928,-4.27198098773755,0.581747072964178) q[14];
cx q[14],q[6];
u1(0.684235447011793) q[6];
u3(-1.43850425482282,0.0,0.0) q[14];
cx q[6],q[14];
u3(2.19991641055756,0.0,0.0) q[14];
cx q[14],q[6];
u3(1.71536226675491,0.349488085128089,-2.61720091831210) q[6];
u3(0.648219882293559,-2.19112669060576,-3.30507790675122) q[14];
u3(1.72145018119177,-0.743156578461508,2.94013869883654) q[8];
u3(1.21633201551977,-1.42283154693610,-1.65471648754784) q[3];
cx q[3],q[8];
u1(2.38037577400669) q[8];
u3(-1.78725521826071,0.0,0.0) q[3];
cx q[8],q[3];
u3(3.14374338245936,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.20082834700043,-1.11669733885010,2.62217095325890) q[8];
u3(2.74433699451547,3.72149879104398,-1.84265374613592) q[3];
u3(2.44707897469709,3.63191235239607,-0.938612444333539) q[5];
u3(3.01907276952373,5.04091915069627,0.171653264229565) q[2];
cx q[2],q[5];
u1(4.30493566860259) q[5];
u3(-3.08084333214842,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.322211238086467,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.531375473266521,1.12715053587099,-0.843521462398624) q[5];
u3(1.29840700005698,1.09507417591790,1.73049818660265) q[2];
u3(0.950255883406595,3.36482737283782,-2.73603873093141) q[13];
u3(0.911989771591326,-3.45278120074944,2.17028192357146) q[15];
cx q[15],q[13];
u1(1.78989491693764) q[13];
u3(-2.88179724902222,0.0,0.0) q[15];
cx q[13],q[15];
u3(0.904544285606851,0.0,0.0) q[15];
cx q[15],q[13];
u3(1.97004411614541,4.16291409794926,-2.06193435869265) q[13];
u3(2.05458766292141,-3.70462898498390,-0.880700260231160) q[15];
u3(2.78421977053925,-3.59907819161843,0.985601429779762) q[11];
u3(1.76204601339140,0.0773794114509118,1.93731319583860) q[12];
cx q[12],q[11];
u1(4.60795619304121) q[11];
u3(-3.44656221140404,0.0,0.0) q[12];
cx q[11],q[12];
u3(-0.347550070267594,0.0,0.0) q[12];
cx q[12],q[11];
u3(1.64012998830884,-2.70938102906542,0.768497501362105) q[11];
u3(0.299907274355949,1.39972524728436,2.44698106303062) q[12];
u3(1.66115979924856,3.21092409883444,-2.22435423831668) q[1];
u3(0.771596057683218,2.71102906162046,-2.71945439252349) q[9];
cx q[9],q[1];
u1(2.32830755585251) q[1];
u3(-2.70646834164841,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.10138919753433,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.21139637954217,-1.66926893792141,4.31939067738985) q[1];
u3(1.02320254344859,0.523819484703736,4.56245751736072) q[9];
u3(1.69778802096454,0.521336733340395,1.48341952795449) q[7];
u3(0.925316004790626,-2.11787243500685,-2.84893156511016) q[4];
cx q[4],q[7];
u1(0.0293256094233738) q[7];
u3(-1.27231995082054,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.12416180757282,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.795184524540416,2.93247514356366,-2.67944808026213) q[7];
u3(0.313980639914606,-2.92336798501774,-1.89723789092099) q[4];
u3(2.87777892105934,0.771052282837065,-1.82403602549914) q[10];
u3(1.97658059627391,1.89402108157188,-4.31178344232357) q[0];
cx q[0],q[10];
u1(2.36064350616513) q[10];
u3(0.402427800225537,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.91604338399793,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.790174605498226,-0.503332598819395,0.227037629178746) q[10];
u3(0.525820395773640,-1.37579206023742,3.32033746448065) q[0];
u3(2.23890877386338,-1.42176422641731,0.338686739802440) q[7];
u3(1.72360734049412,-4.00357576609006,1.43265875760489) q[11];
cx q[11],q[7];
u1(0.155803708617796) q[7];
u3(-0.885053530359482,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.62283280182743,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.06741702065621,-3.40640429483485,-0.364156298798783) q[7];
u3(0.974822628815564,-4.05981294018215,-2.20855687214950) q[11];
u3(0.901516419506772,2.55056349966024,-3.11149353975296) q[8];
u3(0.842794200426441,-3.05626175288248,2.43388400624098) q[1];
cx q[1],q[8];
u1(1.30510951468556) q[8];
u3(-0.0253159728770069,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.89165417659185,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.28505815259073,-1.30406584446846,1.82155350382861) q[8];
u3(1.98483632970816,-1.15217854652565,-0.353324928780376) q[1];
u3(0.624852603712448,1.72327256122661,-1.13783116228375) q[6];
u3(1.58167776499951,0.0705963543444263,-3.52266165250996) q[12];
cx q[12],q[6];
u1(2.72107700901480) q[6];
u3(-2.19811936626999,0.0,0.0) q[12];
cx q[6],q[12];
u3(1.44950016376220,0.0,0.0) q[12];
cx q[12],q[6];
u3(1.65584384463094,4.35289498658955,-0.799146605373752) q[6];
u3(1.52746786759727,0.0311650289872167,5.27418798318840) q[12];
u3(2.02408527150841,2.78419740097722,-0.929569109631295) q[9];
u3(1.55398964792844,1.42321687068691,-1.67909439367486) q[15];
cx q[15],q[9];
u1(3.61687056238990) q[9];
u3(-1.37522479130319,0.0,0.0) q[15];
cx q[9],q[15];
u3(2.28937597833411,0.0,0.0) q[15];
cx q[15],q[9];
u3(1.85381299541404,-0.184274118373315,-3.79551847846207) q[9];
u3(1.84849725264031,4.78122860312878,0.0684482038276775) q[15];
u3(1.35846301152315,-0.940788646839223,1.55592598179385) q[4];
u3(1.92141818939527,-1.25365127491940,-1.92633162909383) q[5];
cx q[5],q[4];
u1(2.30902145542051) q[4];
u3(-0.0545151504847183,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.49421080428911,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.54387968616666,0.255342323012584,2.67087246032343) q[4];
u3(1.34491524194307,1.44511071586630,-3.89136183775513) q[5];
u3(2.21360580032819,-0.268410206720035,-1.13141321970023) q[0];
u3(0.245480911141532,0.381331311518508,-5.46639179693046) q[10];
cx q[10],q[0];
u1(-0.975129899235514) q[0];
u3(0.396223470801109,0.0,0.0) q[10];
cx q[0],q[10];
u3(3.82121711374354,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.13460381816161,-3.82535479274959,1.43952092669028) q[0];
u3(0.990150745913675,0.550120599987318,1.23890869958287) q[10];
u3(0.797349125029382,-0.968015822801437,1.43867284548958) q[14];
u3(1.34644698797027,-1.24774338064230,-1.60086545020264) q[13];
cx q[13],q[14];
u1(2.98532951017292) q[14];
u3(-1.80132347326181,0.0,0.0) q[13];
cx q[14],q[13];
u3(1.14231420576782,0.0,0.0) q[13];
cx q[13],q[14];
u3(1.92014030683492,3.91276022312478,-2.32051376783538) q[14];
u3(0.773054599429539,-1.18164057081851,-3.96640279845862) q[13];
u3(2.81413737484773,-0.668626830052841,3.72432969247447) q[2];
u3(1.56135902668738,0.617893132852769,1.77434697278968) q[3];
cx q[3],q[2];
u1(1.84230812411467) q[2];
u3(0.281516170264065,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.848048293539197,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.48651781645817,-0.111919381242826,-0.433072476169268) q[2];
u3(1.21576901832794,2.19223350973510,-0.571035922340056) q[3];
u3(1.82999772405827,0.0625911219404812,1.81534540782594) q[5];
u3(1.44894532395769,-0.335508975237552,-1.32151237839702) q[13];
cx q[13],q[5];
u1(3.13741634508264) q[5];
u3(-2.29562376873149,0.0,0.0) q[13];
cx q[5],q[13];
u3(1.24546980629073,0.0,0.0) q[13];
cx q[13],q[5];
u3(1.00361277320904,2.11788458207895,-0.195034588828298) q[5];
u3(1.50939552169590,0.108292337689535,1.15328660390370) q[13];
u3(2.17530952087162,3.99617519422380,-2.22624384786908) q[10];
u3(0.239477077646155,0.911585110735129,0.0929426497945267) q[11];
cx q[11],q[10];
u1(4.37762323860161) q[10];
u3(-3.69712841705922,0.0,0.0) q[11];
cx q[10],q[11];
u3(-0.389794792156923,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.21481768068069,0.918571738852368,0.121401282243106) q[10];
u3(1.60162925007075,-0.777417144132248,-5.03162414167917) q[11];
u3(2.94648627493307,-3.07770584816155,0.322400752368514) q[7];
u3(2.37778717061962,0.879556211599613,3.23846657219308) q[14];
cx q[14],q[7];
u1(2.69907003018853) q[7];
u3(-1.90634747115888,0.0,0.0) q[14];
cx q[7],q[14];
u3(1.61496381807126,0.0,0.0) q[14];
cx q[14],q[7];
u3(2.85148024732035,0.506874172060849,-2.33897527242695) q[7];
u3(1.45483536711399,0.989503369454094,-4.42032704614907) q[14];
u3(1.90314102311490,1.46672494757597,-2.43656150584925) q[9];
u3(2.88589076694933,1.46927936441124,-3.51651634028489) q[15];
cx q[15],q[9];
u1(3.37807806896575) q[9];
u3(-1.60000362693969,0.0,0.0) q[15];
cx q[9],q[15];
u3(2.36733839736796,0.0,0.0) q[15];
cx q[15],q[9];
u3(2.34914973642561,0.482092080553042,-3.92283252478284) q[9];
u3(1.45782776140121,1.59175260140288,3.57960576190636) q[15];
u3(0.617277781739870,1.80814195310637,-0.128254346228231) q[6];
u3(1.03606205212083,-0.928609189875419,-3.32025300980263) q[0];
cx q[0],q[6];
u1(1.84545029748338) q[6];
u3(-2.71070524814283,0.0,0.0) q[0];
cx q[6],q[0];
u3(3.27712912056131,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.40647398984519,-0.119464011014257,1.66644006093477) q[6];
u3(0.653411508366198,-0.766375069506882,2.34350457351699) q[0];
u3(2.33753474103264,1.80130831830611,0.0328717967374486) q[1];
u3(2.92845630227171,-0.546069145233858,-5.12558598536282) q[12];
cx q[12],q[1];
u1(2.98112960117337) q[1];
u3(-1.54686660120433,0.0,0.0) q[12];
cx q[1],q[12];
u3(0.328192196665416,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.99169842414167,-1.03454718364170,1.63632876700203) q[1];
u3(1.02125942230331,0.764697282038680,0.159040988305264) q[12];
u3(1.68415224489898,-4.02139340937723,2.06234348107942) q[4];
u3(0.348879901788592,0.385982259702543,0.639344367122458) q[2];
cx q[2],q[4];
u1(3.19861508311721) q[4];
u3(-1.18911448755073,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.85723284025434,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.780384449182208,2.17619841840937,-3.53631255217903) q[4];
u3(2.29329236025986,1.26226006479608,-0.431775057976428) q[2];
u3(1.88830193577873,2.08818820753334,0.150890578383365) q[3];
u3(0.785028679264363,-1.19854750182616,-1.64711475499108) q[8];
cx q[8],q[3];
u1(0.408264015412738) q[3];
u3(-1.14163554898620,0.0,0.0) q[8];
cx q[3],q[8];
u3(3.18129084069870,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.40218878200868,3.15946765032525,-1.70149448614383) q[3];
u3(2.24451463631759,-3.04690553866185,0.171108470061066) q[8];
u3(2.34488106225578,-2.20118016167659,1.15389903672116) q[2];
u3(1.86655857155851,-1.97294288682628,-0.353759568407176) q[10];
cx q[10],q[2];
u1(2.86019925025087) q[2];
u3(-2.14493551861830,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.439761046613798,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.43674278939752,-2.92083025365419,1.64056531703097) q[2];
u3(1.69596054311117,1.24708720880863,2.73815516837961) q[10];
u3(2.85129989073158,-1.10616369845849,-1.57439926710200) q[4];
u3(1.11668240813747,-3.97002054813003,-0.0934810107647321) q[14];
cx q[14],q[4];
u1(1.13408681212717) q[4];
u3(-0.955334349959937,0.0,0.0) q[14];
cx q[4],q[14];
u3(3.15673589129660,0.0,0.0) q[14];
cx q[14],q[4];
u3(0.865489056378598,0.0639578978678386,0.415703780772908) q[4];
u3(0.445173191261908,-0.750080739065324,3.23197931562053) q[14];
u3(0.797970214105273,-2.12312683275228,2.44515537035006) q[5];
u3(0.643053593579999,1.95745924980852,-2.92360551663716) q[8];
cx q[8],q[5];
u1(-0.0719113674595311) q[5];
u3(-1.65703943503903,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.714980620475385,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.63509260644156,0.0997044035254834,0.838353432950298) q[5];
u3(1.41730912479083,0.340909694433594,1.34692847687813) q[8];
u3(1.51158222201961,2.16385328643428,-0.641040505449037) q[6];
u3(0.963626154821869,0.706721449869741,-2.71748979194199) q[11];
cx q[11],q[6];
u1(1.54045004748375) q[6];
u3(-0.0985656782491520,0.0,0.0) q[11];
cx q[6],q[11];
u3(2.74427476603198,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.51886665326377,-0.106374254452591,-2.76329180631611) q[6];
u3(2.42681890879366,2.92327223887294,-2.12122048767333) q[11];
u3(1.84315463928588,0.625113639518462,-3.20232207041562) q[9];
u3(2.04954776016948,-1.41954699276074,4.52264594463700) q[1];
cx q[1],q[9];
u1(0.0398191222409467) q[9];
u3(-1.63273755243170,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.90870189078726,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.25932875894076,-0.00934883218461935,3.08109145141544) q[9];
u3(1.18004198912757,-2.61392835598485,-2.48903872960228) q[1];
u3(2.05235427919045,0.252222998380090,2.74666615990355) q[7];
u3(2.22914748528917,-2.01900824129452,-2.22042901214587) q[13];
cx q[13],q[7];
u1(1.80780662152406) q[7];
u3(-2.16179008100003,0.0,0.0) q[13];
cx q[7],q[13];
u3(3.57131075230575,0.0,0.0) q[13];
cx q[13],q[7];
u3(1.83086160840965,-2.28108113877494,2.74776485224812) q[7];
u3(1.16241465199506,-0.915528390452793,0.495651743175900) q[13];
u3(0.464621895348375,-2.50462582543862,1.79686067804812) q[12];
u3(0.688499113486848,1.70166114961542,-4.01037680035562) q[15];
cx q[15],q[12];
u1(2.78567048770344) q[12];
u3(-1.76562633843385,0.0,0.0) q[15];
cx q[12],q[15];
u3(3.11362400219520,0.0,0.0) q[15];
cx q[15],q[12];
u3(1.26332479639381,0.448104603940831,-1.74227029416298) q[12];
u3(0.603055133589497,0.618750143613246,4.81592044587495) q[15];
u3(1.46316462388880,-1.30197630191234,2.74879664048269) q[3];
u3(0.465417094660841,-1.22397388475439,0.166254143040273) q[0];
cx q[0],q[3];
u1(1.93253294000539) q[3];
u3(-2.99499967577523,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.786220979646450,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.949075768281908,-2.82197173022514,1.54367619812540) q[3];
u3(1.96599009740203,2.61159948351320,2.61567201890476) q[0];
u3(1.25535781405800,1.71176925632969,-2.28415578963806) q[10];
u3(1.25543734810474,-2.16547118127114,2.57417130887469) q[7];
cx q[7],q[10];
u1(0.934526123341319) q[10];
u3(-1.52372255110507,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.52781666566772,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.421573134333356,-1.74915224990547,-0.513354684197210) q[10];
u3(1.85113248728137,-0.0974026906721250,-3.32791326848950) q[7];
u3(1.04708887845846,3.05880303056327,-2.21547607719796) q[4];
u3(1.38045243016521,2.21000831343849,-2.24891404879282) q[5];
cx q[5],q[4];
u1(1.39796307157828) q[4];
u3(-0.187596907567198,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.38664820233980,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.06871837992810,-0.291715789223215,-0.729585557644270) q[4];
u3(2.80551614988070,2.39034505816525,2.19841394613308) q[5];
u3(2.00088328129178,-2.41227551855580,-0.237484336669790) q[14];
u3(1.99633706953348,-2.36056685907131,-0.165836847167979) q[12];
cx q[12],q[14];
u1(2.41141824042679) q[14];
u3(-1.69732945086017,0.0,0.0) q[12];
cx q[14],q[12];
u3(1.17731988008958,0.0,0.0) q[12];
cx q[12],q[14];
u3(1.85084490775942,0.0998870995200980,3.19249906780358) q[14];
u3(1.41706020490267,-0.382320262995805,0.894252450645723) q[12];
u3(2.17445221445354,1.55777561813812,-3.07083004631973) q[15];
u3(1.93051557150592,-1.84975238531302,2.87264800559007) q[3];
cx q[3],q[15];
u1(2.02075350839480) q[15];
u3(-1.46575253767952,0.0,0.0) q[3];
cx q[15],q[3];
u3(3.58557632654243,0.0,0.0) q[3];
cx q[3],q[15];
u3(2.00714306094796,1.64028185005982,-1.61647521590139) q[15];
u3(1.11284500301863,-4.55209071154697,-0.0660259649603123) q[3];
u3(1.94721760113717,2.16605703509408,-1.77906904425111) q[6];
u3(1.60516149781079,1.41900853626611,-2.78675652269268) q[13];
cx q[13],q[6];
u1(2.02703672276780) q[6];
u3(-2.28920719311270,0.0,0.0) q[13];
cx q[6],q[13];
u3(0.353403409850088,0.0,0.0) q[13];
cx q[13],q[6];
u3(1.99897947288114,-0.000640692618702698,1.01926226429275) q[6];
u3(1.68699439271277,-4.65202672342553,0.0696511606700700) q[13];
u3(1.36559427490996,3.92319653123413,-1.80163652671782) q[2];
u3(0.475649445440399,1.83567560705593,-2.50850677023689) q[1];
cx q[1],q[2];
u1(3.86211805087726) q[2];
u3(-3.41481380940783,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.901778391727282,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.408190428326485,1.33998276806121,-0.740066066520843) q[2];
u3(1.83296412585590,1.42168044345994,-4.77019950989142) q[1];
u3(2.62415044430934,1.32952371607897,0.752381162191832) q[11];
u3(0.660212299415517,-4.98246833105963,0.684321831524916) q[0];
cx q[0],q[11];
u1(2.13540184923684) q[11];
u3(0.420319884792387,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.42449133298742,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.987662406338083,1.63656551324536,-3.64743677150659) q[11];
u3(1.33843409452978,-2.72385217485898,-0.844272562882039) q[0];
u3(1.77656389875318,-0.538428151626630,0.940994200318721) q[8];
u3(2.04066388053329,-0.849751911272095,-1.40245405966144) q[9];
cx q[9],q[8];
u1(1.11902056793463) q[8];
u3(-0.827597277902025,0.0,0.0) q[9];
cx q[8],q[9];
u3(-0.0942033943720759,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.68152675739087,-2.91824646560639,0.924124068718531) q[8];
u3(2.18800197029854,-2.83945104581998,1.46374077528949) q[9];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[14] -> c[14];
measure q[15] -> c[15];
