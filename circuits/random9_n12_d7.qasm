OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(0.895742625182354,2.14140340331781,-2.99737500322541) q[10];
u3(1.53960493259269,-2.61048900842477,2.80163733591129) q[7];
cx q[7],q[10];
u1(3.07772993483307) q[10];
u3(-2.19561598636994,0.0,0.0) q[7];
cx q[10],q[7];
u3(1.05765810183123,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.27001988580384,0.903784720877567,-0.00708521854281979) q[10];
u3(1.02704130169086,1.10778582506017,-2.76581136308909) q[7];
u3(1.17877220755452,2.12926848474267,-2.10444520241321) q[0];
u3(0.306374177439699,1.23860467679966,-2.12099112957548) q[6];
cx q[6],q[0];
u1(3.19518380999051) q[0];
u3(-2.21423724591289,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.29473579908927,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.44620196852501,1.87469956856663,-0.395930094458715) q[0];
u3(0.766619866582313,2.77610845309857,-1.26895945623149) q[6];
u3(1.14953494903159,1.23517928698246,-3.45989502984544) q[11];
u3(1.43401495350540,3.86836963397809,-2.37931315191147) q[3];
cx q[3],q[11];
u1(2.14397473937981) q[11];
u3(0.544915417816733,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.33417701098837,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.10603536183032,-0.870695130773601,1.56486173056600) q[11];
u3(0.857276499797201,0.0237647658575404,2.78689623149435) q[3];
u3(1.34488799300547,-2.91619250898272,1.12271723634591) q[2];
u3(1.34452150203775,-3.33528222053641,0.0200105863615609) q[5];
cx q[5],q[2];
u1(3.44906234842891) q[2];
u3(-4.22253694713406,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.191620945127764,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.93679016270335,0.159165833969948,-2.51937178716376) q[2];
u3(2.08577237507052,4.78491047320852,1.02306594918464) q[5];
u3(1.11603988467582,-0.343054445215756,1.29623409366306) q[9];
u3(0.760651286480754,-1.29744606439433,-0.363233344616263) q[8];
cx q[8],q[9];
u1(3.35968784856318) q[9];
u3(-2.14377440436614,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.61722509998378,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.50195427897364,-1.59635665263800,2.87709460001037) q[9];
u3(1.61466781933952,2.20804960505612,-3.62193467052832) q[8];
u3(0.130452493619816,2.34761272147515,-1.95399632386907) q[4];
u3(0.748120355802376,-2.88378581615968,1.36519048493579) q[1];
cx q[1],q[4];
u1(1.61976453353207) q[4];
u3(-0.515944958298976,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.34176051493228,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.31272899363452,1.37990715274420,0.614137053130481) q[4];
u3(0.605401684374986,1.90509235481951,-0.211265291414180) q[1];
u3(0.741518853442220,-0.194095175115283,-1.12087298827372) q[7];
u3(1.89423551328277,-3.41455052175339,1.87500893550684) q[1];
cx q[1],q[7];
u1(2.00655545866173) q[7];
u3(-2.23289630609605,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.200549410028184,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.63461275395788,1.00995604015159,-0.0721528575780555) q[7];
u3(1.19757400643907,5.84585202242295,0.295319505527013) q[1];
u3(1.31032780464395,-0.150209813563097,2.26334211645234) q[8];
u3(1.91966889060676,-2.72541749841757,-2.10268882084547) q[6];
cx q[6],q[8];
u1(-0.469711942419823) q[8];
u3(-1.72480089013823,0.0,0.0) q[6];
cx q[8],q[6];
u3(0.835005281667338,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.01495231597633,2.51917262815327,-3.48794499739303) q[8];
u3(0.172004415261370,0.752572417732794,-3.03500387325576) q[6];
u3(0.877076108061353,1.74683352727066,-3.16136859679603) q[11];
u3(1.81302098701451,-2.30638521896591,3.63729630547301) q[9];
cx q[9],q[11];
u1(1.74363941610740) q[11];
u3(-2.31025112382995,0.0,0.0) q[9];
cx q[11],q[9];
u3(3.46010084529672,0.0,0.0) q[9];
cx q[9],q[11];
u3(0.469693541739410,-1.06083448189047,1.97659873083635) q[11];
u3(0.378419940471835,2.41050064824890,3.81388778630499) q[9];
u3(1.47116514376009,-1.47655252539409,-0.970545780947895) q[4];
u3(0.933010166305851,-4.23444266655973,0.578422660928133) q[2];
cx q[2],q[4];
u1(-0.495185453438724) q[4];
u3(-1.65665763523464,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.962276739421456,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.20309666190074,-3.34145541658031,2.73110861797428) q[4];
u3(2.07651091773915,2.60446819438425,2.63554430660285) q[2];
u3(1.29286298080322,1.86383439640031,-1.52382492833218) q[10];
u3(0.493033672080136,-1.19901635336521,-0.325973774856467) q[0];
cx q[0],q[10];
u1(3.04765713450599) q[10];
u3(-2.03406611106534,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.701531061281072,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.30728287795579,1.19803512875077,-4.80698305528461) q[10];
u3(0.737025762734040,1.66200000284964,2.04397477291437) q[0];
u3(2.81817228124929,-2.31019958474411,1.43666568178331) q[3];
u3(2.57780306371953,1.60823177611800,2.78197356804787) q[5];
cx q[5],q[3];
u1(-0.419299257567320) q[3];
u3(-1.72948241237338,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.977099662206709,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.42174427546267,0.0913402878062450,-2.84552411512404) q[3];
u3(0.454317914316878,2.80185015124562,-0.558815938473328) q[5];
u3(2.31340463846039,-0.844645884037669,-0.406188597297486) q[6];
u3(1.23325389489489,-2.93358689875719,-0.551605934357024) q[11];
cx q[11],q[6];
u1(0.985589557986951) q[6];
u3(-0.457226101735456,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.99140358721532,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.90714999588668,1.44580702481516,-2.29784717740298) q[6];
u3(0.785717693971268,4.39620854533494,0.269501404130209) q[11];
u3(2.47361110809579,-1.12957207886343,1.37477722551029) q[7];
u3(2.09105660451685,-3.29632315055272,-0.377886262850463) q[9];
cx q[9],q[7];
u1(-0.167956940970356) q[7];
u3(0.773720835117580,0.0,0.0) q[9];
cx q[7],q[9];
u3(3.76642596132523,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.67585568765006,-0.139442274709670,-2.15974540376841) q[7];
u3(0.959935971075141,-0.148867887260763,-5.69051622849864) q[9];
u3(1.74928636838416,-2.55319991456596,3.60641788549336) q[8];
u3(1.83515306140282,1.61599833054628,-2.04533653572156) q[1];
cx q[1],q[8];
u1(1.78196872476286) q[8];
u3(-2.99619014460132,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.273685207876580,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.03045640032308,-1.52344898693637,-0.237959774406344) q[8];
u3(1.53543433675361,4.61652297101881,0.389236772969058) q[1];
u3(1.46297288036029,0.263630655470558,-0.546484150780437) q[3];
u3(1.71635537545505,-3.94160588259025,1.58936020593889) q[2];
cx q[2],q[3];
u1(0.656163791419108) q[3];
u3(-1.40656827559501,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.73892333007555,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.207273434959826,-1.23053311973681,-3.14487508349348) q[3];
u3(1.07978082193646,-1.06669275487640,3.89022404519456) q[2];
u3(0.961527900169002,2.22408059776456,-1.80763930737799) q[5];
u3(0.987470006710706,1.38969295084858,-1.64645919726351) q[4];
cx q[4],q[5];
u1(3.09043079795612) q[5];
u3(-2.14749886410064,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.324195461721424,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.57982909822838,0.145001575219398,1.00161963225909) q[5];
u3(1.21761523059259,-4.79533084286982,-1.11737544238270) q[4];
u3(1.55204164705674,-1.49431877673554,0.123557729895673) q[0];
u3(0.284128160078442,-1.48083963808404,-0.707850432417577) q[10];
cx q[10],q[0];
u1(0.704993679322050) q[0];
u3(-1.59712912279221,0.0,0.0) q[10];
cx q[0],q[10];
u3(2.73644955747298,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.42210570060514,-0.622833226593028,1.58342004144645) q[0];
u3(1.43255628701817,2.42940792863075,3.01612799749928) q[10];
u3(1.73815243119871,-1.88835562043051,0.705042539339079) q[7];
u3(1.71449300148100,-3.41232548993325,0.780940883530775) q[11];
cx q[11],q[7];
u1(3.68258025381545) q[7];
u3(-3.36624667029120,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.827624440626834,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.83772917705657,-1.29739002195501,2.86613043770731) q[7];
u3(0.0585279792717048,-0.168464568838780,2.18175067743438) q[11];
u3(1.97588418825161,3.04191783182192,-0.133957600105175) q[2];
u3(2.63216296955148,2.55458804632327,-0.954580422323528) q[6];
cx q[6],q[2];
u1(0.142554124332973) q[2];
u3(-1.56397612154434,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.56788128946621,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.26093719841461,-2.84915690440638,1.23312984486143) q[2];
u3(1.63333395938979,1.11127531176998,0.655836505960847) q[6];
u3(0.926792632127337,2.27852189449201,0.169852196160500) q[8];
u3(1.27817754827797,0.783783636436207,-3.40274472878915) q[3];
cx q[3],q[8];
u1(1.45418558324735) q[8];
u3(-0.430293009726671,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.09958918669604,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.29172063443564,-1.54863064396980,3.04465521534378) q[8];
u3(1.81773477409734,-5.11892272464022,0.303983559904674) q[3];
u3(1.98921804761521,-1.62929359342028,-1.49541278121270) q[9];
u3(0.189808378237153,-1.08651791855152,-4.71826294851379) q[10];
cx q[10],q[9];
u1(1.27204614430615) q[9];
u3(-0.780434727808763,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.76764499767309,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.03603495545888,1.97967480084744,-1.13632217467658) q[9];
u3(1.89030645213794,1.34457765292034,-2.24336863901337) q[10];
u3(1.67661670768291,-1.68075985431831,-0.906402222604582) q[5];
u3(0.911856265126666,-3.00262150472088,0.0766379657223721) q[0];
cx q[0],q[5];
u1(2.88899391066134) q[5];
u3(-1.92605980141887,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.690086347883353,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.59708380853791,-2.53412072273734,-0.712418733265807) q[5];
u3(2.85377162727662,-1.62197709012013,-0.290583430206807) q[0];
u3(1.20216832147690,-0.229298444841586,-1.34731116025556) q[1];
u3(1.47989226513091,-3.76548491725623,1.79094128701728) q[4];
cx q[4],q[1];
u1(2.02517936720893) q[1];
u3(-3.02517037074151,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.56015999566104,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.18382264252026,-1.58834364852385,4.25829298983989) q[1];
u3(0.830731170964085,-0.761028008115444,-4.67318349564207) q[4];
u3(1.40809608668122,2.81941202964559,-3.20060017269293) q[0];
u3(1.75735016951468,-3.49032176777615,2.15630338966446) q[7];
cx q[7],q[0];
u1(1.93146396787696) q[0];
u3(0.256203226691708,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.50020950759024,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.59510549691782,-3.39253495865181,0.389636438391492) q[0];
u3(0.920007818985837,2.09790553294345,0.419062454525898) q[7];
u3(2.57702374066238,-0.671534223489656,2.92369843547435) q[5];
u3(2.18480514116195,1.01871349343345,3.21038472808133) q[6];
cx q[6],q[5];
u1(1.08342189401321) q[5];
u3(-3.47900528671823,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.57201528636805,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.61987392854324,0.896155409972905,-1.64407493364785) q[5];
u3(2.29168642155633,-0.204844537742202,-4.31171229566589) q[6];
u3(1.39188972639091,1.65489500489828,-3.13807737799018) q[11];
u3(0.487951199083321,-2.56120485075569,3.23939638101879) q[10];
cx q[10],q[11];
u1(0.228520776117132) q[11];
u3(-1.57197983974259,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.38041287128743,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.73435761218357,-0.389764365467822,2.83724018890034) q[11];
u3(2.83131284910989,-0.242324759097708,-3.61843357049811) q[10];
u3(2.27784102926139,0.0372586430253270,-1.63495321750156) q[1];
u3(2.03151892867160,-3.91628456833016,1.03380475576607) q[2];
cx q[2],q[1];
u1(1.25043319709744) q[1];
u3(-2.97590613780926,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.59349791896167,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.73094878829754,2.03176110651224,-3.37264152184558) q[1];
u3(2.39475093908946,0.304515130907019,-1.94150073025441) q[2];
u3(2.55641474243742,-2.04800949901061,1.91136932318983) q[3];
u3(2.11449986366404,-3.49079863547016,-2.69037694159211) q[8];
cx q[8],q[3];
u1(0.111364404282682) q[3];
u3(-1.29823103989025,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.64621579837386,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.03735502305840,-0.916030535175159,1.16529595348991) q[3];
u3(2.68383433808417,-3.49802256745881,-1.14737150378328) q[8];
u3(1.05541203127718,-1.63933819085480,-0.755813796342483) q[9];
u3(0.860751193117890,-4.09334321897490,0.565388138120916) q[4];
cx q[4],q[9];
u1(2.74052188598028) q[9];
u3(-1.84962255049347,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.0305679352126882,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.61996267942139,1.32138821305580,-2.35709534627108) q[9];
u3(2.52591235435740,0.00582728527271792,4.70961985317259) q[4];
u3(1.73517999722680,0.203218046968082,-3.34188895974617) q[11];
u3(2.04230673780722,3.56696579589945,-1.58165685383153) q[3];
cx q[3],q[11];
u1(2.19597359204471) q[11];
u3(-2.97598834915652,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.31554257195033,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.35738430837465,-0.991496245048665,-0.940708734477684) q[11];
u3(1.85199082131584,-3.80353567494746,-1.50831238532685) q[3];
u3(1.79711170343165,0.614280173503162,-2.67409123704960) q[9];
u3(1.01725608432911,-2.87466181651613,2.71048804363816) q[4];
cx q[4],q[9];
u1(3.27313184454485) q[9];
u3(-2.05247885000071,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.52480099275203,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.14131891222785,3.06891347333416,-1.90447937836833) q[9];
u3(0.478388996344464,4.57748572834166,-0.128046798272952) q[4];
u3(1.38471502049287,2.10043943816041,-3.00580166176057) q[7];
u3(2.12377079635601,-3.07817960464268,3.04710987393952) q[8];
cx q[8],q[7];
u1(3.49387083716470) q[7];
u3(-0.863144076797661,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.61725154427802,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.848570404665825,2.03557788281061,-0.562255993869413) q[7];
u3(2.00625972319072,-0.592471805334160,4.33152006243795) q[8];
u3(2.02611032982163,-1.96797902131991,-0.909314462710312) q[10];
u3(1.55851913838225,-4.58600315208637,-0.251134494018318) q[5];
cx q[5],q[10];
u1(0.924937733309905) q[10];
u3(0.0103713328849984,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.40993824157130,0.0,0.0) q[5];
cx q[5],q[10];
u3(2.20439648295343,1.27321438928539,-0.506027623414655) q[10];
u3(0.662490452392610,-5.81810258399350,-0.358435069272376) q[5];
u3(1.24135633809667,0.165938990991431,0.546986687840103) q[1];
u3(1.16563464759769,-2.00269236032437,-1.90411190062968) q[6];
cx q[6],q[1];
u1(-0.118354031227684) q[1];
u3(-1.64561729452240,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.24045046157326,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.06182565038276,2.21627542681628,-0.914376169768750) q[1];
u3(0.643842107755775,-2.72920326179771,1.70625006970613) q[6];
u3(1.54373595445532,1.34180375556511,-2.38840809384246) q[0];
u3(1.64781704442658,1.65202022371389,-4.39956647449947) q[2];
cx q[2],q[0];
u1(3.16589633210915) q[0];
u3(-1.84992949906583,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.648427852348485,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.461351973462622,1.80505794392926,0.618958879192973) q[0];
u3(1.08571364545365,3.59432952682193,1.40400127333416) q[2];
u3(1.18185296762195,2.47482338971659,-2.33952120736200) q[3];
u3(0.877781639803415,1.16347816010194,-1.89965596994787) q[10];
cx q[10],q[3];
u1(1.11774222370798) q[3];
u3(-0.629584406569081,0.0,0.0) q[10];
cx q[3],q[10];
u3(0.0200045899774881,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.607879039095641,0.932082469778279,-4.39086880249173) q[3];
u3(0.821667154297580,-4.74071339145627,1.35667388069888) q[10];
u3(1.04790347009186,1.31056764588891,-3.97552349779612) q[6];
u3(1.27948589081653,2.08372344458313,-3.07777621077560) q[8];
cx q[8],q[6];
u1(-0.502606915418589) q[6];
u3(-0.0530103021314678,0.0,0.0) q[8];
cx q[6],q[8];
u3(4.03445376631678,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.87743769434425,1.70305949736646,-1.31946643080385) q[6];
u3(1.59522178610342,4.73006718235096,-0.595486835943137) q[8];
u3(0.770648319220984,-3.31074906598232,2.35649121020517) q[11];
u3(0.663142431007823,0.946281425496188,-3.20778334610063) q[5];
cx q[5],q[11];
u1(1.94172014674345) q[11];
u3(-0.153044406900144,0.0,0.0) q[5];
cx q[11],q[5];
u3(0.370007616893850,0.0,0.0) q[5];
cx q[5],q[11];
u3(2.08947553430042,2.09397552340465,-1.54846664256921) q[11];
u3(2.88086421264875,-2.72270707861289,-2.04121636047933) q[5];
u3(1.12457057779448,-0.183750243011917,-1.32914319990038) q[4];
u3(1.46201551429016,0.691638589239981,-5.10856824382695) q[9];
cx q[9],q[4];
u1(0.827622274790819) q[4];
u3(-1.11357776026015,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.35130104735658,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.69687151497552,-3.59215817552260,0.873874448053696) q[4];
u3(2.11910617811199,0.586212759939534,-3.40865128331211) q[9];
u3(0.213631773833858,2.24487680806604,-3.24675070714216) q[7];
u3(0.554673206481997,2.15765790185710,-3.58762042047528) q[2];
cx q[2],q[7];
u1(2.08936051423776) q[7];
u3(-3.20457531123777,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.439330761895129,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.43568798397323,-0.563585133227365,2.39604944896818) q[7];
u3(2.64300145351796,-3.27811729074127,-0.676126015088869) q[2];
u3(0.871924567575331,1.32189306001708,-2.41243813507007) q[1];
u3(0.437198054160613,-0.828580423505703,-0.410774714719117) q[0];
cx q[0],q[1];
u1(1.90714485598436) q[1];
u3(-3.13015594218361,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.910521241518355,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.66614892746930,-0.201291386208700,-0.398325869124997) q[1];
u3(2.44508158349964,3.47274449387664,-2.64443491790702) q[0];
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
