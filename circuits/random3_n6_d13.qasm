OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.03063879616370,-1.49817692998028,1.36167351313310) q[3];
u3(0.579947297205646,-0.544002401824284,-1.06182355044226) q[5];
cx q[5],q[3];
u1(-0.577316184199975) q[3];
u3(1.28667592013747,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.53434792391835,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.20520701872380,-0.440794731553415,-0.354949236531270) q[3];
u3(0.804238491985837,2.75691652693942,-3.36927663031045) q[5];
u3(0.612363633141781,2.60638662438519,-0.499545656871176) q[0];
u3(2.13892907920407,0.286097683485793,-1.58078189473794) q[1];
cx q[1],q[0];
u1(2.72263741358506) q[0];
u3(-1.96604806545556,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.811046778399414,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.74985391198237,3.91331607582380,-2.01620934933618) q[0];
u3(0.429786399908232,3.38387695034569,-2.46534128659272) q[1];
u3(2.66835883959266,-1.38397179876570,4.30790506111825) q[4];
u3(1.20658548954035,0.750304616187724,0.239997741190035) q[2];
cx q[2],q[4];
u1(2.94319935141545) q[4];
u3(-1.87034009909610,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.812762131939457,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.813657819018101,-1.10514414914825,2.15430128855886) q[4];
u3(0.990887689712163,1.56052209488135,3.07097150027719) q[2];
u3(1.18119667964802,1.37892251488211,0.208175288105476) q[2];
u3(0.187723574429124,0.444139816470182,-2.06110679503671) q[1];
cx q[1],q[2];
u1(0.170576953328509) q[2];
u3(-1.53220890163320,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.33425846433274,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.603909673885362,-1.38082698795577,0.0184571939320753) q[2];
u3(1.85052540774879,-4.97112553431527,0.0664089493262288) q[1];
u3(1.19111033445491,1.02034754399377,-2.30710821159391) q[0];
u3(1.33319393862272,-2.27387766891486,2.50938968859485) q[5];
cx q[5],q[0];
u1(3.29068987869561) q[0];
u3(-0.754569126201199,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.84539929523267,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.775303914325097,-2.17053050965642,1.16636083771548) q[0];
u3(1.87358613126330,1.73768925984850,4.47206653660862) q[5];
u3(2.12039723427956,-0.745304710365083,2.34396258230098) q[3];
u3(1.31024691255384,-1.62178790760416,-2.03284930240104) q[4];
cx q[4],q[3];
u1(0.384802363225496) q[3];
u3(-1.91567276632015,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.0264196139720039,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.65785228566077,0.827097534388636,0.0453425234341368) q[3];
u3(1.65418163620141,-0.309251626056493,-2.04122506282855) q[4];
u3(0.665092341185723,2.58247482948404,-0.617908880184192) q[2];
u3(1.61671244630284,-0.561240451714428,-3.57928018483738) q[5];
cx q[5],q[2];
u1(1.18759223454438) q[2];
u3(-0.699675730421454,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.02812516541452,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.611833382488006,-1.66058765386402,-0.145030202319466) q[2];
u3(2.73316557752221,3.51422708833342,2.09892484939140) q[5];
u3(0.806641549339944,1.86491024189580,-1.29244103970209) q[3];
u3(0.866309170759578,0.277346069823786,-2.72748284155838) q[4];
cx q[4],q[3];
u1(2.50683263230069) q[3];
u3(-1.88325950488677,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.14509671146015,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.68133686110097,1.34389780137306,0.218136825570174) q[3];
u3(1.20292728411799,-1.06677955542220,3.21274154226169) q[4];
u3(0.546262155985750,-0.785786231768347,0.798646829941245) q[0];
u3(0.327736111105775,-1.32171172210506,-0.622614871529620) q[1];
cx q[1],q[0];
u1(0.634941999136818) q[0];
u3(0.0295758043399421,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.96920420547255,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.16003237083309,1.91523706132393,1.13964882168697) q[0];
u3(0.981051687140009,-1.65835192423265,-2.77359531837120) q[1];
u3(1.69128792934416,-0.169319116763247,1.29577230996914) q[1];
u3(1.67033904528405,-1.71810568784574,-1.02984564750978) q[2];
cx q[2],q[1];
u1(3.13752683731280) q[1];
u3(-0.804147689979207,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.90228394564080,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.54141939490764,-0.252982881820162,3.16987178137559) q[1];
u3(1.88769946950870,0.239356443475459,-1.35532803197470) q[2];
u3(1.39988447644634,0.118663007347513,-0.444544796907107) q[3];
u3(0.224085862111651,1.44596580404114,-4.39463627005382) q[5];
cx q[5],q[3];
u1(1.26932802777217) q[3];
u3(-3.28347455637338,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.01822157637837,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.633069532363908,1.37253573686845,-2.39528103832299) q[3];
u3(1.33907979272995,-5.17128887860475,0.00546908451196337) q[5];
u3(2.95601842189375,-0.801220376785500,1.44858289128621) q[0];
u3(2.40039627610025,-0.807151000946992,0.627664697518232) q[4];
cx q[4],q[0];
u1(-0.648729954903527) q[0];
u3(0.330281813276885,0.0,0.0) q[4];
cx q[0],q[4];
u3(4.20513329343863,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.08539071990596,-2.46935646368475,2.34264570609748) q[0];
u3(1.64603284410959,-1.33442332994878,1.01201372324567) q[4];
u3(1.53421601457319,-0.226880010088181,-1.08704381316408) q[0];
u3(1.83167958569949,-4.73243165403343,1.13519286188129) q[3];
cx q[3],q[0];
u1(1.29698924375786) q[0];
u3(-0.958354812974608,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.172128578870669,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.78210664227272,-0.333851651707314,2.47647318797849) q[0];
u3(0.215128736277030,-1.78907696283973,-4.29005714531347) q[3];
u3(1.50590248224188,0.851073689138941,-0.407941948645564) q[5];
u3(2.63192998146232,-4.16729815133352,1.64272421156710) q[2];
cx q[2],q[5];
u1(2.17400227709847) q[5];
u3(-2.99716433859945,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.17068721443366,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.611346088555345,3.54495566293322,-1.87232369740174) q[5];
u3(0.689213768616554,1.42695813186412,2.34697292010718) q[2];
u3(2.43503848602132,1.60969302459410,-2.19778106987456) q[1];
u3(1.66650153124938,1.81415122035536,-2.58760523448675) q[4];
cx q[4],q[1];
u1(2.69155506598367) q[1];
u3(-2.05537701573680,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.189173945468128,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.59582582420556,-0.520229307209546,4.02183100479389) q[1];
u3(2.71951894661608,2.54989048324451,-2.53385908500874) q[4];
u3(2.30639228216848,2.10128166830712,-2.52889948520416) q[4];
u3(1.84434696646126,2.69723837192145,-3.30422150061227) q[5];
cx q[5],q[4];
u1(-0.0712725666658560) q[4];
u3(-1.18123797208836,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.53159356206810,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.02175147721528,-1.19047549241610,-1.96223716804364) q[4];
u3(2.06892195441146,3.56573171962211,1.11172284488174) q[5];
u3(1.36037893087204,-2.09941675343137,1.70531316690845) q[1];
u3(0.324632701716033,1.54306167921429,-3.37808275824773) q[3];
cx q[3],q[1];
u1(2.96286201604420) q[1];
u3(-1.28739917487025,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.76795327488988,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.31413201705175,-0.433407555535941,0.534639288154380) q[1];
u3(2.41261807879581,-3.96353896041086,-0.243639579002982) q[3];
u3(2.16353040259719,0.775635522364342,0.336257239785612) q[0];
u3(0.379750259094736,-3.46716371331111,-1.64619064880248) q[2];
cx q[2],q[0];
u1(1.68416893832854) q[0];
u3(0.271580954315506,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.966662533955419,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.70636013040880,0.216538901824107,-1.22571560613126) q[0];
u3(2.11057422892622,-1.69919801085356,-4.27878975408220) q[2];
u3(2.31271487487610,0.831345090667768,1.54400428564243) q[3];
u3(2.12200922345952,-2.02770933497448,-1.61076188169286) q[1];
cx q[1],q[3];
u1(0.0301112792150071) q[3];
u3(-2.22083847017723,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.57450968424322,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.36928371621295,1.42054748856121,-3.02251509286971) q[3];
u3(0.814587623598650,-1.85725587737379,3.10971033820100) q[1];
u3(0.735286500606721,0.210591427287138,-1.95077365466462) q[0];
u3(1.73736180909073,-2.69347458664769,2.54721154442205) q[2];
cx q[2],q[0];
u1(4.56397544797220) q[0];
u3(-3.27251274970551,0.0,0.0) q[2];
cx q[0],q[2];
u3(-0.0865052951492591,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.77776397724557,0.227795995781001,-3.58301253348097) q[0];
u3(1.49630901295680,-1.86619450394433,2.40491701878301) q[2];
u3(2.63116536944995,0.680966596873176,-0.660295251115321) q[4];
u3(1.84186277365857,-4.22914226681385,1.74840196085961) q[5];
cx q[5],q[4];
u1(-0.472722451324552) q[4];
u3(-1.75971297624566,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.16758426934906,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.92851334438043,0.761131684271152,2.44115363436541) q[4];
u3(1.49370411030297,-2.50393162194104,-2.29450897561270) q[5];
u3(1.84345347542652,1.06477972936279,0.375543784105952) q[2];
u3(2.10747188805177,0.162107622949606,-3.00496844043017) q[1];
cx q[1],q[2];
u1(-0.0967873492712479) q[2];
u3(-2.41716968877425,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.27024622975169,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.38773991108353,-1.22539389884159,-0.297079500860706) q[2];
u3(1.36486778069339,-1.28308217776736,4.75516465570202) q[1];
u3(2.03925188062832,-3.38168153637824,0.465212084100064) q[3];
u3(2.67726223434792,-3.41657327700218,-1.68268254588350) q[0];
cx q[0],q[3];
u1(2.80026536058782) q[3];
u3(-2.44045984831750,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.24633532865559,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.782461906191968,0.351688027314799,-2.71169064665185) q[3];
u3(1.93195861724795,-4.25370006523957,-0.0571704333690466) q[0];
u3(0.243176076309662,3.00888322215122,-2.95774579861431) q[5];
u3(0.453856033228379,-4.01064430803183,1.82587559501300) q[4];
cx q[4],q[5];
u1(2.69671320418415) q[5];
u3(-2.49308131304619,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.30697561643326,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.32596198481272,2.87329491851688,-3.28573374942646) q[5];
u3(2.63705378535613,4.38690177364888,1.22361821808451) q[4];
u3(2.49064479479189,1.41768410375311,-0.374559928517660) q[1];
u3(1.86349259249718,-0.819085730978120,-2.31228770260874) q[0];
cx q[0],q[1];
u1(1.37444543742200) q[1];
u3(-3.17074155691021,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.67824308732033,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.32666735092467,-2.63871212813050,-0.902137425780071) q[1];
u3(1.51451603838179,-0.687709367287001,-3.48995346350900) q[0];
u3(2.25756290907156,0.819976449689760,-2.59312547250190) q[5];
u3(2.72999129006969,0.604033258421969,-3.53794794798338) q[2];
cx q[2],q[5];
u1(2.48574540786126) q[5];
u3(-1.68628094404864,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.325872052222514,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.22152528087252,0.390544180216762,-0.510080320457903) q[5];
u3(1.41576186588269,-2.18010109657717,0.938759223937874) q[2];
u3(2.88455900463927,2.23730327627659,-1.43830842718339) q[3];
u3(2.02560761054770,5.81803164532961,0.377203817507759) q[4];
cx q[4],q[3];
u1(3.59028974526802) q[3];
u3(-0.831736020253621,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.84909121452321,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.24658331525336,-0.616094483522906,-0.851751549426337) q[3];
u3(1.72913957552334,3.46377818380353,-1.79102680544098) q[4];
u3(2.39693709269794,-2.15999324283368,-0.0343785473004332) q[0];
u3(2.54234121834838,-0.947597537216836,0.180818150382600) q[2];
cx q[2],q[0];
u1(1.33576899464619) q[0];
u3(-0.124137544438912,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.30283958155382,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.08119935684130,4.25112588257044,-0.604287775497271) q[0];
u3(1.76809093697544,-2.11579132365580,-0.336751930759146) q[2];
u3(1.57872287842450,1.27601630823999,0.327306158781250) q[4];
u3(1.47043637582081,-1.10904235025102,-1.94800461008497) q[1];
cx q[1],q[4];
u1(1.94411793852769) q[4];
u3(0.616894924519023,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.54677104921213,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.41008693797653,-2.01401621272700,3.74254921173920) q[4];
u3(2.55497871423500,2.72357861530910,3.38538648316549) q[1];
u3(2.73648181937715,0.103694018149219,-1.06610145434214) q[5];
u3(1.22661035251335,-4.52853780413003,1.10072038749557) q[3];
cx q[3],q[5];
u1(0.813587758978021) q[5];
u3(-3.76716128849019,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.49364542202839,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.58486522417825,1.67854836140999,1.41247526536595) q[5];
u3(0.642083255346589,0.649894907442671,-2.97811651684835) q[3];
u3(2.29664643715570,0.606238085684985,1.00883859489314) q[3];
u3(1.51978471103933,-2.04258074782828,-2.25075331970641) q[4];
cx q[4],q[3];
u1(-0.119432915333657) q[3];
u3(-1.54120114647216,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.34960662988460,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.758060368659887,3.33918717324669,-1.32792950981733) q[3];
u3(2.10530548907020,1.41620413297285,-3.14293754224173) q[4];
u3(1.52145128601346,-0.318222459547830,-2.22971124440963) q[1];
u3(1.88598037744869,-4.82638113634347,1.45520764662170) q[5];
cx q[5],q[1];
u1(2.82772377524095) q[1];
u3(-1.88736620209039,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.435405242126769,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.606532811303204,-0.835791858127700,-2.35971897128556) q[1];
u3(0.894942580127106,0.795604559943305,-5.32932837316365) q[5];
u3(1.14027150504680,0.638154708271432,1.37894463119207) q[2];
u3(1.28219251600077,-1.76701298959563,-0.883932487479256) q[0];
cx q[0],q[2];
u1(1.57253546786315) q[2];
u3(-2.38622492613238,0.0,0.0) q[0];
cx q[2],q[0];
u3(3.61274239106019,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.35751659438174,-2.60831673905111,1.49285855497762) q[2];
u3(1.50664375264694,-2.74867585596583,-3.51116073516020) q[0];
u3(0.892967579270897,0.891818270672039,-0.754625819805585) q[1];
u3(0.779107371241787,-0.893716874173158,-0.731598224577290) q[3];
cx q[3],q[1];
u1(2.75844597660299) q[1];
u3(-1.72786967601619,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.724815902785648,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.67170143215410,-0.993112757117689,3.08598796810283) q[1];
u3(0.853324611141607,-0.00623360379696458,1.61075229622353) q[3];
u3(1.10012807243428,1.67731799795997,-2.89389597547853) q[2];
u3(1.85676387203925,1.81448785970725,-4.20120254932954) q[4];
cx q[4],q[2];
u1(0.316533938931109) q[2];
u3(-0.803489552095623,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.94804291482244,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.846428476191429,-0.772395020515125,-0.576310705116149) q[2];
u3(2.30547271874751,-0.464811430573465,-1.09031605541760) q[4];
u3(2.24559009314830,-1.68905291512170,3.80922977280439) q[0];
u3(1.78535455787042,1.39954712040575,1.29542340199708) q[5];
cx q[5],q[0];
u1(2.98014264205524) q[0];
u3(-2.67365350172724,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.09061187699226,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.55859298853450,0.246199789090507,-2.56409457406741) q[0];
u3(1.45585917224860,-1.46693990235758,3.21153562954309) q[5];
u3(0.0576021954348316,1.63550691499212,-2.61338839763243) q[4];
u3(1.59973698801713,-2.74883060927941,2.93408943371778) q[1];
cx q[1],q[4];
u1(0.461316330108869) q[4];
u3(-0.801322591856301,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.63500280101546,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.174798945698361,-3.88190259031433,2.10852369339065) q[4];
u3(2.06571779813494,0.695935563673103,0.527808963222786) q[1];
u3(1.56242255503746,-1.07470542818276,1.45345337400999) q[0];
u3(2.21379196134935,-1.50624281677896,-2.49859210982732) q[3];
cx q[3],q[0];
u1(1.34167180407548) q[0];
u3(-3.25451865166827,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.04675806960254,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.668717054005546,1.22631592502763,1.18270441696311) q[0];
u3(2.53980461231159,-1.98540320387022,-2.39108779396539) q[3];
u3(1.24091309671453,0.886361942404789,-1.26824680178988) q[5];
u3(0.242582522531146,-4.57709568466752,1.64346725784772) q[2];
cx q[2],q[5];
u1(1.74156487324543) q[5];
u3(-2.61527381834421,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.690490050965695,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.10214995556997,0.852115203250198,-2.13224790227300) q[5];
u3(1.84070835154839,-0.502963321350241,0.308184636506328) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
