OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0.0725498862285003,2.43027753711774,-2.66045473703913) q[13];
u3(0.822722545518041,-1.26113905491063,-0.814550248158256) q[4];
cx q[4],q[13];
u1(0.976771946146318) q[13];
u3(-3.07113133949348,0.0,0.0) q[4];
cx q[13],q[4];
u3(1.99885206639770,0.0,0.0) q[4];
cx q[4],q[13];
u3(1.64491994837961,-2.50113175402754,1.49951906030307) q[13];
u3(0.786380715079279,-3.67806687618441,-0.695714157059719) q[4];
u3(0.584041526685607,0.350804551036941,-2.02482231671986) q[8];
u3(1.58287035585125,-3.75918891216640,1.76862713815123) q[7];
cx q[7],q[8];
u1(2.33478714041152) q[8];
u3(-0.00289112278542714,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.50662846213753,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.520169217900182,-0.657947181177980,1.61193467757673) q[8];
u3(1.39091999341888,4.37338024471304,-0.779876420154217) q[7];
u3(1.73604753692523,1.88567462366334,1.04004718592291) q[9];
u3(2.04919317204908,1.11839101131268,-2.90465215191429) q[1];
cx q[1],q[9];
u1(-0.205966996488485) q[9];
u3(-2.06714048407593,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.17508586696931,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.912698373819766,-0.444982093517660,4.11305596823936) q[9];
u3(2.60758905225596,3.55093176268282,1.30454624775181) q[1];
u3(1.02123566096734,-0.605219252754876,0.839957116702233) q[5];
u3(0.933982354104035,-2.60465618613511,0.120314809755931) q[10];
cx q[10],q[5];
u1(2.53692220575773) q[5];
u3(0.174792985613897,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.84953822683083,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.82512034787168,-3.70231145875810,1.98915614698333) q[5];
u3(0.339437840549598,-1.41787352206117,-2.89562074437417) q[10];
u3(2.45191612717996,-0.645315769493646,-1.98677350268999) q[2];
u3(2.02189706007702,-5.02391384508243,1.07170738661079) q[6];
cx q[6],q[2];
u1(3.52392697271032) q[2];
u3(-1.51236376766474,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.19368236897243,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.428022135687657,1.01551061466342,-3.39064545019518) q[2];
u3(0.856394915488511,-0.760864022467197,-4.11111541772706) q[6];
u3(2.78539180631900,1.92392618692671,1.19542589647619) q[3];
u3(1.49372074643357,-0.0608907025984120,-3.53859662672452) q[12];
cx q[12],q[3];
u1(0.867983574432919) q[3];
u3(-0.334813170048800,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.58697138489813,0.0,0.0) q[12];
cx q[12],q[3];
u3(0.685185094141937,3.58887358111238,-2.27849804989509) q[3];
u3(0.366412973734337,4.59842400803082,0.568654506313695) q[12];
u3(1.15707032959465,-0.502063771116105,-2.06618158215499) q[15];
u3(0.420628974981539,-4.40701460158792,1.23052389279444) q[11];
cx q[11],q[15];
u1(0.346809086237694) q[15];
u3(-1.34607774396564,0.0,0.0) q[11];
cx q[15],q[11];
u3(2.34287707766049,0.0,0.0) q[11];
cx q[11],q[15];
u3(0.294938269261664,3.19289235788430,-1.76949978233969) q[15];
u3(1.45802022446057,-4.38871087775701,0.283315803087175) q[11];
u3(1.35143008054344,3.77652599279888,-1.36600040166847) q[14];
u3(2.28013251778211,2.73445015304530,0.129780956486109) q[0];
cx q[0],q[14];
u1(2.07651500662598) q[14];
u3(-2.48301840599855,0.0,0.0) q[0];
cx q[14],q[0];
u3(1.38900275652366,0.0,0.0) q[0];
cx q[0],q[14];
u3(1.61561113440595,-3.09947527405308,1.07037390448708) q[14];
u3(1.71189843222291,1.17581871198568,-0.273722223706431) q[0];
u3(2.19159605306242,3.66017931932003,-0.548583647241302) q[4];
u3(2.15463171209135,2.63804622993173,-0.237153236926950) q[9];
cx q[9],q[4];
u1(1.80822755251031) q[4];
u3(-2.38336798014526,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.11790705324169,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.20072707765292,-0.915968896041004,1.54020758564128) q[4];
u3(0.947194775455521,2.43988910497004,-2.27810540253422) q[9];
u3(0.848179392988911,2.79377213905280,-3.41544099391738) q[5];
u3(0.866533842866187,0.257517883282231,-2.23234609778646) q[2];
cx q[2],q[5];
u1(3.00338525016494) q[5];
u3(-2.33811075986475,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.68433381271015,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.46423113276761,-3.26662531903635,2.24099418591862) q[5];
u3(1.29173621481006,3.05289116707248,1.19072285372783) q[2];
u3(1.41931475046504,1.43466203197688,1.26526787982774) q[14];
u3(1.86780986010006,-1.75100865749434,-0.867018978806250) q[8];
cx q[8],q[14];
u1(1.23025302216675) q[14];
u3(-0.750616092521363,0.0,0.0) q[8];
cx q[14],q[8];
u3(1.87043668398705,0.0,0.0) q[8];
cx q[8],q[14];
u3(1.98171130059454,-1.87843102826615,1.80163817783772) q[14];
u3(0.505434989535598,0.691222847799070,2.51387863684633) q[8];
u3(0.874651074285146,-0.00675089064049195,1.07826695941081) q[13];
u3(1.32165250754978,-0.593614196882880,-1.75231461320858) q[6];
cx q[6],q[13];
u1(-0.224896339826810) q[13];
u3(0.0865028831218995,0.0,0.0) q[6];
cx q[13],q[6];
u3(4.20027400782562,0.0,0.0) q[6];
cx q[6],q[13];
u3(2.04751715771948,-3.07473855747758,0.684943820691822) q[13];
u3(2.00292190998340,-3.89229013905635,-0.516496143060253) q[6];
u3(2.16528615568830,-1.21331679449751,3.89133097324501) q[1];
u3(1.17100624566165,1.42098564823241,1.20466561743090) q[0];
cx q[0],q[1];
u1(3.06225197090089) q[1];
u3(-1.62458015596767,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.948905955396218,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.33479778891338,0.841313814647387,2.14133695961946) q[1];
u3(0.585021413349184,2.37860885794214,2.65242686353010) q[0];
u3(1.04568406940685,2.16671802430894,-1.27089351150421) q[12];
u3(1.37642208664006,1.28263830244768,-0.454785343020083) q[7];
cx q[7],q[12];
u1(2.74160056720557) q[12];
u3(-2.42177150937152,0.0,0.0) q[7];
cx q[12],q[7];
u3(1.49216156862691,0.0,0.0) q[7];
cx q[7],q[12];
u3(1.77484508244803,0.199062004993799,-1.15427150118554) q[12];
u3(2.05674223314459,3.66527173011611,2.55858775294041) q[7];
u3(1.12619008778060,-1.41309960851141,0.493417021299445) q[11];
u3(0.623744627321040,-2.00399561013920,0.0617327088751878) q[15];
cx q[15],q[11];
u1(0.325969875430769) q[11];
u3(-1.13461214631940,0.0,0.0) q[15];
cx q[11],q[15];
u3(1.69479626609036,0.0,0.0) q[15];
cx q[15],q[11];
u3(2.16469859765009,2.86305186234268,-0.742705214410286) q[11];
u3(2.46898769109283,0.0876482586434029,4.54289753022469) q[15];
u3(1.65969455298400,-0.667612306404493,1.70525235534560) q[10];
u3(1.36438306571165,-2.33135630042966,-2.44854427707827) q[3];
cx q[3],q[10];
u1(1.32610431956949) q[10];
u3(-3.13342376128507,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.75351318492228,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.09949331449306,-0.108555692767959,2.79267504365864) q[10];
u3(1.35981706827140,2.72776355252718,1.95193406690104) q[3];
u3(1.61375628898521,-0.102066760966109,1.79999897809788) q[7];
u3(1.95560940282452,-1.77088583953368,-0.413937501388616) q[15];
cx q[15],q[7];
u1(1.81237874787532) q[7];
u3(-2.78354338672098,0.0,0.0) q[15];
cx q[7],q[15];
u3(0.879370762984449,0.0,0.0) q[15];
cx q[15],q[7];
u3(1.25208159017357,-0.988289435720067,0.999609459858193) q[7];
u3(1.50856166089317,0.102310136253726,-5.00032804571408) q[15];
u3(0.875306458668812,2.71288679663863,-2.60956330260000) q[5];
u3(0.850459412334790,0.724282749644439,-2.38428047222374) q[4];
cx q[4],q[5];
u1(2.49203013155050) q[5];
u3(-1.83301457938280,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.16779774856326,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.345957840059587,2.95676751915305,-1.46408276875217) q[5];
u3(1.47300565510306,-2.18689926031912,-4.02088289033867) q[4];
u3(0.677828910715681,2.32447609133191,-2.24772728975604) q[6];
u3(1.00906697017336,-3.30722122500301,1.73329389613593) q[3];
cx q[3],q[6];
u1(-0.146376786214500) q[6];
u3(-1.53631223214197,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.24563778869088,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.25992855087614,-1.75050971738270,1.55472542516841) q[6];
u3(2.96191795323245,-2.35055732091969,0.152148468773224) q[3];
u3(2.17108827678948,0.130297703190348,2.96604566089782) q[8];
u3(2.27476747182743,-2.42145506421972,-2.29742625070471) q[14];
cx q[14],q[8];
u1(3.36762085082162) q[8];
u3(-1.73616771751322,0.0,0.0) q[14];
cx q[8],q[14];
u3(1.50326169241136,0.0,0.0) q[14];
cx q[14],q[8];
u3(2.90001530114541,-1.94224271445539,4.02995474211199) q[8];
u3(1.96907659882145,-1.06564051905526,3.36958757506656) q[14];
u3(1.02894698962646,-0.140141739718825,1.73258024470982) q[1];
u3(1.85827569753705,-0.446603774681632,-1.85212440654788) q[10];
cx q[10],q[1];
u1(0.692493997034070) q[1];
u3(-1.53285393147369,0.0,0.0) q[10];
cx q[1],q[10];
u3(-0.603243586039452,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.92916910445270,-3.00473313821537,2.32019837813280) q[1];
u3(1.10680581525964,-0.779819506903533,2.67510520027825) q[10];
u3(2.43243173158686,3.57997678342136,-1.27981659078779) q[11];
u3(0.947790784977903,1.63202525261122,-1.14211493894954) q[13];
cx q[13],q[11];
u1(0.405960051392852) q[11];
u3(-0.836870427552844,0.0,0.0) q[13];
cx q[11],q[13];
u3(2.01628190773220,0.0,0.0) q[13];
cx q[13],q[11];
u3(1.68706121202853,2.88068664481222,-1.27179359048774) q[11];
u3(0.311975103599633,-3.98348024723468,0.361935177528555) q[13];
u3(2.84978373210370,-0.136107244916739,-1.88277751120341) q[12];
u3(1.64479983109393,-4.88330287681739,0.930013093979574) q[9];
cx q[9],q[12];
u1(-0.365482012875999) q[12];
u3(-1.69017557883321,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.03289593939880,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.86842657482079,-0.609133003431339,-0.274536105620518) q[12];
u3(2.67720557240515,1.14341504879336,0.489414044735747) q[9];
u3(1.26407341861914,2.60747638301790,-1.08361327454605) q[0];
u3(1.48322159549063,1.74860154416243,-0.565278591645363) q[2];
cx q[2],q[0];
u1(1.63399205096384) q[0];
u3(0.294939105749761,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.718722730873118,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.02812506206420,0.647074540826095,2.29630416144847) q[0];
u3(1.52980735912794,-3.63215391838078,0.0840943216082801) q[2];
u3(2.83137281717506,1.18091500330016,1.48373799033314) q[14];
u3(1.41229806484095,-2.94296515648814,-2.97796707627798) q[12];
cx q[12],q[14];
u1(3.04188089183285) q[14];
u3(-1.24345287729073,0.0,0.0) q[12];
cx q[14],q[12];
u3(2.24956637568719,0.0,0.0) q[12];
cx q[12],q[14];
u3(2.68271732762616,1.37004846415805,-1.16790126067292) q[14];
u3(2.88603785552937,5.72133880584052,-0.117624409929127) q[12];
u3(2.50624039536075,1.23552872538149,-0.443476208639151) q[1];
u3(1.51652776404745,-0.207907437575805,-4.03645195719710) q[15];
cx q[15],q[1];
u1(1.62991252517951) q[1];
u3(-3.23581109263614,0.0,0.0) q[15];
cx q[1],q[15];
u3(2.74442722040435,0.0,0.0) q[15];
cx q[15],q[1];
u3(1.27816787910010,0.610848976896699,1.76627276666886) q[1];
u3(1.65173905848686,3.30153980202789,-0.241002618070314) q[15];
u3(2.01375097898309,-0.114513782660636,0.521119982629609) q[7];
u3(1.92793399372421,-1.09845774395881,-1.38914807926577) q[13];
cx q[13],q[7];
u1(-1.23138639555639) q[7];
u3(0.773521942317025,0.0,0.0) q[13];
cx q[7],q[13];
u3(4.00439359484769,0.0,0.0) q[13];
cx q[13],q[7];
u3(0.401849324759197,1.35141757457536,2.40220431918319) q[7];
u3(0.487639357857643,-0.258868654221970,0.546095652385211) q[13];
u3(1.65286762131538,0.191188869947077,1.95684323325240) q[6];
u3(1.57720535824092,-2.42050737197855,-2.36908983357503) q[8];
cx q[8],q[6];
u1(1.69398126160222) q[6];
u3(-0.197187316818856,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.728579873932224,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.44503135316242,4.57892484717074,-1.47580746769843) q[6];
u3(2.19813674387634,-4.34393433697929,1.45618270765328) q[8];
u3(1.02296754686854,-3.33662689442606,2.90769684046756) q[11];
u3(1.49066821338563,-2.97125049289705,2.38935113830506) q[10];
cx q[10],q[11];
u1(3.05994375259920) q[11];
u3(-2.56742546320163,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.25645242215495,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.23280548288758,-1.91217266566022,-0.667976167771257) q[11];
u3(1.25473675687783,0.366910711732672,5.87381064065409) q[10];
u3(0.670463312824913,2.83915444350626,-2.28139659484695) q[5];
u3(0.367298391449669,0.783998683644464,-1.77928356452301) q[9];
cx q[9],q[5];
u1(2.28783308624184) q[5];
u3(0.399805316818963,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.44633212168389,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.29903687572956,-1.35310081874995,-3.13329249134977) q[5];
u3(2.40549833472142,-2.75677368269142,1.71767710745466) q[9];
u3(1.65649046677560,-2.11987364359176,-0.973809976153254) q[0];
u3(1.26040050497625,-3.96010941799183,0.284261358026242) q[2];
cx q[2],q[0];
u1(1.29654362071106) q[0];
u3(-3.22225701692706,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.24763941510178,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.93866461871083,-0.837718209772987,-0.0118675083418709) q[0];
u3(1.53271335276716,0.731260981252975,2.77311292528952) q[2];
u3(1.89121826155431,-0.225035810917771,2.45578392546021) q[3];
u3(2.63149299240089,-2.98191572804201,-1.71464455542598) q[4];
cx q[4],q[3];
u1(2.05613935905255) q[3];
u3(-2.48568389490337,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.178289630883965,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.823609808535451,0.360289399723016,0.234363939918692) q[3];
u3(2.38818151528009,3.10317879358450,2.36068552844656) q[4];
u3(0.857696905898136,0.605332860701512,-1.38022199519833) q[0];
u3(1.69847521515312,-4.33439564205319,1.57068903471723) q[11];
cx q[11],q[0];
u1(3.71814721025924) q[0];
u3(-3.35876234245450,0.0,0.0) q[11];
cx q[0],q[11];
u3(-1.11682547148430,0.0,0.0) q[11];
cx q[11],q[0];
u3(0.502395270326148,0.161957789029180,-1.75687943581443) q[0];
u3(2.25491591893328,1.35223739778435,-2.62487390023728) q[11];
u3(1.51399964546530,-0.413462576184036,-1.28926892677114) q[7];
u3(2.41882813712116,0.958667151273797,-4.86309301884724) q[2];
cx q[2],q[7];
u1(2.91244209559463) q[7];
u3(-2.76811566712582,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.03298203446969,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.11180332858967,-2.52579224893084,3.11131786625228) q[7];
u3(0.702053559547278,-5.67746554248239,-0.406961073509285) q[2];
u3(1.94066493188456,0.485053478288504,-1.10238880492268) q[10];
u3(2.60509037675344,-4.87909531674626,1.17497719920284) q[5];
cx q[5],q[10];
u1(3.43262013040269) q[10];
u3(-0.833218778287648,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.08905998560824,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.609701906329386,-2.10461980432466,-0.467030169429049) q[10];
u3(1.70925047464406,-0.970090648897779,-0.0643527190946054) q[5];
u3(1.92590867430818,-0.628774956316226,1.86489076304706) q[13];
u3(1.96883654170736,-2.23832769963846,-0.314881898160473) q[3];
cx q[3],q[13];
u1(3.95845531017623) q[13];
u3(-4.18144124260652,0.0,0.0) q[3];
cx q[13],q[3];
u3(-0.637538333673901,0.0,0.0) q[3];
cx q[3],q[13];
u3(0.810087487003517,1.03140356239705,-3.58589427104088) q[13];
u3(1.34860651310257,4.01970335043123,-1.29350125781678) q[3];
u3(2.71984205231028,0.0193548529103067,2.84880900902598) q[15];
u3(2.89775846775016,-3.42071496086826,-2.43180548237181) q[12];
cx q[12],q[15];
u1(0.723315008438059) q[15];
u3(-0.330882037321779,0.0,0.0) q[12];
cx q[15],q[12];
u3(1.74385154803643,0.0,0.0) q[12];
cx q[12],q[15];
u3(2.43805505544834,0.543805530039838,-0.348116589723089) q[15];
u3(1.77875063080164,-3.98539721781728,-0.696835165897050) q[12];
u3(1.57382330768690,-0.639394207937169,1.42244101908561) q[1];
u3(2.07804578957540,-2.15331174444838,-0.821905472302067) q[8];
cx q[8],q[1];
u1(0.0578049278841577) q[1];
u3(-1.21877321046819,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.99425011759737,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.762359643095942,0.404163219518942,3.58724314660565) q[1];
u3(0.576949418255406,-2.51058585653965,-2.20336141701838) q[8];
u3(1.12261293723640,-1.09354730033164,0.769518152842520) q[14];
u3(1.31876497518623,-1.40551593883457,-1.15331880667394) q[4];
cx q[4],q[14];
u1(-1.20864797234525) q[14];
u3(0.539108493604088,0.0,0.0) q[4];
cx q[14],q[4];
u3(3.39940193034640,0.0,0.0) q[4];
cx q[4],q[14];
u3(0.881104839713664,1.57173627262925,-2.62993348699729) q[14];
u3(0.126648235137014,-2.07949835628533,-1.55760727400301) q[4];
u3(1.38648660398523,-2.05368497706969,2.13085182353807) q[9];
u3(0.868302832442037,1.02204384472167,-2.69249910015156) q[6];
cx q[6],q[9];
u1(1.57153440817221) q[9];
u3(-2.35459934875357,0.0,0.0) q[6];
cx q[9],q[6];
u3(3.06459977384395,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.503241987433921,-3.51921071797054,1.06714428016251) q[9];
u3(1.22433529704703,-2.98361132864365,1.31465518316370) q[6];
u3(2.27232667986537,4.60220868036963,-1.61504435446139) q[4];
u3(0.410367779091802,0.912553178403365,0.769822184756385) q[12];
cx q[12],q[4];
u1(3.11159751776317) q[4];
u3(-2.36656374041182,0.0,0.0) q[12];
cx q[4],q[12];
u3(0.308151009220435,0.0,0.0) q[12];
cx q[12],q[4];
u3(2.29152177583283,-2.18722838829705,3.12785723583298) q[4];
u3(1.81731763280721,4.02691671088633,1.16155618939874) q[12];
u3(2.04312096476262,1.41141275744645,-3.49387283247313) q[11];
u3(1.42413635414152,-2.27822274562040,2.75132692667239) q[7];
cx q[7],q[11];
u1(1.53963920333214) q[11];
u3(-3.64035624409746,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.84130974071854,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.00888022396317,0.817059613631073,2.04601035973570) q[11];
u3(1.30999095292777,-1.25473467629034,3.96208215533049) q[7];
u3(0.524619525946521,2.52565533230857,-0.994427656503978) q[2];
u3(1.43884971106608,1.31721852753589,-0.409172036351470) q[6];
cx q[6],q[2];
u1(3.64029786013007) q[2];
u3(-3.19724834595000,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.883071070352162,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.728440677910540,2.35781675861919,-3.01501268719282) q[2];
u3(0.310207918104299,-1.93611712058909,2.92078945632446) q[6];
u3(1.49604352953835,1.19017738767524,0.654915903065537) q[14];
u3(0.642981684136308,-0.726207875707293,-1.91842079311666) q[8];
cx q[8],q[14];
u1(3.04709660392987) q[14];
u3(-1.70816111357344,0.0,0.0) q[8];
cx q[14],q[8];
u3(0.775312665447325,0.0,0.0) q[8];
cx q[8],q[14];
u3(1.64766291971373,-0.477815617898170,3.86552788587267) q[14];
u3(0.578701803965963,2.72878888095520,-2.48778245242615) q[8];
u3(1.25074782512557,-0.0281648036405707,-0.182941873734475) q[15];
u3(1.10260425818830,-3.18570791439192,1.24046381407199) q[9];
cx q[9],q[15];
u1(2.47290035298269) q[15];
u3(-1.86639367788885,0.0,0.0) q[9];
cx q[15],q[9];
u3(3.46280172567612,0.0,0.0) q[9];
cx q[9],q[15];
u3(1.34438442446285,-0.143818168874858,1.79955816865394) q[15];
u3(1.21617071387861,0.921117932975414,3.42119019457716) q[9];
u3(1.29622133595372,0.321851380352753,-3.19632568603786) q[0];
u3(1.49388767016846,-0.887541552446921,4.74406196624813) q[13];
cx q[13],q[0];
u1(0.699846274501343) q[0];
u3(-1.43470899679580,0.0,0.0) q[13];
cx q[0],q[13];
u3(3.21208476192947,0.0,0.0) q[13];
cx q[13],q[0];
u3(1.56979856281946,2.36072927527871,-1.05800488930942) q[0];
u3(1.66614670624664,-1.54474054618804,2.84320504699893) q[13];
u3(2.01604026177139,-0.855305023721852,0.472972995410355) q[1];
u3(1.97921882876693,-2.35341001401928,0.633565201807646) q[3];
cx q[3],q[1];
u1(1.16135088980887) q[1];
u3(-3.30014354080219,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.30509924628464,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.22624002205534,0.534197155191449,0.855478409694099) q[1];
u3(0.296936803326452,-1.90163119348815,-0.358810469100098) q[3];
u3(2.69667232343258,-0.354904729015739,0.000690078007708567) q[5];
u3(1.18826655023554,-2.75089635904565,-1.46706653241377) q[10];
cx q[10],q[5];
u1(-0.0417414528740956) q[5];
u3(-1.89537325518516,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.639412390206738,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.610219270058742,2.09070193096799,1.32236211792812) q[5];
u3(1.87474913552618,-0.140565164287792,1.28332906692647) q[10];
u3(0.861561910090108,-1.09877537162828,1.34398729592488) q[14];
u3(0.419454319669924,-0.0788034859123435,-0.847577753109930) q[2];
cx q[2],q[14];
u1(1.44600528109663) q[14];
u3(-2.80706245454600,0.0,0.0) q[2];
cx q[14],q[2];
u3(-0.0436987446812906,0.0,0.0) q[2];
cx q[2],q[14];
u3(1.98480322705910,-2.39926118853466,1.19439930101529) q[14];
u3(2.33231713009922,-5.33095213030471,0.485047346338635) q[2];
u3(1.74768709236424,-0.728951787956653,1.49848014427118) q[11];
u3(1.71475979185351,-0.979617872469935,-1.17534757169488) q[9];
cx q[9],q[11];
u1(1.01180216166891) q[11];
u3(-3.42430932885874,0.0,0.0) q[9];
cx q[11],q[9];
u3(1.66868770524648,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.14575720726300,-1.52496705278049,-0.188037087898390) q[11];
u3(2.05258368846110,-5.39489615811830,0.00444024986713165) q[9];
u3(1.48818718102479,1.78996940834197,-3.11237655423952) q[8];
u3(0.695471329562844,2.20612125938295,-2.84566251558015) q[3];
cx q[3],q[8];
u1(2.89958716821931) q[8];
u3(-2.64807684326721,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.25010602017193,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.55210673900584,-2.06092741680844,1.44478437105420) q[8];
u3(0.937465963902551,2.68212326476165,2.54444069377999) q[3];
u3(1.06886283145616,1.72620447231432,-3.14456642389416) q[7];
u3(0.498295403020934,3.61980737679531,-2.61488545394713) q[1];
cx q[1],q[7];
u1(1.73166293044906) q[7];
u3(0.0248601522732645,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.949231230221295,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.23177542126844,3.58796127700860,-0.742222852301904) q[7];
u3(1.60721181173006,-0.109549746276894,0.0760943469393853) q[1];
u3(1.39954477217580,0.328358399146364,1.26836381552784) q[4];
u3(2.77234431161092,-1.52138276410600,-0.106265721960736) q[13];
cx q[13],q[4];
u1(1.08569375631256) q[4];
u3(-0.388711054575702,0.0,0.0) q[13];
cx q[4],q[13];
u3(2.64865591907338,0.0,0.0) q[13];
cx q[13],q[4];
u3(1.51214345988041,1.54195028686854,-4.53095045643258) q[4];
u3(1.98837419928609,1.64845179996067,3.15857025186774) q[13];
u3(2.56308239331574,1.75544371438245,-1.86072136491643) q[5];
u3(2.22480075645457,1.43265904512147,-3.10081022329375) q[15];
cx q[15],q[5];
u1(2.19297228015298) q[5];
u3(0.112806426341539,0.0,0.0) q[15];
cx q[5],q[15];
u3(1.45649150148713,0.0,0.0) q[15];
cx q[15],q[5];
u3(1.52685420312472,0.479698843830779,-1.52278309950792) q[5];
u3(1.72357135685016,1.29065567607311,-4.31536518087626) q[15];
u3(2.09011637263332,2.68922983526985,-3.08987384985716) q[6];
u3(2.45739495048090,-4.26085954202978,1.95592403614145) q[10];
cx q[10],q[6];
u1(1.06557154478658) q[6];
u3(-0.713291289114392,0.0,0.0) q[10];
cx q[6],q[10];
u3(3.12562121265735,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.91509553840631,0.814310854407769,-0.637580958081602) q[6];
u3(2.46219408391307,-2.82151534929462,-1.45395403175826) q[10];
u3(1.99653802975881,0.411666360698289,-2.83872516582542) q[0];
u3(1.74738519551910,2.82279706288963,-3.33123619313289) q[12];
cx q[12],q[0];
u1(1.45127347233926) q[0];
u3(-3.03552305265959,0.0,0.0) q[12];
cx q[0],q[12];
u3(2.58947118354731,0.0,0.0) q[12];
cx q[12],q[0];
u3(2.14907037242832,-3.35684236768122,0.295078231252316) q[0];
u3(1.94793872660489,-1.16215771787430,-0.108877706412997) q[12];
u3(1.12151451078183,-0.0645797439770601,2.34262607588887) q[12];
u3(1.56713614711561,-1.38171815334973,-1.51610301778233) q[3];
cx q[3],q[12];
u1(1.46779813050050) q[12];
u3(0.136761689860879,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.319747010837091,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.32274657313456,0.688244730187503,-0.411079405235247) q[12];
u3(0.736474680114371,2.20642876283046,-1.02788799679553) q[3];
u3(0.910440727379714,-1.37528165723759,-0.699065161563159) q[5];
u3(1.71296169525112,-4.91147646468463,0.936484721527278) q[4];
cx q[4],q[5];
u1(1.65294116648861) q[5];
u3(-2.38482963391219,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.10562712143311,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.684068983093425,3.43079230190661,-1.73958079072040) q[5];
u3(3.02672376075061,-1.03672647845394,5.16426071052368) q[4];
u3(2.14614276096369,3.53926595532472,-0.641387549951381) q[7];
u3(2.69359093387442,-0.736877659214794,-5.51876940487367) q[8];
cx q[8],q[7];
u1(-0.815590448194226) q[7];
u3(-1.61021610422738,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.894363214208981,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.77773020487205,-1.24311571329475,-0.604569649358720) q[7];
u3(1.34130056551646,0.557869946263205,5.37842403219932) q[8];
u3(1.41023624586572,2.46301741432854,-3.53738024967735) q[6];
u3(0.871265672098015,2.79210496925595,-2.31731380730028) q[13];
cx q[13],q[6];
u1(1.79827427445763) q[6];
u3(-0.497982741259428,0.0,0.0) q[13];
cx q[6],q[13];
u3(0.0599336484681698,0.0,0.0) q[13];
cx q[13],q[6];
u3(0.969038588436943,-1.98376736418531,1.61189008355359) q[6];
u3(2.33169419650645,1.49147508567251,-4.57927613730726) q[13];
u3(2.22798775199447,0.138478734838704,2.06948564281066) q[11];
u3(2.12246301826790,-1.09106000860249,-1.81707901592858) q[0];
cx q[0],q[11];
u1(1.40103026481786) q[11];
u3(-2.94955039429135,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.443280680109836,0.0,0.0) q[0];
cx q[0],q[11];
u3(1.20683594145615,0.849211832146724,1.61560348279430) q[11];
u3(1.63229633438064,0.103346108582529,-3.45066467344037) q[0];
u3(2.15227271040520,-2.03216609171552,3.56956285711520) q[14];
u3(0.964084672148582,-1.32383586012626,2.59931918823858) q[15];
cx q[15],q[14];
u1(1.44431089989386) q[14];
u3(-0.836278142495644,0.0,0.0) q[15];
cx q[14],q[15];
u3(2.76020590157180,0.0,0.0) q[15];
cx q[15],q[14];
u3(0.839846743950301,2.41819390102445,-0.581885358342013) q[14];
u3(1.60395524169648,1.23361671957435,1.00300135901100) q[15];
u3(1.89078935100010,2.06192508206996,-3.99826261720611) q[9];
u3(2.04984121251986,3.28795398126897,-2.88491032227177) q[1];
cx q[1],q[9];
u1(0.691993743897889) q[9];
u3(-1.74188035991215,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.68238208599250,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.21415273692577,3.12931770162880,-1.46602601751180) q[9];
u3(2.18635050800358,1.64820843579108,4.61239434722458) q[1];
u3(0.892759700771481,0.685394140171024,-1.34485987558548) q[2];
u3(0.754273421349240,1.30693877254892,-4.45638058311884) q[10];
cx q[10],q[2];
u1(0.554952007029828) q[2];
u3(-1.34068895823950,0.0,0.0) q[10];
cx q[2],q[10];
u3(2.48251061087621,0.0,0.0) q[10];
cx q[10],q[2];
u3(0.342929515685581,-2.29821082852851,2.79988934939840) q[2];
u3(1.83402692756583,-0.680699445108536,4.33777797294208) q[10];
u3(2.37599611054389,-2.19998775808838,3.45511961796255) q[0];
u3(1.02919825448479,-1.46851314221950,2.39021381636430) q[10];
cx q[10],q[0];
u1(4.26566377581125) q[0];
u3(-3.17259544426072,0.0,0.0) q[10];
cx q[0],q[10];
u3(-0.398001490111930,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.33797940592521,-0.666306455427504,3.19031983932212) q[0];
u3(2.70994904429932,-3.07038033625755,0.990380104181775) q[10];
u3(1.94411982587102,1.55300053888264,-3.39592666685480) q[15];
u3(1.18460826442965,-2.46471754795063,3.18371309788230) q[3];
cx q[3],q[15];
u1(1.38875910665821) q[15];
u3(-0.0621534442048113,0.0,0.0) q[3];
cx q[15],q[3];
u3(1.48630604906775,0.0,0.0) q[3];
cx q[3],q[15];
u3(2.06712537761123,-2.80765435787797,3.15944204004993) q[15];
u3(0.348993573926005,1.42385388674930,-4.26682466713863) q[3];
u3(2.18753860399851,-0.806757838643319,1.22268906374655) q[7];
u3(2.17924949122678,-1.95132138213778,0.00144812353587331) q[12];
cx q[12],q[7];
u1(1.96893009546148) q[7];
u3(-2.42022922048962,0.0,0.0) q[12];
cx q[7],q[12];
u3(-0.0568230873444275,0.0,0.0) q[12];
cx q[12],q[7];
u3(1.47913126272986,3.76248674943028,-1.79266028855538) q[7];
u3(0.700767037300231,2.59589862683366,-3.38786923762307) q[12];
u3(0.592421809902137,1.67484394254389,-1.19416077069410) q[9];
u3(0.602345061846953,0.315376737924186,-3.01104356408304) q[1];
cx q[1],q[9];
u1(2.25725875346631) q[9];
u3(-1.74171179772915,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.63414533596456,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.80340129297824,2.39173820097160,1.00071964567535) q[9];
u3(1.01273869246933,1.28314335350501,-1.38690054078836) q[1];
u3(1.09319638398121,1.60872254350627,-1.21771278006270) q[6];
u3(0.443165116451779,1.83281028558696,-4.02722648246613) q[2];
cx q[2],q[6];
u1(-0.415759835472601) q[6];
u3(-1.89918473787553,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.809637175023691,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.55047531947413,-3.69340444103145,-0.293211946551522) q[6];
u3(2.38484929053117,-2.24519652972274,0.312742431341282) q[2];
u3(2.34510728688473,2.81717427195159,-2.93552250847311) q[13];
u3(1.88229662133059,1.65214641306310,-1.40100524348816) q[8];
cx q[8],q[13];
u1(1.40717670866772) q[13];
u3(-0.420966611816642,0.0,0.0) q[8];
cx q[13],q[8];
u3(3.20040057231325,0.0,0.0) q[8];
cx q[8],q[13];
u3(1.71877445873987,2.38448978183121,-2.47960306688967) q[13];
u3(1.52590282990294,1.17097068165971,3.29659615727484) q[8];
u3(1.84740975347187,1.50025665582460,-3.10907208149129) q[14];
u3(1.99445881012469,1.65476782283613,-3.37230268722159) q[5];
cx q[5],q[14];
u1(0.574204004969367) q[14];
u3(-1.01009976046280,0.0,0.0) q[5];
cx q[14],q[5];
u3(1.78430772446198,0.0,0.0) q[5];
cx q[5],q[14];
u3(2.28373898307175,0.815898506406051,-2.38146425494619) q[14];
u3(1.80356422537544,-2.26729895686093,-0.837677851339216) q[5];
u3(0.469298123139362,2.17211904959020,-2.25098235298238) q[4];
u3(0.735502834827367,0.251312674444368,-1.53967141856522) q[11];
cx q[11],q[4];
u1(-0.0384604363489658) q[4];
u3(1.00376885949695,0.0,0.0) q[11];
cx q[4],q[11];
u3(3.69821664681358,0.0,0.0) q[11];
cx q[11],q[4];
u3(0.920392338185737,1.76247806153583,-0.930784921154798) q[4];
u3(1.14743701333713,2.47889828679474,2.69387489942302) q[11];
u3(1.48327497696471,0.826417877916408,-3.69726922447976) q[14];
u3(2.24399137777467,2.82529827723304,-2.94722622859760) q[10];
cx q[10],q[14];
u1(1.66926850707948) q[14];
u3(-3.01875791609795,0.0,0.0) q[10];
cx q[14],q[10];
u3(2.28477824415374,0.0,0.0) q[10];
cx q[10],q[14];
u3(1.64213202206192,1.80063057275173,-4.35860683167785) q[14];
u3(1.78231569517471,2.04220685768615,-0.626505721197428) q[10];
u3(0.373313092509292,-2.52877800597269,1.39439438528651) q[12];
u3(0.551171209247311,-2.22197542080264,0.460473863390943) q[3];
cx q[3],q[12];
u1(3.93461778710603) q[12];
u3(-4.52489195089274,0.0,0.0) q[3];
cx q[12],q[3];
u3(-0.902851851755513,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.37842572675390,0.785701942820774,-0.772928017430022) q[12];
u3(0.614153655694867,-5.41034201735113,-0.133537921736044) q[3];
u3(0.342825813504102,0.287881330505310,2.79400690843635) q[9];
u3(1.43846208308372,3.01196758659572,3.22923689575356) q[15];
cx q[15],q[9];
u1(3.54204702410571) q[9];
u3(-1.31083224508502,0.0,0.0) q[15];
cx q[9],q[15];
u3(1.98246311335514,0.0,0.0) q[15];
cx q[15],q[9];
u3(0.800271218623828,2.73331882441437,-1.39766225983635) q[9];
u3(0.994231082636466,-0.559637678134613,-1.77899814895256) q[15];
u3(1.64188956300613,0.659317230340325,0.811376253912382) q[7];
u3(1.35202265098232,-1.35920143828638,-1.72699970826195) q[6];
cx q[6],q[7];
u1(2.38889071743424) q[7];
u3(-2.18741500348487,0.0,0.0) q[6];
cx q[7],q[6];
u3(-0.281160542656428,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.761378376193263,0.892095825272436,1.14507239354865) q[7];
u3(2.26619308563617,-1.71105211006248,2.65239190534633) q[6];
u3(1.38023393703173,1.60925196019027,-0.193180190037179) q[5];
u3(2.68977416235771,0.835712766382755,-1.89508334166897) q[2];
cx q[2],q[5];
u1(-0.148253470527465) q[5];
u3(0.330977572941368,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.18532446118511,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.65293349042728,2.33099113496527,-2.97012039116125) q[5];
u3(2.55617279249846,2.57053357563552,3.60701206600054) q[2];
u3(0.785982502098445,1.18743602479499,1.25252939694714) q[13];
u3(2.83678564974445,-0.671910779471283,-2.11727997280646) q[11];
cx q[11],q[13];
u1(3.41157242043439) q[13];
u3(-4.20532916234722,0.0,0.0) q[11];
cx q[13],q[11];
u3(-0.530978261469328,0.0,0.0) q[11];
cx q[11],q[13];
u3(2.54791602878173,2.20934685982624,-2.46768927727881) q[13];
u3(2.08359804296744,-0.291553021010512,-5.41748169847217) q[11];
u3(1.53126981636790,1.29664307819609,-0.682758995488924) q[4];
u3(2.21301042138335,-4.43575097731715,0.693570169363929) q[1];
cx q[1],q[4];
u1(1.70643071018549) q[4];
u3(-3.24477790538295,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.334222311305503,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.29432042782529,3.11853100292372,-2.27147637783456) q[4];
u3(1.27013508205759,-0.931761332027884,-0.936907678142818) q[1];
u3(0.691860212143059,1.59635170735887,-3.17980724034012) q[0];
u3(1.44608065074379,2.73708567979710,-3.26566961572829) q[8];
cx q[8],q[0];
u1(1.09756507290924) q[0];
u3(-1.46568580718628,0.0,0.0) q[8];
cx q[0],q[8];
u3(3.53839510446777,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.26450386900752,-1.73111263616935,1.81229939943772) q[0];
u3(1.63054285724877,2.34004854261865,3.29191230591147) q[8];
u3(1.37593862929957,1.40518137778148,-2.80009014311286) q[11];
u3(2.31367961082331,-1.61582849468817,3.35502576943272) q[14];
cx q[14],q[11];
u1(3.58577980857303) q[11];
u3(-1.18322557705357,0.0,0.0) q[14];
cx q[11],q[14];
u3(2.24035585091653,0.0,0.0) q[14];
cx q[14],q[11];
u3(0.133036410286938,3.16738965293312,-2.53228183687752) q[11];
u3(2.98041814222143,-1.36188183977154,-1.50337469239737) q[14];
u3(1.13028367985544,2.54426634901751,-3.63528547239054) q[15];
u3(2.02096999104321,2.92997319890957,-2.40300779393850) q[8];
cx q[8],q[15];
u1(-0.121220987014162) q[15];
u3(0.524337122083856,0.0,0.0) q[8];
cx q[15],q[8];
u3(4.00458279136099,0.0,0.0) q[8];
cx q[8],q[15];
u3(2.24322802500004,-1.01179211565464,1.98126934976446) q[15];
u3(1.70052523166836,-0.234342729680016,-4.00696751381580) q[8];
u3(2.03167498051766,-0.302713805567445,2.53329157029452) q[13];
u3(1.84926073377782,-1.31958652114628,-1.05554232451439) q[5];
cx q[5],q[13];
u1(1.45570931882382) q[13];
u3(-0.880425450319648,0.0,0.0) q[5];
cx q[13],q[5];
u3(-0.372765110766525,0.0,0.0) q[5];
cx q[5],q[13];
u3(2.51435539078387,-2.15374000187125,1.54304909839402) q[13];
u3(1.75466606190431,2.73365041845202,3.17288128971953) q[5];
u3(1.16308144550719,-1.65070107612210,-1.31491123702421) q[6];
u3(2.01167267731035,-2.25678092623150,0.0656650009019464) q[3];
cx q[3],q[6];
u1(0.746377023065708) q[6];
u3(-1.27831326090469,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.84503397554421,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.637316784053095,1.53476326943670,-0.906837071431845) q[6];
u3(1.14386132976119,-3.07273031933700,-0.303651133607426) q[3];
u3(2.02109602561849,-2.35336428351613,-0.187887583488020) q[0];
u3(2.19004911011304,-3.11005982746616,-1.72073881949544) q[2];
cx q[2],q[0];
u1(1.71879327527380) q[0];
u3(0.320487760816350,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.06824718503349,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.46712891606310,-4.32967831015464,0.851874856709125) q[0];
u3(0.164718043663762,4.24933563971437,-0.716853332949238) q[2];
u3(1.32312669783644,-2.45200959067461,3.48015165109449) q[4];
u3(2.20702424471286,1.77096408426865,-1.64703659590572) q[7];
cx q[7],q[4];
u1(0.252420825483191) q[4];
u3(-1.24747753619849,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.17142835702186,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.66583585643803,3.36318761639618,-2.17145956183454) q[4];
u3(1.19014785324123,4.40005964081140,1.85711341589040) q[7];
u3(2.74882650076526,-2.01891156118194,3.54249289998781) q[1];
u3(1.21336439117477,2.08717846319288,0.291097146660468) q[12];
cx q[12],q[1];
u1(2.09963506147480) q[1];
u3(-2.89889920663633,0.0,0.0) q[12];
cx q[1],q[12];
u3(0.740817944060782,0.0,0.0) q[12];
cx q[12],q[1];
u3(0.568572457220606,1.60463324205359,0.541732368977454) q[1];
u3(1.87423206311684,-0.693150141728863,4.21335714717378) q[12];
u3(1.99371206610826,0.995394684950960,0.911044921127980) q[10];
u3(2.14534905988438,0.574068702602453,-2.50344592402216) q[9];
cx q[9],q[10];
u1(1.80204522727470) q[10];
u3(-2.44417235296005,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.282561938821859,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.54090747139650,-1.04734458831492,0.133461683885191) q[10];
u3(1.75195886230724,1.57921028827569,3.77703601746174) q[9];
u3(1.52580420892619,0.810860031320506,-3.00259356678572) q[2];
u3(1.19321534510675,-3.53548234739937,1.96980365885659) q[7];
cx q[7],q[2];
u1(1.01784606608051) q[2];
u3(0.0717236198757332,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.01197099116211,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.23834769729821,-4.13819383162497,1.37537246653362) q[2];
u3(0.778223200903737,-1.26023089165420,-3.06861298214057) q[7];
u3(0.748442717663183,-0.281567951580710,0.0903013492995124) q[12];
u3(1.19490030674596,-2.98998523643958,1.32162802415492) q[13];
cx q[13],q[12];
u1(1.54965896864990) q[12];
u3(0.361886590011672,0.0,0.0) q[13];
cx q[12],q[13];
u3(0.632207703498966,0.0,0.0) q[13];
cx q[13],q[12];
u3(0.542159027065526,0.513539024368793,0.448440949354995) q[12];
u3(1.53887022194383,-0.187742049483948,-5.48528405780349) q[13];
u3(2.34403794199055,0.841941977225817,-3.73022946179614) q[0];
u3(1.01505858896362,3.04163500731579,-2.97856171030318) q[8];
cx q[8],q[0];
u1(1.59513286159707) q[0];
u3(0.00237357897557700,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.398697744415864,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.862418088887214,1.63909556279066,-0.928746012160283) q[0];
u3(0.701386534746489,3.42766696402300,0.477190503049338) q[8];
u3(2.19047104922250,-0.726924423076212,1.33692096096087) q[9];
u3(1.67363149618260,-2.01410728190356,-2.31951861450853) q[1];
cx q[1],q[9];
u1(1.95958715916281) q[9];
u3(0.465541044843167,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.911420850427896,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.814890030683226,-1.69752097885942,3.38757726937660) q[9];
u3(2.91002617811549,-3.95180575616683,0.501288826612512) q[1];
u3(1.98190534772165,-1.84901909970299,-0.786392534296950) q[3];
u3(1.33478402335048,-3.86495642846186,0.211329414362126) q[14];
cx q[14],q[3];
u1(2.81648250592814) q[3];
u3(-1.62779156723802,0.0,0.0) q[14];
cx q[3],q[14];
u3(0.686637841454061,0.0,0.0) q[14];
cx q[14],q[3];
u3(1.54643825995497,1.00834311488392,0.450264865236819) q[3];
u3(1.65407890575260,1.76533276071917,0.219454072414963) q[14];
u3(1.80564286554193,0.306039527135413,-2.93481936321454) q[4];
u3(1.91472810237971,3.23257099684827,-2.33761365878108) q[10];
cx q[10],q[4];
u1(-0.446763730596824) q[4];
u3(1.29635410935897,0.0,0.0) q[10];
cx q[4],q[10];
u3(3.50967979910567,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.308291921432708,-0.110819538928659,2.66323397720445) q[4];
u3(2.71611202136980,5.49046776778209,0.371648068856508) q[10];
u3(1.68584967834113,2.24151876672376,0.319186346319551) q[11];
u3(2.36555807301887,0.771170757926098,-3.13456699777552) q[15];
cx q[15],q[11];
u1(1.47216487270181) q[11];
u3(-3.16196918130958,0.0,0.0) q[15];
cx q[11],q[15];
u3(2.93102937545707,0.0,0.0) q[15];
cx q[15],q[11];
u3(2.22439064080317,3.70073455270931,-1.31304196404373) q[11];
u3(0.191028281059131,-0.344026650175182,5.44035586839434) q[15];
u3(2.73130357539693,1.16337584159942,-1.75277412447235) q[5];
u3(1.80013068738242,1.29237611853299,-4.59228153738658) q[6];
cx q[6],q[5];
u1(0.940094979481819) q[5];
u3(-0.0326972572820730,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.77138723424973,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.292333795089897,-0.0712614868440151,-4.18712094566067) q[5];
u3(1.89171687951569,5.17912833748095,0.741964004110537) q[6];
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