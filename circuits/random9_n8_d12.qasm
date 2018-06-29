OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(0.505247423771709,-1.15968837222959,-0.0166296227768253) q[3];
u3(1.38966996646412,-4.30012262390232,1.71119780339958) q[0];
cx q[0],q[3];
u1(2.17471269213222) q[3];
u3(-2.75565426779251,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.866378433298655,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.75614041600987,-2.08840118613702,-1.19511890340141) q[3];
u3(2.40439781606727,-1.21674449449822,-2.64671723940262) q[0];
u3(1.34019738608319,0.609615510868442,0.943181006272543) q[2];
u3(1.77371993073650,-0.520198848033387,-2.27824636086882) q[5];
cx q[5],q[2];
u1(2.56686847078167) q[2];
u3(-2.19346828719368,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.88643259693749,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.20524143326645,1.85475283745839,-1.81361069224151) q[2];
u3(0.730007892167093,0.0542372901727666,1.35476402614934) q[5];
u3(0.710791136945060,0.927646490623265,-1.68799067036841) q[1];
u3(0.858907531001732,-2.83576393175752,1.83582487675537) q[6];
cx q[6],q[1];
u1(1.52351970186282) q[1];
u3(-3.36180991546266,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.539006024657299,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.949784937363018,2.09769117510716,0.110758686090853) q[1];
u3(2.25015164481557,-1.59119285931823,3.20327215459758) q[6];
u3(1.04639045790310,3.80270728985646,-1.47170672138870) q[7];
u3(0.892517868237820,0.680500372754571,-0.631625666502725) q[4];
cx q[4],q[7];
u1(0.121877362112442) q[7];
u3(-1.26185657702609,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.82979554911761,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.14677155875361,3.03327272295932,-2.18694511610593) q[7];
u3(2.21573942646156,0.706481779402419,-0.789271925913535) q[4];
u3(1.16109320103723,3.30590028240423,-1.04587639950307) q[3];
u3(1.42743678743063,2.03722269027915,-1.07116589788620) q[7];
cx q[7],q[3];
u1(0.0215686598727971) q[3];
u3(-2.40424823100140,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.07834835794974,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.30336922361865,2.53173299715341,1.44807967684491) q[3];
u3(1.18084077048439,4.23740147185190,-0.178183772677882) q[7];
u3(1.70638714657018,1.37700337048614,-3.12936365283323) q[5];
u3(1.47586581708435,-2.68663075369458,2.90918464424051) q[4];
cx q[4],q[5];
u1(3.15604530060534) q[5];
u3(-1.57051066041751,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.63203859712805,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.52755798936537,-0.0199539867395788,-1.75702328499208) q[5];
u3(1.89184940248716,-0.850091865640289,-1.66094187011730) q[4];
u3(2.10321124790882,-0.345193275804332,1.84936531302303) q[2];
u3(1.86090016666036,-2.02300229383730,-0.594617331722060) q[1];
cx q[1],q[2];
u1(1.78902538800912) q[2];
u3(-2.42021466531021,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.123666387054986,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.43815060131353,-3.04017178298575,0.558214707918695) q[2];
u3(1.33000147595470,1.41527818126343,4.24211224554736) q[1];
u3(1.58261798945071,0.0902354764449383,1.99684292816335) q[6];
u3(2.00713031715126,-2.08115220137001,-1.18995144021522) q[0];
cx q[0],q[6];
u1(2.38429075286925) q[6];
u3(0.260922245108103,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.31646160251116,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.31057782144745,-0.624518837220515,-1.06818989797211) q[6];
u3(0.846571659723269,1.87945615597488,0.125502993661139) q[0];
u3(0.705092334343501,0.562506514247249,0.372327757881330) q[4];
u3(1.33066326305626,-0.683041449380616,-2.69678415709522) q[2];
cx q[2],q[4];
u1(1.98347726923497) q[4];
u3(-2.55251382332797,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.0916381259590873,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.07984583229680,2.83343624194430,-2.18849292473755) q[4];
u3(2.30795984426185,-0.0556128609312919,1.37757945016058) q[2];
u3(1.72544804957933,0.0110078478153528,1.48302809526720) q[6];
u3(1.86766461799382,-1.02561942249588,-2.38057231868964) q[1];
cx q[1],q[6];
u1(3.26163186474865) q[6];
u3(-1.09093054687327,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.81478705308576,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.33961894919832,-0.203772443264178,3.83135480398440) q[6];
u3(1.16623488897137,0.567889824791248,-3.46418946739760) q[1];
u3(2.07959499000392,-2.62227625569324,2.30830269588532) q[5];
u3(2.37664501294460,-0.599720889071157,0.729605535210958) q[7];
cx q[7],q[5];
u1(2.10482849047983) q[5];
u3(-1.59206900299394,0.0,0.0) q[7];
cx q[5],q[7];
u3(3.91294790981330,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.79489739303937,1.54049589635671,-4.25393734468150) q[5];
u3(1.71159066251326,-1.09048942354935,-4.78123518522752) q[7];
u3(1.35855273177089,-0.425404957096353,-0.984733521287497) q[3];
u3(2.26081106537587,0.517741982716861,-5.39578489295163) q[0];
cx q[0],q[3];
u1(1.54469915736910) q[3];
u3(-0.604133910685828,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.14599276895018,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.08528389536765,-3.57450799252725,2.06380785122028) q[3];
u3(1.20993704726799,-3.48696681667682,2.62429930077672) q[0];
u3(1.19639709168630,-1.65943405426449,0.00483907482378398) q[1];
u3(0.671308044579109,-1.96695769546431,0.641069860463076) q[7];
cx q[7],q[1];
u1(0.220268443460101) q[1];
u3(-0.692129502626976,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.34650307301665,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.53461742540946,2.41292525751960,-1.65188593020607) q[1];
u3(0.758750610765331,-2.06372663875744,3.73788594031225) q[7];
u3(2.39158839396305,-2.14923269786657,0.00133516698223413) q[5];
u3(1.59715977240325,-3.89004189393288,1.09924429034985) q[0];
cx q[0],q[5];
u1(2.08016722475239) q[5];
u3(-2.41002561159281,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.438596049444916,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.893341425219863,-1.53197776734379,0.436884815141649) q[5];
u3(0.245909251352874,-2.01658067651126,-2.80281369173902) q[0];
u3(2.17958457985922,1.35014358905489,0.360488956859121) q[4];
u3(2.04777740134454,-1.06497431091958,-3.33600535576171) q[3];
cx q[3],q[4];
u1(4.22111682645297) q[4];
u3(-3.82005796272390,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.394317199297382,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.97405722982903,-3.20345732257184,2.68453938785877) q[4];
u3(1.23376181322824,3.26587265643506,1.53423110921418) q[3];
u3(1.38765367448994,2.34766343561966,-3.65058189671278) q[6];
u3(1.88367116955611,-2.41144865080982,3.68436951639074) q[2];
cx q[2],q[6];
u1(3.33001741580570) q[6];
u3(-0.847703712508684,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.70317844178991,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.06033687680500,-1.49552861789358,0.961538000630934) q[6];
u3(1.03162463292421,-0.680319006788606,-0.388084262110653) q[2];
u3(0.776898127078406,1.64863663878416,-1.77765603036287) q[0];
u3(0.354590315797135,1.27052364719875,-1.55114441012376) q[3];
cx q[3],q[0];
u1(0.142231380554823) q[0];
u3(-1.81975346693500,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.37371198346602,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.644040923143365,2.94554850473954,-1.80943621148195) q[0];
u3(2.16783528647981,-5.39634920477501,0.840092626321806) q[3];
u3(1.53029803027241,-0.331536501676580,1.84755702413534) q[1];
u3(1.49459685831942,-1.29759382679216,-2.08175669937508) q[6];
cx q[6],q[1];
u1(4.37995998685040) q[1];
u3(-3.29349109061768,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.475431393042688,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.56552644211359,-1.22478656449393,0.227782392466050) q[1];
u3(1.32092308449193,3.27915788018877,-0.280421549143635) q[6];
u3(0.917083219466465,-0.436203729138739,0.780024412119289) q[4];
u3(0.377651439363570,-1.98240556703466,1.57588351457425) q[7];
cx q[7],q[4];
u1(-0.318055357882101) q[4];
u3(-2.27783479899239,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.32288763475102,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.46858217505228,3.50217046108417,-0.860322465164837) q[4];
u3(2.10299695698028,-4.82014287969258,1.38263671973787) q[7];
u3(2.30483859922596,-0.382508725201205,1.60859644704833) q[2];
u3(2.10578587242002,-2.44606208149062,-0.585291597396646) q[5];
cx q[5],q[2];
u1(0.165547805352812) q[2];
u3(-1.40609713641545,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.23212999715814,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.46729647083232,2.51766696365047,-0.526687374501340) q[2];
u3(1.68212349841918,2.35386958571963,1.96182585892408) q[5];
u3(1.32738018764631,0.407838648719445,-2.41016049466139) q[6];
u3(0.389013268695396,-3.39678515792718,2.05008406208770) q[0];
cx q[0],q[6];
u1(2.61146358014873) q[6];
u3(-1.76905758198100,0.0,0.0) q[0];
cx q[6],q[0];
u3(-0.105511835069891,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.60631114337474,-2.34177480610965,1.88550664211525) q[6];
u3(1.03397023482790,-1.36323897914792,3.18630027411623) q[0];
u3(1.69244529685161,-1.60442307113381,-0.115474974669483) q[2];
u3(1.43838362927230,-2.82625097040904,-0.713383184568117) q[5];
cx q[5],q[2];
u1(1.00404569687850) q[2];
u3(-3.41431894712630,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.08035211370438,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.320337384498996,0.944238445888378,2.47025794943548) q[2];
u3(1.99010692106599,0.770684763074423,-1.79042399076662) q[5];
u3(2.66010153959228,1.80136006877487,0.806284887767387) q[7];
u3(1.98066393949293,-0.101576982981174,-3.08093935912059) q[4];
cx q[4],q[7];
u1(1.97299056928134) q[7];
u3(-1.70238049426068,0.0,0.0) q[4];
cx q[7],q[4];
u3(3.27529971349508,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.666252811047224,3.60652861780483,-2.25190932091557) q[7];
u3(0.574106292744875,-1.63403290775407,1.57669681351683) q[4];
u3(1.95058438202746,-1.00288531628867,-0.745873820046038) q[1];
u3(0.482399356248496,0.276364549061088,-5.53666144514669) q[3];
cx q[3],q[1];
u1(1.07163990547964) q[1];
u3(-0.129002641928224,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.42685262311195,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.00718859846388,4.40411658200399,-1.60961479574306) q[1];
u3(0.214033273318308,0.237767729973114,3.28490915925862) q[3];
u3(1.73447342071084,-1.39095187287571,-0.500141748353379) q[0];
u3(1.84506251768414,-2.56995568200536,0.0108579938213025) q[5];
cx q[5],q[0];
u1(1.98692303590914) q[0];
u3(-2.75900358973569,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.738073309463650,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.79662867043277,2.73804730513185,-1.14379201745597) q[0];
u3(0.263422694436479,3.95115352127044,-2.05305860815573) q[5];
u3(0.667679824577110,1.50063789709454,-2.46921994671268) q[3];
u3(1.57799506101581,2.43745568276602,-3.46707276321187) q[6];
cx q[6],q[3];
u1(1.29544083885093) q[3];
u3(-0.554576536576481,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.48320311118585,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.50348769688854,0.788609431068701,1.67988056353017) q[3];
u3(0.950194045840308,1.02964434674279,-3.60400809093430) q[6];
u3(1.93620410078115,-3.83846323090260,1.38450815379615) q[4];
u3(1.56694145382578,0.166702739596118,2.92606971806066) q[1];
cx q[1],q[4];
u1(0.00764537039004676) q[4];
u3(-1.63306798263892,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.41260898329878,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.64494933137351,0.415683420032581,2.59022435139959) q[4];
u3(1.73442308491034,1.53256181723544,-3.18105542126092) q[1];
u3(1.76058850973047,0.731335141943279,-3.06046262567269) q[7];
u3(2.62775610697036,3.84387387720836,-2.13150109262434) q[2];
cx q[2],q[7];
u1(-0.0481725800194777) q[7];
u3(-1.67889735589905,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.772266556284300,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.945352583618469,-2.29043834855809,3.18478939468937) q[7];
u3(0.988312862412270,-0.204598630769294,2.39067161408985) q[2];
u3(1.45261562313336,0.327919845289110,1.93346281051068) q[3];
u3(1.63582209983474,-0.835534929187609,-2.72213608461117) q[2];
cx q[2],q[3];
u1(1.29370414662701) q[3];
u3(-0.769915112395503,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.88138445764787,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.22718141805839,2.36270941093804,0.0492412819013945) q[3];
u3(0.760037136559883,0.386847654815879,-5.00497716594348) q[2];
u3(2.00700202169151,-2.59885598781032,1.22485631183109) q[6];
u3(2.61677604120028,-4.05589732743236,-1.99701311354549) q[7];
cx q[7],q[6];
u1(-1.09044216971710) q[6];
u3(0.518078728383633,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.86304738092791,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.902979946162249,-3.39205671697211,0.747613511929685) q[6];
u3(0.911695561906027,-2.31708845691070,1.07717210341899) q[7];
u3(1.92682231967525,-2.48045621080932,0.0953965846952809) q[4];
u3(1.93914551767900,-2.88072031511316,0.603833033768730) q[1];
cx q[1],q[4];
u1(2.97644760545810) q[4];
u3(-1.71312272966353,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.573280942113592,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.836322372104430,-2.42319846973701,0.00431982146078447) q[4];
u3(2.41679104920174,-0.936506293182934,-1.27950062321855) q[1];
u3(0.343057946962979,-1.51981979446758,1.96803302842834) q[5];
u3(1.28861240775811,-2.18557984529727,1.42887042489477) q[0];
cx q[0],q[5];
u1(3.28781087059895) q[5];
u3(-0.641984188593509,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.00837418822914,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.61880231778835,-2.28134998164562,-0.566823768342949) q[5];
u3(2.27315843741578,2.41469561110035,-0.183507629722778) q[0];
u3(1.32325930945466,1.05104407810595,1.53004339145735) q[6];
u3(1.09579265430101,-1.14565410990406,-1.96560417846982) q[5];
cx q[5],q[6];
u1(1.95204068644028) q[6];
u3(-2.41310410005724,0.0,0.0) q[5];
cx q[6],q[5];
u3(-0.102452078648932,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.01655949732374,-3.25502804067216,2.52908831149147) q[6];
u3(1.83076218066486,-0.778446552920811,2.11528773088728) q[5];
u3(0.579453658130203,-1.07684576462510,0.764942185011780) q[0];
u3(0.430434601640877,-2.21727591856610,-0.110944848757717) q[1];
cx q[1],q[0];
u1(1.51570852100932) q[0];
u3(-3.32147099791632,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.61833261958027,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.806442945039366,-0.387888911115807,-1.89625207257739) q[0];
u3(2.19466482325449,1.18251280586724,-1.22292423618542) q[1];
u3(2.07672473670673,1.50484715305277,1.08816150365777) q[7];
u3(0.663624594126335,-5.16510340671791,0.0826662680748953) q[4];
cx q[4],q[7];
u1(1.66075790993273) q[7];
u3(-2.63835578665981,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.07457268052436,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.02417847586404,1.98628583385037,-4.07689568203375) q[7];
u3(1.71708528898910,1.68313028419901,-3.78718601891850) q[4];
u3(2.47568910979150,-1.22906619946423,0.239724910869776) q[3];
u3(2.35783433610053,-2.35334731035494,-0.763401396113521) q[2];
cx q[2],q[3];
u1(1.62746269095023) q[3];
u3(-0.0896015306811178,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.21346907601939,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.90518614801816,0.747838549212934,-2.78153435109358) q[3];
u3(1.79857177739154,-2.32292452559248,-3.06401579338463) q[2];
u3(0.981536942015493,0.191736804142372,1.65578139822703) q[6];
u3(1.51033086555275,-0.248810719233987,-0.683115742387791) q[3];
cx q[3],q[6];
u1(0.915782191325076) q[6];
u3(-3.11212719375100,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.85885177869155,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.37844219097234,-1.77047591171040,1.99045678786280) q[6];
u3(1.40325468997949,-4.03425524892201,-0.236754297168617) q[3];
u3(2.52326765430443,-1.38850562995769,-1.45877909847025) q[4];
u3(1.26140310704749,-4.61950049437552,1.08027355161253) q[5];
cx q[5],q[4];
u1(1.08500295004726) q[4];
u3(-0.160217360717550,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.68899915713043,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.632749947046000,-1.56820107582741,2.33430697333372) q[4];
u3(2.94102482414720,-4.76452307210237,0.447278243863403) q[5];
u3(0.525053797801359,0.0552246930619187,0.840615365711535) q[0];
u3(0.159592223744032,0.290094674660019,-1.26805957744600) q[1];
cx q[1],q[0];
u1(1.36065799174519) q[0];
u3(0.120646230736750,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.23836977647690,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.68366903318413,3.97411615602393,-1.96275467952739) q[0];
u3(2.84958363874908,-0.249412663662090,4.45511659673541) q[1];
u3(0.382032589561286,2.08652924966182,-2.96417022985913) q[7];
u3(1.16358778851465,0.795103720062827,-1.67091142509060) q[2];
cx q[2],q[7];
u1(1.54143230384635) q[7];
u3(-2.59788237735929,0.0,0.0) q[2];
cx q[7],q[2];
u3(3.28960676609227,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.28354003510294,0.526979503128334,-0.932988375633564) q[7];
u3(2.24134125153518,0.767203471441511,4.23616275895562) q[2];
u3(1.95356600130085,-0.191629214707586,-1.07175547238912) q[1];
u3(1.67941638029353,-3.39336823460107,0.976452402352279) q[7];
cx q[7],q[1];
u1(1.68326672112884) q[1];
u3(-2.99933630610399,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.712283391474100,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.50319772860440,-1.06007903620054,2.38591341603791) q[1];
u3(2.13965441054236,0.568138771659147,-4.63937126691838) q[7];
u3(2.16526440535228,-0.248296286655046,-1.19387971478168) q[2];
u3(0.820585437637124,-4.77051536819680,1.27557624060070) q[0];
cx q[0],q[2];
u1(1.08196702738965) q[2];
u3(-0.204294625499204,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.42485526874680,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.74760758996403,-0.0994221875774734,-1.05075888518311) q[2];
u3(0.962539779906334,-5.30573802322733,0.161314089056246) q[0];
u3(1.29316770687667,1.51230222769885,0.787650612474630) q[6];
u3(1.81603088301927,1.10775078159492,-3.30092594390043) q[3];
cx q[3],q[6];
u1(1.43769479853649) q[6];
u3(-0.484068530677395,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.15220114256761,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.58339575417958,2.34906332059797,-1.97428115259898) q[6];
u3(2.29831145274132,0.493565920142245,-1.72867038507853) q[3];
u3(1.26315603262739,0.263215641357182,-2.13143083719272) q[5];
u3(2.26549018862308,-3.37865979515755,2.74876632549649) q[4];
cx q[4],q[5];
u1(0.937037754403982) q[5];
u3(-0.212595066131246,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.55708547674082,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.48758846978244,-0.118546535218528,2.11809385543143) q[5];
u3(1.56695552875446,0.628843476860176,-3.38506253967280) q[4];
u3(2.60521095189791,0.683891549871977,-3.16233635253899) q[0];
u3(2.55938590595061,1.69911014089457,-2.66415848696916) q[7];
cx q[7],q[0];
u1(2.72218305739646) q[0];
u3(-1.87741417706953,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.211669494625742,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.51133656748112,3.50455152746510,-1.55390824379875) q[0];
u3(2.65068212092446,2.68008533999789,2.52593759565618) q[7];
u3(1.21760393030920,-0.133394119874954,1.63403747265423) q[6];
u3(1.37285950585703,-0.903098789544168,-2.40347999346925) q[2];
cx q[2],q[6];
u1(3.10794134786719) q[6];
u3(-0.981598843666721,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.84127742129655,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.37481201694791,3.11416628689908,-1.52028255234190) q[6];
u3(2.23213149994102,-3.20035373618556,-1.04699993770926) q[2];
u3(2.32900421503783,-0.321620792283976,2.30178988770792) q[3];
u3(2.65861804200738,-2.03383296219690,-0.678514254441817) q[5];
cx q[5],q[3];
u1(0.596343186625954) q[3];
u3(-0.886987671439127,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.06018927372290,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.64884078338663,1.47566674753440,1.14141663717642) q[3];
u3(2.18445405448252,3.08867446254835,-1.19144286793035) q[5];
u3(2.32811520334738,-0.181708752926750,-1.62466723656811) q[4];
u3(1.56160338807425,0.627347569467218,-5.36809046885522) q[1];
cx q[1],q[4];
u1(2.14222073458225) q[4];
u3(-2.50659629064303,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.15507881471788,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.27741827793914,2.25125023197458,-0.316015622682270) q[4];
u3(2.32833504550402,1.74258049056211,1.99138045699018) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
