OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.15358143252011,-3.76524972170014,1.86831315549386) q[7];
u3(2.53622228319144,-2.50321277214543,3.10345203112235) q[3];
cx q[3],q[7];
u1(3.30252698691532) q[7];
u3(-1.18128381708592,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.51652980721062,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.42333031588663,-2.64131440040354,0.956338561687275) q[7];
u3(2.60764618647671,-3.54407008813852,-0.253057439261992) q[3];
u3(2.33117755184522,1.81237720847231,-3.53401110958752) q[0];
u3(1.50808907243522,2.05818034237034,-2.80210273707893) q[4];
cx q[4],q[0];
u1(3.33848526301809) q[0];
u3(-2.13722166309917,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.67130018595157,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.696380404565340,-0.392738637492667,-1.21208139889398) q[0];
u3(1.75901337483126,-1.79187980405115,3.56768516755572) q[4];
u3(1.51091213525653,2.50696277768737,-0.497875513418084) q[5];
u3(2.42742259717564,0.614478383053154,-2.64019688446464) q[2];
cx q[2],q[5];
u1(3.34319217010511) q[5];
u3(-0.778070297231054,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.89298127852722,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.39305033693880,-0.907844268231633,-0.422423385118203) q[5];
u3(2.59569199013923,-2.64791712164645,3.26795573585937) q[2];
u3(1.83546793954466,-2.76052351652771,0.794623550326773) q[6];
u3(1.57892425496570,-3.27703543009429,-0.163987141786319) q[1];
cx q[1],q[6];
u1(1.17854180003065) q[6];
u3(-0.791591393314225,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.92213936580871,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.87731535946066,0.666342436639987,0.216323816982503) q[6];
u3(2.34707389521640,-1.95492856742987,4.20569247140077) q[1];
u3(0.660952766622066,-1.35491232292817,0.900202539846766) q[0];
u3(0.540794752279976,-2.69576443242817,0.704439948754579) q[3];
cx q[3],q[0];
u1(3.27584694142133) q[0];
u3(-1.19264554343577,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.33890755105241,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.997228776786771,-0.673875256677553,3.46225609337478) q[0];
u3(0.981068047162437,1.33587248349720,1.82802261255123) q[3];
u3(0.239096604028420,2.42696958162023,-0.505177886821654) q[2];
u3(1.29255955524648,-0.0470823749358291,-3.74413103666825) q[4];
cx q[4],q[2];
u1(2.68134962643450) q[2];
u3(0.0311982364804138,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.67165417466643,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.67069389317973,-4.48329512508695,0.833865670995238) q[2];
u3(2.72696797997029,1.16554915605768,-2.80348860918610) q[4];
u3(1.90960182152324,0.304957801202025,2.57554188695212) q[6];
u3(1.88111897934441,-0.441150094188136,-1.36997848170299) q[7];
cx q[7],q[6];
u1(2.32265848868076) q[6];
u3(0.308887931026312,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.47350691542882,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.56770773521089,-0.525027181386571,-1.40794184173100) q[6];
u3(0.800495283502511,4.18317784293912,-0.919383975792741) q[7];
u3(1.67371421930739,0.0419885206718799,0.839570435635071) q[5];
u3(1.21513495791279,-0.467417943905987,-1.58906090466673) q[1];
cx q[1],q[5];
u1(3.62183914024355) q[5];
u3(-3.30515843435773,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.696260144544835,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.407288117532151,1.71970386024444,-1.84913763792836) q[5];
u3(0.217395747134892,3.88922752582353,0.171376519336147) q[1];
u3(0.824087382471036,-1.18170737454105,2.12332687291188) q[6];
u3(1.08894286792255,-1.28531550205597,-0.271818379885107) q[7];
cx q[7],q[6];
u1(3.32270989652505) q[6];
u3(-2.15580029736381,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.18445552295521,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.38462960768035,2.03811567555396,1.24075762247005) q[6];
u3(1.42047230366339,4.14270074108534,0.854102192408123) q[7];
u3(1.36979899304101,1.56044011176873,-3.20694002370647) q[4];
u3(1.02485363989858,-2.38783633020267,2.76893352228053) q[0];
cx q[0],q[4];
u1(2.64972604454681) q[4];
u3(-0.0460784283075573,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.34886214232435,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.57780047641191,1.79108702597812,-1.91559613519419) q[4];
u3(0.700358463736960,1.34117723299207,-0.406075394964853) q[0];
u3(1.84969544014748,2.07414896607059,-0.0148146387631420) q[2];
u3(2.71168527639596,-0.117173070539213,-3.97958182997708) q[1];
cx q[1],q[2];
u1(-0.00848402944696924) q[2];
u3(-1.90071705667340,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.809643444833019,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.330436418222458,1.01239020742640,-0.458296559759083) q[2];
u3(2.07254612054924,-1.15092836543445,-1.70664652340330) q[1];
u3(1.42034891102306,0.535282721096519,1.10325871492985) q[3];
u3(1.16553853972354,-1.26462480399719,-0.680393137903667) q[5];
cx q[5],q[3];
u1(-0.362437915825045) q[3];
u3(1.32777856337599,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.68133274910867,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.57303226476176,1.32517184208148,-2.77205578593159) q[3];
u3(1.38843644530731,0.899697989606570,-3.10667591672228) q[5];
u3(1.58309268317145,2.01064654306817,0.120150803261335) q[6];
u3(2.86602634363016,0.298089730906158,-4.06965610256827) q[2];
cx q[2],q[6];
u1(-1.00899953983878) q[6];
u3(0.811152596451329,0.0,0.0) q[2];
cx q[6],q[2];
u3(3.70932205148148,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.622391646012781,1.03526142703740,-3.22021881489680) q[6];
u3(2.26551729901517,4.98424684352698,1.17435857444515) q[2];
u3(1.62977865159610,1.36811949238497,-3.64171301267131) q[0];
u3(2.16825667014848,2.71037305220339,-2.70832525743625) q[5];
cx q[5],q[0];
u1(2.72832290091320) q[0];
u3(-2.38172003132156,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.03942525751696,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.02792531905027,0.677083387855098,1.15129044973713) q[0];
u3(1.62977268450146,4.19507389861908,0.231963168924959) q[5];
u3(2.63059024910124,-3.39943996202772,2.80258256295991) q[4];
u3(0.589329307807487,-0.00703340038117684,1.57758381334270) q[1];
cx q[1],q[4];
u1(-1.19326613799898) q[4];
u3(0.629628300359978,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.98997483535033,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.889102239902077,-1.34471823044703,0.990715762476554) q[4];
u3(2.48097140439319,5.15893612490304,0.242423560781999) q[1];
u3(2.32226411051120,-3.89985791292034,2.19506013913958) q[7];
u3(0.500552852364669,-2.26756135235287,3.04197819090270) q[3];
cx q[3],q[7];
u1(3.43580356328617) q[7];
u3(-1.16737497630695,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.23709827992663,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.77758622998121,1.67706961372651,-3.98583614532811) q[7];
u3(1.15957260925815,-0.711650590161277,4.91276993406328) q[3];
u3(1.15049090978774,1.78177967126619,-3.41133769012993) q[2];
u3(2.60859750942506,-1.98102862116387,3.13126087775119) q[4];
cx q[4],q[2];
u1(1.06150902488160) q[2];
u3(-0.623634538681065,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.90838998061834,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.20386920349230,-2.90356740430910,2.92832009449819) q[2];
u3(1.90197119682060,-0.368904807446722,-1.76629888516828) q[4];
u3(0.481591579640204,2.90295897708694,-2.71793970546995) q[3];
u3(1.17743948526901,1.01679402975964,-1.70340341542279) q[6];
cx q[6],q[3];
u1(2.44706624932531) q[3];
u3(-1.88354794501399,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.0294602067135430,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.70576142698677,3.16424945500622,-1.61169303208150) q[3];
u3(1.32881069787371,1.13400988217901,2.79515973443727) q[6];
u3(1.66152379478127,-0.715463295900538,-0.596084916820844) q[1];
u3(0.475537279096580,-2.35764390441520,-2.80122629183663) q[7];
cx q[7],q[1];
u1(1.65175398494862) q[1];
u3(0.235129990543973,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.09225606613857,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.23768231858396,3.08511157871743,-3.06619734275491) q[1];
u3(1.40873411631503,-0.166306333725159,4.70442067434846) q[7];
u3(1.55427817092874,2.33038247541353,-3.67285432869162) q[5];
u3(0.448730377388682,2.65011272092342,-1.44767729285218) q[0];
cx q[0],q[5];
u1(1.13962366349918) q[5];
u3(-0.394294355690741,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.56928570228591,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.12467468256413,-0.0462607836627347,-1.78281855945962) q[5];
u3(0.961578196220199,3.08002447575298,-0.503801460166409) q[0];
u3(2.36420830437891,3.48217764944701,-0.857166816828540) q[7];
u3(1.01032683933837,1.96369554218432,-1.48918769311811) q[2];
cx q[2],q[7];
u1(1.01996796424591) q[7];
u3(-0.540951653342307,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.83415120595710,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.532981165903635,0.639434377154724,-0.630236784136037) q[7];
u3(0.550430436545713,4.52865898335963,-1.66869915890240) q[2];
u3(3.04390898488362,-0.498123872943982,0.637514593640884) q[3];
u3(1.29253972397230,-1.78151869088029,-3.01755351763313) q[0];
cx q[0],q[3];
u1(0.0631580165470531) q[3];
u3(-1.31765732476864,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.35157246552370,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.29326598807522,3.33908430130023,0.257102117244151) q[3];
u3(0.538671230407873,2.96628664756983,-1.27010930747632) q[0];
u3(1.41502826940791,1.87957743460750,0.258109183802273) q[5];
u3(0.929900300983742,0.779209478413754,-4.59176546872433) q[6];
cx q[6],q[5];
u1(3.13837092811981) q[5];
u3(-1.18162353685499,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.58776147393841,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.37078684386340,0.606166465414705,-1.43794649922200) q[5];
u3(2.03118131517981,-5.25310975607423,0.221261007470424) q[6];
u3(1.28537511208055,1.25545862775358,-2.17814785708570) q[4];
u3(1.73287796260703,-4.20056982231807,1.92587600992009) q[1];
cx q[1],q[4];
u1(3.27851662932106) q[4];
u3(-0.904386397376415,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.85940947253455,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.15552388226228,1.38925937845781,-1.76177574498472) q[4];
u3(2.62283768476988,-4.81607062837383,-1.00219665813985) q[1];
u3(1.19157323277522,2.89565013479935,-2.56543270531857) q[4];
u3(0.343694077479635,0.269123611778049,-0.963610304329358) q[2];
cx q[2],q[4];
u1(1.71393661727897) q[4];
u3(-0.0159430300051937,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.93359177918465,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.29247591137926,4.77644998787581,-1.40014378804631) q[4];
u3(2.62957233013024,-1.98593586236933,0.737000002225174) q[2];
u3(1.72379435624509,-0.476590201687745,1.35962079410261) q[3];
u3(1.65362216445653,-1.94100583443959,-2.73209248927465) q[1];
cx q[1],q[3];
u1(0.107002154766917) q[3];
u3(-0.398560120617011,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.09855401307834,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.88977491955348,0.771067564862500,-2.42972401534958) q[3];
u3(1.78063259659446,-4.18002358286455,2.00427484879954) q[1];
u3(1.83854904976624,-0.708357689982088,-1.73350368690167) q[6];
u3(1.20168017580198,-4.63106932186178,1.09867851937457) q[5];
cx q[5],q[6];
u1(3.70895207805620) q[6];
u3(-3.38866611349267,0.0,0.0) q[5];
cx q[6],q[5];
u3(-1.18419608654905,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.62843322018554,2.74188954077810,-1.96777979586246) q[6];
u3(0.923704786622902,2.67591945772205,1.65293818318275) q[5];
u3(1.38128170651959,-2.20129253896183,0.328051033459968) q[0];
u3(0.440248158909852,-2.43018854578658,0.374108862634537) q[7];
cx q[7],q[0];
u1(4.17459243834608) q[0];
u3(-3.83872181071398,0.0,0.0) q[7];
cx q[0],q[7];
u3(-0.555378425359408,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.62861331365878,-3.86197553509747,1.31932839104944) q[0];
u3(2.41881822757840,1.35591692445489,-0.578124583554211) q[7];
u3(0.783783281043688,0.265898599993379,0.316488843873526) q[7];
u3(0.863969286238494,-0.663187622378969,-0.758139180842125) q[4];
cx q[4],q[7];
u1(1.24342516350282) q[7];
u3(-1.58634346626270,0.0,0.0) q[4];
cx q[7],q[4];
u3(3.40766968775527,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.79774523184259,1.21023034147223,2.53381628115314) q[7];
u3(1.63390359055687,-1.26944034879610,-0.0139936514252200) q[4];
u3(2.29068303959029,1.91090508336325,-3.88359779114726) q[6];
u3(1.12873118345502,-2.17254621884849,3.64234192272692) q[3];
cx q[3],q[6];
u1(1.62506276213736) q[6];
u3(-2.06559213026693,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.345418583650686,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.28743299193390,0.108535776639381,-0.959881976636942) q[6];
u3(1.46546649620843,-3.47081352950409,1.16030300189607) q[3];
u3(2.49467727042015,2.06776329906071,-1.82681250675493) q[5];
u3(2.56005446136547,-0.0543765765423392,-4.25044924257170) q[0];
cx q[0],q[5];
u1(3.96526559553265) q[5];
u3(-1.42928392262142,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.03097889237208,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.96442735722031,0.563220441618504,-4.34814406811012) q[5];
u3(0.688269500862807,0.590231460640292,1.02588868111143) q[0];
u3(1.71522803379858,0.888583552989761,-2.74284716202658) q[2];
u3(2.28576094651414,2.43950547048202,-3.82079737960829) q[1];
cx q[1],q[2];
u1(-0.0767552770658413) q[2];
u3(-1.84407861477458,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.20106960719289,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.366701287657358,-2.11664581429433,2.47506114488467) q[2];
u3(2.37121363026403,-2.88044092766284,-2.02152410942388) q[1];
u3(0.455066836337381,1.83803015292279,-3.07852139704359) q[6];
u3(1.82178583215682,-2.40435583557332,3.71093503526891) q[1];
cx q[1],q[6];
u1(0.233658747690764) q[6];
u3(-0.659804748183892,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.50204924545194,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.60127486232761,-0.0482208626005813,2.16771332755676) q[6];
u3(0.988039047216483,1.16351973324535,-0.339839863971353) q[1];
u3(1.65002450647854,2.91110941877959,-0.850244945878732) q[3];
u3(2.24533746821791,2.87262620599631,-0.541486782148324) q[7];
cx q[7],q[3];
u1(0.311221154647220) q[3];
u3(-1.18477905013999,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.81256231114383,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.16661220072855,0.615912790166629,0.218754374596183) q[3];
u3(1.63657701503500,1.76841190032063,3.91229722336993) q[7];
u3(1.54880847616033,1.85709983981694,-3.90964045268984) q[5];
u3(1.89838290827086,2.86148570142335,-2.48408342617976) q[4];
cx q[4],q[5];
u1(1.74105715222607) q[5];
u3(-2.92770955169433,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.0584469327362485,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.25013041798077,-2.46367935046011,1.30107247876986) q[5];
u3(1.22774938062684,0.853245844349308,1.52483277972337) q[4];
u3(1.76771929394341,-1.51449414585597,1.61710233843049) q[0];
u3(2.53555887296276,-2.03902843055989,-2.34643490757082) q[2];
cx q[2],q[0];
u1(1.80897464931178) q[0];
u3(-2.52973130433353,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.14846297681652,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.03175691720956,0.354781768339317,-1.38725480806871) q[0];
u3(1.69523390841834,1.04183154047492,3.55780065931012) q[2];
u3(1.22064643211007,-0.0238708679663140,0.532393615472334) q[7];
u3(1.16011158680645,-2.10248411156264,-1.20689876646303) q[5];
cx q[5],q[7];
u1(4.62897240387093) q[7];
u3(-3.83538875863347,0.0,0.0) q[5];
cx q[7],q[5];
u3(-0.642987528168572,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.71354246161948,1.34869957458574,-1.97416726437689) q[7];
u3(0.325029616456487,-1.36791633693175,0.373877541535582) q[5];
u3(1.97481708244446,3.95822816118153,-1.94564317694784) q[6];
u3(0.493247170839379,-2.17645569767841,3.70700287157086) q[0];
cx q[0],q[6];
u1(-0.0559072636160658) q[6];
u3(-1.90784632163891,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.946406271072859,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.802099831416594,-2.61901594806798,-0.229774161495535) q[6];
u3(1.96604793071094,-0.384538708556763,-2.37008962367675) q[0];
u3(0.0807857470249941,2.68769249850018,-2.48739425923731) q[1];
u3(0.783750764600230,-1.62400926820690,1.12260565473514) q[4];
cx q[4],q[1];
u1(0.353975311167429) q[1];
u3(-1.51018623325510,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.47732759801909,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.41227284588025,-1.76302605091646,0.357254063294410) q[1];
u3(0.933401422164306,0.893678051604634,2.33406543804834) q[4];
u3(0.281824320813648,0.773140326314208,-0.593657716770428) q[3];
u3(0.886417543839823,-3.23326607033164,0.636931394615484) q[2];
cx q[2],q[3];
u1(3.23793982418070) q[3];
u3(-3.63072668369330,0.0,0.0) q[2];
cx q[3],q[2];
u3(-1.00217095042423,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.81543756697308,-2.55917059715498,3.60183156970346) q[3];
u3(0.495483989219616,3.92865035587609,-1.62895481581865) q[2];
u3(1.55080479919673,2.02143151865688,-1.27392850981264) q[6];
u3(1.48171409144019,1.27076842135520,-1.06269531117144) q[2];
cx q[2],q[6];
u1(0.0614007665607130) q[6];
u3(-1.40630396851933,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.44238362509187,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.67756917475828,2.01421423427376,-2.75951316519757) q[6];
u3(2.05725970403949,4.71190449245790,-0.635839005050135) q[2];
u3(1.78406448345904,2.65932075932322,-1.39088007720900) q[1];
u3(2.24847394977755,2.30794134926777,-0.637686112470492) q[7];
cx q[7],q[1];
u1(-0.222735086700070) q[1];
u3(-2.25304902979466,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.36567272715577,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.03110879808526,0.111384663299795,-1.39592352880446) q[1];
u3(1.23514438970819,-1.10477493769775,2.63680936103552) q[7];
u3(2.57879796485318,-2.16354446180575,-0.413339504080594) q[3];
u3(2.13925578524999,1.97467613791003,3.81302360659244) q[0];
cx q[0],q[3];
u1(0.765525392125813) q[3];
u3(-1.19108961628761,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.79268665542135,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.60356099167065,1.63389475655699,0.965561230059551) q[3];
u3(1.55167752920415,-2.35535261240672,2.88075716951085) q[0];
u3(2.15097018263875,-2.75750993497694,-0.0557606842064067) q[5];
u3(1.92343490960481,-3.04227494379676,0.140743913605434) q[4];
cx q[4],q[5];
u1(0.903830590233569) q[5];
u3(-1.34520932171290,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.357108289419164,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.71791407325339,0.0697306406993503,0.151546213072920) q[5];
u3(2.68361405285116,-0.961065485718479,-3.62623328558645) q[4];
u3(1.11815713376719,2.38262947958700,-3.43387942624609) q[0];
u3(0.327818576399950,-0.118970437299787,-1.92368172164097) q[4];
cx q[4],q[0];
u1(2.68357752168328) q[0];
u3(-1.48334512072060,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.631247335143271,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.73361772705472,-3.12284290502479,0.412581198161299) q[0];
u3(0.780286172709804,-2.12424100532017,-1.38359737487311) q[4];
u3(0.793496022127880,0.978779975844107,1.38131147591090) q[7];
u3(0.834416834205280,-1.41583544899152,-2.01965422691981) q[2];
cx q[2],q[7];
u1(1.12391114192551) q[7];
u3(-3.27399668954390,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.75716404764076,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.873769255026825,-0.0656020642241391,1.90116084072857) q[7];
u3(1.55419040436369,-2.37483972815944,-0.692762467506663) q[2];
u3(1.50049920216280,-0.431048852520601,-1.45405034486664) q[6];
u3(1.42383823434508,-3.82352640296401,0.449019353545660) q[1];
cx q[1],q[6];
u1(3.62049423375802) q[6];
u3(-1.20475799866305,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.14884662850136,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.98504711441415,-1.93965592076556,1.39557625499806) q[6];
u3(2.58408072844459,4.26474251330563,-0.275596811355175) q[1];
u3(1.55224872865427,-1.66733987334257,0.390797303597696) q[5];
u3(1.16580692377280,-1.92714516618984,-0.177832541274527) q[3];
cx q[3],q[5];
u1(2.89497939213915) q[5];
u3(-1.58607077820364,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.17014642925840,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.691216858356725,-0.403353268157550,-0.351352242618453) q[5];
u3(1.64297393366375,-1.32819679577172,-0.702369364266456) q[3];
u3(1.38358936955686,-0.616700977201126,1.95534549256184) q[3];
u3(0.919547230546290,-0.703746783356679,-0.210997638762746) q[4];
cx q[4],q[3];
u1(1.60142861447625) q[3];
u3(0.171060446092238,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.737547800134138,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.45676905359064,4.15018720466078,-1.93690267691204) q[3];
u3(1.15637368810267,-3.01673439625142,-0.246779223343215) q[4];
u3(1.13656908077787,-2.70574605841059,2.12770214419559) q[0];
u3(0.468477394093601,0.847073251839033,-3.37297438266198) q[6];
cx q[6],q[0];
u1(1.65151323126273) q[0];
u3(-2.46167771856805,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.10108824004781,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.478676608203675,-1.37642341787695,-1.20330755133052) q[0];
u3(1.98392077978664,-1.25701718574154,0.126680860447860) q[6];
u3(1.26258399489553,-0.558958328106446,1.50281500095769) q[5];
u3(1.46608196375655,-0.964596685770223,-1.83590038713238) q[7];
cx q[7],q[5];
u1(2.17265097871542) q[5];
u3(-1.72246543226485,0.0,0.0) q[7];
cx q[5],q[7];
u3(3.53675036657470,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.436802271702908,-0.271876656668821,1.22505324798493) q[5];
u3(1.93463020906467,0.114098002570480,-5.32763699834300) q[7];
u3(0.370259884505336,-0.0776146744862782,-0.754332452938676) q[2];
u3(0.345116148569789,-1.47582062797381,-0.465870604823815) q[1];
cx q[1],q[2];
u1(1.51593346707643) q[2];
u3(-0.726847952389399,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.239619171841923,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.0781539333503886,2.10382369702925,-0.756392280418729) q[2];
u3(2.19291315326840,-3.54588700818146,1.38854088857851) q[1];
u3(1.72381867190231,2.71016362991811,-2.90442489349544) q[1];
u3(2.13170600283412,2.82087683238688,-3.45796854831533) q[0];
cx q[0],q[1];
u1(-0.112670944559723) q[1];
u3(-1.79563482869847,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.30700822765892,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.90398532688283,0.303020440045657,-1.22852788432590) q[1];
u3(2.27138354732477,-4.42884766759953,1.15090845740117) q[0];
u3(0.860873502854746,0.990913087652379,-2.66225292720858) q[4];
u3(1.22848059391097,2.57799131610096,-3.23761357960360) q[5];
cx q[5],q[4];
u1(2.35626063537680) q[4];
u3(-1.74790852609358,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.10694529615088,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.22860073017454,0.156720430950231,-0.0428323453834358) q[4];
u3(2.62982396968122,1.43195760247619,-1.02852384648227) q[5];
u3(1.94947688474666,-0.289365782585003,1.37291967994614) q[6];
u3(1.62614940286425,-1.99052885660773,-1.53328631627182) q[7];
cx q[7],q[6];
u1(2.82344781591089) q[6];
u3(-1.47286578514040,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.626138858505187,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.03934471947869,0.438816342901531,3.03479314962066) q[6];
u3(1.65820155791474,0.731018908553720,-5.49898487613752) q[7];
u3(1.51937710644828,1.41659346927375,-3.07249723292082) q[3];
u3(2.85336191247044,2.52916173996821,-2.85843974782810) q[2];
cx q[2],q[3];
u1(0.572718117576644) q[3];
u3(-1.18364891618148,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.03887752552890,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.952787449453743,1.86136951140491,-0.00313188996380975) q[3];
u3(1.41737192350070,1.89362075435180,4.08170835455830) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
