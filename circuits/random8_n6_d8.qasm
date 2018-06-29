OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.48949015803309,-0.696834381053006,1.49012778668688) q[5];
u3(1.37882808062495,-2.06382753473904,-0.702100586462075) q[1];
cx q[1],q[5];
u1(1.88624214912259) q[5];
u3(-1.74214782271991,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.02584639361178,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.17315811466830,-1.17297874223477,-0.636063486583977) q[5];
u3(1.25789531328603,1.53073895508270,-3.06274707277878) q[1];
u3(0.339958410391770,1.87005127492539,-1.50847433733467) q[0];
u3(1.41873030682574,1.46942103946407,-1.60968129248822) q[4];
cx q[4],q[0];
u1(1.69911471329974) q[0];
u3(0.688449111047390,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.02802660657344,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.37467791642140,1.35275154049809,1.19652929365193) q[0];
u3(2.10905455235284,0.0423233000539924,-3.81816865500724) q[4];
u3(1.38266857117439,-0.605356027969033,-0.747018082356360) q[3];
u3(0.149358923113268,-1.91396792654303,-3.16542207826740) q[2];
cx q[2],q[3];
u1(2.74187280585219) q[3];
u3(-1.74095292344438,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.783626405285329,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.23966338914254,-0.259647671222456,3.39518156988411) q[3];
u3(1.45791750325330,1.19179937033653,4.19536863254056) q[2];
u3(1.70950939801235,2.58648331575832,-1.22715891935144) q[3];
u3(0.476391755379352,2.77471991511517,-3.22878734678209) q[0];
cx q[0],q[3];
u1(2.27183226346305) q[3];
u3(0.169394352343396,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.31741983682945,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.05721452807569,-2.06907693245867,0.487891324327302) q[3];
u3(0.660032825570139,-3.73276664513830,0.721558597464749) q[0];
u3(1.33443951683310,0.736700610382545,-1.46205401091544) q[4];
u3(0.309968111917548,0.709735748507385,-3.70126201755105) q[2];
cx q[2],q[4];
u1(0.990619469590418) q[4];
u3(-0.458148179377732,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.59819726757483,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.922228409973458,0.533758560541243,-0.512256451166667) q[4];
u3(2.11391771186331,1.42987905270405,-0.438675590706033) q[2];
u3(1.88382240300994,-2.27902245843446,-0.772424966052021) q[1];
u3(1.56011495660613,-3.11488959378822,-0.855677219341198) q[5];
cx q[5],q[1];
u1(3.34536979933528) q[1];
u3(-4.21470359368532,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.507835592548685,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.26986240188431,3.69192788380079,-0.935292588975347) q[1];
u3(1.67562883191989,1.86311682418327,-4.29901811318914) q[5];
u3(2.21865098551332,2.18585457946152,-0.942050179411070) q[2];
u3(2.78399095436138,2.01276296833160,-1.61774629959806) q[5];
cx q[5],q[2];
u1(2.03831548732957) q[2];
u3(-2.71549171382135,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.0291017100267352,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.31513727990874,-0.506487591508672,0.937979224738080) q[2];
u3(0.823474171561440,0.389960526836951,-4.80454681175972) q[5];
u3(1.40221213171832,1.72316614398496,-3.31826215762806) q[1];
u3(2.54695571367640,-2.57765849530966,2.98855563357953) q[0];
cx q[0],q[1];
u1(3.24902642145373) q[1];
u3(-2.07431052698343,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.346616051364035,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.70949853672051,-1.21771963534481,3.58657642556906) q[1];
u3(2.57360116627769,1.00183688589365,2.21542994046060) q[0];
u3(1.61109689845487,1.67237649811481,0.578763959950528) q[4];
u3(1.98803239946612,0.110889205233135,-3.72391279760416) q[3];
cx q[3],q[4];
u1(1.15100673090757) q[4];
u3(-0.839110568847341,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.200944847345966,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.11357145854554,2.13491212054873,-2.25970621610680) q[4];
u3(3.10122704580719,-3.48031878160920,2.79708279827183) q[3];
u3(0.788154804641336,0.237870492853613,0.00958776810841788) q[0];
u3(0.580661803764688,0.0308527562984174,-1.35542772293100) q[5];
cx q[5],q[0];
u1(2.02902719676478) q[0];
u3(-3.34702140857125,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.33242705975409,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.76387050480821,0.184264982886788,0.746553616474220) q[0];
u3(1.62299906045163,0.798140962055565,-0.0186469463366229) q[5];
u3(0.881019487975038,1.16702874322659,-0.253532219839027) q[2];
u3(0.605319728337206,0.769921007426854,-1.78693992109172) q[1];
cx q[1],q[2];
u1(-0.219584402022233) q[2];
u3(-2.33079049954917,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.01020132028230,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.675567947178524,-0.205379732024583,-0.805683322531075) q[2];
u3(1.51923182831840,0.648658708178673,0.226574450506745) q[1];
u3(2.32880944185981,0.326755428081617,2.17830066566616) q[3];
u3(1.99517909259827,-2.93591503780353,-3.15574069860295) q[4];
cx q[4],q[3];
u1(0.834938943529248) q[3];
u3(-1.65604008467215,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.28237533742407,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.504785851580175,-0.742188574234002,2.82139220519715) q[3];
u3(1.22345738313863,0.0382111966482515,-4.24268688560564) q[4];
u3(0.366895314737003,-0.968692827182891,0.769142929567987) q[5];
u3(1.03521513285007,-3.01343537672312,1.67947942959424) q[3];
cx q[3],q[5];
u1(0.305836598145479) q[5];
u3(-1.58467874616440,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.71806280201460,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.92128366439260,1.94500452864225,-0.227252947360838) q[5];
u3(1.47545821063661,0.541495193485461,4.51311230316143) q[3];
u3(1.42626060008711,2.49276666162124,-1.60639310954524) q[0];
u3(1.06749042597831,0.777395111157198,-1.04611619277007) q[1];
cx q[1],q[0];
u1(-0.557665512597662) q[0];
u3(0.858866597891123,0.0,0.0) q[1];
cx q[0],q[1];
u3(4.10970806425439,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.971855390607066,1.23439256464568,-0.838925404324395) q[0];
u3(1.77501116845218,-3.08217515265833,-0.0154284855286553) q[1];
u3(2.23977715166788,0.458715027691788,2.42360520648648) q[2];
u3(2.22367976507259,-1.85226690803716,-2.56374506720252) q[4];
cx q[4],q[2];
u1(2.46913401581009) q[2];
u3(-0.114255256312025,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.59314600063834,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.22664942232927,2.26430558073057,-0.751489518295824) q[2];
u3(1.92083027399391,-1.91123587780439,-3.45991399240981) q[4];
u3(2.01985823587256,0.939458660984706,1.38479262058523) q[0];
u3(1.55901572359113,-1.23281989667861,-2.61185529349050) q[5];
cx q[5],q[0];
u1(2.30804513728571) q[0];
u3(-3.05728212215799,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.51110960226602,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.80092044937906,-1.68699326833937,3.06872977113722) q[0];
u3(1.25332245570880,1.13804736531568,0.307726451234058) q[5];
u3(1.22982336260191,-1.31377869697820,2.07693298881285) q[2];
u3(1.10451918142771,-1.67253539609922,-2.91136750254792) q[1];
cx q[1],q[2];
u1(0.673903373784879) q[2];
u3(-1.59720504986550,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.0572110504001846,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.94297008845750,-1.72456437840437,-0.407359744484671) q[2];
u3(2.29757934704354,-6.10860608518456,-0.164551462587450) q[1];
u3(1.36537949169899,-1.01323506426078,2.21344735123006) q[4];
u3(1.16826236215811,-2.05138084766698,-1.27223951402517) q[3];
cx q[3],q[4];
u1(0.930507772052960) q[4];
u3(-1.57429278276243,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.83364485710201,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.07369531032167,0.919535230409626,-2.95571982878011) q[4];
u3(1.44024918248006,-0.929586669987311,-4.77824752172458) q[3];
u3(1.66571697796019,-1.88647945171388,0.180969191448265) q[5];
u3(1.61763317772056,-3.63340909065168,-0.000533404823991335) q[4];
cx q[4],q[5];
u1(1.69207151438543) q[5];
u3(-3.50661091961845,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.10549422373937,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.573037464270462,-1.51644708996443,-0.957184040750739) q[5];
u3(2.94506939529910,-1.14614211002390,-4.50832632908907) q[4];
u3(0.527467863678543,0.967048759824714,0.0712792793210599) q[3];
u3(1.01772306739401,0.302136318874365,-1.69680546004468) q[2];
cx q[2],q[3];
u1(1.77797024949347) q[3];
u3(-3.07468840856267,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.821055732011280,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.63999543113966,-1.46792397372328,0.644511084882625) q[3];
u3(1.19135417751955,-1.34468738611471,-0.984875415589948) q[2];
u3(0.753207201770039,-1.64338799037518,-0.580651172852014) q[1];
u3(0.766972954541078,-2.84275205853749,-0.947675807192763) q[0];
cx q[0],q[1];
u1(0.794670701984513) q[1];
u3(-3.51293352861490,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.91191465308742,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.987526297904488,0.142249608467848,1.54160209176085) q[1];
u3(0.854640546739736,0.318117003653014,1.75129650713469) q[0];
u3(2.04839374755901,-2.18127251216553,-0.530991987842327) q[3];
u3(1.09249911899446,-3.92710006447654,-0.955871583172238) q[0];
cx q[0],q[3];
u1(1.48556772671342) q[3];
u3(-0.255389786370770,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.89602300123639,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.24675249892247,3.10033108645000,-1.94650041794984) q[3];
u3(0.792713565916507,0.841908212149425,-4.15071547758331) q[0];
u3(2.33173314516557,-1.21246072299733,0.756490378604393) q[1];
u3(2.33134720548604,-1.21045719546299,0.0589756732528925) q[2];
cx q[2],q[1];
u1(1.43020546679651) q[1];
u3(-0.532139600100275,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.70764947921544,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.06466212232830,2.19516789361722,-3.50634599646324) q[1];
u3(0.865244054297699,2.41104034044724,-3.14492134634878) q[2];
u3(1.38242508854465,-0.344940929175902,-1.68130189354760) q[4];
u3(2.77418718120481,-4.55061091166407,0.514399914710927) q[5];
cx q[5],q[4];
u1(2.71252676259855) q[4];
u3(-2.09251858800124,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.474139270029742,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.62467006708231,-0.181773790214274,-1.78492280591056) q[4];
u3(2.45629927867053,-0.870181933022660,-3.69937278133240) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
