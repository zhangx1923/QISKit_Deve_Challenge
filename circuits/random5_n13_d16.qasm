OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(0.740227404205154,2.41032226843333,-0.627745670559100) q[12];
u3(1.86513122317393,2.13860924169117,-1.54197722293090) q[9];
cx q[9],q[12];
u1(3.21725038243881) q[12];
u3(-1.23206404243958,0.0,0.0) q[9];
cx q[12],q[9];
u3(2.29130216191984,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.83253473758990,4.72874549548012,-1.34841141251652) q[12];
u3(1.90588098051211,0.460966912756585,-5.39462768322042) q[9];
u3(0.966343713943062,-0.720284517903263,3.05592253958153) q[5];
u3(2.41758542893510,2.93232201617125,1.33224660695988) q[6];
cx q[6],q[5];
u1(1.39316266835162) q[5];
u3(-0.330227460726542,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.37036712133891,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.20104263739958,-2.61822037599684,1.27366659439067) q[5];
u3(1.23278523762825,-1.97713894630138,1.62679786797927) q[6];
u3(2.33023424272521,2.77811697279912,0.302156285613611) q[3];
u3(1.99720167041598,3.18286079266346,-1.56286558464461) q[2];
cx q[2],q[3];
u1(0.772418311246483) q[3];
u3(-0.118322468895616,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.42530231339718,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.54611582903870,-2.31383196729258,0.821851671058592) q[3];
u3(1.77832908782290,-1.55058163394484,-3.04929881645864) q[2];
u3(0.146397213995954,-0.733188326554656,0.334217690902378) q[8];
u3(0.293940890131205,-2.25245529469402,0.627097876068601) q[0];
cx q[0],q[8];
u1(1.74038634875504) q[8];
u3(-3.35193628182791,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.479502060080286,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.44863088515938,-0.268895787181256,2.74104094204321) q[8];
u3(1.69445755997556,-3.89985636993377,-1.13014215442710) q[0];
u3(1.43658612317736,0.143021618818827,-1.71334374292185) q[10];
u3(0.938222854794475,-3.23865395936474,1.20499144969268) q[4];
cx q[4],q[10];
u1(-0.0875729057078989) q[10];
u3(-2.53827104317955,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.44547858164035,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.30549821062843,0.105467575766764,2.70776101851106) q[10];
u3(0.845161440717181,-1.34195964708229,-3.84976631557468) q[4];
u3(2.59221079872654,-1.45184040278363,2.90696669573356) q[7];
u3(2.34259234042803,-1.45116998294312,-0.271082384705054) q[11];
cx q[11],q[7];
u1(0.403232840080144) q[7];
u3(-1.40098013197687,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.233080288056870,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.16121695233277,1.43816669753119,-4.31079108066659) q[7];
u3(1.78157157712675,-0.986501227633040,-0.765275657497343) q[11];
u3(2.02323366416443,1.51321474623718,1.51739361286890) q[5];
u3(2.07193505357535,-0.363010489419073,-3.32912577757822) q[0];
cx q[0],q[5];
u1(1.97095916654447) q[5];
u3(-2.15707242988242,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.181129600852525,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.25292223798747,0.835793829924703,-1.81697719927948) q[5];
u3(2.83232154770262,-0.623737699195137,-2.44525551898490) q[0];
u3(2.17203076875701,1.37539422192815,-3.22704051717618) q[1];
u3(1.44346576294714,2.78551747093175,-2.92682528019214) q[6];
cx q[6],q[1];
u1(-0.308961022177735) q[1];
u3(1.18116005329271,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.81007292501236,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.778643894031923,1.97123009953117,-0.795962181687042) q[1];
u3(2.31312520796242,2.21711595054249,-3.91838384604361) q[6];
u3(2.88039533271046,-0.0471718933028937,2.97561318528897) q[3];
u3(2.92001165161330,0.256777872550167,1.30657041043725) q[10];
cx q[10],q[3];
u1(2.20832734821291) q[3];
u3(-1.71977641367614,0.0,0.0) q[10];
cx q[3],q[10];
u3(0.398462784637767,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.91387947388979,-1.89353324356580,3.88804637050981) q[3];
u3(1.36347920053476,-3.83297149963238,-0.878912183577484) q[10];
u3(1.76975549116790,0.555021674844497,-0.770713198601050) q[11];
u3(1.73610792265774,-0.662773138595797,-3.74837986419301) q[2];
cx q[2],q[11];
u1(0.416013907706233) q[11];
u3(-1.39686765783177,0.0,0.0) q[2];
cx q[11],q[2];
u3(3.06517891631577,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.88456945201521,-0.339008597938384,1.67967931650390) q[11];
u3(1.75782418921722,3.66083365285683,0.844959055557422) q[2];
u3(1.65502042970024,2.92390711976799,-2.47678747523854) q[9];
u3(1.19790682090064,3.29956269922641,-2.85392596923251) q[8];
cx q[8],q[9];
u1(1.79149201887371) q[9];
u3(-3.20222865189618,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.725073560624015,0.0,0.0) q[8];
cx q[8],q[9];
u3(0.801441464310797,-2.47658621119951,2.51103706371552) q[9];
u3(1.48504914413781,3.94735309680991,0.559977867519194) q[8];
u3(2.97933934043709,3.61465213293765,-2.47655422738102) q[12];
u3(1.55618053887892,-2.50694757720588,3.56680895530621) q[4];
cx q[4],q[12];
u1(2.93350985281578) q[12];
u3(-1.84039301378487,0.0,0.0) q[4];
cx q[12],q[4];
u3(0.836780225134433,0.0,0.0) q[4];
cx q[4],q[12];
u3(3.02931211645550,-0.179308924577244,-1.69094414162559) q[12];
u3(2.81559911092442,-5.58812024023696,0.197990778353558) q[4];
u3(1.46455716386085,1.47270007431267,-0.400060458620088) q[4];
u3(1.17306499572952,-0.107912378781987,-3.25338842666417) q[7];
cx q[7],q[4];
u1(2.38168357858232) q[4];
u3(0.280598391353100,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.10151820806088,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.590770691018411,-0.780729050443663,-0.551260478534690) q[4];
u3(1.92467240356239,-3.26170937129282,-1.05246377148661) q[7];
u3(1.44728595396158,1.57187581352250,-0.0955276194049015) q[10];
u3(0.900574837597030,-0.692268777256971,-2.85393577416965) q[2];
cx q[2],q[10];
u1(1.71201507404398) q[10];
u3(-0.542476776293299,0.0,0.0) q[2];
cx q[10],q[2];
u3(2.63940577677858,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.18893920689071,3.53670853805121,-1.86890059972960) q[10];
u3(2.29394387225207,2.30487323030579,1.15336827294992) q[2];
u3(1.51831057424306,3.22579834432000,-2.08732993022055) q[6];
u3(1.42227548217671,3.31347197109435,-2.48241456447038) q[1];
cx q[1],q[6];
u1(3.39087422183752) q[6];
u3(-1.49119490670776,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.46638115686144,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.26887649472650,0.191113421461406,3.92718455991515) q[6];
u3(0.411327338654564,3.56232992504699,-1.20036563923477) q[1];
u3(1.02825342215288,-2.02903025551759,0.284426342651148) q[5];
u3(1.26608924485055,-3.18621099149399,-0.893662601045897) q[3];
cx q[3],q[5];
u1(1.95352231758799) q[5];
u3(-2.23959438164280,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.299838194032601,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.41670424334078,-0.573424937057367,-1.14151076935364) q[5];
u3(1.17353560094763,1.53199192443302,-3.81333471841357) q[3];
u3(2.14483652343221,0.419145776326730,1.94825542943822) q[11];
u3(2.09756304053579,-0.766206019166192,-2.08106248817002) q[12];
cx q[12],q[11];
u1(1.54569036827477) q[11];
u3(-1.07134704457187,0.0,0.0) q[12];
cx q[11],q[12];
u3(3.93065172713381,0.0,0.0) q[12];
cx q[12],q[11];
u3(1.59492795118375,1.83287129958898,-2.31718990095734) q[11];
u3(1.87009111554444,0.946267628818689,-3.95028655721432) q[12];
u3(2.57109867061972,1.42594954063936,-0.188696383544610) q[9];
u3(1.72820029192127,-0.284106388594686,-2.76675172677707) q[0];
cx q[0],q[9];
u1(0.108909819900918) q[9];
u3(-1.29607003248247,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.46391222853806,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.00276949876730,-2.37150646848864,1.28068248919938) q[9];
u3(1.55252790966024,-2.32684848656563,-3.04407131387730) q[0];
u3(1.22665743817713,-0.332435887879099,1.45302422534448) q[5];
u3(1.06774325431426,-2.37926538509390,-1.85213070667661) q[4];
cx q[4],q[5];
u1(-0.436316396149039) q[5];
u3(0.134575569953363,0.0,0.0) q[4];
cx q[5],q[4];
u3(4.11970150681197,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.81405419144904,1.17722468689746,0.995649546448128) q[5];
u3(1.43734736639569,-1.01335896028847,2.81523242243873) q[4];
u3(2.32396945692732,-1.11810410401794,2.42117307891732) q[8];
u3(2.77178035164772,-0.735899445451226,-0.242793517770363) q[9];
cx q[9],q[8];
u1(3.31457274392268) q[8];
u3(-1.64834847945386,0.0,0.0) q[9];
cx q[8],q[9];
u3(2.75775270665719,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.740112323551881,-0.122338620348570,-2.70006821352482) q[8];
u3(0.731318685609827,1.44331621643340,4.00053550335143) q[9];
u3(2.09128579513646,-0.912434628966698,-1.22352170254396) q[2];
u3(1.53442235410199,-2.91397574856163,0.0398922606854448) q[7];
cx q[7],q[2];
u1(4.57253992253274) q[2];
u3(-3.70786788915674,0.0,0.0) q[7];
cx q[2],q[7];
u3(-0.477098729691711,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.27173331966435,0.748729848659519,-3.71400445507622) q[2];
u3(1.67170571938553,-1.64963327398068,3.35373672850589) q[7];
u3(1.54645064516263,-1.13616410617096,1.29585171529325) q[11];
u3(2.43572465094542,-1.69220415711055,-2.01209618891241) q[6];
cx q[6],q[11];
u1(-0.302149332328213) q[11];
u3(-2.08730261965195,0.0,0.0) q[6];
cx q[11],q[6];
u3(1.36698933910519,0.0,0.0) q[6];
cx q[6],q[11];
u3(2.47930428533403,1.25627337445863,2.09756755135853) q[11];
u3(1.28142384587686,-1.36492171257779,1.90920410023562) q[6];
u3(2.36048777567578,0.117756568918080,2.16472435137484) q[10];
u3(1.54830057867511,-0.392571513251187,-1.59307637893966) q[0];
cx q[0],q[10];
u1(2.89454975353020) q[10];
u3(-1.60031776098759,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.22042259073453,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.880413464663037,-0.321673748554260,-0.356901872281729) q[10];
u3(0.333972976613275,0.647038026071176,-4.30680932578057) q[0];
u3(2.05449607674607,-0.895849087457213,1.30304632865013) q[1];
u3(1.91840089125506,-2.02030022223892,-1.69850129599276) q[12];
cx q[12],q[1];
u1(-0.332814980813202) q[1];
u3(0.991516606010811,0.0,0.0) q[12];
cx q[1],q[12];
u3(3.20125518531356,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.36163493267134,1.91542497272789,-4.02167856441046) q[1];
u3(0.590369150589635,0.189476955887997,-1.93501123414353) q[12];
u3(1.66303636608883,0.389466890161817,-1.88149536295367) q[12];
u3(1.08125048827340,-3.74859722434198,1.68505652830014) q[10];
cx q[10],q[12];
u1(2.76781028839381) q[12];
u3(-1.66907624398322,0.0,0.0) q[10];
cx q[12],q[10];
u3(3.19421696796203,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.02629447419111,1.46501876437720,2.39990895881091) q[12];
u3(1.92097717573963,2.33358523805469,3.66763205311222) q[10];
u3(2.24728949788786,-0.393675291588123,-1.68552361389607) q[4];
u3(1.37546896092597,0.478073875879977,-3.94608564669056) q[11];
cx q[11],q[4];
u1(1.49112152654502) q[4];
u3(-0.633355820437386,0.0,0.0) q[11];
cx q[4],q[11];
u3(3.31591707501782,0.0,0.0) q[11];
cx q[11],q[4];
u3(0.287840037113122,-2.40344099208470,-0.0822032521183926) q[4];
u3(0.995455365354594,-6.03177944086118,-0.213585440900422) q[11];
u3(1.50171983305139,-0.917141410337497,-0.564874814594627) q[0];
u3(1.47867486385866,-2.58038427170861,-1.03552222827452) q[2];
cx q[2],q[0];
u1(0.899102835327232) q[0];
u3(-3.54346669198506,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.66863454919160,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.89028444279975,2.67242163981299,-2.58454246780316) q[0];
u3(2.38198390889850,-0.546838834440963,3.11738384032045) q[2];
u3(0.281568021677865,-1.48165614849122,2.52778365624768) q[7];
u3(0.485176710181324,-1.42092619225683,0.373434131852604) q[8];
cx q[8],q[7];
u1(0.331824410068939) q[7];
u3(-0.779348271686525,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.16614534659892,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.62593654996187,1.65888193335600,-2.71095757524011) q[7];
u3(1.48164427826663,5.83533528275254,0.433457613565211) q[8];
u3(0.801091431366964,1.64637271664396,-2.62281285997005) q[6];
u3(1.45423636941253,1.66934326035709,-3.97141042712093) q[3];
cx q[3],q[6];
u1(0.0176858246500808) q[6];
u3(-0.575737733104419,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.64051478169689,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.78505508095119,-1.10805691752211,-2.11437783580455) q[6];
u3(1.41715201793105,-4.31856623257308,1.17300753353562) q[3];
u3(0.950284784724912,1.56972487766667,-0.855679695845666) q[9];
u3(0.358689659646870,-2.11089583756576,0.836394170935807) q[1];
cx q[1],q[9];
u1(1.45103195535062) q[9];
u3(-0.817209494674092,0.0,0.0) q[1];
cx q[9],q[1];
u3(-0.259225778689839,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.01810497512869,2.92423527714597,-3.02999148071135) q[9];
u3(0.642558475108191,-5.66330465752184,-0.303306832973276) q[1];
u3(1.81047963131741,-3.79247414309505,2.44648823267960) q[3];
u3(1.41225380226164,3.60768442617768,-2.30835401717714) q[5];
cx q[5],q[3];
u1(0.660393500083436) q[3];
u3(0.148004944728924,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.90196586074646,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.16419445333303,2.40896635685030,-1.02697860109632) q[3];
u3(2.60002278293939,-3.13614024005195,2.58241641814515) q[5];
u3(1.17418975751007,1.34189923827725,-3.60990055749286) q[12];
u3(2.10330040904846,-1.96568123373988,3.60709902146255) q[11];
cx q[11],q[12];
u1(2.94951908961582) q[12];
u3(-2.02457160434531,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.01539873297710,0.0,0.0) q[11];
cx q[11],q[12];
u3(2.57502449972107,0.884243860714335,-0.0335097272150356) q[12];
u3(1.71440479880353,-2.52366187954767,0.221032584327413) q[11];
u3(0.526900184567938,1.54461681081508,-2.69238149445668) q[7];
u3(0.456938922101728,1.43649497959093,-3.43352905567832) q[6];
cx q[6],q[7];
u1(1.40379866411488) q[7];
u3(0.469229787147835,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.85511291283758,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.99745731599507,-0.429848334541490,-1.51666884275113) q[7];
u3(0.800123275969272,-1.70917354200733,1.91692106132867) q[6];
u3(0.495622156917053,-2.91139574437448,-0.170647215160229) q[0];
u3(2.11033872664651,-3.62142295940523,-0.405708575873998) q[4];
cx q[4],q[0];
u1(-0.109082684508546) q[0];
u3(-2.28252848225677,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.32360805167458,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.993846028034856,2.27271367686798,1.21704834968081) q[0];
u3(0.744232430985821,3.72651349463775,-1.44099695147871) q[4];
u3(1.23398840441414,1.44572508390685,-2.82356146255191) q[2];
u3(2.39961502119123,-3.64533740069491,2.58147934002701) q[9];
cx q[9],q[2];
u1(0.896466577109423) q[2];
u3(-1.18115314820130,0.0,0.0) q[9];
cx q[2],q[9];
u3(3.33176559355094,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.26243879258058,0.878060029454199,0.254568192985015) q[2];
u3(1.12744134760553,-3.66483037270424,1.62620741073059) q[9];
u3(0.761423581037872,1.49179022182000,0.327374600614941) q[8];
u3(1.05527348781006,-0.0643921567183845,-2.96686744330399) q[10];
cx q[10],q[8];
u1(1.66142246211610) q[8];
u3(-3.16845114558591,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.317026100942468,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.33127407733908,2.21557859265527,2.09605634299844) q[8];
u3(1.42934715460309,1.34171121660404,1.20553537344778) q[10];
u3(1.03594114032651,0.0100882363551731,-1.21275635945493) q[9];
u3(1.74507948336220,1.74433189927269,-3.74993427082036) q[4];
cx q[4],q[9];
u1(0.442541335188504) q[9];
u3(-1.36526098138649,0.0,0.0) q[4];
cx q[9],q[4];
u3(2.40437490668202,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.266273825372437,-0.355554429523859,1.80155307844509) q[9];
u3(2.34671520699016,2.34978226994708,3.26622195775482) q[4];
u3(2.17682787801643,-3.02963490615912,2.75007647626298) q[10];
u3(1.13459327160782,2.98252428825004,-1.61834613934183) q[6];
cx q[6],q[10];
u1(2.23891019593344) q[10];
u3(-2.56828016595680,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.55893517054858,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.94733268016641,0.459493176216876,0.864098898339008) q[10];
u3(1.57134823326197,5.55201507775418,0.435698857649748) q[6];
u3(2.73612268223983,-3.49692169591514,2.34075256341952) q[12];
u3(0.863369435456100,3.02067604361213,-1.30178974321010) q[11];
cx q[11],q[12];
u1(0.234376503894639) q[12];
u3(-1.62699130083351,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.38667899255996,0.0,0.0) q[11];
cx q[11],q[12];
u3(2.00815146421502,3.52497323090314,0.471946280590047) q[12];
u3(1.84522769048538,0.632863500350048,2.18697668522202) q[11];
u3(1.41205748490008,1.19882720734825,-2.62140436904870) q[0];
u3(2.24770612053581,-2.38171506167911,3.14119865270251) q[3];
cx q[3],q[0];
u1(3.07041946420560) q[0];
u3(-1.61214636131962,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.787581953058070,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.41701250070364,-2.30131430793517,1.12727600417434) q[0];
u3(1.54314027757544,-1.48368787225356,2.00117996833512) q[3];
u3(1.65098006964953,2.44655328003695,-1.73873587590887) q[5];
u3(2.36624521604983,1.38717519348185,-0.302500862239025) q[2];
cx q[2],q[5];
u1(0.872182145720309) q[5];
u3(-1.46343197589370,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.81474623241292,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.917744347458721,-2.70796003915191,1.86762902425230) q[5];
u3(1.94850718482882,1.70253861503337,-4.52777214725115) q[2];
u3(0.998940138902071,-0.191748212652949,0.460044744733742) q[8];
u3(0.350328204847734,-1.83885560373803,1.42164079038390) q[1];
cx q[1],q[8];
u1(1.11040517529148) q[8];
u3(0.0178605424300651,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.59455466106954,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.22267567872685,-3.79505508780432,2.08467646493967) q[8];
u3(0.296147308682214,1.62553053994108,2.45341841335749) q[1];
u3(1.06085122984726,1.98890616810121,0.762191397357303) q[6];
u3(2.15727554607432,0.0424856742293502,-3.20690421474791) q[12];
cx q[12],q[6];
u1(2.02155666545704) q[6];
u3(-2.52031238082159,0.0,0.0) q[12];
cx q[6],q[12];
u3(0.272166421916909,0.0,0.0) q[12];
cx q[12],q[6];
u3(2.30986292709882,-4.01708149037918,1.79829618061813) q[6];
u3(1.82661542831884,-2.92520544561950,0.569055070440027) q[12];
u3(2.30864012211321,0.827914439711501,1.25021646758331) q[2];
u3(1.77752349037199,-1.49113866583314,-1.87236682710792) q[10];
cx q[10],q[2];
u1(1.70468598006232) q[2];
u3(-1.96216785559901,0.0,0.0) q[10];
cx q[2],q[10];
u3(-0.108074059420958,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.16944553564574,1.34014121053853,0.403071864652224) q[2];
u3(1.59246125163944,-0.189650375160779,-3.25884174167666) q[10];
u3(1.20182058291144,0.155829149915495,0.839012007988254) q[3];
u3(1.25358473907199,-0.444894144470823,-1.70562847130217) q[5];
cx q[5],q[3];
u1(2.33920552157294) q[3];
u3(-0.0840509476112348,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.49162421405760,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.994370079476238,0.905923070328381,-0.0753040377304178) q[3];
u3(1.87003996100430,0.0469729378193067,3.01427006200760) q[5];
u3(1.55818688824128,-0.753652888582746,1.88660177161743) q[11];
u3(1.77225120508591,-1.81560386101184,-0.489283987495920) q[7];
cx q[7],q[11];
u1(2.15014418799397) q[11];
u3(-3.02604650508546,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.378595564570207,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.27937430540182,2.14454498532910,-3.03456856306553) q[11];
u3(1.58813075776217,-2.07354964989281,4.04366395396432) q[7];
u3(1.71582676491721,-1.64323689703224,0.628244818713799) q[1];
u3(1.77193288422208,-4.03736601684083,-0.572339122994455) q[0];
cx q[0],q[1];
u1(3.50162466131784) q[1];
u3(-1.36527874238713,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.31500074014888,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.45044452678484,0.826911938182137,-1.24359538219230) q[1];
u3(2.12244296204229,-4.60697346118301,-0.828307127072795) q[0];
u3(1.73490764029940,-0.187670788032569,-1.82334997808259) q[4];
u3(2.31878489512630,1.15481867489561,-3.82073413802149) q[9];
cx q[9],q[4];
u1(-0.389252185346650) q[4];
u3(1.19914843097976,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.25625723762643,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.48390974781523,1.88626499322222,-3.03598030627545) q[4];
u3(1.68156190198111,3.68785772635929,-2.35142758620357) q[9];
u3(1.39238792281299,2.57339491022780,-0.110659411507511) q[0];
u3(1.98978378009744,-0.00607978782759799,-4.75186039206954) q[7];
cx q[7],q[0];
u1(0.660906172649663) q[0];
u3(-3.28403150816752,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.80077762104223,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.227984058200741,-2.32240179543591,3.88741930576022) q[0];
u3(1.50538507944701,0.676081582444656,-2.98705614685003) q[7];
u3(2.93085377391616,1.32552908136128,-3.43493763921300) q[11];
u3(1.61248441124607,-2.88676602911409,2.63070000371581) q[12];
cx q[12],q[11];
u1(1.12437477890970) q[11];
u3(-3.73784397795288,0.0,0.0) q[12];
cx q[11],q[12];
u3(1.51101569908903,0.0,0.0) q[12];
cx q[12],q[11];
u3(2.05928014756159,1.07020060060909,-1.02178400923816) q[11];
u3(1.95687490848386,-0.664032706916467,1.19983890131396) q[12];
u3(1.30004193207624,-0.137601229434338,0.631861343788183) q[1];
u3(1.33506284215860,-1.60099815882291,-1.56172729167785) q[2];
cx q[2],q[1];
u1(1.48562466099656) q[1];
u3(-0.753573741101336,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.0894203258874560,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.04899207790656,-2.08561984776149,-0.702733374046061) q[1];
u3(1.58246492291000,-1.69149036134234,-2.41817702836326) q[2];
u3(0.176366141912496,-1.30315258067211,1.55961996651323) q[4];
u3(0.652209446796211,0.854710023957799,-1.40366740754994) q[9];
cx q[9],q[4];
u1(0.709532319789359) q[4];
u3(-0.227154233427233,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.57477696189928,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.80239277964961,1.60549457637802,-2.59575871695585) q[4];
u3(0.165986859315985,-1.02122302659217,-1.58591022359475) q[9];
u3(1.67139054287759,0.445897968877416,1.85513724606299) q[3];
u3(2.08655548001620,-2.79364920838229,-1.93742423631601) q[10];
cx q[10],q[3];
u1(3.19869390030514) q[3];
u3(-1.84361683155632,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.78279380534198,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.37397681576981,-3.62208992067372,1.47721858477095) q[3];
u3(2.45184497182587,-0.717657379641036,-3.01214280607473) q[10];
u3(2.08998998999806,1.83285511804083,-3.36547807812792) q[8];
u3(0.927176638773957,-2.38193278670368,2.96437141078810) q[5];
cx q[5],q[8];
u1(-0.409995657281249) q[8];
u3(-2.34511319104732,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.83540594195287,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.47479160558965,2.79990929772607,-3.22086235866672) q[8];
u3(0.471571626495912,4.42406120329618,-1.50634309562620) q[5];
u3(1.08737423656344,3.68417119095696,-1.13396708398307) q[10];
u3(0.831273245817131,0.913298299365936,-1.13054475714079) q[1];
cx q[1],q[10];
u1(2.36835852705132) q[10];
u3(-2.86081514721088,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.650278838660518,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.80445474252282,1.88609216353451,-0.321433697054679) q[10];
u3(0.711379650195475,0.104364453150493,-2.93169497779079) q[1];
u3(2.71561325841887,-0.416452212120091,-1.67791204716111) q[7];
u3(1.67126089252734,-3.61622496951494,1.20441478558621) q[3];
cx q[3],q[7];
u1(2.12054938764306) q[7];
u3(-2.54903904922890,0.0,0.0) q[3];
cx q[7],q[3];
u3(-0.0422284223557428,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.78251876993633,1.04971188013943,-4.09776608251355) q[7];
u3(0.801952845889915,-3.54353322342654,-1.28752311204703) q[3];
u3(2.52377335149304,2.05710295887261,-3.70777395127558) q[8];
u3(0.733427568828568,3.58418592764174,-2.59986624034793) q[2];
cx q[2],q[8];
u1(-0.233164914403848) q[8];
u3(-2.26268474076679,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.34930856825511,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.48092549783947,-1.49420209915347,2.22481150304254) q[8];
u3(1.86369285957688,0.400131041021767,-5.85068803419809) q[2];
u3(2.34747574776358,1.23716581028273,0.00887911909486860) q[11];
u3(2.07568074417185,-0.186884738280666,-3.41710184865314) q[9];
cx q[9],q[11];
u1(1.53944239526605) q[11];
u3(-2.85195874718157,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.751211317393087,0.0,0.0) q[9];
cx q[9],q[11];
u3(0.490357126648061,-2.31219124423840,3.44109568420055) q[11];
u3(1.02153234730227,-0.244779329728497,1.83818652172006) q[9];
u3(2.51288212940076,0.342999681073979,1.09312939547763) q[0];
u3(1.08730843353621,-2.85312631501101,-1.85910103676160) q[4];
cx q[4],q[0];
u1(1.08733069619777) q[0];
u3(-1.50722714446937,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.87179005003104,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.15948930711164,-3.64119935419756,0.804971998489938) q[0];
u3(2.39269454828004,-1.45692142557325,0.453846390782978) q[4];
u3(1.35355482994820,-0.334552306127608,-1.61095651730810) q[12];
u3(1.22700224912689,0.548732693492213,-4.45525134734613) q[5];
cx q[5],q[12];
u1(2.33940104139343) q[12];
u3(-1.44929006355055,0.0,0.0) q[5];
cx q[12],q[5];
u3(3.43711783623254,0.0,0.0) q[5];
cx q[5],q[12];
u3(0.0398863282865395,-3.79064204749633,1.99668023691781) q[12];
u3(2.42491359265634,1.78469120418917,4.24298573853187) q[5];
u3(1.26302283308256,0.535062038722100,0.461501453630488) q[12];
u3(0.909691027656403,-1.57594514993326,-1.62838515224030) q[1];
cx q[1],q[12];
u1(3.10874641737966) q[12];
u3(-1.13413872138041,0.0,0.0) q[1];
cx q[12],q[1];
u3(2.66504988034102,0.0,0.0) q[1];
cx q[1],q[12];
u3(2.74353920520234,1.76868272052527,-0.897003613071036) q[12];
u3(1.43664644972997,-1.94852277730810,1.42023170727267) q[1];
u3(2.87861236014638,4.26627355518110,-1.27950996469249) q[8];
u3(0.902612282936768,2.60919606619594,-0.580317562821016) q[6];
cx q[6],q[8];
u1(3.91492504402005) q[8];
u3(-3.56561387591655,0.0,0.0) q[6];
cx q[8],q[6];
u3(-0.627613807215559,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.89761093454578,1.84801698652631,-1.72503657833323) q[8];
u3(1.15277702294519,-0.991323019360249,-3.25479270403011) q[6];
u3(1.89686012556887,-1.06645875805689,-0.454497779652948) q[7];
u3(1.69976322325263,-2.17809855561274,0.877182505133797) q[9];
cx q[9],q[7];
u1(3.04269677125451) q[7];
u3(-2.33596529038726,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.55855883493696,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.738356291391983,0.544175741276350,-1.73261661936203) q[7];
u3(1.46985483415104,3.79701755556146,1.18962673795148) q[9];
u3(0.692923512088337,3.03159663321686,-2.80114951565441) q[5];
u3(0.666833443171728,0.875990749896397,-1.71452234076007) q[2];
cx q[2],q[5];
u1(2.02929081637498) q[5];
u3(-2.40941706718800,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.41373712722917,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.91209716204221,-1.83615943593937,3.71950537889659) q[5];
u3(1.80571341868038,1.13101318160478,0.953217856641280) q[2];
u3(1.67357105675182,-0.886610847923031,-0.0758349275636697) q[11];
u3(0.918682478062795,-3.16336271498046,-0.544230474741095) q[3];
cx q[3],q[11];
u1(0.416008064428673) q[11];
u3(-1.32613735052983,0.0,0.0) q[3];
cx q[11],q[3];
u3(2.35722886805150,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.04434528222905,2.64492582794624,-3.12490561945549) q[11];
u3(0.987316855692751,-4.48767616875936,-0.798917746549737) q[3];
u3(1.67718925886157,0.837407511528639,0.283070531821615) q[10];
u3(2.19668661074220,-0.480962819120621,-3.37633101207076) q[4];
cx q[4],q[10];
u1(3.26121418512366) q[10];
u3(-1.30925948616629,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.15994218316233,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.980775998030398,-0.543567738070814,1.99298596317070) q[10];
u3(1.61805800909657,-1.48924357247702,-2.10508212335678) q[4];
u3(1.28630857007191,2.75689937544136,-1.71931330799563) q[5];
u3(2.14133437912358,2.35137529980980,-0.561652163078292) q[8];
cx q[8],q[5];
u1(2.34072498714700) q[5];
u3(-1.68774841661191,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.03574463583644,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.549455130178076,0.776367936888068,0.535182425952900) q[5];
u3(1.21156312450863,-0.526633547007582,0.176983980393866) q[8];
u3(1.80957556148754,-1.10213632497851,0.739351222957856) q[3];
u3(2.09109011962554,-4.20495180728060,0.0445380043740662) q[4];
cx q[4],q[3];
u1(1.67401781334148) q[3];
u3(0.688548500581672,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.05270708424300,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.63399299077535,-1.65483363745929,4.48404991988300) q[3];
u3(2.28199829553468,1.08749680012042,3.23465883729111) q[4];
u3(1.24629669767903,-0.509505721763955,-1.62827766031940) q[10];
u3(1.88807272238807,0.700129486745292,-5.06396006632961) q[12];
cx q[12],q[10];
u1(3.00693924154965) q[10];
u3(-2.08935751628041,0.0,0.0) q[12];
cx q[10],q[12];
u3(1.67315391816852,0.0,0.0) q[12];
cx q[12],q[10];
u3(1.50463922329021,1.61453551035080,-1.05256883471285) q[10];
u3(2.19536082357169,-2.08910498403120,3.56684608322409) q[12];
u3(1.65396322797891,0.124203321332206,-2.20225178057142) q[0];
u3(2.06116055466275,1.93493180629295,-4.00595351493998) q[1];
cx q[1],q[0];
u1(-1.03039784280430) q[0];
u3(0.103993341322544,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.59667251415126,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.10373742851305,-3.74407826267134,1.59676805555120) q[0];
u3(1.33771602603190,3.35401269105349,-2.80996694179090) q[1];
u3(1.58892596812155,-1.56657878919491,3.90421984249718) q[7];
u3(0.126246493593175,-0.483637204128628,0.496334520181824) q[9];
cx q[9],q[7];
u1(0.611842559241169) q[7];
u3(-1.58226858860464,0.0,0.0) q[9];
cx q[7],q[9];
u3(-0.389245468905540,0.0,0.0) q[9];
cx q[9],q[7];
u3(3.03904524598214,-1.88008426844016,3.71416725640355) q[7];
u3(2.14879332443929,-0.0249216442018771,1.07521240867119) q[9];
u3(2.47269010226070,0.232853787799911,-3.29325240336264) q[11];
u3(2.17439191266734,4.36016352996990,-1.14803434380771) q[2];
cx q[2],q[11];
u1(1.58694317493364) q[11];
u3(0.206348545841731,0.0,0.0) q[2];
cx q[11],q[2];
u3(0.740151415582037,0.0,0.0) q[2];
cx q[2],q[11];
u3(2.71438946130888,1.84960605379906,-2.59831999666802) q[11];
u3(1.50812979289151,1.62625512005787,-0.743902544128320) q[2];
u3(1.26085328726634,1.08397052547719,-3.78471791783195) q[9];
u3(1.62716433175134,2.95986790303247,-2.58067843073710) q[4];
cx q[4],q[9];
u1(2.65356879847837) q[9];
u3(0.128076379034487,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.48937948971762,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.715234406460075,2.43309483617862,-0.762847565187919) q[9];
u3(1.77962579459568,-0.327342976159646,2.85937962136678) q[4];
u3(2.15399012178502,3.14635031561846,-0.985902768779689) q[6];
u3(2.51005361382498,4.06976666376858,-0.848421495796896) q[5];
cx q[5],q[6];
u1(0.296834803386824) q[6];
u3(-1.18826335881223,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.46831925102358,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.05019086943720,-2.02947799057965,1.78693035916545) q[6];
u3(1.09862210555828,-1.64140873856104,-3.24186753964956) q[5];
u3(2.35360059133594,0.575550427051746,1.99590918468915) q[1];
u3(1.98747342955074,-2.02318709571642,-2.80790349598197) q[3];
cx q[3],q[1];
u1(1.83952951835696) q[1];
u3(-2.15220285792894,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.427948628147091,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.75181104019982,-0.630562242263086,-0.345172731021200) q[1];
u3(2.11859687335546,0.0842624170978201,1.49618839964281) q[3];
u3(1.77434617378601,0.323790237099630,-0.748095902376619) q[0];
u3(1.49021364005506,0.864577851396917,-4.03295021521394) q[2];
cx q[2],q[0];
u1(3.27191283605888) q[0];
u3(-1.18131977768675,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.45692174811639,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.43661276270666,0.630905851053632,-2.46693651303334) q[0];
u3(1.41376291894740,1.89113595286257,-2.87832911503739) q[2];
u3(0.253695339808316,-1.06574241763646,-1.20095986187017) q[8];
u3(1.62987139897515,-4.76868620209431,1.36138810099837) q[7];
cx q[7],q[8];
u1(0.715125554021204) q[8];
u3(-1.17855967099202,0.0,0.0) q[7];
cx q[8],q[7];
u3(-0.0238014999200564,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.06918930532008,-1.09884028288378,1.74358494004837) q[8];
u3(1.26252086912996,-1.24544994414002,1.13663803181421) q[7];
u3(2.09596032743865,-0.639502056209948,1.95767142546258) q[12];
u3(2.19390044032708,-2.23103210637790,-1.86075482599068) q[10];
cx q[10],q[12];
u1(2.37222310488233) q[12];
u3(-1.57859658247972,0.0,0.0) q[10];
cx q[12],q[10];
u3(0.425615301424791,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.58833125867765,3.44428679968540,0.680840534937220) q[12];
u3(0.666069287401727,-1.16153679540456,1.66225941302519) q[10];
u3(2.43636500559432,1.15233167650416,-2.74230059921914) q[3];
u3(1.42783246577407,-2.11627632075849,2.27146099861263) q[4];
cx q[4],q[3];
u1(1.75894010660387) q[3];
u3(-2.90903714889227,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.396394361225573,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.07145318664086,-0.709762300919045,-2.26114885396744) q[3];
u3(1.52596413406857,1.23084662601716,2.84215811711177) q[4];
u3(0.703845496384966,1.90037166301394,-3.21780791231471) q[7];
u3(2.15429727382900,3.21591094564007,-3.04974346599005) q[6];
cx q[6],q[7];
u1(1.40163538133339) q[7];
u3(-0.728494852686977,0.0,0.0) q[6];
cx q[7],q[6];
u3(-0.209884126801771,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.443800315216757,2.02351996235313,-1.88806390875108) q[7];
u3(1.48151871705084,-0.134591471466957,2.82702944726411) q[6];
u3(2.05777271169553,0.138490527446940,1.11590288463367) q[0];
u3(2.22901108856221,-0.805151832263891,-1.22170375302820) q[5];
cx q[5],q[0];
u1(0.700744890698263) q[0];
u3(-1.56773366204729,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.97142467508167,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.903958402404350,0.558295425074294,-1.57310001812723) q[0];
u3(2.14392442525818,1.17130579508582,-3.69472497625564) q[5];
u3(2.65375119572240,3.92165850024520,-0.919983527498501) q[8];
u3(1.68854248050234,1.98824728176235,-0.774706838087781) q[11];
cx q[11],q[8];
u1(0.791304916489008) q[8];
u3(-3.42373766271912,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.56434097405150,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.29250439767337,1.09874469702404,0.882198242728919) q[8];
u3(1.19959740170318,-0.926522104671491,3.88264897544807) q[11];
u3(2.44054107299552,-1.68358572745177,3.78160393697525) q[1];
u3(0.388624975365501,3.35985049363329,-1.49590996235232) q[9];
cx q[9],q[1];
u1(1.85524874273964) q[1];
u3(-2.42827919789169,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.46399181048850,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.09345006407053,-0.0784530674642795,-1.71804387746174) q[1];
u3(2.72266289329850,-1.92708341116945,1.90932969755649) q[9];
u3(1.15211263312802,0.500846460140026,-1.13205154651058) q[10];
u3(0.855293582918595,-4.80764769293272,1.12010717290915) q[12];
cx q[12],q[10];
u1(0.699657427902474) q[10];
u3(-0.104145898987746,0.0,0.0) q[12];
cx q[10],q[12];
u3(2.22240744894757,0.0,0.0) q[12];
cx q[12],q[10];
u3(0.507465639225579,-1.46941920711345,0.537103292832593) q[10];
u3(2.30651766271653,-0.774832895584215,-2.78927010184080) q[12];
u3(1.90912830822233,-0.979313758374825,0.292879024093221) q[1];
u3(1.09885361229349,-2.40957149045347,-1.24155172404257) q[6];
cx q[6],q[1];
u1(1.49424811126157) q[1];
u3(-0.146449507991116,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.66954310121739,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.62698401509186,0.635640768599113,0.175286100273106) q[1];
u3(0.911151878427322,-2.40454850758830,0.749268704636993) q[6];
u3(0.256470963168037,-0.999941111026216,0.0989743765770406) q[4];
u3(0.743867723459669,-0.894525728154480,-1.18360721114442) q[10];
cx q[10],q[4];
u1(0.515439586679588) q[4];
u3(-1.72452727524591,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.89055700341264,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.12901696841180,2.28522118117524,-1.10887811047219) q[4];
u3(1.41695062035769,3.78749563478181,0.316627820455164) q[10];
u3(0.320412712979096,-2.37979209630530,3.33231516892619) q[11];
u3(0.432402784105834,0.0868767052803963,-1.01663880101803) q[5];
cx q[5],q[11];
u1(2.50209695689214) q[11];
u3(-1.40484017275163,0.0,0.0) q[5];
cx q[11],q[5];
u3(0.177242918173269,0.0,0.0) q[5];
cx q[5],q[11];
u3(2.11807278835574,-1.31794413313956,-0.985865063988464) q[11];
u3(0.966877484234623,-0.831571225232916,3.62548332053068) q[5];
u3(3.03768205703495,-0.437217184599507,-0.136393149731855) q[12];
u3(1.71403299613191,-2.27684502495867,-2.09449366994923) q[0];
cx q[0],q[12];
u1(3.43145104989682) q[12];
u3(-1.59149998043548,0.0,0.0) q[0];
cx q[12],q[0];
u3(2.37787783967241,0.0,0.0) q[0];
cx q[0],q[12];
u3(0.784473012300356,-1.65894010288767,3.22668709366677) q[12];
u3(2.11482212043098,0.925433992455821,-2.49954549583450) q[0];
u3(0.927372476549375,0.232826115822933,1.01865687802385) q[9];
u3(1.92252496396939,-0.0853199434339381,-2.06222878284175) q[2];
cx q[2],q[9];
u1(1.43167776184827) q[9];
u3(-3.32207554001191,0.0,0.0) q[2];
cx q[9],q[2];
u3(2.11385919089006,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.44994494641585,-1.06487688901585,1.89406380473325) q[9];
u3(1.90499368037672,-1.71867321494328,0.0767384178687182) q[2];
u3(3.05202381973212,0.898194766419817,-1.83245622413928) q[8];
u3(1.68480474525421,0.138803763581991,-3.28328684877363) q[7];
cx q[7],q[8];
u1(3.38562744834001) q[8];
u3(-3.82566824153077,0.0,0.0) q[7];
cx q[8],q[7];
u3(-0.987662559067174,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.848336145595686,2.66387501432375,-2.35842415739674) q[8];
u3(0.809270789819318,-2.06284761943668,3.39109420946635) q[7];
u3(0.661017232540484,0.977523859970212,2.03116769679543) q[8];
u3(1.61188272611561,-1.70426005435345,-1.53726464382509) q[6];
cx q[6],q[8];
u1(2.30583097404353) q[8];
u3(-3.12151984207338,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.62686670603409,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.54605476905618,2.76786714918041,0.186830564777873) q[8];
u3(0.711947080671141,2.85390050171830,-2.62955493543454) q[6];
u3(1.68605510145853,1.41085195003996,-3.65803466244681) q[10];
u3(0.789858757516159,1.97143773931014,-1.85233191120153) q[2];
cx q[2],q[10];
u1(1.55827523902112) q[10];
u3(-0.0714102936319814,0.0,0.0) q[2];
cx q[10],q[2];
u3(0.410909563126635,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.94243684192460,0.526088874851637,2.99132588103309) q[10];
u3(0.290422033077743,-3.20056505915728,0.440755160977515) q[2];
u3(0.445798534845202,-3.41465444033684,2.67286281811157) q[9];
u3(0.844298266158773,0.395405472689212,-1.67420030927684) q[3];
cx q[3],q[9];
u1(1.18012989054076) q[9];
u3(-3.13655272974240,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.61626708093093,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.17646578271215,0.861424740169701,-1.94563039843043) q[9];
u3(1.03152275597841,2.98662300697559,1.85377625500260) q[3];
u3(0.901568301612885,-0.533305238898799,0.944805660032054) q[1];
u3(0.866111812448151,-0.876786714726321,-1.94643480417499) q[4];
cx q[4],q[1];
u1(3.02186022134552) q[1];
u3(-1.04655971436745,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.82697253063319,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.47733395419093,-4.41800672802248,1.51075218294162) q[1];
u3(1.51876216451951,-1.42372390607901,2.90330613485095) q[4];
u3(2.13154445470953,3.74623788580834,-1.26009463719554) q[12];
u3(1.52882012764228,2.56710998720266,-0.365751383408868) q[11];
cx q[11],q[12];
u1(3.44691937588734) q[12];
u3(-1.48282180824506,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.89695767495246,0.0,0.0) q[11];
cx q[11],q[12];
u3(2.13220796020580,-0.305407587014787,-1.91691441990155) q[12];
u3(0.851882874406841,-1.54273254579309,0.721974430477417) q[11];
u3(0.393071729282406,2.10473286664340,-2.08373381351206) q[7];
u3(0.541937424650846,1.41420433981081,-2.36577634413739) q[0];
cx q[0],q[7];
u1(1.65924841675440) q[7];
u3(-2.98177082160987,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.732127500195520,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.65056239686283,-1.09983911045075,1.17711650918422) q[7];
u3(1.08353379555732,2.81140571687880,-2.30143821892074) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12];
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