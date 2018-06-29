OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(0.945980781603226,0.568776598659022,-1.36995591570411) q[3];
u3(0.443279709297155,-1.68726388793712,-0.144315620069711) q[6];
cx q[6],q[3];
u1(0.574306602929879) q[3];
u3(-1.48605817020027,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.107313694123787,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.02483794111100,0.586869539124740,-2.83171623353803) q[3];
u3(0.799342076515446,2.17448780302621,-1.53023088120183) q[6];
u3(0.765961529316478,1.90786451644183,-2.61543820095903) q[5];
u3(1.57171652739703,-2.88346808373753,3.05245799073327) q[4];
cx q[4],q[5];
u1(2.16204701411739) q[5];
u3(-2.85859028447890,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.856038799053701,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.130798239939511,2.62930865083970,-1.37078384684487) q[5];
u3(2.06351926323721,1.95135735569849,-1.15069968178816) q[4];
u3(1.59399624926459,0.570875274343426,-1.08317307032510) q[2];
u3(0.991405195107199,0.00443043721072889,-3.93311903648040) q[7];
cx q[7],q[2];
u1(0.573231779490672) q[2];
u3(-1.16514329970039,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.06538000146260,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.960011408729822,2.26694730356779,-3.13605135461473) q[2];
u3(1.06965611330458,3.66180042639752,-1.17999616270394) q[7];
u3(1.23212846078971,-3.91560241613571,2.18942610946344) q[0];
u3(1.86309812623577,2.74775257817176,-3.33147179911751) q[1];
cx q[1],q[0];
u1(1.18398200632792) q[0];
u3(-0.637981011520792,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.72473599681705,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.00300673584800,-2.11546461119701,4.14776871129944) q[0];
u3(2.56962163294424,-2.32053587324729,0.248873412022699) q[1];
u3(1.35590137116107,0.554880014047826,-3.41810208332121) q[7];
u3(0.537521111179584,3.04394961172435,-2.60823153143529) q[2];
cx q[2],q[7];
u1(0.202815264255167) q[7];
u3(-1.65913938268582,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.75505417411995,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.98092829292001,1.10002148433619,-0.450306948091931) q[7];
u3(2.26841886692881,-4.87919478543262,0.361265155204338) q[2];
u3(0.570953641087509,1.84903426152635,-1.81440448031922) q[6];
u3(0.485775456735690,0.616259609303228,-2.92445247793138) q[0];
cx q[0],q[6];
u1(2.22224907848754) q[6];
u3(-0.0442861216147983,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.59768625846484,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.07605606670257,1.87866677249514,-2.79966926275063) q[6];
u3(2.05274854920432,-2.09282602146916,-3.64249246657710) q[0];
u3(2.48629814346032,2.77163980579806,0.344390722574508) q[3];
u3(2.45059682527071,2.49629632318609,-2.20675348256817) q[5];
cx q[5],q[3];
u1(-0.244441426418488) q[3];
u3(-1.64190716536058,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.790108859128233,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.164788802091438,-0.872723453224360,5.11631805214677) q[3];
u3(2.07393350384395,-1.61270500429981,-2.63154481643669) q[5];
u3(1.40979038985670,-0.212236697033780,-1.74200077518469) q[1];
u3(2.50532822389672,-4.11976191804466,2.02363517989570) q[4];
cx q[4],q[1];
u1(1.22586303894333) q[1];
u3(-3.46070230040031,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.45039810930214,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.66315025173948,3.22107822434485,-1.50580025246068) q[1];
u3(1.07854431442948,2.22128863642232,0.530092363136173) q[4];
u3(1.06088052224951,-1.00074875861658,-1.71739640993694) q[6];
u3(2.52548001096857,0.380507834249903,-5.45757369408385) q[2];
cx q[2],q[6];
u1(3.09135573622920) q[6];
u3(-2.01452179815466,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.740333425697416,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.33694054629479,1.01008645494875,-3.75994559079779) q[6];
u3(1.73791418462652,-0.938142325482236,-4.64176968619652) q[2];
u3(1.62632489893194,1.64640357926728,-0.101179694407367) q[3];
u3(0.882647948565736,-0.00190651801789032,-4.58021434349944) q[7];
cx q[7],q[3];
u1(-0.280729075561747) q[3];
u3(-1.74731414824872,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.13402641998821,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.54274079565657,1.80261918425234,-3.40495422082029) q[3];
u3(0.194568588806224,4.84103438671318,-0.0546658967813718) q[7];
u3(1.25041545458211,2.09454134597513,-3.10698923767098) q[4];
u3(2.50321746094716,-1.78519506692230,3.53506131781907) q[1];
cx q[1],q[4];
u1(2.33186651477552) q[4];
u3(-2.93393818538309,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.36802939250502,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.33907660642712,0.526269987886300,-2.08893162185152) q[4];
u3(1.17080470380243,2.12276103358420,-0.760189011252191) q[1];
u3(0.682680468111151,-2.42020082016494,3.33094844111020) q[0];
u3(0.905611454669107,0.657653653268709,-2.49912532484943) q[5];
cx q[5],q[0];
u1(1.48791417715658) q[0];
u3(-3.25612498782021,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.48710263819663,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.43292155792789,-0.902642835212890,0.516770150153346) q[0];
u3(0.509820338727377,-2.68253304572668,2.69786460502065) q[5];
u3(2.58063047241099,-2.90198401488379,0.00132739981895180) q[6];
u3(3.02367076012656,-4.15724045200531,-2.10342839215387) q[7];
cx q[7],q[6];
u1(-0.488455392874168) q[6];
u3(-2.02020763791973,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.01016705146807,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.73574472292439,3.86818279702164,-0.0837786705880061) q[6];
u3(2.03697471368018,4.72232100880145,-1.50261472586206) q[7];
u3(2.08524190621997,0.203901845249749,1.16562346551709) q[5];
u3(2.09889272941472,-1.34932308670513,-1.19801768832459) q[4];
cx q[4],q[5];
u1(2.57197281115995) q[5];
u3(-3.38209441797762,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.57044055625385,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.24176433330466,4.31591274976069,-1.35372991515342) q[5];
u3(2.22158337128560,1.90547129218966,3.23071417357716) q[4];
u3(1.76285888172895,0.333566862188853,1.72443614802077) q[2];
u3(1.77352384432247,-1.06236963992135,-0.926670207573673) q[1];
cx q[1],q[2];
u1(0.0792563002180862) q[2];
u3(-1.76940668725510,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.650691899333633,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.67319676403566,0.343486749890651,0.783088887876025) q[2];
u3(2.21198961302131,-0.151875570098350,-5.16314005365289) q[1];
u3(1.24794189299905,-2.02029877690862,0.389674969429208) q[0];
u3(1.30428703799619,-4.18828444309395,0.406371759150838) q[3];
cx q[3],q[0];
u1(2.71080973964995) q[0];
u3(-1.54997538646290,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.165228503258674,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.70583719398838,0.535476320618346,3.30504305948051) q[0];
u3(0.949742961418930,1.09728884461747,2.28957914149075) q[3];
u3(2.11452427186080,-0.697506413149879,1.37450877422140) q[4];
u3(2.09058050583250,-2.10081380976290,-0.914644781596706) q[7];
cx q[7],q[4];
u1(1.10254317190322) q[4];
u3(-3.48068170613555,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.63739383078458,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.399965474600224,1.95955775059177,-2.32359877243566) q[4];
u3(1.58771875774666,-4.56435283062496,-1.00376501069845) q[7];
u3(1.26712568978021,-2.70299694445325,2.17218829883450) q[6];
u3(0.801186534933180,0.749257643477919,-3.04281852720770) q[0];
cx q[0],q[6];
u1(-0.366445622961049) q[6];
u3(-1.50626874150209,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.90639809927368,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.790184118071880,1.49423729221525,-0.132185358167855) q[6];
u3(0.692921478959410,3.06015708514492,2.18507515471747) q[0];
u3(1.85064428362672,3.56706312534453,-1.73658012155195) q[1];
u3(0.261755498277513,1.00345268318275,0.727481550534651) q[5];
cx q[5],q[1];
u1(1.72993082534251) q[1];
u3(0.544877100699476,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.02489163258744,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.66973679863611,-0.334981400100507,1.92075253687206) q[1];
u3(0.644672403952031,-0.596899829278212,5.58354719067870) q[5];
u3(2.16646693244846,3.91971335243359,-0.841486980320264) q[2];
u3(2.00114231001337,2.91848185334524,0.451476859564034) q[3];
cx q[3],q[2];
u1(1.82279417992419) q[2];
u3(0.281949370452261,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.664898093079852,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.06044527764831,-1.31905607113912,-0.0202275384984367) q[2];
u3(2.04049401929774,-2.25897228675910,-1.05033585894833) q[3];
u3(2.79146864313486,0.195701005944116,-1.34359241823886) q[6];
u3(1.72490902752538,0.887209703132306,-3.54889174654550) q[2];
cx q[2],q[6];
u1(0.558684252606956) q[6];
u3(-1.32222454641660,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.83553880612888,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.12472738003250,2.88912690366245,-0.882827061513417) q[6];
u3(1.33895822060384,-2.46056262241288,2.53738671981479) q[2];
u3(0.572826690377013,2.75480540951871,-2.58191049656649) q[3];
u3(0.975658262752319,1.34432274571850,-2.07996997526377) q[1];
cx q[1],q[3];
u1(1.84359719826217) q[3];
u3(-0.255786604873299,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.60094770986745,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.40214539386215,1.59497387349105,-3.10468491523376) q[3];
u3(1.97962951794882,-2.15713402995364,-2.44477050635511) q[1];
u3(0.971643846478819,-2.42977487490319,0.721918585114970) q[5];
u3(0.542870948825611,-2.01498069741550,0.142063584459764) q[4];
cx q[4],q[5];
u1(1.79252857153370) q[5];
u3(-2.43229023507935,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.93705429346851,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.46442305204696,-3.64936089981051,1.24589142477094) q[5];
u3(1.00109810110539,1.64664031334811,1.89786001117420) q[4];
u3(1.83210028495365,3.04561473651432,-1.85624396852495) q[7];
u3(0.841730065823898,2.51106493801183,-1.91098200946042) q[0];
cx q[0],q[7];
u1(1.48432752552412) q[7];
u3(-0.364163912193227,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.33893235881139,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.55653247770714,0.116161792703927,1.04505398047870) q[7];
u3(2.05366881539031,-1.33935906294837,-2.61574132134910) q[0];
u3(2.39355396844211,2.67634111525427,-3.54342743722282) q[0];
u3(0.269053101218699,-0.0435160183440177,2.95875881878622) q[5];
cx q[5],q[0];
u1(1.18741942980833) q[0];
u3(-0.453730452851398,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.24126714702866,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.960284478240778,-2.41594345480565,2.73276594854066) q[0];
u3(1.77646078736488,0.872086283668057,-4.89645221183824) q[5];
u3(1.53615131012759,2.10481716922788,-2.57746570059138) q[6];
u3(1.31360909545857,1.51679029811350,-2.18514244557703) q[3];
cx q[3],q[6];
u1(1.78655647646223) q[6];
u3(-3.16654114567102,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.865530919747317,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.73156613751278,1.30250813073986,-4.21886739498695) q[6];
u3(0.930045615668316,4.37419568848516,0.601489275059661) q[3];
u3(0.737567511021452,-0.219651938578072,-0.931459172581073) q[2];
u3(0.964703764975979,-3.71233019508835,1.36486768049795) q[4];
cx q[4],q[2];
u1(1.13983127066107) q[2];
u3(-0.277240403332567,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.58570759554823,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.524543024967918,-3.43519270574663,2.07096336330011) q[2];
u3(0.504402004444783,3.44836976483774,2.41616964416058) q[4];
u3(1.54043612535993,-0.688042741185766,-1.16797940263149) q[1];
u3(1.19149581064104,-4.55395718712997,1.16641227826717) q[7];
cx q[7],q[1];
u1(2.08187226857913) q[1];
u3(-2.57774209595529,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.47557214980952,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.73043942883386,-1.34727998952084,-1.00550926640601) q[1];
u3(0.998510302785625,-3.37765355053051,1.69588210069263) q[7];
u3(0.164657710546121,0.339358726578554,-1.58941781771619) q[2];
u3(0.828736754244267,-3.39850771813852,1.73252956636298) q[3];
cx q[3],q[2];
u1(0.301020249846070) q[2];
u3(-1.32189439407205,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.06201931318360,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.68259033648879,0.355323366546419,0.247761755739096) q[2];
u3(0.992789915401833,-4.51711491754990,-0.948050730297610) q[3];
u3(0.554464686039269,2.98982202248713,-3.13553867855987) q[4];
u3(1.28170333073662,0.715307739821885,-1.51960140669737) q[1];
cx q[1],q[4];
u1(3.22221454237227) q[4];
u3(-0.585778706609443,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.88732195861765,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.09808462641661,-0.831293593328274,0.584087840717997) q[4];
u3(1.53673657342413,-3.79652119684420,1.76303827837899) q[1];
u3(1.90938330126979,2.20754221857159,-3.47844331256136) q[0];
u3(2.20461080608650,-2.56627988165466,3.23778136164119) q[6];
cx q[6],q[0];
u1(2.42966539220847) q[0];
u3(-2.81203890480673,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.12545462568054,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.14833655652975,0.953505791227470,-0.0261791261937000) q[0];
u3(1.53143950667271,-5.63369620726795,0.288697279137928) q[6];
u3(0.321240179941964,3.31960034848571,-2.36654974325172) q[5];
u3(1.09152017545921,-3.42201715805344,1.97816930902872) q[7];
cx q[7],q[5];
u1(1.06982414714226) q[5];
u3(-2.89598154680589,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.46687228087357,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.279517526855797,0.755227761910982,0.129770011619812) q[5];
u3(0.571483862267621,1.36095016255796,0.962752918369837) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];