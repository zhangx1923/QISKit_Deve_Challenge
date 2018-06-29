OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.55024176345278,-0.0111328544358973,-2.03330161797162) q[7];
u3(0.349292569080017,1.06152299120930,-3.76138134500726) q[6];
cx q[6],q[7];
u1(2.15736709853268) q[7];
u3(-2.40040801327451,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.45777341133839,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.93008759254517,1.89288183688628,-3.88705665800680) q[7];
u3(1.74257397810509,-1.39713284067490,-1.60514974030831) q[6];
u3(0.121513789528781,2.24439351496796,-2.28209208107822) q[2];
u3(1.19126215877261,-3.31645112890847,1.17617788106468) q[5];
cx q[5],q[2];
u1(0.960322570729395) q[2];
u3(-1.60579109831875,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.75553592802276,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.19884475879576,0.661801097072490,-5.00754466195099) q[2];
u3(2.31293716644107,-2.07225190372873,-0.315177755808376) q[5];
u3(0.894158537324198,0.533007710194768,-0.920034154110340) q[8];
u3(0.854325953132625,-1.28388017549725,0.213423114360138) q[3];
cx q[3],q[8];
u1(0.809021104166361) q[8];
u3(-1.15962379857462,0.0,0.0) q[3];
cx q[8],q[3];
u3(-0.0228326155191845,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.58555245658641,-0.343588048720551,0.985210998846316) q[8];
u3(1.80325185374803,3.16176711984043,2.90663884585405) q[3];
u3(2.46419992778062,0.591664585982511,-0.838290841913108) q[0];
u3(2.02988412704770,-4.15146684976602,0.716870024640796) q[4];
cx q[4],q[0];
u1(-0.256169711780239) q[0];
u3(-2.30413830904346,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.49837935313616,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.96815074423621,4.02214047856729,-1.47585981758427) q[0];
u3(1.65527381569959,2.46597894514772,-2.79194553935497) q[4];
u3(2.69613788784552,2.39292551887096,-0.824966466716313) q[9];
u3(2.07176518522489,1.22781753740964,-2.01016631448606) q[1];
cx q[1],q[9];
u1(-0.196859731473173) q[9];
u3(-1.62673474840363,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.973750001412768,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.880489982169615,2.37897561869954,-2.61029735509391) q[9];
u3(2.21998872418169,-1.68016425065381,-1.64548270689767) q[1];
u3(2.54660797600368,1.53636712531289,-2.21902165932294) q[9];
u3(1.88322771474396,2.46070209533138,-2.53749882386289) q[1];
cx q[1],q[9];
u1(0.452042584887559) q[9];
u3(-1.34127239451131,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.92741964684709,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.661870272745253,-2.94381196837126,2.67040591334974) q[9];
u3(1.51839360245671,0.0554611631568875,1.07018039786498) q[1];
u3(2.62886196185637,3.30444490847574,-0.341238267504191) q[0];
u3(2.79840622517968,4.72252338309958,0.0742547614210531) q[7];
cx q[7],q[0];
u1(0.674851022719710) q[0];
u3(-0.318221410043499,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.08797863889573,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.12594297484616,1.84758795220519,-0.827023337613529) q[0];
u3(0.737578328546087,2.83397255545883,0.0718397243875324) q[7];
u3(0.878138293005708,-0.507068910221206,-1.60613584179220) q[8];
u3(0.812417528287767,-4.07097693765937,1.19967690318715) q[2];
cx q[2],q[8];
u1(2.53823202465632) q[8];
u3(-1.62642469346105,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.05370824210235,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.25472379152269,-2.92327535499219,2.40100838644281) q[8];
u3(1.05858082173105,-0.434887865011827,0.827332311968602) q[2];
u3(2.28181500821832,0.974207421633986,-3.78221258370012) q[4];
u3(1.42988901453393,2.75685119289810,-2.72122299588051) q[6];
cx q[6],q[4];
u1(3.12287843947310) q[4];
u3(-2.24387881371398,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.518704543975987,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.503546503625248,-2.31276731993579,2.40713517184952) q[4];
u3(1.51094159644157,-2.00218364490525,-1.68219921978007) q[6];
u3(2.55296206531200,0.943679809003667,-1.30581816761323) q[3];
u3(1.43748218113217,-4.33160164402508,1.85818777714018) q[5];
cx q[5],q[3];
u1(1.36263595474269) q[3];
u3(-3.16938190627951,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.64315638171980,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.56592605196605,2.63778492360515,0.0364343333800567) q[3];
u3(0.650818972655690,-0.942286753176065,-3.24346075079384) q[5];
u3(0.892879644166044,-0.179382325514830,-1.29986801390912) q[7];
u3(1.54748719016997,-3.54945184407560,2.62622034713324) q[6];
cx q[6],q[7];
u1(-0.183193666305184) q[7];
u3(-1.50181379007973,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.21680271653578,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.88831704714040,1.53252494480015,0.173012395208579) q[7];
u3(1.71970016741560,0.874141278262333,5.30471941018070) q[6];
u3(1.12061878570929,-0.777230955128477,0.794390968150754) q[1];
u3(0.462353224351320,2.15478729234434,-2.66158641378246) q[4];
cx q[4],q[1];
u1(0.965953765348476) q[1];
u3(-0.187026945629740,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.83201808666122,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.818571120337071,3.23794897995633,-0.337312211600071) q[1];
u3(2.10098151850451,3.27866695739937,-0.632547424366571) q[4];
u3(1.76205328574729,0.490603873574187,1.82934912465112) q[2];
u3(1.92959712059056,-0.569782376200435,-0.523023676586525) q[0];
cx q[0],q[2];
u1(2.24117173699506) q[2];
u3(0.201909922312856,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.56281163758353,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.367516293172185,0.750807219913840,-3.63989259721139) q[2];
u3(0.411242581635292,-0.697357346079948,-2.81656353067572) q[0];
u3(1.46013066222737,-0.800529671166827,1.05603624994104) q[5];
u3(2.31700681850315,-2.01867693505162,-2.40615380296957) q[9];
cx q[9],q[5];
u1(3.47129222304688) q[5];
u3(-1.44129766957854,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.57333725172407,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.929286560168379,1.24025434502258,0.482022160422204) q[5];
u3(1.07702217951006,-0.975145582902174,-4.04302951032939) q[9];
u3(1.43329290880813,3.72133454141729,-2.32662530735769) q[3];
u3(2.15267014517819,0.737522781836735,-2.89717732904441) q[8];
cx q[8],q[3];
u1(0.460193352566034) q[3];
u3(-0.928803309109925,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.04271213695784,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.23727460347489,1.13264748142510,-1.82145450483597) q[3];
u3(0.876966179921255,-2.47614860644157,-3.40066460720259) q[8];
u3(2.01254549543221,1.34497570049675,-4.05818023002732) q[4];
u3(1.27036625069996,3.76462691238897,-2.50650100118885) q[0];
cx q[0],q[4];
u1(2.75599968709258) q[4];
u3(-2.34896343789185,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.69380004463554,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.21212051697175,1.03831605354624,-2.26619877601441) q[4];
u3(1.35759471928044,1.32489267015938,-2.09668009252813) q[0];
u3(0.746089393010312,-1.74599349539472,-1.02565829556471) q[7];
u3(0.977757221765130,-3.09209608051448,-0.297746175523764) q[3];
cx q[3],q[7];
u1(-0.298386724020550) q[7];
u3(-2.59954683681325,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.38349487058712,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.28167282664082,2.29423702677287,-2.50180770142535) q[7];
u3(1.30110146296829,2.90609446709309,1.01599339982979) q[3];
u3(1.96911568505386,-0.380936347105786,0.367847760258125) q[1];
u3(1.59784834685689,-2.69585377874488,-0.487628434801031) q[9];
cx q[9],q[1];
u1(2.33729925307789) q[1];
u3(-1.79606930037476,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.964058330347386,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.590564757204536,1.24409878683478,1.38613296882805) q[1];
u3(0.815196143203161,2.11121616298990,2.61775435331139) q[9];
u3(2.43439685584478,1.51577785240433,-3.97533624249826) q[2];
u3(1.53829131009711,-2.05905912039045,4.03415140602236) q[6];
cx q[6],q[2];
u1(2.52371036078530) q[2];
u3(-1.49963811684152,0.0,0.0) q[6];
cx q[2],q[6];
u3(3.44723072900077,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.879541135979481,2.16256472419547,0.883807486359892) q[2];
u3(1.98479872029079,0.495360002456906,-5.03684738030796) q[6];
u3(1.98042166439719,0.555071785082430,-1.93792350369157) q[5];
u3(2.15276040948583,2.68185833999749,-3.09136354403875) q[8];
cx q[8],q[5];
u1(-0.164040865847806) q[5];
u3(1.00792586759088,0.0,0.0) q[8];
cx q[5],q[8];
u3(3.35914848085363,0.0,0.0) q[8];
cx q[8],q[5];
u3(0.585648367115348,-0.216370192881721,-2.79648413818160) q[5];
u3(1.44639852857211,1.42580487398937,-2.84048220220568) q[8];
u3(1.92698245580514,-0.398851951858878,-2.70000926637714) q[2];
u3(2.74608384041039,-0.0261249397445296,-5.00084145319237) q[5];
cx q[5],q[2];
u1(0.262421303535226) q[2];
u3(-1.52605920653406,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.11719462942291,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.14904126537473,-1.82574319882172,1.36788201220067) q[2];
u3(1.52530504527683,-1.42349872441033,-0.628442670125751) q[5];
u3(2.15044940771456,2.74557234846422,-1.76909310336816) q[3];
u3(2.37109118536532,1.41912926962420,-2.27529341336155) q[0];
cx q[0],q[3];
u1(3.84162762513672) q[3];
u3(-3.72377053017798,0.0,0.0) q[0];
cx q[3],q[0];
u3(-1.07627405200454,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.28421824558434,-2.55965541172493,1.08445863375565) q[3];
u3(0.328861477006846,0.777309056553776,-4.88189863220566) q[0];
u3(0.761541668013163,0.568891578112788,-2.64846858854433) q[9];
u3(1.30713108449591,1.34912223220254,-4.81768715217618) q[4];
cx q[4],q[9];
u1(1.96726599087465) q[9];
u3(0.283943292260177,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.14187917332550,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.63519331857352,-1.68152717701729,2.03662629652939) q[9];
u3(0.132683260288961,2.71028230658723,2.25980724428202) q[4];
u3(0.388036485868234,-3.21999309324267,2.85482672451728) q[6];
u3(0.658002933108950,-3.48592048227451,2.74221532955339) q[8];
cx q[8],q[6];
u1(1.21374782331835) q[6];
u3(-1.63456460309528,0.0,0.0) q[8];
cx q[6],q[8];
u3(3.49031629415922,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.85898825288559,0.790094858478499,1.41332377226618) q[6];
u3(1.34017439006827,4.22998092991758,1.67122173905823) q[8];
u3(0.655145468390417,1.03122315463705,-2.97569939721326) q[7];
u3(1.31850274319963,3.74410571841054,-2.12348462859405) q[1];
cx q[1],q[7];
u1(2.97708600321866) q[7];
u3(-1.41877061881174,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.14924368054376,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.795932345604482,2.69997394340020,0.313742164020399) q[7];
u3(1.09905864051556,2.30366227436699,0.393944908354556) q[1];
u3(1.23091961342260,1.52178389290722,-0.0877813231527894) q[7];
u3(1.18720711025546,0.164760956629691,-2.10059423184525) q[5];
cx q[5],q[7];
u1(2.33500127995265) q[7];
u3(-1.97118982953240,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.40583736928899,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.85821422186009,-1.54249810062504,-1.49255257097701) q[7];
u3(1.53340683293979,-4.87917209538032,0.922176859732010) q[5];
u3(1.15823778995805,1.09992337522246,-3.51426200815268) q[6];
u3(2.05231577610422,2.55250007991100,-2.50162746450013) q[4];
cx q[4],q[6];
u1(-0.482100108867106) q[6];
u3(1.37087762057413,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.44533634135273,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.17628606774257,4.06893274432020,-1.58478427004364) q[6];
u3(1.27198468760422,-1.42783244377465,-2.32235994331187) q[4];
u3(1.31108819654384,-1.85877738968554,1.65177308891364) q[1];
u3(1.63097242079529,-1.64647986650609,2.72031335633346) q[0];
cx q[0],q[1];
u1(2.72723869712638) q[1];
u3(-1.87499402943829,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.913899462160434,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.946250571216900,-3.20032745764487,2.42306415022641) q[1];
u3(1.04003783381840,1.03045927589804,5.14045228015880) q[0];
u3(2.65622142976285,0.308252603886742,-0.977824455532793) q[3];
u3(1.55773128665460,-3.81183120790132,1.30528099249575) q[8];
cx q[8],q[3];
u1(0.623401481391765) q[3];
u3(-1.56850223174743,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.14316456922455,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.844162317405775,-3.44094983755747,0.850335643506283) q[3];
u3(1.46044197194633,3.18654010830965,-1.57311547969784) q[8];
u3(1.39052461364298,-0.845977155984890,-0.0165650266291892) q[2];
u3(0.444415060539772,-1.97495695759102,1.12350940531126) q[9];
cx q[9],q[2];
u1(1.37562906424641) q[2];
u3(-0.603487538323776,0.0,0.0) q[9];
cx q[2],q[9];
u3(-0.237586387590029,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.69284375136225,2.36713025435308,-0.765813736453521) q[2];
u3(1.31657698711297,1.92239309218539,2.07210432262755) q[9];
u3(0.346038355163176,1.25840530913847,-0.667565771630780) q[9];
u3(0.315166440234381,-1.70391033967361,1.02508078757852) q[6];
cx q[6],q[9];
u1(0.629936934283136) q[9];
u3(-1.90888767678132,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.231026227989948,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.191583937360290,-0.588793610560452,2.17216133388342) q[9];
u3(2.07996704593504,0.866324807811902,-2.18765941310644) q[6];
u3(0.403704081197863,-1.80444738272917,1.42824012178141) q[3];
u3(0.915158893885960,2.01381774377344,-3.10318912983968) q[0];
cx q[0],q[3];
u1(0.0843139838735645) q[3];
u3(-2.47452226549670,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.845572156235350,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.06549833468329,-4.40071700701581,1.50256645178125) q[3];
u3(2.63900437961181,-3.80784667568186,-0.801194894259383) q[0];
u3(1.41605148998046,2.37642498878404,-1.13222796917218) q[1];
u3(1.46082596324212,1.55451373368591,-0.456150878400389) q[5];
cx q[5],q[1];
u1(1.85047434566208) q[1];
u3(-2.77086998710531,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.616433769839178,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.39720084286017,-3.45080851210023,0.967667839984020) q[1];
u3(1.70338675573005,-1.36996974790884,2.09201386989315) q[5];
u3(2.56950680929766,2.50352364864703,-1.23956247192604) q[7];
u3(2.33328989608657,-0.163004451574956,-5.47943568935774) q[2];
cx q[2],q[7];
u1(2.95333649791934) q[7];
u3(-2.70074318083820,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.99028186917617,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.28981076012914,-1.71235272762259,3.23290785482295) q[7];
u3(2.04784674638806,-1.45996138036916,-3.65376867790814) q[2];
u3(2.98509895104281,0.424096393593699,1.20147819719550) q[4];
u3(1.32035776269134,-3.18296653314244,-1.49285639521586) q[8];
cx q[8],q[4];
u1(3.05263506800523) q[4];
u3(-1.73366316675177,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.696009070783117,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.82618592171993,1.02117876810861,-2.77935078341937) q[4];
u3(0.234858638028964,0.775567468870631,3.03737266936409) q[8];
u3(0.733646523577749,-2.98927399249176,3.13448047803109) q[6];
u3(1.37438424332234,0.0607555613169187,-1.57331260099666) q[9];
cx q[9],q[6];
u1(-0.171952808563527) q[6];
u3(-1.37677886683348,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.17191575167335,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.58668443379101,0.830567977285118,-3.91866562574758) q[6];
u3(1.77229354368124,-3.33865009525499,0.448396008439329) q[9];
u3(2.12971957340957,0.327151948415026,0.750151253919665) q[5];
u3(0.942051929498809,-2.45754556325796,-2.17891911618196) q[0];
cx q[0],q[5];
u1(2.05158200562414) q[5];
u3(-1.70835427530512,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.94671220664410,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.70605912957619,-1.02496140035119,-1.04439856307209) q[5];
u3(2.77081576214005,0.750120397178103,-3.92917298561660) q[0];
u3(1.28697358291768,-1.22568720997734,1.98929944266368) q[4];
u3(2.01630774655508,-1.65327189182952,-1.87623966953812) q[1];
cx q[1],q[4];
u1(0.713467644851414) q[4];
u3(0.0524101580704972,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.18176776652836,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.06699638603274,0.802857990646837,0.539296171450694) q[4];
u3(2.18307043296987,1.01819078901425,-1.28557629210691) q[1];
u3(1.13058693601413,0.437118176640173,2.17997437838978) q[2];
u3(1.81195335532948,-1.62346749243913,-1.07179097948640) q[3];
cx q[3],q[2];
u1(1.79770629402340) q[2];
u3(-3.17698255260769,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.77356427441400,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.66959658796446,-1.66016184501016,-0.529142230861426) q[2];
u3(1.83778636943481,3.09730672560565,0.277834338150405) q[3];
u3(1.98566669449799,-0.434968906242928,1.12234363809649) q[8];
u3(1.99123638124118,-2.79126846648467,-1.59665946013727) q[7];
cx q[7],q[8];
u1(1.59804815883049) q[8];
u3(-2.60461806895737,0.0,0.0) q[7];
cx q[8],q[7];
u3(-0.207388114136021,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.0494143264238296,-2.30308630095055,0.656643273729137) q[8];
u3(1.37596144144913,-2.88641511484182,0.483294180649822) q[7];
u3(2.32347447656848,1.42689263454761,-2.48920519738943) q[4];
u3(1.26936519349587,-2.58542504547118,2.01947347998082) q[6];
cx q[6],q[4];
u1(1.78845860142255) q[4];
u3(-2.10513936518598,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.129292638778252,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.25057591230680,3.19994472561025,-2.33215052520873) q[4];
u3(1.25116836633448,-3.41802177903873,2.22498675297979) q[6];
u3(0.527510592047717,-0.140061037111940,-2.34304770480612) q[7];
u3(2.00855069101766,1.09841727150775,-5.02247311967316) q[0];
cx q[0],q[7];
u1(3.20547558737030) q[7];
u3(-1.16180672766236,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.26969342468453,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.47245042469727,3.61223369935606,0.250526387717254) q[7];
u3(2.71583049595541,3.54712479091454,-0.967015800673453) q[0];
u3(0.670431428727312,-0.833686414903095,1.11673436256356) q[8];
u3(0.0563893566688816,0.0897661209419048,-0.513326386408687) q[5];
cx q[5],q[8];
u1(0.772122127771956) q[8];
u3(-0.0585104814868949,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.75627100534545,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.42587531827041,1.78739523263226,-1.50791842599951) q[8];
u3(2.01589329425502,-1.37960071785896,-0.179555606477424) q[5];
u3(2.02818279274618,-0.915352895624181,0.00596370477976160) q[3];
u3(2.24651470462260,-1.91269920474763,0.943679149025580) q[9];
cx q[9],q[3];
u1(-0.467388661240140) q[3];
u3(0.0394379204139346,0.0,0.0) q[9];
cx q[3],q[9];
u3(3.96739488441892,0.0,0.0) q[9];
cx q[9],q[3];
u3(2.09768980611888,2.34517001532375,-2.81339159196100) q[3];
u3(3.03330683038684,2.55034137023803,0.536805363622581) q[9];
u3(1.92971231639003,-2.91006313311134,2.31131035894090) q[2];
u3(2.56082377851931,-2.92181861573472,1.14180895286681) q[1];
cx q[1],q[2];
u1(1.80281042908037) q[2];
u3(-2.47569087158763,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.226269261302481,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.79325278544171,0.114936887156349,-0.699723397497548) q[2];
u3(2.32080532840934,-0.758550197017432,4.91040607103212) q[1];
u3(2.69435111336882,1.70547210024244,-1.42528178506501) q[1];
u3(2.42709575484700,1.72481495795820,-3.98065722586686) q[4];
cx q[4],q[1];
u1(0.686827999375792) q[1];
u3(-0.214548388456791,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.98843457767050,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.49653790130077,3.27128707119943,-0.0752143111915287) q[1];
u3(1.82343159326322,1.30349785548313,-0.651399953379009) q[4];
u3(2.17153075752026,-0.00347348424253635,2.28089453072734) q[3];
u3(2.54810416638320,-2.08386301980107,-1.68836908800795) q[7];
cx q[7],q[3];
u1(0.653594030524143) q[3];
u3(-0.111935052144215,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.35900801247384,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.16867393551958,0.333678593318884,-0.224999922224830) q[3];
u3(1.20586170894237,-0.350295572607808,0.484960580755716) q[7];
u3(1.40964862989412,0.0729913511949585,-2.04270417139230) q[6];
u3(0.915358534211202,1.83524040200561,-3.92863255397632) q[9];
cx q[9],q[6];
u1(1.27160886559247) q[6];
u3(-0.0508623919563895,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.17656463002340,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.998713768259263,-0.916946971804858,2.86584851821137) q[6];
u3(0.979693974010174,0.102300350103725,-2.58437794130272) q[9];
u3(1.37501054936781,3.48398453178228,-1.92958729036182) q[5];
u3(0.590211926501189,1.43444446737724,-0.364256004881191) q[0];
cx q[0],q[5];
u1(1.33468802189739) q[5];
u3(-0.887883191760736,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.87506673963366,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.71044607857814,-2.69338847879691,2.78802227034084) q[5];
u3(1.35753304870661,2.29608426487457,2.97491448864299) q[0];
u3(1.37079188093080,3.36160313136500,-1.38423940676823) q[8];
u3(1.27917562406140,1.74341432799363,-2.19988246906342) q[2];
cx q[2],q[8];
u1(1.84737143128289) q[8];
u3(-0.00424903765372076,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.331752887772745,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.52142045789663,2.21892988538419,-1.21766458198323) q[8];
u3(0.955829878329704,-0.728174446011565,1.15511803687085) q[2];
u3(1.28647124396634,-0.463413604860871,1.85838263210922) q[0];
u3(1.38960215770487,-2.32130195238628,-2.53245041727370) q[3];
cx q[3],q[0];
u1(1.59420442770951) q[0];
u3(-3.07563189011580,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.55178559400666,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.14492228365137,-1.31014807015883,3.39734856514113) q[0];
u3(1.24957176697511,-0.764165075920642,-4.85218602126802) q[3];
u3(2.65621803523852,-1.01176200275692,-1.43812687025971) q[5];
u3(0.820602171326279,-3.46603134330103,-1.06764975682431) q[4];
cx q[4],q[5];
u1(1.59999679408978) q[5];
u3(-2.38036146260029,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.312758607338921,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.36030996348688,3.09275850951345,-2.44599075251406) q[5];
u3(1.35490253329161,0.438044495060320,-4.82855645009856) q[4];
u3(1.61028718960809,1.60550737968095,-3.70792093069165) q[2];
u3(1.02068602357870,2.62353221607357,-2.27445489562162) q[6];
cx q[6],q[2];
u1(0.984536675865105) q[2];
u3(-1.42359640504219,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.0170128807090113,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.559291426960843,-2.92436898860081,0.295927883237549) q[2];
u3(1.30797117346704,-0.199918925776774,-1.16411325051613) q[6];
u3(2.41287415075397,-1.22115457667514,0.657611758617641) q[7];
u3(2.06702646821760,-1.87380446121076,-0.263859498272495) q[9];
cx q[9],q[7];
u1(1.44514549369904) q[7];
u3(0.544050679234651,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.63710771591080,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.80824622747020,0.368353102189499,1.30760506873982) q[7];
u3(0.529796087127667,-3.46857399488059,1.51746462811748) q[9];
u3(1.57587770522008,-2.36544677986691,-0.300842355127275) q[8];
u3(2.04771166173623,-3.31670022716761,0.191912309740913) q[1];
cx q[1],q[8];
u1(0.636914131448556) q[8];
u3(-1.29716835948355,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.73486756726228,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.73287678177115,0.160961112533074,1.57996197618682) q[8];
u3(1.34254628611074,-3.93351561703629,0.613555323905908) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
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
