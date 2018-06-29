OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(2.24800776511031,0.810922566182099,2.15211321579449) q[4];
u3(2.25443348029124,-1.51646345437916,-1.19376433044769) q[2];
cx q[2],q[4];
u1(2.95602905275803) q[4];
u3(-1.37863716968885,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.92021213318205,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.37472391558757,-0.304635443531529,-2.74054918423765) q[4];
u3(1.44467601290255,-5.05688868328120,0.299704428665252) q[2];
u3(1.74039810159806,2.08834196048948,-2.99713758277732) q[5];
u3(2.32398697223164,3.76611183720570,-2.51484474229858) q[0];
cx q[0],q[5];
u1(2.41658531035267) q[5];
u3(-1.72590059773692,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.324398707854054,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.44775297360796,-0.194678786477492,1.65250936687498) q[5];
u3(1.41247277919415,4.27559093467467,-1.11529414348183) q[0];
u3(2.00498573924239,0.729893347599176,0.515132384308235) q[3];
u3(2.13760256009076,-1.85012847545207,-1.20415825220930) q[1];
cx q[1],q[3];
u1(1.55074789401636) q[3];
u3(-3.23155830991469,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.309708304653603,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.09398707961909,0.438745840758101,-0.0677855271207244) q[3];
u3(0.780754261872143,2.67947561266724,-1.62381581075142) q[1];
u3(0.705800645970273,-0.195234219384022,-1.07231420963113) q[3];
u3(0.812636728401654,-0.295197472447513,-1.15654962854614) q[5];
cx q[5],q[3];
u1(-0.0769380954195182) q[3];
u3(-1.81152557356142,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.12828625132014,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.78875878654578,1.31091706694397,2.33617842132888) q[3];
u3(1.70261038802075,-0.954086061787135,-1.08842299732265) q[5];
u3(0.598612088268449,2.27387192745814,-1.27976409673554) q[1];
u3(1.20240824833327,0.612958646738800,-2.58389986648757) q[2];
cx q[2],q[1];
u1(2.78283206884199) q[1];
u3(-2.01320168709515,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.809505735440535,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.31814853808584,1.03667855023561,0.300449284094800) q[1];
u3(2.58249276872336,2.93382369206965,3.18707866802845) q[2];
u3(0.515503727979038,-3.49897647246023,2.02128031375106) q[4];
u3(1.49358118994181,2.93603912144167,-2.73714602381391) q[0];
cx q[0],q[4];
u1(3.34054013185084) q[4];
u3(-1.55147460313050,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.867718448261518,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.90125926910633,0.148523056608277,0.838300106239008) q[4];
u3(1.77225777974983,-1.27415351115761,-1.18091313810453) q[0];
u3(2.54014033382314,0.467646307120193,-2.41784023482895) q[1];
u3(2.26183807597196,1.60530575567126,-3.82078326692525) q[5];
cx q[5],q[1];
u1(1.19340189891053) q[1];
u3(-3.14466262612168,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.10051032754927,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.72705657632372,-0.745047483219131,0.534560294005916) q[1];
u3(1.18929481334882,-1.71703051239051,-0.672233886230985) q[5];
u3(1.94090296683439,-0.194283335333825,2.62917750356840) q[2];
u3(2.41651914136863,-2.11738138266833,-1.90261636804165) q[4];
cx q[4],q[2];
u1(1.21933074384180) q[2];
u3(0.00902645485072462,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.36570420922842,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.87319437144435,-1.85277188550550,2.77770202263745) q[2];
u3(1.19739688221866,-0.119080534507272,-2.82943226678249) q[4];
u3(2.67676399937389,0.564253016160588,-2.55620653335384) q[0];
u3(2.78129851225683,-0.524223633901276,-4.35939793387819) q[3];
cx q[3],q[0];
u1(2.65718696460937) q[0];
u3(-1.58339266831746,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.23326017555490,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.883760673896671,-0.336215219097792,4.11005964224213) q[0];
u3(2.53538464209021,-1.43959668961510,1.74401261179437) q[3];
u3(1.69737666704722,0.691208090620054,-0.851268880764758) q[0];
u3(2.78215700538894,-3.62035462566446,1.80839689316679) q[1];
cx q[1],q[0];
u1(2.96393439263436) q[0];
u3(-2.41243199297641,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.17581116317716,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.72921105273772,0.457026604210417,2.72524263744421) q[0];
u3(0.414137400742671,1.42238374194934,0.616728357597147) q[1];
u3(2.08079174435553,-0.904813944703466,1.07451508727912) q[2];
u3(1.59935482302116,-1.85241253703540,-0.561510889210135) q[3];
cx q[3],q[2];
u1(2.25019568449146) q[2];
u3(-0.117884412733765,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.55485593548549,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.44890243306569,1.13821945948160,1.83753060578828) q[2];
u3(1.61165534907469,5.35417195390048,-0.525876310176494) q[3];
u3(2.69313485411580,0.947242996913850,1.29241320743726) q[4];
u3(0.392733984244616,-3.78123915155483,-1.34628024251326) q[5];
cx q[5],q[4];
u1(2.21900796390355) q[4];
u3(-2.71166540029166,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.21528386518474,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.52053796070905,-1.10305522957476,4.59407822797977) q[4];
u3(1.63579710833595,-2.42955016883882,-1.97306344877639) q[5];
u3(2.08710626347664,2.28055225071433,-1.17736264332771) q[4];
u3(1.54844042548534,1.00326394270429,-1.93714181308278) q[5];
cx q[5],q[4];
u1(1.70689595062221) q[4];
u3(-2.82121946805189,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.471805591297327,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.33012851079360,0.966584112902849,-2.81688328433560) q[4];
u3(1.33830328182588,-2.85926134914439,-2.88439775153857) q[5];
u3(1.89208480489006,1.05857195941889,0.365001508224475) q[0];
u3(0.818853070559008,-1.38647447866380,-2.13956723595298) q[1];
cx q[1],q[0];
u1(3.00872828696827) q[0];
u3(-1.34509160359657,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.44429897823226,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.42647749648863,0.738801599648894,-4.39579469693544) q[0];
u3(2.44453312127983,-4.03766361464188,-0.161369338628011) q[1];
u3(2.57443070515093,-0.914198117878449,2.44172577406787) q[3];
u3(2.33450405890726,-1.67729147096330,-0.889353665909511) q[2];
cx q[2],q[3];
u1(2.47913144289108) q[3];
u3(0.403730006375729,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.52613278620307,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.495017166486044,-2.53083808795901,3.10098212882890) q[3];
u3(0.480325127477314,4.03076032471536,-1.42553326104962) q[2];
u3(1.01820861573527,-2.62101161877074,0.0138195604613911) q[1];
u3(1.31332940559571,-2.61238171847305,-0.969475679866329) q[3];
cx q[3],q[1];
u1(1.89286487574084) q[1];
u3(-0.0121735506466021,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.709710928612025,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.82113202119972,-0.779293366967485,3.96733044589113) q[1];
u3(0.506396725918228,0.372258403969031,-3.36475511306256) q[3];
u3(2.67736533133273,-2.74663540785615,-0.0645714193090847) q[5];
u3(2.10226218029747,-3.61899226692681,-2.26479010944631) q[2];
cx q[2],q[5];
u1(0.782150972334519) q[5];
u3(-1.32970166833285,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.91981275162706,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.61046969822789,-0.623329649956551,-2.37130136852205) q[5];
u3(1.10698447851917,-4.43818023278400,0.206450231955531) q[2];
u3(2.46641694545427,1.73634350140463,0.290150050988923) q[4];
u3(1.30795657003134,-1.45987335908306,-2.85863676201215) q[0];
cx q[0],q[4];
u1(1.85999735170544) q[4];
u3(0.503713416119351,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.892402864426386,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.39750630599955,-2.75273073756204,1.37923998975826) q[4];
u3(1.01728676044835,2.56658885428924,1.43229772020759) q[0];
u3(2.15561925515489,-2.35050422587534,1.19013707340706) q[4];
u3(1.30250014294880,-1.23784746821869,1.01507979511738) q[5];
cx q[5],q[4];
u1(3.52151834752338) q[4];
u3(-4.18886825022704,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.756111179520641,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.490271509749607,-2.97474528071698,2.31684024452434) q[4];
u3(1.83429534453974,2.08559442901925,0.364577608162598) q[5];
u3(1.91220567620476,0.737777566706973,2.32413979233763) q[1];
u3(1.88691539609549,-0.939821572212373,-1.18944285103750) q[3];
cx q[3],q[1];
u1(0.623366699781428) q[1];
u3(-0.766771220050575,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.62084694724773,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.301642003694733,-1.72329775266541,0.101710200307304) q[1];
u3(2.02462370532299,2.61656593858500,-1.47054432022935) q[3];
u3(1.39324011148204,3.70564885272183,-0.865559250570403) q[0];
u3(2.13676876088222,3.27596826187772,0.405934839283900) q[2];
cx q[2],q[0];
u1(-0.720793603539613) q[0];
u3(0.275690568796056,0.0,0.0) q[2];
cx q[0],q[2];
u3(4.09843883275362,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.56908117703945,2.23570302510273,-3.50047435912573) q[0];
u3(1.59128762216407,-0.310523076655587,1.87044108832004) q[2];
u3(0.912001529687794,-2.97873626895779,2.32258792137256) q[1];
u3(0.867756131554073,-3.98511071454195,2.23223848621733) q[3];
cx q[3],q[1];
u1(1.55853531459688) q[1];
u3(-0.368616490130324,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.62234071089687,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.440580270274555,-2.66108847095144,2.51493395447111) q[1];
u3(1.43490206296020,3.27474523697181,-2.59048160414843) q[3];
u3(1.00026530187378,0.136896622790307,-0.512091248519749) q[4];
u3(0.938856493916160,-2.91057546942440,1.26696545541994) q[0];
cx q[0],q[4];
u1(-0.119882063323315) q[4];
u3(0.739985050169235,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.54158870013661,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.515126171900598,-2.17961092691912,3.38874694838725) q[4];
u3(1.77885124963373,0.488468026978141,2.14859835106062) q[0];
u3(2.42903773827262,-1.26197667286990,-0.880741385073809) q[2];
u3(0.854569252104788,-3.07359733619296,0.0847769973986503) q[5];
cx q[5],q[2];
u1(-1.07732977423888) q[2];
u3(0.439757015292920,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.91631756810501,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.319254880974127,-0.587502189741041,-1.03878543529009) q[2];
u3(1.79530774829860,-0.782533299560356,2.00923292341575) q[5];
u3(1.86739133819667,0.921001149170116,0.824642034059284) q[4];
u3(1.75671942003944,-1.46127319905325,-1.03904943812993) q[2];
cx q[2],q[4];
u1(3.00952697665387) q[4];
u3(-2.73111724246848,0.0,0.0) q[2];
cx q[4],q[2];
u3(-1.24199189033285,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.27472980031905,-2.69440577656912,2.82561414009713) q[4];
u3(0.644260441295143,2.29626197719104,0.561698137199916) q[2];
u3(2.60831976024830,-3.26680826279791,1.06501071743479) q[0];
u3(1.35952898770169,-0.0760659959905121,4.03627335090466) q[3];
cx q[3],q[0];
u1(0.484651480742079) q[0];
u3(-3.14338062692354,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.66398313438503,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.569001367609998,1.27912075584225,2.00880627876496) q[0];
u3(1.82329249049037,-2.03130808025888,2.55167816241838) q[3];
u3(1.88430835597482,-2.64811765117839,0.231136664962120) q[1];
u3(1.81658904994053,-3.32570851503723,-0.0636633216417564) q[5];
cx q[5],q[1];
u1(3.55006595763047) q[1];
u3(-0.471988326551921,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.58242399124098,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.66216364978815,0.237047020681465,1.76029370622696) q[1];
u3(2.02647263376462,0.859183814851479,-0.143318163007944) q[5];
u3(0.823941290463317,2.45332656946447,-1.38442831740147) q[3];
u3(0.936003118815958,-3.21598811434057,0.966484126560413) q[1];
cx q[1],q[3];
u1(2.86020392644223) q[3];
u3(-2.15349943699706,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.49263586517275,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.83096750683444,0.717736724802548,3.56986002920195) q[3];
u3(1.66809088228065,2.70657957588062,-0.636035426396161) q[1];
u3(2.17907713761179,2.22428229724497,-0.551197279523586) q[5];
u3(2.71944659442256,3.00260184172369,-0.857368978452666) q[4];
cx q[4],q[5];
u1(4.32835130921126) q[5];
u3(-3.37155759384254,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.556272674133485,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.60553675564234,-4.08655526931597,-0.409835786442030) q[5];
u3(1.70991183418601,-2.62834991762528,2.77720387776864) q[4];
u3(1.17313209156679,2.60939702766040,-2.01769063412788) q[2];
u3(0.655963315046175,1.41005262087611,-2.25451731541576) q[0];
cx q[0],q[2];
u1(2.80122613051371) q[2];
u3(-2.46646598770986,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.53192036508254,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.91623256010156,0.799202438791772,1.86542000091475) q[2];
u3(1.79684466231158,4.24492072360033,0.850661501315905) q[0];
u3(2.21856715137415,-2.20560553848280,2.38467647431884) q[0];
u3(2.36261322262694,-2.05493678595976,2.50486649351299) q[4];
cx q[4],q[0];
u1(1.25845965424034) q[0];
u3(-0.387392636138758,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.44705127129268,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.96979993001812,-1.42616958145212,3.74500870896758) q[0];
u3(2.40110475750336,2.84368506698309,0.0637775240811040) q[4];
u3(2.49139582306576,2.53524334225924,-3.68203581706758) q[3];
u3(0.353597648594973,-1.91091932024365,4.15850897449974) q[2];
cx q[2],q[3];
u1(0.291360321142492) q[3];
u3(-1.16876178660061,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.77688100933956,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.85753332093896,3.16791132579836,-2.97850629488953) q[3];
u3(2.74710190280283,-0.643829110494164,-0.545171318724699) q[2];
u3(2.32561432349086,1.38194417424869,-2.53333807972667) q[5];
u3(2.28448655778431,2.09792371455046,-2.61367942040384) q[1];
cx q[1],q[5];
u1(0.687331402711764) q[5];
u3(-1.58277534899972,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.432584699217626,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.79107686812165,2.61051600541272,-2.45454509363519) q[5];
u3(1.43164308427101,-0.0953166863214230,-2.67048452718533) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
