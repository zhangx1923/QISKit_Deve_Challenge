OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(2.00841057121025,0.202994550844849,1.75854522055423) q[0];
u3(1.57463396224245,-3.08820092837536,-2.87694540108090) q[2];
cx q[2],q[0];
u1(1.29638472484449) q[0];
u3(-0.470760452205125,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.92176010446867,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.886360011121902,0.907218200969244,0.473018884043819) q[0];
u3(2.39540379128643,-1.39830276552259,-1.40940695691125) q[2];
u3(1.91501911449252,0.103227423501650,-0.136904437623701) q[7];
u3(0.770411421527715,-0.384523751048520,-5.14323446623155) q[3];
cx q[3],q[7];
u1(3.02120600861813) q[7];
u3(-2.59160657172485,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.20490282313909,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.26136091300442,1.29258664072439,-2.68680130795345) q[7];
u3(1.65396498731086,0.0937861211196358,5.86157046156340) q[3];
u3(1.97799436756619,-2.66564831679016,2.10786768904875) q[6];
u3(2.72318776087095,-3.01493680948134,1.77679418896398) q[4];
cx q[4],q[6];
u1(-0.166902649442258) q[6];
u3(-1.73146595462631,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.373316009365028,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.59929233785415,-1.49331611918343,1.54824796878119) q[6];
u3(0.927094995512277,2.28204517503583,-0.815475491321479) q[4];
u3(2.01625624645715,1.24448372033280,-3.91147328744183) q[5];
u3(0.506567787075685,-2.81916553720729,3.27321068352296) q[1];
cx q[1],q[5];
u1(2.70757391129692) q[5];
u3(-2.17708561380895,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.232502271235784,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.70284220148424,-2.11978407979845,-1.29717747275133) q[5];
u3(2.16947945067618,0.240583814335261,-1.94948391417253) q[1];
u3(1.56657355182589,-0.419780423366421,2.35458919523039) q[1];
u3(1.57803458099012,-1.60059130146378,-1.65590816310104) q[2];
cx q[2],q[1];
u1(1.61074514520747) q[1];
u3(-0.989019867143625,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.524864757052467,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.04306425855786,-2.56161025789836,-0.289660823114483) q[1];
u3(1.66825901779077,-3.60693711519702,-0.692922453543757) q[2];
u3(2.29053434535364,-0.426232633377727,-1.10389466768584) q[7];
u3(1.38991788485824,0.274893578339078,-5.67502595387469) q[6];
cx q[6],q[7];
u1(3.39014371637511) q[7];
u3(-1.30768754946538,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.60033491691509,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.547845139545978,2.17851918748687,-0.540601098777732) q[7];
u3(1.76346472727655,-5.24816691283469,0.302237121010272) q[6];
u3(1.67619610025582,-0.521715147874619,0.405655132775585) q[5];
u3(1.68183139951018,-2.61818293419963,-1.65555304610907) q[4];
cx q[4],q[5];
u1(1.95692056117760) q[5];
u3(0.143084294571792,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.905446828597681,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.02744894389169,0.547979898254142,1.89617607202135) q[5];
u3(2.49456690723121,1.77016815204954,-0.150739819879123) q[4];
u3(1.15047495003398,-1.50243521156886,-0.0695703760964814) q[3];
u3(1.08172114998949,-2.49737806749924,-0.119220108417567) q[0];
cx q[0],q[3];
u1(2.35372418761593) q[3];
u3(-1.67501953010825,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.206774667943258,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.93979103814872,-1.31261015762134,3.99834129280310) q[3];
u3(1.63537800858051,-1.29337842271750,-3.85921929312775) q[0];
u3(2.01292285749538,-3.81968761651235,0.732754854547933) q[6];
u3(2.15108303445954,0.703660001890656,3.01295141705572) q[7];
cx q[7],q[6];
u1(0.0931561043777125) q[6];
u3(-1.16124482174396,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.11195858499784,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.65030560470688,1.83003412260010,-3.18174185290471) q[6];
u3(1.50003043711174,-0.766296787532461,-1.60322581156713) q[7];
u3(2.65261549667063,-2.51113237287095,-0.357903055640188) q[5];
u3(2.67039716734635,-0.767729586918331,0.372195168449059) q[0];
cx q[0],q[5];
u1(1.80652841912816) q[5];
u3(-2.30515352370518,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.10886803356779,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.22841347621488,-3.10003914378413,0.797078023316312) q[5];
u3(0.762433750135815,0.542783650699957,2.54446413793832) q[0];
u3(2.69935594142452,-1.19729015819358,3.08635017184617) q[1];
u3(2.01529195995284,-0.542603401300021,0.625500660805759) q[4];
cx q[4],q[1];
u1(-0.468753683662337) q[1];
u3(1.12406008447219,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.73672748674651,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.07960265928736,-0.348808346408768,1.75142508720965) q[1];
u3(1.55014234157412,3.55954643140077,1.86333536350994) q[4];
u3(1.77362064893070,-2.58216596486249,1.11807832873246) q[2];
u3(2.13388877958884,-3.25003507180202,0.547949347085008) q[3];
cx q[3],q[2];
u1(1.46756330088855) q[2];
u3(-0.646230987042030,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.02064481967034,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.69673402599268,-1.50092582184609,2.74465262761718) q[2];
u3(1.45726435971960,-2.67246554750774,-0.902601829027878) q[3];
u3(2.19257963541619,2.09424364406523,-3.80780946227491) q[3];
u3(0.506400265268825,1.17567579052019,-0.229109389273217) q[7];
cx q[7],q[3];
u1(1.46838210956275) q[3];
u3(-0.490064389271653,0.0,0.0) q[7];
cx q[3],q[7];
u3(2.30742270237043,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.85733414310483,0.174073593692921,3.42953930987335) q[3];
u3(1.96914138495763,-0.591133769598515,1.22595205852012) q[7];
u3(1.99476368987884,-0.244475998293594,-0.0302489742274997) q[6];
u3(1.11988809439578,-0.124840530287456,-5.26532745080645) q[0];
cx q[0],q[6];
u1(0.530019732287829) q[6];
u3(-1.12668397430085,0.0,0.0) q[0];
cx q[6],q[0];
u3(3.24135548003689,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.19646334304135,-1.43446134177955,2.57684571164239) q[6];
u3(1.91980694360081,-3.61472443840452,-0.480398323471306) q[0];
u3(1.78701805149351,-0.913276314222641,0.230685383471930) q[2];
u3(1.95741343757828,-3.82890668878305,-1.23748883483732) q[4];
cx q[4],q[2];
u1(1.70246584562908) q[2];
u3(0.606364077143700,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.06537078149805,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.44266338451709,-0.0948699255981543,2.00365106155098) q[2];
u3(0.962212558151572,2.08865485031954,2.33928814332562) q[4];
u3(2.21594209109554,1.83913886399922,-0.135868129248852) q[5];
u3(1.73790401437557,-0.150352485820850,-3.74488636277646) q[1];
cx q[1],q[5];
u1(0.460934593481917) q[5];
u3(-1.27532298497547,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.02265419754973,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.897128059292676,-1.82407128254040,3.75443917910292) q[5];
u3(2.49697467603726,-2.23531970242168,2.98714602358679) q[1];
u3(1.24727420367198,0.151328032016104,1.94245992172366) q[0];
u3(1.46998627273014,-0.401114896268911,-1.89903405622595) q[7];
cx q[7],q[0];
u1(2.68094087178056) q[0];
u3(-1.16200480091263,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.643774704486374,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.91646429107887,-1.63359757091612,2.26445077402866) q[0];
u3(0.946518056033521,-1.04787705308522,4.90089171817279) q[7];
u3(1.75340539392555,0.547245591699748,-1.05389775338698) q[2];
u3(2.40368327355783,-4.95383185379966,0.783187300845263) q[3];
cx q[3],q[2];
u1(0.641841512997940) q[2];
u3(-0.244956438129422,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.69696222384249,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.459881041724461,0.757845371459057,-0.952400548148287) q[2];
u3(0.852393975010162,-0.00279705919958095,-3.92340468092010) q[3];
u3(1.14404578557357,-1.58053898185868,-0.224398347917199) q[5];
u3(1.77362968593434,-3.88259606008643,0.652566823213292) q[4];
cx q[4],q[5];
u1(1.95151113918999) q[5];
u3(0.348600540822778,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.04244696772442,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.61078924258122,-0.0698852452568931,-0.278643810909693) q[5];
u3(2.33623573312306,2.85757570642590,0.475734381272865) q[4];
u3(0.803191686567166,1.65811772719440,-0.906075049592583) q[1];
u3(1.23726300293520,0.0181743242164942,-3.49620305784503) q[6];
cx q[6],q[1];
u1(0.514863576428093) q[1];
u3(-1.33565314151884,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.17682842236354,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.10739372354597,1.23220141857957,-5.02510724415801) q[1];
u3(1.26967949869857,0.0227590703275622,-5.93814445892437) q[6];
u3(1.42401720094149,1.39965175109844,0.789387177868398) q[0];
u3(0.906064649157725,-1.45891906309903,-1.59123696220260) q[2];
cx q[2],q[0];
u1(1.35781625997224) q[0];
u3(-0.828943342330138,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.11241286485294,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.92392996738228,-1.34489023454627,-1.76842673935699) q[0];
u3(1.36745296293107,-0.539504662559047,2.09157046293195) q[2];
u3(1.70363944934126,2.03125264848481,-2.91336112749738) q[7];
u3(0.685513899152175,2.57066349361748,-1.67524952407830) q[1];
cx q[1],q[7];
u1(1.46652113662971) q[7];
u3(-0.185048219735603,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.54667861472353,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.04605734139449,-2.41949009602683,0.404541849370975) q[7];
u3(1.13212767227888,1.49739851229625,-2.74942106711950) q[1];
u3(0.681007053340943,-1.67011334697396,1.26390235330407) q[3];
u3(0.552198293176241,-2.82485927078033,1.96942873212374) q[5];
cx q[5],q[3];
u1(1.59231754414884) q[3];
u3(-3.54694426533296,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.59656335501445,0.0,0.0) q[5];
cx q[5],q[3];
u3(2.12881183882565,3.20509476236048,-1.86451086191370) q[3];
u3(1.95057926678229,-2.56792951312870,-2.54279802423264) q[5];
u3(0.678795180378472,0.582882895105869,0.731150578191201) q[4];
u3(1.42024673131911,0.677613389163923,-2.49536489191633) q[6];
cx q[6],q[4];
u1(-0.377543155246126) q[4];
u3(1.16218270626898,0.0,0.0) q[6];
cx q[4],q[6];
u3(3.33089461907542,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.88593438825888,-0.649307715987218,1.46904444800533) q[4];
u3(2.39939523561758,-2.07360202214400,0.238257121755015) q[6];
u3(2.42716162777293,0.574918486958952,-1.94237600507820) q[4];
u3(1.98255636476063,1.48898598516171,-3.84655259421115) q[3];
cx q[3],q[4];
u1(4.19518225034653) q[4];
u3(-3.59755509478191,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.241497215795195,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.00026007673493,-1.89369914342997,4.19888306547238) q[4];
u3(2.25199956633385,-3.17586928229363,0.526196275761275) q[3];
u3(2.43087955150343,2.21070501057367,-0.931949453087305) q[6];
u3(2.69937463495139,3.45747106430330,-0.0576512005519145) q[7];
cx q[7],q[6];
u1(1.37547203807426) q[6];
u3(-3.32700392239676,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.26480383284354,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.404077548574269,0.632716871360820,-1.09867962189684) q[6];
u3(0.561211238735841,-0.235450395844812,-6.00956310444762) q[7];
u3(2.61072423739828,0.301790892131910,-0.855189226272651) q[1];
u3(1.41004974526704,0.584768876597462,-4.41688608947620) q[2];
cx q[2],q[1];
u1(1.90615624123324) q[1];
u3(-2.69526759203704,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.199117762868864,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.25009859490816,-2.09396713740907,4.12060806110787) q[1];
u3(0.690795350615353,-0.0442877318197095,4.63748358527378) q[2];
u3(2.03549952206334,-1.61589402666966,4.64796908528361) q[0];
u3(0.485330531860143,2.80208637160488,-1.04555535476087) q[5];
cx q[5],q[0];
u1(1.68834003256062) q[0];
u3(-0.626955879902206,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.93058928380486,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.307369046602513,1.28919736923224,-1.59817892170719) q[0];
u3(1.74629250740658,0.853834438652125,-5.11080176809406) q[5];
u3(1.60994557795059,2.26893185719668,0.0403326646893709) q[5];
u3(2.40460935626743,0.0796582024390446,-3.47133565434966) q[4];
cx q[4],q[5];
u1(2.11505074349151) q[5];
u3(-2.84090719235941,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.12418989793830,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.89605860302002,-0.672036390558163,4.05659549838424) q[5];
u3(0.830117144016006,-0.00935151441615933,-2.47507650488810) q[4];
u3(1.05900023505128,1.40022671336593,-3.68335833628788) q[1];
u3(1.75731300324268,-2.36754081028839,3.78431876377600) q[3];
cx q[3],q[1];
u1(1.53820070828083) q[1];
u3(-0.490030673048899,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.88159676743885,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.26071494846569,0.256679505855280,-0.273006687123117) q[1];
u3(0.617274502821400,-2.65108859285557,-2.28981961157451) q[3];
u3(1.94989596600700,1.08256057793346,-2.66987256285920) q[6];
u3(0.529486914713592,-2.94407835467715,2.25507582073582) q[0];
cx q[0],q[6];
u1(1.37716938282781) q[6];
u3(0.0648041158818322,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.930108916831068,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.284620577916934,-1.04439698194645,0.168832207164345) q[6];
u3(0.787837137818429,-3.56736065852036,-0.921117394386810) q[0];
u3(2.09706228646165,-0.296676490896073,1.29962923338064) q[7];
u3(2.29926096201962,-0.321919549403752,-1.82321993409949) q[2];
cx q[2],q[7];
u1(1.14921034156608) q[7];
u3(0.218889604571932,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.60149797025397,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.39237187264205,0.839734778821068,-1.19239348674785) q[7];
u3(2.73120845409958,1.06301150693723,1.86962758414386) q[2];
u3(1.62810146157294,0.324845081311043,-1.67682845751727) q[4];
u3(0.625466234709411,-3.63522540611108,0.906681018433552) q[0];
cx q[0],q[4];
u1(0.195957064082556) q[4];
u3(-1.02753094938759,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.52881468330702,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.85637103083696,2.94137666682884,0.749481096105362) q[4];
u3(0.692559415362388,5.33968583924953,-0.371383289333243) q[0];
u3(2.35555969959568,1.99080724219843,-3.03214830909277) q[1];
u3(0.974534781860062,-2.86742746354993,2.74835393890953) q[5];
cx q[5],q[1];
u1(1.65626793109692) q[1];
u3(-2.08629688078208,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.898799196590889,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.91463834727599,-1.05735209295343,2.14425946794329) q[1];
u3(0.772148107425594,-2.56253188205436,-0.402211730184691) q[5];
u3(2.06171984044422,-2.59740248015036,-0.352995251287182) q[2];
u3(1.59574658590054,1.21977031728381,4.60889355016236) q[6];
cx q[6],q[2];
u1(1.41718273437842) q[2];
u3(-0.709860549881201,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.212486906227587,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.977695742987986,3.98236980533735,-1.59306227520275) q[2];
u3(1.13218727715861,2.21967629051550,-0.251920852923439) q[6];
u3(1.78237914076915,-1.64451600500727,0.866803275954752) q[3];
u3(2.13743710817684,-1.58785714728122,-0.961668795879531) q[7];
cx q[7],q[3];
u1(2.72128123039686) q[3];
u3(-1.99284161795502,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.222903872993176,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.49183015002747,-2.76574290269645,1.23440940505362) q[3];
u3(1.70709752552449,0.164041697928408,4.93477476030864) q[7];
u3(2.27797303840969,-4.46479830648040,1.44758949629769) q[3];
u3(0.265021119354542,3.87985296059609,-1.52527663040216) q[0];
cx q[0],q[3];
u1(0.689467009320524) q[3];
u3(-0.974714007411256,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.92131555696816,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.76633588891174,0.417120861650498,1.01442147483053) q[3];
u3(2.10653412839302,-0.580982848947794,-3.96290483467156) q[0];
u3(2.36591403783253,-2.44648024006366,0.112147657904966) q[6];
u3(2.28183460885769,-3.01775862110748,0.164715513768332) q[7];
cx q[7],q[6];
u1(0.456723137148581) q[6];
u3(-1.36197385561805,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.40605414168642,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.69551186264807,2.89878554072318,-0.306744227795294) q[6];
u3(1.69336805297600,-1.01369146995630,-0.633970132563354) q[7];
u3(1.97585751389381,-2.25067389350153,-0.312097239636602) q[4];
u3(2.29946475386912,-4.76217508503365,-0.896994562144266) q[5];
cx q[5],q[4];
u1(0.989380865060298) q[4];
u3(-0.200750541747222,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.77197707046147,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.89619192630216,0.743262705837246,1.02019424592756) q[4];
u3(2.24519228280472,-0.958232985759763,5.04750496301224) q[5];
u3(2.01664015034343,0.243225328734254,-3.19153742676694) q[2];
u3(2.16647024037648,3.76149234684256,-2.17843279942339) q[1];
cx q[1],q[2];
u1(-0.284186663866821) q[2];
u3(-1.53496930701760,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.773157378936559,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.42818700736384,3.68436846258636,-1.46926446933943) q[2];
u3(1.85070311352252,2.49555055496041,-2.02199324986227) q[1];
u3(2.17743770336845,0.464274221948416,1.23372376893590) q[3];
u3(1.76698081298878,-2.04536409778413,-2.69107658120435) q[6];
cx q[6],q[3];
u1(2.37856977871824) q[3];
u3(0.0480049355474090,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.29031143851050,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.04231285512931,-1.54132270062724,0.343732516606098) q[3];
u3(0.229739833253535,0.215848302815290,0.929430768007673) q[6];
u3(1.18952697946700,3.16103404543836,-0.531510474186822) q[2];
u3(2.43377156254843,2.58597661092918,-0.828468778802077) q[4];
cx q[4],q[2];
u1(0.427941797375848) q[2];
u3(-3.33537062481255,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.77949264129477,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.273501745580879,2.53357098098131,-1.14320335168003) q[2];
u3(2.01362404194946,-0.368852378092849,-1.71292517329933) q[4];
u3(1.58630874346207,1.24000724189726,-3.13550822816220) q[7];
u3(0.921644523123463,2.73213427104583,-2.96312956234440) q[0];
cx q[0],q[7];
u1(0.879906623476374) q[7];
u3(-0.441302222471303,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.85729013781148,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.659785919045338,0.410460085442507,-0.628981703797532) q[7];
u3(1.61844286937684,-2.81381996419831,0.279223139131239) q[0];
u3(0.252265098426313,-2.23262137417202,1.66080213198965) q[5];
u3(0.171928357584541,-0.0821566251136736,-1.50271726456609) q[1];
cx q[1],q[5];
u1(0.300044638605105) q[5];
u3(-1.65601508540497,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.94054976021798,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.200702183709906,-1.30681239043138,0.374173870429043) q[5];
u3(1.99860997193113,2.46174518442829,-3.00373535415501) q[1];
u3(1.80772799035772,1.71542408670544,-0.437095697190960) q[6];
u3(0.717207295293451,-0.0502594614361689,-3.63574756838185) q[2];
cx q[2],q[6];
u1(1.89193872625733) q[6];
u3(0.222195464416297,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.13014766135367,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.39893355224606,-0.283870672793940,-3.71045785600392) q[6];
u3(2.22896714940567,0.930455960037927,5.34195494190719) q[2];
u3(1.53588241484549,0.531570154723035,-1.48974228567538) q[7];
u3(2.45172631732440,-3.93207455282375,1.23713216109965) q[0];
cx q[0],q[7];
u1(3.38406381386497) q[7];
u3(-0.879661703984869,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.36771634920398,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.24223093578339,-0.618006975688528,-3.26474985100411) q[7];
u3(1.84857897409743,4.28805350379949,0.0210655865968139) q[0];
u3(1.33379617206916,0.0563976632863297,-2.14607559310668) q[1];
u3(2.12257827469739,-3.73356234182477,2.34965352118185) q[3];
cx q[3],q[1];
u1(1.68009446495927) q[1];
u3(0.428193923345233,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.613975698624707,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.13003203674512,-0.157519915642962,2.69779696800181) q[1];
u3(1.45488092871687,2.05047128808283,-4.13589227770361) q[3];
u3(0.957749186828309,-1.09486715715257,0.789520874246880) q[4];
u3(1.32980944703549,-1.36417528888584,-1.27686588996031) q[5];
cx q[5],q[4];
u1(0.178098999366974) q[4];
u3(-1.67651309662982,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.19593366484080,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.820657005996930,-2.01398185332190,0.728528106240619) q[4];
u3(1.09907799352393,-2.11274119636980,-1.49760340634968) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
