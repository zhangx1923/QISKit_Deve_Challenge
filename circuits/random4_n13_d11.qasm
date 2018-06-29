OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.28668390198032,1.33493526909984,-2.57853382958651) q[5];
u3(0.920618259559524,-3.19947559586731,2.42842420979198) q[11];
cx q[11],q[5];
u1(2.37612581277194) q[5];
u3(-1.54286579892561,0.0,0.0) q[11];
cx q[5],q[11];
u3(3.34153818236365,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.79714610131150,2.65194243692974,0.186609605833561) q[5];
u3(1.95827196894131,-0.0578032396400232,-4.95025281665279) q[11];
u3(1.81432732964665,1.78152643333854,0.779292937642204) q[1];
u3(0.930067841584909,0.271054443722326,-3.29666768198974) q[10];
cx q[10],q[1];
u1(2.76927939023440) q[1];
u3(-1.66325860198569,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.13948049946609,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.621430854214772,1.96516821667796,-0.666377554020538) q[1];
u3(1.16307706748820,-1.79178260501665,-1.44949877957777) q[10];
u3(0.623772997215545,1.51879267428245,-1.16518070550011) q[0];
u3(0.370021201434800,0.185820856514230,-2.58813146927972) q[3];
cx q[3],q[0];
u1(1.82329826404486) q[0];
u3(-3.25751975126654,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.07250468627228,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.35297660457733,2.18160284214381,-1.47295254609548) q[0];
u3(2.23364856456036,2.47526705471678,-0.429140059350114) q[3];
u3(0.320644382794883,-1.27352094378595,2.25301364762321) q[8];
u3(1.18986949368313,-2.98032011301486,2.04267573468209) q[12];
cx q[12],q[8];
u1(0.921056685657654) q[8];
u3(-3.31926763886807,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.93909278243047,0.0,0.0) q[12];
cx q[12],q[8];
u3(0.973261171740451,1.67591115703730,-2.48708410466598) q[8];
u3(1.68900506413817,-0.375258644998280,-4.18394457307400) q[12];
u3(1.92681897425677,0.668235274084686,-0.746919826606766) q[2];
u3(0.971389160409490,-4.42967904146235,1.28493330041800) q[7];
cx q[7],q[2];
u1(1.80860945186926) q[2];
u3(-2.21127986782197,0.0,0.0) q[7];
cx q[2],q[7];
u3(3.49759862879895,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.71078102464763,-2.30417591926890,0.459619727922600) q[2];
u3(1.63620038104427,4.06135854746170,0.643604062032148) q[7];
u3(2.04889295396592,1.23856229372527,1.40147043302865) q[4];
u3(1.52910142744299,-1.74163561945330,-2.11232900043822) q[9];
cx q[9],q[4];
u1(-0.192590522280992) q[4];
u3(-2.46515671739633,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.54243660926656,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.23623645507187,-0.563634136883136,-0.181108467465393) q[4];
u3(1.69991492718793,-3.25962354868956,-0.426788662530159) q[9];
u3(1.23697479186983,-1.64385704252508,2.26717854063684) q[3];
u3(0.661489639820439,1.11041415839235,-2.19114268536382) q[1];
cx q[1],q[3];
u1(3.66800162055018) q[3];
u3(-3.31118218660915,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.959167075928000,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.0311762869900405,-2.08879690832952,1.77204001664564) q[3];
u3(1.75561011398935,-1.23128088362755,3.05258527520579) q[1];
u3(2.88388511858463,0.873011382232369,-2.70305205477004) q[11];
u3(2.69620001616889,1.08387667136606,-3.50527236621207) q[9];
cx q[9],q[11];
u1(2.63267064230537) q[11];
u3(-3.08681434793256,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.969349564784284,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.63131185203397,-0.852827214896737,1.35775303021109) q[11];
u3(1.02380551724121,0.230876257148660,0.355427490995778) q[9];
u3(1.42249725391004,-1.14255944724546,0.328303758563904) q[5];
u3(1.09950659036529,-2.94488372342657,0.400542120778274) q[7];
cx q[7],q[5];
u1(3.24048431394344) q[5];
u3(-4.42057759800678,0.0,0.0) q[7];
cx q[5],q[7];
u3(-0.313557382216134,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.49637580723449,-0.477952404669509,0.522040064646253) q[5];
u3(0.612045480587381,-4.11114779072355,-1.25312477017725) q[7];
u3(0.900680885006035,0.561842952650706,-1.96727653555239) q[8];
u3(2.16688109486633,2.31208498790841,-3.65753359627661) q[6];
cx q[6],q[8];
u1(-0.690612620787856) q[8];
u3(1.02720839509423,0.0,0.0) q[6];
cx q[8],q[6];
u3(4.24647077160393,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.467610178479973,-1.87275665533936,1.37027326920508) q[8];
u3(1.41213862015734,1.90724288552375,3.13594232990067) q[6];
u3(1.12398953864826,-1.18018976688338,-1.14270838375098) q[10];
u3(1.27087167852344,1.59697561805395,-3.89960168087501) q[12];
cx q[12],q[10];
u1(0.462569573762150) q[10];
u3(-1.30582297328568,0.0,0.0) q[12];
cx q[10],q[12];
u3(2.14432164394245,0.0,0.0) q[12];
cx q[12],q[10];
u3(1.89787129036847,-2.67560439395310,3.52971421895105) q[10];
u3(0.840730495224280,0.0298348347722061,0.707551091469683) q[12];
u3(2.01738572250403,-3.71751865246646,2.21793802165411) q[4];
u3(0.666552052281014,1.53159112293028,0.363887852815938) q[2];
cx q[2],q[4];
u1(3.04612797509994) q[4];
u3(-1.83181985690331,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.743803421429376,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.16985154462932,-0.835389519159343,0.814700503126351) q[4];
u3(1.97387518826389,-2.38960884161313,0.885613569510062) q[2];
u3(2.00979077345473,1.63348725490444,-3.65347793662929) q[11];
u3(0.933667784784250,3.63156481262555,-2.39436274802848) q[0];
cx q[0],q[11];
u1(1.67975146970212) q[11];
u3(0.0835313853102597,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.709865869672048,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.480634878578720,-0.621083088362693,2.59935628941612) q[11];
u3(1.45003591023036,-0.453640582560912,5.07934988741518) q[0];
u3(2.00833756950168,1.02271211206232,-3.58584902442474) q[5];
u3(1.55429363456412,2.03405414985060,-3.04422879767061) q[1];
cx q[1],q[5];
u1(2.96117097455198) q[5];
u3(-1.26181061770096,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.359619169218137,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.46521406620624,-3.32849659384078,0.378790576347771) q[5];
u3(2.05928488249584,-3.81384414061181,-0.0221264240124337) q[1];
u3(0.405401062557098,-0.710571240322661,0.802314114523652) q[2];
u3(0.745356611131192,-1.42213190277480,-1.59218989212344) q[9];
cx q[9],q[2];
u1(1.99849956153977) q[2];
u3(-0.145720389803102,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.22158708157094,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.54860229068421,1.84282625524755,-1.13958715530107) q[2];
u3(0.733535088118711,1.37446900350807,-3.33356240361247) q[9];
u3(0.869938341919773,2.54820227997443,-0.921313789811406) q[6];
u3(1.38681008462329,0.517362363155166,-2.50530385481267) q[4];
cx q[4],q[6];
u1(1.76335342276258) q[6];
u3(-2.79071573854820,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.640528656114002,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.73885764030345,-4.40710262961559,1.57395712584879) q[6];
u3(2.41120996511298,0.914481135777416,2.84173789685394) q[4];
u3(1.47037042452629,2.34224339869066,-3.29671030996312) q[12];
u3(0.492970993397086,-2.50234484594098,3.22224439011653) q[3];
cx q[3],q[12];
u1(0.904613591263995) q[12];
u3(-1.24386572028665,0.0,0.0) q[3];
cx q[12],q[3];
u3(0.220404525005611,0.0,0.0) q[3];
cx q[3],q[12];
u3(0.603588224191060,0.375660943190454,-4.27579330929256) q[12];
u3(1.18463352746249,-4.51909834567405,0.645660917054640) q[3];
u3(2.23118846510335,-0.253059535384784,1.59480571960645) q[10];
u3(2.14201631878697,-0.198677437539256,-1.48188271797880) q[7];
cx q[7],q[10];
u1(1.54300101674207) q[10];
u3(-3.66884516793618,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.35353368180120,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.644310202896011,-0.0735791831319703,-1.01893211461691) q[10];
u3(0.940081174813998,-1.41831949846660,-2.92291979029630) q[7];
u3(1.00956357912482,0.760196766819667,-1.97965923317969) q[12];
u3(0.967200925329910,1.17404754654938,-3.68025441405583) q[2];
cx q[2],q[12];
u1(1.90437746793308) q[12];
u3(-3.09466356448453,0.0,0.0) q[2];
cx q[12],q[2];
u3(1.15343599448247,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.65853061512141,-1.90348905738709,2.78891559525001) q[12];
u3(0.990732967596575,-0.335228611506140,4.21925345583711) q[2];
u3(1.59250952658644,1.71685862232002,0.410799006899750) q[10];
u3(0.862936930438051,0.429396024796312,-2.95380958913854) q[8];
cx q[8],q[10];
u1(0.509635753212774) q[10];
u3(-1.16809446917838,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.23461410191831,0.0,0.0) q[8];
cx q[8],q[10];
u3(0.823239625798494,-2.37464675056537,2.44285506649846) q[10];
u3(2.20878594578286,3.14149557018758,0.251473945793214) q[8];
u3(0.715996297509692,-1.63608544935877,1.49643767990455) q[3];
u3(1.24664674915584,1.89275056503036,-2.59601782254866) q[4];
cx q[4],q[3];
u1(3.97971097088639) q[3];
u3(-4.16500373801914,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.749159543598565,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.51527491135108,-3.26902841697633,2.23823145134708) q[3];
u3(0.821542732978940,-0.306180070642898,-1.70226752066563) q[4];
u3(0.642826686898609,-0.560829487419528,0.964461362105787) q[6];
u3(1.35056708624972,-2.63162229343275,-0.0613198700449893) q[7];
cx q[7],q[6];
u1(1.90134309947729) q[6];
u3(0.175749279157202,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.54405377373557,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.673591380992092,-1.11361641199192,2.44720332504579) q[6];
u3(2.53384210878733,-1.21301698333004,-0.219456099060986) q[7];
u3(1.60687673980830,-1.20612911239501,-0.215172513325985) q[0];
u3(2.10353409797292,-3.13040468434287,0.664933428090082) q[11];
cx q[11],q[0];
u1(-0.761104546461006) q[0];
u3(0.393279300150652,0.0,0.0) q[11];
cx q[0],q[11];
u3(4.08659501575078,0.0,0.0) q[11];
cx q[11],q[0];
u3(2.10367484403270,0.856048697128201,1.34090060697421) q[0];
u3(0.227852556400472,2.73823243632434,-0.740481937429944) q[11];
u3(2.60158796344927,0.180624742825858,-1.70760861032147) q[9];
u3(1.67251439593855,-4.23763901498677,0.582401723711718) q[1];
cx q[1],q[9];
u1(0.151599535032164) q[9];
u3(-1.44638800001626,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.57370143122503,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.65719041756001,0.395074827006312,3.44011365542038) q[9];
u3(2.52595614925866,-2.21561711066859,-0.550901193218282) q[1];
u3(1.99725136511515,-0.898164988246034,-0.333112013718832) q[7];
u3(1.74427865651584,-2.88995774313205,0.667703848312755) q[2];
cx q[2],q[7];
u1(2.15902458273665) q[7];
u3(-3.32369223568314,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.365690651033676,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.35240955014510,3.02586896781397,-0.567171545473534) q[7];
u3(2.97279251341615,1.43481224975954,3.44801645869106) q[2];
u3(1.60653078253844,-1.20559954522121,-0.416771077938956) q[10];
u3(1.74768431993225,1.27601738404190,-4.48124492966115) q[11];
cx q[11],q[10];
u1(0.792639830634731) q[10];
u3(-0.371563580411885,0.0,0.0) q[11];
cx q[10],q[11];
u3(2.00276361286249,0.0,0.0) q[11];
cx q[11],q[10];
u3(0.744303083323763,2.69351387532042,-2.78862120471821) q[10];
u3(2.64090969050133,2.55343363975297,1.91310398780852) q[11];
u3(1.03404385530661,0.850142217378298,-1.85031066434040) q[4];
u3(2.18178029037687,-4.67071262543532,1.32078706249705) q[6];
cx q[6],q[4];
u1(0.892473980388780) q[4];
u3(-0.189943291529278,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.33719026861206,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.03381354981612,-1.43622722217432,-0.288411203201545) q[4];
u3(1.41467226014226,1.88301765406545,-3.02694613487338) q[6];
u3(1.63786301846130,-0.597713354063553,1.06697786047105) q[8];
u3(1.94988115266694,-2.25902041303532,-1.57387739613808) q[3];
cx q[3],q[8];
u1(1.43274741655414) q[8];
u3(-0.774093996248486,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.60226647584707,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.838545945755397,4.68207210612241,-1.00153664815237) q[8];
u3(1.26771511685512,-0.644054772428937,-3.00649537120468) q[3];
u3(2.38866526959017,0.0208789870660342,-1.28305234785023) q[1];
u3(1.35266649001344,-0.0682060614338149,-4.26654743607202) q[12];
cx q[12],q[1];
u1(0.416009992763219) q[1];
u3(-1.32267346241242,0.0,0.0) q[12];
cx q[1],q[12];
u3(2.93191437892603,0.0,0.0) q[12];
cx q[12],q[1];
u3(2.38031787347840,-2.61336297769819,1.70252875092037) q[1];
u3(0.711434241778262,-0.428941082107195,-4.16158652871667) q[12];
u3(1.78712055082911,-1.80659288760759,1.00223122721717) q[0];
u3(2.06335212628596,-3.83035181265571,-0.685187538592673) q[5];
cx q[5],q[0];
u1(3.99564261168104) q[0];
u3(-3.69539909429013,0.0,0.0) q[5];
cx q[0],q[5];
u3(-1.03147083483107,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.50848028592514,1.04024090280530,-4.45945728051357) q[0];
u3(1.69294760284220,1.39235744001528,3.52455198423689) q[5];
u3(0.864951255905740,1.01187853993472,-4.03417515151123) q[6];
u3(1.41089801035249,2.54198714209352,-2.46032269799290) q[10];
cx q[10],q[6];
u1(1.42628997565763) q[6];
u3(-0.996608075308401,0.0,0.0) q[10];
cx q[6],q[10];
u3(-0.415250160605799,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.61429093056118,2.75121612731999,-1.05770302838407) q[6];
u3(1.68908575905325,0.546864990836430,0.0246648550114298) q[10];
u3(1.38163279605964,0.251533371108788,0.823651938779643) q[0];
u3(1.76128734855412,-0.315661471547408,-2.63796852028273) q[3];
cx q[3],q[0];
u1(3.05824736272990) q[0];
u3(-1.76888438355133,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.14403383194819,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.25218508850883,-0.232306621879801,-1.47595764856987) q[0];
u3(2.61327511976190,0.0729202906452716,0.812750636629831) q[3];
u3(1.74837371675139,1.70167398451158,0.741925549363753) q[9];
u3(0.756747732235758,-0.525324219242074,-2.90471702376504) q[1];
cx q[1],q[9];
u1(1.98092651739280) q[9];
u3(-2.25271463395977,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.229230799552135,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.48814602506686,-1.23464234851179,2.37454317213346) q[9];
u3(1.53167180034318,-1.39259653593837,1.07221172825975) q[1];
u3(0.298787618851115,-0.404931767893375,0.377308568373371) q[8];
u3(0.878526854863992,-0.298906577534575,-1.71558039565689) q[5];
cx q[5],q[8];
u1(3.02386511924920) q[8];
u3(-1.50753645915860,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.47486769788544,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.701326845738695,1.32128911917367,-3.83103333008133) q[8];
u3(1.16491505993618,-1.14384272684396,3.04096923797792) q[5];
u3(1.46524976304715,3.03615892907915,-1.02389148893819) q[2];
u3(0.950864228250970,0.792402178173517,-0.288936874875703) q[12];
cx q[12],q[2];
u1(3.09036775854134) q[2];
u3(-2.42254285280106,0.0,0.0) q[12];
cx q[2],q[12];
u3(1.57419858294142,0.0,0.0) q[12];
cx q[12],q[2];
u3(2.50379094014906,-2.08036360636995,1.49402107103038) q[2];
u3(0.969447417419861,-0.228263584460870,-2.82355250872807) q[12];
u3(1.36229409741242,0.530389113418541,0.600056120736545) q[4];
u3(1.66485680719788,-0.0746825282168884,-2.52330911865971) q[7];
cx q[7],q[4];
u1(-0.179648104627664) q[4];
u3(-2.17730177065323,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.45349609943992,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.15656323833338,1.21897255075453,-0.120639017991638) q[4];
u3(2.02387695440583,3.41755989303952,1.18549113110333) q[7];
u3(1.13817850180228,3.30258607924270,-0.909255618923927) q[9];
u3(2.25873455039858,2.07773053547885,-1.30789175407006) q[6];
cx q[6],q[9];
u1(4.34330365062991) q[9];
u3(-3.62790940315129,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.844356778101277,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.0544806107358452,-1.56903224474469,2.47893148232196) q[9];
u3(1.19669577098074,1.95804368022361,-2.07918453802376) q[6];
u3(2.02620850149418,0.120340644077331,0.323027773177226) q[7];
u3(1.72101985595490,0.0930701280736173,-2.46970604094509) q[0];
cx q[0],q[7];
u1(0.880992030092227) q[7];
u3(-0.181710767791758,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.49897660439071,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.796872818100607,-2.38113246806754,-1.44750152534538) q[7];
u3(2.26687465901954,0.223549707796291,4.89462239112727) q[0];
u3(2.91716897654813,1.98329840458112,-0.663668442691328) q[8];
u3(1.90672102736631,-0.202320148472098,-3.49191919154982) q[2];
cx q[2],q[8];
u1(0.912942475331729) q[8];
u3(-1.49793275155282,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.82783427980899,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.09978119389914,1.52962798652853,-3.04655639191078) q[8];
u3(1.28975486123086,3.82935761003933,1.74734768701373) q[2];
u3(2.22736131341571,1.28754930618232,-3.83627219689741) q[3];
u3(0.782032468999738,-1.45933483084950,3.46411224986931) q[11];
cx q[11],q[3];
u1(1.00406827025017) q[3];
u3(-0.141474541782329,0.0,0.0) q[11];
cx q[3],q[11];
u3(1.99576651620761,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.10639549509217,-4.24729347903959,0.493549116656832) q[3];
u3(1.30842249647784,0.216624200833166,4.37896632299266) q[11];
u3(0.869253437249242,1.96136947931802,-2.44200190637518) q[4];
u3(1.44346588084967,-2.67899009711150,2.95414002289808) q[10];
cx q[10],q[4];
u1(3.49647016799730) q[4];
u3(-3.88646127307404,0.0,0.0) q[10];
cx q[4],q[10];
u3(-1.03876367458714,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.38434504624212,1.79582804663759,-1.86217193210148) q[4];
u3(0.373159125967331,-3.63165640947833,0.563971944202418) q[10];
u3(1.05560942625896,-0.315494002644997,-2.22822323357520) q[5];
u3(2.23069406624385,-4.22530138920488,2.00569616896001) q[1];
cx q[1],q[5];
u1(2.74846606796379) q[5];
u3(-1.78822387367827,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.741873845493736,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.11554542477769,2.83224676360488,-2.08914568281662) q[5];
u3(1.93515785171322,-5.78325767082986,0.427204460419708) q[1];
u3(1.65340378739510,2.73848412082432,-1.83433606814562) q[6];
u3(1.69642495412058,1.05567963044499,-2.59808642293150) q[11];
cx q[11],q[6];
u1(1.73872102642619) q[6];
u3(-2.78657730568583,0.0,0.0) q[11];
cx q[6],q[11];
u3(3.05540028935535,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.834691831156208,2.84163607986123,-1.60394830095154) q[6];
u3(1.06614496341567,2.65502515508521,-1.16630715718427) q[11];
u3(1.25279167012277,2.75563816980555,-2.43401907687400) q[4];
u3(0.602966497135297,2.41956976349421,-2.79783088040422) q[10];
cx q[10],q[4];
u1(-0.148354874432087) q[4];
u3(-2.68140027634148,0.0,0.0) q[10];
cx q[4],q[10];
u3(1.37158941722235,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.62804703978159,-2.52163621861352,3.12503400975141) q[4];
u3(0.893433834603276,-4.74176716305745,0.649488904847231) q[10];
u3(1.50008830789415,1.14137551746703,-2.62712081021455) q[0];
u3(0.798874561638046,-2.78619560092339,2.13940017297318) q[12];
cx q[12],q[0];
u1(1.26834164690489) q[0];
u3(-0.854748780708198,0.0,0.0) q[12];
cx q[0],q[12];
u3(2.80564650029809,0.0,0.0) q[12];
cx q[12],q[0];
u3(1.62261254012686,0.719991048382644,-1.45384376348053) q[0];
u3(0.0638719889594994,4.28121306095315,1.39576754034060) q[12];
u3(1.09540821481375,-1.40333377919393,1.47234773657013) q[3];
u3(0.810790719742546,1.52997586236171,-2.68567917461281) q[7];
cx q[7],q[3];
u1(3.01692428795326) q[3];
u3(-1.87579348057328,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.797248235524030,0.0,0.0) q[7];
cx q[7],q[3];
u3(2.66139561447592,-2.38977708546734,2.24802032235463) q[3];
u3(0.227382025214346,-2.21907063636096,-1.43155505824450) q[7];
u3(1.51307450602484,2.50142790762797,-2.78879022988026) q[5];
u3(1.83686216251536,-3.36373520773475,2.35067760491836) q[9];
cx q[9],q[5];
u1(1.52679513046155) q[5];
u3(-3.08297673498864,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.959547007651665,0.0,0.0) q[9];
cx q[9],q[5];
u3(0.943531676709221,-2.21622137968911,0.778632955179469) q[5];
u3(1.27788456048163,-2.20642171904262,-1.55784547019755) q[9];
u3(1.10503614070744,3.45363346241432,-1.18509037727442) q[8];
u3(2.04759795302590,2.32067233965491,-0.842295239144398) q[1];
cx q[1],q[8];
u1(0.789655108370488) q[8];
u3(-1.21944523992256,0.0,0.0) q[1];
cx q[8],q[1];
u3(-0.203521426371299,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.55590105085659,-2.91711171068095,1.85111207961616) q[8];
u3(1.26827645197235,-4.82843093367076,0.769928716318778) q[1];
u3(2.10714313976155,1.78682919750673,0.564110432900012) q[0];
u3(1.99303101085702,0.636550946405982,-2.20625716903077) q[5];
cx q[5],q[0];
u1(2.48920930542742) q[0];
u3(-3.02409638748298,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.54095843495322,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.27210453280643,3.51804123577682,-1.84176748705634) q[0];
u3(1.93390643419797,0.991045205508216,-3.22804262290225) q[5];
u3(0.631910834138213,0.249648311014833,-2.05889387068283) q[10];
u3(1.81109481331990,-3.01250814412816,2.47280920713767) q[6];
cx q[6],q[10];
u1(2.53112523139096) q[10];
u3(-1.78013669063622,0.0,0.0) q[6];
cx q[10],q[6];
u3(0.829507612351795,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.406674793635225,-0.882978822112216,1.90080826441936) q[10];
u3(0.803674250970415,3.65892581685395,-2.30383902809028) q[6];
u3(1.70262855240344,-4.08218893084254,1.45612373996355) q[3];
u3(0.202782992768405,-1.01712280498067,2.42310222776604) q[4];
cx q[4],q[3];
u1(1.60736345929140) q[3];
u3(-0.878231128612421,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.44131509475874,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.20870635597878,3.04960722067477,-2.24015863787882) q[3];
u3(2.68168964825693,0.0802189336027038,0.588052412535648) q[4];
u3(2.23577510141875,-1.91841560648945,-0.0751832668059719) q[8];
u3(1.75403150333700,-4.58131869153734,-1.36640155744104) q[7];
cx q[7],q[8];
u1(1.65142492215636) q[8];
u3(-2.53262590552204,0.0,0.0) q[7];
cx q[8],q[7];
u3(3.46398558092834,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.08258421143207,-1.08824363980183,0.730531146753341) q[8];
u3(1.86650168909923,-2.94176843749772,-0.786347215335244) q[7];
u3(0.696721943314993,1.29567169282801,-0.248775853779011) q[9];
u3(0.351560977949694,1.13660561390956,-2.60192569278410) q[11];
cx q[11],q[9];
u1(3.38545461781023) q[9];
u3(-4.31383098102414,0.0,0.0) q[11];
cx q[9],q[11];
u3(-0.588576219556109,0.0,0.0) q[11];
cx q[11],q[9];
u3(2.68745671751876,-3.76243287724861,0.777318486124169) q[9];
u3(0.642953193468271,1.08610591465207,5.18231159324097) q[11];
u3(1.01885494039334,1.58083359917300,-3.30805652163749) q[2];
u3(1.11572442057585,2.04938813597884,-3.25672295958449) q[12];
cx q[12],q[2];
u1(1.28819346064161) q[2];
u3(0.215876845535254,0.0,0.0) q[12];
cx q[2],q[12];
u3(2.37794538346264,0.0,0.0) q[12];
cx q[12],q[2];
u3(1.10136193811029,0.279801001190740,-2.58637249782818) q[2];
u3(1.55969454258977,-5.31809319649690,0.220352264229617) q[12];
u3(2.79825095476643,-1.30031496603332,1.47169042659878) q[9];
u3(1.31692869905807,-2.07217473239279,-0.991760224157625) q[4];
cx q[4],q[9];
u1(1.86054682321946) q[9];
u3(0.756710081169233,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.05068581635839,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.81877655658032,-0.249038532146333,-1.15986739176781) q[9];
u3(0.914040141450100,-0.917069558850146,2.08395952620855) q[4];
u3(1.87029509821506,3.54517600303974,-1.48086937950602) q[7];
u3(2.64088507505291,2.47877966342130,-1.41621427624773) q[0];
cx q[0],q[7];
u1(1.97540007466317) q[7];
u3(-2.67724378702464,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.791706654180889,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.988917605748611,0.757212035549459,0.933911180522021) q[7];
u3(2.45215359368018,-3.10652467057150,1.08127800786689) q[0];
u3(0.128996126324003,-2.03925878866553,2.57963162070299) q[5];
u3(0.0724164427418655,1.83204532253621,-3.50716031706630) q[2];
cx q[2],q[5];
u1(1.29639401684688) q[5];
u3(-0.0737287855040347,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.504021671924620,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.17001323070141,0.273926805402930,-0.911464881919059) q[5];
u3(2.68695228794713,3.63453764293478,0.380376560081005) q[2];
u3(1.90859519917526,0.849403585159779,-2.49536575529343) q[11];
u3(2.76330037310324,0.727784107853258,-4.43151734702698) q[10];
cx q[10],q[11];
u1(1.83581453516446) q[11];
u3(0.0804713916138273,0.0,0.0) q[10];
cx q[11],q[10];
u3(0.676243449145845,0.0,0.0) q[10];
cx q[10],q[11];
u3(3.05945243453342,1.96492407971539,-0.258035394624913) q[11];
u3(2.63282821493149,-2.83255543773122,2.22660915260976) q[10];
u3(0.772016318370280,-0.0191001433791902,-0.876202595145841) q[8];
u3(1.44280541577815,-3.76066975376523,2.30284769175159) q[3];
cx q[3],q[8];
u1(0.350037851746230) q[8];
u3(-0.749718314302319,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.22901859006180,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.978296328796933,-1.36143245150551,0.720231997583354) q[8];
u3(0.621923793690162,-3.77665424652754,-1.30821330428808) q[3];
u3(2.22127778759356,3.21264131493967,-2.13998351001324) q[1];
u3(0.902849869360368,2.45442050690811,-2.11995225448322) q[6];
cx q[6],q[1];
u1(1.54780466347726) q[1];
u3(-3.37594001062492,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.21860330122400,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.58310060402801,1.06602989483731,-4.24051432375477) q[1];
u3(1.37163428870495,5.83477092041658,0.294817134600703) q[6];
u3(2.50454172107315,3.32933790715329,-0.744061988596369) q[0];
u3(0.994956485691008,0.779961436506328,-1.75951316557406) q[7];
cx q[7],q[0];
u1(-1.15168468313982) q[0];
u3(0.348715903692868,0.0,0.0) q[7];
cx q[0],q[7];
u3(3.41040535698911,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.685715599877605,-0.186519149211490,-1.13553126823287) q[0];
u3(1.24016146804174,-1.16941012192958,1.17281608651050) q[7];
u3(1.52154087206627,2.40321035614387,-1.22327544708699) q[3];
u3(1.97256638446389,1.53958459664205,-2.10445595258901) q[11];
cx q[11],q[3];
u1(-1.26207892417941) q[3];
u3(0.961042386686947,0.0,0.0) q[11];
cx q[3],q[11];
u3(4.15257448219991,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.40827617421573,-0.707563309475524,3.32119616419863) q[3];
u3(2.57929703568031,-3.72487587974895,1.28183548865076) q[11];
u3(2.47631564027340,3.50979986724563,-0.897252416132510) q[4];
u3(1.95897968620316,1.51752078536134,-1.70051810233460) q[8];
cx q[8],q[4];
u1(2.95418384483041) q[4];
u3(-2.34760081500903,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.25579011145070,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.76117327991855,1.82281317467757,-0.939717897663125) q[4];
u3(3.06612001668087,-1.57893490702368,-3.06310152490966) q[8];
u3(1.90509312668621,-1.88753219566726,-0.165218008309024) q[5];
u3(2.93597601360067,-1.06884550042533,0.687597575662760) q[6];
cx q[6],q[5];
u1(2.69582891693907) q[5];
u3(-1.63657688146538,0.0,0.0) q[6];
cx q[5],q[6];
u3(3.27369656451256,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.36361560113619,-3.91386889387122,2.25701287419914) q[5];
u3(0.331125072285845,2.17002396159132,2.20919201587132) q[6];
u3(2.22064751854246,2.51682001663777,-1.99829087515341) q[10];
u3(1.55648137738875,2.50191044718585,-3.59480295036385) q[12];
cx q[12],q[10];
u1(2.85708917479299) q[10];
u3(-1.63447924428798,0.0,0.0) q[12];
cx q[10],q[12];
u3(0.639549555632978,0.0,0.0) q[12];
cx q[12],q[10];
u3(2.03919503066342,-1.67134326533343,3.60517064187187) q[10];
u3(2.18054694060001,-0.0427016596931797,-5.35566906943190) q[12];
u3(1.21000741173599,-0.341898260683630,-2.43934962137977) q[1];
u3(0.543079091411015,-4.53815871891240,1.17805407145313) q[9];
cx q[9],q[1];
u1(0.578263347301290) q[1];
u3(-1.45185297305190,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.12327789797534,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.497243171577560,1.91655783822933,-4.16596387728875) q[1];
u3(1.47696690471723,2.78822516241422,-2.30811104163685) q[9];
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
