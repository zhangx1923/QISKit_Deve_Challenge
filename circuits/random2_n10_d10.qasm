OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.64596259481487,2.76413127338345,-1.89447131788529) q[5];
u3(0.321235334018136,2.85779960769838,-2.55814825953557) q[2];
cx q[2],q[5];
u1(3.01517938966285) q[5];
u3(-2.34271934046621,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.15103437650423,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.963747843219478,0.992876796977377,2.12205760202236) q[5];
u3(2.24033310178656,-0.256689940962392,-3.89313267558949) q[2];
u3(0.861202631197065,2.93859699390384,-0.502718033344697) q[0];
u3(1.55283002054791,0.394647600146662,-0.986938972493354) q[6];
cx q[6],q[0];
u1(4.49285499695934) q[0];
u3(-3.80084489054666,0.0,0.0) q[6];
cx q[0],q[6];
u3(-0.811127225924769,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.43777804337709,-1.50123570017376,3.90069483577624) q[0];
u3(0.796059029562240,2.42143741118754,-2.87365527170154) q[6];
u3(0.955488256453383,0.158411064697552,1.65829023839377) q[9];
u3(1.48838861758564,-1.38219303089041,-0.962196910618296) q[4];
cx q[4],q[9];
u1(1.34039205098664) q[9];
u3(-0.300679151521482,0.0,0.0) q[4];
cx q[9],q[4];
u3(2.44410309644016,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.31600804161535,-2.17218099856369,1.10495374267656) q[9];
u3(1.28126598620955,-2.10692304516506,-2.00699868208129) q[4];
u3(1.47360210327055,1.21234860954024,-3.76409647426684) q[1];
u3(1.96987243650858,2.19805877735042,-2.46186617394388) q[7];
cx q[7],q[1];
u1(2.62089684079647) q[1];
u3(-1.84706489038583,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.0253588945204721,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.73976537375909,-0.319964298511846,1.36974570464231) q[1];
u3(1.22523015234269,-2.61929028265839,3.18643292502251) q[7];
u3(2.56888768201917,-0.504685775386977,-0.152184761957851) q[8];
u3(0.589610790292170,-0.223678591431891,-4.60342029501977) q[3];
cx q[3],q[8];
u1(2.99563243978890) q[8];
u3(-1.37375739881144,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.86962917887149,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.35104501468771,3.25667939479949,-1.70870923368419) q[8];
u3(1.55498014754424,0.992733352736735,3.87165986110888) q[3];
u3(0.648429424263530,1.59811495463034,-3.67417859442838) q[6];
u3(1.82421822342231,2.01493335559584,-2.55664074157955) q[5];
cx q[5],q[6];
u1(1.75042195013424) q[6];
u3(-2.49125710575672,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.53789138181052,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.31341579959373,4.02667232506935,-2.15059658974202) q[6];
u3(1.07476757656562,0.323681950477823,-3.73855970068660) q[5];
u3(2.16305522481690,1.26405534730542,-4.10055729299537) q[7];
u3(1.19761684819292,-1.86287699959251,4.04944952268407) q[4];
cx q[4],q[7];
u1(-0.281853475131887) q[7];
u3(-1.70384729657033,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.05813070696536,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.52663356513659,2.00674684079220,-2.25388553084186) q[7];
u3(0.917475649427358,0.747391835620644,-5.52700817989893) q[4];
u3(1.68660628345634,1.83561019994180,-0.265460427555595) q[1];
u3(1.01477625146605,0.298923732731561,-3.59877462393023) q[3];
cx q[3],q[1];
u1(-0.161996456109290) q[1];
u3(-1.73118217711696,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.05139148912990,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.61112324044234,2.61105296062119,-1.91563154910195) q[1];
u3(2.34339584240230,2.58051523578072,1.07344541976432) q[3];
u3(0.533085618007206,1.39364806771510,-3.29281563103393) q[8];
u3(1.60583588406232,3.63835932048310,-2.40966273061672) q[2];
cx q[2],q[8];
u1(1.24799736017499) q[8];
u3(-0.941615496780686,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.01022295325667,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.767191069344681,-2.73459199136302,2.87713634893111) q[8];
u3(1.36716299955656,-3.02334525851131,3.18390918159378) q[2];
u3(0.896321681121935,3.79802372861179,-1.77745765738779) q[0];
u3(1.59403714966211,2.12669661654221,-1.11243770763672) q[9];
cx q[9],q[0];
u1(2.71782549176237) q[0];
u3(-2.01289332656480,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.35555109286207,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.15312006706872,0.907134461910723,0.224948886436201) q[0];
u3(0.831090390022458,-0.0205908952515923,3.72864247389039) q[9];
u3(0.418296014015637,0.828103026126308,-1.28877637788016) q[9];
u3(0.916282608694739,-3.72316385555679,1.69709426221934) q[3];
cx q[3],q[9];
u1(1.81141722702714) q[9];
u3(-2.48377581374128,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.879411012248256,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.96033657173426,1.28092703283625,-1.72358116813949) q[9];
u3(1.20952370350283,-2.31821532639811,3.21175010942448) q[3];
u3(1.56535361199883,1.33551966605029,-4.18441128376671) q[2];
u3(2.43833263545011,-1.79748746397158,4.39696081608807) q[5];
cx q[5],q[2];
u1(0.00474834984531980) q[2];
u3(-1.61467538229827,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.708598300241724,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.40140241660877,2.61972221417620,-3.61548920589546) q[2];
u3(0.960904301639774,1.89603161837232,3.06436279281330) q[5];
u3(1.10858453054400,1.26842135200598,-0.601437961593684) q[1];
u3(1.79299502285622,-0.568000448010298,-3.63223594770255) q[8];
cx q[8],q[1];
u1(2.50785374447935) q[1];
u3(-2.28448664434140,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.30535648929972,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.84763484800003,-0.406342058651239,0.104092553607669) q[1];
u3(0.279253912355581,-4.50741521812977,1.60648910016461) q[8];
u3(1.45563030160505,0.938488783862369,-3.44941240788629) q[4];
u3(1.31121854409182,2.84634154100895,-3.15478165404587) q[6];
cx q[6],q[4];
u1(0.146677140540658) q[4];
u3(-1.83203105486034,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.36270745556705,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.96071555365982,-2.56159023605332,3.56012492863842) q[4];
u3(0.914677916546033,0.145145432143549,-6.04431257555974) q[6];
u3(2.36045611222035,0.696881240519561,-2.99594274410526) q[0];
u3(2.37031028057938,0.910855092512890,-4.00915105474018) q[7];
cx q[7],q[0];
u1(0.988437023087885) q[0];
u3(0.00433592366847591,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.52741472069832,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.323768212607934,0.879909910527472,-2.28058675092786) q[0];
u3(2.77333945428366,1.34735867996102,-1.06677633401063) q[7];
u3(0.741729509687150,1.86182845424977,-1.62957480772496) q[5];
u3(1.18791807903147,0.596516489970512,-2.47985309782292) q[0];
cx q[0],q[5];
u1(3.27109502733713) q[5];
u3(-1.27886527185953,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.50152497878093,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.30688212993010,1.75661111769269,-0.207581996026123) q[5];
u3(1.68155360693341,0.343304660905723,-0.777463591688247) q[0];
u3(2.04464073808182,2.33343796714210,-2.89301744248684) q[7];
u3(1.41156069657747,2.21844801558517,-2.62831065465527) q[1];
cx q[1],q[7];
u1(1.65301811907646) q[7];
u3(-2.32395425681658,0.0,0.0) q[1];
cx q[7],q[1];
u3(3.16172216369640,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.88376729110995,-1.67130144666183,0.501741525186604) q[7];
u3(1.20252244662582,-0.371078901099820,-0.544584858932762) q[1];
u3(0.583989132300313,-1.96370197648204,2.36961710350774) q[8];
u3(0.495166180228115,-3.19232438015321,2.02661268022913) q[2];
cx q[2],q[8];
u1(2.57042101983964) q[8];
u3(-1.78030702275833,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.164723315766047,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.723443174436432,1.42916148821714,-1.84622478557928) q[8];
u3(0.182219123956146,-2.13774149294226,4.05793611745881) q[2];
u3(2.29969239503138,0.703759485457638,0.403359882755343) q[4];
u3(1.18016933697537,0.389264343315480,-5.59559720129459) q[9];
cx q[9],q[4];
u1(-0.590616119300587) q[4];
u3(1.02670254708513,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.44617659571850,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.41591767562671,2.49411557732962,-2.64593512434333) q[4];
u3(1.39942776966454,-3.34816785818439,-1.05683382900324) q[9];
u3(0.353557980613364,1.34785550113561,-0.929146067598044) q[6];
u3(0.103041514843317,-0.459035389345610,-0.689330133570295) q[3];
cx q[3],q[6];
u1(3.11754517457706) q[6];
u3(-0.580999061740499,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.84833145303285,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.508274805403042,-0.907206773172190,0.290451872384689) q[6];
u3(2.53410977684749,0.947967058688592,0.882820999887142) q[3];
u3(2.08703109157145,-2.61339911111803,0.343189747587617) q[3];
u3(3.06960897088801,0.0980521640585131,0.588293143065159) q[4];
cx q[4],q[3];
u1(1.35376825989032) q[3];
u3(-1.24421127216358,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.99591878390852,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.777035906145073,-1.46086992031936,2.41350559547115) q[3];
u3(1.82433790928454,4.03221784960759,1.80912773384236) q[4];
u3(2.84365451274355,-3.65026569317506,0.597011391759648) q[2];
u3(2.21844919985591,-1.86363535405931,-0.0491732311787806) q[1];
cx q[1],q[2];
u1(0.405515424251371) q[2];
u3(-0.0417879942313426,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.26939049881321,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.61152582043760,2.61469525438649,-2.32131266857457) q[2];
u3(1.74303500589032,1.22402272119220,2.63265198073020) q[1];
u3(2.94842166426909,2.32462595918324,-2.28112734993340) q[9];
u3(1.90603084126647,2.53035589995186,-2.45554423043642) q[0];
cx q[0],q[9];
u1(2.08014793541183) q[9];
u3(-1.78620624287898,0.0,0.0) q[0];
cx q[9],q[0];
u3(0.168163476905859,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.03912270536717,1.61133770294255,-0.766665900153625) q[9];
u3(1.89785103755550,4.46734570633484,1.06374432768336) q[0];
u3(2.26329178336783,0.615049007145436,-3.71683454849962) q[5];
u3(1.99102137588851,2.40985289154802,-3.00387317980722) q[7];
cx q[7],q[5];
u1(1.52643966364144) q[5];
u3(-2.94903959855209,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.23910992625775,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.41620476146619,-0.853758068472759,1.31015353584112) q[5];
u3(0.928015050133632,-4.25067012613721,-1.44701220282128) q[7];
u3(1.59042322544886,2.97265004876512,-2.34467380191034) q[6];
u3(0.142782100866489,0.681821710815690,-0.522578733402687) q[8];
cx q[8],q[6];
u1(1.20740729052418) q[6];
u3(-0.362628175689506,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.58190441952694,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.36935315276124,0.496953962845983,2.73423243345859) q[6];
u3(1.79948520212378,-1.56683401481738,0.694316774820205) q[8];
u3(0.852338724839160,-3.06052405140573,1.98432256538535) q[9];
u3(1.04389259737058,0.330798812575126,-2.12033368812979) q[7];
cx q[7],q[9];
u1(2.11651259513292) q[9];
u3(-3.63724523509621,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.19383212105720,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.24493241038612,0.796556215945754,0.534649606359011) q[9];
u3(1.26345342028636,0.389619570409650,2.76858622993198) q[7];
u3(2.01450707435393,-0.243891537092512,1.88704110311034) q[1];
u3(1.94244882432904,-2.10570777347532,-2.23018597762712) q[4];
cx q[4],q[1];
u1(1.62164969734154) q[1];
u3(-2.16432621147551,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.52351857515899,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.40529320443519,1.56175044907054,-1.54733054764718) q[1];
u3(2.83814710287548,-2.53609748582175,2.23580964818010) q[4];
u3(1.28215450697259,1.43824643149703,0.693833150921991) q[3];
u3(1.04923182719246,-1.58008851698826,-1.29392626693634) q[2];
cx q[2],q[3];
u1(1.21225773526161) q[3];
u3(-0.0604146749639491,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.84839867359608,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.72192869638179,-3.63475055315820,1.72613239461104) q[3];
u3(2.48267211418091,2.00831682618079,-1.79424470247767) q[2];
u3(1.60407605090326,0.768181957188176,1.59772969873888) q[8];
u3(1.89045633105626,-1.28663356524575,-0.119963323197300) q[5];
cx q[5],q[8];
u1(1.75467093552259) q[8];
u3(-2.06662708414399,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.190254830584046,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.25338114494499,0.288371692320285,-1.95235447410069) q[8];
u3(1.67065876187804,2.75778973276896,3.39622831081985) q[5];
u3(0.880087114935426,1.16730387701888,0.487364030304857) q[0];
u3(0.855206712835357,-0.151028585317722,-1.70217778115233) q[6];
cx q[6],q[0];
u1(0.353185996483771) q[0];
u3(-1.54515884806695,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.42678785852432,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.40354284268300,-2.76643282092922,0.778700774061708) q[0];
u3(0.398649862525722,0.652449917453655,-4.40398618747501) q[6];
u3(1.23692909101672,0.458482307697997,-1.21369516594603) q[5];
u3(0.561346666834985,1.65957159183326,-4.60037090298180) q[8];
cx q[8],q[5];
u1(-0.539154645006787) q[5];
u3(-1.62647091485862,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.675983312476923,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.27587924571456,-0.199406101598733,0.153883486344721) q[5];
u3(0.597195314457742,-4.04965857005128,1.58475129023294) q[8];
u3(1.97316396835068,0.0832291169026972,1.37134110918769) q[4];
u3(1.90518507684271,-1.48979333385032,-1.79952495595272) q[0];
cx q[0],q[4];
u1(1.55220846620071) q[4];
u3(0.717198829815580,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.35841227179147,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.331314602092446,1.90737006173135,2.39237967849079) q[4];
u3(1.14050947740399,-3.39141762329234,-0.249591373230192) q[0];
u3(1.32126725136103,-1.70405299867718,0.559413654720883) q[3];
u3(1.46873264399463,-3.87859953466529,-1.02938742956086) q[6];
cx q[6],q[3];
u1(0.369710140085446) q[3];
u3(-1.61677836355982,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.57118850615765,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.66190421794147,-0.000573490746173411,-1.47998893374734) q[3];
u3(1.40814753415845,1.58356407390275,3.97535518120694) q[6];
u3(2.55847348239590,-0.0904116321755060,-1.75168294160108) q[2];
u3(2.03405329503117,0.531476389694689,-4.46520437811974) q[1];
cx q[1],q[2];
u1(0.449064382607528) q[2];
u3(-1.59685762564790,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.99703262719560,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.39942801350179,-2.34278215227315,1.51934895709246) q[2];
u3(0.471116569916113,-0.896478834768405,3.57925213279808) q[1];
u3(1.22400727211885,0.457266538990383,-2.08885757683735) q[9];
u3(1.87482412944658,-4.20484515411281,1.46361164761930) q[7];
cx q[7],q[9];
u1(2.77173389820167) q[9];
u3(-1.86142005528050,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.879511807714304,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.527579379751029,2.66864836450974,-1.24908393714893) q[9];
u3(0.381547182412012,1.00449034930670,3.12315544185025) q[7];
u3(1.65888462151270,1.26701055527118,-2.71255211142270) q[8];
u3(1.54422158122749,-2.59442362112628,2.68658882508704) q[9];
cx q[9],q[8];
u1(2.36338285703508) q[8];
u3(-1.54723760273638,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.378521881385186,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.05958875809005,-1.55114877823758,0.439152177287077) q[8];
u3(1.69337329444762,1.29914051713415,2.24027109168373) q[9];
u3(2.15303334927843,-2.28108185640682,-0.443702283535273) q[6];
u3(1.17997267141577,-3.67082272893917,-0.476235947711861) q[1];
cx q[1],q[6];
u1(1.53389886297547) q[6];
u3(-0.139142380930497,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.33599463851065,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.514871743511777,4.08622363472967,-1.88279124047729) q[6];
u3(2.64580208862265,0.425799198860057,-3.45978101509766) q[1];
u3(2.91582734955091,-2.45674345650708,3.76495565361816) q[0];
u3(0.403443364105549,-0.341555208368264,1.75634134157873) q[7];
cx q[7],q[0];
u1(3.51205495362889) q[0];
u3(-0.896824743446610,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.09189780812941,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.87851073078898,0.193695332674283,2.70536781194901) q[0];
u3(2.33678290665547,1.59804939719179,-3.65889862218971) q[7];
u3(2.49213057361659,-2.19056845217577,-0.766010426287069) q[2];
u3(1.72825135266443,-5.09223408331942,-0.902719598742109) q[4];
cx q[4],q[2];
u1(1.51422458987103) q[2];
u3(-3.59047953083015,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.10231352596206,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.76903668935698,-2.93543583829488,2.17756416277683) q[2];
u3(0.326110787419276,-5.36768139746075,0.207830039608725) q[4];
u3(1.47175546064235,2.97563208644222,-0.679118057220808) q[3];
u3(0.984738598603794,0.865230677336842,-0.863040934809413) q[5];
cx q[5],q[3];
u1(0.262590866913579) q[3];
u3(-1.19650356587315,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.40742004399421,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.68442298025349,-2.69075992363490,1.08721109354006) q[3];
u3(1.34189765457585,1.59020652618953,-2.56653460632099) q[5];
u3(2.22664094562346,0.233438689868483,-2.78148562402870) q[6];
u3(2.53701309615866,-0.530611091844932,-4.32703909279704) q[0];
cx q[0],q[6];
u1(2.70913071516456) q[6];
u3(-1.75550421731078,0.0,0.0) q[0];
cx q[6],q[0];
u3(0.0250206475154264,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.81360212664443,-0.184980682546101,2.75151551832550) q[6];
u3(1.35677057605202,3.79472203721916,-1.09098458888542) q[0];
u3(2.75024632709418,-1.74638066413642,-1.19026920436248) q[9];
u3(1.00230366750267,0.409290391477463,-5.57244961791995) q[7];
cx q[7],q[9];
u1(2.50005177544884) q[9];
u3(-1.90774373231773,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.92016756630099,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.46136575414158,-2.51222426219806,1.69622601987181) q[9];
u3(2.07470975334184,-4.00199517879245,0.0575223486491767) q[7];
u3(1.56766134493278,0.552121215110416,1.77505667719575) q[3];
u3(1.33929156526615,-0.815637193008889,-0.515015897702768) q[1];
cx q[1],q[3];
u1(1.50401163929420) q[3];
u3(-0.198205517116885,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.27046100146749,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.44780165658988,3.19899405266795,-0.403109791318352) q[3];
u3(1.52452374779586,-0.560721950362486,0.939902463364856) q[1];
u3(1.53345241669185,0.447708905781960,-2.08830216666033) q[4];
u3(2.02094864696741,-2.90764264960543,3.02792238450366) q[8];
cx q[8],q[4];
u1(1.69866037026489) q[4];
u3(0.0835382886423637,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.590821022318696,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.30617242079631,-0.557527421225973,-1.34851093417898) q[4];
u3(1.57551413255429,3.78483428449426,1.57527324794340) q[8];
u3(1.44665948837894,1.00455331333836,-3.40865728609554) q[5];
u3(2.34602232149587,-1.03907635592888,4.38360212585872) q[2];
cx q[2],q[5];
u1(0.896604674434372) q[5];
u3(-1.52134734862675,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.14026127346247,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.12577896936567,-1.57157896851855,4.12357533260082) q[5];
u3(1.14412991746784,2.27659287936426,-3.68739864191755) q[2];
u3(1.30427039508535,1.46775004488733,-2.30857229490794) q[1];
u3(2.12273126398469,1.74590433350007,-4.41925015338246) q[4];
cx q[4],q[1];
u1(3.55900258928459) q[1];
u3(-1.32965728230791,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.64110063781853,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.98169367588810,1.32668161385149,-3.43396935186931) q[1];
u3(0.722598079722381,-3.16592886979363,-2.44405559005134) q[4];
u3(0.785414921699582,1.08549678369171,-0.301495723078096) q[7];
u3(1.54621516181072,0.337475134460540,-2.20126562523647) q[0];
cx q[0],q[7];
u1(-0.955098847352058) q[7];
u3(0.387199080770951,0.0,0.0) q[0];
cx q[7],q[0];
u3(3.43485695769203,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.00888060859059,-1.49245550773658,3.57740153078912) q[7];
u3(1.43158025430467,-0.825706953934239,2.07586548468534) q[0];
u3(1.95552351589820,-2.09248238249149,0.504785570238623) q[9];
u3(1.55379308399222,-3.42340035316591,0.277834972598821) q[8];
cx q[8],q[9];
u1(2.20850965331230) q[9];
u3(0.0552538867893722,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.799694182054370,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.89819454584553,-1.74817440811004,-0.228461366153994) q[9];
u3(2.09615321307417,-1.71490800717785,-1.66291521479221) q[8];
u3(1.29550729328315,-0.341745238683018,-2.15025848795248) q[5];
u3(2.47668575911357,0.188181407044615,-5.67154983177732) q[6];
cx q[6],q[5];
u1(1.84750822360603) q[5];
u3(-2.84234777369270,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.113435062218946,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.802654368748163,-1.16691210677476,3.06984036975946) q[5];
u3(1.99116415241171,-1.93346964055245,0.504181140567464) q[6];
u3(2.17158888588937,1.20722285913812,-3.24532282484631) q[3];
u3(1.61436713824919,3.09026900247056,-3.09848245756714) q[2];
cx q[2],q[3];
u1(1.73097717666128) q[3];
u3(-3.03537606155267,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.868599428188341,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.24777972509262,4.14517249874420,-1.73870141702387) q[3];
u3(1.65290150173849,-1.67421075185736,2.80768166547177) q[2];
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