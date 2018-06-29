OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.48364324545409,-2.85418915031714,0.923229611375816) q[3];
u3(2.60541515053412,-2.99463544212728,-1.92557444502074) q[4];
cx q[4],q[3];
u1(2.12384187116966) q[3];
u3(-1.79352188345183,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.523804171973406,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.648844103933000,3.66429810946550,-2.12397868525307) q[3];
u3(1.92224712374164,-3.41971640191879,0.893481382955740) q[4];
u3(1.46522122476837,3.07904383792952,-0.600198077541586) q[0];
u3(2.03028578925795,2.43426382436462,-1.87305255748996) q[1];
cx q[1],q[0];
u1(0.0142438864457692) q[0];
u3(-1.21410144616926,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.39270685974537,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.58856105013456,1.42891789010898,-0.855501271573735) q[0];
u3(0.805134042949345,0.812639657436163,1.63534898759956) q[1];
u3(1.12217307172301,2.29236788659637,-1.58469125143273) q[4];
u3(0.539603465963068,2.23001781722288,-2.64528828222357) q[3];
cx q[3],q[4];
u1(0.560025170227938) q[4];
u3(-1.28394397301630,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.94691498588635,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.11525383666008,3.46900432265226,-2.11946348859717) q[4];
u3(2.03924664067328,2.23330860836886,-2.21330229512080) q[3];
u3(1.85003118794666,-0.410067973576752,-2.32312656238643) q[1];
u3(0.948259360648011,1.52312814752406,-3.80804425995908) q[2];
cx q[2],q[1];
u1(1.09114509848861) q[1];
u3(-0.501039850141187,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.17588142151336,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.94385035473163,-3.40533799622641,2.81722369816263) q[1];
u3(1.54370863201726,0.762771318946412,2.23877099808624) q[2];
u3(1.69498135707432,-0.121315003648048,0.789488371904231) q[1];
u3(1.76907728618424,-2.64699223162064,-1.80131608531971) q[3];
cx q[3],q[1];
u1(0.713758270637218) q[1];
u3(-1.28270761141568,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.02224203461643,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.556020415603492,0.973607589241115,-1.18103057626622) q[1];
u3(1.94792960587288,1.05231404424674,-1.60641018978623) q[3];
u3(1.89116328339090,1.67503993765993,-2.86523852992374) q[2];
u3(1.34838081134025,2.50527810909560,-3.65252121935001) q[0];
cx q[0],q[2];
u1(2.08259514948018) q[2];
u3(-3.27868910781162,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.19983729212550,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.707880953298178,0.139042813892527,-4.10337416195053) q[2];
u3(2.10618539492364,4.09491534946784,-1.48165957901258) q[0];
u3(1.31113783085892,0.238818937630246,1.98180895908504) q[0];
u3(1.53608945026808,-1.16325495566098,-0.690466653779366) q[1];
cx q[1],q[0];
u1(2.10857970228092) q[0];
u3(-2.84845138089733,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.75502882382246,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.48830196712220,2.57593278348596,-2.29557228657995) q[0];
u3(2.01376419546862,0.314920108200935,-3.91345250319815) q[1];
u3(1.13101095686682,-0.0520398513892063,-0.127197156475570) q[4];
u3(0.145027329408531,-2.12791701588872,-0.457644374558738) q[3];
cx q[3],q[4];
u1(-1.04783438304941) q[4];
u3(0.251670267953626,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.85686524759516,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.51058630100815,1.20089851094822,-0.525722690146971) q[4];
u3(1.94119301236391,-3.11165118366182,0.903462491700986) q[3];
u3(0.900550377229916,-1.24050473002222,-1.00615359018068) q[2];
u3(0.703203153964279,-3.66556574346061,-0.197363838884172) q[4];
cx q[4],q[2];
u1(2.46070949273790) q[2];
u3(0.246585524706320,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.43519891657165,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.50734166200232,1.18254766634746,2.52919219371193) q[2];
u3(1.40794875434440,4.92379077345058,-0.138425887739119) q[4];
u3(0.748030564597743,0.379524694887145,2.14154932014631) q[3];
u3(1.13247462028743,-1.85175281046540,-1.49459527952407) q[1];
cx q[1],q[3];
u1(2.36974766262578) q[3];
u3(-3.09806097566546,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.39235322807531,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.694047464174384,-1.71568320568154,1.44213928626752) q[3];
u3(0.913568959633645,-1.28133396939775,0.358119401663595) q[1];
u3(1.20184255679255,-0.232610662078583,1.12209906285569) q[2];
u3(1.47332877717070,-0.742451222026960,-1.41943138158888) q[4];
cx q[4],q[2];
u1(3.02447925286665) q[2];
u3(-2.00707938507147,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.507053810624868,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.61311885368128,-0.283602224690469,0.994980194614794) q[2];
u3(1.18617733906073,3.99219025559203,-0.928558601604089) q[4];
u3(1.63072786417860,-1.13383952860081,1.05480255500331) q[3];
u3(0.825596755828764,-3.23181888639778,0.706225176378301) q[1];
cx q[1],q[3];
u1(3.70661440862032) q[3];
u3(-4.37993304794182,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.821731211282217,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.58381482923537,2.52713792539839,-2.45745455128733) q[3];
u3(1.52853325160330,0.421083215044825,2.41237865872055) q[1];
u3(1.74275729730257,1.29432673462037,-3.24225477528729) q[2];
u3(1.19985175524697,2.91319624599806,-2.84876401074312) q[1];
cx q[1],q[2];
u1(1.36619279134814) q[2];
u3(-0.563284030710679,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.210656515917739,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.694716484225102,0.355946953432225,-3.59020010558178) q[2];
u3(1.38802189815660,2.92796617177739,-2.25314472885089) q[1];
u3(2.89957259110274,-4.10914830922657,2.05433934657513) q[0];
u3(0.787507530545119,0.720996771057959,1.24086788139022) q[4];
cx q[4],q[0];
u1(-0.360692642952444) q[0];
u3(-1.74596928929375,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.706897514383544,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.89789677174537,2.79333413058520,0.557761005655240) q[0];
u3(1.78816600646966,-0.565057118905052,3.72602499371766) q[4];
u3(1.81975959955415,0.934356979599503,1.44175728890684) q[1];
u3(0.980543972076299,-1.21125025448292,-2.51873924233244) q[0];
cx q[0],q[1];
u1(1.61336384025636) q[1];
u3(-2.34005661521948,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.458475773705338,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.95400695721621,0.715123938913989,2.08760553380559) q[1];
u3(2.33631292591889,0.471831056437138,-4.09516766328879) q[0];
u3(2.57915697675762,-0.155732020283978,-0.629458210253435) q[4];
u3(0.650008579694134,-5.04879597278097,0.707099991183709) q[2];
cx q[2],q[4];
u1(3.70885692282218) q[4];
u3(-4.23197511783589,0.0,0.0) q[2];
cx q[4],q[2];
u3(-0.181420037145422,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.46835120532364,0.802029381756059,-1.94138526972325) q[4];
u3(0.671842337813702,1.12682870748944,1.66902815071277) q[2];
u3(1.47179189887291,0.861207825251864,-0.0351097288652165) q[1];
u3(0.439329442308623,-0.0105631877668151,-3.50052815114494) q[4];
cx q[4],q[1];
u1(2.06929939142444) q[1];
u3(0.219267827400029,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.53523051724920,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.97254308416310,-2.20903736914621,2.98530068217568) q[1];
u3(1.38766759346208,-2.50131059490463,-1.91304739722538) q[4];
u3(0.950408166455530,1.39908500854500,-2.59519853733080) q[0];
u3(1.47633919986503,-4.87549420754205,0.961704251677470) q[2];
cx q[2],q[0];
u1(1.15794468944577) q[0];
u3(-0.142178407843097,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.57247313868780,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.25710532392134,-1.46553631289029,-0.143807361264849) q[0];
u3(1.38472790839638,-0.798677835475222,2.75214383168778) q[2];
u3(1.57427507089676,0.885070574037230,-2.96793937476583) q[4];
u3(1.58304114179195,-2.29544502539173,3.49348459872246) q[3];
cx q[3],q[4];
u1(-0.213284316961909) q[4];
u3(-0.907294051837680,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.82427233871384,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.755717731549453,0.639611150061902,-0.890215513745849) q[4];
u3(2.82342129584553,-0.107617360585340,-4.13296983531894) q[3];
u3(1.47139811698910,1.05938384380443,-0.731540176577456) q[1];
u3(2.50203624555054,-1.09695229416482,-3.93169685695669) q[0];
cx q[0],q[1];
u1(0.0814322689465237) q[1];
u3(-1.75916310467483,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.03674394570090,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.65237863886683,2.84010256638587,-2.29545772728653) q[1];
u3(0.892920621841391,-0.767457480335585,-3.82360517850840) q[0];
u3(1.61191950248736,0.897727901145339,-3.72494216479625) q[3];
u3(1.53062045079229,-2.06284539084005,4.19313116619424) q[0];
cx q[0],q[3];
u1(2.84913447664743) q[3];
u3(-1.88834847660568,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.613366007685103,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.90564904509423,0.703703472444143,0.856288394635535) q[3];
u3(2.42894800387128,-0.364455814564227,-2.44722239773827) q[0];
u3(1.94292973789678,1.17864761630816,-3.29514639682441) q[2];
u3(1.15620071664962,-2.08725788153147,2.48984598527968) q[1];
cx q[1],q[2];
u1(2.70458827393688) q[2];
u3(0.0876447126102313,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.72045247863130,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.97538873157533,0.526588338624315,-1.43502917106121) q[2];
u3(1.89144387239529,0.555721344816205,3.69790669622570) q[1];
u3(0.660759159899284,2.01349601139643,-2.73294373059001) q[4];
u3(1.21555761684047,2.46689565501984,-3.75619151363601) q[1];
cx q[1],q[4];
u1(3.14984975780147) q[4];
u3(-2.64104919501405,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.40331423523560,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.19290612150529,0.881753533829433,-0.777755374337837) q[4];
u3(1.32890873594932,-3.05860153075961,1.98410604695643) q[1];
u3(2.23866393277033,0.0235645495544237,1.75289504630280) q[2];
u3(1.91162146160197,-1.99550044104962,-3.10961154284499) q[3];
cx q[3],q[2];
u1(3.57590592495524) q[2];
u3(-0.537598625914087,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.60221992287404,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.19998278334559,-0.773189874143126,4.79448253993865) q[2];
u3(1.64428075780697,1.19108224938911,1.26179138554515) q[3];
u3(2.62668210157141,2.51511230102186,-0.734077237954563) q[2];
u3(2.05937442439472,1.24169184104270,-4.15375671491107) q[1];
cx q[1],q[2];
u1(-0.311142122185184) q[2];
u3(-2.46749734507461,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.28112065826021,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.46770574590453,0.369207586254615,2.84312581169571) q[2];
u3(2.34009469687415,-1.37503237940527,-3.85801009976988) q[1];
u3(1.71853175258672,2.92781694028323,-2.04823477512153) q[3];
u3(1.43477720479015,2.50661471493240,-2.41302907714418) q[0];
cx q[0],q[3];
u1(3.20486338480293) q[3];
u3(-1.00659011632243,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.07072196625137,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.753918448429643,1.39584876233175,-2.12416952804704) q[3];
u3(1.41756501517012,4.38697060471538,-0.663706074916307) q[0];
u3(1.12536108934192,-1.27737540172099,0.129175254989752) q[4];
u3(2.02650408718463,-2.67237926193488,0.667841450710047) q[2];
cx q[2],q[4];
u1(1.47464950986868) q[4];
u3(-0.273791326031012,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.23413283263085,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.552985797099001,-0.920980225790722,-1.18541745323309) q[4];
u3(1.65552471896538,-1.59059027315890,3.55126951542662) q[2];
u3(1.73936031824538,-2.27764580265825,-0.153004996550755) q[3];
u3(2.01916100810183,-3.80055833048366,-1.22380872458827) q[1];
cx q[1],q[3];
u1(0.688665140116863) q[3];
u3(-1.11742530196542,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.30549266310548,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.568756720095808,-0.919129257130499,-0.529045018513904) q[3];
u3(1.91164464795159,1.70410560572311,-3.47280904688684) q[1];
u3(1.61683887407064,1.88761063728742,-2.95193091115012) q[1];
u3(1.39743408097588,-3.04326111387861,3.00502308816816) q[0];
cx q[0],q[1];
u1(2.22270218752369) q[1];
u3(-1.66776134953679,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.617985454538994,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.77897979324775,1.41200400212075,0.866579843616204) q[1];
u3(2.09311182396389,0.499784086159580,3.49694477350931) q[0];
u3(2.83554290030139,-1.07737902123263,2.22809565278252) q[3];
u3(2.50253647329466,1.00923071147342,3.66935550309724) q[2];
cx q[2],q[3];
u1(3.23314348562687) q[3];
u3(-0.772629553733662,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.81607971531133,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.25175102956205,0.992188214713916,-1.41999326107497) q[3];
u3(0.582518673213165,-1.43608605508719,-3.34106998987234) q[2];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
