OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.48673553293714,1.24607238673067,1.45605708218660) q[7];
u3(1.44513238902520,-1.34334737717408,-2.44573065694640) q[1];
cx q[1],q[7];
u1(3.19675331345997) q[7];
u3(-1.50639936721693,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.45673870353071,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.47827952317171,-4.61332246613162,1.57946917850744) q[7];
u3(2.28980766349967,1.50506738839557,-2.63727507414782) q[1];
u3(1.44326036957887,-2.69050669238518,3.40754745616169) q[9];
u3(1.62073725163052,2.44984931939571,-1.43664708521529) q[4];
cx q[4],q[9];
u1(1.21284002659493) q[9];
u3(-0.334461889109398,0.0,0.0) q[4];
cx q[9],q[4];
u3(2.70779557661722,0.0,0.0) q[4];
cx q[4],q[9];
u3(2.53023888135139,0.609332492989691,3.85132527988656) q[9];
u3(1.67228259978030,-4.38276723434782,-1.25080190103041) q[4];
u3(0.572802293855033,-1.41483050133453,0.602825317654661) q[8];
u3(0.116721705871746,-2.01569603472634,-0.253205218065123) q[6];
cx q[6],q[8];
u1(1.63554279108248) q[8];
u3(0.00970978789592158,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.20087947250669,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.01291212656573,3.03592510697776,-2.53876858708758) q[8];
u3(1.17315823810207,0.139254474400167,5.90022747153887) q[6];
u3(1.59420323534756,3.62817930677556,-2.03611209762710) q[0];
u3(1.52992180714850,2.46964006518778,-2.10210138185666) q[5];
cx q[5],q[0];
u1(3.12672575920369) q[0];
u3(-2.27422109585962,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.22146903313738,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.79079288809210,0.841645227666342,1.06797084963330) q[0];
u3(1.75705181560453,0.298596282228788,0.734858998448267) q[5];
u3(2.17358994224709,1.81937410661398,-2.92917685920887) q[2];
u3(0.933625305149385,-2.76679877729030,3.11925525196796) q[3];
cx q[3],q[2];
u1(1.90728598643867) q[2];
u3(-1.60950542974396,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.79538974179082,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.28642690595890,-3.29498988002444,2.16022821313061) q[2];
u3(0.328408337973651,2.98358382522316,0.109671515478937) q[3];
u3(2.78992709354649,2.56141212516254,0.283806669400849) q[0];
u3(1.87831371848069,0.306743160371254,-3.74626806195188) q[8];
cx q[8],q[0];
u1(3.67277592756641) q[0];
u3(-1.45755476381507,0.0,0.0) q[8];
cx q[0],q[8];
u3(2.23043716753813,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.23420428457663,1.84033159034066,-1.57534665679993) q[0];
u3(1.81259707269152,-1.62133503984839,-2.58399745152370) q[8];
u3(1.41105941434163,2.85354284399177,-3.03568640066926) q[3];
u3(1.21186490197881,2.77708333614908,-2.07811708281632) q[4];
cx q[4],q[3];
u1(3.23355052997381) q[3];
u3(-4.17415003849788,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.563773674543584,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.06658999704059,-2.51955561735389,2.09406906977491) q[3];
u3(1.66419038429924,-0.0720740142262273,2.66484534631094) q[4];
u3(1.17108810923388,2.91177980746730,-0.635513647265462) q[2];
u3(2.23666653567059,1.51577322652179,-2.20578373978875) q[6];
cx q[6],q[2];
u1(0.0777625946889373) q[2];
u3(0.874361660503798,0.0,0.0) q[6];
cx q[2],q[6];
u3(3.56036965861020,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.34352196331991,-0.452657115290433,-2.17100067685512) q[2];
u3(1.51811526569789,-0.967461793028747,0.963067183209666) q[6];
u3(2.02216823594833,-2.57142113660774,0.612050861535451) q[5];
u3(1.83067441118916,-3.42213487976913,-0.115812831644257) q[7];
cx q[7],q[5];
u1(0.454264511409155) q[5];
u3(-0.0955180327105736,0.0,0.0) q[7];
cx q[5],q[7];
u3(4.12796211211909,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.10328678118193,-1.21474953264900,2.46535178216647) q[5];
u3(0.563934040238388,-5.16234617219127,0.797756477565379) q[7];
u3(0.946703034997015,-0.547596436316154,-0.0550756446263194) q[9];
u3(0.558265094369893,-1.13199194095963,0.609185451621364) q[1];
cx q[1],q[9];
u1(-1.09008304850168) q[9];
u3(-0.281791970291230,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.15976180384966,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.34243842743070,1.38700220201118,0.457291719791275) q[9];
u3(1.24327700828325,2.75111490230317,-3.50382368838717) q[1];
u3(2.68556578374382,2.46974428551961,-1.07028589475717) q[6];
u3(1.93612750316943,1.37551228543503,-3.86235262492285) q[8];
cx q[8],q[6];
u1(2.68458626842209) q[6];
u3(-1.94952474368050,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.528677041782078,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.532564753183751,-0.880149695280225,-2.40600828585253) q[6];
u3(1.97678689001849,-5.10036411805994,-0.139906702834657) q[8];
u3(2.31288499915060,2.48280694914560,-0.150523065956832) q[2];
u3(2.67595320299067,0.0565130952376034,-4.70115120927820) q[7];
cx q[7],q[2];
u1(0.933690205263731) q[2];
u3(-1.22295678519953,0.0,0.0) q[7];
cx q[2],q[7];
u3(-0.293872161786284,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.863333419381578,1.99061179392488,0.889052666830809) q[2];
u3(0.804916366241335,-0.451625968293919,1.43055307644402) q[7];
u3(1.20695700448235,-2.56546341899040,-0.439170913481628) q[5];
u3(1.40595926273657,-2.73054689522983,-0.0772971349334801) q[1];
cx q[1],q[5];
u1(0.876394743697176) q[5];
u3(-3.13408176808196,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.12152022573671,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.47411563161876,1.31349132677799,-0.730415921484271) q[5];
u3(2.52007409485864,3.01977968719304,-0.532656128769121) q[1];
u3(2.43043770787198,2.86214262139165,-1.73340636456314) q[9];
u3(1.56640707124971,1.02688876973220,-2.82039210242254) q[0];
cx q[0],q[9];
u1(-0.582979067348780) q[9];
u3(-1.56850981525749,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.12371702359017,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.49886640232278,0.461404987965954,1.72170177533619) q[9];
u3(1.91576282060007,3.37193749065874,-1.85068567997423) q[0];
u3(1.82952982030310,-0.598045656714647,-0.533364782090656) q[3];
u3(1.16324595531049,0.935232445142475,-4.87597524588546) q[4];
cx q[4],q[3];
u1(0.894410644890235) q[3];
u3(-0.136681294893886,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.80646044171022,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.46462732951354,-0.603600729283969,-1.13029304954975) q[3];
u3(0.713140953524786,0.461919076778473,2.48713799459631) q[4];
u3(2.13483830206664,1.10211301850973,-0.714959973421957) q[6];
u3(1.94650600723254,-0.0846004922489649,-3.05240066656578) q[2];
cx q[2],q[6];
u1(1.78118461924232) q[6];
u3(0.348334537841498,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.597120127666601,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.36568124527365,-0.851716864699069,3.52769169173793) q[6];
u3(1.85344539922772,0.924111899252079,4.76793704229090) q[2];
u3(1.15167434580878,-1.22536702968851,-0.767381891277817) q[3];
u3(1.94243717135037,-3.86860917601585,0.804905574965578) q[8];
cx q[8],q[3];
u1(1.21935564159155) q[3];
u3(-0.218175654937055,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.47311848694447,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.744766454232964,-3.24678371762334,1.24874675651861) q[3];
u3(0.506509353927762,-2.34396838294767,-3.19171017114491) q[8];
u3(1.02346580485831,0.301178063293527,2.67645191522457) q[4];
u3(1.51311805659081,-0.338076791894432,-1.02377646537910) q[5];
cx q[5],q[4];
u1(1.43419554200121) q[4];
u3(-3.43005464196040,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.27143189061678,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.80506601323510,2.85934150602293,-1.35517089207450) q[4];
u3(1.40737940239146,-0.514429841815231,5.30162590546815) q[5];
u3(1.29807811077712,-0.0840770248128648,-2.19680468938724) q[7];
u3(1.67413986206819,1.05446492847603,-3.99518319488506) q[0];
cx q[0],q[7];
u1(1.58752966095738) q[7];
u3(-4.25326421542582e-5,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.22139232929089,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.80052719942798,-1.42214764466017,-0.248203350386602) q[7];
u3(2.01139933826913,-2.24563642577952,2.68707927496605) q[0];
u3(1.52110977458666,3.25252769440432,-0.666086290749243) q[9];
u3(2.33306852311586,2.85105732288421,-0.406926047321163) q[1];
cx q[1],q[9];
u1(2.63739187799335) q[9];
u3(-2.23373798325468,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.140680206230221,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.763888025217520,2.11394687146798,-2.18467766394031) q[9];
u3(1.77692655685668,-0.00728708363826791,-0.800276609534606) q[1];
u3(0.495670072150783,-0.229025204281658,0.00169877331648711) q[6];
u3(1.31933936135601,-3.91805552997988,1.21699141202711) q[5];
cx q[5],q[6];
u1(-0.565852427901979) q[6];
u3(0.138618034381585,0.0,0.0) q[5];
cx q[6],q[5];
u3(4.08859211027263,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.15430021690609,-1.97582092429184,-0.550668643631770) q[6];
u3(2.16937536379461,-1.03899740644119,-2.75801602624898) q[5];
u3(0.929182623301586,2.15278102625632,-2.46305227485824) q[9];
u3(0.186423735331765,-3.71011504484069,1.28091111411938) q[3];
cx q[3],q[9];
u1(3.48397699337973) q[9];
u3(-4.50898918752004,0.0,0.0) q[3];
cx q[9],q[3];
u3(-0.400436327662275,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.67210775436334,-1.55186607432392,2.43465066503730) q[9];
u3(2.41076745614818,-1.98949371759776,-3.92367957992130) q[3];
u3(1.03417275758237,0.527349997341078,1.89273110812881) q[8];
u3(1.01228922254431,-1.46627155324877,-1.61946929137495) q[1];
cx q[1],q[8];
u1(3.46631418325059) q[8];
u3(-0.867836689000601,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.68714090443022,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.07605399808598,4.46579481962759,-1.64624756671997) q[8];
u3(1.17813465651343,1.36259305043148,3.47470837673918) q[1];
u3(0.841928119374549,-3.44441560538338,2.79394417290544) q[0];
u3(0.702802712760567,0.861438108598537,-1.96361541424676) q[2];
cx q[2],q[0];
u1(0.110481539553817) q[0];
u3(-1.36482216990831,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.529235160458067,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.36514840555182,-3.47441845923042,2.64525914883117) q[0];
u3(1.53934210553386,4.40268014293923,-1.40970731796623) q[2];
u3(1.26373298777847,-2.28680117930643,2.38496727692669) q[4];
u3(1.43259747443501,-1.72807231446374,1.15174305272014) q[7];
cx q[7],q[4];
u1(-0.0794803248657383) q[4];
u3(0.946833489576221,0.0,0.0) q[7];
cx q[4],q[7];
u3(3.82232831530583,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.18412312146925,2.75024804223781,-3.11432976441529) q[4];
u3(1.64937792294891,0.820432786072211,0.175827451413380) q[7];
u3(2.75930088430965,-0.489238303703691,0.0691120813567416) q[9];
u3(1.07661479964622,-2.89723996148311,-1.21114677406468) q[7];
cx q[7],q[9];
u1(1.50107858975278) q[9];
u3(-3.33967691612424,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.34801522041526,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.03596476435873,3.10733067479492,1.23087567755015) q[9];
u3(1.24415516837249,-4.68851176516176,0.841297320274942) q[7];
u3(1.81897016094371,2.14542157762638,-0.666050583411157) q[4];
u3(2.17985357499550,-0.996447167433296,-5.07214864193182) q[3];
cx q[3],q[4];
u1(0.492637480367712) q[4];
u3(-0.838409320637348,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.58419987167157,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.33226595111502,4.41666750359858,-1.17570152278023) q[4];
u3(2.75743632810809,3.60339615239357,-2.50060876095730) q[3];
u3(1.15860878170238,1.19488887808423,-3.25737646900578) q[8];
u3(2.04792540875098,1.84396795803717,-4.19513288301210) q[1];
cx q[1],q[8];
u1(0.288306271966672) q[8];
u3(-0.683540955248646,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.95137790545376,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.50017440998430,-3.31381862986050,1.09377486013472) q[8];
u3(1.97462478759732,1.36752831233710,-4.70332361929606) q[1];
u3(0.725678539707140,-0.625672309143295,0.266311985399001) q[5];
u3(0.292201611582187,-1.65766766295740,1.06027285955688) q[2];
cx q[2],q[5];
u1(3.40099792755277) q[5];
u3(-4.30577324215927,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.662152683475281,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.90799918852640,3.97635385630791,-1.72721528885987) q[5];
u3(1.79544028235327,3.52221411409351,2.19143716921122) q[2];
u3(1.48584020400882,-0.351003210631873,-1.17568822566812) q[0];
u3(2.88368585787194,0.510807325686159,-4.91847975224224) q[6];
cx q[6],q[0];
u1(0.0788044490405391) q[0];
u3(-0.696945888143164,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.65874286691215,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.64558735299382,0.110928869100387,0.875604313194435) q[0];
u3(0.890124536061860,1.05979567521223,-4.45335276144248) q[6];
u3(0.629821798078365,3.00426092971458,-2.84350953582799) q[6];
u3(1.36897366616974,0.431555904043377,-2.03736776495469) q[4];
cx q[4],q[6];
u1(1.08770619365978) q[6];
u3(-0.952855628021677,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.72481480172297,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.52174553648755,0.503037517537747,-3.82120159133818) q[6];
u3(2.17678729190204,-2.85106134274971,-0.968016360837436) q[4];
u3(1.39051070997783,2.80528870352880,-1.35899204708726) q[8];
u3(1.40691364843183,1.72964594957825,-0.323130744109934) q[0];
cx q[0],q[8];
u1(-0.134612434150557) q[8];
u3(-2.17288423573885,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.46806671765374,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.18757428351910,2.12062160231258,-3.54253758801035) q[8];
u3(1.56687634062856,-2.72943005472930,-2.88013165701320) q[0];
u3(1.05840385515206,0.246836087262130,-1.50279125453685) q[7];
u3(2.12962060884598,-4.86681839318308,1.33218133999002) q[5];
cx q[5],q[7];
u1(0.922427999482461) q[7];
u3(-0.987859979325412,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.92903259647581,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.86888299769503,3.83860206272925,-2.32968613960308) q[7];
u3(1.17777624030613,-3.66445621618379,1.59134439418002) q[5];
u3(0.939556462547744,0.449679378465677,-2.62134811909560) q[2];
u3(1.83232926473876,-3.04435636283627,3.07310738328207) q[3];
cx q[3],q[2];
u1(0.449347009681858) q[2];
u3(-0.238251556376179,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.12324359074056,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.39539508540823,1.03976023709097,2.43415152031282) q[2];
u3(2.07658589586184,2.90181795592902,-3.29669714625975) q[3];
u3(1.69661003773332,-1.69654754316131,0.775499881180759) q[9];
u3(1.86914277812101,-3.34792976865234,0.194577872885434) q[1];
cx q[1],q[9];
u1(2.48536461518996) q[9];
u3(-2.93321666980264,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.37152820027893,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.11916249591995,-0.524314665134607,-0.324655555153414) q[9];
u3(1.94983750679526,0.500400248499642,0.0163582915846429) q[1];
u3(1.28591698001009,0.757770148244415,-2.70000686286486) q[5];
u3(1.64246301219475,-2.60263104409754,3.09839996681126) q[6];
cx q[6],q[5];
u1(0.337590357500958) q[5];
u3(-0.586887878510353,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.94041728431137,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.454091459424675,-0.0732107809699059,2.71707683105405) q[5];
u3(1.81155915332777,0.735669494698725,-4.51053237686834) q[6];
u3(2.17977560871441,-0.307197631823342,-1.26357470052395) q[1];
u3(1.88215509555780,0.878118042083141,-5.14601762435963) q[0];
cx q[0],q[1];
u1(1.36057506027106) q[1];
u3(-3.45346479824772,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.29546096136760,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.26792164726142,-0.00599734494632265,2.85120347210605) q[1];
u3(1.87565401786040,-3.25299464564087,-2.24084176205386) q[0];
u3(1.79351315645683,2.12161487816548,-3.01332899434974) q[8];
u3(1.95268378424939,2.08431652216331,-3.59832873365272) q[7];
cx q[7],q[8];
u1(-0.0395407270253836) q[8];
u3(-1.31308565009500,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.97714056313639,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.53081791011534,2.81233807119911,-1.29374624059637) q[8];
u3(2.22962186128041,1.40140582310214,-2.68509373914664) q[7];
u3(2.12922035673647,-2.84033599437492,2.15820865901754) q[4];
u3(0.399010735912203,-1.76849817784662,3.15320492997879) q[2];
cx q[2],q[4];
u1(0.618954617972139) q[4];
u3(-1.21125202991942,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.27000782195999,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.911733892602946,-0.195625334421086,4.05247287223989) q[4];
u3(1.27902311993888,2.26585410226054,3.20596729481968) q[2];
u3(1.57568875953260,1.37355208734642,0.0847950594412649) q[9];
u3(1.96633819503100,-0.234728596706742,-4.15766679786354) q[3];
cx q[3],q[9];
u1(1.15398454905992) q[9];
u3(-0.688772533445612,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.0187757760600158,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.29656848082174,-2.84083706417254,-0.563315227511293) q[9];
u3(1.26556303363449,-0.0651818890658959,-5.14454844026926) q[3];
u3(1.12469671964760,0.638846973561944,0.494040945791151) q[9];
u3(1.20076206391277,-0.665502103862365,-2.74015597713668) q[3];
cx q[3],q[9];
u1(1.62571787434176) q[9];
u3(-0.929908784621052,0.0,0.0) q[3];
cx q[9],q[3];
u3(-0.755015307604344,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.28748806685849,1.59071015395437,-1.73946321170460) q[9];
u3(2.31010434933251,0.102421944424428,4.39755608001027) q[3];
u3(2.31473039823684,0.790156074111245,1.01424260822745) q[2];
u3(0.820930315460647,-2.49364053308658,-2.32184949992576) q[1];
cx q[1],q[2];
u1(1.69593036420083) q[2];
u3(-0.413241340433116,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.46175347901054,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.25834821224485,1.75301724012338,-1.84386634451482) q[2];
u3(1.88975502463662,-1.37023853485294,3.48109449965942) q[1];
u3(2.95681070537492,-2.43515063019125,-0.673638514007208) q[8];
u3(2.40705491703785,2.28151295674046,3.83614724716534) q[7];
cx q[7],q[8];
u1(0.0616818265842729) q[8];
u3(-0.807425209484925,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.95570333536541,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.91894706008449,2.87293758053378,-2.34812817627321) q[8];
u3(1.01139196501119,-5.79035023233936,0.0219529829694585) q[7];
u3(2.45742172737985,3.11263450818027,-0.724182859519881) q[5];
u3(1.26468277487225,1.77477259507140,-1.64886311883664) q[6];
cx q[6],q[5];
u1(1.22497967118235) q[5];
u3(-0.565873642227012,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.14925679601574,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.158099223378009,-2.05196074812870,0.148899648278788) q[5];
u3(0.444084211945251,0.0313606819239467,-2.51562930129114) q[6];
u3(1.97036520847125,0.258762639551258,-1.31130206603146) q[0];
u3(1.70991189088340,-0.0621510703609254,-3.27492373346346) q[4];
cx q[4],q[0];
u1(1.49106443060823) q[0];
u3(-0.811953522787666,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.05126165402533,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.836963323762200,-0.379950190039715,1.02860367788262) q[0];
u3(2.01663703553932,2.18596402804682,2.44076680558550) q[4];
u3(0.998430798419712,-2.07051633261693,-0.0945690847124404) q[0];
u3(1.59980161499089,-2.11223548239874,-0.130789708148175) q[4];
cx q[4],q[0];
u1(4.08903765946013) q[0];
u3(-3.17421792095277,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.520560770327963,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.83562330024002,2.99758105311441,0.181138609503998) q[0];
u3(2.16178962233061,1.85614516670823,0.233599496453599) q[4];
u3(1.04359893799888,-3.80401801339533,0.933764027500502) q[1];
u3(2.19242566235580,0.413792770766584,5.81809024448845) q[5];
cx q[5],q[1];
u1(0.777230325761192) q[1];
u3(-3.35259269169473,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.90960863594185,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.22550873462088,0.444213767688693,-1.14647931759425) q[1];
u3(1.50538096205931,0.525157917467844,-3.43440243669699) q[5];
u3(1.71135572941415,0.755110949309258,-0.182060230207088) q[2];
u3(0.476154206991054,0.350468867232470,-4.35913409098053) q[9];
cx q[9],q[2];
u1(0.0244642028831306) q[2];
u3(-0.588594581825643,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.88669102318733,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.07414049725915,3.84178223739524,-1.61771577956512) q[2];
u3(1.72294748838921,-0.467284165452882,0.686093313507025) q[9];
u3(1.84424232790836,1.00788309533345,1.09800041346516) q[8];
u3(0.0679202424388775,-3.64823059439297,-2.15782714296028) q[3];
cx q[3],q[8];
u1(-1.01682578092099) q[8];
u3(0.439802136988691,0.0,0.0) q[3];
cx q[8],q[3];
u3(3.54442933149893,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.09944064751588,0.883346228273039,-2.79813262858107) q[8];
u3(2.30454416083463,-2.83635932995655,-0.657009146813206) q[3];
u3(0.461201351612216,-1.76962178567592,1.36066101178704) q[7];
u3(0.692862731967883,-2.87585452195166,1.35793523683323) q[6];
cx q[6],q[7];
u1(2.52373914240281) q[7];
u3(-1.97028072002843,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.482294058312951,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.25393817917733,0.285651321380231,-1.51200744047343) q[7];
u3(2.04814709167193,2.22521919281985,-3.20819390531194) q[6];
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
