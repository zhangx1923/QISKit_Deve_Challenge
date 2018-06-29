OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(2.15063907965890,-1.09166118090025,-0.582514777000805) q[1];
u3(0.305195660004344,-2.03453901738221,-3.12027376202938) q[5];
cx q[5],q[1];
u1(0.504134500416019) q[1];
u3(-1.30566851474345,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.24105161629157,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.742228212166524,-1.74413278399362,4.53777284214616) q[1];
u3(2.20893096938020,1.52817540111204,-4.12633461327928) q[5];
u3(0.550350397406095,2.21264142574633,-1.36566068807370) q[12];
u3(1.01534941570731,1.58667443909132,-3.21751095041276) q[0];
cx q[0],q[12];
u1(1.44953749283077) q[12];
u3(-2.64592142884463,0.0,0.0) q[0];
cx q[12],q[0];
u3(0.209880345938430,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.23849442189243,-2.12618442636572,3.79539705652521) q[12];
u3(1.45383199345573,-1.61205414663825,3.04635080411097) q[0];
u3(2.81637798888299,-4.23480831917999,2.00199407255236) q[3];
u3(1.30054512197630,3.45386002130894,-2.08651159733046) q[2];
cx q[2],q[3];
u1(1.33370513667401) q[3];
u3(-0.811029906765256,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.90574986206786,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.15212031031875,-2.97622787574478,2.81389477321211) q[3];
u3(1.43761141126188,-0.621724846918737,-0.189886313046564) q[2];
u3(2.15205910167525,0.417039317137912,-2.09292842589529) q[10];
u3(2.38116045156892,-0.404766342983800,-5.37867422739752) q[7];
cx q[7],q[10];
u1(3.30358717508644) q[10];
u3(-1.36383319213979,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.61115255632820,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.11874340230231,1.27124204375553,0.0890699598180099) q[10];
u3(1.69337637436365,0.439566592371235,-0.569271560583121) q[7];
u3(2.37482318836039,-1.37186612974151,1.69877156795550) q[9];
u3(2.57007028380838,1.31360774820992,1.77614161398327) q[11];
cx q[11],q[9];
u1(0.696833353772568) q[9];
u3(-1.47215779434108,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.83784800664876,0.0,0.0) q[11];
cx q[11],q[9];
u3(2.66467710978383,-1.73654797327966,-0.538949314793336) q[9];
u3(0.507067626871036,-1.61702705269431,-4.30417940642383) q[11];
u3(0.578985046594528,-2.12816505715813,1.64746272895732) q[8];
u3(0.679893917983417,-2.25537953701337,-0.244947416691063) q[6];
cx q[6],q[8];
u1(1.35677117033351) q[8];
u3(-3.52110065668807,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.62997463443041,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.51006900162484,4.36310914811037,-1.77657725559374) q[8];
u3(0.890839447536141,-0.419978955244395,-5.25763643530921) q[6];
u3(1.94887677129037,0.608787730489005,-3.52136723141999) q[11];
u3(2.52923850348182,3.61345560962532,-2.16677375109947) q[4];
cx q[4],q[11];
u1(1.13708986935928) q[11];
u3(-1.29249242393455,0.0,0.0) q[4];
cx q[11],q[4];
u3(0.330282268095635,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.51273484910014,-3.81676311642282,1.49999445202072) q[11];
u3(0.781463436222503,0.355190002914282,-5.79747237085422) q[4];
u3(1.27747442225121,-1.32666472570224,-0.799979178644684) q[6];
u3(1.19658062994322,-2.34575917147478,0.0911737523857199) q[10];
cx q[10],q[6];
u1(2.07301309154632) q[6];
u3(0.0699743489837656,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.35480503181038,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.21988852335284,1.99353145606512,0.390040351010674) q[6];
u3(2.24198483476528,0.858371742896338,-3.59812185191561) q[10];
u3(0.515612800638294,0.975908083984871,-0.384678562610016) q[2];
u3(1.38634866099933,-0.185314132049724,-2.69702131243055) q[5];
cx q[5],q[2];
u1(1.44901341085570) q[2];
u3(-0.859559614389879,0.0,0.0) q[5];
cx q[2],q[5];
u3(-0.203458634929757,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.76593103251810,2.67471661142960,-3.37251222901582) q[2];
u3(1.23470712221794,-2.11351796554354,1.85374115664990) q[5];
u3(1.94681491274993,1.95900273901047,0.890016398065076) q[8];
u3(1.57862675420293,0.316183691485012,-3.53592851171464) q[9];
cx q[9],q[8];
u1(0.824802469542350) q[8];
u3(-0.280023169260197,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.86015917847454,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.540001359513471,3.59272070343714,-2.01224541802333) q[8];
u3(2.33883579472331,-1.19802997858558,0.357812329055646) q[9];
u3(0.609189654906331,1.68668357652121,-3.00383819890805) q[3];
u3(1.56754033655248,-2.47898073728575,3.62310653877025) q[12];
cx q[12],q[3];
u1(-0.0832593785422908) q[3];
u3(-2.05320817066064,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.46957954301682,0.0,0.0) q[12];
cx q[12],q[3];
u3(0.579322904643154,1.13314165382223,-3.07033637187225) q[3];
u3(0.839829100997579,1.04234093650860,-2.68905400504936) q[12];
u3(1.06284528158824,-0.335598000754834,-1.75628728091956) q[1];
u3(2.69632740822847,-3.06417838630826,1.92870749121017) q[0];
cx q[0],q[1];
u1(1.21486499508899) q[1];
u3(-0.425228891701167,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.11995051294662,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.27633570422773,2.29146761516483,-0.755363599600893) q[1];
u3(1.59189583439431,5.28545453039739,0.787955117764660) q[0];
u3(1.64311765336288,2.90952817631710,-1.76370193718081) q[4];
u3(1.08210401790311,0.869341634702641,-0.741254666330111) q[8];
cx q[8],q[4];
u1(2.33912053119984) q[4];
u3(-2.87613441988875,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.38779063117987,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.18771271146625,3.10177421552851,-1.87640590448390) q[4];
u3(0.547879982015979,-2.42329978390204,-0.911208644234173) q[8];
u3(2.23903020183556,-1.89531142466114,-0.390964885381846) q[10];
u3(2.09020312334481,-3.98257001064467,-0.771960167000262) q[12];
cx q[12],q[10];
u1(3.20547369030028) q[10];
u3(-1.29543020024526,0.0,0.0) q[12];
cx q[10],q[12];
u3(2.74393481381826,0.0,0.0) q[12];
cx q[12],q[10];
u3(1.12620531856679,-0.0427520108267048,0.229404562300327) q[10];
u3(1.71673662300863,1.35143430935202,-3.37678046040163) q[12];
u3(1.85779139674920,-1.47516071319758,-0.969519120561599) q[0];
u3(0.558326122151374,-3.84626524778224,0.320355479854105) q[6];
cx q[6],q[0];
u1(0.658819139798550) q[0];
u3(-1.21791998480944,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.89876184724528,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.35998962705956,-2.75537181173452,3.14222248927971) q[0];
u3(1.84450124032351,-0.301627215591001,-3.34587457485643) q[6];
u3(2.84885458738091,-2.63067694902275,0.388963059094153) q[11];
u3(2.79784310032894,-2.53035562326971,-2.03120039036767) q[9];
cx q[9],q[11];
u1(1.15413152119233) q[11];
u3(-3.43802222171373,0.0,0.0) q[9];
cx q[11],q[9];
u3(1.78833742713809,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.52733722291130,-2.03725879835129,3.21276947217321) q[11];
u3(2.63607449624432,-0.933004569009394,5.12233848524975) q[9];
u3(1.27381840054606,1.20589662598596,0.635771620427939) q[3];
u3(1.93015792672672,1.07106062840451,-2.48349284642542) q[2];
cx q[2],q[3];
u1(1.45088339557218) q[3];
u3(-0.430545212376464,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.45324017275897,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.87607451922059,2.09468178035089,0.494418499318003) q[3];
u3(2.26903840818467,-2.90498822248401,1.93789838052394) q[2];
u3(0.245489538382221,-2.31385655064306,2.14009926775109) q[5];
u3(0.678903525858886,-0.663263992595935,-1.29195491387951) q[7];
cx q[7],q[5];
u1(1.62044231548462) q[5];
u3(-3.34286398069836,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.69908372894085,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.927079762590194,0.0294478948743728,-2.65929899671728) q[5];
u3(1.35200024179083,-3.19718519326291,1.70023171193171) q[7];
u3(1.97071758012774,-0.604646053886425,-1.89886542172056) q[8];
u3(1.29768556459716,1.25330581035212,-4.73382257749907) q[3];
cx q[3],q[8];
u1(0.399367103322073) q[8];
u3(-1.18064712107444,0.0,0.0) q[3];
cx q[8],q[3];
u3(2.36807599123917,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.971815333581558,1.60341093206181,-1.64985737629129) q[8];
u3(0.777653571821472,-2.12213066840889,3.40465112801561) q[3];
u3(1.60430735804404,3.58336146424139,-2.03179628828366) q[4];
u3(0.952149417923558,2.82159396735390,-1.81360319643028) q[11];
cx q[11],q[4];
u1(3.36785233177917) q[4];
u3(-3.85141803501958,0.0,0.0) q[11];
cx q[4],q[11];
u3(-0.603227926181694,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.97638317997535,0.521756060774654,1.40023432468672) q[4];
u3(1.35839080287971,-3.55530784546286,0.464335723786749) q[11];
u3(2.47963952309444,-1.87139577747965,0.835458976493457) q[12];
u3(2.85681676482885,2.25819207681935,3.59289415996091) q[9];
cx q[9],q[12];
u1(1.40997186731377) q[12];
u3(-0.116797342397399,0.0,0.0) q[9];
cx q[12],q[9];
u3(2.12439115834373,0.0,0.0) q[9];
cx q[9],q[12];
u3(0.620926952725410,-0.368995010715443,0.407965527476937) q[12];
u3(2.40645135686633,-3.31029441037231,-0.953819814350587) q[9];
u3(1.06062879259071,1.01398517296295,-1.41989833804494) q[5];
u3(0.450378868049955,-1.81661912479147,0.217559349696419) q[2];
cx q[2],q[5];
u1(0.107033140706268) q[5];
u3(-1.51876137466496,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.47278968204702,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.550839491048444,0.170985857839401,2.24947532215768) q[5];
u3(1.29992842469886,0.401977504650279,-0.0770402020324451) q[2];
u3(1.36643889261538,2.31606448200640,-1.15017989951382) q[10];
u3(0.552517292993850,1.39157901148897,-0.0713480625984282) q[0];
cx q[0],q[10];
u1(-0.634830997098557) q[10];
u3(1.33538819082454,0.0,0.0) q[0];
cx q[10],q[0];
u3(3.89618642930024,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.94546696935168,3.88366075576924,-2.22149927750178) q[10];
u3(0.779001642536037,4.56181188532402,-1.56621407233513) q[0];
u3(0.668809651348728,-0.694437197221766,0.857186771411399) q[6];
u3(0.560237866206874,-1.61344462182642,0.535261393289319) q[1];
cx q[1],q[6];
u1(-0.231927667543821) q[6];
u3(-1.61148114016504,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.942810188507364,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.32969932251201,-1.37033076064677,-2.63936454885140) q[6];
u3(0.927080000131601,-1.44678257882103,1.98373193431601) q[1];
u3(0.672468701447347,-0.466221154305683,1.22519314537326) q[5];
u3(1.51850896776907,-1.47586185394745,-2.57666754041129) q[2];
cx q[2],q[5];
u1(2.21635220165252) q[5];
u3(-2.57479507681847,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.26241943579895,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.578778713863869,1.45308578418738,-2.05166963440137) q[5];
u3(0.933405580545906,-0.0272458984238475,5.34826993892906) q[2];
u3(1.98968845715789,0.770248324709039,-2.02691664642283) q[12];
u3(2.86581477857790,1.37806159900543,-4.73322027584388) q[4];
cx q[4],q[12];
u1(1.14297748480193) q[12];
u3(-0.223545165039878,0.0,0.0) q[4];
cx q[12],q[4];
u3(2.06677424736222,0.0,0.0) q[4];
cx q[4],q[12];
u3(2.51053683816869,3.54819767474505,-0.224641432575461) q[12];
u3(1.34972852947388,4.05934227282186,1.34921061867192) q[4];
u3(2.94934680167661,-2.60239474602656,3.10897492189420) q[3];
u3(1.48210790019956,3.88403024557932,-2.35546798061232) q[10];
cx q[10],q[3];
u1(1.24358285921014) q[3];
u3(-3.45614536608473,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.28826282626260,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.952718356057534,2.22041705044543,-2.80654062479960) q[3];
u3(0.781772452528792,5.35932224284993,0.681309056710806) q[10];
u3(1.08476685857752,0.288549492089853,2.23651204712611) q[1];
u3(1.48952822733107,-2.40324054817120,-1.70050524257096) q[0];
cx q[0],q[1];
u1(-0.219523831241086) q[1];
u3(-1.51939182721288,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.866624284960961,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.69113234875108,5.26251007070886,-0.815811037594713) q[1];
u3(1.87420965751719,-0.167362651510974,6.07343901754436) q[0];
u3(1.83368575837434,0.0829045736513634,2.69037829950735) q[6];
u3(2.87564557236086,-2.94092935097100,-2.08506084010027) q[7];
cx q[7],q[6];
u1(-0.439757319105222) q[6];
u3(1.13482040684290,0.0,0.0) q[7];
cx q[6],q[7];
u3(3.65163310598173,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.44959563364437,2.68608830850686,-0.651928385009810) q[6];
u3(1.90905060483798,-2.27018970336847,3.88183029518725) q[7];
u3(1.49581015310385,-0.462586733666371,-1.01973643714323) q[8];
u3(2.38098892354952,-5.09440812404785,0.647068952394650) q[9];
cx q[9],q[8];
u1(3.14631777782597) q[8];
u3(-1.73113757950808,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.414763310701771,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.65518714623221,-2.97431521602708,-0.0747884148026381) q[8];
u3(2.22737222542118,-2.82818994092359,0.845745331079671) q[9];
u3(1.04697073346299,0.0335625871707282,-1.06806929775790) q[2];
u3(1.97508613002656,1.63062116841899,-4.15452422061580) q[0];
cx q[0],q[2];
u1(4.18678484586444) q[2];
u3(-3.48175728521070,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.577494535991104,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.24339875836353,1.25889188481455,0.639994540578326) q[2];
u3(2.91850368210333,-0.727417556216234,-1.58192819005530) q[0];
u3(1.75975293552408,-0.226651272578853,-0.364787192017966) q[1];
u3(0.332620018180955,0.481665199695430,-4.91526118999027) q[8];
cx q[8],q[1];
u1(1.77390303334604) q[1];
u3(0.324309548685812,0.0,0.0) q[8];
cx q[1],q[8];
u3(0.985176772114906,0.0,0.0) q[8];
cx q[8],q[1];
u3(2.48156734288484,-1.63926612216156,0.553935602286690) q[1];
u3(1.83786487267469,1.83937838848280,-3.05291949939619) q[8];
u3(0.829317546128607,0.721649070004880,-3.44819752707217) q[11];
u3(2.24807811963090,-3.03119654665912,3.21896728801908) q[10];
cx q[10],q[11];
u1(0.0839950534622376) q[11];
u3(-1.74989484761753,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.10664949466445,0.0,0.0) q[10];
cx q[10],q[11];
u3(0.519779563740181,2.44138435697559,-0.824875977930767) q[11];
u3(1.76214337865716,1.73484896580049,-1.59197753664859) q[10];
u3(1.13206192035647,0.310474003338107,2.56488424830146) q[7];
u3(1.70742308513221,-2.58860890293386,-1.63305756196855) q[4];
cx q[4],q[7];
u1(3.35415433038148) q[7];
u3(-0.622766669662961,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.29730009978443,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.18622150634426,4.01435328760269,-1.84949796465516) q[7];
u3(1.59165371626823,-3.25325649234157,-1.14675984878589) q[4];
u3(2.72215561611237,1.45997855296859,-2.71629513106474) q[12];
u3(2.65963611869978,0.229656331355364,-4.63341849717576) q[5];
cx q[5],q[12];
u1(2.35244438186777) q[12];
u3(-2.72667453293634,0.0,0.0) q[5];
cx q[12],q[5];
u3(1.38763806862856,0.0,0.0) q[5];
cx q[5],q[12];
u3(1.05214300695261,3.11184992146787,-3.16658143794605) q[12];
u3(0.134370605238746,-1.40403600400898,-0.859182380085781) q[5];
u3(0.683146995605357,-0.481207031662250,0.595235978757059) q[6];
u3(1.97349893924676,-1.29373639061715,-1.70095453970665) q[9];
cx q[9],q[6];
u1(1.57050311350694) q[6];
u3(-0.445459673319101,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.44112036703785,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.746006137054906,0.587271950877883,-3.30179598095585) q[6];
u3(2.43561674924574,-0.0138210933585601,0.0474317181144277) q[9];
u3(1.79635244090027,1.66811875904462,-2.10411890504602) q[1];
u3(2.53992777415788,5.64028803804646,0.230308845479676) q[6];
cx q[6],q[1];
u1(2.47422823686153) q[1];
u3(-1.68392346289221,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.129903700370669,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.42825011520268,-1.32294587671429,-1.33368963818669) q[1];
u3(0.687876890755211,3.85234844165296,2.01322840085978) q[6];
u3(3.06270835727713,1.36077769612217,1.09016014613206) q[2];
u3(0.942516505658679,0.511399558572289,-4.25037448230084) q[11];
cx q[11],q[2];
u1(-0.921667515480160) q[2];
u3(0.273930840675475,0.0,0.0) q[11];
cx q[2],q[11];
u3(3.47878892814553,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.50704114523837,0.656104081267000,2.89816966824401) q[2];
u3(2.10491258134380,3.61491706567664,0.613816471012165) q[11];
u3(1.51254167540422,1.27039486128324,-3.26574106863227) q[7];
u3(2.00045658762392,3.38052214790402,-2.87802956387824) q[12];
cx q[12],q[7];
u1(-0.355022286974023) q[7];
u3(-2.09356905266389,0.0,0.0) q[12];
cx q[7],q[12];
u3(1.40768770156405,0.0,0.0) q[12];
cx q[12],q[7];
u3(0.974206875716861,0.647177476848129,1.52863419189644) q[7];
u3(1.40006833749467,-1.46017100129546,2.75715653243683) q[12];
u3(0.695287634947258,-2.17807552566299,2.81687689753571) q[10];
u3(0.532013403774745,1.82400769802267,-3.62039129414196) q[9];
cx q[9],q[10];
u1(1.80490866014041) q[10];
u3(-0.621699390572063,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.106921370437174,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.87177013050854,-0.846925147483388,-1.62186964982766) q[10];
u3(3.10507182714282,2.35502162956922,2.25682933532683) q[9];
u3(0.532441494615955,0.264630480435367,-1.18627783684546) q[3];
u3(0.989689875798411,0.0502060339116281,-1.37347976682659) q[8];
cx q[8],q[3];
u1(1.55375146070799) q[3];
u3(-0.872772387505209,0.0,0.0) q[8];
cx q[3],q[8];
u3(-0.420950379618833,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.04519933608995,-0.975326585686825,-0.584225314964426) q[3];
u3(2.00390694377731,-1.02084670587182,-2.66342693990906) q[8];
u3(1.45135935473916,0.0934048542415542,2.12466523708464) q[5];
u3(1.64272249964366,-0.918214822528402,-1.15463852148197) q[4];
cx q[4],q[5];
u1(0.771361324180429) q[5];
u3(-0.129344626781854,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.21389667151180,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.48278547923202,-2.89150474095770,2.49156022126605) q[5];
u3(0.793038014588203,2.43789711694828,-1.75888392785346) q[4];
u3(1.39228057581623,0.156775784498641,-1.94022618104037) q[5];
u3(1.32859165861722,-3.56806694428833,1.88460541600408) q[6];
cx q[6],q[5];
u1(1.79626213156757) q[5];
u3(-0.443415224902214,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.44617294074546,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.27953836466489,1.54883465151843,-1.57285203094758) q[5];
u3(2.96195088376332,-0.346634684933974,4.76896734353011) q[6];
u3(2.07924444270215,-0.305633960444939,0.258271013277361) q[10];
u3(1.72484196700356,-3.18219541134829,-0.883529291508016) q[3];
cx q[3],q[10];
u1(1.94466824254948) q[10];
u3(-2.08397853783578,0.0,0.0) q[3];
cx q[10],q[3];
u3(-0.152629250264707,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.680052721273186,-1.48202992074247,1.27659349305198) q[10];
u3(1.75119675470958,5.75329353643237,-0.144226485650495) q[3];
u3(1.03409872277246,2.67705607389882,-0.904311509531441) q[0];
u3(2.16329580472176,1.23093096465821,-2.68980619474168) q[2];
cx q[2],q[0];
u1(1.31058741238840) q[0];
u3(-3.51986288033861,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.39833128219000,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.63333824628305,-0.858950870193044,-1.54740787611340) q[0];
u3(1.72522365336063,-2.10597673864503,1.73616683146362) q[2];
u3(0.916074672493646,2.50991109700951,-2.88072023331831) q[8];
u3(0.947605045966130,1.93306710815339,-3.95560498545052) q[9];
cx q[9],q[8];
u1(1.38743036550616) q[8];
u3(-0.609236798909697,0.0,0.0) q[9];
cx q[8],q[9];
u3(3.06219712249541,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.65401497640103,0.499661171307040,0.172531873775365) q[8];
u3(1.37355446538620,-0.983634478330929,3.48051124574399) q[9];
u3(2.04220622975684,1.90756824823488,0.154437150655828) q[12];
u3(2.37222636777229,0.148862254515700,-3.74850953384731) q[11];
cx q[11],q[12];
u1(1.45424890822576) q[12];
u3(-1.05682980833345,0.0,0.0) q[11];
cx q[12],q[11];
u3(-0.179569038339271,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.70130557257366,3.86752473015818,-1.90159508869715) q[12];
u3(1.50579291419149,-3.41873522244222,2.71534127866061) q[11];
u3(1.52402820425694,2.19860338910901,-3.87107837037683) q[4];
u3(1.68128228725201,-2.70838659397349,3.50087113141041) q[1];
cx q[1],q[4];
u1(1.09577662959781) q[4];
u3(-3.09936681692853,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.79726939829098,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.03421170505873,-0.673655281217790,2.74142422053132) q[4];
u3(1.20143741876331,3.93444651002489,0.763139799041410) q[1];
u3(1.57859943138816,0.208681276527796,-0.673998921200462) q[8];
u3(2.36607570568763,-5.21191436575901,0.645425239409411) q[2];
cx q[2],q[8];
u1(0.574023536500569) q[8];
u3(-1.05832002710315,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.18959841329216,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.28193067423830,2.65929949428675,-2.76991091289516) q[8];
u3(1.18887063201602,2.42480617937929,2.77858120450291) q[2];
u3(1.60504888873828,1.59346153269436,-0.493724370844830) q[3];
u3(1.26395523134681,1.41868266537325,-3.95003382170354) q[4];
cx q[4],q[3];
u1(2.02781569947972) q[3];
u3(-2.91589281370910,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.841034104356095,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.855613927866389,0.651983981831943,1.27911458713402) q[3];
u3(1.41599771809615,-5.90506395250407,-0.150899717634768) q[4];
u3(1.24177992071267,0.109336211980917,1.98725011087195) q[0];
u3(1.39461760770224,-0.505775476133039,-1.41868682741185) q[9];
cx q[9],q[0];
u1(1.39080464672313) q[0];
u3(0.212841871809147,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.61758239625429,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.25867095818969,-0.884468387770990,2.87908300951489) q[0];
u3(1.73643730623108,0.410028012427032,3.74093844209853) q[9];
u3(0.164530379071007,2.57350889530775,-2.67697410323742) q[1];
u3(0.236382192135704,-3.63863105148076,1.67514100296270) q[10];
cx q[10],q[1];
u1(1.35213259826680) q[1];
u3(-2.47069618798179,0.0,0.0) q[10];
cx q[1],q[10];
u3(3.02325814298322,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.618842438474004,-2.42994617887627,2.66041223259964) q[1];
u3(0.991739297713522,-2.08271078627097,-3.23845512775990) q[10];
u3(1.36002713398242,2.88512398895470,-0.553059909286877) q[6];
u3(1.30069689582816,1.22441133355160,-1.62355545512957) q[7];
cx q[7],q[6];
u1(1.37839403343301) q[6];
u3(-0.803024418883102,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.57971472170261,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.69215429999189,-0.272948684455507,0.0317831382020836) q[6];
u3(1.13373122821030,-2.09308994360053,3.45628168583300) q[7];
u3(0.957999288674160,-1.32387769914343,1.06675720845813) q[5];
u3(1.46780103310860,-2.12741182016539,0.453829148709253) q[12];
cx q[12],q[5];
u1(1.93411635617726) q[5];
u3(-2.50057112692704,0.0,0.0) q[12];
cx q[5],q[12];
u3(0.182671450626316,0.0,0.0) q[12];
cx q[12],q[5];
u3(2.78346063585785,-2.80756779358263,-1.03346848963722) q[5];
u3(1.99458618201069,-4.05376336207367,0.817416439003623) q[12];
u3(2.08033470002425,-2.96861421733709,2.07625019165640) q[7];
u3(1.75917917774325,-3.07504274074892,2.04723150622184) q[1];
cx q[1],q[7];
u1(1.96038745318228) q[7];
u3(-2.38114883645207,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.949539965092608,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.84266235153643,2.99235466489766,-3.25336249200097) q[7];
u3(0.792668113527204,0.139296324426421,3.11880690131247) q[1];
u3(1.71460793671287,-1.35427199513721,-0.534372285119596) q[11];
u3(2.16248632793264,-1.94800162846850,1.40753355804050) q[6];
cx q[6],q[11];
u1(1.69729966982777) q[11];
u3(-3.28034014896305,0.0,0.0) q[6];
cx q[11],q[6];
u3(1.11417651775350,0.0,0.0) q[6];
cx q[6],q[11];
u3(2.02338724779385,-2.13362503413890,-0.959809820817943) q[11];
u3(1.54403833861749,-3.86917445778502,-0.546723877747207) q[6];
u3(3.09314791601505,2.33448609200823,-0.786484455908628) q[12];
u3(1.92099501948820,1.06942820515769,-4.09117630615124) q[0];
cx q[0],q[12];
u1(2.36386008557462) q[12];
u3(-1.73548421739165,0.0,0.0) q[0];
cx q[12],q[0];
u3(3.55921367780357,0.0,0.0) q[0];
cx q[0],q[12];
u3(2.00825499521554,-2.10010555235468,3.51797047814751) q[12];
u3(0.949110625011924,-2.13957610647736,3.78336127929575) q[0];
u3(2.26236803141986,-3.47869873921528,2.79442266536751) q[3];
u3(1.19901248101141,3.21068354117272,-2.01760210229005) q[4];
cx q[4],q[3];
u1(2.30436432940490) q[3];
u3(-3.09658876704575,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.51385486168891,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.49809669213655,1.62391014746104,2.64694072248997) q[3];
u3(0.315828431075110,-2.40287938857387,-1.85089394359599) q[4];
u3(1.25420551584286,2.55121033415821,-0.860781284613804) q[8];
u3(1.70067712987528,1.36978826064347,-0.761653888687836) q[2];
cx q[2],q[8];
u1(3.71299383889128) q[8];
u3(-4.40339345256183,0.0,0.0) q[2];
cx q[8],q[2];
u3(-0.754704996806645,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.18552076684593,1.07818174177764,-1.13504182246853) q[8];
u3(2.63250673400468,-0.240755937896426,2.69873360959592) q[2];
u3(0.647496199875663,1.50734742960421,-1.12984436138853) q[9];
u3(0.0904362416808176,-1.40802780210935,0.472989309791316) q[5];
cx q[5],q[9];
u1(1.53092565810176) q[9];
u3(-2.66708259911209,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.04136840845691,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.35247797111727,-1.64362998638217,-0.318285663920627) q[9];
u3(1.38251917087535,5.08028337091169,1.19899807386167) q[5];
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
