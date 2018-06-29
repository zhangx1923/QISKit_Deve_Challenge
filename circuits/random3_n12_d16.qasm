OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.16671970261739,-3.57541708698993,1.38771452899409) q[7];
u3(1.17518167218370,-0.501185272918311,1.53108572727068) q[9];
cx q[9],q[7];
u1(-0.862349528115798) q[7];
u3(0.188881986198106,0.0,0.0) q[9];
cx q[7],q[9];
u3(3.36686812025147,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.449747802311396,3.30085384080054,-2.66963326675536) q[7];
u3(1.17407231895637,3.53627269567267,-2.47292205064710) q[9];
u3(2.19009848905773,-1.58799503900645,1.07653813870362) q[3];
u3(2.60672490701322,2.52820703217253,3.20781700515791) q[0];
cx q[0],q[3];
u1(2.99906190995428) q[3];
u3(-2.37914495687752,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.47405374274018,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.138807983110285,2.62043729889597,-1.29265413921476) q[3];
u3(1.32248179493423,-3.20192036546454,2.99744333179538) q[0];
u3(2.08536766601597,1.78252070994187,0.178046523815691) q[6];
u3(2.43888598266928,0.798646265381507,-3.04200198670511) q[10];
cx q[10],q[6];
u1(3.21962116579452) q[6];
u3(-2.35328135736628,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.24064778026209,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.842443729755016,-1.91356995785691,0.765838644556608) q[6];
u3(2.30108407448670,1.37241902416936,-4.01772950877639) q[10];
u3(1.59190935108885,1.73743255774666,-0.595663129303592) q[4];
u3(0.149416982197816,-0.358550173378811,-1.88992192075546) q[2];
cx q[2],q[4];
u1(0.794283795269736) q[4];
u3(-1.34762289876261,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.94281895404750,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.00409212876526,-0.0486494204860335,0.469113065014915) q[4];
u3(1.35198979423335,-4.19349003872714,-1.54833935919067) q[2];
u3(1.73227970717872,0.217736163031186,1.10058068462725) q[11];
u3(2.26200594190491,-1.64756279392123,-2.38862525857781) q[1];
cx q[1],q[11];
u1(1.56086359115505) q[11];
u3(0.213754459966785,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.10861751922395,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.16799044945204,-0.730875196801341,-1.23590153926481) q[11];
u3(0.995394424028173,2.56968497267199,-2.51794209154620) q[1];
u3(1.79960615656034,-0.361209266713342,-1.75442244940285) q[5];
u3(1.60073679973933,1.03217234861552,-4.48942431447659) q[8];
cx q[8],q[5];
u1(2.96840299701953) q[5];
u3(-2.33716379343684,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.13550424802168,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.88022275170187,-0.0578895876552791,1.94517367763153) q[5];
u3(0.755607324137755,-2.45344620958348,-0.263315973141245) q[8];
u3(1.46550160995207,2.82027871341166,-1.92385097288947) q[0];
u3(2.01236643951022,2.40982122814495,-0.434102498537177) q[10];
cx q[10],q[0];
u1(3.19338985746251) q[0];
u3(-1.36065953975150,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.81327077593582,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.62844440087630,-1.39993697866705,0.0826597783402729) q[0];
u3(2.93686539868907,1.33729060462689,0.155033139669947) q[10];
u3(0.918689677337013,-1.18153949835231,1.82310863168573) q[7];
u3(1.42865678658967,-1.52676135397284,-2.29703918413436) q[1];
cx q[1],q[7];
u1(1.84529619476031) q[7];
u3(-2.63750896815513,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.84643264911420,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.72102541403016,-3.50439199635521,1.25813105548527) q[7];
u3(2.62865990914215,-6.01066655237580,0.0639476133337169) q[1];
u3(1.26381506998451,0.460420034872945,-2.43500411204888) q[6];
u3(2.54073646054330,-4.77130345634011,0.996125944139420) q[3];
cx q[3],q[6];
u1(0.0927833500526176) q[6];
u3(-0.902757017017797,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.45295890683092,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.91323544590026,-0.306839861245278,1.29882097914237) q[6];
u3(2.10643092412528,0.183722907711472,-1.53800616911101) q[3];
u3(1.65964701779158,2.33784107800822,-3.69409977476196) q[11];
u3(1.00493483694522,2.51530526926942,-1.79867485777946) q[4];
cx q[4],q[11];
u1(1.89419567853682) q[11];
u3(-3.05195191821137,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.09749771517799,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.54201396679295,-2.25187158903815,2.32091005950868) q[11];
u3(1.23402399516603,1.54279509980219,-3.47699914139568) q[4];
u3(1.97232345966184,-1.77786243844520,0.701689930777610) q[9];
u3(2.50851463743901,-1.50904324725562,-0.265912752649020) q[5];
cx q[5],q[9];
u1(2.13129688712367) q[9];
u3(0.0598335576150748,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.23419718608571,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.03770397486243,0.0107058084284902,-1.00092288512165) q[9];
u3(2.97034038178779,0.681429928015337,0.145150303462984) q[5];
u3(2.39144795968676,-2.48573850121394,3.68991036352170) q[8];
u3(0.719689945562641,-0.278405155152257,1.60333959443149) q[2];
cx q[2],q[8];
u1(1.51293093201058) q[8];
u3(-0.720188841976436,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.00022619925163,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.830185782007207,1.68411246620225,-1.07651868287245) q[8];
u3(2.26759544790636,4.91181765358187,0.636207641350645) q[2];
u3(1.68013446619432,-1.11603494296542,0.205452279525622) q[7];
u3(1.75865668688816,-3.16572765141768,-1.09668745683079) q[3];
cx q[3],q[7];
u1(0.0299338048527873) q[7];
u3(-1.21057321787896,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.75307850525489,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.11798974312782,-0.0884024721202246,3.30817169572365) q[7];
u3(1.16845094410635,-1.83605423435187,4.11407585796996) q[3];
u3(1.39964850706202,1.18779136734728,-1.88445927526828) q[10];
u3(1.67446361087887,1.49000652481124,-4.65148555581653) q[2];
cx q[2],q[10];
u1(3.09213173281350) q[10];
u3(-2.66296929850664,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.35538757969042,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.92560334068580,0.332196906958058,1.04877438648221) q[10];
u3(0.922001079096214,0.418434205251433,5.12451361365661) q[2];
u3(2.24530353600530,-1.33205928856036,-1.61231591631405) q[0];
u3(1.04598567927409,-1.57465758166821,-3.64414714069279) q[5];
cx q[5],q[0];
u1(1.44711105085387) q[0];
u3(0.0535394519113142,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.91189776227583,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.11051170576174,-3.25249035077579,2.35827709547451) q[0];
u3(2.38344408969256,-3.47088891398762,2.73869965212926) q[5];
u3(1.17520404868916,0.208152244432764,0.802831920350482) q[1];
u3(1.93992494105251,-0.573551074432056,-2.27173594603487) q[6];
cx q[6],q[1];
u1(-0.538065231193204) q[1];
u3(0.149672288811871,0.0,0.0) q[6];
cx q[1],q[6];
u3(4.29894894077981,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.55214288003343,-3.22539317072394,0.0557746760110691) q[1];
u3(1.24855506341971,-1.55883066445210,4.42108657052189) q[6];
u3(1.22101540760977,1.17449545358604,-1.24077827269774) q[11];
u3(1.25547346439195,-0.804700913089276,-3.15390658432887) q[4];
cx q[4],q[11];
u1(0.287131172135158) q[11];
u3(-1.31872110186594,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.28269295148688,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.16098167943514,1.02782273305907,0.0732929576825506) q[11];
u3(2.01447157389751,-2.35992253920957,2.85347387818169) q[4];
u3(1.93834949823172,3.58003115459382,-1.13951178302083) q[8];
u3(0.642501857579483,1.80014523084358,-2.88171472469473) q[9];
cx q[9],q[8];
u1(-0.281074376647961) q[8];
u3(-1.65500623585210,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.07979496025419,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.31586400092905,-2.48402461840619,2.28922789456823) q[8];
u3(1.75352532472813,6.14327891175867,-0.0845030466184324) q[9];
u3(1.68041451769216,-0.245026875257985,0.0815591978535384) q[11];
u3(0.701289854779535,-2.28196647387234,-1.30951359829938) q[8];
cx q[8],q[11];
u1(1.78500401380471) q[11];
u3(-2.57499116179303,0.0,0.0) q[8];
cx q[11],q[8];
u3(0.985500172132154,0.0,0.0) q[8];
cx q[8],q[11];
u3(0.472773118777643,0.339012450887098,-3.22603546984967) q[11];
u3(1.35323077518177,-5.14885861414896,-1.06716221965970) q[8];
u3(1.46537876295521,1.30863277235235,-3.55736577857889) q[7];
u3(0.787901505955728,1.90439448514855,-1.87781276015240) q[2];
cx q[2],q[7];
u1(0.00281475021277888) q[7];
u3(-0.890081369298882,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.88719659778446,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.636026643092978,-1.57267683847931,4.68170447769065) q[7];
u3(2.33931269080382,-4.71558012812501,-0.249027991222680) q[2];
u3(1.88214769824446,0.0562554332892595,2.15409506318661) q[5];
u3(3.07429568984266,-0.364364731149234,0.883277163796308) q[10];
cx q[10],q[5];
u1(1.61733762037980) q[5];
u3(0.276520457241729,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.501139580053304,0.0,0.0) q[10];
cx q[10],q[5];
u3(2.50827926734923,0.905762662554335,2.51276712162132) q[5];
u3(1.27689771898598,-0.215380083420980,-3.93702557234704) q[10];
u3(2.18601906171284,0.472290891736653,2.64454592782328) q[3];
u3(1.19794755573406,-0.584199358853267,-1.36778294372399) q[0];
cx q[0],q[3];
u1(1.08109844242032) q[3];
u3(-0.101142435691224,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.53212070046844,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.77862143292307,-2.02167910685679,3.33547339295102) q[3];
u3(2.18868816377990,2.16715300998570,1.87857842476380) q[0];
u3(2.64834003685958,-1.75472579447152,-0.997027704767693) q[9];
u3(0.996583060227053,-1.77516043739544,-3.48219172397498) q[1];
cx q[1],q[9];
u1(2.36777618036392) q[9];
u3(-1.89974888353908,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.50640060153979,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.96045575999036,-1.71525912591159,1.36587054071199) q[9];
u3(0.913179078115355,-0.551511836020808,-4.51436859269442) q[1];
u3(2.57428526880503,0.0760457169905182,-0.472630475623758) q[6];
u3(1.09869039507197,0.216732606188617,-5.29294106872154) q[4];
cx q[4],q[6];
u1(0.154889780288489) q[6];
u3(-0.508218943469633,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.23030930655081,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.32605600983135,4.20345978172494,-1.26737399671951) q[6];
u3(2.20395608653413,0.503633633204628,4.94551917580988) q[4];
u3(1.27508971215500,1.48335306427863,0.102325553889317) q[6];
u3(0.436218415883533,0.212892942276185,-4.09961078843248) q[9];
cx q[9],q[6];
u1(-0.182964403060948) q[6];
u3(-1.70130525525062,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.914566774230742,0.0,0.0) q[9];
cx q[9],q[6];
u3(1.49932675525244,2.32568063575848,0.0947417963463366) q[6];
u3(1.00098674592650,-2.13286492859025,4.10733011311929) q[9];
u3(1.81689927619486,0.959103548974816,1.07864209660455) q[10];
u3(1.15301138311629,-1.14890452710537,-2.22492934644829) q[4];
cx q[4],q[10];
u1(1.63527511091377) q[10];
u3(-2.40692341943977,0.0,0.0) q[4];
cx q[10],q[4];
u3(-0.0339783945190000,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.78000315074431,1.80066051276400,-3.39369810529304) q[10];
u3(2.37022342089408,3.88002051029846,0.0544457570093930) q[4];
u3(0.766977125067993,3.73297523226140,-1.56431859659812) q[2];
u3(1.44692367933320,0.757304656404556,-1.80519507151903) q[0];
cx q[0],q[2];
u1(1.35226193501704) q[2];
u3(-0.536497473612671,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.34405202594559,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.93322723971154,-0.698959542050851,1.64002423420675) q[2];
u3(1.12388686284826,-1.40005678234748,-3.82896702442083) q[0];
u3(1.95690956827837,1.95542898451315,1.04701828929408) q[1];
u3(1.83487053009118,-0.148824877589018,-3.28884916842594) q[8];
cx q[8],q[1];
u1(1.89472821611491) q[1];
u3(-2.45020213569298,0.0,0.0) q[8];
cx q[1],q[8];
u3(3.12881839571500,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.56162189854454,-1.84292121004084,0.682194845614186) q[1];
u3(2.19201804812702,-5.11036611610110,-0.608273476850844) q[8];
u3(2.88933528172539,1.52893100907473,-1.62141203421352) q[3];
u3(2.78712276591425,1.59377035770950,-4.63115413293013) q[7];
cx q[7],q[3];
u1(2.51381290627722) q[3];
u3(-1.71705523028830,0.0,0.0) q[7];
cx q[3],q[7];
u3(3.30647638176965,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.86544320614823,-2.03423561098672,1.11596838234797) q[3];
u3(1.35683065534111,-2.95509936060169,-1.38736729159054) q[7];
u3(1.69138557519364,-1.50608649118385,0.718137434414729) q[5];
u3(0.535595437650516,-2.21144236225928,-0.382835862528945) q[11];
cx q[11],q[5];
u1(4.25465997266601) q[5];
u3(-3.15039686075728,0.0,0.0) q[11];
cx q[5],q[11];
u3(-0.432636387630129,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.48504760193270,-0.890926096778236,3.78821622943311) q[5];
u3(1.48933933568054,4.65160466614987,1.51471288492791) q[11];
u3(2.11567576596190,3.28676883567447,-0.198120385836594) q[7];
u3(1.94590345791339,2.13787618442584,-1.16945408988980) q[9];
cx q[9],q[7];
u1(1.89318145945611) q[7];
u3(-2.27044276983446,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.495706707945694,0.0,0.0) q[9];
cx q[9],q[7];
u3(2.60640015850527,-0.908554593378256,0.326102003402525) q[7];
u3(2.75476873577718,-2.48366415940233,-3.37802459189553) q[9];
u3(1.08627915312089,2.58445742970978,-1.36091429041865) q[0];
u3(1.95903703376876,1.84220818932705,-1.07134133537110) q[4];
cx q[4],q[0];
u1(1.41759847116980) q[0];
u3(-2.56224793317035,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.50458454901603,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.06725988972755,-0.679678770142405,0.763022944598390) q[0];
u3(1.81406472779768,1.50322313754700,3.64483632065628) q[4];
u3(1.32721782952783,2.42317493713406,-3.38908402036946) q[11];
u3(1.85679507968491,3.34525036472364,-2.56939911952159) q[10];
cx q[10],q[11];
u1(2.34204215974720) q[11];
u3(-1.56141893347014,0.0,0.0) q[10];
cx q[11],q[10];
u3(3.33877577447426,0.0,0.0) q[10];
cx q[10],q[11];
u3(0.579692138132977,1.51196535151854,-2.02883883973013) q[11];
u3(0.603954328460648,2.54351593766536,3.01091876942435) q[10];
u3(2.00080310402272,-1.26410677583854,-1.46904881330433) q[8];
u3(0.628405507783745,-3.79625617173335,0.0375488964436319) q[5];
cx q[5],q[8];
u1(-1.39622591948677) q[8];
u3(0.562274262224493,0.0,0.0) q[5];
cx q[8],q[5];
u3(3.84259335111138,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.10321323880563,0.0855609107300178,1.64586539661242) q[8];
u3(1.13731820386690,1.93436170460051,-3.75965246057258) q[5];
u3(2.46176557372527,1.84770715555652,-3.91622413146271) q[6];
u3(1.28863850953159,1.60987541568769,-0.811133886432247) q[3];
cx q[3],q[6];
u1(3.07020554930731) q[6];
u3(-1.73602070969975,0.0,0.0) q[3];
cx q[6],q[3];
u3(0.765084346021329,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.22670310511107,-0.0275841472552261,-2.07464975845963) q[6];
u3(1.29677756715416,0.111491988221536,2.82767735135800) q[3];
u3(1.38288323612530,0.566196755312945,-1.20285494955111) q[1];
u3(0.769882620479157,-4.65481871312993,1.45875786021983) q[2];
cx q[2],q[1];
u1(1.92021531289567) q[1];
u3(-2.42860675029268,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.888543649534507,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.86230263290190,-1.45355174819792,-0.205793871903611) q[1];
u3(1.84048967838276,1.32097653958185,-4.50423782306777) q[2];
u3(2.37740271543036,-3.45554042735611,2.71327294545393) q[9];
u3(1.34761155950117,3.43517834120979,-2.60283765164119) q[3];
cx q[3],q[9];
u1(2.53790274966288) q[9];
u3(-1.54141664845680,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.157354027598642,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.02828068728716,-0.739674416823013,1.77810997428696) q[9];
u3(1.71941181852768,2.78750814216939,-2.60514781001731) q[3];
u3(1.50366824211020,1.30902430046861,-0.236927937061425) q[0];
u3(1.73100862357074,-0.717132635840491,-4.37359759370746) q[4];
cx q[4],q[0];
u1(2.30322395333619) q[0];
u3(-1.91301537570245,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.154083555654009,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.30465123820745,4.19525037609787,-1.91862811659512) q[0];
u3(0.921498512026135,-2.79153753388001,-3.26666218716551) q[4];
u3(1.63508756303178,-1.50124701696707,0.532135671876602) q[7];
u3(2.07084848576096,-3.72293188241527,0.0743125665075841) q[8];
cx q[8],q[7];
u1(0.647204768918450) q[7];
u3(-1.51011399942641,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.45468617312733,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.07170245089936,0.604780100547811,-3.47484740749825) q[7];
u3(0.696251998954737,-0.532262468160138,3.68319531250808) q[8];
u3(1.93426344557917,2.67990466511704,-0.725796433172361) q[6];
u3(2.99100342004788,5.09581542684904,-0.0244028353445529) q[5];
cx q[5],q[6];
u1(0.448662142116509) q[6];
u3(-1.51305035375160,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.90111036472991,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.16636731733033,3.99550097884867,-0.694415177012972) q[6];
u3(1.37266747821910,-0.0914331281469182,-1.72666157798771) q[5];
u3(2.46232772362208,0.908369263894798,2.11646294422598) q[10];
u3(1.09144742110446,-2.15250956727309,-3.18830477404341) q[1];
cx q[1],q[10];
u1(2.78820230192856) q[10];
u3(-2.22680047612684,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.03455784970505,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.60319720695783,-0.0241950010154472,0.803135272497446) q[10];
u3(2.60063687767299,-1.37755919150415,4.57446017627573) q[1];
u3(1.12984745191541,3.46128965639269,-2.14448503891461) q[11];
u3(1.37140776685166,1.13904160956991,-2.20024012342033) q[2];
cx q[2],q[11];
u1(1.75108828065265) q[11];
u3(0.184559515595984,0.0,0.0) q[2];
cx q[11],q[2];
u3(0.771579490455821,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.499960823115814,0.169217924578063,-1.71392770649536) q[11];
u3(0.559515210336538,0.281749418517415,-3.79419699838608) q[2];
u3(1.33465261342768,0.194288902218464,-1.92106155593139) q[11];
u3(1.77900293739533,3.19095788250094,-2.90754331158650) q[8];
cx q[8],q[11];
u1(-0.168678978822143) q[11];
u3(-2.12763193900016,0.0,0.0) q[8];
cx q[11],q[8];
u3(1.19589971063842,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.22537372498762,-3.48464707308123,2.76820295089663) q[11];
u3(1.47405613888075,-0.726160298131066,3.30810095201655) q[8];
u3(1.82805077863037,3.49545835594721,-0.422019329971992) q[10];
u3(2.51604768765508,2.40062739234105,-1.00055883571345) q[0];
cx q[0],q[10];
u1(1.99583524890819) q[10];
u3(-2.82820962115796,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.47452070078992,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.11353358205918,2.10280733487388,-1.04431593430639) q[10];
u3(0.758515260013446,3.12334453778237,1.66552081833830) q[0];
u3(2.23569566683923,-1.17552251194816,0.400002391010682) q[3];
u3(1.90019045029179,-2.17605955761758,-0.386510058653806) q[6];
cx q[6],q[3];
u1(2.86300373159427) q[3];
u3(-1.41887837828082,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.486321020373513,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.578790826731211,1.84687587328992,-4.12161078514040) q[3];
u3(1.86580436612897,0.730130750719093,-0.553884906639929) q[6];
u3(1.20221645716781,0.416071348446713,-0.854008089163236) q[4];
u3(0.834798342411723,-3.45790954343513,0.497897465147486) q[2];
cx q[2],q[4];
u1(-0.370704536507718) q[4];
u3(-2.38715288529355,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.86844416776031,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.84334265987989,-1.57660627135762,0.209418033999650) q[4];
u3(1.05437168214733,1.22389154745219,2.29398137277623) q[2];
u3(0.948372205361209,3.48245122933607,-1.65258481622009) q[7];
u3(1.43274547958507,2.18590586469531,-1.48602417138813) q[9];
cx q[9],q[7];
u1(3.59944795101427) q[7];
u3(-4.11119375801745,0.0,0.0) q[9];
cx q[7],q[9];
u3(-0.824287696271887,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.87498893530857,-0.534023252694126,1.10845527831436) q[7];
u3(2.26632819033164,-0.752210742949480,4.33353418002640) q[9];
u3(2.51744508201505,-0.285779100304393,1.98303722745429) q[5];
u3(1.67274904524957,-2.51003441852340,-2.49631431467424) q[1];
cx q[1],q[5];
u1(-0.300921465268045) q[5];
u3(0.792208256586122,0.0,0.0) q[1];
cx q[5],q[1];
u3(4.12674642097801,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.973467852146225,-1.03514544696533,-1.33129049484754) q[5];
u3(2.29251878871595,0.735093777423170,-0.601657912393723) q[1];
u3(1.34480342844322,1.87412832816860,-3.62756920971465) q[5];
u3(1.35302254447994,3.32200860962871,-2.65462977805983) q[3];
cx q[3],q[5];
u1(1.33869391933760) q[5];
u3(-3.04134955321546,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.577149572772624,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.932606814131380,1.31397496580584,-1.43876205852601) q[5];
u3(2.24481954284329,0.210104917008120,5.64387237809472) q[3];
u3(1.20708439208954,-0.640757569937249,2.14821389447717) q[8];
u3(0.915355112015082,-2.32150067231886,-1.01500399482982) q[11];
cx q[11],q[8];
u1(0.214233531832915) q[8];
u3(-1.57484799317067,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.0175163297114000,0.0,0.0) q[11];
cx q[11],q[8];
u3(0.473321175923281,1.55660777207491,-0.0673819678970606) q[8];
u3(1.44957678810011,-0.528554320830210,-5.23790941582727) q[11];
u3(0.221336630587558,2.72151015436284,-2.58594342312486) q[9];
u3(0.434702814073789,1.24095675080145,-2.39239008387495) q[10];
cx q[10],q[9];
u1(1.61330815049657) q[9];
u3(-2.85990600310844,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.829721670731142,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.85214780344006,-2.44502444569805,-0.321400837523588) q[9];
u3(1.46797290903552,-2.21149374028321,1.72752732120474) q[10];
u3(1.95724829189288,0.0740453703012092,-2.95540629897744) q[6];
u3(2.07852915863206,-0.323461571751359,-4.63821814176830) q[2];
cx q[2],q[6];
u1(0.263737396480046) q[6];
u3(-1.18015292241447,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.26490169173346,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.80664711691227,1.39576660439434,-2.87276391789849) q[6];
u3(0.816928348579610,-0.195613849457287,-5.55209194988589) q[2];
u3(1.76658882523320,0.908785730922794,-2.48513092817563) q[1];
u3(1.76199606962669,-4.11641024931780,2.04262897577050) q[0];
cx q[0],q[1];
u1(1.87983826788491) q[1];
u3(-2.94534067718197,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.32210603787118,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.82291970557091,-2.39806391619639,-1.31257579261697) q[1];
u3(0.494160858145623,-2.54412751948416,-1.22917201468910) q[0];
u3(0.673647466814607,2.44370803425819,-2.65567001536439) q[4];
u3(0.775113323662645,0.595017053705714,-1.81990071630341) q[7];
cx q[7],q[4];
u1(0.0884725994768296) q[4];
u3(-1.03029225716650,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.83046357166447,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.51416305010021,-1.57116543016372,4.06152097155606) q[4];
u3(1.13003906343743,5.18381240580169,0.392102451332435) q[7];
u3(1.35417653745544,2.30668039174096,0.686778735800243) q[11];
u3(1.76318843671114,0.481586997975154,-2.59258670744278) q[7];
cx q[7],q[11];
u1(2.56552255272029) q[11];
u3(-1.88990499123056,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.58142786919309,0.0,0.0) q[7];
cx q[7],q[11];
u3(0.814075189167284,-1.32368468200694,2.34998679687496) q[11];
u3(2.80963558518353,2.76438813192845,1.78752268286810) q[7];
u3(1.69307885938055,2.15915677615276,-1.85500962970516) q[6];
u3(0.829279826744732,2.08632874585444,-2.94873812456485) q[5];
cx q[5],q[6];
u1(1.57732758717021) q[6];
u3(-2.42113670675994,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.66886948879295,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.81071332999191,-0.249504724573463,4.41046278096593) q[6];
u3(1.15620123020781,3.65760101768748,0.0900589347424532) q[5];
u3(2.22219772121837,1.01031973558042,-0.599008298764276) q[2];
u3(1.80094255251520,1.46127329502764,-4.56033960090905) q[8];
cx q[8],q[2];
u1(3.51429094951590) q[2];
u3(-1.30948469558702,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.43038180335499,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.62711930156082,2.82995811708849,-1.80720960059645) q[2];
u3(1.63401733432562,4.81265185307403,-0.869865447526357) q[8];
u3(0.598482357954904,-3.18083293258659,2.38189753146862) q[10];
u3(1.82761512272453,-2.14591601713744,2.41365089718831) q[0];
cx q[0],q[10];
u1(1.56968280128614) q[10];
u3(-3.37210820435034,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.52557636242539,0.0,0.0) q[0];
cx q[0],q[10];
u3(2.62540767630304,1.35707859228751,-0.942037617679674) q[10];
u3(1.83425300447881,1.86437885158820,1.75729875799110) q[0];
u3(0.634435772972602,-1.23011546144839,0.904070619313997) q[3];
u3(1.27850329063854,-2.89522706047480,0.448396352603673) q[1];
cx q[1],q[3];
u1(3.38295839153699) q[3];
u3(-1.02528544643149,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.51113956468885,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.00074695424499,-3.04690904409595,2.50322617334205) q[3];
u3(1.29612630526241,2.32137467467684,-0.835925071564992) q[1];
u3(2.78152775067517,-2.73405185564054,1.07411661458975) q[9];
u3(2.38083098186649,-1.92047479318775,-0.455516445930328) q[4];
cx q[4],q[9];
u1(1.01551323616163) q[9];
u3(-1.26071373212385,0.0,0.0) q[4];
cx q[9],q[4];
u3(-0.191902086652816,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.38654067192323,2.29898589244513,-2.32020455730425) q[9];
u3(1.35896245366329,2.50221428477025,1.40416795712802) q[4];
u3(1.16506986676248,3.27293631229058,-1.79633939502726) q[9];
u3(1.89757757004785,0.652077495384337,-0.194835797438056) q[10];
cx q[10],q[9];
u1(1.66680536533854) q[9];
u3(-2.56205693307070,0.0,0.0) q[10];
cx q[9],q[10];
u3(3.49200107596492,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.58980729146059,-3.43027576139867,1.03347372732911) q[9];
u3(1.09789776910683,-3.35528124768712,-1.98623175420475) q[10];
u3(2.22886022767397,-3.24790480303050,0.462778772594428) q[8];
u3(2.77898037957336,-3.80697163274038,-2.15711550269869) q[2];
cx q[2],q[8];
u1(0.633266157683338) q[8];
u3(-3.41845644961199,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.88538452488118,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.780459693455800,1.63511705938748,-4.05327314333610) q[8];
u3(1.92466631448828,4.35803887161218,-1.46434738742296) q[2];
u3(1.81287622111744,2.77141150034222,-0.914583584240188) q[5];
u3(2.01640821919894,1.47149137607874,-0.871838425542487) q[1];
cx q[1],q[5];
u1(1.99314458544124) q[5];
u3(0.391882235031815,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.47198128833239,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.31906410858940,-2.34711979663475,3.12131588755754) q[5];
u3(1.12114793674278,-3.79602423213238,-2.01007865651490) q[1];
u3(2.02684802643197,1.35513667653410,0.767170419423802) q[11];
u3(1.79844127207300,-0.106961181860177,-3.79391485042089) q[4];
cx q[4],q[11];
u1(2.35617338759679) q[11];
u3(-1.80108483712378,0.0,0.0) q[4];
cx q[11],q[4];
u3(1.54244684632975,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.21921838119613,3.88526799216593,0.240613973753710) q[11];
u3(2.10234139634262,-1.42231367200463,1.00010492876081) q[4];
u3(0.627404218706089,0.880779480932384,1.81945253337650) q[6];
u3(1.45070084661602,-1.45222744082661,-0.749838408534309) q[3];
cx q[3],q[6];
u1(1.41793579851135) q[6];
u3(-3.59733709352035,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.18336250433051,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.891662277733888,1.83835017612037,-3.63033957311457) q[6];
u3(2.14505462942187,-5.63406323834400,-0.213536139133116) q[3];
u3(2.25673657756794,0.406594751620418,-2.09553065704472) q[7];
u3(1.20917755694233,-3.98212402028050,1.77205336472698) q[0];
cx q[0],q[7];
u1(0.661916941990159) q[7];
u3(-3.55615745107907,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.36652184018941,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.10736496227123,2.71923028242839,-2.26781811215487) q[7];
u3(2.88503591505312,-3.69820108871586,1.89851373895380) q[0];
u3(1.86137370333805,1.00090907536840,-3.05633264792109) q[5];
u3(2.42359051117230,2.14007154164102,-3.27334198926864) q[8];
cx q[8],q[5];
u1(0.606214523708135) q[5];
u3(-1.61798205704934,0.0,0.0) q[8];
cx q[5],q[8];
u3(3.06461511024179,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.39384767836241,-2.40710199671580,-0.105606329386661) q[5];
u3(1.75938937801683,2.03474420973838,-2.00797531722098) q[8];
u3(1.86585777226669,3.30832542819195,-1.81909230906186) q[3];
u3(1.97702684879553,1.72253823670368,-2.84527398815008) q[10];
cx q[10],q[3];
u1(0.914324132254947) q[3];
u3(-1.07566359991033,0.0,0.0) q[10];
cx q[3],q[10];
u3(-0.317239803419679,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.23036986630566,-0.392530368070056,-0.123270696157183) q[3];
u3(2.09527904942099,0.245354343986917,-1.89765849080377) q[10];
u3(1.29068966783607,3.75800094352289,-1.50951213894729) q[0];
u3(2.02608760414925,2.05503273480325,-0.0788209836583488) q[6];
cx q[6],q[0];
u1(-0.572870538972937) q[0];
u3(1.32801698895216,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.59930299442793,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.812743904357731,-0.0257910731922077,-2.11803524475529) q[0];
u3(1.51300733047937,-2.16851843493939,-2.58395601283769) q[6];
u3(2.49301448187234,-0.489275075879338,3.05875943144991) q[9];
u3(2.76409660245228,-0.329695396914974,1.36995055190336) q[4];
cx q[4],q[9];
u1(2.01421702872663) q[9];
u3(-3.12547418964047,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.535889717912592,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.737709708827694,2.04989259505425,-1.60191241836847) q[9];
u3(1.22968320127443,0.582601028852952,-3.12194045326717) q[4];
u3(0.776949768716528,-0.529985452468345,1.00138892677548) q[11];
u3(0.606759454169975,-1.44470973340115,-1.39641907307385) q[1];
cx q[1],q[11];
u1(2.51301086473551) q[11];
u3(-0.112289404232094,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.50881174038722,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.13389829682069,-0.352263958608768,1.18695165797654) q[11];
u3(1.12421787880745,4.11379213301954,-1.65253327463519) q[1];
u3(1.52288129948375,0.0125425146268436,-1.72889181748407) q[2];
u3(1.97520182670112,0.409858181851863,-4.87530957451482) q[7];
cx q[7],q[2];
u1(2.15582545581764) q[2];
u3(-3.14411075127002,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.691967304962070,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.83691341289779,2.49522614797450,-2.59973125948009) q[2];
u3(1.60853992874431,1.82566690666274,-0.0958138141478901) q[7];
u3(1.84068341852158,-1.10073265330064,0.521641427716817) q[2];
u3(1.75651318928306,-2.94725697399165,-0.718093176590804) q[3];
cx q[3],q[2];
u1(3.46057676921187) q[2];
u3(-4.36441163966057,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.0169518273478668,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.75274915326133,1.25987549667003,-4.88745882844781) q[2];
u3(2.93446614345142,0.974192232114935,-3.44244702982013) q[3];
u3(2.23837628445008,0.263779657263399,1.99466331253051) q[9];
u3(1.55527986241600,-0.552863882848022,-1.62920491865616) q[0];
cx q[0],q[9];
u1(2.78388867118865) q[9];
u3(-1.73799884659154,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.31213741679938,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.66731326667276,-1.57404231894432,-2.68817049239535) q[9];
u3(0.698925044905140,-3.28001874054078,-0.830781402999360) q[0];
u3(1.58138251679051,0.0500246223027281,1.29878137941669) q[11];
u3(1.24391371628357,-1.87722542861688,-2.29854548551979) q[7];
cx q[7],q[11];
u1(4.09961747612095) q[11];
u3(-3.89853892871326,0.0,0.0) q[7];
cx q[11],q[7];
u3(-0.262348369120533,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.65067102383624,-1.89626064860472,1.76733831231774) q[11];
u3(1.65356530019962,0.159428454341778,-1.12805706103745) q[7];
u3(1.80678677280425,1.48521811202123,-0.715017422874494) q[4];
u3(0.229235332461422,1.32964720063194,-3.67602258326742) q[10];
cx q[10],q[4];
u1(1.66305717695238) q[4];
u3(-0.431016851365319,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.12998848025115,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.44037028757087,0.422775676842747,1.26838753557469) q[4];
u3(0.250380658489530,1.56907931258743,2.65857346433262) q[10];
u3(0.0778871818123771,-0.408480478403712,0.995539345614717) q[8];
u3(0.606594159873373,-2.39529192345785,0.536132868145138) q[5];
cx q[5],q[8];
u1(-0.314822292329167) q[8];
u3(-1.72283287840458,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.881221492144524,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.49183355430405,-1.40370067966169,4.77605212617130) q[8];
u3(2.28470291391945,0.663951584815523,4.20965220476630) q[5];
u3(1.44883787121980,-1.59837384747185,0.138012506701273) q[6];
u3(2.50032201568653,-4.03930345687577,0.811840886751166) q[1];
cx q[1],q[6];
u1(0.480330353331715) q[6];
u3(-1.42801137625651,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.0144428391021576,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.84462281727382,-1.72069066517156,-0.437191685287136) q[6];
u3(2.71554216003595,1.10314117608995,4.97356908414600) q[1];
u3(1.40570016619072,0.930739924686630,-0.501330353504848) q[8];
u3(1.05314436559356,1.29267672168035,-4.51804186010279) q[11];
cx q[11],q[8];
u1(1.91343641543949) q[8];
u3(-2.57632904976473,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.09584028423332,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.76539162652928,3.91150534903224,-2.24618163772391) q[8];
u3(1.84482445859404,-2.43154347293740,2.80158875994045) q[11];
u3(2.13843558198031,0.335702638890753,1.45086848253111) q[0];
u3(1.83318750535571,-1.34299767831820,-1.95356996161863) q[3];
cx q[3],q[0];
u1(2.31046654624410) q[0];
u3(-1.89513998505761,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.58999985813665,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.06012048329353,0.653416794530376,0.926344185507817) q[0];
u3(1.33345790101581,3.98933004855070,-2.06623760078919) q[3];
u3(0.790822899709925,-2.28926582034090,1.59251877361821) q[2];
u3(0.738472619365618,2.43790330209776,-3.20671192251103) q[4];
cx q[4],q[2];
u1(2.29395226098450) q[2];
u3(0.248971706043256,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.61378467563188,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.16794355271254,0.612174135477791,2.66511897031397) q[2];
u3(0.991393549548993,-0.565465054424736,-0.656765483511344) q[4];
u3(2.21452992769961,1.81569937023517,-3.48284910760565) q[10];
u3(0.823437985087033,3.32044384349942,-2.15324746590363) q[5];
cx q[5],q[10];
u1(1.24717134259671) q[10];
u3(-0.468572989196943,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.45881696862430,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.873963091535023,-1.93892618867926,-0.666803613920323) q[10];
u3(1.74772474538336,-3.65104494029863,2.41879621568423) q[5];
u3(1.31263752429735,2.16027849503035,-3.53538747780304) q[1];
u3(0.741475255627252,2.99588790045568,-1.78043854135908) q[7];
cx q[7],q[1];
u1(1.02837662671207) q[1];
u3(-0.384649440184639,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.48305845225938,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.453130047660127,-1.37611179678833,4.23396359367549) q[1];
u3(1.77721417426756,-3.99308315709564,-0.553356326474886) q[7];
u3(0.746345907203647,1.93947506434309,-1.60500659366343) q[9];
u3(0.515283053280781,0.731328598088924,-3.09130142354343) q[6];
cx q[6],q[9];
u1(0.582527614378773) q[9];
u3(-0.329483639531839,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.18543406681631,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.19307455432585,2.37510267669604,-0.642596550254028) q[9];
u3(1.70621149954696,0.753362433282359,3.31332661211622) q[6];
u3(0.617723286690257,2.88561576392524,-2.03725677250864) q[10];
u3(1.10112726528785,-2.79337748276576,1.33214250540556) q[3];
cx q[3],q[10];
u1(0.168593949164525) q[10];
u3(-0.867020374403188,0.0,0.0) q[3];
cx q[10],q[3];
u3(3.07232706900891,0.0,0.0) q[3];
cx q[3],q[10];
u3(0.834013767973695,-1.42027314842693,1.87136050159173) q[10];
u3(1.23219640551178,0.497346752023994,-5.57150100257030) q[3];
u3(1.86282637151220,0.536157526178067,1.65141977776777) q[0];
u3(2.57634724727754,-2.41020695550556,-1.50629600858353) q[6];
cx q[6],q[0];
u1(3.04956675698320) q[0];
u3(-2.12358014945555,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.880074333597216,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.78480789632878,-0.784345144999526,-3.04478090167872) q[0];
u3(0.629708057794156,2.25857285416135,2.57526361909566) q[6];
u3(0.341540473405835,-1.67300285205016,1.94452400461476) q[5];
u3(1.24718834501138,-2.66911296203833,1.49437571156776) q[4];
cx q[4],q[5];
u1(1.52387103194519) q[5];
u3(-2.18646213867497,0.0,0.0) q[4];
cx q[5],q[4];
u3(3.76748615486645,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.70096072975084,-1.04901192467710,-2.74565382355579) q[5];
u3(1.04206571019268,-1.87940995844571,3.60120147656464) q[4];
u3(1.64972565757247,2.28024114457745,0.363687138695098) q[1];
u3(2.20033752576449,0.340762014616071,-2.42910453795874) q[7];
cx q[7],q[1];
u1(1.57669493391304) q[1];
u3(-2.87285078452264,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.890649300839539,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.41178459834078,-1.68877082353678,2.22273820918416) q[1];
u3(1.56395691351998,0.154023363940147,-2.66703446092287) q[7];
u3(1.64886138065801,2.72676388644225,-0.777186936101382) q[2];
u3(1.18586231689133,1.13992033209841,-0.945716397764724) q[9];
cx q[9],q[2];
u1(1.40514977031225) q[2];
u3(-3.31609791560671,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.50949289779343,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.97532748119700,2.10718626367665,-0.127363476590936) q[2];
u3(1.53969404985908,0.747102650812963,-4.03657732746366) q[9];
u3(1.37929020539212,1.97568582524977,-1.90904728013859) q[8];
u3(0.156610075905513,-1.76101112935882,0.956348426597314) q[11];
cx q[11],q[8];
u1(2.72724873101666) q[8];
u3(-1.62116361992562,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.157446231901104,0.0,0.0) q[11];
cx q[11],q[8];
u3(0.963055043939053,0.224605899384492,0.922865403562497) q[8];
u3(0.304898048106192,2.30056507703003,0.973465227848626) q[11];
u3(0.579476692309348,-2.38387221439978,1.48953074703294) q[1];
u3(0.475660488087521,-3.46706770188905,2.19106180461796) q[6];
cx q[6],q[1];
u1(2.50604129015423) q[1];
u3(-1.77907247914435,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.350039535317892,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.30261132360429,-1.90453663369669,3.29637101856889) q[1];
u3(0.379140020297716,0.346179786143343,0.437890183828268) q[6];
u3(2.09414581669877,0.0317069675903325,2.56463085083459) q[2];
u3(1.99069968145078,-1.47242356500730,-2.08482821620123) q[8];
cx q[8],q[2];
u1(0.857787567320442) q[2];
u3(-3.00946797939161,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.83973740168107,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.12240352009127,-1.37877188795129,4.29033749995268) q[2];
u3(1.43756676208523,-0.992755299775779,-4.71815776933496) q[8];
u3(1.36287311808647,0.761573298915978,-2.37575574779374) q[11];
u3(1.05815677877166,2.29881992358239,-3.32053230211258) q[4];
cx q[4],q[11];
u1(-0.512569881699701) q[11];
u3(1.21696273909645,0.0,0.0) q[4];
cx q[11],q[4];
u3(3.56608035828294,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.47968164840754,-1.63532931696502,2.26386226459055) q[11];
u3(1.03465730219970,-2.67999860165684,0.0459990912066017) q[4];
u3(2.38396540005332,1.60349192174991,-1.00480740613338) q[5];
u3(2.18236828956271,5.07350446605119,0.703902793717199) q[7];
cx q[7],q[5];
u1(1.22692685787089) q[5];
u3(-0.317960789803003,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.77426013962888,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.925854925834323,0.465453492205171,0.902803674735269) q[5];
u3(2.56030025373391,-0.463815630872332,1.02878340705461) q[7];
u3(2.09762529360459,0.602603355900080,2.52984335652308) q[9];
u3(1.70911258416223,-0.407323670054569,-1.58164210440614) q[3];
cx q[3],q[9];
u1(3.09991822222060) q[9];
u3(-1.60274311375902,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.862465752868200,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.86201814081190,0.252307792332672,-1.94175698344764) q[9];
u3(2.06560530604198,-1.68692169976399,-2.01784449521450) q[3];
u3(0.901780132167146,0.852006073487805,0.0955499635201025) q[0];
u3(1.04522509161033,0.382075990970372,-2.01432222681549) q[10];
cx q[10],q[0];
u1(0.689156096086039) q[0];
u3(-3.35365838747661,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.84856827894037,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.30300362520256,1.31129894791348,-2.30041959565151) q[0];
u3(1.41660883360414,1.48330358848804,-4.07309062891190) q[10];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
