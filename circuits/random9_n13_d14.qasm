OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.22984384298753,0.503466409181399,-3.64211741240840) q[12];
u3(1.02730614123189,-1.29948887501473,4.65161409664516) q[3];
cx q[3],q[12];
u1(1.52994873899904) q[12];
u3(-0.349765169432933,0.0,0.0) q[3];
cx q[12],q[3];
u3(3.24526943114800,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.43590288644678,0.955892615716757,2.47612266887859) q[12];
u3(1.41358952927731,-2.29199956543991,-2.70367495245021) q[3];
u3(2.01970778773992,2.24090969063459,-3.32331078640634) q[5];
u3(1.13544705740411,3.45446331980298,-2.51221543425859) q[8];
cx q[8],q[5];
u1(0.409530674499184) q[5];
u3(-1.57554136277437,0.0,0.0) q[8];
cx q[5],q[8];
u3(2.26405920426460,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.14016323424258,-3.91205002238501,1.62217183128384) q[5];
u3(1.67653819753823,-0.690630661351826,2.91180748210817) q[8];
u3(1.30120606896954,1.48355292006474,0.948863016262829) q[10];
u3(2.16683573703231,-1.49584809481546,-0.987611018670962) q[0];
cx q[0],q[10];
u1(0.547174300019613) q[10];
u3(-1.61524813963017,0.0,0.0) q[0];
cx q[10],q[0];
u3(-0.118101457060728,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.62454094607946,0.0922918821455287,-3.86367444769896) q[10];
u3(2.68550046164245,1.04493288719575,-0.605552944316695) q[0];
u3(1.91692251336841,-0.346157757455629,1.86013874148556) q[11];
u3(2.30230717756123,-1.99796214746401,-0.640346517945745) q[6];
cx q[6],q[11];
u1(2.06631024632606) q[11];
u3(0.0306551344164712,0.0,0.0) q[6];
cx q[11],q[6];
u3(1.25014599471060,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.34845918356866,-0.662268717200681,-1.48633677162143) q[11];
u3(1.41279540152563,-0.407564050122649,-2.95071372393628) q[6];
u3(1.40097719958072,0.00402870832218583,1.42716131671515) q[9];
u3(0.823211312397000,-1.44473595881205,-0.534490676918044) q[2];
cx q[2],q[9];
u1(3.97610759798864) q[9];
u3(-3.39022080685440,0.0,0.0) q[2];
cx q[9],q[2];
u3(-0.395310358126835,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.26549556105764,3.72616845273415,-0.660036218307185) q[9];
u3(1.83703082463513,2.82712307265503,-2.25732107445277) q[2];
u3(1.07969576223744,-1.08080918607047,-0.626654766831150) q[7];
u3(2.39769062748428,-4.03113598512075,0.324523363826007) q[4];
cx q[4],q[7];
u1(1.00232777573888) q[7];
u3(-0.155106609134305,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.35528077712943,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.67404626902718,0.148785133674088,-3.80987872642125) q[7];
u3(1.71649764348536,-1.34963014946843,3.82721173468260) q[4];
u3(0.932114876655021,3.71404229493346,-2.23665565584766) q[5];
u3(1.20017532556435,0.118724072735850,0.132945612180276) q[0];
cx q[0],q[5];
u1(3.55194900406709) q[5];
u3(-4.25077047618100,0.0,0.0) q[0];
cx q[5],q[0];
u3(-0.611247098656788,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.795352186863180,-2.04508196742000,-0.248781329115050) q[5];
u3(1.29676498698658,-0.862076437797792,2.04813349094822) q[0];
u3(2.20559660958612,0.0406523127098338,1.43565970660985) q[4];
u3(1.21469060114197,-2.36314013174003,-2.40222803343834) q[10];
cx q[10],q[4];
u1(-0.200026282859140) q[4];
u3(-1.64287716768525,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.628264261266298,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.04471571084512,-0.946409018603948,-3.02820175450843) q[4];
u3(1.23460635381084,0.461217700934675,-2.24418151727266) q[10];
u3(1.09975223968773,0.974063661425209,-3.56574279703888) q[11];
u3(0.281011012824927,2.30754408647048,-1.83427119021279) q[3];
cx q[3],q[11];
u1(0.373300844612225) q[11];
u3(-1.04281898205448,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.91265900992557,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.39582972919552,-1.53273730705489,0.596139432567971) q[11];
u3(2.39566344223137,3.72287308096859,0.342988590434562) q[3];
u3(1.07548353437069,3.37349126882316,-2.38829894825100) q[2];
u3(0.945739108778319,0.357440040253462,-0.648381716762178) q[1];
cx q[1],q[2];
u1(2.93961880598670) q[2];
u3(-1.75683467387052,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.10598514341509,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.64666630994621,-2.46285299069478,1.97949278999037) q[2];
u3(0.965082403944161,-1.43509781883975,1.21049483608279) q[1];
u3(1.68833559781870,0.724225582485872,1.19119283903231) q[6];
u3(2.29245098108778,-0.972469318418159,-0.839099846032557) q[9];
cx q[9],q[6];
u1(1.85116363393277) q[6];
u3(-0.0840949931229074,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.927391227387856,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.598450305335350,-1.74471214283309,3.32189621909028) q[6];
u3(0.947527953422458,-3.44990099826239,0.791948257121187) q[9];
u3(1.11124665531192,1.63048588303460,-3.50434367730593) q[7];
u3(1.49684204023560,-3.77782045566482,2.50189968640396) q[8];
cx q[8],q[7];
u1(0.534918067863452) q[7];
u3(-1.25818516934417,0.0,0.0) q[8];
cx q[7],q[8];
u3(-0.198343795561686,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.741879508536981,-1.30690923866169,-1.28286931951384) q[7];
u3(2.55038294994402,-2.11305234809981,-3.03042293104834) q[8];
u3(2.23183037907600,1.88347997310528,-1.29651451454403) q[12];
u3(1.98439246180721,0.0575019659987706,-3.08963465994048) q[0];
cx q[0],q[12];
u1(-0.783935275338971) q[12];
u3(-2.02689708297672,0.0,0.0) q[0];
cx q[12],q[0];
u3(1.57353993364476,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.40070438664457,-2.34566411605880,1.81974376777554) q[12];
u3(1.27891801709898,1.15607389907313,-1.99725072501092) q[0];
u3(2.09562463540080,-1.91433590629996,0.570505664766944) q[5];
u3(2.03907587146666,-1.69803617160885,-1.37233359475939) q[8];
cx q[8],q[5];
u1(2.17496306676417) q[5];
u3(-2.58848401829349,0.0,0.0) q[8];
cx q[5],q[8];
u3(0.879216799431063,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.41915042387763,-2.70207581093251,3.04763074554776) q[5];
u3(0.525888781196795,-2.08552752922395,0.880962668193159) q[8];
u3(1.67254766399920,0.770398851731489,-3.54602078950569) q[7];
u3(2.01791863595391,3.96534829270733,-2.06518793228347) q[11];
cx q[11],q[7];
u1(-0.203377328803290) q[7];
u3(-1.33199982917969,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.18449303326800,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.94372222024150,2.04044253838977,-2.47396985821861) q[7];
u3(2.82112265361280,0.106399738952762,0.0791831389948092) q[11];
u3(1.35953861109870,-3.05086270494109,2.36277005787595) q[4];
u3(1.78298857459336,-3.28826066240670,2.62991133545402) q[10];
cx q[10],q[4];
u1(1.56723175708613) q[4];
u3(-0.125640405085801,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.78315773675749,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.861254974032030,1.42881781808305,-3.56774578287537) q[4];
u3(1.46100148142551,-0.230026949404263,-5.90574745896756) q[10];
u3(0.595934453838368,1.01134441140174,-1.32351187583572) q[2];
u3(0.891313392368402,-0.0986515796940524,-1.95709737606040) q[9];
cx q[9],q[2];
u1(1.49623903397449) q[2];
u3(-0.717997504177758,0.0,0.0) q[9];
cx q[2],q[9];
u3(3.18206591269632,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.06181699924660,-1.33063850872280,2.19158446376407) q[2];
u3(0.911492018531679,-0.713953040625584,-4.55397983169023) q[9];
u3(2.53322752472150,0.543999565188672,-2.46502555182100) q[1];
u3(1.93625824579962,-3.76026089888669,2.26775358727680) q[3];
cx q[3],q[1];
u1(0.624965105480316) q[1];
u3(-1.47911305007424,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.128004499180850,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.53719942896299,1.59976429912767,0.161438761547687) q[1];
u3(2.72387304050503,2.51394235014718,-2.43023778365823) q[3];
u3(0.585765000468123,-0.977904764521235,0.948439729600743) q[2];
u3(0.582164673925165,-1.70056658194627,1.58075291925902) q[1];
cx q[1],q[2];
u1(-0.0173949245781864) q[2];
u3(-1.19232015501714,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.42249136554783,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.74394738553836,3.62017135272709,-1.58316017896790) q[2];
u3(1.53178647727841,1.27584451545943,0.0505228766341923) q[1];
u3(2.06465533219431,-0.663805181208612,1.33682816239800) q[9];
u3(1.97108616783182,-1.41726420767059,-0.969610526796377) q[10];
cx q[10],q[9];
u1(0.936895524042613) q[9];
u3(-0.0765599422460013,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.55652631901548,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.18590983950174,0.605489303603653,1.87318000739908) q[9];
u3(0.653137755902425,-2.82458184968252,0.216183889108485) q[10];
u3(2.11118421474872,0.101101395716301,0.983534218012000) q[0];
u3(1.96665135329016,-1.88492975249344,-1.53724061575740) q[4];
cx q[4],q[0];
u1(1.56566777020615) q[0];
u3(-1.92643375007654,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.375623154581970,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.31110847907735,-3.58178294009824,1.38513391731236) q[0];
u3(2.13748212453879,1.86799423439091,0.393561095106172) q[4];
u3(2.19523787170656,3.82651160711284,-1.48666380022921) q[3];
u3(1.19640947571768,2.16850808034011,-0.328772971580646) q[8];
cx q[8],q[3];
u1(0.328346542784324) q[3];
u3(-1.08061754586999,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.24434437262304,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.662862174864347,-2.61384818588104,3.55885093843940) q[3];
u3(1.25268079562226,-1.13436022478959,1.22995843665340) q[8];
u3(0.882799515930270,0.969206475124831,-0.424663947867149) q[12];
u3(2.21072498398451,-0.626417694742479,-3.05149838045078) q[11];
cx q[11],q[12];
u1(3.03434538990073) q[12];
u3(-2.59772853269657,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.859216002030004,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.31555181746878,3.72276078585547,-1.33126182162411) q[12];
u3(2.18647929760852,-1.00081125273652,-1.04707379570968) q[11];
u3(1.48843709330270,-0.431296400722654,1.25376763432356) q[6];
u3(1.09448588642291,-1.85367328886597,-2.20407363274078) q[7];
cx q[7],q[6];
u1(2.19970436653294) q[6];
u3(0.0622015797418785,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.784645984567287,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.09611436589394,-4.17237107395093,0.111772282824240) q[6];
u3(2.19551976071608,3.73289696418738,1.84523522717561) q[7];
u3(2.04589018096504,-2.77206460707011,0.0881378597938329) q[12];
u3(1.81468596254773,-3.30838237990065,-0.113878899017335) q[9];
cx q[9],q[12];
u1(2.19241651272236) q[12];
u3(-0.333681601837078,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.24006887817901,0.0,0.0) q[9];
cx q[9],q[12];
u3(0.786748469199277,-1.97747251758857,-0.793046740053348) q[12];
u3(1.73627759591576,-1.38033504804450,-1.22275402783594) q[9];
u3(0.493935428971151,2.25363023736487,-0.585684003675304) q[1];
u3(1.64704232581014,0.824964312745474,-1.21713223633731) q[0];
cx q[0],q[1];
u1(-0.435852580556374) q[1];
u3(-1.83780559453515,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.912279411631190,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.47266030170597,-0.173781624753014,0.720764090455515) q[1];
u3(2.58184463385160,2.86997954309220,1.54664233653277) q[0];
u3(2.50045925790492,3.40460500486259,-2.42088869174716) q[7];
u3(0.858259645494393,-0.818216020846407,3.18949912674052) q[10];
cx q[10],q[7];
u1(3.81422636105119) q[7];
u3(-4.19574780089583,0.0,0.0) q[10];
cx q[7],q[10];
u3(-0.963674730079223,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.376187169416019,1.36033878716827,-2.54004979121720) q[7];
u3(2.21290454782576,1.33705974038180,-0.413302592816810) q[10];
u3(0.662154071219396,3.73630931515021,-2.17910184591933) q[8];
u3(1.89166559278972,1.27253010305297,-2.33673828189681) q[2];
cx q[2],q[8];
u1(2.88798761067133) q[8];
u3(-1.75425313055671,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.780137225931596,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.830831567851100,1.79213207994809,1.14132351574683) q[8];
u3(2.74806290061920,2.21998955049325,-2.31182181647582) q[2];
u3(2.65428357209770,0.661344797428661,1.29288605375993) q[11];
u3(1.94527118728522,-2.31531613915817,-2.44148183593275) q[4];
cx q[4],q[11];
u1(0.981277108393012) q[11];
u3(-1.07543393238571,0.0,0.0) q[4];
cx q[11],q[4];
u3(-0.345059986291567,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.70541478335415,1.94554007492582,-4.11861552512702) q[11];
u3(2.07472362564151,1.78600633370057,-1.28756002355862) q[4];
u3(1.87665334507208,-0.0771183883145761,1.42930532947453) q[5];
u3(1.66706652221710,-2.40443920404521,-2.11513166935184) q[6];
cx q[6],q[5];
u1(1.79094638089487) q[5];
u3(0.224397391079447,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.24482141810624,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.76857539843623,-0.912349985490849,3.69920651572244) q[5];
u3(1.81564979388499,1.52654440066080,-1.98584842109849) q[6];
u3(2.08318575589579,1.82997878324462,-4.44697688451354) q[6];
u3(0.582186701738124,-1.55028602359444,3.18955170893723) q[8];
cx q[8],q[6];
u1(-0.174547360399308) q[6];
u3(0.963611842806089,0.0,0.0) q[8];
cx q[6],q[8];
u3(3.71461724676598,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.915604542392177,-2.64053378393285,0.245868818133599) q[6];
u3(2.40155729072681,4.63307152477292,-0.867509511066747) q[8];
u3(1.09797728715554,1.21991425142793,-0.196680305207892) q[5];
u3(1.46080584495225,-0.00236075231607091,-2.18862433381124) q[10];
cx q[10],q[5];
u1(0.526337402655690) q[5];
u3(-1.14050204040207,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.93013938699477,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.17049747215857,2.07845270457551,-2.04708863427301) q[5];
u3(2.13301168163853,-2.34581966178021,-1.37178200230651) q[10];
u3(2.52413213959606,0.798639658280732,-0.901688034138525) q[2];
u3(2.20863081826908,0.902700069112985,-4.19914752652719) q[12];
cx q[12],q[2];
u1(2.02437939546531) q[2];
u3(-2.66149853959460,0.0,0.0) q[12];
cx q[2],q[12];
u3(-0.0666147607995144,0.0,0.0) q[12];
cx q[12],q[2];
u3(2.63407239414057,-2.66340009932666,3.50276975511838) q[2];
u3(1.51304007173838,-1.56547769599572,0.693236052905303) q[12];
u3(2.02580490066855,2.26424658846141,-1.16206752775636) q[9];
u3(1.45961419329294,0.176539971530283,-3.19488367035832) q[11];
cx q[11],q[9];
u1(1.83496908431258) q[9];
u3(-2.76444351801814,0.0,0.0) q[11];
cx q[9],q[11];
u3(0.596127108243264,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.52194502122338,0.613944831837099,0.187090788427681) q[9];
u3(2.34656025766051,1.32272512012184,-0.909324122629961) q[11];
u3(1.89979285799916,-3.18225462665255,0.876820646893656) q[7];
u3(1.66723241422158,-2.99724079831311,0.654134705016252) q[1];
cx q[1],q[7];
u1(1.38795135186812) q[7];
u3(-0.814048231813420,0.0,0.0) q[1];
cx q[7],q[1];
u3(-0.428605342988614,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.74434966390526,-0.415082769650331,3.22104121267155) q[7];
u3(0.933458759459365,-3.69814313502695,-2.34462016725112) q[1];
u3(1.93322144147225,-0.168611307928236,1.87987401730153) q[4];
u3(1.18007404922943,-2.55467630631374,-2.55648263165912) q[3];
cx q[3],q[4];
u1(-0.193052530199112) q[4];
u3(-1.54525071792716,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.18297764544586,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.660001591051672,-2.04996830577146,3.81155642793287) q[4];
u3(0.875712026909576,-1.89415488625727,0.508141122878869) q[3];
u3(1.67188092635432,0.288759546229563,1.85947313274861) q[3];
u3(1.28875606079144,-0.677933296264728,-1.51792209804294) q[5];
cx q[5],q[3];
u1(0.00740646987916005) q[3];
u3(-1.64997149276185,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.05114378684066,0.0,0.0) q[5];
cx q[5],q[3];
u3(0.866260584509044,-2.90666516659589,2.51742839516581) q[3];
u3(2.06553888431225,-0.250634538459214,-2.50681883602704) q[5];
u3(2.42593233655928,-1.18644834101727,-0.274522836200539) q[4];
u3(1.66396299294402,-2.08916243251058,0.956778711114324) q[9];
cx q[9],q[4];
u1(2.81015774876820) q[4];
u3(-2.18757586333095,0.0,0.0) q[9];
cx q[4],q[9];
u3(1.56733966570380,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.472722131855914,2.57320083731160,-1.35813029564236) q[4];
u3(1.50263381504528,1.90537371218788,3.30549241180238) q[9];
u3(0.972015804784292,0.934735540481936,-1.35381136759910) q[1];
u3(0.951510919484121,-0.269079470705218,-0.271437296374163) q[6];
cx q[6],q[1];
u1(0.00238035293068894) q[1];
u3(-1.53629074343398,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.920785288964799,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.27350560974711,0.0570668045495287,-1.84160499450030) q[1];
u3(1.49548377310094,1.78807266213991,-1.92509332335261) q[6];
u3(2.47506811580720,-1.97747250995692,1.70302792665406) q[11];
u3(2.00809956850464,-2.23654896670954,-0.333838698009760) q[2];
cx q[2],q[11];
u1(1.07887210772158) q[11];
u3(-0.270496559330242,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.73315693427704,0.0,0.0) q[2];
cx q[2],q[11];
u3(0.430915902206336,-4.48834048831638,1.59389027279326) q[11];
u3(2.48772868812779,0.657091202176498,5.27904562009223) q[2];
u3(2.50638591632635,-0.710031803994548,2.31489837841637) q[7];
u3(2.07629519099262,-1.81034271629782,-1.75093627794564) q[12];
cx q[12],q[7];
u1(1.18592315514511) q[7];
u3(-3.20487276537897,0.0,0.0) q[12];
cx q[7],q[12];
u3(1.70354651280554,0.0,0.0) q[12];
cx q[12],q[7];
u3(0.894596092763034,0.0433875682770331,1.44045690438934) q[7];
u3(1.62123313516394,-0.0853766496554078,-3.32408866073557) q[12];
u3(1.30927151607057,2.79054062474204,-0.0626683178242369) q[10];
u3(1.02313488234812,0.914103356594465,-1.39522782727213) q[0];
cx q[0],q[10];
u1(2.16762383703636) q[10];
u3(-2.89876649379486,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.48500526826065,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.348613244660654,-2.09991501645008,0.0536601992948700) q[10];
u3(0.522575022848659,2.02298560900089,-1.24778484886169) q[0];
u3(2.11586383639605,1.50716294013251,1.39786769752567) q[2];
u3(0.470052644472666,0.380003667956787,-5.18651050085906) q[5];
cx q[5],q[2];
u1(1.92317352077458) q[2];
u3(-2.88779880888232,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.0152179221787536,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.69573333646401,0.221753847268533,2.38363221262018) q[2];
u3(2.33032717814134,2.96649381143771,-0.406390339853473) q[5];
u3(1.31332913668751,1.35108315470610,-3.05376313756764) q[8];
u3(0.721246117876020,2.61513021471853,-3.58092378631183) q[9];
cx q[9],q[8];
u1(0.712494509149496) q[8];
u3(-3.17779664158789,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.94487643219785,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.97622222786234,3.42526186224128,-0.373391946345089) q[8];
u3(1.96051846307051,1.24069117673236,-0.319234872274036) q[9];
u3(0.524730736703851,-2.77909826169014,3.03335639207829) q[0];
u3(0.651876695171755,0.592637241444334,-0.990336934223693) q[10];
cx q[10],q[0];
u1(-0.137488259994888) q[0];
u3(-2.03305349875868,0.0,0.0) q[10];
cx q[0],q[10];
u3(1.49477793018322,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.23708989237068,4.33452239072763,-1.56855709422401) q[0];
u3(0.920896381297682,3.70128355019825,1.53117310645900) q[10];
u3(1.90996796914080,-1.38699012726606,0.407806456086727) q[11];
u3(1.76537248468389,-2.69433665761589,-0.906168273707395) q[7];
cx q[7],q[11];
u1(2.91770836111739) q[11];
u3(-2.62270242758319,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.15079027997027,0.0,0.0) q[7];
cx q[7],q[11];
u3(2.46820934876283,1.39780485592769,-2.19207501538129) q[11];
u3(2.05150799872426,-3.27933744743734,-0.254588407190568) q[7];
u3(1.77166773733748,-1.47241055642059,-0.0559560636910598) q[6];
u3(1.64838542597420,-1.96590775609726,0.929672484010177) q[1];
cx q[1],q[6];
u1(0.624281551209716) q[6];
u3(-1.50810601366810,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.93305836139368,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.38904634613197,0.175786794043274,0.484104768567740) q[6];
u3(0.860173429682216,4.38558153271977,0.0200898226887616) q[1];
u3(2.63896643594805,-3.30867758813713,1.63753369335334) q[3];
u3(1.56902540919215,-0.475899556419815,2.94625468418128) q[12];
cx q[12],q[3];
u1(1.63309244141311) q[3];
u3(-0.188661524793115,0.0,0.0) q[12];
cx q[3],q[12];
u3(2.26328934978402,0.0,0.0) q[12];
cx q[12],q[3];
u3(1.04269431282087,-3.05564873564502,0.276548425542579) q[3];
u3(1.64911322205262,2.58413455986382,-2.53723178005375) q[12];
u3(2.06305707610799,2.68219287506314,-0.406494119340532) q[7];
u3(2.59117450097223,1.89057889946573,-1.60930624186512) q[10];
cx q[10],q[7];
u1(2.06379159007223) q[7];
u3(-2.91248193100350,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.424630809281998,0.0,0.0) q[10];
cx q[10],q[7];
u3(0.711878374737788,0.989964338078064,-2.62361055148655) q[7];
u3(1.10809489381871,4.43739691164729,-1.19990583177489) q[10];
u3(1.33274174365737,1.76153479706208,-3.47604938714786) q[12];
u3(2.01815911283220,2.22201018321074,-3.01862680713727) q[0];
cx q[0],q[12];
u1(3.27529860654335) q[12];
u3(-1.22673897176561,0.0,0.0) q[0];
cx q[12],q[0];
u3(2.39516159490717,0.0,0.0) q[0];
cx q[0],q[12];
u3(1.51956301530821,0.267928727094031,3.20399414196747) q[12];
u3(1.70845424275385,0.491927356246094,-5.63537970075059) q[0];
u3(1.65548144049066,1.12693237112944,-1.23733032863603) q[6];
u3(1.69652073991745,-0.358126643845811,-3.40775522024727) q[8];
cx q[8],q[6];
u1(1.25260784464279) q[6];
u3(-0.262774258980768,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.00849370617752,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.786393904697419,1.80819693139707,-1.73549999984673) q[6];
u3(0.395942417105670,-1.10297461605693,2.21941736709150) q[8];
u3(2.13980096062483,-0.0659400415366502,2.52352631401874) q[11];
u3(1.22669536698114,3.32258851289399,2.78224978735344) q[4];
cx q[4],q[11];
u1(1.58507564551045) q[11];
u3(-3.27290857055218,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.72885119015876,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.46462638485436,-2.71325294002883,0.716295723042497) q[11];
u3(1.21811990081556,-0.760801076701570,5.21921634900154) q[4];
u3(1.82307902451417,-1.30695168173626,3.91519403296290) q[5];
u3(1.27865879604520,1.54914602355128,1.77500714817984) q[3];
cx q[3],q[5];
u1(1.83989535507722) q[5];
u3(-3.16879547249358,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.907262053266952,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.828042200046039,-0.132305466504869,1.99202049296556) q[5];
u3(1.78001345002323,-0.0187966782183879,-3.67723417271501) q[3];
u3(1.89159227149422,1.05037121061571,-4.09464529315615) q[2];
u3(2.12005245826527,-1.81548622755743,4.11944453381407) q[1];
cx q[1],q[2];
u1(3.28063664113713) q[2];
u3(-4.29160293812809,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.515701617303971,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.68023919585259,-4.39857717827909,1.31482928792383) q[2];
u3(0.487975449763676,3.93876849959137,-0.802406975970501) q[1];
u3(2.44731073205360,3.52860316337000,-1.39175485852127) q[10];
u3(2.03564918288036,1.92324892203840,-0.118310702476183) q[0];
cx q[0],q[10];
u1(0.824793226551289) q[10];
u3(-0.316408253564936,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.66284334284573,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.85047391580553,0.720579865928220,1.65464068834554) q[10];
u3(1.03539228137371,-2.43356931260711,-1.69804676028467) q[0];
u3(2.51808773455172,3.42916361825747,-0.882381550338008) q[9];
u3(1.00803224977254,0.751030512200953,-1.88094290383337) q[8];
cx q[8],q[9];
u1(-0.577016980507618) q[9];
u3(-1.99060386075751,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.06575775405540,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.66616091196845,-1.09983039682690,4.55226757643650) q[9];
u3(0.735694373992051,-4.75162244129178,1.38726774821715) q[8];
u3(2.11433995877864,-0.643273956318450,0.395952119527299) q[7];
u3(1.70420332838206,-2.34764119009529,-1.00812599671099) q[5];
cx q[5],q[7];
u1(3.14450043105539) q[7];
u3(-0.924065535363090,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.60711618845042,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.17132091686272,-0.172397678197358,-0.890710135309970) q[7];
u3(1.84607180268622,-3.16234463288226,-2.07746652269657) q[5];
u3(1.43823782198401,-0.950862460469996,1.75717664483616) q[12];
u3(2.35602611695571,-1.51492375360314,-2.37244147752821) q[3];
cx q[3],q[12];
u1(1.93558615273049) q[12];
u3(-3.27696090834541,0.0,0.0) q[3];
cx q[12],q[3];
u3(1.14380876541899,0.0,0.0) q[3];
cx q[3],q[12];
u3(2.20119172435495,-3.72395229590878,0.636709993484350) q[12];
u3(2.22687018212259,-0.0191656336124053,-2.78303745625130) q[3];
u3(0.961653102506073,-1.66728290459088,-0.225795649423925) q[6];
u3(2.05403805306755,-3.64441943314636,0.778830399915496) q[1];
cx q[1],q[6];
u1(3.53119499395840) q[6];
u3(-4.36480734160452,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.617992392077435,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.50499101958784,-4.70062226579032,1.33963164462713) q[6];
u3(1.03804149376836,-3.02826605125902,1.74643130360121) q[1];
u3(1.34209285851882,-0.434233008906686,1.04168705482554) q[11];
u3(1.81231390422876,-1.55493689542854,-2.20965947125578) q[2];
cx q[2],q[11];
u1(-0.115807304910839) q[11];
u3(-2.58658223250638,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.72246884561537,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.97156273482576,1.08677999149221,1.23405945232324) q[11];
u3(1.78087275784340,-3.64493960996215,2.31155777540479) q[2];
u3(0.854537745839257,-0.124448449412333,-2.06400336814141) q[10];
u3(1.22226702573431,0.686980026078811,-4.48615827383259) q[4];
cx q[4],q[10];
u1(1.87713927196570) q[10];
u3(-2.23142931688129,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.187837377011142,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.01777744519500,-2.24955394965532,2.08040780155003) q[10];
u3(1.10707934608332,4.40263091452618,0.248177435053560) q[4];
u3(0.522301859514675,-2.23783882702315,0.444907399041196) q[12];
u3(1.47768639366645,-3.16537942771209,1.39334098015321) q[3];
cx q[3],q[12];
u1(3.27005960025760) q[12];
u3(-4.06455856490856,0.0,0.0) q[3];
cx q[12],q[3];
u3(-0.573417072415130,0.0,0.0) q[3];
cx q[3],q[12];
u3(1.29585783662765,2.35747363823780,-1.10572579692747) q[12];
u3(2.01601254106538,-2.35449239186217,1.19727433140586) q[3];
u3(1.69895796606178,1.54050242965091,-0.741167574336726) q[9];
u3(0.787825890292630,0.318942431900263,-3.30593191355834) q[0];
cx q[0],q[9];
u1(0.135831085771903) q[9];
u3(-1.00799072964121,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.99938826379747,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.03063490643482,-2.26025240331109,1.28071762441284) q[9];
u3(1.09242611977526,-0.871233648450945,2.04267851799803) q[0];
u3(1.78024336964398,2.36276177437375,-2.77931857346958) q[11];
u3(0.789130399252486,-3.00542958386891,2.58982571530121) q[7];
cx q[7],q[11];
u1(3.01746625316602) q[11];
u3(-2.13613372544044,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.12278670848983,0.0,0.0) q[7];
cx q[7],q[11];
u3(2.87708336975424,0.365982292088967,2.36479810438200) q[11];
u3(2.47719338098335,4.26875593061594,-1.08280008902659) q[7];
u3(1.38539998583150,1.42605203743397,-3.05946414557861) q[5];
u3(1.23140866121577,2.46714784692576,-3.15291823276241) q[6];
cx q[6],q[5];
u1(1.21916533828502) q[5];
u3(-0.327654234565429,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.0628357387076153,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.24655733644551,-0.888081301792047,2.55387944761263) q[5];
u3(1.91258774189604,0.784844869059393,5.49671250919122) q[6];
u3(0.722531270042786,-1.78777736729582,2.14803616113604) q[8];
u3(0.481427850056311,-0.846909075029979,-1.61378559434860) q[2];
cx q[2],q[8];
u1(2.66658978846039) q[8];
u3(-1.48350754321924,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.237253095563190,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.61100507654540,-0.682078453785109,3.85803692489162) q[8];
u3(0.684868878204576,-0.864922909732301,-1.75307939143549) q[2];
u3(2.62493826562470,1.34546337889128,-1.57653327033208) q[6];
u3(2.54227217756376,0.272118027066015,-4.35665827202718) q[2];
cx q[2],q[6];
u1(-1.09097486190730) q[6];
u3(0.599905230010033,0.0,0.0) q[2];
cx q[6],q[2];
u3(3.44966523967638,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.59415058923228,1.41098858233065,-3.90852395365997) q[6];
u3(0.540334600164567,-0.418519190655094,-3.94502418438074) q[2];
u3(0.700435848973380,2.26051617460378,-3.36876229648886) q[3];
u3(0.764133093082538,1.33461919522284,-3.40788549542175) q[12];
cx q[12],q[3];
u1(3.20529558331674) q[3];
u3(-0.826484915464992,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.27299661332540,0.0,0.0) q[12];
cx q[12],q[3];
u3(2.20874934561914,3.53731613426667,-0.518088548692467) q[3];
u3(1.55248615403926,-1.71255089704712,1.89225545839192) q[12];
u3(0.881499937231795,1.00258324982205,-1.02071099332972) q[8];
u3(0.563881721956303,1.76015242271728,-4.38510547908120) q[4];
cx q[4],q[8];
u1(1.60397004270731) q[8];
u3(0.330116693556856,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.612787019741211,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.72722610519555,-1.50310610330808,2.32571880605765) q[8];
u3(2.90473301977832,-4.75998514259522,-0.0622945056213324) q[4];
u3(2.35908643855544,2.25554369112598,-0.839221651685597) q[11];
u3(2.45846408032196,1.36289319554940,-4.26135118297794) q[1];
cx q[1],q[11];
u1(1.41024397665882) q[11];
u3(-0.000692114015063749,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.45114131487310,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.44498502599848,1.39763627639507,-2.66181700565857) q[11];
u3(2.03311956506032,-1.29606650693089,0.354859451210625) q[1];
u3(1.43066573830197,1.58176021996912,-1.36055376573438) q[0];
u3(1.22116358434228,0.0163575025123885,-3.11796653913531) q[5];
cx q[5],q[0];
u1(3.41844773051295) q[0];
u3(-1.34764125845459,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.30352772412110,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.507151083003832,3.06627802004446,-1.05673519953127) q[0];
u3(1.96729597520037,-4.41422917271171,0.542515211482496) q[5];
u3(1.73019899077355,2.86527404199899,-2.19552550373581) q[9];
u3(0.974375634153474,2.10617962716818,-2.67445206795872) q[7];
cx q[7],q[9];
u1(-0.230013921006244) q[9];
u3(-0.852961247096376,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.80139554844010,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.846951890607824,2.09694558534364,-3.95444916223623) q[9];
u3(1.65133926145610,2.78866627849412,-1.18926137212048) q[7];
u3(1.76316836504366,0.738887361098638,0.0819536088628293) q[4];
u3(0.247663485582281,-3.08819359035296,-1.88089525868410) q[6];
cx q[6],q[4];
u1(-1.25198753010182) q[4];
u3(0.349111287171485,0.0,0.0) q[6];
cx q[4],q[6];
u3(3.59581404045289,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.834429657472281,3.06419040716802,-1.65697735483797) q[4];
u3(2.84592531455296,4.12840946577411,1.64750871377886) q[6];
u3(2.21922595948645,1.21852447592350,-3.01051194423050) q[3];
u3(1.33699854814581,-3.55175187874058,2.58068845055660) q[10];
cx q[10],q[3];
u1(0.101803693897692) q[3];
u3(-1.52819992054025,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.00164169136310,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.16341452803207,-3.30838649675195,1.48609849296046) q[3];
u3(1.86477841213676,5.06112344897067,0.727336185280488) q[10];
u3(2.31173243901157,-2.55857893170690,-0.430865804861881) q[8];
u3(1.67625628941091,-4.67682149086034,-1.36644720057076) q[2];
cx q[2],q[8];
u1(2.48169343550147) q[8];
u3(-2.90188051174454,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.06381061188982,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.99977225219929,0.797911585534126,2.89225400211357) q[8];
u3(1.05731320724742,-1.13325105177333,0.232583942990881) q[2];
u3(2.51954303391283,-0.350393392961226,-1.59022135324124) q[5];
u3(1.68817680207862,0.559536263831832,-3.75655128557454) q[9];
cx q[9],q[5];
u1(1.62226096284796) q[5];
u3(-2.47042020673453,0.0,0.0) q[9];
cx q[5],q[9];
u3(0.226123968012178,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.29682023524975,-2.16452736694770,-1.59321141616258) q[5];
u3(1.97700692352753,-3.84528529767111,2.00184871107083) q[9];
u3(2.45337749635304,2.49532434628351,0.303344001117216) q[12];
u3(2.62048902448213,4.02050139726156,-1.28012265635879) q[7];
cx q[7],q[12];
u1(2.58459628694421) q[12];
u3(-1.88582541574948,0.0,0.0) q[7];
cx q[12],q[7];
u3(0.742923862341806,0.0,0.0) q[7];
cx q[7],q[12];
u3(0.705188643098766,0.584893181154899,-1.26482752141444) q[12];
u3(1.56322067540756,3.25358486054330,0.416949079474046) q[7];
u3(1.24973975216509,1.92869074647537,-2.86788589058342) q[0];
u3(1.01871007669667,1.03490769617425,-1.90243581779583) q[11];
cx q[11],q[0];
u1(1.55836132126674) q[0];
u3(-1.03376170715666,0.0,0.0) q[11];
cx q[0],q[11];
u3(-0.532176447485099,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.81395260837598,0.0860910527813774,2.93173229960379) q[0];
u3(2.42989741373244,0.317448102539698,-2.86376128586035) q[11];
u3(1.56903539702408,1.48109765538868,-3.73308833469242) q[8];
u3(1.49609672933796,-2.32474123138882,3.71612162088192) q[5];
cx q[5],q[8];
u1(2.14609552180116) q[8];
u3(-2.99269340552462,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.67390326266189,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.48630837005612,-1.98046983578638,0.521879197154617) q[8];
u3(0.251151310436842,0.109804469309557,-3.49749749137384) q[5];
u3(1.22033276437853,2.37645898696943,-1.41667914567667) q[7];
u3(1.58827063575323,1.31502413006371,-0.123671607471635) q[1];
cx q[1],q[7];
u1(1.47954580062297) q[7];
u3(-0.313553262923408,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.19441031287077,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.272302779851774,-3.21341220272133,-0.954841613900366) q[7];
u3(1.83392400217332,-0.160194791758356,0.797375520010730) q[1];
u3(1.16778675507827,-3.54356774069445,2.19377621812382) q[4];
u3(1.36503561408179,2.88751361033821,-2.52388138502289) q[6];
cx q[6],q[4];
u1(1.47020960211351) q[4];
u3(0.0325848735029433,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.746906036828929,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.76674494926585,-1.05825474265590,3.19756726814382) q[4];
u3(1.36699363494134,-0.362267179823436,-1.97670080832520) q[6];
u3(2.38628876171187,-1.05933317308511,0.413070134722960) q[12];
u3(1.71624001013203,-3.63415720321246,1.16004934692172) q[11];
cx q[11],q[12];
u1(3.06266868100261) q[12];
u3(-2.19079156453123,0.0,0.0) q[11];
cx q[12],q[11];
u3(0.441406959190378,0.0,0.0) q[11];
cx q[11],q[12];
u3(1.62232765986861,-0.660847885504991,-1.50135694019988) q[12];
u3(2.03438782555534,-1.09595242923834,-0.486236450152121) q[11];
u3(0.478596249911836,-0.201636509504957,0.252243615626220) q[2];
u3(0.596575391111761,1.88088300067391,-2.26766936482645) q[3];
cx q[3],q[2];
u1(0.593498724082653) q[2];
u3(-1.47216886124287,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.48286753068881,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.50968312108207,-3.19935923531615,0.265909906086376) q[2];
u3(1.86785981882549,-2.71733409060176,0.394751961955157) q[3];
u3(2.29049103184708,-1.43252039479037,4.26401320721354) q[10];
u3(1.30303025094090,1.13330374822401,0.678570375704454) q[0];
cx q[0],q[10];
u1(0.270299307410849) q[10];
u3(-1.44303456225345,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.75869190064981,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.856749724328542,-1.86169527905107,0.387968141918916) q[10];
u3(0.614643576589400,-3.77214661300331,-2.35578997007587) q[0];
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