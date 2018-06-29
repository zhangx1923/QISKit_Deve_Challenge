OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(2.06208592619961,1.85688385062932,-0.286778660300714) q[5];
u3(1.46726234535306,0.334866362940915,-2.19466333509160) q[7];
cx q[7],q[5];
u1(2.20655508022759) q[5];
u3(-2.65066536955809,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.0789554058191240,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.22003826106266,0.866358488216535,-4.23045683014738) q[5];
u3(1.87976857513652,-4.08559211981495,-1.73199705136364) q[7];
u3(2.39089113056381,2.78723490528330,-0.656548337434726) q[12];
u3(2.42267306692498,1.23054205673805,-4.25228791264642) q[8];
cx q[8],q[12];
u1(3.36630071191332) q[12];
u3(-0.791957223499842,0.0,0.0) q[8];
cx q[12],q[8];
u3(1.81533478987044,0.0,0.0) q[8];
cx q[8],q[12];
u3(2.45537962211554,-1.93975448867924,0.753895553911251) q[12];
u3(2.07049584359764,-0.634651765753632,5.20817808348643) q[8];
u3(0.213224181994545,-0.421104372314268,0.913541692060361) q[4];
u3(0.140093394783462,-2.94018913110202,1.10970975617059) q[3];
cx q[3],q[4];
u1(1.46940305279112) q[4];
u3(-0.952247628137490,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.293334745721583,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.31965617503112,0.266462232694181,0.580273955266189) q[4];
u3(1.19530422157324,3.11746005426526,-0.654853591356327) q[3];
u3(2.34226847441760,-0.0947010832144654,-1.15245490149871) q[11];
u3(1.70035687078920,-4.59378067860418,0.997849239316905) q[6];
cx q[6],q[11];
u1(2.09904059736158) q[11];
u3(-2.52890008168890,0.0,0.0) q[6];
cx q[11],q[6];
u3(1.45886987543183,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.50858481388219,-1.46068345793954,-1.15401753084304) q[11];
u3(0.487080028547120,1.41389224151756,2.76676929980590) q[6];
u3(0.654747578115396,0.232941689819998,-0.174376606314697) q[13];
u3(0.569130451769560,-0.512170649630690,-1.67242821048043) q[0];
cx q[0],q[13];
u1(-0.371641858701169) q[13];
u3(0.911966872730103,0.0,0.0) q[0];
cx q[13],q[0];
u3(3.86377539945690,0.0,0.0) q[0];
cx q[0],q[13];
u3(0.803540115167950,1.80103090823421,-2.59946829980060) q[13];
u3(1.90770048706885,0.730802789979377,-5.02949135692794) q[0];
u3(1.20997869247378,-0.411646408401206,1.06070893281278) q[9];
u3(0.251454845157399,-2.04269302520692,0.756430744212881) q[14];
cx q[14],q[9];
u1(3.23990234084377) q[9];
u3(-2.16706105385462,0.0,0.0) q[14];
cx q[9],q[14];
u3(1.35967494524380,0.0,0.0) q[14];
cx q[14],q[9];
u3(1.09216937551531,2.56692662451047,0.937921393983097) q[9];
u3(1.20551895243833,3.11353702476713,2.04238098374495) q[14];
u3(1.07510226379717,-1.57193801215453,0.0351327932277513) q[10];
u3(0.755100415617976,-2.12242600135904,0.137324224823366) q[1];
cx q[1],q[10];
u1(1.76787833045854) q[10];
u3(-3.12460409715418,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.06211244714317,0.0,0.0) q[1];
cx q[1],q[10];
u3(2.55916772139384,0.564127460492658,-3.26743490797060) q[10];
u3(1.30161015385761,-2.17253922217038,-1.54572361365152) q[1];
u3(1.64409840063710,3.48832555001271,-1.75921460625063) q[2];
u3(0.666441040295966,0.935167003680486,-0.00685484642905987) q[15];
cx q[15],q[2];
u1(-0.351461335636891) q[2];
u3(0.512154522678232,0.0,0.0) q[15];
cx q[2],q[15];
u3(4.42542086947050,0.0,0.0) q[15];
cx q[15],q[2];
u3(1.00004208398221,0.930184040203861,-0.636275087520999) q[2];
u3(1.58964241962120,-3.33860827173856,0.921902762259754) q[15];
u3(1.24635313116690,-0.416508688373389,-1.90574548699186) q[2];
u3(1.82888241610784,2.73752119276209,-3.45303071363966) q[5];
cx q[5],q[2];
u1(1.82991560319715) q[2];
u3(-2.39260175671520,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.483944302052076,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.561570490116197,-2.27278495277233,0.465391069490966) q[2];
u3(1.47503872143734,-0.396640844850299,5.16669038377053) q[5];
u3(1.44252200254507,2.13193524708467,-0.201374548266233) q[6];
u3(1.97643087438429,0.364841325362875,-3.91354743975887) q[11];
cx q[11],q[6];
u1(2.04976981923371) q[6];
u3(-2.27866288135841,0.0,0.0) q[11];
cx q[6],q[11];
u3(3.22479898492880,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.992679161529916,-2.84196620468757,2.48096956669311) q[6];
u3(0.469734413496303,-0.404738562348480,-0.671598066269034) q[11];
u3(1.01863239149742,-2.83733353083497,0.971162423912446) q[9];
u3(1.76056214440215,-3.77272410433189,0.438486197750675) q[3];
cx q[3],q[9];
u1(0.806530637618498) q[9];
u3(-3.42640776564277,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.96544723925035,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.34122921318379,0.136428635375004,-1.56945073535797) q[9];
u3(1.14085814441515,-5.29091370590455,0.0353560031873781) q[3];
u3(2.10003023836082,-0.739239816210370,-1.34098322724280) q[15];
u3(2.38456286179052,0.876252862156913,-5.19835777247012) q[14];
cx q[14],q[15];
u1(1.65330086076673) q[15];
u3(-2.23971024936027,0.0,0.0) q[14];
cx q[15],q[14];
u3(3.80943569101440,0.0,0.0) q[14];
cx q[14],q[15];
u3(1.17769477119569,1.64791623399398,-2.15571440089045) q[15];
u3(1.41507190268224,4.04166240128066,-1.28541602615813) q[14];
u3(1.26560008598844,-0.169013896740178,1.37458266611767) q[4];
u3(1.22906577012387,-2.67557370705457,-0.744533100890056) q[8];
cx q[8],q[4];
u1(-0.124099221905019) q[4];
u3(-1.06047197747472,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.44560114840135,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.54976066149836,2.56357527473788,1.26842968847642) q[4];
u3(2.55507130026141,2.11384304048615,0.345404758769974) q[8];
u3(2.07332111012794,0.435400755449616,-0.451466235694887) q[1];
u3(0.450127323054290,0.159189685371104,-4.75384511184808) q[10];
cx q[10],q[1];
u1(1.80791392036519) q[1];
u3(0.247122677887207,0.0,0.0) q[10];
cx q[1],q[10];
u3(0.776478466398277,0.0,0.0) q[10];
cx q[10],q[1];
u3(2.07461159553549,0.395757727550240,2.48685417167085) q[1];
u3(1.16537182169998,-1.42215199097788,4.70740930555798) q[10];
u3(0.728086751892459,-0.785500365143834,0.669658595097637) q[13];
u3(0.353590377069722,-1.43246942423891,-0.445246443169422) q[7];
cx q[7],q[13];
u1(1.83641744531296) q[13];
u3(0.462261794049890,0.0,0.0) q[7];
cx q[13],q[7];
u3(0.943608550273042,0.0,0.0) q[7];
cx q[7],q[13];
u3(0.335708285469959,-0.110993552890565,-3.22573427242651) q[13];
u3(0.130455580002758,2.59529109683767,1.16286943288218) q[7];
u3(1.08139953726962,0.600366932300984,-0.630455667173452) q[0];
u3(0.863800784340037,-3.44501843420341,0.752639077548384) q[12];
cx q[12],q[0];
u1(1.57643783647300) q[0];
u3(-2.17721251483545,0.0,0.0) q[12];
cx q[0],q[12];
u3(-0.0693815455246658,0.0,0.0) q[12];
cx q[12],q[0];
u3(0.958565468434040,-2.02398425972275,0.898908760208195) q[0];
u3(1.49482554957706,3.56293610365039,-2.31723699017280) q[12];
u3(2.36926943228760,-2.33707966918361,0.645615589188815) q[0];
u3(1.91300473323127,-2.67455294521696,0.0854598966215951) q[12];
cx q[12],q[0];
u1(1.31722788180969) q[0];
u3(-0.703962015946018,0.0,0.0) q[12];
cx q[0],q[12];
u3(0.0721646071604032,0.0,0.0) q[12];
cx q[12],q[0];
u3(0.357707790681252,-1.89898830455465,2.86530672054683) q[0];
u3(2.71438986264014,4.19740435119379,1.76266800249319) q[12];
u3(2.08416625228303,-0.543023490179450,-0.217427831260548) q[2];
u3(0.238525980396047,-0.528743753256952,-4.91841391474506) q[15];
cx q[15],q[2];
u1(2.61231900674392) q[2];
u3(-1.74161469502470,0.0,0.0) q[15];
cx q[2],q[15];
u3(0.772453135289323,0.0,0.0) q[15];
cx q[15],q[2];
u3(2.25876336227664,-1.91555365915850,0.996689462648218) q[2];
u3(2.03561156533391,2.90891670268561,-0.424460614431463) q[15];
u3(1.93902347835352,-1.29829275614943,0.740708245035685) q[5];
u3(1.99117446753127,-3.68205978612444,0.541748678122037) q[10];
cx q[10],q[5];
u1(1.20370758691104) q[5];
u3(-3.49946438436090,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.57862415188693,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.970225550186662,0.117253497440751,2.73345789238436) q[5];
u3(0.480997279311047,-2.08990646645717,3.11308070943871) q[10];
u3(1.45240756583216,1.04069318573470,-1.70178530590539) q[8];
u3(2.49289267762734,0.827765273338357,-5.22091526320237) q[4];
cx q[4],q[8];
u1(2.33026146296392) q[8];
u3(-1.91096248540839,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.146168561710368,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.90767806148475,0.382341984348072,-4.40779998970424) q[8];
u3(1.76074882018969,2.18598533602260,1.45698157471386) q[4];
u3(1.69325925496109,-3.35216646540506,1.51797856760833) q[9];
u3(2.06291205414610,0.0950103944548157,2.84135827272453) q[7];
cx q[7],q[9];
u1(1.65184943959476) q[9];
u3(-0.623926488631380,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.66853191559590,0.0,0.0) q[7];
cx q[7],q[9];
u3(0.817649226771064,-2.61344082709191,2.72416434869276) q[9];
u3(1.49569698943558,-1.66284759975075,1.19842079009949) q[7];
u3(1.73176330362206,0.825358001303952,0.443781663755088) q[3];
u3(0.920594844688935,0.877516989951066,-5.08947225498125) q[13];
cx q[13],q[3];
u1(1.34830062273635) q[3];
u3(-0.353397163591945,0.0,0.0) q[13];
cx q[3],q[13];
u3(1.93689127107850,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.81755850911385,-0.805848626186498,1.19977442291605) q[3];
u3(1.49038903003203,2.47702220925111,2.26618972000499) q[13];
u3(1.15534557005380,-0.870151632933300,0.515149785330576) q[1];
u3(1.58438627422533,-1.17294163396304,-1.68240134504254) q[14];
cx q[14],q[1];
u1(1.63594393785426) q[1];
u3(0.0888128372578543,0.0,0.0) q[14];
cx q[1],q[14];
u3(2.72824077214260,0.0,0.0) q[14];
cx q[14],q[1];
u3(0.483748568540171,0.805468531423256,-1.15883022380718) q[1];
u3(2.53851749891833,-1.83477207322902,1.18992222687406) q[14];
u3(1.53610243426471,0.390101002438257,-3.35112890095658) q[11];
u3(1.92159272479547,3.27836364243083,-2.47921903547245) q[6];
cx q[6],q[11];
u1(1.51793697947123) q[11];
u3(-0.729447501030331,0.0,0.0) q[6];
cx q[11],q[6];
u3(3.43172198024537,0.0,0.0) q[6];
cx q[6],q[11];
u3(1.57271821050391,-2.43527047806999,0.360588425922582) q[11];
u3(2.03465085325785,-0.890866499975859,1.19123695146040) q[6];
u3(2.07175415111291,-1.29392286993314,-0.998939987819433) q[1];
u3(1.83047621198436,-2.88641392423762,-0.0725128942965880) q[9];
cx q[9],q[1];
u1(2.87951558012946) q[1];
u3(-1.86789352982473,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.851493683445245,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.16650435975822,1.34095780422685,-2.31811925380714) q[1];
u3(2.02138500846544,5.02088696622826,0.920362621470795) q[9];
u3(0.528024685049308,2.26337604156036,-0.465784271841669) q[11];
u3(1.53818152128110,1.87543437766615,-1.06440626546054) q[10];
cx q[10],q[11];
u1(1.00576858810834) q[11];
u3(-1.34765190005051,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.32837424478951,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.25818532449806,-1.28572515915206,2.46231956674595) q[11];
u3(2.00849351630362,-2.03205831313093,-4.23047598897145) q[10];
u3(1.44014490822956,-0.187668662007326,2.00168827353478) q[2];
u3(1.44802761059390,-2.62925632929323,-1.36017973205515) q[15];
cx q[15],q[2];
u1(0.798226556295760) q[2];
u3(-1.56784356787909,0.0,0.0) q[15];
cx q[2],q[15];
u3(2.73762612520746,0.0,0.0) q[15];
cx q[15],q[2];
u3(1.26029322606242,2.08169734871057,-3.31597578791385) q[2];
u3(0.817081926409668,1.69666195040785,1.11932170221553) q[15];
u3(1.02623534188572,1.42056162018447,-1.96431259397020) q[0];
u3(0.163987449344703,1.35836064722165,-3.05105626478116) q[5];
cx q[5],q[0];
u1(2.29931131551736) q[0];
u3(-1.72443175905610,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.303562363420070,0.0,0.0) q[5];
cx q[5],q[0];
u3(0.692277316314997,-0.508148001322869,3.65516712636922) q[0];
u3(0.848533080920093,1.43098958838189,-4.15587538072138) q[5];
u3(1.52331661084276,-0.422038010908154,0.990221703345656) q[7];
u3(1.90755634492948,-2.46032953013776,-2.01748094270229) q[12];
cx q[12],q[7];
u1(4.35315504205112) q[7];
u3(-3.69145716624802,0.0,0.0) q[12];
cx q[7],q[12];
u3(-0.860449117528969,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.02364183055286,-0.268453329199009,0.449976075310088) q[7];
u3(1.76173163368828,-0.882857213165601,-2.75731265967112) q[12];
u3(2.71739417980087,-0.157212748660708,1.08656844721631) q[3];
u3(1.04210323549878,-2.92829036908314,-1.77495462291705) q[13];
cx q[13],q[3];
u1(2.71374982399962) q[3];
u3(-1.75248105289459,0.0,0.0) q[13];
cx q[3],q[13];
u3(0.833579382292438,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.07081817356476,2.13864092905996,-4.01603931540344) q[3];
u3(1.62066591900529,-1.14115163726187,3.97874400942304) q[13];
u3(1.53802758098869,2.66224865385903,-1.62131741440592) q[8];
u3(1.69164081362474,1.54066412545552,-2.63278917519644) q[4];
cx q[4],q[8];
u1(3.60928295904435) q[8];
u3(-1.15385575111830,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.76526200955712,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.99879630105192,3.31838690060933,-2.94050223738712) q[8];
u3(2.83753674971190,-0.715722839782718,2.59655188878159) q[4];
u3(2.21060791233329,2.43340990909307,-2.64442096899356) q[6];
u3(1.43001347391119,0.777822037460854,-1.91191742681228) q[14];
cx q[14],q[6];
u1(-0.345267747648160) q[6];
u3(-2.18209361019663,0.0,0.0) q[14];
cx q[6],q[14];
u3(1.28885464339273,0.0,0.0) q[14];
cx q[14],q[6];
u3(1.51258774114533,-0.189282711367144,-3.50025250016415) q[6];
u3(1.95450300994853,1.28270971020813,2.97152943763825) q[14];
u3(2.61976337630496,2.55445122120162,0.0635812435567074) q[12];
u3(3.00369881464501,1.15729484372349,-3.21383331608562) q[4];
cx q[4],q[12];
u1(0.766709974696437) q[12];
u3(-1.21656647489164,0.0,0.0) q[4];
cx q[12],q[4];
u3(3.11599415861751,0.0,0.0) q[4];
cx q[4],q[12];
u3(1.12394449280112,0.467500720915347,-2.96837172872239) q[12];
u3(0.802888040717762,-4.12309957871823,-0.660987776796494) q[4];
u3(0.943285004922291,-1.60147210306989,0.192986608535745) q[15];
u3(0.982471500326930,-4.19892253680328,-0.00390829243966317) q[7];
cx q[7],q[15];
u1(2.19270116955382) q[15];
u3(-1.80863956761480,0.0,0.0) q[7];
cx q[15],q[7];
u3(-0.0127211591561203,0.0,0.0) q[7];
cx q[7],q[15];
u3(1.64205281788958,2.49596736173462,-1.04864320227179) q[15];
u3(1.79728289145992,0.623332026024363,4.95246802359754) q[7];
u3(0.201227295238455,-2.14491427333846,1.49799130809282) q[6];
u3(0.897964273589964,-2.99395679271278,1.27920239336740) q[0];
cx q[0],q[6];
u1(3.49956918172288) q[6];
u3(-1.34494047562127,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.14716727185315,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.29069201933716,-1.69825099015599,3.29443024380744) q[6];
u3(2.46440255323486,-3.02860954287969,2.74770794432126) q[0];
u3(0.880475091465765,-2.26979838252324,1.07957225007114) q[8];
u3(0.474989599911282,-1.62186770372314,-0.499473739892164) q[13];
cx q[13],q[8];
u1(0.949876377663325) q[8];
u3(-1.46726392215393,0.0,0.0) q[13];
cx q[8],q[13];
u3(-0.159103646222038,0.0,0.0) q[13];
cx q[13],q[8];
u3(1.67147426628762,1.05774299024262,-3.13464093330983) q[8];
u3(1.31465454838060,-0.880018117578527,-2.55759604189031) q[13];
u3(1.85336485850225,-0.317946186447133,-2.56304136970250) q[11];
u3(2.34599200136656,0.242564487929245,-5.02712926125502) q[3];
cx q[3],q[11];
u1(0.0917993344604595) q[11];
u3(-1.26511004441605,0.0,0.0) q[3];
cx q[11],q[3];
u3(2.43008055988708,0.0,0.0) q[3];
cx q[3],q[11];
u3(2.82043861613865,-1.75512221923279,3.82037501647018) q[11];
u3(0.527468733570005,3.44336896080719,0.979123612461040) q[3];
u3(2.02423538565114,-0.461587684357139,0.294332105956529) q[2];
u3(2.22043646125208,-0.929785553113093,-1.70826098752261) q[1];
cx q[1],q[2];
u1(1.88547306998840) q[2];
u3(-0.0215110528453979,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.653286598142290,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.50758775730821,1.77633934439061,-1.89222397271589) q[2];
u3(1.39207027215121,-1.53626881200892,1.27925772274688) q[1];
u3(2.73082466526610,1.06312517084886,0.555273091812954) q[10];
u3(0.943043181064994,-4.92589428329392,0.599781038512561) q[14];
cx q[14],q[10];
u1(0.808213032354103) q[10];
u3(-0.589574704233441,0.0,0.0) q[14];
cx q[10],q[14];
u3(3.22766243731629,0.0,0.0) q[14];
cx q[14],q[10];
u3(1.66721044715836,0.113304929486091,-1.11876250500394) q[10];
u3(1.90052026303186,2.76249056736458,2.41960094539610) q[14];
u3(1.84480871943456,3.69577740555342,-1.49055602957697) q[5];
u3(0.343852457091741,1.72285903225644,-1.47158518578170) q[9];
cx q[9],q[5];
u1(0.944433723898646) q[5];
u3(-3.71156535993427,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.54311414079818,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.37109100940834,-1.28395495692229,0.764231593766641) q[5];
u3(2.36468641100210,2.72064351997921,-0.768618589897818) q[9];
u3(2.02754725242289,0.815989209703433,-1.33712837537213) q[10];
u3(0.730520760626407,-4.44699219073671,1.30786510187655) q[8];
cx q[8],q[10];
u1(1.95437864582598) q[10];
u3(-3.18234427168140,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.487246665939045,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.40352925136015,-1.68079685941508,0.674423385948881) q[10];
u3(0.856143372039266,1.28382574581696,-1.01657718198391) q[8];
u3(1.50804183517683,0.711603212887450,-3.47211578080419) q[2];
u3(2.79960321159783,3.66866338976229,-2.27685709567200) q[12];
cx q[12],q[2];
u1(3.13874127471147) q[2];
u3(-1.11536708526975,0.0,0.0) q[12];
cx q[2],q[12];
u3(2.27576556702710,0.0,0.0) q[12];
cx q[12],q[2];
u3(0.154148119556318,-0.296153577521222,0.878792529046136) q[2];
u3(2.70594346208545,-4.26083910870865,-0.991485845190578) q[12];
u3(1.85356979801872,0.329375498042973,-2.20597560630317) q[7];
u3(1.98153110550113,-3.14935126990233,2.39537865870987) q[1];
cx q[1],q[7];
u1(2.93819911540802) q[7];
u3(-2.00660668480511,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.39091781047172,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.24210665871827,1.29899675302828,-2.02726792450004) q[7];
u3(1.88078021080278,0.908830031152916,3.14811869325544) q[1];
u3(2.08066373352033,0.177159095969295,1.71850258628118) q[5];
u3(2.10063178732207,-1.26694315552456,-0.170196329534310) q[3];
cx q[3],q[5];
u1(-0.476150053903403) q[5];
u3(-1.60254260052812,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.75536853250882,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.785473639122989,-0.334749444431353,-1.49895821774557) q[5];
u3(0.900464774755783,-5.01018860392948,0.502905886560639) q[3];
u3(1.47528618038854,-1.16604820903409,2.15047505352693) q[15];
u3(0.767840135031576,-1.32372645350441,-2.09335129456781) q[11];
cx q[11],q[15];
u1(1.61097163111646) q[15];
u3(-0.797328891793500,0.0,0.0) q[11];
cx q[15],q[11];
u3(-0.494315870414859,0.0,0.0) q[11];
cx q[11],q[15];
u3(2.49352568453047,-0.187355475766578,-2.98401317725748) q[15];
u3(2.04737427157823,-4.03036151940992,1.97236022542081) q[11];
u3(2.30195442059768,4.21149002840544,-1.28332625402074) q[4];
u3(1.12729943512563,1.38580353998159,0.0869289162480202) q[0];
cx q[0],q[4];
u1(-0.122775319740196) q[4];
u3(-1.65196121672047,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.947530787873167,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.41635445880114,2.93131747243972,-1.07858995058714) q[4];
u3(2.01482777193540,1.31517931022114,-0.144213874347003) q[0];
u3(2.22131323988353,1.69204345395243,-2.37985757582936) q[14];
u3(2.10911110912533,1.77525670438685,-3.39738863856459) q[9];
cx q[9],q[14];
u1(2.23484274045871) q[14];
u3(-1.60448889059272,0.0,0.0) q[9];
cx q[14],q[9];
u3(-0.00209715553319012,0.0,0.0) q[9];
cx q[9],q[14];
u3(1.90234482142122,-1.87589406120689,1.16663029892970) q[14];
u3(1.84969191003937,5.51856643212802,-0.709623302874605) q[9];
u3(1.41361282080998,3.07958078298629,-2.34920736610249) q[6];
u3(2.58271037956481,0.957152702293177,-2.06879224051788) q[13];
cx q[13],q[6];
u1(3.44601197557823) q[6];
u3(-1.56814806989116,0.0,0.0) q[13];
cx q[6],q[13];
u3(1.92804848220724,0.0,0.0) q[13];
cx q[13],q[6];
u3(1.65831761876456,-1.04428274888054,-2.43463648510856) q[6];
u3(1.65106846400929,-2.61613499577090,-1.44019848928494) q[13];
u3(2.32394908353142,0.619789966489130,-1.63820073062624) q[12];
u3(1.94077164984790,4.39276593561913,-0.00377308425817180) q[10];
cx q[10],q[12];
u1(-0.347482183090788) q[12];
u3(1.27832971680573,0.0,0.0) q[10];
cx q[12],q[10];
u3(3.74250765375163,0.0,0.0) q[10];
cx q[10],q[12];
u3(1.14757542270972,4.43199597854356,-0.636257053096553) q[12];
u3(1.21005060944505,4.77307080341907,1.46617299024714) q[10];
u3(0.778210610147418,0.0318438438009285,-2.68561974838951) q[2];
u3(2.20100401925149,-3.68432390037721,2.09181832176487) q[3];
cx q[3],q[2];
u1(-0.452231274861679) q[2];
u3(-1.76922180626624,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.758594905308250,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.21041686404694,-2.26825283709143,-1.56162190351101) q[2];
u3(0.751168028208974,-1.69115243971860,2.99004090937742) q[3];
u3(2.19491442019599,-2.82538776251788,-0.193218873751137) q[11];
u3(2.43319870997522,0.344653521582806,1.22898014890204) q[4];
cx q[4],q[11];
u1(2.10119522570897) q[11];
u3(-1.60144126788801,0.0,0.0) q[4];
cx q[11],q[4];
u3(4.13127140725107,0.0,0.0) q[4];
cx q[4],q[11];
u3(1.88944388178078,-2.97559866413079,2.14132147897171) q[11];
u3(0.859772954021523,3.76822768640711,2.23064918937789) q[4];
u3(1.24441232535841,-2.47583463211516,-0.0515565259787460) q[6];
u3(1.38548851224473,-3.44283907161537,0.799207832633693) q[8];
cx q[8],q[6];
u1(2.21492411838357) q[6];
u3(-3.08774391768504,0.0,0.0) q[8];
cx q[6],q[8];
u3(1.09132638690478,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.30372652310525,1.71344152147118,1.10077137195524) q[6];
u3(2.42526402573326,1.90970622575615,2.70045487557256) q[8];
u3(2.51277566797670,-1.93924707473991,4.06878447678671) q[1];
u3(1.72067440030205,0.721583722436619,0.792512820045777) q[5];
cx q[5],q[1];
u1(0.503921606854984) q[1];
u3(-1.25897738766975,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.68442813548414,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.33507614728121,-1.24665931466154,1.88036181242713) q[1];
u3(2.02959104594855,-2.10988225388617,2.34035804525413) q[5];
u3(1.48046710328194,-0.0762788249522459,-2.25097349023929) q[0];
u3(1.19553908909289,0.371233098825865,-3.88461650005272) q[9];
cx q[9],q[0];
u1(1.06194423458329) q[0];
u3(-3.08837358495638,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.81542981015493,0.0,0.0) q[9];
cx q[9],q[0];
u3(0.503443913956683,0.435625941475682,1.08589928853395) q[0];
u3(2.09538613964306,2.81753545722710,-2.18266448434486) q[9];
u3(1.62894898693274,1.79007455857639,-3.55824522086679) q[15];
u3(1.37021375298351,-2.42616917047622,3.52207507448219) q[14];
cx q[14],q[15];
u1(3.23375415102961) q[15];
u3(-1.34028386469318,0.0,0.0) q[14];
cx q[15],q[14];
u3(0.984571220880107,0.0,0.0) q[14];
cx q[14],q[15];
u3(1.89426652203771,3.05370920747480,-0.295420065373222) q[15];
u3(0.565621310631560,3.89273689643717,-1.08418981604656) q[14];
u3(1.30675108495889,-0.578129986387937,-0.245053904003632) q[13];
u3(1.27276119979971,-3.03682499360863,-0.156244084670778) q[7];
cx q[7],q[13];
u1(-0.712429572091979) q[13];
u3(1.15898980743861,0.0,0.0) q[7];
cx q[13],q[7];
u3(3.60909364224692,0.0,0.0) q[7];
cx q[7],q[13];
u3(1.16886132045541,0.386931545082378,2.55192182333635) q[13];
u3(1.86676682254786,4.11678361455168,0.774015446047362) q[7];
u3(2.20860076775447,-0.0927574811094914,1.17087727853813) q[13];
u3(2.22756954579167,-0.771658358177802,-2.14151699337837) q[4];
cx q[4],q[13];
u1(1.59643791131855) q[13];
u3(-2.35630474559883,0.0,0.0) q[4];
cx q[13],q[4];
u3(3.23547611090121,0.0,0.0) q[4];
cx q[4],q[13];
u3(1.19891360410526,-1.13887214969132,-1.95703735803869) q[13];
u3(0.100366320891139,1.63263883598429,-1.76651234116448) q[4];
u3(1.63896624240768,2.70272496189294,-0.620522404179142) q[2];
u3(1.04189624926898,0.859763792788755,-1.61850219135750) q[9];
cx q[9],q[2];
u1(1.00958352707663) q[2];
u3(-1.33220413187765,0.0,0.0) q[9];
cx q[2],q[9];
u3(3.01606337689808,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.96349486222748,-2.49379257169390,-0.705861428637589) q[2];
u3(0.743957550441662,-1.75510933133008,0.617836062220129) q[9];
u3(0.665939582452303,0.872263807715617,-0.0216227275699200) q[7];
u3(0.786569743518701,-1.67316101612704,-0.305558950726916) q[8];
cx q[8],q[7];
u1(1.32825274233235) q[7];
u3(0.0449564175827086,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.99993400836289,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.51375502982647,4.03055018023947,-1.37002390199378) q[7];
u3(2.87933220167120,3.41201055187476,-2.46714542079472) q[8];
u3(1.41996396724072,0.810804562291891,-2.11602257498089) q[5];
u3(0.937492864376354,2.34130359808273,-3.33274772903503) q[10];
cx q[10],q[5];
u1(2.97773674187875) q[5];
u3(-1.25630412302778,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.265881372532450,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.24550288032007,0.474951469477216,-3.87771783751317) q[5];
u3(2.11705351159148,3.90745241752810,-0.462291484063316) q[10];
u3(2.15021506417777,0.967842790620239,-3.59749869145198) q[3];
u3(1.61776736760511,2.50436683576526,-3.24530582522569) q[11];
cx q[11],q[3];
u1(1.34659411593510) q[3];
u3(-3.54966237254919,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.14672925402741,0.0,0.0) q[11];
cx q[11],q[3];
u3(1.23004826193950,-3.42367567995360,1.57610320766338) q[3];
u3(1.30341774538331,1.48939691299256,1.04723284605997) q[11];
u3(2.00644293551816,2.23883260947309,-3.67753127662900) q[0];
u3(0.149479246060799,2.39572757221276,-1.22253595429518) q[1];
cx q[1],q[0];
u1(0.752233673840247) q[0];
u3(-3.32739230705421,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.38358400739952,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.09110931664783,3.75033250373281,0.121879053776980) q[0];
u3(1.51491380621380,-3.74000399981020,-2.20106312876742) q[1];
u3(0.369447140517663,-0.777977822043639,-0.213213947589187) q[15];
u3(0.342262873911891,-1.20764575075352,-0.631157298209564) q[12];
cx q[12],q[15];
u1(0.244894060768196) q[15];
u3(-1.43428555398611,0.0,0.0) q[12];
cx q[15],q[12];
u3(2.58526766124460,0.0,0.0) q[12];
cx q[12],q[15];
u3(1.21992297120717,0.185332538265015,1.23855446899759) q[15];
u3(0.484801287767157,1.52277752328342,-2.42026802670811) q[12];
u3(2.25791365533011,1.39098894822606,-2.00524696387255) q[14];
u3(1.84551840286187,2.01444907321213,-3.21154847107543) q[6];
cx q[6],q[14];
u1(-0.125103777724426) q[14];
u3(-1.58733870494141,0.0,0.0) q[6];
cx q[14],q[6];
u3(0.720409592021847,0.0,0.0) q[6];
cx q[6],q[14];
u3(2.11676086212641,-2.11808008436404,3.44537167336369) q[14];
u3(1.76320272838912,-1.12879721994923,3.40911360898655) q[6];
u3(2.17885299668785,0.563742425341942,2.21551366656171) q[4];
u3(1.40955260765503,-2.77943754904274,-2.58395892831092) q[13];
cx q[13],q[4];
u1(2.38402395365201) q[4];
u3(0.312612119322767,0.0,0.0) q[13];
cx q[4],q[13];
u3(1.52702968882404,0.0,0.0) q[13];
cx q[13],q[4];
u3(2.67652185888247,-1.34718007515598,-0.569302410405714) q[4];
u3(0.440273825695754,1.85430934125052,0.840584858445329) q[13];
u3(1.71527625148835,-2.14739804332120,-0.0593720929881341) q[1];
u3(2.25097743812514,-3.59114555934379,-1.05133003158040) q[9];
cx q[9],q[1];
u1(2.43269003987686) q[1];
u3(-2.06282005363003,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.24165390763227,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.35088236628628,1.13199420118232,-1.83206892500775) q[1];
u3(0.701559050274218,4.34478493396770,1.32732745928978) q[9];
u3(1.36502036640440,-1.64847236235654,1.96421366126776) q[5];
u3(0.465387408938518,1.57284149658966,-2.07893770050036) q[0];
cx q[0],q[5];
u1(0.358151540270652) q[5];
u3(-1.11346022848292,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.69597540687092,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.94895908066556,0.391703687170187,-1.77430875773634) q[5];
u3(2.65230660402458,-4.69215543244336,-0.896777735028874) q[0];
u3(2.68566013089507,-1.26426783623182,-1.34516974732560) q[15];
u3(1.17502263483577,-5.00692597216164,0.698850524879507) q[6];
cx q[6],q[15];
u1(-0.575814950571117) q[15];
u3(-2.09080258265705,0.0,0.0) q[6];
cx q[15],q[6];
u3(1.54208347433629,0.0,0.0) q[6];
cx q[6],q[15];
u3(0.410354130408359,0.530969304520478,2.48617551088410) q[15];
u3(0.889364974809331,1.00401509379811,1.72963580508825) q[6];
u3(2.52795192259500,-2.69235243805530,3.40795430951511) q[8];
u3(0.742742530899216,-0.908733581118755,2.01638705141625) q[2];
cx q[2],q[8];
u1(-0.150698401097256) q[8];
u3(-2.10026945099493,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.13862552087465,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.846193498489982,-2.75593048434352,2.06565404074818) q[8];
u3(1.30281516352107,2.45722362369193,0.720230743600787) q[2];
u3(2.07929463103360,-1.10601624475803,-1.10992723640047) q[3];
u3(1.30893290694925,-2.68443320403889,0.206496433690134) q[14];
cx q[14],q[3];
u1(1.94113652927892) q[3];
u3(-3.14523889793337,0.0,0.0) q[14];
cx q[3],q[14];
u3(0.881736647847894,0.0,0.0) q[14];
cx q[14],q[3];
u3(1.97009429256567,-2.46639957220021,1.35852108673329) q[3];
u3(1.47959627604188,0.786557757585650,1.69877827668230) q[14];
u3(1.10187286909270,0.234479865145657,-2.61172236268073) q[12];
u3(1.91099927315166,-3.30953155089586,2.30346600288998) q[10];
cx q[10],q[12];
u1(1.75863541662287) q[12];
u3(-3.35089058912115,0.0,0.0) q[10];
cx q[12],q[10];
u3(1.21403983896910,0.0,0.0) q[10];
cx q[10],q[12];
u3(1.52344078254238,0.352059946083386,-1.70001283834112) q[12];
u3(0.574772619527355,0.698041534891849,0.660317672361850) q[10];
u3(2.23130736781448,-2.20100602235755,-0.925871550087442) q[7];
u3(0.668394801240746,-4.26638670590362,-0.299243646523710) q[11];
cx q[11],q[7];
u1(1.53960535913567) q[7];
u3(-0.0980612727151577,0.0,0.0) q[11];
cx q[7],q[11];
u3(2.56602661717605,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.10569449522009,3.25610896967270,-2.52379828420736) q[7];
u3(1.58963578717736,3.80311196294450,0.376089075330034) q[11];
u3(1.79673303978779,0.0244576290037282,-2.33243152608077) q[11];
u3(1.80606164519703,0.848780824851961,-4.10671017458176) q[15];
cx q[15],q[11];
u1(1.31069482051037) q[11];
u3(-2.77489670916267,0.0,0.0) q[15];
cx q[11],q[15];
u3(2.21517658446282,0.0,0.0) q[15];
cx q[15],q[11];
u3(2.77796434856264,-2.92425597847433,0.230499250421408) q[11];
u3(1.82462389611566,-3.43358081017337,1.99670419622529) q[15];
u3(1.03631939158991,-0.669113624690038,-0.695182306452566) q[13];
u3(1.19250390475045,-3.24471688404818,1.04147834631473) q[9];
cx q[9],q[13];
u1(0.714904266511970) q[13];
u3(-3.07941418290605,0.0,0.0) q[9];
cx q[13],q[9];
u3(1.95830566739183,0.0,0.0) q[9];
cx q[9],q[13];
u3(0.493218145328038,-1.08557086834155,1.90923580759480) q[13];
u3(1.79442395134798,-0.915970698833958,-3.44811861053775) q[9];
u3(1.34607666651305,0.416962698672930,1.60041194061568) q[5];
u3(1.66075212045391,-2.10463122908377,-1.02082056446741) q[6];
cx q[6],q[5];
u1(2.56340303694595) q[5];
u3(-1.61111166943589,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.309359015803704,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.93203720347844,-1.75029522219658,2.94830423431059) q[5];
u3(2.14726329020534,4.36777074659470,-1.61670990821634) q[6];
u3(0.284173552245695,3.30549893158999,-2.62420946023962) q[14];
u3(1.36265169346079,0.571285332524669,-1.17213795466699) q[12];
cx q[12],q[14];
u1(1.61275363570557) q[14];
u3(0.0552109925651314,0.0,0.0) q[12];
cx q[14],q[12];
u3(0.502066806051536,0.0,0.0) q[12];
cx q[12],q[14];
u3(1.91576623515280,-0.266961435018389,-0.0535120201249802) q[14];
u3(2.48964783980030,-1.53539515744997,3.59129188583800) q[12];
u3(1.22900948564558,-1.39335136323739,-0.466063659698575) q[0];
u3(0.875289364745225,-1.95343829757496,0.526535104379832) q[2];
cx q[2],q[0];
u1(1.65387035890239) q[0];
u3(-0.360077280498109,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.09067909965225,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.08420602331062,1.11355204131981,-0.807340992895937) q[0];
u3(1.41440213186459,-5.11188398577529,0.375163217884771) q[2];
u3(1.43224950207258,-1.29005575016273,0.815729551614985) q[4];
u3(1.06227367251284,-2.02728357856853,-0.153783567550464) q[1];
cx q[1],q[4];
u1(0.900062762128977) q[4];
u3(-1.57080498508526,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.389195078338710,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.691595505938887,0.972455214391625,0.572783684154273) q[4];
u3(0.821256474309036,0.654315184542395,1.24787537788143) q[1];
u3(1.97407642002519,2.13148145306526,-3.22907457763190) q[3];
u3(0.864010651650215,2.81062841904146,-2.68483068265147) q[10];
cx q[10],q[3];
u1(-0.475614046078818) q[3];
u3(-1.99783608922563,0.0,0.0) q[10];
cx q[3],q[10];
u3(1.55826151015947,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.22710085570752,1.76068574495482,-4.31109612697750) q[3];
u3(0.470366248635001,0.549136720982264,4.37283839597025) q[10];
u3(1.08528849230763,-0.976158686901298,1.73285518965639) q[7];
u3(1.06467009716326,-1.47704120286596,-1.98952332129455) q[8];
cx q[8],q[7];
u1(-0.530628216595169) q[7];
u3(0.978113582114752,0.0,0.0) q[8];
cx q[7],q[8];
u3(3.42074010022101,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.27865451936606,-0.0132477750159383,-3.37256109355525) q[7];
u3(0.998728865462287,0.430084350766873,5.43390490085105) q[8];
u3(1.54509895972924,-0.254738877062691,-1.32247090570168) q[15];
u3(2.46261833110497,0.709170655261727,-5.09004828412729) q[3];
cx q[3],q[15];
u1(0.643736673624062) q[15];
u3(-3.22653727543603,0.0,0.0) q[3];
cx q[15],q[3];
u3(1.76991718010998,0.0,0.0) q[3];
cx q[3],q[15];
u3(1.47462128760414,-1.63844973598760,2.42986006139909) q[15];
u3(0.425390927592186,-1.31900533877559,-2.16908857225239) q[3];
u3(0.677707372048317,-2.53733101551643,1.48168426002958) q[9];
u3(0.427686940588420,1.63294606366696,-3.26614163800760) q[7];
cx q[7],q[9];
u1(0.913048339660207) q[9];
u3(-0.711695096813771,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.83551308002296,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.80740012216942,0.603345927790458,-3.25483053521165) q[9];
u3(1.83615757042645,0.0873866204406341,0.468123912418362) q[7];
u3(0.761834468969683,1.46090680435521,-3.73206283627187) q[10];
u3(1.16735297056640,2.31949296493389,-2.40787159099690) q[5];
cx q[5],q[10];
u1(2.69412513400242) q[10];
u3(-1.54021795643048,0.0,0.0) q[5];
cx q[10],q[5];
u3(3.39705442957601,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.04417249302637,-3.66856533048485,1.89753059725931) q[10];
u3(2.05929188039785,-0.417376478305918,3.81797665280021) q[5];
u3(0.853278716860885,1.70519941503666,-0.871580729508099) q[11];
u3(1.39792489203668,0.176918566956974,-2.92451134824045) q[14];
cx q[14],q[11];
u1(3.24735802460557) q[11];
u3(-3.54247654725994,0.0,0.0) q[14];
cx q[11],q[14];
u3(-0.916135587388474,0.0,0.0) q[14];
cx q[14],q[11];
u3(2.13486740068019,-2.45250448515250,1.12266869140256) q[11];
u3(0.608492094702085,-2.97047683110475,2.97822845202913) q[14];
u3(0.837348343440346,-1.37363810408418,0.460832210300136) q[2];
u3(1.09350929753289,-1.62129635196810,-0.693594616870348) q[1];
cx q[1],q[2];
u1(0.991767065572296) q[2];
u3(0.278295201679958,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.77969151455718,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.40387418189939,-0.303007787960817,-1.35813525272969) q[2];
u3(0.732186639335318,-2.90466667289025,-3.02506880786963) q[1];
u3(1.09122197638241,0.888722737017361,-0.662598665277419) q[0];
u3(0.211722991041117,-3.83357273457898,1.40586587599769) q[6];
cx q[6],q[0];
u1(3.04339893304574) q[0];
u3(-1.59378401985835,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.70664589561475,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.613635325705737,-4.22473622517930,2.05631773127184) q[0];
u3(1.13209045977843,3.67395600973122,0.375375044688607) q[6];
u3(1.76452440669024,-0.591234923384009,2.26411656955032) q[8];
u3(2.53740764132906,-1.28873137534800,-1.21627018225934) q[12];
cx q[12],q[8];
u1(0.914601009084229) q[8];
u3(0.0181648919904960,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.84676387335682,0.0,0.0) q[12];
cx q[12],q[8];
u3(1.90623992968076,0.313125700298555,-3.81500589168172) q[8];
u3(0.741198246188282,1.68242046597785,-2.59152419929966) q[12];
u3(2.10914552367125,3.11922994539198,-2.58319416694986) q[13];
u3(1.96371323941919,1.43547166503288,-2.00027125168820) q[4];
cx q[4],q[13];
u1(0.582755902451172) q[13];
u3(-0.225013275825888,0.0,0.0) q[4];
cx q[13],q[4];
u3(1.63571935528580,0.0,0.0) q[4];
cx q[4],q[13];
u3(1.63678142904377,-1.32480662246980,-1.00764852445707) q[13];
u3(1.18440715078968,1.79850421463607,3.94911200207505) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
