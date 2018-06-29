OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.07241842766082,3.17442204939775,-3.02643557959715) q[0];
u3(0.618656991816202,0.255766326276967,1.36785735932551) q[3];
cx q[3],q[0];
u1(2.23947263659983) q[0];
u3(0.458407180881179,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.45888935806245,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.05411610732190,2.96069950008294,-1.14553738120497) q[0];
u3(2.11124295221501,-4.60626976104599,-0.199479961353931) q[3];
u3(1.37855307079145,-0.464336891970095,1.21064239656164) q[8];
u3(1.19906053698270,-2.83121382627884,-0.435086421500787) q[9];
cx q[9],q[8];
u1(0.365454514155481) q[8];
u3(-1.40613838381277,0.0,0.0) q[9];
cx q[8],q[9];
u3(2.23122120391214,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.894736085833451,-3.22607217806242,2.81859761201215) q[8];
u3(2.39722122194627,2.57139826760024,0.120423886783060) q[9];
u3(2.36849247648699,1.93953841326027,-2.18143501462908) q[5];
u3(2.11433758585629,2.34504731501558,-3.38658540773955) q[2];
cx q[2],q[5];
u1(0.806008110810858) q[5];
u3(-1.61173336839468,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.15633665007990,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.968685752094265,0.323327994878182,-0.995753493026589) q[5];
u3(2.19446559742989,0.0477373382838495,5.29511576157461) q[2];
u3(2.40166391976759,-0.398002396815928,0.0314333580503025) q[6];
u3(1.46947753085000,-2.84064645075012,-1.50616888078917) q[1];
cx q[1],q[6];
u1(0.0878139666651860) q[6];
u3(-1.54134312929610,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.21770531036568,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.64802616912807,-0.841796808130358,1.13154920086448) q[6];
u3(1.19115515503791,1.87804012710725,2.75452144149892) q[1];
u3(2.29516707245496,-0.319454256059904,0.637498530864971) q[4];
u3(1.12898507233904,-2.21730796664113,-1.86157802397979) q[7];
cx q[7],q[4];
u1(1.40115345604471) q[4];
u3(-3.55867677490630,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.26272280536617,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.16334787202929,-1.50176135169761,0.500867965248459) q[4];
u3(2.00443044345185,3.50615043686057,-0.694258158620696) q[7];
u3(0.865853317224867,3.00981772376027,-2.80291105023597) q[10];
u3(0.835840511221014,-3.85856255591116,2.32101877084826) q[2];
cx q[2],q[10];
u1(1.47062057893732) q[10];
u3(-3.34588605642969,0.0,0.0) q[2];
cx q[10],q[2];
u3(2.57479081269463,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.04162874398607,-2.10109578369464,3.19784286174574) q[10];
u3(2.21713670377124,1.97281614042403,-3.22798583352199) q[2];
u3(2.04364746707093,-1.17917467276058,-1.20282294197070) q[1];
u3(2.14520439639726,-3.55496753080777,-0.302814006303876) q[5];
cx q[5],q[1];
u1(2.35461897884507) q[1];
u3(-1.91176193991446,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.49986037569126,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.70491666303326,1.73464581500715,2.33849234850754) q[1];
u3(0.418602439381801,0.769720270503014,-4.52319458717621) q[5];
u3(1.51507938238853,0.679512726009088,-1.80236941099059) q[4];
u3(2.86846794236972,-4.58693316227737,1.16204265352627) q[8];
cx q[8],q[4];
u1(4.14097388818591) q[4];
u3(-3.36455958244296,0.0,0.0) q[8];
cx q[4],q[8];
u3(-0.0830401424901155,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.92406035341684,-0.0627532796494687,-0.718070097168115) q[4];
u3(1.01714442521834,-0.686354919917883,0.105320935477082) q[8];
u3(2.38004871952316,2.25305361196555,-3.17993624737489) q[7];
u3(1.07521056751688,0.190863565860064,1.20562100449999) q[6];
cx q[6],q[7];
u1(2.84511117059509) q[7];
u3(-2.35124532538349,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.63766495322253,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.47276607487265,0.157510660543970,2.42894903386240) q[7];
u3(2.29061035986456,5.38079764518590,0.511440224983117) q[6];
u3(2.27500822851573,1.56061343047623,-2.09913233588614) q[3];
u3(2.29679058140031,-0.0632547214373269,-5.79173420798098) q[0];
cx q[0],q[3];
u1(1.28116340604639) q[3];
u3(-0.864882052505814,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.75673882019946,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.50905826013410,1.97608209117610,-1.85743650151413) q[3];
u3(1.07204424711981,2.57754992188773,0.182654452853872) q[0];
u3(2.44826801119869,0.383425275850870,-1.65227392406078) q[3];
u3(1.49982789337928,-5.10429647859528,0.793957676386103) q[10];
cx q[10],q[3];
u1(2.26086709615618) q[3];
u3(-1.61420744377559,0.0,0.0) q[10];
cx q[3],q[10];
u3(3.86940207706692,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.08024671866943,-1.18387038106173,4.52463429621355) q[3];
u3(2.02789520031098,3.28890408355136,-1.81517173707339) q[10];
u3(1.48980562401643,1.20861526730956,0.363288905177269) q[6];
u3(2.61726929392152,0.522019905501530,-2.79545120627569) q[5];
cx q[5],q[6];
u1(1.55984489616214) q[6];
u3(-2.89429218139183,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.08048984860830,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.65397089817048,-0.498577534091041,-1.28571889727888) q[6];
u3(0.766140191975601,-5.52042634345297,-0.180373431896658) q[5];
u3(0.955039986879041,-2.25286066044586,-0.354113954265456) q[0];
u3(1.34159798782759,-3.60557066894332,-0.176643593720595) q[4];
cx q[4],q[0];
u1(2.83704218022839) q[0];
u3(-1.77192000119869,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.624072202704871,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.94048564546475,1.80084021501086,-2.15963385308312) q[0];
u3(2.21914848207078,-0.0978433787068611,2.83749292555259) q[4];
u3(2.21022598274737,-0.437122855016660,1.48890011407233) q[7];
u3(1.39179038112430,-1.83092603649455,-0.391336831967673) q[8];
cx q[8],q[7];
u1(1.46729132372940) q[7];
u3(-2.75174005488491,0.0,0.0) q[8];
cx q[7],q[8];
u3(3.21375424616477,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.30916629631813,3.20737693630770,-2.74369436637417) q[7];
u3(1.33994139716236,3.34845203839362,0.556038609220483) q[8];
u3(2.31673833050026,-0.682504597428622,0.119150356107895) q[9];
u3(1.50923993507399,-2.93865957561023,-0.400891126513029) q[1];
cx q[1],q[9];
u1(0.880106679187553) q[9];
u3(-1.37786060311340,0.0,0.0) q[1];
cx q[9],q[1];
u3(-0.424338298409600,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.88594741891902,-3.45818586289481,1.95576085909717) q[9];
u3(1.10661657519852,0.811643721334198,4.67459314798177) q[1];
u3(1.35400980612355,-2.19681756599555,1.87197519501528) q[10];
u3(0.317421379275485,1.41045275846125,-3.34246409625139) q[3];
cx q[3],q[10];
u1(1.43947522061078) q[10];
u3(-2.95721880433417,0.0,0.0) q[3];
cx q[10],q[3];
u3(0.466047485979779,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.18195647526400,-1.10502321760776,-1.74588618764457) q[10];
u3(0.903201328068703,3.15824862670264,1.46913744912906) q[3];
u3(2.29782004100513,-2.34545986677484,-0.0187489874085893) q[9];
u3(1.57688187320454,-2.93435249809967,0.701940168190659) q[8];
cx q[8],q[9];
u1(2.62122350111691) q[9];
u3(-1.87289692261090,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.403135444405315,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.26911762014598,-0.948531057318434,-1.05398491620296) q[9];
u3(1.68058580337465,-1.91741602003752,1.37100111624369) q[8];
u3(1.20485763597647,2.65645630775628,-0.815100455227093) q[2];
u3(0.943937649817680,0.218403816203629,-2.35729119090601) q[6];
cx q[6],q[2];
u1(1.73005971997268) q[2];
u3(-3.01054245211476,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.738579294878543,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.737292028056463,0.294358081262836,0.171344130990488) q[2];
u3(1.18317451252156,1.06399677263777,1.09010554563832) q[6];
u3(1.62035911546050,2.14481102436393,0.295419396173664) q[4];
u3(1.26965445613449,0.623196455459605,-3.74484118223478) q[5];
cx q[5],q[4];
u1(2.58451444987027) q[4];
u3(-2.10528775272424,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.355782178875697,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.259355383236640,-4.22904354332379,1.42779006331382) q[4];
u3(1.94456555326872,-2.26978496485625,3.87121902972607) q[5];
u3(2.49531046511040,1.49507334060703,-1.24598193436844) q[7];
u3(2.63358614901251,1.92190080508621,-3.17129032349841) q[0];
cx q[0],q[7];
u1(3.70580989467563) q[7];
u3(-4.57593966456604,0.0,0.0) q[0];
cx q[7],q[0];
u3(-0.630955480821170,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.763912152405150,1.63175886657770,0.0382086125971138) q[7];
u3(2.09245304585599,3.65622917486406,1.56549455908440) q[0];
u3(0.896953547889083,0.871118076653741,-0.154151888233360) q[5];
u3(2.00034669548127,-0.393212257600731,-3.33105397081704) q[10];
cx q[10],q[5];
u1(0.342029144915553) q[5];
u3(-0.673883617009606,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.40205417292388,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.748173334219226,-0.569133080686047,2.84007645716764) q[5];
u3(0.322786344261715,-1.15594395169551,-0.0654248156305525) q[10];
u3(2.25334784916053,3.35065671751379,-2.15981833033972) q[0];
u3(2.41723360283547,1.80186988486635,-1.62453653097904) q[4];
cx q[4],q[0];
u1(-0.409058999413300) q[0];
u3(-1.85684492950165,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.12818002395820,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.704947150614736,-0.740330664593599,-1.54645016576912) q[0];
u3(2.28291960013574,1.17546073769329,-3.80133988553293) q[4];
u3(1.32910470357932,-1.91377738870663,-0.270984049828033) q[9];
u3(1.46000068700774,-3.80482449197230,-0.588297653894838) q[1];
cx q[1],q[9];
u1(0.879374400618538) q[9];
u3(-1.17644848791881,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.92964735538325,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.603910774529578,1.94043625813069,2.52972690822391) q[9];
u3(2.28070239457245,-0.823432070500768,0.765757200066177) q[1];
u3(1.47970038825206,0.718260202485958,-3.27277394959456) q[2];
u3(1.21507881476484,2.83399006898595,-2.92990092915937) q[3];
cx q[3],q[2];
u1(2.32412362031953) q[2];
u3(0.217531541775201,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.27687956401770,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.516048437522779,-2.21724619871261,-1.30247347713477) q[2];
u3(1.50949257272474,-1.62277217242848,-1.80491896949229) q[3];
u3(1.97611286352430,-1.57112624166360,4.26646118809713) q[8];
u3(1.68458410739191,1.66882825789242,1.92258241841570) q[7];
cx q[7],q[8];
u1(1.21878178221670) q[8];
u3(-3.22256887120956,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.71891942079378,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.27397702683038,3.95394402560963,-1.21571553460812) q[8];
u3(1.67089938882942,1.15093149047071,-3.30250890472289) q[7];
u3(0.718405967851299,0.685701622639106,1.11334457949779) q[9];
u3(1.87195700345782,0.199784632423848,-3.35673378281022) q[7];
cx q[7],q[9];
u1(-0.170626738292579) q[9];
u3(-2.05563391761810,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.43185092153822,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.73707611510212,-3.37466048904146,0.960693778943168) q[9];
u3(1.47830269483980,2.55161615802290,-3.18316231917999) q[7];
u3(1.35384997887447,0.360376953381490,-2.21952475240315) q[0];
u3(2.01098772875337,1.71325419845550,-4.18164734251035) q[1];
cx q[1],q[0];
u1(0.899763207849404) q[0];
u3(-0.185474054427777,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.78716128902188,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.36607321268885,0.378734903555159,0.758690842235658) q[0];
u3(1.54461182092608,-2.91289828368775,0.682496833427870) q[1];
u3(0.510254139773609,-0.315398592071682,0.0972869840293401) q[2];
u3(1.05870942250797,-0.0850401700489877,-0.814849925225267) q[10];
cx q[10],q[2];
u1(2.59174129623994) q[2];
u3(-2.05706014936175,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.64557370904980,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.38612444997990,-0.298281352217735,1.97077021307700) q[2];
u3(2.64031906373796,-2.66371596335276,1.72823364109808) q[10];
u3(2.52022796290815,2.76357371902492,-2.30741173776536) q[3];
u3(1.29119534612111,2.46105731577338,-1.99136906851886) q[6];
cx q[6],q[3];
u1(2.86047432478021) q[3];
u3(-2.07514554784243,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.35526313001689,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.70086296285993,0.570885946131127,-2.35875408601760) q[3];
u3(1.40150633743253,2.22204213324963,2.18227965688207) q[6];
u3(1.80556306101708,0.784867301781339,1.83236509930052) q[8];
u3(1.61022684832793,-1.61712802499280,-2.11935230282519) q[5];
cx q[5],q[8];
u1(-0.0254511950478273) q[8];
u3(-0.615919826359786,0.0,0.0) q[5];
cx q[8],q[5];
u3(1.50435968418088,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.89706086195929,3.05845201588636,-1.77217752650475) q[8];
u3(1.61646811351970,-1.56698751635938,-4.52593445614755) q[5];
u3(1.90309715406292,0.472532514408234,1.16720158594622) q[7];
u3(1.95742057491999,-1.37311516059635,-1.91466937893926) q[8];
cx q[8],q[7];
u1(1.59918418447109) q[7];
u3(-3.08458961136642,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.09394227856704,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.810699011651758,1.96692063098727,-2.28120525225853) q[7];
u3(2.34408552644616,-1.60291287363063,-2.70901233594210) q[8];
u3(1.46537294387192,3.04204222558420,-2.36616301365517) q[9];
u3(1.82809229149589,1.15558079674422,-1.96717685552241) q[1];
cx q[1],q[9];
u1(0.705435180282966) q[9];
u3(-1.25547289556358,0.0,0.0) q[1];
cx q[9],q[1];
u3(-0.0253333813621781,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.93265113004004,2.10695360803088,-1.13990558750948) q[9];
u3(2.18040302835376,-0.960225948024343,-1.39778875086571) q[1];
u3(1.67176017737583,0.173423957725629,1.64860781479626) q[6];
u3(1.94452564635724,-1.86303847078259,-1.17889841136658) q[5];
cx q[5],q[6];
u1(2.95615209911845) q[6];
u3(-1.19737697336856,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.384091168851592,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.47323730367781,-1.12702057402027,-1.50832705801429) q[6];
u3(1.58278616323339,0.415225718601358,1.63319364946601) q[5];
u3(0.740270274644157,-1.42385533541723,0.423034745788200) q[10];
u3(0.115568917076613,0.358571875542945,-2.01448157356798) q[0];
cx q[0],q[10];
u1(0.280235974801189) q[10];
u3(-1.58633676141552,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.33001006009528,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.32088646721241,-0.101104419516617,-1.45088068685683) q[10];
u3(0.885762168311374,0.192071129375587,2.94099808432172) q[0];
u3(1.52099773489292,3.81092007211012,-1.37303638208505) q[4];
u3(1.54045428554531,1.77007556331842,-2.36831203221249) q[2];
cx q[2],q[4];
u1(1.15726176221925) q[4];
u3(-0.341339702693867,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.01113372231860,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.79139365509612,0.684328802020794,0.0504574089222519) q[4];
u3(1.15154092341509,1.05945801707553,-3.53068047053101) q[2];
u3(1.50144471529334,1.21181655401355,-0.818856341669282) q[7];
u3(1.38758449625891,1.32883000271934,-4.63054036299046) q[5];
cx q[5],q[7];
u1(-0.424360040134662) q[7];
u3(1.25240875796891,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.84210214053401,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.34212561071698,3.89574694933160,-1.98901976084452) q[7];
u3(1.19889451844822,-0.933756798184092,-3.64265018931242) q[5];
u3(0.700811310457782,0.972579669919495,0.212848000185727) q[2];
u3(0.883548880675187,0.0479344846525847,-2.21959274001676) q[4];
cx q[4],q[2];
u1(2.10463411302107) q[2];
u3(-2.43642344559245,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.181360790034486,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.93569556594221,2.30123665991878,-3.27766720774338) q[2];
u3(2.04087194537588,-1.85057904735505,4.07507936782906) q[4];
u3(2.18444262590943,-0.865822974192221,-1.98762654763384) q[6];
u3(1.84844166866897,0.846306331316784,-5.03425475344266) q[1];
cx q[1],q[6];
u1(2.29134188883308) q[6];
u3(-2.01720270725618,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.971671708249018,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.34572153972929,1.49587323860340,-2.41999783646292) q[6];
u3(1.09883215696661,1.98410457317426,-0.840404452856095) q[1];
u3(1.86493218352258,-2.41161245951773,-0.0381104121435005) q[9];
u3(1.27987829579750,-3.52161717477463,0.928455646512947) q[10];
cx q[10],q[9];
u1(1.47312850713613) q[9];
u3(-0.110357068634173,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.48082292085142,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.12614345205735,-1.81184649949087,3.13298757311458) q[9];
u3(2.80114661929741,0.707332679944836,1.06935348180041) q[10];
u3(0.851060619327910,1.42377684392567,-2.06201223284052) q[8];
u3(0.499308193421602,-0.327451995689877,-0.358473513306054) q[0];
cx q[0],q[8];
u1(0.727170988050348) q[8];
u3(-1.43005885226753,0.0,0.0) q[0];
cx q[8],q[0];
u3(0.291941193759026,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.45265876301377,3.89210597736012,-0.299238112367409) q[8];
u3(1.36575887485182,2.65217808840321,-0.366577023038641) q[0];
u3(0.805228830587744,0.375376196418909,2.03610409981173) q[6];
u3(1.18704287769781,-2.30122938865947,-1.82305356243701) q[2];
cx q[2],q[6];
u1(3.03615336841323) q[6];
u3(-2.43238288165579,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.36932884232025,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.97234003645958,-0.318067406177142,-2.48249451729351) q[6];
u3(1.09555668159205,1.37738760696804,-3.49062078739614) q[2];
u3(0.886613075130575,2.00565006013596,-0.790259948201681) q[10];
u3(1.48212779678077,-0.733743999890883,-3.22255044653188) q[0];
cx q[0],q[10];
u1(-0.182854999244244) q[10];
u3(-2.34471470690609,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.21857600665397,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.60438568068983,-0.558363181448716,-0.686466391899597) q[10];
u3(0.800446114839311,-1.05689724068269,-4.90474405441335) q[0];
u3(2.24745056557700,-0.125578182491476,0.521199339741991) q[5];
u3(1.05173444517192,0.128326044503565,-5.10823636492125) q[7];
cx q[7],q[5];
u1(2.03490003522729) q[5];
u3(-2.80237662122240,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.722645396794497,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.499261799444569,-1.03198541118174,2.77562589815444) q[5];
u3(2.83912413096403,-4.35001140406832,1.44327349480626) q[7];
u3(1.37562008335136,1.82644832209883,-0.423159701125425) q[1];
u3(2.58158150708592,-0.477189384514725,-3.74981936833071) q[9];
cx q[9],q[1];
u1(-0.0523120626792590) q[1];
u3(-2.13972889676935,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.36682816411978,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.883708150384068,-1.77735675829927,3.60655233987069) q[1];
u3(0.241739128085368,2.37235483781140,-0.369680174609948) q[9];
u3(1.61184324374517,3.45586750133976,-1.76023889342198) q[8];
u3(2.06372620183087,2.20092511095818,-0.746717643758326) q[4];
cx q[4],q[8];
u1(0.920704057910744) q[8];
u3(-1.30535867286052,0.0,0.0) q[4];
cx q[8],q[4];
u3(3.11511465227944,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.702946791228454,-2.24557712875699,2.74852632545219) q[8];
u3(2.36612226081118,4.34737450334527,1.75359450225438) q[4];
u3(1.39549686183243,-1.72972780734029,-0.648580659892383) q[8];
u3(1.43395640690261,-3.91222197284537,0.658052784533970) q[9];
cx q[9],q[8];
u1(2.34130884612621) q[8];
u3(-1.86890064284942,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.509078929394276,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.35610149621790,-1.00185273437114,0.575204631355997) q[8];
u3(0.959929883814717,0.125798081215322,3.52271955344680) q[9];
u3(2.41577545769723,2.07707179247468,-0.101518750011483) q[7];
u3(2.07993643098807,0.781092755827568,-2.86967226962725) q[5];
cx q[5],q[7];
u1(3.16606651353622) q[7];
u3(-0.936978399069004,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.65318969046944,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.14259355290341,-4.70262881224930,1.51807548356098) q[7];
u3(2.18585416922790,-2.58068124650482,2.82291468554185) q[5];
u3(1.87588508162309,-0.627337080215898,2.03458382321206) q[2];
u3(1.53388012001101,-1.34094942660231,-1.39260751974862) q[1];
cx q[1],q[2];
u1(2.28218353980111) q[2];
u3(-1.99283957455478,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.479562770058892,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.45911789647542,-0.710436999146816,1.32233360370070) q[2];
u3(0.722915707931140,0.235121772064752,-0.860718001947689) q[1];
u3(1.43134307256109,-2.17716572474908,4.04013069831236) q[0];
u3(2.21709301357191,1.12085079081489,-2.98101875884760) q[6];
cx q[6],q[0];
u1(-1.25475429509534) q[0];
u3(-0.266205894233822,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.61171156593769,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.45865636755419,0.221133305896579,-4.13202965743995) q[0];
u3(2.47642402383065,-2.78594432692294,-1.10810012766792) q[6];
u3(0.831275296473332,0.663347036315236,1.17865057451516) q[4];
u3(1.23947247478145,-0.550146726694555,-3.58319967311694) q[3];
cx q[3],q[4];
u1(1.21105705211605) q[4];
u3(-0.352179608191500,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.38011273412869,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.75998026419091,0.763768228338770,-1.18811993374937) q[4];
u3(2.07646713125931,-3.95011240482799,-1.11349148125030) q[3];
u3(0.285705112436837,-2.23357910150466,1.84759074025466) q[5];
u3(0.906596462707913,-3.77026148339675,1.96329178027202) q[2];
cx q[2],q[5];
u1(3.13409696242186) q[5];
u3(-2.26602498032129,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.16609423676169,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.60512859673948,1.75151965937089,-0.616827934999296) q[5];
u3(1.46530699304478,1.22293803037198,4.27358595864860) q[2];
u3(1.80683330966495,0.144786710240013,1.97296703134615) q[1];
u3(1.96801557520079,-1.85691619408819,-1.07095423701920) q[3];
cx q[3],q[1];
u1(2.71958066928955) q[1];
u3(-2.56370522388464,0.0,0.0) q[3];
cx q[1],q[3];
u3(-1.27255955338296,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.305987763336286,2.15096198278325,-2.53839865743002) q[1];
u3(1.37265054036045,0.313519140046614,4.88255028904939) q[3];
u3(2.08726143947340,2.65950334990683,-3.04521181017768) q[0];
u3(2.47497527213567,2.08674064413948,-4.11776405704386) q[9];
cx q[9],q[0];
u1(1.50170276144584) q[0];
u3(-2.00077425549600,0.0,0.0) q[9];
cx q[0],q[9];
u3(3.67938238914083,0.0,0.0) q[9];
cx q[9],q[0];
u3(0.862968260793491,3.15940718362929,-2.69924521273106) q[0];
u3(1.54698365247556,0.278778660547443,-0.960169057485328) q[9];
u3(2.40837516581756,2.82657477992964,-3.04444745265900) q[8];
u3(0.206098886531040,0.0799203684464764,1.98857373596149) q[7];
cx q[7],q[8];
u1(2.35598419465723) q[8];
u3(-2.05237671868095,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.17378525218999,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.08234294327710,1.30982919114986,-1.96155406124508) q[8];
u3(0.826727644477588,1.64606699170918,3.57541925921954) q[7];
u3(2.08175192848972,-0.539206680015916,-0.583522014115993) q[4];
u3(1.61150727401223,-3.23707933864239,0.330052443952351) q[6];
cx q[6],q[4];
u1(2.40374862689234) q[4];
u3(-1.64796843759519,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.196802588149172,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.35822731566520,0.617792767927114,2.38224223717900) q[4];
u3(1.65332451686234,-0.985428996283063,3.33711457376025) q[6];
u3(0.998783515353325,0.0337627360094156,1.75497245126922) q[1];
u3(1.60731390618491,-2.44533186929665,-1.51076364923863) q[0];
cx q[0],q[1];
u1(0.286319740403267) q[1];
u3(-0.942069688382274,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.69517360842731,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.30937101946160,-0.458751951337121,2.54967899246064) q[1];
u3(0.802790262394865,-0.220328848477965,-3.92289762852407) q[0];
u3(1.53798640600551,0.390204433486603,-1.90981751814792) q[7];
u3(1.36078926971497,1.34195618990805,-4.27581808670803) q[6];
cx q[6],q[7];
u1(-0.494676455413735) q[7];
u3(0.249116118666024,0.0,0.0) q[6];
cx q[7],q[6];
u3(4.19947587058066,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.10582800390252,1.23048215750737,1.53329812409112) q[7];
u3(0.697901344578645,-0.262082324238229,3.12630948182847) q[6];
u3(1.72964559537160,2.25485576740592,-0.0227181047144309) q[8];
u3(2.12935682300577,0.274455662005289,-2.00749980716892) q[2];
cx q[2],q[8];
u1(0.312325753462799) q[8];
u3(-1.21315124263109,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.0483384651102405,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.08860367914025,-1.66037331250481,1.50387725496488) q[8];
u3(1.42330274345571,2.64195937818753,-2.14784428882775) q[2];
u3(1.27727508586862,1.15798847456181,1.72539879492819) q[10];
u3(2.24934849005021,-1.74598844979647,-1.16144083342314) q[5];
cx q[5],q[10];
u1(0.0806412864849435) q[10];
u3(-1.50975188272696,0.0,0.0) q[5];
cx q[10],q[5];
u3(2.12662489252752,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.35569973550272,1.24476020310273,-3.31682301217470) q[10];
u3(1.11494627188884,-0.287733823773289,0.128140134595032) q[5];
u3(1.40610795808031,1.22476588991369,-2.36478251877270) q[9];
u3(2.81002832089593,2.41692386844420,-3.81551293270840) q[4];
cx q[4],q[9];
u1(0.346569329304096) q[9];
u3(-0.915328363083012,0.0,0.0) q[4];
cx q[9],q[4];
u3(1.54159638110359,0.0,0.0) q[4];
cx q[4],q[9];
u3(0.453902242788348,2.15456405855307,-0.271949939443764) q[9];
u3(0.888918964687596,-0.695418631483775,-5.47120194082585) q[4];
u3(2.19873616678844,3.05785675705358,-1.76618830064332) q[4];
u3(1.99781594831699,2.07117266308942,-0.464191030754219) q[0];
cx q[0],q[4];
u1(-0.0384972248338997) q[4];
u3(-1.84344085969089,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.548877353180826,0.0,0.0) q[0];
cx q[0],q[4];
u3(3.00162857474351,2.04488816932804,-0.311838375731249) q[4];
u3(2.50868662713593,3.14116776282623,-1.47943622315835) q[0];
u3(0.812924537881472,1.03467191454322,2.09328955204121) q[10];
u3(1.44024815449414,2.05511177647316,3.87294596425594) q[1];
cx q[1],q[10];
u1(3.29434726079080) q[10];
u3(-1.13841829417233,0.0,0.0) q[1];
cx q[10],q[1];
u3(2.52168845516532,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.01342909558261,1.90213885840834,-0.590355937094616) q[10];
u3(1.14614369883596,1.52272283409969,-0.557773975570313) q[1];
u3(0.539263836045565,1.83345110994636,-3.78174948514186) q[8];
u3(1.26025516430002,2.77699617080029,-2.81328121254505) q[7];
cx q[7],q[8];
u1(0.00815495895148977) q[8];
u3(-0.337768307979821,0.0,0.0) q[7];
cx q[8],q[7];
u3(4.05208143334114,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.62052634496283,0.863714107887728,2.82246115897658) q[8];
u3(1.80584676864807,-1.17221130175489,-2.61787727101259) q[7];
u3(2.02422998853257,1.82756387784564,-3.24143139794588) q[3];
u3(1.61984358371027,-2.09469732016056,2.33238813050062) q[2];
cx q[2],q[3];
u1(1.54436619480144) q[3];
u3(-3.52749341489765,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.18988417807962,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.17912256989286,-1.71330613967578,4.04225610973495) q[3];
u3(1.81217803582621,2.35439556517524,0.381793863392222) q[2];
u3(0.759996692330068,0.239355673711850,0.861553356423206) q[6];
u3(1.69214299352889,-0.515334562255345,-2.67433874145918) q[5];
cx q[5],q[6];
u1(1.01696953280741) q[6];
u3(-1.44965126023638,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.46486501648138,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.74265707847899,-3.71553220721609,0.0872246665988439) q[6];
u3(1.73826328303737,1.08188259869929,3.56559617954652) q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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
