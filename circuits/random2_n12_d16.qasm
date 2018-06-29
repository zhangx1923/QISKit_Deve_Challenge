OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(0.369841055349706,1.17017244565248,-3.96590800678195) q[5];
u3(1.37895119012999,2.59293195735367,-2.70458199926135) q[4];
cx q[4],q[5];
u1(1.54599126783704) q[5];
u3(-2.50038132699743,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.261249596890287,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.01310769426337,2.87803573848609,-2.12939351296999) q[5];
u3(1.99182520829678,1.83314569057087,-3.51580163536906) q[4];
u3(1.51550435803691,-0.448360492834866,1.59228033805319) q[10];
u3(1.35137210590556,-1.21520335664046,-3.13614000633136) q[6];
cx q[6],q[10];
u1(2.79171349460186) q[10];
u3(-1.96077250700080,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.35545748900899,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.37591816018518,-1.99891067554790,0.910497573924928) q[10];
u3(1.92170470557341,-2.03295097691389,-4.02887644577502) q[6];
u3(2.10870198943066,-3.09128195586420,0.880563911249957) q[2];
u3(2.64730570475027,-3.74409859388389,-2.14842132731446) q[0];
cx q[0],q[2];
u1(2.12343387881236) q[2];
u3(-2.98206072822041,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.52349408564113,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.19376934766958,0.0621169020038214,-0.257067156096251) q[2];
u3(1.41311551401728,0.530784064217066,-0.264345675157557) q[0];
u3(0.258439064209699,0.960845291532242,0.00385886725592233) q[3];
u3(1.52930332566280,1.16804053787636,-2.63615967045228) q[8];
cx q[8],q[3];
u1(1.64589100109937) q[3];
u3(-0.262064196709662,0.0,0.0) q[8];
cx q[3],q[8];
u3(-0.0727849717408762,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.00300526932312,-0.468161086924254,2.23029601133510) q[3];
u3(1.73920628680840,2.16776812247339,4.07964200562525) q[8];
u3(2.48377943783389,0.0402581395866015,-0.364176710298778) q[7];
u3(0.509911927005823,-2.14347696106862,-2.10104029163364) q[11];
cx q[11],q[7];
u1(-0.444296830660484) q[7];
u3(-1.77661198759852,0.0,0.0) q[11];
cx q[7],q[11];
u3(0.814437413962150,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.44432589353097,0.317319812249098,3.02838831095939) q[7];
u3(1.11481142798689,4.36143233513361,-0.240474713561266) q[11];
u3(0.629832456025571,0.790088050399467,-3.23902575606241) q[9];
u3(1.53367931646549,3.01917740096596,-3.06900402469211) q[1];
cx q[1],q[9];
u1(0.00638397270986890) q[9];
u3(0.420500319019188,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.89730565621329,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.933063480835622,-0.477078352419069,0.678858424329026) q[9];
u3(1.10673019758927,-0.919073838743230,-2.78171304388783) q[1];
u3(1.81303576350954,2.27085678105846,-0.510194833313859) q[6];
u3(1.23302748906866,0.166848157251173,-3.83541772055999) q[3];
cx q[3],q[6];
u1(1.34193787161940) q[6];
u3(-0.289560894062779,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.02551134004580,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.13372036296189,0.676983046908171,1.26338697758167) q[6];
u3(1.72252580918071,2.82637432550178,-2.85276908706524) q[3];
u3(2.14646920266963,2.54833728529615,-1.98552811317569) q[2];
u3(1.07492573887714,2.07842752361061,-2.98774847624674) q[4];
cx q[4],q[2];
u1(2.90768112983499) q[2];
u3(-2.23933064289469,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.65491411936008,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.19398902412393,-0.807188878408176,2.85340746988648) q[2];
u3(1.25708525234355,-1.98691075304524,3.79506678079161) q[4];
u3(2.09046975977168,0.137637135122505,-0.478873717985914) q[9];
u3(0.809579622092253,0.490902425779895,-5.41575060578842) q[1];
cx q[1],q[9];
u1(2.68064975648813) q[9];
u3(-2.11572950115049,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.961233963598281,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.22788420536278,0.574135458429301,0.372475566678878) q[9];
u3(2.10275298047906,-0.353722703615873,2.70404638750686) q[1];
u3(1.50862403393641,2.31986076009038,-3.68114598724687) q[10];
u3(1.77472023291126,2.95020741964618,-2.54712598959481) q[0];
cx q[0],q[10];
u1(2.57449880521764) q[10];
u3(-3.04298033149866,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.17447827862235,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.02272692027860,3.86605563124097,-1.10224157490853) q[10];
u3(2.09206965657985,0.194165123640956,-6.00359991002429) q[0];
u3(1.22030803019337,3.73985794051640,-1.25165623540605) q[8];
u3(1.19124209264170,0.871875505224751,-0.753697857212736) q[11];
cx q[11],q[8];
u1(0.448097047548110) q[8];
u3(-1.03291830464902,0.0,0.0) q[11];
cx q[8],q[11];
u3(2.31660401382844,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.77922207120533,2.38803704618089,-0.613838096815219) q[8];
u3(1.82575528438024,3.96008258288864,-0.781733227397511) q[11];
u3(1.60307156865847,-1.32945396351440,-1.24205308523309) q[5];
u3(2.09290089400495,1.56318452001061,-3.66730535966084) q[7];
cx q[7],q[5];
u1(3.24828984612831) q[5];
u3(-1.53663990239781,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.24352737966092,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.42487412313073,-3.09516483813019,-0.856260607058657) q[5];
u3(0.587173469524928,-2.50832339057479,1.62113844170093) q[7];
u3(1.70535661324944,1.94422374498007,-3.34316301501310) q[8];
u3(2.26498948251939,2.02584214447397,-2.89675578162777) q[10];
cx q[10],q[8];
u1(-1.17468512928939) q[8];
u3(0.107073015622796,0.0,0.0) q[10];
cx q[8],q[10];
u3(3.63116410019699,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.56214298695399,1.95067129400082,-3.35802188801265) q[8];
u3(1.22912070995000,0.721132746186588,-4.23755120045553) q[10];
u3(0.294607374384655,2.40570371127658,-3.46834893177034) q[7];
u3(0.554099580966584,0.912216351936375,-2.15616835582126) q[11];
cx q[11],q[7];
u1(3.86222491705349) q[7];
u3(-3.41710184594735,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.994169273319736,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.19074393374068,-1.68897095913101,4.26867699430024) q[7];
u3(1.29345899568345,-3.47914355532795,-1.50022415089566) q[11];
u3(1.17036081790036,1.88545538030082,-2.20545007874487) q[9];
u3(0.503759884158570,-3.04452436994715,2.11054389834632) q[5];
cx q[5],q[9];
u1(1.32720327259822) q[9];
u3(-3.44382507675198,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.45230446825040,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.50498454210235,-0.0943606203909262,-0.785366678725257) q[9];
u3(2.24461628909308,-0.300202048115140,-5.78193074147976) q[5];
u3(1.08292556514677,-1.12483271180547,0.932509382738700) q[2];
u3(1.32379794869369,-2.67844322312018,-0.0556089820804095) q[4];
cx q[4],q[2];
u1(0.834611387565512) q[2];
u3(-1.33018148381575,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.0631058302216359,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.695330750039991,-1.39146954364586,2.45310888102481) q[2];
u3(2.64202848188363,-1.44670370752319,-0.826197494903482) q[4];
u3(2.51152991850994,1.49316663989803,-3.30556600906746) q[6];
u3(1.16508715181355,2.67593086033604,-2.80765901895646) q[1];
cx q[1],q[6];
u1(0.231120566111923) q[6];
u3(-1.64876577430295,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.514816380506685,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.12484238413629,-3.15571676951366,0.431733836347092) q[6];
u3(1.76323067952035,2.98053276037263,2.28051041011419) q[1];
u3(0.523191769938182,2.88601837966719,-0.756751701709230) q[3];
u3(1.59088273362973,0.491121571252974,-3.53867474206982) q[0];
cx q[0],q[3];
u1(1.29545635966750) q[3];
u3(-3.22336742121033,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.46872056900769,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.51904632908987,3.94629596201798,-0.430321413210544) q[3];
u3(2.43325768583498,-1.76478432168148,-0.987979085826004) q[0];
u3(2.64359067552122,0.417779093824075,-0.204440606381355) q[2];
u3(1.50932391744501,1.12440364580724,-4.67041729661917) q[9];
cx q[9],q[2];
u1(1.08938018953670) q[2];
u3(-3.37712654097333,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.94591337491652,0.0,0.0) q[9];
cx q[9],q[2];
u3(0.928581755468651,-0.282475264628110,0.509699592150940) q[2];
u3(1.49874312578086,-1.73644597183421,2.57714555448274) q[9];
u3(0.489857732953643,0.597315014210310,-2.85337441472088) q[6];
u3(1.24914842601541,3.41341018506818,-2.76418346536455) q[0];
cx q[0],q[6];
u1(0.453506359338946) q[6];
u3(-0.249690088515139,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.09920781761221,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.02067246960357,0.244306024641073,-3.83015191147673) q[6];
u3(1.44281437987819,0.680622447157121,-3.96316619265090) q[0];
u3(1.76547635038728,3.64663769128841,-2.56206949712195) q[4];
u3(1.91575488821595,2.14396588059417,-2.10734266643841) q[8];
cx q[8],q[4];
u1(2.63465083847174) q[4];
u3(-2.89890949019840,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.06277335062976,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.85867495704071,-0.887724389771705,0.816623535667752) q[4];
u3(2.87604597172745,0.782239078834951,-3.23046951653939) q[8];
u3(1.42184486589304,1.63613894179073,0.136682156712853) q[7];
u3(0.354232243145160,-0.526257109137182,-2.60888019484976) q[11];
cx q[11],q[7];
u1(-1.23311633916797) q[7];
u3(0.889038953883815,0.0,0.0) q[11];
cx q[7],q[11];
u3(3.85053637083898,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.41663790420969,-0.251707976487829,2.40372594435500) q[7];
u3(0.298715563773724,-2.23283446990073,3.98485106645094) q[11];
u3(0.378426387571106,-0.224831528668629,-1.97705116200116) q[10];
u3(1.50670687716567,1.39594175485254,-4.54653207706145) q[3];
cx q[3],q[10];
u1(1.27898618998300) q[10];
u3(-0.470253379107813,0.0,0.0) q[3];
cx q[10],q[3];
u3(2.08460364142533,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.36020081378110,0.357599093927279,-1.11752973894376) q[10];
u3(1.28151187093273,3.09856544930619,-3.01121739243567) q[3];
u3(2.79179632798212,-0.193576458211732,1.51636393980199) q[5];
u3(2.53939433412990,-1.39213104173503,0.171026620075973) q[1];
cx q[1],q[5];
u1(1.05072202424734) q[5];
u3(-1.56688470717090,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.601809184143436,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.00148497400439,5.10994980131713,-0.889261371335091) q[5];
u3(1.64878856402736,1.38419465650650,3.24503680339643) q[1];
u3(1.60485615094669,0.904842905182135,-2.45846842972565) q[11];
u3(1.31252674542067,-4.00679057364745,2.18771656877199) q[5];
cx q[5],q[11];
u1(0.863693190121204) q[11];
u3(-3.13109537844801,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.82316722847790,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.57501319662453,-2.15867319825189,0.128152085345842) q[11];
u3(1.06357366546394,4.74202960601167,-1.15757076664460) q[5];
u3(0.799952564989482,0.421288670273934,1.04736670965109) q[7];
u3(1.56582237556182,-0.394193082175325,-1.98101327740948) q[9];
cx q[9],q[7];
u1(1.49411663920614) q[7];
u3(-0.0188296496267712,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.51023262500347,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.75844852485197,2.01816808552206,-2.12959963679320) q[7];
u3(2.25080046794491,-0.825180374193518,-2.20395675005018) q[9];
u3(1.70576676667627,-3.64163186653490,1.89947084811916) q[4];
u3(0.338474963157898,1.65196406882780,0.306837675025874) q[1];
cx q[1],q[4];
u1(1.12735563679439) q[4];
u3(-0.686017754619787,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.79007556928228,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.934043146771777,-1.30710264170985,2.48387664151141) q[4];
u3(1.27548200752737,-3.11976351973256,0.802792578926399) q[1];
u3(0.990040474779899,0.391762875192316,-1.42428016269742) q[8];
u3(0.643827149577300,-0.569620082929952,-0.974155211368043) q[10];
cx q[10],q[8];
u1(1.64796386174281) q[8];
u3(-2.47181242832949,0.0,0.0) q[10];
cx q[8],q[10];
u3(3.30223395417671,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.52744396678760,-3.23608393973493,0.723413354355844) q[8];
u3(1.57260973766449,3.08073316040276,-0.259685237828911) q[10];
u3(1.53517399381905,1.60201705464605,-0.219867849793384) q[6];
u3(2.61413040285825,0.322584430336289,-4.19667791588432) q[0];
cx q[0],q[6];
u1(0.431742655350795) q[6];
u3(-0.980866712117311,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.87781108871713,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.29277914850354,-1.48413812509172,-2.94714602124905) q[6];
u3(2.32823620706992,4.74681761801360,0.988848684308023) q[0];
u3(1.23672131556518,-1.47205686588227,0.642413630490792) q[3];
u3(1.49312811590288,-3.88473582249999,-0.280732571262807) q[2];
cx q[2],q[3];
u1(0.881674546189801) q[3];
u3(-3.35920021777075,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.55062918645073,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.728374190001628,3.34751095809571,-2.63561504300454) q[3];
u3(2.46930969462072,-0.857385717922090,-0.836007713471972) q[2];
u3(2.11700435588798,-2.77981918224175,3.47682576751666) q[4];
u3(0.847964036831412,-0.379689223601950,2.06378230125313) q[3];
cx q[3],q[4];
u1(3.25392340647528) q[4];
u3(-0.924736876915791,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.74616338101640,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.84661436437763,-2.06913931184684,2.39582869047296) q[4];
u3(0.358239953737174,-1.06437984427140,-4.37265655247072) q[3];
u3(1.58385559038479,-0.786936023248473,1.30031817602316) q[9];
u3(1.32292456915019,-0.989484985427375,-1.10507936725244) q[5];
cx q[5],q[9];
u1(1.65542898435257) q[9];
u3(-3.71895944886400,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.18648021796804,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.00895937504450,-0.168346492378511,-2.36041065420191) q[9];
u3(2.80150448031097,-3.98469515988660,1.46721612011543) q[5];
u3(2.79079641567751,0.809730461845429,-0.511121001342888) q[0];
u3(1.56371208360941,0.347106438122755,-4.63478527686596) q[7];
cx q[7],q[0];
u1(1.50018905281144) q[0];
u3(-0.829681582531359,0.0,0.0) q[7];
cx q[0],q[7];
u3(-0.337708071035213,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.53799904816928,-2.17883314186832,2.35935646421508) q[0];
u3(2.39735251437602,2.06291934030187,2.97204144965790) q[7];
u3(2.17427419533337,1.65106723370200,-4.34989614629903) q[2];
u3(0.576392692712431,2.69754222049196,-1.74796283569198) q[1];
cx q[1],q[2];
u1(3.08124920479503) q[2];
u3(-2.04981194441822,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.49112123898829,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.00522568487960,-2.06895719473570,2.99401216788206) q[2];
u3(1.15493287226626,0.342021383889359,-5.08303809494959) q[1];
u3(2.31101212685476,0.347297140780359,2.37608310312980) q[8];
u3(1.10098968688753,-0.487894361568235,-1.27896255905573) q[10];
cx q[10],q[8];
u1(1.66641810789015) q[8];
u3(-3.62166862456228,0.0,0.0) q[10];
cx q[8],q[10];
u3(1.96487119358573,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.57538678051326,-1.55207806755500,2.05197796150463) q[8];
u3(1.15590344180136,-0.323363315581481,-0.430414423427864) q[10];
u3(2.48273709499635,-2.05589902492756,1.42658359625566) q[11];
u3(2.60013490624402,-1.16589309631742,0.933038535261291) q[6];
cx q[6],q[11];
u1(0.372394822886265) q[11];
u3(-1.70846119884324,0.0,0.0) q[6];
cx q[11],q[6];
u3(-0.221247649146682,0.0,0.0) q[6];
cx q[6],q[11];
u3(0.921352610855271,-0.372833945701787,-1.00117415219680) q[11];
u3(0.556785101849511,-2.92088387366333,-1.99045485869243) q[6];
u3(2.31878459103330,1.12232055232423,-3.93422441140795) q[3];
u3(0.998068841575302,2.96455270794761,-2.45629990339332) q[4];
cx q[4],q[3];
u1(0.903859832793645) q[3];
u3(-0.0778936099633594,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.31828121762469,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.85063092524104,-3.80017765109409,1.76559696740390) q[3];
u3(0.982865624222217,-3.16372547822208,-1.82867357939742) q[4];
u3(1.21221694036614,2.42291460625951,-2.12783218624308) q[2];
u3(0.686924115300705,1.89563596845541,-2.42563268144132) q[6];
cx q[6],q[2];
u1(1.68102577032646) q[2];
u3(-2.91317329458781,0.0,0.0) q[6];
cx q[2],q[6];
u3(0.970245641952316,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.65467436184038,0.395071581867425,-2.23378432052363) q[2];
u3(2.18866599312717,3.50218736008518,-0.135503667805731) q[6];
u3(2.20118800903046,-0.0604450817365569,0.711823081033495) q[5];
u3(1.71167703902528,-2.17296037612851,-1.21633054596840) q[8];
cx q[8],q[5];
u1(0.130007145013397) q[5];
u3(-2.29311753107450,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.64989826120185,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.34956350030486,-0.775902745109470,1.77691449016111) q[5];
u3(1.59915234384147,3.17407646776491,0.473425925162942) q[8];
u3(2.70760278295993,4.31856194171552,-1.89448393541861) q[7];
u3(1.32182260428638,-1.12353235311502,2.85978819836553) q[0];
cx q[0],q[7];
u1(2.14315557369663) q[7];
u3(-2.56133341211456,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.0530930697244421,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.65326682530542,1.09790749429638,-2.90371119174135) q[7];
u3(1.17267193948344,-0.301566580588041,0.361838151761866) q[0];
u3(1.97651651627471,0.104992882144521,0.875578037665333) q[9];
u3(1.79951172671785,-0.956501898035981,-1.78787218446875) q[1];
cx q[1],q[9];
u1(1.73253734259018) q[9];
u3(-3.05920941381217,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.762547375536998,0.0,0.0) q[1];
cx q[1],q[9];
u3(0.652180682029329,0.277378953523597,0.408805229546554) q[9];
u3(1.50793780268701,-0.572028881866013,-5.43002188279786) q[1];
u3(1.24933434348974,1.41521239368795,0.878310457743401) q[10];
u3(1.29355870514178,0.0267684974148865,-2.65837376729501) q[11];
cx q[11],q[10];
u1(2.67156093295083) q[10];
u3(-2.52815869114244,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.40152582697893,0.0,0.0) q[11];
cx q[11],q[10];
u3(2.62921861685514,-2.73847180889870,0.0155452450448601) q[10];
u3(2.11626260109289,-1.30093650559035,1.78953949726540) q[11];
u3(0.621309579387515,1.72073760188170,-1.46544206212200) q[8];
u3(0.852619185192134,0.194941976563006,-2.24533520495023) q[0];
cx q[0],q[8];
u1(1.60173298406924) q[8];
u3(-0.0156971196108389,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.08665696317455,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.41402351590949,-1.94279332810658,1.59086227432195) q[8];
u3(1.41289665240982,0.159734376339482,-2.97636763389922) q[0];
u3(1.72087019968410,2.23987671752621,-0.385089694942091) q[1];
u3(1.10536828440035,0.848628814137444,-3.69229984220476) q[11];
cx q[11],q[1];
u1(1.19620996366669) q[1];
u3(-0.995361004018017,0.0,0.0) q[11];
cx q[1],q[11];
u3(2.95637728115635,0.0,0.0) q[11];
cx q[11],q[1];
u3(2.28396779507699,-2.42511790468197,3.05147927633888) q[1];
u3(2.30033855231634,-0.542910236524752,3.96642219895797) q[11];
u3(0.928931686114574,2.25438669865635,-2.35401420387884) q[2];
u3(1.16871285768543,1.55892111873519,-1.38134211547867) q[6];
cx q[6],q[2];
u1(3.18263585741355) q[2];
u3(-1.19533461195835,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.33181221294760,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.45987793300079,-0.672975895104020,-0.802519645371431) q[2];
u3(0.981786890571989,1.64533956337042,-0.209825867995913) q[6];
u3(1.98832853420392,2.89178366494782,-2.75312615004187) q[9];
u3(0.635845919097827,3.39812609314843,-2.37411192712795) q[5];
cx q[5],q[9];
u1(3.28322583421260) q[9];
u3(-1.04207064926817,0.0,0.0) q[5];
cx q[9],q[5];
u3(2.02508354253860,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.29982695538111,-1.74607825764495,-0.246827896351759) q[9];
u3(1.29371866814839,-3.00066645681618,-1.68263037157943) q[5];
u3(2.50516871224936,2.55362400137861,-0.372723475776207) q[7];
u3(2.77540611309179,5.08063235858021,0.714984685863912) q[4];
cx q[4],q[7];
u1(3.47682317950697) q[7];
u3(-0.932068888321488,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.59538496139665,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.11688381101173,1.55127279931751,-0.626962540898582) q[7];
u3(2.05555249888291,-2.95086470576750,1.98366235469432) q[4];
u3(0.841091492527412,-1.81668917859142,0.878010427001047) q[10];
u3(1.48450662680322,-3.20191729998817,0.0766645533388524) q[3];
cx q[3],q[10];
u1(3.60085018299747) q[10];
u3(-4.24659086757319,0.0,0.0) q[3];
cx q[10],q[3];
u3(-0.712173674262430,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.54784784687688,0.567570272140326,-1.97730065006117) q[10];
u3(1.46947895555022,-3.99420534065984,1.51238929385350) q[3];
u3(2.66086199353121,-0.202351086792586,-0.144699792340706) q[1];
u3(1.30259763988453,-3.63875918049529,-0.320383105857376) q[5];
cx q[5],q[1];
u1(1.24087918203817) q[1];
u3(-0.232651631258907,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.37929211734235,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.34251404020617,0.189492128917894,-1.13830560418198) q[1];
u3(1.35167539144908,-2.45501074513554,0.848505512364050) q[5];
u3(1.19731396660742,3.99658039428363,-1.50839450275837) q[2];
u3(1.10945910275512,1.84593462453062,-1.18120718328967) q[8];
cx q[8],q[2];
u1(3.25538418379904) q[2];
u3(-1.94382692089500,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.922434767162000,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.48341287060381,1.05653384747392,-2.56907971522559) q[2];
u3(1.78201122281100,-0.241326844312521,5.92928933714462) q[8];
u3(2.48746383446694,0.193146773209633,1.83219356005423) q[9];
u3(1.46464009427027,-0.0919076192885226,-1.44541920674997) q[0];
cx q[0],q[9];
u1(4.14629533213420) q[9];
u3(-3.45142239202752,0.0,0.0) q[0];
cx q[9],q[0];
u3(-0.728147132820087,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.92264657046461,-0.176290876644961,-1.94476787290497) q[9];
u3(2.59354204238200,5.34098152871902,0.222021875113155) q[0];
u3(2.69365341586665,-3.45641142009865,2.74911837895852) q[6];
u3(0.854793958080886,1.79948335474044,-0.575120275053629) q[3];
cx q[3],q[6];
u1(0.159344358107883) q[6];
u3(-0.776465852215216,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.05453232405861,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.466443914859487,-1.70161806922253,1.91275075799097) q[6];
u3(1.97999376406322,-2.14014283984263,2.03849949882880) q[3];
u3(0.595631657834932,0.388729231730272,-0.681454537632943) q[10];
u3(1.28920178128781,-3.45109318569153,0.791823871219257) q[4];
cx q[4],q[10];
u1(-0.379559018396951) q[10];
u3(-1.55154120019075,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.07424390976812,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.75901865423136,2.08830424155820,-0.886727523346616) q[10];
u3(1.51004538558352,-2.01054756401650,-2.46163300854960) q[4];
u3(2.16113257503776,1.44300771104296,-0.383048315374928) q[11];
u3(1.36696250657686,0.330896504403025,-2.40205092037517) q[7];
cx q[7],q[11];
u1(1.85122314966116) q[11];
u3(-2.77243946016939,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.541277085924394,0.0,0.0) q[7];
cx q[7],q[11];
u3(0.284309812628876,3.77077791087983,-0.392939593649817) q[11];
u3(1.97243956464456,-1.62814487201728,2.76971829239976) q[7];
u3(1.62202537618568,2.19429969508556,-0.221138139310807) q[1];
u3(2.74476161817261,0.803323248951059,-2.20424190287595) q[7];
cx q[7],q[1];
u1(1.01948533648869) q[1];
u3(-0.0488936396819104,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.47824188031496,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.04090674593084,-0.275705957465991,2.86260894599169) q[1];
u3(2.34726588712075,-0.583558056741698,-0.374879732255739) q[7];
u3(1.83525465439546,0.285141675029584,-1.35016797070205) q[9];
u3(1.96783043045758,-3.35967858286813,1.40509346471928) q[10];
cx q[10],q[9];
u1(1.33639666382133) q[9];
u3(-0.509089998519531,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.39728194997501,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.864342824982679,-0.414576265268700,-0.654907835158454) q[9];
u3(0.244729687496463,-0.684554616867652,4.48116501141494) q[10];
u3(2.31341360666354,-2.37140582097499,0.0561034265674110) q[6];
u3(2.33667813955654,-3.46482543889707,-0.504867667152354) q[5];
cx q[5],q[6];
u1(1.30002441932772) q[6];
u3(-2.97670033607042,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.40277946708262,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.44409955057109,0.759258580194678,-1.54573860831877) q[6];
u3(2.34010841255908,-0.114504653708559,3.35185497286247) q[5];
u3(1.15133008979954,1.27591872079511,-3.39592591852567) q[11];
u3(0.898676308398254,-2.64390303789358,3.44624032569473) q[3];
cx q[3],q[11];
u1(0.0198009829773136) q[11];
u3(-0.860730982136727,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.90008432632409,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.51178421439740,1.40155713342061,1.24917523827028) q[11];
u3(1.85818749824004,0.459963830405124,-1.54786322206627) q[3];
u3(2.26329626652160,4.09991302325560,-1.64237223336660) q[2];
u3(1.01449591585796,2.53658607212736,-0.155467413413140) q[0];
cx q[0],q[2];
u1(1.76711255690634) q[2];
u3(-2.96296079955984,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.651071360426361,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.04946022885612,1.57932622946791,-1.33924356255966) q[2];
u3(3.05389523494963,1.25806742282830,2.40975893834136) q[0];
u3(1.91997661460601,1.97750457079139,-2.43286121979711) q[4];
u3(1.67852490943478,-3.10816413475075,2.59938529947913) q[8];
cx q[8],q[4];
u1(0.0249471283083467) q[4];
u3(-2.26827466282081,0.0,0.0) q[8];
cx q[4],q[8];
u3(1.18422602130516,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.84360065377923,1.55629516995777,2.06738658745895) q[4];
u3(2.52914102023295,-3.86671234014118,-1.41892768382476) q[8];
u3(1.54586199380826,2.60758444750461,-2.15186535819437) q[2];
u3(2.65073349959817,2.15467758635471,-0.576701956499182) q[9];
cx q[9],q[2];
u1(1.40410625761670) q[2];
u3(-3.37072029728003,0.0,0.0) q[9];
cx q[2],q[9];
u3(2.25074835266246,0.0,0.0) q[9];
cx q[9],q[2];
u3(0.501173018722908,2.17826207870741,0.339890668796664) q[2];
u3(1.89615218210180,2.75113911202984,-0.834019282695640) q[9];
u3(0.511722948591661,2.35177711731508,-2.15017570461094) q[1];
u3(0.732526253283221,0.289805298735325,-1.33821631980204) q[11];
cx q[11],q[1];
u1(1.50652484613347) q[1];
u3(-1.10966077679537,0.0,0.0) q[11];
cx q[1],q[11];
u3(3.11443612320772,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.06820628197971,2.31512703184637,0.656963933222341) q[1];
u3(0.619364103347500,-0.130321631567304,2.27830985984883) q[11];
u3(1.19701486087989,3.83602807342711,-1.09915464390385) q[8];
u3(1.61535887047890,2.20642094096517,-1.44268519520595) q[5];
cx q[5],q[8];
u1(-0.0977184684824783) q[8];
u3(-1.35548180818773,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.22703754641737,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.731213068719987,0.737060349455471,-3.76952939367366) q[8];
u3(1.42004556506471,-3.39316581443944,-0.235935528411348) q[5];
u3(0.978922005356707,-1.45506344666663,1.12077031701090) q[3];
u3(0.0930568558082549,0.0110136537392207,-2.38220212008073) q[10];
cx q[10],q[3];
u1(1.53257417472222) q[3];
u3(-0.326617784143941,0.0,0.0) q[10];
cx q[3],q[10];
u3(2.70702185467236,0.0,0.0) q[10];
cx q[10],q[3];
u3(2.14686201259783,-0.294415574362522,3.99876637450436) q[3];
u3(0.744870274334599,1.87623735993414,2.15005488725274) q[10];
u3(2.94380729130093,-2.38453316156181,3.65900708524690) q[4];
u3(1.20487533999698,-0.756390677370064,2.23138201314570) q[0];
cx q[0],q[4];
u1(1.63136188434963) q[4];
u3(0.353725072153583,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.991756834775024,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.45036863223108,-0.503979122388498,4.06208651478572) q[4];
u3(0.586400452976620,-0.677503216494755,-0.431445523989203) q[0];
u3(1.20679912199789,0.747314876400808,-2.76913028807260) q[6];
u3(1.90544186412911,-3.18708960702150,2.85126017006112) q[7];
cx q[7],q[6];
u1(-0.577794955116725) q[6];
u3(0.0990063659545808,0.0,0.0) q[7];
cx q[6],q[7];
u3(4.04137742395699,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.06679922996808,-2.57655827879898,1.19119401031385) q[6];
u3(2.78411260846251,2.21832859914992,-1.62007780603093) q[7];
u3(0.383862744594720,0.572788437360542,-0.125524570411923) q[8];
u3(0.986999851955482,0.311586695228056,-1.02524112245399) q[7];
cx q[7],q[8];
u1(0.327418603431150) q[8];
u3(-0.459690018810718,0.0,0.0) q[7];
cx q[8],q[7];
u3(2.09302808959456,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.41134622595073,-1.20861667797954,-0.440092176934436) q[8];
u3(2.38875521160163,-2.04901290866060,4.01602210634638) q[7];
u3(1.90172601141044,3.35011822426601,-1.05783075402648) q[10];
u3(1.48778215021277,3.32938722401466,-0.247345690660128) q[2];
cx q[2],q[10];
u1(1.03931885789564) q[10];
u3(-1.42204480303125,0.0,0.0) q[2];
cx q[10],q[2];
u3(-0.848514212608229,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.575713166292530,-4.69399235436636,1.20792235089179) q[10];
u3(1.02547518626545,-1.56436924736302,0.0983072039644912) q[2];
u3(0.185742084711487,-2.72630276900124,2.53678401872246) q[1];
u3(0.927734031501718,-3.03773408948265,0.651269730956756) q[0];
cx q[0],q[1];
u1(1.88162297614325) q[1];
u3(-2.88493913546642,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.865060081771304,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.40668068267597,1.50599027686198,2.14514629724348) q[1];
u3(0.993172372486006,-2.75410447832986,-3.09817962983139) q[0];
u3(0.889844441052776,-0.976020310211271,1.80287476513129) q[3];
u3(0.422574881330431,-2.15208440485890,0.439950965870081) q[9];
cx q[9],q[3];
u1(0.466459441915885) q[3];
u3(-0.279780576226268,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.66807283308357,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.90279111885174,1.05810241829103,-2.20264270443620) q[3];
u3(1.13632727014616,2.99907378864924,-2.22540583697017) q[9];
u3(2.95643362403775,-0.789897095576133,3.39119123680191) q[11];
u3(1.47703363109707,-0.610921849751583,-0.289962226168800) q[5];
cx q[5],q[11];
u1(0.323572584579885) q[11];
u3(-0.981233901970967,0.0,0.0) q[5];
cx q[11],q[5];
u3(2.57865436773351,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.51582714977259,0.0426261263783700,0.0514243530087953) q[11];
u3(1.17417204050902,-0.0729974647626211,-1.05885347580344) q[5];
u3(1.69140558161977,1.46909681358418,-0.685397373666232) q[4];
u3(0.455623925241576,-0.185842821908783,-2.87392471806161) q[6];
cx q[6],q[4];
u1(-0.00720453340403804) q[4];
u3(-1.34073800075882,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.09020105808956,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.879729330599558,-1.46302664546283,3.61897433225263) q[4];
u3(0.848377890295439,2.06128559517274,-1.25137276504409) q[6];
u3(0.177538806729071,-0.147963889577377,-0.476719762654884) q[1];
u3(0.751330061152086,-1.05558978897006,-0.928440372931707) q[6];
cx q[6],q[1];
u1(2.06751921050815) q[1];
u3(-2.92944721790878,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.38409135516262,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.67987973980714,1.84196789301035,-0.993006743428668) q[1];
u3(2.29189696029774,-0.940815180223541,-1.95866286443709) q[6];
u3(2.44772430074980,-0.853651456395752,-1.43830334391512) q[5];
u3(1.45088158921318,-4.74016102010642,0.652793359042582) q[9];
cx q[9],q[5];
u1(-0.213101681458691) q[5];
u3(-2.26127327424932,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.40955042051512,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.76129196283675,-3.68984419753170,0.173059990950574) q[5];
u3(1.87351567531311,0.872991077982075,-2.87477484955041) q[9];
u3(1.94758757343980,-2.42596548504288,0.307765721055232) q[4];
u3(0.633076657350304,-4.19760318705461,1.25625704327269) q[0];
cx q[0],q[4];
u1(1.99817630658632) q[4];
u3(-2.70361384920009,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.642231188634561,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.97695356516504,-3.65763538404886,0.0785196130592640) q[4];
u3(0.932778198014409,-0.806353346038391,-0.222025727765456) q[0];
u3(1.80344989504413,1.53941304103709,-3.72757782696732) q[3];
u3(1.47237896464714,2.01289694057475,-2.26617272213108) q[11];
cx q[11],q[3];
u1(1.28879799616063) q[3];
u3(-3.16384051754011,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.56783824584335,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.18590536470100,0.720359536572040,1.92183552839802) q[3];
u3(0.828026102148782,-3.58948827440088,1.13497650900461) q[11];
u3(1.70454596312396,-0.779681916567801,2.16578351101605) q[7];
u3(1.24925618908620,-2.02055847131637,-1.13763886741457) q[2];
cx q[2],q[7];
u1(2.91233549309605) q[7];
u3(-1.60750060372945,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.85097143863184,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.86817233939406,-0.160555185932843,-1.86343595174109) q[7];
u3(1.65778768842080,4.93330488939677,-0.994488555309768) q[2];
u3(1.42449661689285,-0.496161670451557,0.593327253486458) q[10];
u3(1.60592003786631,-0.699832708443597,-1.36562675038809) q[8];
cx q[8],q[10];
u1(-0.124267828356826) q[10];
u3(-0.887757095361635,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.57329694995939,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.98291430412727,0.951778991063669,1.59766745514325) q[10];
u3(2.33087844614200,-2.00986221565669,2.43630130218850) q[8];
u3(2.75543769062179,1.97989325473852,-1.98300719519452) q[0];
u3(1.92847159312582,1.69489105393616,-3.58319693552643) q[3];
cx q[3],q[0];
u1(2.18252795403106) q[0];
u3(-2.90645698923463,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.21104858954062,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.31480421202774,-0.593331556791463,-0.996045210249825) q[0];
u3(1.63050208001225,4.29215523624550,1.14277429918272) q[3];
u3(1.92043325673588,1.61339131512060,0.707647038385733) q[7];
u3(0.331782444104672,-4.63344188431885,0.150000230876437) q[5];
cx q[5],q[7];
u1(1.92098594442363) q[7];
u3(-2.78530138977647,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.61333098645258,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.981744458268783,-2.02295939605495,-1.82207113692022) q[7];
u3(0.244128417377273,4.41734863090828,1.75566489432381) q[5];
u3(2.59155172302586,-0.251640830503899,-0.103743659067458) q[4];
u3(0.973944713342707,-2.86899703120654,-1.33680607070797) q[1];
cx q[1],q[4];
u1(1.48144803312863) q[4];
u3(-3.51440363758366,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.93229730626767,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.32382269749595,-1.09135179809950,4.04043785901228) q[4];
u3(1.83841057714695,-5.36265804967105,-0.469463399678603) q[1];
u3(1.52383603650026,-1.31652799404512,2.03473362304716) q[11];
u3(1.12417833808260,-1.38464194415631,-1.31446871944711) q[10];
cx q[10],q[11];
u1(0.0364682007682444) q[11];
u3(-0.912257475157217,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.88484065500063,0.0,0.0) q[10];
cx q[10],q[11];
u3(0.779461197037104,-0.644915096093736,2.67960145709845) q[11];
u3(2.90665648129091,-1.58645379263566,-2.44418380784450) q[10];
u3(0.877254132709595,-0.454192545757148,1.69539751860415) q[2];
u3(0.598165670739294,-0.647245925708244,-0.829269595850511) q[8];
cx q[8],q[2];
u1(1.56713152815622) q[2];
u3(-0.345229580873564,0.0,0.0) q[8];
cx q[2],q[8];
u3(3.06910935691458,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.18746566062142,-0.715092026501832,3.18648310156144) q[2];
u3(0.360699436798063,-1.49822191470256,-1.88612960213249) q[8];
u3(1.48466162402748,1.48844913139631,-3.46985562424377) q[6];
u3(0.956422862735222,2.23490732459478,-2.58038472186297) q[9];
cx q[9],q[6];
u1(1.67345976609628) q[6];
u3(0.289888731295447,0.0,0.0) q[9];
cx q[6],q[9];
u3(2.03256567616315,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.517227712838438,0.0213417853810167,2.09198736494445) q[6];
u3(1.15821716819786,-1.31082441235687,-1.76443057977100) q[9];
u3(1.04308878057562,2.04098552206247,-3.77060862816131) q[0];
u3(1.44080849197547,-2.17970236990861,3.63456578185250) q[5];
cx q[5],q[0];
u1(1.61858108908740) q[0];
u3(-2.66910256202294,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.01742333169172,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.84693463381949,-0.0162744474469543,4.26746784661991) q[0];
u3(0.841466384268364,-2.08740933472169,0.0314463128494100) q[5];
u3(1.85174153540238,-2.41339799065527,-0.0392655803629924) q[8];
u3(1.27651138196221,-2.95110456816278,1.22340297854492) q[1];
cx q[1],q[8];
u1(1.06949281369414) q[8];
u3(-0.558231950678642,0.0,0.0) q[1];
cx q[8],q[1];
u3(2.89229354946333,0.0,0.0) q[1];
cx q[1],q[8];
u3(2.24519679305997,-0.847761992668128,4.77354120461314) q[8];
u3(2.64596519469924,-3.35515561813509,2.59909611011576) q[1];
u3(2.91843232893127,0.458862915602098,-3.20319875688022) q[7];
u3(2.65959516848489,4.04742952778867,-0.932613482756443) q[2];
cx q[2],q[7];
u1(2.35278014451644) q[7];
u3(-2.57837521540563,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.981507259532784,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.93559587296643,2.90243571006271,-2.40706190983690) q[7];
u3(0.150591911731519,1.31786632434960,2.31662203737652) q[2];
u3(1.60452019579898,2.40027398902558,-2.57520055974238) q[4];
u3(0.735494966868341,1.00616176767389,-1.83977443848263) q[10];
cx q[10],q[4];
u1(0.0335236850296390) q[4];
u3(-2.58135310117287,0.0,0.0) q[10];
cx q[4],q[10];
u3(1.05736685165347,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.27648541555283,0.495531719016478,-0.617147648071563) q[4];
u3(1.79478664644615,-5.14254125840814,-0.534811678711157) q[10];
u3(0.844200903738689,3.35206662870140,-0.634067900041421) q[6];
u3(0.657094347095457,1.55319555951593,-1.48190263487252) q[9];
cx q[9],q[6];
u1(2.75581485840458) q[6];
u3(-1.78988638200534,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.774012107186827,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.97243571365361,-3.27096919019620,2.88591982782720) q[6];
u3(0.850472131911376,-0.237657574417491,-5.00272859393088) q[9];
u3(2.47067038337673,0.477248584235429,-1.45067931367507) q[11];
u3(1.60843518103930,0.480728733046758,-2.95776854796352) q[3];
cx q[3],q[11];
u1(1.55814169130475) q[11];
u3(0.00372539452487364,0.0,0.0) q[3];
cx q[11],q[3];
u3(0.691279723654535,0.0,0.0) q[3];
cx q[3],q[11];
u3(0.266088414727119,0.566666076354743,-1.14842236677700) q[11];
u3(1.05761089769568,-1.68996785996860,-3.40418000097263) q[3];
u3(1.32087710951151,-2.46395070050408,0.0306433196845868) q[2];
u3(1.27185453328133,-3.61703914223949,0.950916948246558) q[4];
cx q[4],q[2];
u1(0.803970511733623) q[2];
u3(-1.07665867864871,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.26304534920412,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.13211357179068,-3.79214274819076,0.885961692671454) q[2];
u3(1.05576922550725,3.37405740979795,1.18779280522758) q[4];
u3(0.685039433381383,1.87182507462516,-3.01115145659705) q[6];
u3(1.65936677486111,-3.88504767059692,2.20663849234712) q[1];
cx q[1],q[6];
u1(3.26903533565884) q[6];
u3(-1.45135643585052,0.0,0.0) q[1];
cx q[6],q[1];
u3(2.31402707396139,0.0,0.0) q[1];
cx q[1],q[6];
u3(0.810608232614846,-2.45072207636065,-0.583259515195568) q[6];
u3(0.769732402591502,-2.93241847168504,-2.37305826374777) q[1];
u3(1.87565940692217,0.599078177016759,-3.27486663485541) q[3];
u3(2.14252351474148,-2.56344123136616,3.57557271736121) q[11];
cx q[11],q[3];
u1(2.16030506661504) q[3];
u3(-1.38694485887577,0.0,0.0) q[11];
cx q[3],q[11];
u3(0.114938559494689,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.42564255542907,0.462794847010714,-2.67478939663723) q[3];
u3(1.65674823474520,4.02969879860929,-1.92462496020671) q[11];
u3(0.815848968453840,3.16005806434129,-2.36517284894092) q[7];
u3(0.582666417606390,1.72005502352700,-3.45407717627371) q[9];
cx q[9],q[7];
u1(3.15923542519084) q[7];
u3(-1.38663421656356,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.64561789014236,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.29582632450758,-0.740364933196457,0.121206694955649) q[7];
u3(1.33647993606184,-2.37552227290837,-1.62267603750741) q[9];
u3(2.04766504544764,2.23721241508946,-2.99950322392845) q[8];
u3(1.05138843667569,2.37701981145552,-2.43852511249385) q[10];
cx q[10],q[8];
u1(1.21687588772100) q[8];
u3(0.180510239145387,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.613237887845784,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.86013112698721,-0.873117739221276,3.79574495422525) q[8];
u3(2.64865813608759,-0.854944460578198,-0.580424176137169) q[10];
u3(1.76285250925133,0.0162478182228978,0.817710340384833) q[5];
u3(1.52301478204012,-0.784024329657592,-1.30480443057524) q[0];
cx q[0],q[5];
u1(2.32079971570789) q[5];
u3(-3.18651843326173,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.20479768619438,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.70581196182208,0.373187669902807,-2.22008678098306) q[5];
u3(1.67925598613149,-2.54990747365739,-0.547530131781877) q[0];
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