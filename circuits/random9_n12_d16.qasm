OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.37057115650004,1.08348441720378,0.650427724549351) q[0];
u3(1.08996343237569,-1.58667423595171,-1.83383958313070) q[3];
cx q[3],q[0];
u1(-0.183441272739049) q[0];
u3(0.449638370604266,0.0,0.0) q[3];
cx q[0],q[3];
u3(4.37938682257852,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.41265000051293,-2.74064411797712,3.28987516291427) q[0];
u3(2.29534917877985,-3.89323121633153,-1.66343952321245) q[3];
u3(2.34926575390662,-2.01058127525785,0.143373731570989) q[8];
u3(2.77501772663228,-1.85512120524832,-0.295299192978452) q[1];
cx q[1],q[8];
u1(3.57631535170081) q[8];
u3(-1.58120548527879,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.92363961960691,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.06876619329446,-0.340050337152095,0.746522974097537) q[8];
u3(1.15924823389478,-0.922881794731554,-2.47106471654261) q[1];
u3(0.392927853719226,-0.711422419046526,-2.01225871134395) q[5];
u3(1.34995976037209,-3.43902029639493,0.0943649550195638) q[9];
cx q[9],q[5];
u1(3.46857106990889) q[5];
u3(-1.41973054800587,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.26186120157289,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.36275775725552,2.96492898753419,-2.03952992654405) q[5];
u3(1.86505864784920,-4.55081695033052,0.375986679985413) q[9];
u3(2.23915185000522,4.48241899314011,-1.46631322901920) q[7];
u3(0.685041256077241,2.58814484554788,-1.26028796203543) q[4];
cx q[4],q[7];
u1(1.03740692118817) q[7];
u3(-1.41873429865859,0.0,0.0) q[4];
cx q[7],q[4];
u3(-0.335436327737673,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.02807626034931,2.03856808248396,-3.27786291852418) q[7];
u3(1.38628207208311,-0.0322967292794270,2.28813842392506) q[4];
u3(2.66679456595964,0.483879476622790,-2.38566985621728) q[2];
u3(2.71339720281213,4.11068278715352,0.158155853879447) q[6];
cx q[6],q[2];
u1(1.70082602792130) q[2];
u3(-3.10481688504515,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.01135572824867,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.810761482105439,-2.17654550327824,1.97180622415210) q[2];
u3(0.455885332343350,-0.660218645094969,1.00354987120830) q[6];
u3(2.35152630815506,1.35310366985687,0.750282868489748) q[11];
u3(0.998940612784817,-1.63752533309992,-1.76753856140626) q[10];
cx q[10],q[11];
u1(2.03276355772233) q[11];
u3(-1.76670593644241,0.0,0.0) q[10];
cx q[11],q[10];
u3(-0.0100531811519426,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.42714799271605,-1.67947469795766,1.72935281193903) q[11];
u3(0.158757000817373,5.79538871363314,0.460961842678902) q[10];
u3(2.32229574580349,1.36526405214514,-2.55631479323085) q[11];
u3(1.60358529628210,2.31928143493174,-3.46174593376883) q[0];
cx q[0],q[11];
u1(1.03691424687602) q[11];
u3(-1.33512479128470,0.0,0.0) q[0];
cx q[11],q[0];
u3(-0.417216225519083,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.705866015692371,-1.45761808891258,3.18108021196325) q[11];
u3(1.82593809548042,-0.908296102461084,-1.32653866920697) q[0];
u3(1.25001744945288,-1.93522211140672,-0.205702799417079) q[3];
u3(1.73649710082957,-2.19112176848817,0.873918502643263) q[10];
cx q[10],q[3];
u1(3.59038262233843) q[3];
u3(-4.57950536502215,0.0,0.0) q[10];
cx q[3],q[10];
u3(-0.460782165142795,0.0,0.0) q[10];
cx q[10],q[3];
u3(0.413030213742481,2.91604950650750,-2.31920636534406) q[3];
u3(1.29508992120885,1.60228345127620,1.20602230530529) q[10];
u3(0.446260945034335,1.34588507025158,-3.00567490818759) q[2];
u3(1.38474551951447,1.89526184702954,-3.30870972177764) q[9];
cx q[9],q[2];
u1(-0.245204193734910) q[2];
u3(-1.75245974442809,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.801723883787362,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.26867306632452,-0.651126637027773,1.94409497512484) q[2];
u3(0.814274085277608,4.55929905318228,-0.365201262592013) q[9];
u3(0.796551044308609,-1.52739677816127,0.968173001377338) q[4];
u3(0.706718375860405,-1.25531902920637,-0.492566844038368) q[6];
cx q[6],q[4];
u1(1.49596676623115) q[4];
u3(-2.80884140997541,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.642663674436011,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.58374545431621,-0.836108050108963,-0.699496782776383) q[4];
u3(1.73415786010204,-2.30335954611627,-3.12202213317882) q[6];
u3(1.33008776111779,-0.359208565404065,1.44399535334315) q[5];
u3(1.57274327384143,-0.607365910732262,-2.25537140889693) q[1];
cx q[1],q[5];
u1(0.0637027371863566) q[5];
u3(-2.56383727183555,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.44221973058613,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.02156896643823,2.17896716728877,1.29257231309919) q[5];
u3(2.57109904433684,-0.232161122670161,-1.08584273784975) q[1];
u3(1.59144915652140,-0.409916555206406,-0.252601424629809) q[8];
u3(1.68953579862008,-2.62320236821709,1.20416179821742) q[7];
cx q[7],q[8];
u1(1.18744915855928) q[8];
u3(-0.871402304959449,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.181787388749980,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.17534712801460,-2.57287192235690,0.513023003025805) q[8];
u3(1.62706815298532,1.74286330393619,1.55351667286606) q[7];
u3(0.450299144131336,1.33929891681873,-1.76739959577458) q[8];
u3(0.919145786629551,-2.83704475631262,2.10696638043945) q[2];
cx q[2],q[8];
u1(-0.261464121684682) q[8];
u3(-1.54706521205372,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.93030904382899,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.45122007046265,1.09911620833310,-3.67330590905790) q[8];
u3(0.942781219981400,1.58744039443198,-0.860751411468862) q[2];
u3(1.10387633934346,-0.161493671188162,-1.28403498410777) q[4];
u3(1.37257432183098,-3.56704056574680,0.722138202639924) q[10];
cx q[10],q[4];
u1(2.80732326112544) q[4];
u3(-1.27281898816054,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.493297474552502,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.12095596695145,0.468502278044951,3.35639216872379) q[4];
u3(1.81553114223284,-0.212864979415259,-1.06910930997291) q[10];
u3(0.912250513593646,1.88751549437476,-0.676575928787882) q[7];
u3(1.41243383306471,1.97265031140924,-0.669756132989551) q[9];
cx q[9],q[7];
u1(0.412126502962885) q[7];
u3(-1.35849000524669,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.84474746283919,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.26529612436447,0.00602268788898708,3.58550282260316) q[7];
u3(0.868786574944264,-1.31720683516499,3.93608010520328) q[9];
u3(0.345390796415617,-3.04073872571679,2.22728980777744) q[6];
u3(1.24563976103381,-0.634944959057243,-1.61836563963318) q[11];
cx q[11],q[6];
u1(0.795941315889374) q[6];
u3(-0.168493809266961,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.80481909875548,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.72226474936324,-2.06658972717445,1.84243835493599) q[6];
u3(0.246743283810958,0.581058557451345,0.536047140393843) q[11];
u3(2.12617673719238,1.39032921638448,-2.27513262258008) q[5];
u3(0.780554949658817,-2.02909183942394,2.18045529320450) q[3];
cx q[3],q[5];
u1(4.40865410188556) q[5];
u3(-3.63633850466805,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.780178271372118,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.673071998209954,1.75891287574483,-1.71841778266535) q[5];
u3(1.34478106502268,-3.35886390517306,-0.572693031106382) q[3];
u3(1.48519623138766,0.556688642659251,1.19237187415742) q[0];
u3(1.81511640711295,-1.65455826161285,-1.01001277539369) q[1];
cx q[1],q[0];
u1(0.0159104942193786) q[0];
u3(-2.13601163019350,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.48791554415240,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.31325652344365,-1.85613317824387,0.414208226166368) q[0];
u3(1.27536881229180,-1.46079879902921,-2.50349170206728) q[1];
u3(0.391804010895455,0.743100703539808,-1.44297319099311) q[10];
u3(0.614999038376985,-1.72447284063546,-0.0920297460770771) q[0];
cx q[0],q[10];
u1(2.27572839649948) q[10];
u3(0.0105654351384239,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.982525088568110,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.70789597426473,0.00904056924961660,-1.09219403688437) q[10];
u3(0.867142051840371,3.28787490311065,-0.899037506492775) q[0];
u3(1.81682473649904,1.38620844187754,0.546623685535302) q[2];
u3(2.17437439859220,-0.478489910394254,-3.35125292106696) q[1];
cx q[1],q[2];
u1(1.29974813493814) q[2];
u3(-3.93010083250796,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.05795584538796,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.38672513054869,1.27599957011270,-0.323114467733839) q[2];
u3(1.52575727083919,5.22060677670603,0.261005986363708) q[1];
u3(2.64411198494981,1.87904109863511,-3.42396655751326) q[11];
u3(1.51779727315968,1.80651397759967,-1.25578319303824) q[6];
cx q[6],q[11];
u1(1.76813571166072) q[11];
u3(-2.47929942800236,0.0,0.0) q[6];
cx q[11],q[6];
u3(0.330804841621320,0.0,0.0) q[6];
cx q[6],q[11];
u3(2.44508939559912,-2.44278031728906,1.02299602403432) q[11];
u3(2.32887842014991,-1.32744469893296,3.37105616020191) q[6];
u3(2.15060920004955,0.102064095356835,1.99601334197610) q[4];
u3(1.70545590028581,-2.98560334617441,-2.27319346709996) q[9];
cx q[9],q[4];
u1(1.46731392590932) q[4];
u3(0.0187789026122225,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.56763151074088,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.390182388241199,0.0299402516705327,1.72678545215938) q[4];
u3(0.281263739846842,2.80565998504157,-1.13605436156922) q[9];
u3(1.26658783029748,-0.435568991030524,-0.744812939036838) q[5];
u3(1.20631238204915,-3.43890748433991,0.383077474547223) q[3];
cx q[3],q[5];
u1(-0.897595883187417) q[5];
u3(0.334946672231367,0.0,0.0) q[3];
cx q[5],q[3];
u3(3.62113431462710,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.75164000603628,-0.454407811635064,0.357122515325138) q[5];
u3(2.33546976354135,-2.81447042230287,-2.71668858582193) q[3];
u3(0.701913782540316,3.60048683995164,-1.39881020990954) q[8];
u3(1.74102688324569,0.909520469023915,-3.03294627894807) q[7];
cx q[7],q[8];
u1(1.80470949147960) q[8];
u3(-2.97172877335402,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.942611084516219,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.39765985768816,1.09238585436120,-1.37913464066152) q[8];
u3(1.06905137696732,-1.24478289560250,0.601485787854333) q[7];
u3(0.523694606799882,-0.149941916271439,-0.520602946531241) q[4];
u3(1.73230467657810,-3.60218896101467,1.96493998691659) q[3];
cx q[3],q[4];
u1(1.46836030581091) q[4];
u3(-0.795280075206493,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.105498662637589,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.06117404777351,0.604056666256076,-0.589229379833956) q[4];
u3(2.34330747718288,-3.40131710622672,1.47659411226135) q[3];
u3(1.57708831860549,0.930960044375242,0.557347273730547) q[1];
u3(1.48586926906318,-0.686510430588229,-3.09644646398371) q[5];
cx q[5],q[1];
u1(4.44330805054910) q[1];
u3(-3.45844778979912,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.557004523520786,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.96636235237426,1.15285276358134,1.18245106448646) q[1];
u3(1.66880393739773,-0.283880457828498,-2.70593857742214) q[5];
u3(1.90303804280641,3.09736860115907,-1.74028801540549) q[8];
u3(1.64672367166232,1.59892173502716,-1.66090601949731) q[6];
cx q[6],q[8];
u1(2.49642146123370) q[8];
u3(-2.72280141254908,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.43523176206789,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.24882334691013,-0.448864592964521,-1.40871865782249) q[8];
u3(2.22906379304343,-3.89441746196676,-1.86777785586124) q[6];
u3(2.23315075510861,-1.60219827532016,-0.535434240652529) q[10];
u3(1.92782583295643,-4.23884668427599,-0.256778393694591) q[2];
cx q[2],q[10];
u1(4.39811377566026) q[10];
u3(-3.69731997107992,0.0,0.0) q[2];
cx q[10],q[2];
u3(-0.742185800925220,0.0,0.0) q[2];
cx q[2],q[10];
u3(1.70475498043900,3.08847149149428,0.902268535196615) q[10];
u3(2.72627543207748,4.59874995601988,-0.494259706038153) q[2];
u3(1.15840485431819,0.341522750215251,-1.98891200573577) q[9];
u3(2.43386094078154,-3.95649721773548,1.61591134355907) q[11];
cx q[11],q[9];
u1(-1.16647736829820) q[9];
u3(0.543420591703938,0.0,0.0) q[11];
cx q[9],q[11];
u3(4.00606326737972,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.58223013098789,0.0563245691200523,-0.593000464867527) q[9];
u3(1.63898500023053,3.92020440075378,-0.295501678131125) q[11];
u3(1.12273940932030,-0.845351802534918,-1.12032573675388) q[0];
u3(1.14944153466296,1.29877239250916,-4.70060069483001) q[7];
cx q[7],q[0];
u1(2.51422983835337) q[0];
u3(-0.0152612475512539,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.999592224287588,0.0,0.0) q[7];
cx q[7],q[0];
u3(0.872893754741173,-0.727062112078470,1.17996563756259) q[0];
u3(1.94206029750744,-2.24195259151682,2.93645877608621) q[7];
u3(0.408260859898187,0.605170816086510,-0.904398731898347) q[4];
u3(0.881135802021652,-3.01475451593502,0.769999815199132) q[6];
cx q[6],q[4];
u1(1.89598325330635) q[4];
u3(0.490040903454473,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.858300622884447,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.94004114842430,-1.40817040288941,-0.600397111863323) q[4];
u3(0.549509482518084,-1.28667582971991,4.02369981944672) q[6];
u3(1.95481663059087,2.10904762818215,-3.74567468173395) q[1];
u3(0.760049153961264,2.90571416364437,-2.85790745852525) q[9];
cx q[9],q[1];
u1(3.01400703808956) q[1];
u3(-0.684501668762300,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.18520367736142,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.788267761724315,0.463619452294244,-2.08920487957679) q[1];
u3(2.04784935225010,0.724842774657795,-0.253146538474376) q[9];
u3(0.618052960485221,1.08050771968390,-3.31223215809856) q[8];
u3(1.22307367366379,2.95566189032256,-3.22590411255361) q[0];
cx q[0],q[8];
u1(2.36682058176538) q[8];
u3(-0.0877287522965879,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.07680557366178,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.45849341120272,-1.42553352219850,4.14558917925612) q[8];
u3(2.30008890176253,0.610127380638825,1.49721009785325) q[0];
u3(0.895645199225522,0.617932252785644,-2.87730952369760) q[3];
u3(1.16615087131359,-2.93887752680671,3.24713570777028) q[5];
cx q[5],q[3];
u1(1.14226535885756) q[3];
u3(-1.43366621845261,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.177127288406268,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.55741917563094,-2.35250093064285,0.164632664298044) q[3];
u3(1.59746770489603,-3.05886668629531,-0.432777109536755) q[5];
u3(0.910618976498081,-1.91775478171494,1.54195473775478) q[11];
u3(0.0520924600938254,1.73973147664824,-3.35835065781067) q[2];
cx q[2],q[11];
u1(3.11445134329725) q[11];
u3(-0.615594351750239,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.46921344816955,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.41707631522342,-2.71587108120072,2.00408746693304) q[11];
u3(2.33741570872797,-3.95011553783745,-2.06646207457045) q[2];
u3(1.83575818531950,3.70115832618393,-1.25059155652794) q[7];
u3(1.91379588970754,2.56533807998317,-0.367898121579825) q[10];
cx q[10],q[7];
u1(1.04563088507838) q[7];
u3(-0.541062625572951,0.0,0.0) q[10];
cx q[7],q[10];
u3(1.41078835012074,0.0,0.0) q[10];
cx q[10],q[7];
u3(2.69231200064472,-1.95366204556471,0.903248157413184) q[7];
u3(1.01659680438927,3.68665467765423,-1.00807857274941) q[10];
u3(2.76582519351525,3.44470773145063,-2.29125340639253) q[4];
u3(0.814614586846723,3.30754740927414,-2.62404181725764) q[5];
cx q[5],q[4];
u1(-0.259356043567689) q[4];
u3(-2.53928560124494,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.32291130500903,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.66414501068372,0.678977329302519,0.517947696663951) q[4];
u3(1.99231874819668,-0.714331909503707,2.45833190562008) q[5];
u3(1.59737263174337,-0.306936534400419,-0.549812187710960) q[11];
u3(1.44584439299001,-2.79837465241322,0.477544480316195) q[10];
cx q[10],q[11];
u1(1.11531202721456) q[11];
u3(0.189547301315775,0.0,0.0) q[10];
cx q[11],q[10];
u3(1.67270976663006,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.48219503432598,0.0890734583496982,-2.61990697601116) q[11];
u3(0.968935611790147,3.18322828256773,0.914272268689596) q[10];
u3(2.44541195311724,-0.810317046549634,-1.22849781748616) q[0];
u3(0.500652206822796,0.116433088916668,-3.90646546680964) q[2];
cx q[2],q[0];
u1(1.64875615512595) q[0];
u3(-0.550522538862553,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.68040347029576,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.923116790783854,1.21299620819190,-4.07417868647800) q[0];
u3(1.01603117214697,1.78562009782351,-0.153011125452309) q[2];
u3(2.48187157777737,1.25112414986922,-4.34065886197443) q[6];
u3(1.94130113048558,-2.35619359105560,3.67235183819977) q[8];
cx q[8],q[6];
u1(-0.0270108109147871) q[6];
u3(-1.53749838528663,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.818281029850144,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.47482419469051,0.505056960987204,-4.10935995323597) q[6];
u3(2.10760872178814,2.96110297956114,-2.85389371438115) q[8];
u3(1.02202606019613,2.53709743569686,-0.826853130592723) q[3];
u3(1.47923683237966,0.632450257092034,-2.53973716677501) q[1];
cx q[1],q[3];
u1(1.63803431726375) q[3];
u3(-0.203137819373445,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.56677510351672,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.79895790892008,-3.91549704877673,1.36269055134507) q[3];
u3(1.44360906313564,4.87097017485681,0.740356893285027) q[1];
u3(2.60933756465766,0.960470032357731,-0.306508449978553) q[7];
u3(2.11148072901905,0.838563085824060,-3.23250696197489) q[9];
cx q[9],q[7];
u1(1.84531176330054) q[7];
u3(0.295427061985551,0.0,0.0) q[9];
cx q[7],q[9];
u3(0.666422127131183,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.12061903343704,-1.44505044866247,-1.77378288741909) q[7];
u3(2.07563222873408,2.99148345376123,1.54182920358175) q[9];
u3(2.52173634664692,-0.367532861975345,-2.49176915553267) q[4];
u3(2.66215586644307,3.79016753760356,-1.04021572689150) q[3];
cx q[3],q[4];
u1(2.40367433002105) q[4];
u3(-2.90243007112785,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.15784799608986,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.07244946143852,-2.94497315251136,0.713868423410068) q[4];
u3(1.91145742963495,-1.14746691191028,-1.71943213743791) q[3];
u3(1.50978521348925,-1.32029191432752,0.993185634914208) q[11];
u3(2.02967900260883,-4.27964242219910,-0.480043595380432) q[1];
cx q[1],q[11];
u1(-0.194163805332522) q[11];
u3(-1.74422722152347,0.0,0.0) q[1];
cx q[11],q[1];
u3(0.701407874758739,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.50039697332876,-0.609633584213712,-0.403883343898546) q[11];
u3(1.64230255432866,1.22480649169144,-2.10929459496168) q[1];
u3(1.81414339728631,-1.66233708490805,0.439079761760429) q[0];
u3(1.95141977724037,-1.81086566119814,0.740625723938890) q[8];
cx q[8],q[0];
u1(2.13423826389728) q[0];
u3(-1.65654996496721,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.712125775825014,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.193055735436303,-1.96396779237425,-0.371248036634256) q[0];
u3(1.09095784950584,0.214242416115630,-0.330772919748774) q[8];
u3(1.62838655902000,0.391807621390916,-3.09700218016646) q[6];
u3(2.10818372294431,3.54913973819198,-2.34463560078500) q[7];
cx q[7],q[6];
u1(0.796224273217590) q[6];
u3(-1.18865349470586,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.99123918242077,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.37838112177175,-3.44585608748006,2.57595142083815) q[6];
u3(1.28094791356633,-5.37722346282957,-0.791227913158520) q[7];
u3(1.02946912114766,-0.758594777384425,1.78458604497132) q[9];
u3(0.989890332650705,-2.04118013268753,-0.716372844526742) q[10];
cx q[10],q[9];
u1(1.18472921638183) q[9];
u3(-0.0280259154969202,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.36566064819751,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.60640366985389,-0.364014812333187,2.01206296781246) q[9];
u3(1.68250739181850,-1.22896148899663,-1.59954430266920) q[10];
u3(1.69665831596709,-1.20003878642114,2.76201696039813) q[5];
u3(1.32865406882230,-1.38880146309170,-1.66291037252442) q[2];
cx q[2],q[5];
u1(1.62261844186732) q[5];
u3(-0.872134296896729,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.528424523252176,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.32612544835235,-3.46254792608076,-0.0119426080226437) q[5];
u3(1.90239876276350,0.923541095028081,2.44924238850650) q[2];
u3(1.62944153173461,-0.910297787050115,1.95068286580440) q[2];
u3(1.41558165938087,-1.72885373764921,-0.689055593611903) q[7];
cx q[7],q[2];
u1(3.08468542666576) q[2];
u3(-1.71451639475492,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.880337500551416,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.47824778843897,2.11705116100966,-0.882836298909054) q[2];
u3(1.97786532617561,-1.26671283462636,4.56152221301713) q[7];
u3(2.63492333553820,0.199348835024641,-1.85414278654806) q[1];
u3(1.36225305780731,-4.18071132720251,1.62539491096760) q[3];
cx q[3],q[1];
u1(1.32319810651344) q[1];
u3(-3.41211552029322,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.32483551299488,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.899266943375544,0.264158083698407,3.17828251051938) q[1];
u3(1.59637698716149,3.00830028062114,-0.666997162596463) q[3];
u3(1.21197025403007,-2.27298370862323,-0.320225830782535) q[9];
u3(2.39756702805899,-2.67271488954005,0.504845319832451) q[10];
cx q[10],q[9];
u1(1.46407444102230) q[9];
u3(-0.741523202323725,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.76659867562423,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.77270260666168,3.20231625509412,-0.132714649749490) q[9];
u3(3.07784880764747,-0.165837256612408,3.69697931992833) q[10];
u3(2.68994975092447,1.67560792389264,-0.193531732849890) q[6];
u3(1.95969918958695,-0.294051184786723,-2.45593258823870) q[0];
cx q[0],q[6];
u1(3.30918756068766) q[6];
u3(-0.423234617483473,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.80360012017193,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.21221854717273,-2.31519911731989,0.832973658580913) q[6];
u3(0.954116426366318,0.225528149723157,-3.44002739509118) q[0];
u3(1.47399059171609,0.0491937263486985,1.60104709306325) q[4];
u3(1.82153763093295,-2.86795722829591,-1.02305286368851) q[5];
cx q[5],q[4];
u1(1.15570285255407) q[4];
u3(-1.34476959294059,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.372930722113965,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.754411933393326,-1.14774563399945,2.53326507520460) q[4];
u3(0.605059591316873,-3.55720219749817,0.270074076113163) q[5];
u3(0.942408658352364,1.93383788585297,-0.826481075276150) q[11];
u3(0.562665692711395,-2.51590812460742,0.792057423460422) q[8];
cx q[8],q[11];
u1(1.53057940827294) q[11];
u3(-3.14199051524652,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.08027791777494,0.0,0.0) q[8];
cx q[8],q[11];
u3(2.18219281208704,-2.22022493375862,2.36029980066535) q[11];
u3(2.07915298797526,1.35826677269359,-0.0584966658249976) q[8];
u3(1.87126752751086,1.49991115956744,-2.62047273047369) q[10];
u3(1.28834820371271,-1.90218363437723,1.83624908804764) q[7];
cx q[7],q[10];
u1(1.75055056744108) q[10];
u3(0.491948420789678,0.0,0.0) q[7];
cx q[10],q[7];
u3(0.952000523293436,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.48298438988151,1.24676652204256,-0.504458141364285) q[10];
u3(2.10456241001468,-1.18658083386671,1.69518665690073) q[7];
u3(2.26134219877323,2.21665313125924,-0.387600600428141) q[1];
u3(2.44148928432487,0.140135916520389,-5.02494261785260) q[6];
cx q[6],q[1];
u1(2.38401272844607) q[1];
u3(-1.78623019720105,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.465789547225892,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.69573658037507,0.409782436391193,2.79645252524076) q[1];
u3(2.00486258958137,-0.578899582228051,-1.25159711086915) q[6];
u3(2.25528814617336,-0.889918545053552,1.86281403508015) q[2];
u3(2.21338442621694,-2.17983527982795,-0.816704466991607) q[0];
cx q[0],q[2];
u1(0.299714265718615) q[2];
u3(-0.485668660513981,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.68005936523595,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.38779773925373,-1.49299709026421,1.26409298414071) q[2];
u3(2.21699005765247,-2.26036641144618,-1.20561818649282) q[0];
u3(0.385355647053243,2.50454417466980,-1.94321308749720) q[3];
u3(0.857999949616405,-3.54826726944450,1.35048403267882) q[11];
cx q[11],q[3];
u1(3.13140449222482) q[3];
u3(-1.57558920682220,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.43264990458098,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.43129982329899,1.24730829033346,0.623030878223484) q[3];
u3(2.14154009419599,3.56872846846549,-0.319284372966298) q[11];
u3(1.26873141243201,3.05089146831891,-3.06256912450527) q[9];
u3(1.89915923859499,0.914242191143842,-1.65706228803705) q[5];
cx q[5],q[9];
u1(0.0984686817542804) q[9];
u3(-0.498579829056757,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.69361166582721,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.97142274170062,2.11135502067624,-1.07613855340307) q[9];
u3(2.46885631883210,2.17569438570538,-0.913975505155864) q[5];
u3(1.22221757397090,-3.47430499005917,2.70876446707476) q[8];
u3(1.68955299510059,-2.83351544787739,3.24026268275913) q[4];
cx q[4],q[8];
u1(1.74155601281389) q[8];
u3(0.464760787198300,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.956039782184674,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.56110686790400,-0.582577583703330,0.494955612996092) q[8];
u3(0.568563886365048,0.0679125072779665,-3.23318819180677) q[4];
u3(2.38600893309367,1.00831972937682,-2.45072205456279) q[0];
u3(2.50172586522257,4.52956810548655,-0.279949174310663) q[1];
cx q[1],q[0];
u1(-0.202985917848788) q[0];
u3(-2.08610725059213,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.733781196788965,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.68765386968823,-0.733413302408404,2.69275993124971) q[0];
u3(2.58369068629290,-0.363578345862250,-2.05126116186319) q[1];
u3(1.89076442220991,0.326502924920442,1.84355104303808) q[9];
u3(2.30313019634131,-1.21457452937451,-2.82279761438945) q[8];
cx q[8],q[9];
u1(2.99092765936761) q[9];
u3(-2.58531436831040,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.63005209563629,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.78260612402945,3.14128172267136,0.715345122017002) q[9];
u3(1.50084244102638,-4.77515629503639,-1.22154714693348) q[8];
u3(1.91188922850916,-0.957941602452877,0.848374159335225) q[5];
u3(2.44023222226388,-0.744459124637358,-1.80564174771903) q[7];
cx q[7],q[5];
u1(2.55806322863907) q[5];
u3(-1.91316482137598,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.322622238769105,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.778638132379218,-1.76631642456225,2.16072333876997) q[5];
u3(2.29553226763179,-2.23120072497389,0.489807598695439) q[7];
u3(1.19852575678487,2.02315207984248,-3.47762567253777) q[4];
u3(2.11340771309563,2.33169684023029,-3.42110344863070) q[6];
cx q[6],q[4];
u1(-1.17141422809255) q[4];
u3(0.807464623239026,0.0,0.0) q[6];
cx q[4],q[6];
u3(4.23377356528887,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.12764518593539,-0.313616464426502,-0.221159328332630) q[4];
u3(2.10274779626397,5.03226252862997,0.430745447870934) q[6];
u3(1.27337772631726,1.73974001336488,-3.85616498773290) q[3];
u3(1.66848911084790,2.20100314148898,-3.09447531467508) q[11];
cx q[11],q[3];
u1(3.37741172677005) q[3];
u3(-1.59478327048169,0.0,0.0) q[11];
cx q[3],q[11];
u3(2.24102646547113,0.0,0.0) q[11];
cx q[11],q[3];
u3(2.39553289871647,1.91819528392747,-2.37260879283102) q[3];
u3(2.83756849436987,0.850203396108444,-1.05350626737274) q[11];
u3(0.795430764054940,2.10428738458493,-1.82252796162735) q[2];
u3(0.445129469036047,-3.08981088627975,1.47046792580978) q[10];
cx q[10],q[2];
u1(1.90317148168844) q[2];
u3(-3.02584169744597,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.56114198643593,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.55545452505822,-1.63647437517951,0.541103974065521) q[2];
u3(1.18585734366142,1.20969141344656,-1.01799135677394) q[10];
u3(1.19427211027311,0.702148811783454,-1.04400416426178) q[9];
u3(0.655780253619823,-1.26795832557277,0.543903367495576) q[2];
cx q[2],q[9];
u1(0.858771101868784) q[9];
u3(0.140328455483738,0.0,0.0) q[2];
cx q[9],q[2];
u3(1.59589129866617,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.586895207286219,-1.73176444037152,2.36761644487756) q[9];
u3(2.32990636984719,-1.97350918116100,0.533434609019085) q[2];
u3(1.68423225739301,1.35439791511872,-3.46948708867852) q[6];
u3(1.30917594617589,-1.98180755728755,2.66462953334295) q[10];
cx q[10],q[6];
u1(1.53485110414700) q[6];
u3(-3.43509532950406,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.99211692332390,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.975669539285127,1.89155633127773,-0.152102061202802) q[6];
u3(0.968077024896112,-3.89338305125648,-0.195490552401730) q[10];
u3(1.45227653612661,-0.536445154372899,0.924460575379131) q[4];
u3(2.41854795989638,-1.55149782258879,-2.68898360737267) q[7];
cx q[7],q[4];
u1(3.15441542196306) q[4];
u3(-0.521523815148500,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.99640777171892,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.39949668080498,1.30090318835181,-3.70123198296393) q[4];
u3(1.27588738459601,-0.587020347028064,-1.34783930989367) q[7];
u3(1.99212212106790,3.16976966629601,-3.03690293535737) q[5];
u3(0.454943754292876,0.974142669070492,1.29510634582491) q[3];
cx q[3],q[5];
u1(1.81163357273510) q[5];
u3(0.373399268513940,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.08352701969434,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.867369359464168,0.968202251730337,-0.910668761540699) q[5];
u3(2.19045299638444,-3.16504895874754,-1.78794856856319) q[3];
u3(2.57943657351306,-2.15483578532113,3.66376916042483) q[11];
u3(0.849015917530801,0.0579026171997957,1.84121116051867) q[1];
cx q[1],q[11];
u1(3.06479181083485) q[11];
u3(-1.92377151681775,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.16052094470255,0.0,0.0) q[1];
cx q[1],q[11];
u3(2.69071849988366,-2.71031357678203,2.11073399405681) q[11];
u3(2.21945581115949,-4.42973989030151,1.28176773860471) q[1];
u3(1.21798583571669,1.53243704527092,0.0949608159823129) q[0];
u3(1.72024345536332,-0.241375862717388,-2.81128756387618) q[8];
cx q[8],q[0];
u1(-0.280213663221928) q[0];
u3(-1.86125360778476,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.840450104471708,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.25056610336006,3.25682587056484,-1.07125981417718) q[0];
u3(0.985527266554328,-0.510661384857250,4.75935120107119) q[8];
u3(1.40554665528254,-0.605150172344199,1.82930479886528) q[7];
u3(1.21047495459191,-0.831958002739319,-1.18100798751261) q[5];
cx q[5],q[7];
u1(1.80896218678968) q[7];
u3(-2.81856232015498,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.08982994272794,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.41610627748811,-2.89547972008296,2.19974276633011) q[7];
u3(2.17636052812791,4.23432331456167,-1.24766104075337) q[5];
u3(1.29476394041562,2.04642060387333,-2.54740464146262) q[3];
u3(0.353876379542660,-2.43854855664767,2.77573198306089) q[2];
cx q[2],q[3];
u1(0.788813841364750) q[3];
u3(-1.56641853428549,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.83247633233858,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.983078823622303,-1.37031692031428,4.44551539740126) q[3];
u3(3.07593284820097,0.0291609772552663,-1.57581879428570) q[2];
u3(1.88979557482865,0.777309350236280,1.59304440043575) q[8];
u3(1.70892739814067,-0.934582079137926,-2.25622055714365) q[6];
cx q[6],q[8];
u1(3.39171546202372) q[8];
u3(-0.549473970711545,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.79069651667892,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.06898625838295,-0.422939839829184,-0.551966969362329) q[8];
u3(0.899845288213640,-1.01976634625243,-4.96337336392166) q[6];
u3(1.96421665655104,-0.340366071804646,1.44002802552449) q[11];
u3(2.26182321327678,-1.13883135474871,-2.76623987399516) q[4];
cx q[4],q[11];
u1(2.71054386593833) q[11];
u3(-2.34361389062024,0.0,0.0) q[4];
cx q[11],q[4];
u3(1.11782010291177,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.28353804993776,-2.66882524518604,3.30766230501741) q[11];
u3(1.66362823370109,-3.73370757658414,-2.12576138129498) q[4];
u3(1.72912294386882,-0.191452383582237,0.811421090215582) q[9];
u3(1.60218715199829,-2.18299532563016,-1.51054348816744) q[1];
cx q[1],q[9];
u1(1.35123167950858) q[9];
u3(-0.703244158626299,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.114648956336744,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.17169420542389,1.94421755517765,-2.43055452912000) q[9];
u3(0.909630181992107,-2.97764712656085,-3.20281172092941) q[1];
u3(1.45938900406203,0.495647844364056,2.36130737259006) q[10];
u3(0.696364148780630,-3.12386190014484,-3.10612470381918) q[0];
cx q[0],q[10];
u1(-0.679351050628688) q[10];
u3(0.255394173928112,0.0,0.0) q[0];
cx q[10],q[0];
u3(4.25771898541088,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.84147558479985,-2.83172760735930,0.362858073294459) q[10];
u3(1.31636361320958,-2.42135667261286,0.372959138210131) q[0];
u3(1.37862976973761,3.24129577191238,-1.87473277191091) q[8];
u3(2.25357064615139,0.452888665053270,-2.41325536364891) q[4];
cx q[4],q[8];
u1(-0.0400949576717453) q[8];
u3(-1.81880349598800,0.0,0.0) q[4];
cx q[8],q[4];
u3(0.565560294178133,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.59408909338494,-1.69174772243341,4.55202953282475) q[8];
u3(0.329578880778482,-0.449851583994265,-1.85882060153604) q[4];
u3(1.85589229339418,0.405560289717097,0.658859817687097) q[6];
u3(1.38251817908877,-0.758320075488164,-1.59871400225709) q[5];
cx q[5],q[6];
u1(0.892313905748082) q[6];
u3(-3.12309780403278,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.88667108043041,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.69278114790380,-2.39056795585634,2.02719499304686) q[6];
u3(2.28147667263266,-3.65804509786584,-0.509705588997338) q[5];
u3(2.42790682985060,-0.0246059875207030,-2.67883605146685) q[9];
u3(2.31707823134695,1.34060688638341,-3.41469995016263) q[0];
cx q[0],q[9];
u1(1.81440272820764) q[9];
u3(-2.96047881830128,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.06117190701514,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.08063487136945,0.521764025976116,-4.26551481770009) q[9];
u3(2.40989905092643,-1.75588324003785,-2.35728098667132) q[0];
u3(0.996439888480669,2.44619453026806,-3.00949968485698) q[2];
u3(2.54219751108248,-2.03407858060334,3.53422663024037) q[10];
cx q[10],q[2];
u1(1.89921022872341) q[2];
u3(-1.74090059524549,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.800265787720278,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.48871962181102,-3.76112227515859,0.0702975738035139) q[2];
u3(0.397805003414153,2.37801031927076,-1.97068767535005) q[10];
u3(2.65578030313583,-0.709887547719915,-0.0995795951051051) q[3];
u3(0.594630428084770,-3.68614473484867,-0.752362686273093) q[1];
cx q[1],q[3];
u1(-0.531560432066486) q[3];
u3(-1.90308508408625,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.964870062202767,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.69163393576607,-3.92782635092846,1.35847757826782) q[3];
u3(0.383479841384888,3.81814400004130,1.04911909673492) q[1];
u3(2.06728582424482,0.451585963742105,-1.12378796827975) q[11];
u3(2.49587893095817,-4.33536586145291,1.24480058108530) q[7];
cx q[7],q[11];
u1(2.21769242242938) q[11];
u3(-1.93122238662738,0.0,0.0) q[7];
cx q[11],q[7];
u3(0.449065969559743,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.03001979886574,-1.22032613276457,-2.82299439115807) q[11];
u3(1.85054164587654,-1.50982189042408,-0.323698792179785) q[7];
u3(1.13136471689359,-0.414575255332584,1.60494622776391) q[8];
u3(2.00123926770862,-2.04701479868471,-2.62215005705827) q[0];
cx q[0],q[8];
u1(1.58266523088553) q[8];
u3(-0.173472011541763,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.60644441382056,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.01597191218793,-2.75961787578836,3.00699661947695) q[8];
u3(1.29994729477400,2.26616388163546,2.42846899754014) q[0];
u3(0.631026152048659,2.48991368891317,-3.06620867703643) q[11];
u3(1.91274956292923,-2.46232446639008,3.38484505984455) q[6];
cx q[6],q[11];
u1(1.54620489326666) q[11];
u3(-3.16434551562141,0.0,0.0) q[6];
cx q[11],q[6];
u3(2.71673776646589,0.0,0.0) q[6];
cx q[6],q[11];
u3(0.590941798811131,0.0133564490606559,1.69138241439890) q[11];
u3(2.13818730964134,0.555137965877801,4.51273065406958) q[6];
u3(0.692950571889642,-1.72590963535767,0.830775325895325) q[10];
u3(1.54041766552230,-3.10579620721579,0.822253493750399) q[9];
cx q[9],q[10];
u1(0.576943947354043) q[10];
u3(-1.00651695844967,0.0,0.0) q[9];
cx q[10],q[9];
u3(2.27611703672069,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.45764501393601,1.25835406536580,1.08061356294555) q[10];
u3(2.07410404813184,0.485977326418642,-0.549850947907499) q[9];
u3(1.93784180892518,2.18241515390014,-3.23001658150768) q[3];
u3(1.32532458983939,-2.73927406361610,3.42358909353073) q[2];
cx q[2],q[3];
u1(1.39855827994099) q[3];
u3(-3.65945971744007,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.22712366128994,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.44318765560598,0.652047780499795,0.123379536330701) q[3];
u3(2.52297076522468,2.34285829990639,3.65465346348573) q[2];
u3(1.94775540022658,3.31773055337110,-0.495811226710799) q[5];
u3(2.33829549967677,1.86394511225103,-0.983920133277717) q[4];
cx q[4],q[5];
u1(1.43065699526536) q[5];
u3(-1.12560264980892,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.63191478474745,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.48678406778794,-0.663543568629392,2.19396941620780) q[5];
u3(1.14362680231185,-3.20778842522229,-1.11023663309024) q[4];
u3(1.31218579969692,0.700917706337938,-2.53730328126053) q[1];
u3(2.12602576083261,2.27416334595827,-3.91659763968186) q[7];
cx q[7],q[1];
u1(-0.0790764330481923) q[1];
u3(-1.49104677535733,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.532861071263519,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.50925427669039,1.08233434394025,-4.67301618898876) q[1];
u3(0.942527018170617,-4.23932391003961,0.905096241466148) q[7];
u3(1.82861131286379,0.970342705242610,-3.54816244845650) q[1];
u3(2.21226220342122,2.77316009810243,-2.26646070675368) q[7];
cx q[7],q[1];
u1(2.25314532708258) q[1];
u3(-1.85458257727741,0.0,0.0) q[7];
cx q[1],q[7];
u3(3.59674797809113,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.40891338833232,1.36726752037386,-2.78083450296177) q[1];
u3(2.17606375882331,-1.11193727929745,0.756341222813852) q[7];
u3(0.268932015900073,-3.18794412956728,1.93743011275159) q[5];
u3(1.64552111520050,-2.14136672285827,3.37777850361552) q[10];
cx q[10],q[5];
u1(1.91722441084062) q[5];
u3(0.144022162117947,0.0,0.0) q[10];
cx q[5],q[10];
u3(0.654512938025940,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.58428181072468,1.57117310573876,-2.41818638925136) q[5];
u3(1.01711509855525,2.41383986358882,3.15485991017509) q[10];
u3(1.34311832887730,-0.678288310755274,2.28799865339793) q[11];
u3(1.42234321961063,-2.34365353439066,-1.77146759699378) q[9];
cx q[9],q[11];
u1(1.17712989074907) q[11];
u3(-0.770275516988663,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.320460553587716,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.14138337780272,2.00982800020765,-4.07452928504041) q[11];
u3(0.260697726558390,-0.765706958537101,-1.54227735472673) q[9];
u3(2.67479553644754,-2.89775696608317,2.72256840570622) q[3];
u3(0.817551972391628,-1.49737862945339,2.27225337241861) q[8];
cx q[8],q[3];
u1(-0.305523693477339) q[3];
u3(-1.86288647724195,0.0,0.0) q[8];
cx q[3],q[8];
u3(0.899089090145621,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.41623375657796,3.00611146790548,-0.593755590851070) q[3];
u3(1.69458040357507,-2.41386951126973,-0.599885774000754) q[8];
u3(1.81907882432456,-2.06851812395476,4.06203989413792) q[4];
u3(1.03225133055499,-0.954091437041758,2.93650737751329) q[0];
cx q[0],q[4];
u1(0.621293129897669) q[4];
u3(-1.49171638749852,0.0,0.0) q[0];
cx q[4],q[0];
u3(-0.157848572432907,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.78415286447976,2.87480179359027,-1.62126925230040) q[4];
u3(1.63081220316294,1.28376967833388,3.20225071984519) q[0];
u3(2.46062093255826,0.856660214855741,2.02743308562170) q[2];
u3(0.942167042592223,-2.24456725007513,-3.21076144611797) q[6];
cx q[6],q[2];
u1(0.188447957864096) q[2];
u3(-1.05969333939445,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.50766660778280,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.908855602495672,-0.178211369050682,2.75383106140190) q[2];
u3(1.01322227781053,-4.28804366601672,0.939095139073552) q[6];
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
