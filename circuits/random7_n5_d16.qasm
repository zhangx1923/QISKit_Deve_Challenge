OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.40373722639649,3.49548691558260,-1.24725860060826) q[0];
u3(2.20642467949615,2.76194267292077,-0.526639490853636) q[4];
cx q[4],q[0];
u1(2.06919531077690) q[0];
u3(-1.70036827482460,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.95776888507130,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.59046842946780,-2.42703073300382,2.10293811446277) q[0];
u3(0.868729975125790,3.59035930674557,1.56148194302537) q[4];
u3(1.14985977864301,3.25437564506436,-1.72676499390327) q[2];
u3(1.93658641730088,1.23824990271045,-2.48589177563431) q[1];
cx q[1],q[2];
u1(-0.336891570372699) q[2];
u3(-2.70870159723037,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.56367784706275,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.30157065344117,-1.21526197523158,0.810698839789419) q[2];
u3(2.80205270931016,-1.99489405158167,-2.15422558003318) q[1];
u3(2.58116647211699,2.32616757731538,-1.56638314232184) q[4];
u3(2.54873633913338,-0.587046221391334,-5.03062800003972) q[1];
cx q[1],q[4];
u1(2.58486929761932) q[4];
u3(-1.87554958470364,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.642275140190788,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.25372012484050,1.40628913515737,2.27826508467817) q[4];
u3(1.38036020977765,-1.05845225975865,-1.42799549354516) q[1];
u3(1.12493912420959,-1.23224331357007,0.895177904896844) q[0];
u3(0.504100461431254,-1.88160986795366,0.632634881113856) q[3];
cx q[3],q[0];
u1(2.62916406746133) q[0];
u3(-2.92355833775895,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.92183202748027,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.459839825396649,3.26427589804294,-1.94858252550173) q[0];
u3(0.811352030189740,0.739435044365592,4.05605584469890) q[3];
u3(1.58794950501449,2.89638330214503,-1.47045525690284) q[1];
u3(0.859581312820578,1.53504184988000,-1.73881893578944) q[3];
cx q[3],q[1];
u1(1.20667801120150) q[1];
u3(-0.740090140165284,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.65355407633346,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.38726120391156,1.83245109094148,-4.34475273203524) q[1];
u3(2.41319134702065,-2.75927638525200,-2.98779389120107) q[3];
u3(2.56142201351535,1.72804074223638,-1.51180704419452) q[2];
u3(2.17942047993312,4.58897848142421,-1.02175040841901) q[4];
cx q[4],q[2];
u1(-0.452674499585974) q[2];
u3(-1.68481348951443,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.16140487046686,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.41079191260246,-3.37926550871951,0.886385351247973) q[2];
u3(1.34400875931876,0.394659640028950,-4.80108818703496) q[4];
u3(1.30747622381512,2.68372456504915,-3.08742327456655) q[0];
u3(1.88207411479215,-3.15061325294943,2.36377943371206) q[2];
cx q[2],q[0];
u1(2.45856378216439) q[0];
u3(-2.70910327019151,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.00805932481031,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.22750410728375,-0.0934782839347421,-2.05630880144641) q[0];
u3(1.58257491608943,2.24991572208342,-0.0456784023702244) q[2];
u3(2.82211159085766,3.59107961536558,-0.950756452070984) q[3];
u3(2.22895153931989,5.49735840978511,0.565193534385746) q[1];
cx q[1],q[3];
u1(2.55096422856593) q[3];
u3(-3.03988325389483,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.18878434001963,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.03476422810324,-0.344807850416866,2.57706077928197) q[3];
u3(2.36904301579615,-1.55276153029840,-0.651631710944004) q[1];
u3(1.95324801078634,-2.47434090979859,0.114300455050009) q[1];
u3(2.14029722957206,-3.87362948151548,-1.61002757545211) q[0];
cx q[0],q[1];
u1(-0.723649120913029) q[1];
u3(0.517989886565811,0.0,0.0) q[0];
cx q[1],q[0];
u3(4.19722217023660,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.90443960326587,2.17365597543040,-1.32227762517697) q[1];
u3(2.37428373873091,3.45593138445211,-0.582812207758017) q[0];
u3(0.405107296008622,2.26811243698283,-1.40226192192375) q[4];
u3(0.795350965767653,0.482679642210837,-2.11764871276301) q[2];
cx q[2],q[4];
u1(1.92056952933522) q[4];
u3(-3.07810892986860,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.76061023819479,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.17657930245383,-0.439874368406054,-1.19756475697567) q[4];
u3(1.45344979800079,1.74500517851147,1.04227566945250) q[2];
u3(0.516746947736506,1.77730508330958,-0.479944071427714) q[3];
u3(0.543316068534536,0.113183501289445,-1.47682018068793) q[1];
cx q[1],q[3];
u1(0.863776871337839) q[3];
u3(-1.35066195314646,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.44225516029317,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.836424615083662,4.00754649501002,-2.21580684863124) q[3];
u3(1.18078822845143,0.0274596027015517,5.54872961591565) q[1];
u3(1.90209803436285,2.95530583143754,-1.64101761853857) q[4];
u3(1.61012375620238,1.60903549867113,-3.03079494495701) q[0];
cx q[0],q[4];
u1(-0.290509535305062) q[4];
u3(0.116852506358126,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.93438966421287,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.60706668084433,-1.93639793855312,-0.0471013742582667) q[4];
u3(1.26797652523527,4.21127655982446,1.18778580766906) q[0];
u3(2.81686310736500,-1.60443849928392,1.31326925618078) q[0];
u3(1.53472909405941,-1.57101536691495,-0.712233432648921) q[4];
cx q[4],q[0];
u1(-0.244286350210584) q[0];
u3(-1.58559712584934,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.828125229057340,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.583961132363297,0.343185975985048,-4.00958228163416) q[0];
u3(2.13268704622522,2.45256224019085,3.25499165882973) q[4];
u3(0.828159191417828,-1.82594154749053,0.378620711273034) q[3];
u3(1.29130101816852,-4.23780386197766,0.501273901445740) q[1];
cx q[1],q[3];
u1(3.54136612721463) q[3];
u3(-0.659334637268289,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.69602532185641,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.48658031692018,3.71461177224601,-1.87148125937734) q[3];
u3(1.14608190680077,0.0159322570781664,3.63025355086156) q[1];
u3(1.63938727322661,0.668560961402675,-2.65112228209198) q[1];
u3(2.28742971306848,-2.50874259214164,3.32112536754981) q[0];
cx q[0],q[1];
u1(-0.829325293039292) q[1];
u3(0.334291884114106,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.85554071740335,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.722557490117581,2.78256651479143,-1.05649123021042) q[1];
u3(2.21655119952531,1.94738594628152,-1.89082535413339) q[0];
u3(2.01517486496648,-0.576578185312443,-1.90519185749942) q[4];
u3(1.19002734678978,1.31209000081745,-4.25743425534645) q[2];
cx q[2],q[4];
u1(2.77424479274216) q[4];
u3(-1.79712207217712,0.0,0.0) q[2];
cx q[4],q[2];
u3(-0.0473698393259991,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.65011691916854,-1.96347024698157,2.17550238612821) q[4];
u3(2.30444699536718,4.96998420002436,-0.989060196113868) q[2];
u3(2.80158320922183,2.58466381501150,-2.39027128627765) q[2];
u3(1.61689809806235,2.47621992049250,-2.98124638200624) q[1];
cx q[1],q[2];
u1(-0.107090520410549) q[2];
u3(-1.30226718793802,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.55330089185586,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.27332347540295,-2.61161882374468,1.53446846928154) q[2];
u3(2.27647533116321,2.92900130876308,-0.338669201215465) q[1];
u3(0.962297433343842,-0.899702926813956,0.770982019473049) q[0];
u3(2.14849136263530,-1.52409509672830,-1.92747007931082) q[4];
cx q[4],q[0];
u1(2.82317491866759) q[0];
u3(-2.14210976705104,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.144292301108193,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.688938544568575,-1.43805568639020,-1.39385294509855) q[0];
u3(1.18749199307851,-1.04109847548715,0.757120520153263) q[4];
u3(1.82525884696632,-0.0487300041503328,0.790389519829946) q[3];
u3(1.06722007554538,-2.60274396467787,-1.01721656946223) q[0];
cx q[0],q[3];
u1(2.86804283007403) q[3];
u3(-1.46674642860194,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.846504069767810,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.95703615466391,0.550712576487080,-0.113024562570892) q[3];
u3(1.29634214817197,0.152161866936191,0.621673774425517) q[0];
u3(1.20276538631081,3.23709422795691,-2.89830086241989) q[1];
u3(2.40234764027349,0.961328714552558,-2.11992151460398) q[4];
cx q[4],q[1];
u1(0.0354044246646263) q[1];
u3(-1.51933085823857,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.05261962397978,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.667631550649426,-0.266892137545727,1.92824306042274) q[1];
u3(1.32032954185305,-0.549043983504458,-2.72799657491962) q[4];
u3(1.79578097109028,2.06554191722614,-2.47943552907388) q[4];
u3(1.60058103956715,2.25670229584691,-3.73537325954055) q[2];
cx q[2],q[4];
u1(1.55774063281114) q[4];
u3(-0.374512874156169,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.93067148695965,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.22576119478078,0.411628801433300,1.80924956525937) q[4];
u3(2.08693799417652,1.37515431928136,2.46664737813247) q[2];
u3(0.919829379393188,-0.814830933687322,0.287383262721546) q[0];
u3(1.35237835160780,-2.27034760681275,-0.742778695573823) q[1];
cx q[1],q[0];
u1(2.33592431428972) q[0];
u3(-1.72193776025715,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.523462119976860,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.93532582882958,-0.810075039870802,-0.886236062534067) q[0];
u3(0.769677828418992,1.34830464303837,4.11036711542891) q[1];
u3(0.906769282995319,2.38688829544107,-2.59324569826041) q[2];
u3(0.817593278823193,1.02209725645272,-1.33604941894960) q[4];
cx q[4],q[2];
u1(1.79519865012244) q[2];
u3(-2.56717421615603,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.09247652815780,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.65927151005505,-0.300854963710130,-1.84580881400198) q[2];
u3(2.11152940177867,-0.319213398405622,-3.94264566306082) q[4];
u3(2.34367377390503,0.124452097624146,-2.96996097227138) q[0];
u3(2.61726440464912,-1.22554567650972,-4.81544412754963) q[1];
cx q[1],q[0];
u1(1.92489290709777) q[0];
u3(-2.47179743148262,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.24499150145445,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.68451321435553,1.84460694330991,1.21389820515261) q[0];
u3(1.34611012604317,-0.209426576970402,-4.53945032737323) q[1];
u3(0.576441526070426,2.22236237918522,-2.96110685412697) q[1];
u3(0.639269920956120,0.565883521295500,-1.85190514629349) q[0];
cx q[0],q[1];
u1(1.61009480598226) q[1];
u3(-0.805029570927981,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.210888640373444,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.26524406416041,-3.84402119752697,2.03283902726428) q[1];
u3(2.65403032866570,-1.30466785666748,1.34772614228419) q[0];
u3(1.64656346877404,-1.44459778862409,0.0408244292453598) q[2];
u3(1.19039927303732,-2.90291106916512,0.580632905565061) q[4];
cx q[4],q[2];
u1(1.70979435801715) q[2];
u3(0.121721421779802,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.818291537060163,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.80624001984701,0.338119082175450,1.14885991589891) q[2];
u3(1.25419777072845,-0.0971668009294135,-5.11737323930300) q[4];
u3(2.55271911161065,-0.481519164579794,-1.45029214791499) q[0];
u3(1.37997599445513,-3.54670492368043,0.887507661943719) q[2];
cx q[2],q[0];
u1(0.295395856810877) q[0];
u3(-0.997674051060564,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.26419610478745,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.15288187086859,-2.82132367731042,-0.983970867152507) q[0];
u3(1.43927847686278,3.55563930142731,-0.256119147397938) q[2];
u3(2.50148967779803,-3.03940139410309,0.276868385721324) q[1];
u3(2.42981538934604,-3.89385907927003,-2.23133282841328) q[3];
cx q[3],q[1];
u1(0.269027560123994) q[1];
u3(-1.44507632574642,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.22622983044590,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.29931478334158,0.627663185812370,-1.77161070653346) q[1];
u3(0.571803295967370,0.716690637820148,-4.26292368688645) q[3];
u3(1.20563033868031,-0.501514254203372,1.23186646212418) q[2];
u3(1.08281907996262,-1.36665180830712,-1.91555734701026) q[0];
cx q[0],q[2];
u1(1.37652517344790) q[2];
u3(-3.17657855735131,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.99589218397355,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.21394921250991,-1.62232090516806,0.607517854410529) q[2];
u3(0.208223649127764,4.05958280657301,1.21448079927639) q[0];
u3(2.27249377085603,1.51681919445012,-1.59817494191376) q[4];
u3(2.48924046860742,-0.342534582682585,-4.87675751123987) q[3];
cx q[3],q[4];
u1(1.94319007351410) q[4];
u3(-2.95407698525004,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.953782053716408,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.53345790015139,-1.26275038941958,0.636252046091580) q[4];
u3(1.47416764111996,-0.299729489822959,-2.32822971376754) q[3];
u3(2.55534637327434,1.05843644499413,-1.59161814270597) q[0];
u3(2.86204050048856,4.97797424109317,0.585198498624441) q[2];
cx q[2],q[0];
u1(0.812971054220810) q[0];
u3(0.0276462632215284,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.20330659368118,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.20699380607043,2.24058794428271,1.43458490895690) q[0];
u3(1.31388586787175,5.64440326459707,-0.349675079861169) q[2];
u3(1.35201127973410,3.35799236561229,-1.59962600753520) q[3];
u3(2.16808193330557,1.78589889018484,-0.202404649949965) q[1];
cx q[1],q[3];
u1(1.52230531861756) q[3];
u3(0.245970802591048,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.841168398867063,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.08285725778482,1.70563301349713,0.163612316652469) q[3];
u3(2.31971412140238,-2.26177865419327,0.967351654942080) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
