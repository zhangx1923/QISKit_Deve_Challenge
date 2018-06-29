OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.48974707716694,0.690166946627752,-2.56800999454655) q[2];
u3(2.01396814765983,2.07057588045456,-3.87640883175363) q[3];
cx q[3],q[2];
u1(0.706550438305427) q[2];
u3(-1.21016322023568,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.75037119231614,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.93023475724851,-1.06823992127210,-1.26609135690772) q[2];
u3(2.40733058634696,2.98196082875357,1.10845808162785) q[3];
u3(1.05216939865968,-1.38731203172277,2.09541439148484) q[4];
u3(0.0169897997377369,-0.847810415634489,-1.08290804651764) q[0];
cx q[0],q[4];
u1(1.57057983278092) q[4];
u3(-2.22200011444678,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.22295198572267,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.998592092353915,-1.22594940250482,3.88934823764385) q[4];
u3(1.72129494836562,-0.372397810721908,0.110910418760743) q[0];
u3(2.34544681198406,0.110389206151832,-0.215986574574793) q[3];
u3(0.914063270710281,-3.20621127816922,-1.40025606839400) q[1];
cx q[1],q[3];
u1(-0.969605613146296) q[3];
u3(0.708777124865111,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.68093956456970,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.860728992569092,-1.37504708690306,-2.81776057039304) q[3];
u3(0.502206995333196,-3.21054604659661,-2.17370931360562) q[1];
u3(1.53337398225847,1.45423282624786,-3.87502580633928) q[0];
u3(1.16735685944307,-1.74593560314136,3.37595935808863) q[2];
cx q[2],q[0];
u1(2.20797131374346) q[0];
u3(-1.72507383130990,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.12534578660941,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.14851147014201,-0.00682779587698468,1.93383579629446) q[0];
u3(2.18259536202240,-1.36786038083390,-1.68287704037489) q[2];
u3(1.17692584009688,0.511875342767884,1.11873844048122) q[2];
u3(1.30531538566884,-2.62368464890835,-0.836980065422556) q[1];
cx q[1],q[2];
u1(0.0345767560670482) q[2];
u3(-2.68467797494720,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.22093144661658,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.49787812721899,-2.32254814568250,0.491069423021192) q[2];
u3(1.82343673138296,-0.612138945653832,-3.59565104664844) q[1];
u3(1.28572308631237,-2.09232950755164,-0.947472891873315) q[0];
u3(1.94606415200382,-3.52919497393903,0.332414733844064) q[3];
cx q[3],q[0];
u1(0.701855753904353) q[0];
u3(-3.39251947741846,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.95790252642385,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.11926173689274,-3.02533008131021,-1.41440410127492) q[0];
u3(2.33015407667781,3.27144557137100,2.82654289412311) q[3];
u3(2.47711302016170,2.92589513668704,-2.71938677668822) q[2];
u3(1.12927917641397,-2.48889239744960,2.84655214976480) q[1];
cx q[1],q[2];
u1(-0.225771899347924) q[2];
u3(-1.83408097502989,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.781734163613530,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.399176299958913,-0.531857794844874,0.812865442656098) q[2];
u3(2.37755490420620,-2.70719426242042,2.34754125927248) q[1];
u3(1.75347551083090,-1.02251361773946,4.01508033489139) q[4];
u3(0.314025905631818,2.11801016210573,1.93473449630686) q[3];
cx q[3],q[4];
u1(3.63315686020162) q[4];
u3(-1.18926547682084,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.61318622154498,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.38640185246634,-0.154746331362286,-1.71561219960300) q[4];
u3(0.564571708604895,-0.575117107981607,2.63489816454511) q[3];
u3(1.80622335727150,-2.79268410856767,3.08746131105173) q[4];
u3(1.06230871475558,2.87536646106550,-1.82800505597935) q[3];
cx q[3],q[4];
u1(2.36273953649774) q[4];
u3(-3.01207860271200,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.877743072045945,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.73836338205927,0.906163912633251,0.765312681169884) q[4];
u3(0.620120645703155,4.14056388475292,1.84958895312780) q[3];
u3(0.736286720965440,-2.57586418822235,1.66450323631396) q[0];
u3(0.721766853411577,1.06656162721727,-2.65433485423228) q[1];
cx q[1],q[0];
u1(0.617232320731828) q[0];
u3(-1.52375641622005,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.89509073828865,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.749978704232817,2.29352797287075,-0.298288319112386) q[0];
u3(1.75103241340788,-3.42858693343896,2.18998865006776) q[1];
u3(1.03644289793940,0.546099959440607,1.05754365759525) q[3];
u3(1.23609632414255,-0.983644936817271,-1.46825778732138) q[4];
cx q[4],q[3];
u1(-0.354233944954891) q[3];
u3(-2.07454005208199,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.61573717349874,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.969136400929646,-0.823812453818735,2.07955343380830) q[3];
u3(1.13865680020675,0.506308602316293,4.72030445215981) q[4];
u3(1.01897567404878,-0.205034048464856,1.59704662113094) q[2];
u3(1.07190588119731,-2.90455804775253,-1.10474405847050) q[1];
cx q[1],q[2];
u1(1.01509829095616) q[2];
u3(-3.60415540517467,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.67024730412179,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.06071022841077,-2.39225273004100,2.81743616838309) q[2];
u3(2.53569484413460,-1.84676492481716,-0.101559964701122) q[1];
u3(1.32803032176525,-1.56394003818331,0.974060509135545) q[0];
u3(2.14947217588017,-2.09851525846748,-0.0573277617760546) q[1];
cx q[1],q[0];
u1(2.36285310579682) q[0];
u3(-1.78958652470157,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.791745329577860,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.64586121767345,1.97386119827282,-1.33283354434520) q[0];
u3(0.437458227019999,-3.64365300709670,-0.296927618310313) q[1];
u3(1.09912535921067,1.56631713652008,-3.45865844522062) q[3];
u3(2.40835981633137,2.82475550668207,-2.81758085983665) q[2];
cx q[2],q[3];
u1(1.94768432074509) q[3];
u3(-2.47718851276815,0.0,0.0) q[2];
cx q[3],q[2];
u3(3.33995732457784,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.58057408552872,-2.49031266871765,-1.04188735883370) q[3];
u3(0.575517410842824,-5.81812085958640,-0.0249174723217527) q[2];
u3(0.602861789788085,-2.73169406561831,3.49561167755534) q[4];
u3(1.11550506356962,-3.41061086683641,2.16461122645819) q[3];
cx q[3],q[4];
u1(1.39588830273965) q[4];
u3(-3.31795367818052,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.25541064537431,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.576850530639128,0.806966948264200,0.251647981485975) q[4];
u3(1.64232794974473,0.480327319953226,-0.205855398170648) q[3];
u3(1.54713453241766,2.01281538979951,-2.63283255538898) q[0];
u3(0.977897910104120,-2.29808397934128,2.36270763669608) q[1];
cx q[1],q[0];
u1(1.90895476781172) q[0];
u3(-2.97656044590903,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.04243670977369,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.12493424626660,3.79172614285048,-2.07382964098253) q[0];
u3(1.34561564579770,3.34186405583186,2.29865149482563) q[1];
u3(0.894763259998937,-0.470654191182062,0.165079846524473) q[3];
u3(1.65765147100157,-3.44925348402616,0.361203546278340) q[4];
cx q[4],q[3];
u1(0.851807706002828) q[3];
u3(-0.196242460810029,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.60720783201183,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.86690979800637,0.420040352354797,-0.969199231414313) q[3];
u3(0.430811806260666,1.39995147991773,1.48924360253048) q[4];
u3(0.646538629588820,2.73450340004428,-3.28656243725655) q[2];
u3(1.33568822385780,-3.36328930778857,2.48914141069425) q[1];
cx q[1],q[2];
u1(2.14831118951570) q[2];
u3(-0.0448512713462488,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.65459953784738,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.441619244405479,2.66665824393834,-3.10151838053523) q[2];
u3(0.0225325407098567,0.743675876086018,-0.796758576507179) q[1];
u3(2.18856354818681,-0.761602533437483,-1.07382016994995) q[4];
u3(0.567086699434970,-4.08920255066018,-0.451957706241372) q[0];
cx q[0],q[4];
u1(2.19359389662059) q[4];
u3(-2.80248298040949,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.40983671527174,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.551187164071020,4.87514551162007,-0.713075552259925) q[4];
u3(1.06721522300824,3.49876622679207,1.42834410641968) q[0];
u3(2.26163596363363,0.356830986093552,2.38331379745212) q[1];
u3(2.31456993554866,2.15263367208739,3.83055492956797) q[3];
cx q[3],q[1];
u1(3.40878211244528) q[1];
u3(-2.16938444475547,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.0915648552143822,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.916173782378425,1.44305062949084,1.56857500460789) q[1];
u3(0.710355817906671,-0.238794860584002,3.50801141721528) q[3];
u3(2.19430032953655,-1.08490784255161,0.468199045052307) q[1];
u3(1.74166193756067,-3.23131855064731,0.165593310962239) q[4];
cx q[4],q[1];
u1(-0.465419118168291) q[1];
u3(1.09987050872944,0.0,0.0) q[4];
cx q[1],q[4];
u3(4.02388390435260,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.33261566242616,2.62917309127326,-1.44601727730220) q[1];
u3(2.61911952822692,-1.13141855130062,2.36710651394212) q[4];
u3(1.54413812266519,2.56868938877184,-3.67973393463567) q[0];
u3(1.02234018139251,2.82393248482572,-1.95921329627096) q[2];
cx q[2],q[0];
u1(0.667661165220816) q[0];
u3(-3.37358906992336,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.87806985462074,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.90286342110553,1.76076553067631,-0.213278385347249) q[0];
u3(1.28462712600471,1.64839994243021,3.93124482225574) q[2];
u3(2.08258456243194,3.53256279864473,-1.06223715668838) q[4];
u3(2.35020412569443,0.773265938200472,-1.35351672276287) q[3];
cx q[3],q[4];
u1(0.415500292425571) q[4];
u3(-1.16104328122873,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.68855773335114,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.52432526998215,0.376647081061978,-2.42428160944925) q[4];
u3(2.44887003419832,-0.712977372906930,-2.33385826941073) q[3];
u3(1.48360210063668,0.131474060556442,2.62685680873008) q[2];
u3(1.71304441170284,-1.31430837184277,-1.05931918504355) q[1];
cx q[1],q[2];
u1(0.404489786033843) q[2];
u3(-1.13725358206072,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.43940784936267,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.01173627065006,-0.909579160470912,-3.09460358919501) q[2];
u3(2.23669076403448,-2.42881819495701,-2.55186149630058) q[1];
u3(2.35128862155396,-0.111972516567689,-1.41775149204781) q[4];
u3(1.95433674609707,-3.64226389395699,1.14946005973739) q[3];
cx q[3],q[4];
u1(1.78308220744761) q[4];
u3(-2.20256583099613,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.271263073654878,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.81946207487106,-0.998352118920062,-0.944883381993573) q[4];
u3(0.493819728354605,5.13612790343001,0.318092882644435) q[3];
u3(0.168136680117316,-1.88329187463624,0.893264426091986) q[1];
u3(0.600016348198889,-2.54037756089510,1.34393038020793) q[2];
cx q[2],q[1];
u1(2.06814840053675) q[1];
u3(-1.60312117170667,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.73409502997545,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.242827232227526,1.95967488447259,-2.50518116905171) q[1];
u3(1.65311465639079,-0.353993162004737,1.25376102800514) q[2];
u3(1.90225418488305,0.472016019128080,0.00394234187232306) q[2];
u3(2.21976425570350,0.306568716474924,-3.79828758295105) q[1];
cx q[1],q[2];
u1(1.99947295527261) q[2];
u3(-2.44553838869662,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.597822519867684,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.08130469918150,2.09518614895166,-2.45904195504672) q[2];
u3(1.75932969213855,3.06604326077708,2.02391141322681) q[1];
u3(0.0915712134889429,0.516236746436721,0.110308001529629) q[4];
u3(0.713887417462371,0.854965909153778,-1.73792387439445) q[3];
cx q[3],q[4];
u1(2.91273461731202) q[4];
u3(-1.93122699270872,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.790725552538375,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.20761190178950,-3.51805358922986,2.39312258596238) q[4];
u3(1.23562330155252,-0.640323244949760,-1.62212858617233) q[3];
u3(1.21180358083428,-1.10054271384841,1.33901919405693) q[3];
u3(0.197669250353405,0.866804697214088,-1.29387156695785) q[1];
cx q[1],q[3];
u1(1.65907343769891) q[3];
u3(-2.23600735498139,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.05675306833800,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.36379212336196,2.36280560045266,-3.45806085725857) q[3];
u3(1.41689294016382,0.624835922986266,5.28044003633109) q[1];
u3(1.82923272532120,1.97220998109614,-2.74035226129043) q[0];
u3(1.35168860552232,-2.73705939329084,2.23155324713608) q[4];
cx q[4],q[0];
u1(0.216162375839837) q[0];
u3(-1.29415339121609,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.59995150243175,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.781058712199899,-2.31890413638096,2.27224438692157) q[0];
u3(2.31037877406152,-3.28344187850040,1.04135823471482) q[4];
u3(1.61227974611056,1.36935811009960,0.806015500713241) q[2];
u3(2.11824397858860,0.105364909181366,-3.29369272049980) q[4];
cx q[4],q[2];
u1(0.410407988550856) q[2];
u3(-1.13254576314080,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.31901098916769,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.07906836789453,-0.0890131318757781,-2.47271211695030) q[2];
u3(2.46279204404716,1.73269237884041,4.43575444034235) q[4];
u3(0.705660611786557,1.11946626554072,-1.65702896687624) q[1];
u3(0.814137124634134,-0.408247339898366,-0.608885543006131) q[0];
cx q[0],q[1];
u1(2.85689726107899) q[1];
u3(-1.82466924620462,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.52394509140569,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.51498307317960,0.616218093430773,2.90077764692022) q[1];
u3(1.90001671900537,2.21832654820885,-0.439525827838220) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];