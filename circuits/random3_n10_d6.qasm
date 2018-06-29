OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(0.924362631064222,0.944467555478922,-2.17258706156648) q[9];
u3(1.76843045833805,-4.21669613333120,1.12480657992471) q[1];
cx q[1],q[9];
u1(2.78313759453584) q[9];
u3(-2.07457494177272,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.26762398540101,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.07449840337005,-3.02575313296726,2.58380815188688) q[9];
u3(1.71394649913245,0.811100400315232,-0.554097767018179) q[1];
u3(1.69729747507529,0.661467967569080,-3.57952748716355) q[7];
u3(1.08923478531626,-2.69414636042884,2.74796430471831) q[6];
cx q[6],q[7];
u1(2.99561613231985) q[7];
u3(-2.10153181202497,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.47731071104997,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.05558543552392,-1.64761606371892,2.76906595848178) q[7];
u3(0.864105489272334,3.29471711390449,-1.96919062142475) q[6];
u3(2.56843241967369,2.70193820684226,-0.911675140522124) q[8];
u3(2.23925222044244,5.47266214419624,-0.353741030904873) q[3];
cx q[3],q[8];
u1(0.876917265818889) q[8];
u3(-0.159643067372720,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.68673520650896,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.55965304774420,0.358158293030798,-3.29646514428827) q[8];
u3(0.681150732019113,-1.37317404352010,-4.10449938401797) q[3];
u3(1.05389417961652,2.09025961761227,-3.78191313211116) q[0];
u3(0.531629421493551,-1.98730958495824,3.15700059525653) q[2];
cx q[2],q[0];
u1(0.539120820829387) q[0];
u3(-3.38256499733098,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.76403415946693,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.57503625733575,2.85956007148226,0.481427127357524) q[0];
u3(1.71286955492413,1.36379957405230,3.52004210583582) q[2];
u3(0.710564328672450,-0.839773351214760,1.63958656499872) q[4];
u3(0.886714489355502,-2.07184820523943,0.771592559888584) q[5];
cx q[5],q[4];
u1(3.48015089213013) q[4];
u3(-1.35180510139747,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.06057821555887,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.29174841346727,0.165388285479168,-2.07196488059736) q[4];
u3(2.44125455928612,-2.71075924078951,0.256504510620296) q[5];
u3(2.94350398707790,-0.941067307430789,-1.08734552420800) q[5];
u3(1.40761060663014,-1.20692401359121,-4.26204143858888) q[1];
cx q[1],q[5];
u1(1.31551278316797) q[5];
u3(-0.00361072601448265,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.70606376462569,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.93315451267234,-2.18099054359445,1.30949955630611) q[5];
u3(1.23848498126565,-3.56243603578794,0.923080211033733) q[1];
u3(1.72452905671040,1.20364950402437,-0.730468195204802) q[3];
u3(0.937935982937840,-1.16508873670851,-2.56993984818468) q[7];
cx q[7],q[3];
u1(3.49072054661533) q[3];
u3(-1.08426530795735,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.98663390103577,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.689912444375920,-3.13915766140777,2.68743075170971) q[3];
u3(1.45411431884664,-5.06482603380432,-1.18192391116568) q[7];
u3(2.03171142114038,1.29244225471084,0.783510684600971) q[0];
u3(1.21753897809047,-1.51672302014494,-2.25508689729813) q[8];
cx q[8],q[0];
u1(1.50508930789750) q[0];
u3(-0.592262079056616,0.0,0.0) q[8];
cx q[0],q[8];
u3(-0.180279653448450,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.39058746231167,-3.06509848551019,2.09152703714044) q[0];
u3(0.778027505474644,2.00564185709269,-3.11798651460190) q[8];
u3(1.40675072318899,3.42557185840464,-1.56574121529151) q[6];
u3(0.511909406315640,1.36702735781955,-2.83204626384198) q[9];
cx q[9],q[6];
u1(1.77683020254657) q[6];
u3(-2.97023621582451,0.0,0.0) q[9];
cx q[6],q[9];
u3(0.983508639985325,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.612036391275716,3.23923614108254,-1.80420019356301) q[6];
u3(1.24530323473093,2.70597007571052,0.688083196340382) q[9];
u3(1.71837940111474,2.47765527827275,-3.39218160978395) q[2];
u3(1.05190448550675,-2.38713496821744,2.60354536376692) q[4];
cx q[4],q[2];
u1(-0.384116255635390) q[2];
u3(-1.80884242875980,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.691536356493619,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.806656978911312,-0.183844953715264,-1.79098658400029) q[2];
u3(1.17201393021149,3.02538264402818,-1.36485743423132) q[4];
u3(2.68637096423012,-2.50199469614180,0.619301245943185) q[6];
u3(2.13770437108109,1.88910615926114,2.21337104147771) q[3];
cx q[3],q[6];
u1(0.926532224363919) q[6];
u3(-1.49626040075432,0.0,0.0) q[3];
cx q[6],q[3];
u3(-0.432767866794962,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.80900781217750,0.394657428973283,0.177029370880323) q[6];
u3(1.38332689330054,-3.78143826870623,-0.179564265605955) q[3];
u3(0.835629111471917,1.94771475765482,-3.20624595337819) q[2];
u3(0.518215885215934,1.91551858229888,-3.01141732171717) q[9];
cx q[9],q[2];
u1(3.44270508210091) q[2];
u3(-1.10102210355551,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.53887580871466,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.14707863642081,2.69744419018405,-2.99975364806863) q[2];
u3(1.51616626504771,-1.61747239985122,-3.98601762687356) q[9];
u3(1.57741452220692,0.289468781437961,1.29978824996512) q[5];
u3(1.19733102466082,-2.31902710737556,-2.44745055359844) q[0];
cx q[0],q[5];
u1(1.97032384170396) q[5];
u3(-2.96374558756723,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.793325934626043,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.25668385429368,1.46120430568876,-3.32100216226428) q[5];
u3(0.628902752969804,-2.10487687770890,0.793073089562861) q[0];
u3(2.10707071031892,2.54156505904959,-2.99970296429132) q[1];
u3(1.26596593116177,-2.42429444414544,2.65172533581922) q[4];
cx q[4],q[1];
u1(3.65595621736349) q[1];
u3(-3.16908458029696,0.0,0.0) q[4];
cx q[1],q[4];
u3(-1.00982537568234,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.78494670595499,-3.32763025911235,2.10672345450160) q[1];
u3(2.35432545788622,4.36244631545021,0.505215461997174) q[4];
u3(1.29064439357133,-0.0608484882903901,1.94885114193159) q[8];
u3(1.92741928030224,-0.546940340828208,-2.43630009087472) q[7];
cx q[7],q[8];
u1(1.44498082110088) q[8];
u3(-0.0930646878136063,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.968340136342105,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.39312751580006,-2.92173625545317,2.05737433909275) q[8];
u3(2.74887651343697,-3.05002609635876,-2.28996497332973) q[7];
u3(2.60076609148135,-0.424859370528025,3.49565358514278) q[6];
u3(2.07956345892435,-1.17783134618796,0.817034797613744) q[1];
cx q[1],q[6];
u1(1.55733390561843) q[6];
u3(-2.30912385530268,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.04700312079671,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.73377729718846,-2.57125608187197,-0.342492215070024) q[6];
u3(1.14934759357456,1.15149498091220,0.945481422929416) q[1];
u3(1.45559197725283,0.672154399290389,-1.97619675086549) q[0];
u3(2.36393657741021,-5.49932001981767,0.765836916855136) q[9];
cx q[9],q[0];
u1(1.86225357302725) q[0];
u3(0.590284897673810,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.12849374492026,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.96291696790069,-1.35791353134491,4.69205845324757) q[0];
u3(1.60485713922711,0.0213533146994656,-1.50024478920623) q[9];
u3(2.30604837964429,1.14673419288547,-0.207843833017139) q[4];
u3(1.45249949601161,-0.465045633583543,-3.54267408094531) q[2];
cx q[2],q[4];
u1(2.92936515886167) q[4];
u3(-1.03689474119353,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.73492843905169,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.887607278667943,2.61788173640332,-3.06809657316490) q[4];
u3(0.874947903722879,0.547564054142814,-1.37703431861277) q[2];
u3(1.32306768809353,-0.825243322575900,-0.950655130644900) q[7];
u3(2.16900903299362,0.797126083335250,-4.38167865611451) q[3];
cx q[3],q[7];
u1(2.02361023938907) q[7];
u3(0.219821042119462,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.883464455277317,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.735691647733786,1.54275611613887,-4.67495477984876) q[7];
u3(2.14882850436704,0.742107473117591,2.55013609902406) q[3];
u3(1.57866390001741,-2.85521468720452,0.0534066424660928) q[8];
u3(0.555931424095437,-3.48101420751911,0.836179627386238) q[5];
cx q[5],q[8];
u1(2.11361537742139) q[8];
u3(-1.76866144013538,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.254722810314604,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.27476319237229,2.43976258945067,-0.259122661292004) q[8];
u3(1.40359965289918,-0.161682290122413,-4.94134832811341) q[5];
u3(2.09804424045376,2.06190888196018,0.347057570472961) q[8];
u3(1.02524894085267,0.274428765745157,-2.65416887843096) q[4];
cx q[4],q[8];
u1(-1.12573251359208) q[8];
u3(0.481649168467947,0.0,0.0) q[4];
cx q[8],q[4];
u3(3.86425142537311,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.12562422141119,1.02500605033940,-2.97857459374814) q[8];
u3(0.856438667679384,-0.0792466502705971,-3.20272936246063) q[4];
u3(1.41407978117652,-2.92111622458881,1.69715788275965) q[9];
u3(0.357097817847194,0.917719170956715,-2.64881047136247) q[5];
cx q[5],q[9];
u1(1.48656920077072) q[9];
u3(-3.08634789499474,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.23302722434363,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.35439914855692,-2.19748275473144,3.71634650240635) q[9];
u3(1.45905077894479,1.66251586270063,2.34328975890638) q[5];
u3(2.27861997296288,0.279913920603450,2.83159833701080) q[3];
u3(2.57702086171993,-0.950558581143256,-0.159145441494415) q[6];
cx q[6],q[3];
u1(1.46410665188650) q[3];
u3(-3.40352451913495,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.04749212076650,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.339381099913900,-0.304733410490285,-2.26236286928572) q[3];
u3(1.82266034289203,-2.34501432966756,-0.350710610701489) q[6];
u3(2.06052176268702,-0.547277748707461,0.876260377710139) q[7];
u3(2.22523522928574,-2.11131377385548,-0.988035597716753) q[1];
cx q[1],q[7];
u1(1.12570688429544) q[7];
u3(-3.24144889197279,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.32714885261175,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.60502158793627,-0.203853222079743,0.487493938746692) q[7];
u3(1.05195067600970,0.247944963883185,-4.81751011729913) q[1];
u3(2.38358814142692,-3.37127351413216,2.23817554904741) q[2];
u3(1.46350555851768,-0.112712237494747,2.62808375962364) q[0];
cx q[0],q[2];
u1(0.690457217958035) q[2];
u3(-0.0208655755063303,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.85671001125610,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.860974572439487,-1.08656892247029,-0.385298960192094) q[2];
u3(0.783968032867883,-3.14997834090281,-0.0345791163765654) q[0];
u3(1.45347390790864,1.98477954495245,-0.148143630717514) q[3];
u3(0.472765317140166,1.81734600954941,-4.14717970629827) q[4];
cx q[4],q[3];
u1(1.47386965465699) q[3];
u3(0.140249395696336,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.897291392374976,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.33774219929776,-0.298214394378926,-2.34862840380303) q[3];
u3(1.38549252285067,2.79666004336955,-3.02513054610885) q[4];
u3(2.37055262175991,-0.274483144376322,-2.62777030699375) q[7];
u3(2.73170845056393,4.57051188002281,-0.0383898095936552) q[5];
cx q[5],q[7];
u1(2.12410474526505) q[7];
u3(-1.81252570724327,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.227346861278044,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.33203070387444,-2.98260597995722,-0.0535903035267336) q[7];
u3(1.26885296272881,-2.39678795849483,2.53839540901439) q[5];
u3(1.47661684512923,1.46769907077473,0.427507409458033) q[9];
u3(2.08446296699018,0.955559033792938,-2.53197236780506) q[1];
cx q[1],q[9];
u1(1.92849146120766) q[9];
u3(-2.95602011241397,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.859173073064071,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.73853272423051,2.64548890855466,-2.66886488855919) q[9];
u3(2.87202397006250,-2.95834943826061,0.901634643168280) q[1];
u3(1.11991680328415,3.41299328599871,-1.17252869670151) q[8];
u3(1.11484922725124,2.01687685426957,-1.40201852456999) q[2];
cx q[2],q[8];
u1(1.12661447676146) q[8];
u3(-3.58608818828950,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.77192330004706,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.26312854204238,-0.149808555425970,-1.62312342164093) q[8];
u3(1.90682353890527,-2.92444116720546,-2.91730058425567) q[2];
u3(1.07800908384933,-3.18771379754686,2.74223847401280) q[0];
u3(1.86995355805546,-3.00216002416761,2.30071763633911) q[6];
cx q[6],q[0];
u1(1.53838142503710) q[0];
u3(0.318792989517442,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.01792923993305,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.675233571689476,-1.55123087355930,4.02946456572502) q[0];
u3(0.755618526487830,2.01610780292470,2.82359101017050) q[6];
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
