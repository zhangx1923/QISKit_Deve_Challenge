OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.62376861605504,-1.26555280705831,0.373664473863295) q[3];
u3(1.80515668005681,-2.55448645543398,0.690935554454310) q[2];
cx q[2],q[3];
u1(3.58868686296399) q[3];
u3(-3.43583224336970,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.429321541222350,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.31673301415959,-3.45832202919015,-0.950376745015578) q[3];
u3(2.11418959614606,1.38549275563848,4.39504008091512) q[2];
u3(1.97036242259715,-0.863385927653271,-0.220529275543325) q[1];
u3(2.07649591646383,-3.07805934587377,-0.295600923524207) q[0];
cx q[0],q[1];
u1(3.36449458592037) q[1];
u3(-1.69511247856594,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.26832751556980,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.50156660810227,3.27451479197207,0.302401590844478) q[1];
u3(1.29580517682732,-3.44589444634438,2.30079714768078) q[0];
u3(2.09906742401156,0.290240512754029,1.65655327789789) q[4];
u3(2.13358062411455,-0.792275004292168,-1.78130279281289) q[6];
cx q[6],q[4];
u1(3.28215883955665) q[4];
u3(-0.464042048260818,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.27240860504755,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.876659644999930,-2.11274068196752,3.62058509836450) q[4];
u3(0.430132823515186,-0.263775642622530,1.33429995139378) q[6];
u3(1.59850904600142,3.75563673210273,-1.58819193068274) q[5];
u3(2.75511443816313,1.31041894926051,-2.36738883954224) q[7];
cx q[7],q[5];
u1(3.14514713680862) q[5];
u3(-2.07681423897015,0.0,0.0) q[7];
cx q[5],q[7];
u3(0.305246116154303,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.88181231014722,-2.10182215447743,0.0983340066785945) q[5];
u3(1.06371265231857,-1.52171701391194,-3.57987740417431) q[7];
u3(2.44708210119827,-1.60534035881291,0.355595658077748) q[5];
u3(1.89977732053861,-1.33144979347503,0.720099378860970) q[0];
cx q[0],q[5];
u1(1.88970395491250) q[5];
u3(-2.54264698005676,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.324492470946619,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.29758128862967,-0.731120466039429,1.66791148999529) q[5];
u3(1.99997513450337,-3.54109825577792,1.07225173715468) q[0];
u3(2.60844057533536,-1.61721396569879,2.69569018066300) q[3];
u3(2.54313068867441,-3.67128318378397,-2.01482402482604) q[6];
cx q[6],q[3];
u1(0.588253324339345) q[3];
u3(-1.75242713876248,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.03758442443445,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.49518088668914,-2.11462973152132,2.66496185197152) q[3];
u3(1.10940279741954,4.23844760998749,-1.56384228409727) q[6];
u3(1.41454457823259,-0.0414731214728354,1.76480535856236) q[2];
u3(1.92243146040263,-0.486512034808074,-2.29329322608210) q[1];
cx q[1],q[2];
u1(2.22144842999155) q[2];
u3(-1.49477659177781,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.40246803698337,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.79794244961609,1.43644335632898,-1.14951126123352) q[2];
u3(2.21534252205790,-0.853057941041083,-3.43433391615232) q[1];
u3(1.64340376475244,-0.778393019869073,-0.253146680180301) q[4];
u3(1.51798676200794,-2.72778587517351,-0.303659877215476) q[7];
cx q[7],q[4];
u1(2.20652673790505) q[4];
u3(-1.70375596130037,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.337347859670627,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.750872929164361,-4.29910286180786,1.45663674774869) q[4];
u3(2.64948806165357,-2.04620991302214,-1.09422079458571) q[7];
u3(0.993021770626318,-0.576182222296366,0.834287155616600) q[0];
u3(1.61264888693464,-3.16194920816068,-0.0232250884925715) q[1];
cx q[1],q[0];
u1(1.14448822249882) q[0];
u3(-1.53429571113299,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.742763110593522,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.18827771521254,-0.809426984996348,-2.03390195680353) q[0];
u3(2.42437503936321,-0.947652241991299,-1.68479376204197) q[1];
u3(2.16273647890685,2.07857106389198,-3.75850342080990) q[4];
u3(0.356626181934182,3.37153349377513,-2.48185654589964) q[3];
cx q[3],q[4];
u1(1.71370222767056) q[4];
u3(-0.541515163208870,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.02362989052114,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.853981759317293,2.58134999274554,-1.73402396635013) q[4];
u3(2.80019538747000,-0.494045678668883,-5.69102237035097) q[3];
u3(1.21506019802718,1.17020805788039,-1.37213066184471) q[2];
u3(0.353002788466947,-0.0419626453099869,-0.144782445679379) q[5];
cx q[5],q[2];
u1(3.59208667091971) q[2];
u3(-1.50247566933549,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.39516103900169,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.25260000336898,0.558284747613176,0.220382192949311) q[2];
u3(1.45285960411403,-1.13686467848581,3.91500694694850) q[5];
u3(0.889672361107234,1.06174447326224,0.341068028616670) q[6];
u3(0.843213322688471,0.223246513445007,-2.79789752994897) q[7];
cx q[7],q[6];
u1(1.58682761549731) q[6];
u3(-2.93793665160382,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.480964581725802,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.09341212677091,-1.94025051479204,1.24786366650410) q[6];
u3(0.753874461095339,3.32428067939519,2.02151110921418) q[7];
u3(2.87703222334127,-0.920216583144930,-1.39300262711178) q[7];
u3(0.938335741713063,-5.08771446112207,1.15212803389373) q[5];
cx q[5],q[7];
u1(2.89982528073720) q[7];
u3(-1.92399036934604,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.57620694195179,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.29070926085867,-0.670679185268778,2.62697634245321) q[7];
u3(1.55552687857737,-0.345697505305780,2.19235672675093) q[5];
u3(2.73248459363536,-0.778055340697775,1.80372114356525) q[1];
u3(3.01493490175251,-1.65210059885128,-0.509274953460357) q[4];
cx q[4],q[1];
u1(2.52083958002682) q[1];
u3(-1.47375225697096,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.546266900095530,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.24238845430392,-1.52603691629818,-0.483884118556447) q[1];
u3(0.620759175632697,-3.25439989680775,1.45515044309418) q[4];
u3(0.460773741562565,3.36857364866282,-2.79100532176857) q[0];
u3(0.814693608837869,-3.97679677545791,1.36612841629032) q[3];
cx q[3],q[0];
u1(-0.147648115871458) q[0];
u3(1.07662099497766,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.57309149074212,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.410482635957917,-1.57879973843596,0.855323823079969) q[0];
u3(0.756840744987027,-0.870071597969261,2.03216663189889) q[3];
u3(2.13543575198384,-0.986014737717682,0.616378235381840) q[6];
u3(2.29704749692461,-2.41499412891138,-1.19386943024416) q[2];
cx q[2],q[6];
u1(2.28867027008138) q[6];
u3(-0.156312592248927,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.26849934425900,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.95958125084086,2.11071969185898,-1.14145428445253) q[6];
u3(0.317259298906300,-0.606885141958558,-5.11023647768120) q[2];
u3(1.33138286210585,-1.84259558267292,0.690556529582676) q[1];
u3(2.39764003904771,-3.55418856021452,0.191304214185717) q[5];
cx q[5],q[1];
u1(1.49611386533719) q[1];
u3(-0.163524335178442,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.26836194595855,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.86577232024546,-2.90807434639595,2.73526087626781) q[1];
u3(0.470912637993978,1.39714363481291,-1.47096962667295) q[5];
u3(1.42608525596997,0.699152586625541,-2.10811224822255) q[2];
u3(2.41052404984294,1.74486351527098,-4.07601808371782) q[0];
cx q[0],q[2];
u1(1.26492220295366) q[2];
u3(-0.347500274115008,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.83266025889608,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.81066354769081,2.49888166398013,0.215803869140776) q[2];
u3(1.84541344086416,-3.67688129859495,2.43169656218681) q[0];
u3(1.79631369552397,0.861961991252923,-2.56003546099200) q[6];
u3(2.46614384487707,1.47720138464906,-4.54010885516480) q[7];
cx q[7],q[6];
u1(0.874543187061323) q[6];
u3(-3.59903060663832,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.77023136891516,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.76834444300176,-3.49781848631667,2.13812157494817) q[6];
u3(1.44468314747575,-2.25329886885899,1.49361820813386) q[7];
u3(2.62509965648705,2.36854515331932,-3.23983257230342) q[4];
u3(1.95136320199742,-2.78555626610282,2.96442346532128) q[3];
cx q[3],q[4];
u1(1.91764245656144) q[4];
u3(-3.25747600870302,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.512419314431791,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.80461069280880,0.514952316963503,-2.46450911751032) q[4];
u3(1.05250311290921,0.341546588958060,3.11659921889023) q[3];
u3(1.84922989359096,-1.46351734127744,1.89737415453804) q[3];
u3(1.90715868371658,-1.76863301510570,-2.33203694591551) q[4];
cx q[4],q[3];
u1(2.43567910491176) q[3];
u3(0.243042276403591,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.51107712694001,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.35671564848760,-2.92195299468914,2.64739153779991) q[3];
u3(1.51668615936107,5.02319135132212,0.191360901319142) q[4];
u3(2.55419745537837,2.38155401490012,-0.353618479474088) q[0];
u3(2.28350428020956,-0.337896078109622,-3.67647620092479) q[5];
cx q[5],q[0];
u1(1.61313181679353) q[0];
u3(0.639460217485074,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.00654094419408,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.48121231829874,-3.84935346262218,1.22864271642712) q[0];
u3(1.25116068025940,1.40074191553567,3.94975693923934) q[5];
u3(1.18685693851112,-0.736703192058883,0.909732863063676) q[7];
u3(2.11215918209860,-3.48300829399545,0.0492286367169965) q[2];
cx q[2],q[7];
u1(3.86100246438887) q[7];
u3(-3.94153893633676,0.0,0.0) q[2];
cx q[7],q[2];
u3(-0.924676908714992,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.25125025408386,1.63694177199108,-0.223945706024982) q[7];
u3(1.54409764454898,-4.13002607544714,1.75118154691516) q[2];
u3(2.05434348983019,-1.65066846872644,1.64707016789893) q[1];
u3(1.96299755176585,1.53471051153209,3.54528811745283) q[6];
cx q[6],q[1];
u1(-0.416510519013319) q[1];
u3(-1.56173643567389,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.903282996275126,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.13171727356537,1.85412577778597,0.137059855534433) q[1];
u3(1.45055455324595,4.02184170035358,0.871301855469127) q[6];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
