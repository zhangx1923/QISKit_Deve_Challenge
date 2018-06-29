OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.88903446962518,1.23791728062099,1.72992309111820) q[0];
u3(2.00757788357738,-1.46091027332236,-1.59153202815453) q[3];
cx q[3],q[0];
u1(3.32590874636558) q[0];
u3(-0.828378505103096,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.47058678410674,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.711990270150924,0.0235563150761502,2.26198547484973) q[0];
u3(1.71481436642132,-2.87397053494541,2.40846505082127) q[3];
u3(1.17675880747571,1.28227423122123,-0.854204508940435) q[1];
u3(0.131142455330188,-2.79571937021345,0.158013159693245) q[2];
cx q[2],q[1];
u1(2.48528675341564) q[1];
u3(-2.00829771042881,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.0695180627700394,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.05035864171220,4.08500906201681,-1.96851120564114) q[1];
u3(2.27298219108455,1.43142593626453,3.30558180251794) q[2];
u3(1.96988088141723,2.29318181397286,-0.667134517535310) q[1];
u3(2.69399588622230,0.365548635234686,-2.86483179807366) q[3];
cx q[3],q[1];
u1(0.963707096880691) q[1];
u3(-1.36953929519862,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.0786868186538188,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.725245295854124,0.229336214518697,-2.52704001731144) q[1];
u3(2.89001938023716,4.14960381813138,-2.00644264693418) q[3];
u3(2.03804094342940,0.700306202196232,-3.00210824884338) q[0];
u3(2.05208682311290,-2.62172857744105,3.35492339893819) q[2];
cx q[2],q[0];
u1(1.20671901508081) q[0];
u3(-2.63974306912039,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.92098574348788,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.999775964499651,3.12852126604317,-2.84291250409902) q[0];
u3(0.843372048601343,-2.74350195678096,-3.11240817085404) q[2];
u3(1.61670581853583,3.52406780047144,-1.92148970597680) q[4];
u3(1.99278119870094,1.24851279185769,-2.78066919201617) q[0];
cx q[0],q[4];
u1(2.43964208348032) q[4];
u3(-1.72139783622627,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.33102758691121,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.66478170396706,0.366429552978752,-2.24705370280443) q[4];
u3(2.26008680882450,-0.347880701182644,-5.53448559246816) q[0];
u3(1.96043891833759,-2.16439591513396,0.852207935851433) q[2];
u3(1.76678938055263,-2.77138783826173,-0.344399903178754) q[3];
cx q[3],q[2];
u1(3.16538603352130) q[2];
u3(-4.33184040465537,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.297281936350945,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.03204774921638,0.637105702069333,2.76611328281403) q[2];
u3(1.47961492871401,-2.34723476677113,0.889621480573775) q[3];
u3(1.33210400961329,0.669204637685118,-1.41741516233662) q[3];
u3(2.36180973091001,0.880172519726024,-5.25004229805814) q[0];
cx q[0],q[3];
u1(4.02713988132358) q[3];
u3(-4.38234235195303,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.627448431376206,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.82891147537994,-0.192718838988171,-3.88564237266369) q[3];
u3(0.766954282073125,2.72141505455828,-3.14030073784686) q[0];
u3(0.510539442031661,-0.427096407604579,-0.152385396827718) q[4];
u3(0.952068641787415,-2.39293654822826,1.47913104193483) q[2];
cx q[2],q[4];
u1(0.00156292748068343) q[4];
u3(-2.57067657874337,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.41647866093088,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.12069321347775,4.22917968366825,-1.64954080035255) q[4];
u3(0.752716150112908,-3.79801305041574,0.886916370143666) q[2];
u3(1.94857353643828,0.880129924704622,0.597617168695365) q[2];
u3(1.93932101163184,0.0371340021663980,-3.45257034598426) q[1];
cx q[1],q[2];
u1(1.66040531567153) q[2];
u3(-2.31174761144942,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.58003771836485,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.677171769339104,-3.89076261955642,2.22246520504694) q[2];
u3(1.13537819548117,-0.130351981665885,-4.66637321532370) q[1];
u3(1.11222472864179,0.821370574444263,0.140792725501785) q[0];
u3(1.86496130806495,0.0886891602237287,-1.87063064021967) q[3];
cx q[3],q[0];
u1(-0.520050676179973) q[0];
u3(-1.79640801540600,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.05051983208728,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.92982823663353,-0.906522385395762,-0.709893267782343) q[0];
u3(2.06131958643887,1.02590901951301,0.907974010945191) q[3];
u3(2.83058458502398,2.08943987639213,-1.04376272683251) q[1];
u3(1.77110032740308,5.18957489721359,-0.440930542931479) q[3];
cx q[3],q[1];
u1(2.80830517362384) q[1];
u3(-2.47598492848252,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.21648008883860,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.29436077597352,-1.13178502964096,0.130222905006356) q[1];
u3(0.255332751400917,2.51107679499094,-2.86107478167419) q[3];
u3(1.69779407565064,1.09701253081897,1.54775757080619) q[4];
u3(1.75377163424910,-1.77536150697009,-2.29628066225430) q[0];
cx q[0],q[4];
u1(0.427187749637035) q[4];
u3(-1.46395359756083,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.15456160192598,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.432243542118723,1.61271045296625,-1.46525799100350) q[4];
u3(1.42612830180015,-1.68831934780038,4.13259970703466) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
