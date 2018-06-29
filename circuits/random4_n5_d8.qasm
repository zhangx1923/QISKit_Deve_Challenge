OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.957383423370394,1.24742818212308,-2.29858414776707) q[4];
u3(1.23587283022533,-2.30367708503858,2.59796154986578) q[2];
cx q[2],q[4];
u1(-0.224771634971024) q[4];
u3(-2.32105864549501,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.87823170360158,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.33642999582700,1.55703440446957,3.03098563505239) q[4];
u3(2.15413494931528,0.784298062281161,-5.40455440634812) q[2];
u3(2.39719680117826,0.160695672654731,-0.333157039886193) q[3];
u3(0.927681614885770,-3.62853620376339,-1.06905321670990) q[1];
cx q[1],q[3];
u1(3.29935791118663) q[3];
u3(-1.14129747217118,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.35455315374653,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.38493982873285,-2.43486062400750,3.24221111065854) q[3];
u3(0.762041599477892,-0.0405367851202911,-0.0898541921252697) q[1];
u3(1.53579971588666,1.37154940689545,-1.82584485913832) q[2];
u3(0.396999679808217,1.14661536701775,-4.14590035972429) q[4];
cx q[4],q[2];
u1(0.872178555812773) q[2];
u3(-1.23483878995645,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.64729231611147,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.25060217217113,-1.97243452636298,3.99019693177817) q[2];
u3(1.19308773506320,1.43931573992114,-1.65849263004591) q[4];
u3(2.11391010460089,1.29605652184541,0.423696666588790) q[3];
u3(2.03044714090425,0.751265682018651,-2.21378283537227) q[0];
cx q[0],q[3];
u1(1.05445218877071) q[3];
u3(-0.209228287986900,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.31954522910848,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.32809474844541,-2.39702066583311,2.33519530602752) q[3];
u3(2.72110079988111,-4.49786886105605,1.68967235834257) q[0];
u3(2.30869227954098,0.642373909228963,1.43267995679113) q[2];
u3(1.58373417850936,-2.14335357845958,-2.06915267353234) q[4];
cx q[4],q[2];
u1(2.15858911350521) q[2];
u3(-2.98339007068725,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.29577486364081,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.66869268282703,3.32020923658711,-2.24136348820541) q[2];
u3(2.46780377054234,-1.59994282601507,1.92155073648141) q[4];
u3(0.560119363484825,-2.74442436389727,0.725270975387562) q[1];
u3(1.70357798490206,-2.70601039633265,1.30531044521275) q[0];
cx q[0],q[1];
u1(2.91201539177967) q[1];
u3(-2.00016390948566,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.427290690060391,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.95430014434398,0.964121599866670,-0.549670119322005) q[1];
u3(1.63297338461323,-0.419800638391461,5.44849604711739) q[0];
u3(2.18296382246479,3.17163352213926,-0.876402871279948) q[2];
u3(2.40037188401284,1.05397345300198,-2.04226365828794) q[4];
cx q[4],q[2];
u1(3.24922199141140) q[2];
u3(-0.672924940917755,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.00345789434667,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.69593498057857,-1.71545096957732,4.07496085654443) q[2];
u3(0.993236754955694,-2.61543526469441,-2.73039652210199) q[4];
u3(1.51038755581794,0.707619912745745,2.04408941426796) q[0];
u3(0.832352644581321,2.86305988614155,3.27554261361352) q[1];
cx q[1],q[0];
u1(3.70068749779476) q[0];
u3(-1.42256563469844,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.00872925864159,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.55125701104485,-0.802260324618039,-0.368600087945814) q[0];
u3(2.37119082634499,-4.86853025238057,0.179864226658735) q[1];
u3(2.30694528434519,2.29586937312530,-0.231661219298632) q[2];
u3(2.60264733493348,2.09196866218056,-1.80482941757693) q[4];
cx q[4],q[2];
u1(2.85007913425697) q[2];
u3(-1.79297675063585,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.09811122152163,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.15413079502879,1.36169387989830,-1.04007683690668) q[2];
u3(1.88410005745679,-0.599363467896433,-2.42438723724816) q[4];
u3(0.210441215716474,2.23293122309336,-2.64010282373951) q[0];
u3(0.864073707433333,-3.77362340427119,2.21734159096633) q[1];
cx q[1],q[0];
u1(-0.223979846141799) q[0];
u3(-1.54437058880452,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.709022416542284,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.56777944683481,-1.25483462760159,-0.742092822664724) q[0];
u3(1.85259800404567,-3.31286104416571,-2.90310101004776) q[1];
u3(1.34986223572777,0.703743350813834,1.34804668239940) q[2];
u3(0.825839042922763,-1.36362318760034,-3.09757988845196) q[1];
cx q[1],q[2];
u1(1.36096078056867) q[2];
u3(-1.20322299758183,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.268687229673530,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.14376545377430,4.03236843187259,-2.04171595856991) q[2];
u3(0.924211706427835,-2.47910812101865,3.07880280123113) q[1];
u3(2.35521965400982,1.56779008640505,0.609768816921621) q[3];
u3(1.72297532910095,0.00780989275975763,-4.08267113621712) q[4];
cx q[4],q[3];
u1(2.88603823779833) q[3];
u3(-2.35346474983646,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.66110212800933,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.60415664040700,-1.16092949974766,3.57573988496423) q[3];
u3(0.419636471532124,1.54855718379287,-2.19503397676552) q[4];
u3(2.26634982094960,2.18005208489833,0.463507169647365) q[1];
u3(2.23356373913970,0.334256028177244,-3.67989485104738) q[3];
cx q[3],q[1];
u1(3.11673855591077) q[1];
u3(-1.43556894725034,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.24110943453537,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.67557325760167,1.78174826949815,0.452275700769713) q[1];
u3(1.74156049560667,-0.0301414743051285,-4.05213476906763) q[3];
u3(1.57464622077422,-0.170614381655677,1.70293794956237) q[4];
u3(1.45049075972991,-2.63707125406334,-2.23753796264375) q[2];
cx q[2],q[4];
u1(1.86239358754079) q[4];
u3(-2.52491741733362,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.676152347373339,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.779746308650423,-0.740536606455654,3.86354762647276) q[4];
u3(2.32219322350808,0.0374473732730254,1.99359773270694) q[2];
u3(2.28291391616048,2.74018821013845,-3.12720028497237) q[2];
u3(2.47945594723732,-3.64354062401091,2.60100705718395) q[0];
cx q[0],q[2];
u1(2.31247920566237) q[2];
u3(-1.40544658203503,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.235096056435044,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.43334026682213,3.24366855389476,-1.16923427173312) q[2];
u3(0.968824558265063,-1.33783117435898,-0.820198902857993) q[0];
u3(0.488755350686109,2.67212839427006,-1.76647198779615) q[4];
u3(0.657227651885089,2.10336950163753,-3.45035592919944) q[3];
cx q[3],q[4];
u1(2.43436854319966) q[4];
u3(-1.73259220719781,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.641403162502685,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.418283291520452,1.42832474355893,-0.660410470728865) q[4];
u3(2.49176073696555,1.85308242108852,3.48219381269249) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
