OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.07599524715815,0.569022244675786,2.42280961336456) q[4];
u3(2.12266399494427,-2.96109879103142,-2.48226585454941) q[3];
cx q[3],q[4];
u1(3.41751327892523) q[4];
u3(-1.52940897791473,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.12112671762794,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.66999116827645,-0.531479767605430,1.25886979273269) q[4];
u3(1.71010997575797,1.77631244971511,3.11074362947931) q[3];
u3(1.70967973391648,1.04810403442339,1.40066234618811) q[8];
u3(1.73018335978065,-1.21999685823103,-1.69075361604664) q[2];
cx q[2],q[8];
u1(2.88385882732189) q[8];
u3(-2.19183056348288,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.55519692255273,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.48054669336801,-1.03527539049120,2.01738087353393) q[8];
u3(2.05309464455732,-1.54552590556137,-4.63920681132916) q[2];
u3(2.61132374774695,1.56306407901108,0.895831840459716) q[6];
u3(0.901999496609823,-0.151418400517487,-3.78596123489321) q[1];
cx q[1],q[6];
u1(2.22412581767955) q[6];
u3(-1.82938870219836,0.0,0.0) q[1];
cx q[6],q[1];
u3(-0.116229322026550,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.51460104230248,0.787476464324017,-2.76652030687236) q[6];
u3(1.57917284231814,-2.18453537355978,2.97566568290246) q[1];
u3(1.15457244372591,0.900720951927812,-1.07504120655139) q[5];
u3(1.14443200289692,-0.768401465666202,-0.589758950726734) q[9];
cx q[9],q[5];
u1(1.58931780723574) q[5];
u3(-3.46530460839840,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.95020951096007,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.50119385196233,1.79790238172115,-0.694127231225175) q[5];
u3(1.38352895883742,-0.914013958023677,-3.44974821435656) q[9];
u3(1.07709503552242,0.868421497312077,-1.35153923625039) q[7];
u3(0.595454792953708,-1.73351849153994,0.329242629108020) q[0];
cx q[0],q[7];
u1(0.549423688368725) q[7];
u3(-1.02786395636866,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.71253186418751,0.0,0.0) q[0];
cx q[0],q[7];
u3(2.61332888181251,-2.18060757596088,3.80735122066261) q[7];
u3(1.55541968667352,-0.690380078488385,1.99438474616570) q[0];
u3(1.29766948490745,1.75183517346154,-2.91228406650851) q[1];
u3(1.75906301115370,2.06033104581317,-4.13213966887407) q[9];
cx q[9],q[1];
u1(-0.347141859209217) q[1];
u3(1.15346523706842,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.58395786907885,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.59133381003176,-0.245692553262957,2.00853278273603) q[1];
u3(2.32952157415534,1.17084849577536,-4.73550536844021) q[9];
u3(2.22812329320831,-0.739359921400553,1.44728944364889) q[8];
u3(1.72843742905190,-1.73865943923394,-0.197316993816170) q[7];
cx q[7],q[8];
u1(2.67107629711458) q[8];
u3(-2.78471305212418,0.0,0.0) q[7];
cx q[8],q[7];
u3(-1.20960734457677,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.26352313216371,-3.91776928105658,2.19474286636336) q[8];
u3(0.919457268959647,1.69985162212790,-3.21836655915655) q[7];
u3(1.56954866391010,1.65258410238405,-3.18499087400678) q[0];
u3(1.90570662087879,2.05500274903852,-3.84449713381353) q[2];
cx q[2],q[0];
u1(0.354897616304017) q[0];
u3(-1.62070302172860,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.64312427455938,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.807453148350441,0.487836844849841,-3.34252761726508) q[0];
u3(1.74158536165463,-1.83742480070467,4.33128133762121) q[2];
u3(1.71631077321240,0.951428998013814,-3.25294057869632) q[5];
u3(1.23487607603245,3.27576964321245,-2.95675744379138) q[4];
cx q[4],q[5];
u1(3.00104619195868) q[5];
u3(-2.03745447408027,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.540921868977745,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.22043381673778,1.00278169787140,-2.87984122600592) q[5];
u3(2.05074757452152,1.63730811282992,-0.0238834769902226) q[4];
u3(0.697446414029001,0.753015830445418,0.846457165326495) q[3];
u3(0.883825499886591,-1.14361492037977,-3.26193315377890) q[6];
cx q[6],q[3];
u1(1.53984346070764) q[3];
u3(-2.59230391709150,0.0,0.0) q[6];
cx q[3],q[6];
u3(3.08136424519945,0.0,0.0) q[6];
cx q[6],q[3];
u3(1.67073644362859,-3.66883755616183,0.430775634799800) q[3];
u3(1.17214029409134,1.94489027468076,1.45987099818379) q[6];
u3(1.01277423995545,2.46798110093917,-3.36476469557290) q[1];
u3(1.13636505869473,2.90865219791023,-3.19993302443897) q[8];
cx q[8],q[1];
u1(2.18651194634874) q[1];
u3(-0.202417050311976,0.0,0.0) q[8];
cx q[1],q[8];
u3(1.15027376306478,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.54641670323824,-1.22118688941001,3.08887273558752) q[1];
u3(2.10034452370873,-1.42930664292743,-3.16650967150344) q[8];
u3(1.54915645018079,3.43377816039138,-1.25386186710957) q[4];
u3(1.97600427926551,2.95088316311795,0.186769394872400) q[9];
cx q[9],q[4];
u1(-0.471338070094029) q[4];
u3(-1.63729171326767,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.608510331222562,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.51416350889767,-2.45504552362494,2.44391075620975) q[4];
u3(0.300866336133951,1.93302035771044,-2.71072289811027) q[9];
u3(1.45556447654209,-1.91014111444382,1.65177344498735) q[0];
u3(0.380176004279819,-1.63002038926807,-0.183736854101246) q[3];
cx q[3],q[0];
u1(2.17345773133085) q[0];
u3(0.313292252117026,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.65659438728745,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.890550576363826,0.686644418889757,1.71835694491657) q[0];
u3(1.18399937976380,0.386493550146190,2.69409105582406) q[3];
u3(0.544950892311031,-1.07791707381799,1.13714931653210) q[6];
u3(1.03602783119515,0.0887215817157012,-1.11153296759147) q[7];
cx q[7],q[6];
u1(0.0183734257149941) q[6];
u3(-1.50155778195881,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.771088653626367,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.99604190020355,1.70691366823452,-1.32427457095150) q[6];
u3(1.29487800819458,-5.22039997689051,-1.05000432553096) q[7];
u3(1.70527912315453,-2.43028759349855,-0.306735393089768) q[2];
u3(1.85444727187601,-3.85061424568332,-1.04108367603396) q[5];
cx q[5],q[2];
u1(1.27845060376093) q[2];
u3(-0.0919387289590823,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.53539842344222,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.78184707603805,2.80952672928189,-2.87233990494261) q[2];
u3(1.39252468933519,2.34975905537230,-1.15734213177289) q[5];
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
