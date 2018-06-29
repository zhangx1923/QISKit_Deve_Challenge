OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(2.29039479899719,2.89729885128683,-3.14493755100540) q[7];
u3(1.93694777507281,3.01502372779013,-3.00504149136466) q[6];
cx q[6],q[7];
u1(2.61510268887446) q[7];
u3(-2.93249781410574,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.17610495099027,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.18469870644794,-0.127680162005496,0.125027924060199) q[7];
u3(2.45952903025118,0.931696120932002,-0.0809707247803899) q[6];
u3(2.97108451422633,0.123752571908441,2.03112130316975) q[8];
u3(2.60348795416748,-0.805697684220788,0.661770250507201) q[12];
cx q[12],q[8];
u1(-1.30880096417007) q[8];
u3(0.174274182454585,0.0,0.0) q[12];
cx q[8],q[12];
u3(3.50310656471386,0.0,0.0) q[12];
cx q[12],q[8];
u3(2.18448470291005,3.73991535651546,-0.829031592977725) q[8];
u3(0.690972502811650,0.833160333348309,-3.52350960671708) q[12];
u3(1.92810964049836,-1.85602307400224,0.122281864680989) q[2];
u3(1.55665724961835,-4.24767916593753,-1.11876499312606) q[5];
cx q[5],q[2];
u1(2.04755076477633) q[2];
u3(0.0785522264205047,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.978942145981796,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.77521502050786,0.115182144248034,-0.165685858947582) q[2];
u3(2.59894877199190,0.657797030903712,-3.28935736472460) q[5];
u3(1.44849961170769,3.29268989071668,-2.02793435898902) q[0];
u3(1.21720920902646,1.77387331964812,-0.574154453511329) q[9];
cx q[9],q[0];
u1(1.76852450886277) q[0];
u3(-2.77751826603813,0.0,0.0) q[9];
cx q[0],q[9];
u3(0.991166984609606,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.74056704116220,-2.10745563950517,-1.86345757421658) q[0];
u3(2.20199633944387,3.46580382185942,-1.09776705719388) q[9];
u3(0.677782264360596,1.43934784368563,0.590517163422303) q[4];
u3(1.56063908944576,-0.0681933832754635,-3.23381201667229) q[3];
cx q[3],q[4];
u1(1.59526737533007) q[4];
u3(-0.987885733783083,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.64138309984813,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.38248103742890,-0.785980504319881,1.21017984383600) q[4];
u3(1.83084216221797,-0.476297965828363,3.23142670521476) q[3];
u3(0.910229218792582,2.23552520837098,-3.43597368376086) q[11];
u3(1.31654128829389,3.04495364462328,-3.05522724335304) q[10];
cx q[10],q[11];
u1(1.36035344294103) q[11];
u3(-0.0777135369268525,0.0,0.0) q[10];
cx q[11],q[10];
u3(0.696873353865669,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.35347111018908,-1.37814766554805,-1.09246744482905) q[11];
u3(1.85999082467192,2.80277187195327,-1.17095505825607) q[10];
u3(1.07705216086120,-2.62907115953747,3.32604451563164) q[3];
u3(0.571411706346467,0.731810211236365,-0.490390672752818) q[0];
cx q[0],q[3];
u1(1.30721477542181) q[3];
u3(-0.180653200675013,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.21624507292353,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.63252086323701,1.70768650871905,-2.11538343559839) q[3];
u3(0.419478893526102,-4.03728923400310,-1.02991122219019) q[0];
u3(1.56841151596468,-0.673174828636979,-0.123576323266705) q[12];
u3(1.58568537319758,-3.09555095228713,-0.971179804410079) q[9];
cx q[9],q[12];
u1(-1.25663912335095) q[12];
u3(0.691171765413313,0.0,0.0) q[9];
cx q[12],q[9];
u3(3.99341587054799,0.0,0.0) q[9];
cx q[9],q[12];
u3(0.618048490795225,-1.58265607838494,2.49205140383883) q[12];
u3(0.892920986878747,1.97370728028293,-0.182854806558922) q[9];
u3(0.762861588779626,-3.03595427004415,0.0218386357027107) q[10];
u3(1.22866168457113,-2.99737717975834,-0.0231692444416141) q[7];
cx q[7],q[10];
u1(0.550294690041826) q[10];
u3(-3.04747524238560,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.19929536844388,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.62666362639850,2.61725609433743,-0.913647271858380) q[10];
u3(2.86065620924363,0.0614147148925972,3.41535375610330) q[7];
u3(1.37841237976475,-0.169656029666829,1.49683194108128) q[1];
u3(2.08104879452067,-1.14795679869470,-2.99415966960697) q[8];
cx q[8],q[1];
u1(1.08331188144491) q[1];
u3(-1.47101547881060,0.0,0.0) q[8];
cx q[1],q[8];
u3(-0.517804020543571,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.719212758500281,2.19132263837443,-2.34202947868196) q[1];
u3(1.42465857196059,-2.88489459026886,-1.27819199545515) q[8];
u3(1.88962767658788,2.10914763094608,-0.462842614017580) q[4];
u3(2.25409583006639,0.468292223945541,-2.36169887518208) q[6];
cx q[6],q[4];
u1(1.80470937998508) q[4];
u3(-2.53149552516402,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.677944548066795,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.19717301791061,-0.118617337239728,-2.63789277782248) q[4];
u3(1.10655493191679,-4.42632543503618,1.04944416111150) q[6];
u3(2.34514774836802,-2.77722568692767,0.193891287740910) q[5];
u3(1.64196142310074,-2.80767391159307,0.476141675029329) q[11];
cx q[11],q[5];
u1(3.76790307534555) q[5];
u3(-0.984474168176804,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.43728865221787,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.51381550825168,3.34730023298508,-2.53770851081048) q[5];
u3(1.93000320923223,4.50522878503716,1.02396190861263) q[11];
u3(1.98813392634258,-3.61974910513795,2.52609226718371) q[4];
u3(0.790005903770955,2.83051474448125,-1.06205366428989) q[5];
cx q[5],q[4];
u1(0.145954743857873) q[4];
u3(-1.50623136074870,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.47946546300655,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.943291639719466,-1.37590826442022,2.02225130898926) q[4];
u3(1.73374882970667,-2.79352909341259,1.14137178291025) q[5];
u3(2.65323805524504,1.92534201267452,-4.12331979108109) q[3];
u3(1.12714249515843,1.68171752019914,-0.930163966185956) q[2];
cx q[2],q[3];
u1(0.329267528090139) q[3];
u3(-1.45537829734160,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.25676642935835,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.91650693552637,-0.117055872579006,2.46222201566797) q[3];
u3(1.00818099368392,3.72260212558826,-0.476255450662853) q[2];
u3(1.88188653167056,2.17671131601831,-3.74345316372210) q[1];
u3(0.465185918423194,-1.58932587504449,3.07569942358691) q[7];
cx q[7],q[1];
u1(-0.228028561302405) q[1];
u3(-2.30916728897456,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.05809213887093,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.90519139865144,-1.41529444904714,3.49673259375795) q[1];
u3(0.998962219724742,-1.44132860850576,-4.25819054121301) q[7];
u3(2.34003163585012,2.11416483084428,0.00497160836743693) q[8];
u3(1.75963896531475,-0.441465760342482,-1.60821325961599) q[9];
cx q[9],q[8];
u1(3.44601685482964) q[8];
u3(-3.89144238474470,0.0,0.0) q[9];
cx q[8],q[9];
u3(-1.07681369786754,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.91778454552114,-0.833451766267030,-2.14066782944896) q[8];
u3(1.27257763389067,3.97687783641433,-1.63848542097738) q[9];
u3(1.40203450737240,3.03568475216570,-0.410865572446188) q[0];
u3(1.38867904713847,1.67313313982391,-1.87378595064890) q[6];
cx q[6],q[0];
u1(-0.170843396825301) q[0];
u3(-2.05203424778057,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.835006461573068,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.861078219755114,-2.00951932510074,3.21319773456683) q[0];
u3(1.39573935268638,1.91301664532401,2.89570033481498) q[6];
u3(1.84191682764805,-0.274758913679067,0.643750363162957) q[11];
u3(1.90647151965066,-2.42073886374482,-2.29789892391461) q[10];
cx q[10],q[11];
u1(1.80851098639383) q[11];
u3(-3.03143747947328,0.0,0.0) q[10];
cx q[11],q[10];
u3(0.522099454251941,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.70792596178669,1.41222920813844,1.18455370200348) q[11];
u3(1.93639336050118,0.161495992420309,0.427598203701899) q[10];
u3(2.13157551130501,1.09971117738835,-1.71399239263308) q[10];
u3(1.21042579868223,-4.43297583758063,1.53171240046365) q[0];
cx q[0],q[10];
u1(1.19504798419775) q[10];
u3(-3.36710343647863,0.0,0.0) q[0];
cx q[10],q[0];
u3(2.41916734501479,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.898422158301485,3.16819853743137,1.13282468859099) q[10];
u3(1.30547658389600,-0.733773452397376,5.30385707539803) q[0];
u3(2.88057521823533,-0.556983791697599,-0.491493683297714) q[11];
u3(1.38431401293223,-1.00782868477090,-3.58289094585843) q[1];
cx q[1],q[11];
u1(0.934161216896444) q[11];
u3(-0.235533079542776,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.66426796671553,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.01473331895735,-2.16795195398296,-0.702076781987018) q[11];
u3(1.80290069506615,-0.819684609141194,-0.589305294235503) q[1];
u3(2.79956868294523,-1.68117040452814,1.91842527107457) q[5];
u3(2.41762971722709,-2.44642199438192,-0.895961159493786) q[7];
cx q[7],q[5];
u1(2.09233383083786) q[5];
u3(0.524332707238548,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.36343910217046,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.34714949801444,1.84630671733103,0.469112531962202) q[5];
u3(0.471831517399905,-5.21063112490103,-0.499652751360179) q[7];
u3(1.14794652544238,0.142295093993324,0.498271350899289) q[2];
u3(1.39690533654859,-0.952011709592530,-1.49665321772713) q[4];
cx q[4],q[2];
u1(0.531401289025273) q[2];
u3(-1.71076685588743,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.37993718340867,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.33115576086598,4.18368709735710,-1.49906149677106) q[2];
u3(2.17463151012908,-0.952498759233127,3.73687623015641) q[4];
u3(1.13263868471685,-1.48331189161407,0.626364725120880) q[3];
u3(2.27701619844901,-3.90814102712588,-0.449647136377141) q[8];
cx q[8],q[3];
u1(1.65956211538873) q[3];
u3(0.508880975919991,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.04299718520964,0.0,0.0) q[8];
cx q[8],q[3];
u3(2.33182755590574,-1.54078447486878,2.76490492566487) q[3];
u3(1.66203380220342,0.320549009438738,3.24937405712676) q[8];
u3(1.86189966676127,-2.95742376798182,0.725554525102506) q[12];
u3(2.80520779868616,-3.03471096950924,-1.78118993343155) q[9];
cx q[9],q[12];
u1(3.84372155421899) q[12];
u3(-3.58203189159147,0.0,0.0) q[9];
cx q[12],q[9];
u3(-1.10634464864027,0.0,0.0) q[9];
cx q[9],q[12];
u3(2.65177124909173,-2.05221311932677,0.565744434995515) q[12];
u3(2.48353234759487,4.99442782663266,-1.08021219169039) q[9];
u3(2.77807078801643,4.21064979650184,-1.69754756217328) q[3];
u3(1.11483446883544,-0.556354096144600,2.18307651786592) q[12];
cx q[12],q[3];
u1(-0.179210935281571) q[3];
u3(-2.17617663918102,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.30967603627567,0.0,0.0) q[12];
cx q[12],q[3];
u3(1.42681790450560,1.96370513633201,-1.06124393959913) q[3];
u3(1.49497088838789,-0.484016583251331,1.86259423824551) q[12];
u3(1.18230070123489,-0.333295743472049,2.39278211588869) q[11];
u3(1.12576179932377,-1.57778944329540,-2.28688956409658) q[5];
cx q[5],q[11];
u1(0.184253244426045) q[11];
u3(-2.58387798518847,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.80144161608011,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.13207209726677,-1.71840984542800,3.20941332166141) q[11];
u3(2.20306427118056,-2.52325527631932,1.29296859174533) q[5];
u3(1.03948685579693,0.310229796624286,2.81634412717105) q[4];
u3(0.697486466615315,2.94448086410051,3.04256772678423) q[10];
cx q[10],q[4];
u1(3.47099645147206) q[4];
u3(-0.845999555808750,0.0,0.0) q[10];
cx q[4],q[10];
u3(1.75134037857615,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.838588817480560,-3.95357690042665,1.24243318143246) q[4];
u3(2.21254205583512,-2.80565006960131,-2.66451493034544) q[10];
u3(1.57759145593910,0.950686752840386,-2.98521738670868) q[1];
u3(0.472161523709107,2.32450803493700,-2.74554097738600) q[6];
cx q[6],q[1];
u1(1.79366240384948) q[1];
u3(-3.17892164529559,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.72486936876908,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.82288352525121,3.52980933052099,-0.263002032239333) q[1];
u3(1.68125275218646,-0.203610123741267,-2.80242113783230) q[6];
u3(2.81407111786012,1.37905469590092,0.310467076544851) q[2];
u3(1.42835645148093,-0.237689119868717,-2.65231474536523) q[8];
cx q[8],q[2];
u1(4.21769596810030) q[2];
u3(-3.70750637956450,0.0,0.0) q[8];
cx q[2],q[8];
u3(-0.102370598339773,0.0,0.0) q[8];
cx q[8],q[2];
u3(0.376025931801594,0.873864009951606,-1.87288043339134) q[2];
u3(1.87676361110418,4.45743088197653,-0.753364067118867) q[8];
u3(0.804831238253662,2.87947523797011,-0.675119699396467) q[9];
u3(1.56899995004584,0.220779418844538,-1.20737393021638) q[0];
cx q[0],q[9];
u1(1.97036750040937) q[9];
u3(0.534134182538685,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.01778530554584,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.47303122680683,0.979486760278188,-3.68992427528257) q[9];
u3(1.43140326896920,1.98078684852537,1.76172048886920) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12];
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
measure q[12] -> c[12];
