OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(1.53668043080240,0.305584560183670,0.566693423867660) q[11];
u3(1.29194866638233,-1.06719566597671,-1.54569686062901) q[4];
cx q[4],q[11];
u1(1.68793509808363) q[11];
u3(-2.71398677402807,0.0,0.0) q[4];
cx q[11],q[4];
u3(0.753829853764561,0.0,0.0) q[4];
cx q[4],q[11];
u3(2.07208158652559,0.0943640085044912,2.42105470928499) q[11];
u3(2.15176057140123,3.74557817269682,0.137528568337538) q[4];
u3(1.96459181510113,1.46193400549498,-3.19740118995917) q[3];
u3(1.69035288413036,-1.79425112676679,2.77402663042795) q[12];
cx q[12],q[3];
u1(1.78017444410222) q[3];
u3(-0.00673954832222101,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.00122178543803,0.0,0.0) q[12];
cx q[12],q[3];
u3(1.41906362374211,-1.66695412558301,0.852411255511090) q[3];
u3(1.41572534059538,-0.00719152745981466,-5.05385642667149) q[12];
u3(0.845188642669824,1.67934422639486,-3.62185393079197) q[0];
u3(1.74834723891530,-2.85366683111125,3.10906034798966) q[5];
cx q[5],q[0];
u1(1.75626849906796) q[0];
u3(-2.45756812307763,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.62624492358969,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.81358399017720,0.210878882307023,-2.96685869549098) q[0];
u3(0.806157229645170,-5.05320159950299,-0.832392647990563) q[5];
u3(1.39900421052445,-2.62350669038735,0.158591595442386) q[2];
u3(2.06911328943136,-2.80955374331183,-0.335980438603439) q[1];
cx q[1],q[2];
u1(0.182371465635691) q[2];
u3(-0.986904787781363,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.60437331008528,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.15151903463950,3.19585840627456,-0.725795212361402) q[2];
u3(2.34301474595397,-1.25263247544510,-3.70911590383193) q[1];
u3(1.57621583386413,2.51572411172709,-0.218108493698885) q[7];
u3(2.72526221567572,0.484874140743598,-2.58826433794277) q[9];
cx q[9],q[7];
u1(1.63807777187978) q[7];
u3(-3.37415546673035,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.32875402764130,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.41600723305422,-1.13678337018419,1.44464525943872) q[7];
u3(2.29229344071988,0.319821433978624,-0.696736441718465) q[9];
u3(1.81529551157815,1.60138553421172,-3.37520057663363) q[6];
u3(2.25640070523647,-1.53714044367888,4.08512606807846) q[10];
cx q[10],q[6];
u1(0.171399653587341) q[6];
u3(-1.02989455255954,0.0,0.0) q[10];
cx q[6],q[10];
u3(1.50780701409149,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.822144976792549,-2.99465615163210,2.78692851099409) q[6];
u3(0.672547109855882,-2.19420070048861,-3.54810096941494) q[10];
u3(1.32120329646305,1.87473033419756,-0.115712708412715) q[6];
u3(1.37958228452317,0.356725225029886,-3.92884025298330) q[11];
cx q[11],q[6];
u1(2.92458347824003) q[6];
u3(-1.97512626479297,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.60840225795186,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.886460936994896,1.28777592675071,-3.40478995319134) q[6];
u3(2.54267010822857,1.68128584605523,-1.34092181881154) q[11];
u3(0.927629615753651,-1.18714222829186,0.499274907221987) q[2];
u3(0.781680491574521,-1.33500370881337,-1.14359826493996) q[4];
cx q[4],q[2];
u1(3.56632940829598) q[2];
u3(-1.19473677386062,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.33539153021578,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.29644759200525,-1.80395831530786,-2.05984522186525) q[2];
u3(0.550056495840419,-1.27102334514349,2.22268572372343) q[4];
u3(0.413330040577712,2.04814359230547,-2.77036273303948) q[9];
u3(0.564894872459925,0.426712714937351,-1.42433794502327) q[10];
cx q[10],q[9];
u1(0.0822667130921515) q[9];
u3(-2.31362172351004,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.20828715231250,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.75065971232874,1.39801882561834,-3.92604620211943) q[9];
u3(1.79165863374444,0.291232537601217,-0.913784314969192) q[10];
u3(1.42700043908611,-1.81556656291994,-0.100796656719892) q[0];
u3(1.77603560252508,-3.77254734996782,0.757721465650955) q[3];
cx q[3],q[0];
u1(2.78715394398425) q[0];
u3(-2.95636198841295,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.662127987247068,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.72315116885576,-0.214184051492376,2.63323124829560) q[0];
u3(3.01498051685865,-4.56658243248678,-0.884182789642020) q[3];
u3(0.383848806428462,3.15917365306760,-2.84348366604245) q[7];
u3(0.969354395833654,-3.34898483875564,2.61833571409395) q[12];
cx q[12],q[7];
u1(1.70022070707910) q[7];
u3(-2.95995556142562,0.0,0.0) q[12];
cx q[7],q[12];
u3(2.52742788005857,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.83399266732033,-2.32670918195759,2.57948090841275) q[7];
u3(2.07803923760807,0.589252985809310,-5.45930497676599) q[12];
u3(2.70703607384370,-2.47135414122391,0.643927947693934) q[8];
u3(2.85030224904319,-4.28385127710895,-1.62497382218462) q[5];
cx q[5],q[8];
u1(-0.190580172999276) q[8];
u3(-1.89613555117442,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.534931806459997,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.15611271438137,3.57825145084811,-1.59460845792274) q[8];
u3(1.86993035693140,3.30648691872182,-0.154610900595006) q[5];
u3(2.46991272112248,2.60062281486289,-1.18769868754838) q[7];
u3(2.22492086951894,-0.710987488663632,-5.21892883024792) q[12];
cx q[12],q[7];
u1(4.07406645109624) q[7];
u3(-4.33076621095884,0.0,0.0) q[12];
cx q[7],q[12];
u3(-0.795271032292400,0.0,0.0) q[12];
cx q[12],q[7];
u3(2.28716527955763,-0.202423724700498,-2.37282052529190) q[7];
u3(2.15773927682039,-5.40637247530194,0.0469622083160881) q[12];
u3(0.493030359092560,0.0483250296248342,1.18449081357226) q[10];
u3(1.58966410069561,-0.318109987388677,-1.05082037224961) q[5];
cx q[5],q[10];
u1(-0.384415206903476) q[10];
u3(1.29248720171628,0.0,0.0) q[5];
cx q[10],q[5];
u3(3.77424463232245,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.82169153153774,-0.405242078248232,-2.52493621747780) q[10];
u3(1.92992447919369,-0.576376190751540,-1.26376322564453) q[5];
u3(1.87014624637134,-1.41621010717731,-0.366693377463311) q[0];
u3(1.49849558857689,-3.99938007109221,0.118686424450146) q[3];
cx q[3],q[0];
u1(3.61799747216040) q[0];
u3(-4.40250127726593,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.405373944770353,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.0687793012756816,-2.07803839948417,1.60175781580791) q[0];
u3(2.26606337751002,3.76917597733276,2.46873130606393) q[3];
u3(0.896744893481350,-0.320137990523863,1.25918448793428) q[11];
u3(0.336016857471056,-0.829829895638338,-0.280386911431511) q[1];
cx q[1],q[11];
u1(1.11947448272440) q[11];
u3(-0.712867421010667,0.0,0.0) q[1];
cx q[11],q[1];
u3(3.06123042658276,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.55073081591439,-4.19873991276575,0.173924644646601) q[11];
u3(1.54208446213743,-3.25487577846427,0.812255366901082) q[1];
u3(2.21318809040386,-1.74450854122990,4.09161553456043) q[8];
u3(0.536426913030035,3.09418892471852,-1.07660052774521) q[6];
cx q[6],q[8];
u1(4.25581533663854) q[8];
u3(-3.61069059254041,0.0,0.0) q[6];
cx q[8],q[6];
u3(-0.230275290332719,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.88045990409091,3.70000248973893,-2.25740026812708) q[8];
u3(3.03732187315922,1.48311572016502,0.103943432645797) q[6];
u3(1.31402210232352,0.152539047187536,-1.08350935074021) q[2];
u3(0.291893635843641,2.07296123137217,-3.92509096834570) q[4];
cx q[4],q[2];
u1(2.46206755751695) q[2];
u3(-1.96176427482858,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.397428868015314,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.25795686296926,-2.74451821400759,1.90111271774031) q[2];
u3(1.03184354763151,-0.0752476888220825,4.77763095315301) q[4];
u3(1.94301933266594,0.820137053740033,-0.844966993975229) q[3];
u3(1.96235262038055,0.241642239442399,-3.15955256934072) q[7];
cx q[7],q[3];
u1(3.03903759099994) q[3];
u3(-2.39451674629143,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.28146806424071,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.612105370243756,-0.407455437652051,3.24055844499941) q[3];
u3(1.54182191578487,4.58318746933655,-1.30230218317876) q[7];
u3(1.26717245563430,0.872764474428541,-0.875349631781204) q[5];
u3(1.60450271525072,-4.42314089525453,1.08066650636299) q[10];
cx q[10],q[5];
u1(1.63238747468710) q[5];
u3(-0.0936203204162305,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.07854479531951,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.21801097015266,-2.13814068155831,1.63871452185555) q[5];
u3(1.57140277220526,1.53531823640173,-2.97103092420613) q[10];
u3(0.819827438752669,-2.52839723422370,2.69012973159853) q[12];
u3(0.979782125935179,0.589234978589005,-2.44021013177701) q[1];
cx q[1],q[12];
u1(2.98228261045303) q[12];
u3(-1.71759248372452,0.0,0.0) q[1];
cx q[12],q[1];
u3(2.72281883841972,0.0,0.0) q[1];
cx q[1],q[12];
u3(1.85287770936071,1.38537464202007,-0.990389929850264) q[12];
u3(1.96093680065906,-1.03303555110732,-3.07070685883934) q[1];
u3(1.22293375959640,-1.42852559267179,-0.398770104874018) q[2];
u3(1.23053706588126,-3.89262619894174,0.0464017395832941) q[11];
cx q[11],q[2];
u1(2.01270099069996) q[2];
u3(-1.76773082556116,0.0,0.0) q[11];
cx q[2],q[11];
u3(0.360775731906320,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.98173091917900,-0.435121334823655,-2.53164554035835) q[2];
u3(0.981649206981811,3.08886322908347,0.153973784061014) q[11];
u3(0.850396409951933,2.62116569702265,-2.16302010475102) q[9];
u3(0.901704559427343,1.44777179888485,-1.99413307945251) q[0];
cx q[0],q[9];
u1(1.13850020768213) q[9];
u3(-0.179069935541042,0.0,0.0) q[0];
cx q[9],q[0];
u3(2.41391821083616,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.93408838696775,1.91172671010069,-2.37596804436189) q[9];
u3(0.852137671512886,4.83530523626397,0.995112501960160) q[0];
u3(1.60330988903670,2.57883674526569,-1.57212528553132) q[4];
u3(1.28824738146942,0.815142236783468,-0.732739225136551) q[6];
cx q[6],q[4];
u1(0.627179599949644) q[4];
u3(-0.439458238356892,0.0,0.0) q[6];
cx q[4],q[6];
u3(4.22016932993160,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.05109741220014,4.11364969868458,-1.52631858751745) q[4];
u3(0.610996376154537,0.636146171594441,0.744148946089971) q[6];
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
