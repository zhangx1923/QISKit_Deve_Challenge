OPENQASM 2.0;
include "qelib1.inc";
qreg q[13];
creg c[13];
u3(2.38969818100181,3.42361489387142,-1.16281618554022) q[4];
u3(1.85797120678039,1.89970874532896,-0.282594866216465) q[1];
cx q[1],q[4];
u1(1.81989451809360) q[4];
u3(0.624318183094174,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.990439747292628,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.27892905181766,-2.11886675116106,-0.0612239736771192) q[4];
u3(1.16733653444538,4.32953011988760,1.25736648994050) q[1];
u3(1.47358320966595,1.79136130524126,-3.09004539502760) q[8];
u3(2.15526444429223,2.40370986424368,-3.33822684668514) q[0];
cx q[0],q[8];
u1(3.51546444631105) q[8];
u3(-1.33135729786209,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.52079105048130,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.88470690793084,-3.39874153950477,2.40637884062517) q[8];
u3(1.38308798962617,-3.92908845401578,-0.180658673412274) q[0];
u3(0.819268908727141,2.03499693076944,0.780983805526978) q[3];
u3(1.06459790231736,1.17213187316623,-4.19068415289880) q[5];
cx q[5],q[3];
u1(3.24278814577941) q[3];
u3(-1.72069634468667,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.14442500337146,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.65517744835500,0.636303255090354,2.04552400904839) q[3];
u3(1.18520178657136,3.48061351215014,-2.46708206462220) q[5];
u3(1.83067170459720,1.01708151690145,-2.33336542586113) q[6];
u3(2.38901578314752,-2.97256356748760,3.09920953087541) q[12];
cx q[12],q[6];
u1(0.674474508318177) q[6];
u3(-0.259187360913893,0.0,0.0) q[12];
cx q[6],q[12];
u3(3.07623217971225,0.0,0.0) q[12];
cx q[12],q[6];
u3(0.142237862315164,1.73320576621982,-3.46311125512080) q[6];
u3(1.51084279494382,-1.69011311002213,0.328859980025241) q[12];
u3(0.688043944384849,-2.34854766152147,1.21007476285426) q[2];
u3(0.983374813040955,1.48508827071923,-3.18855825721124) q[10];
cx q[10],q[2];
u1(0.503066750910974) q[2];
u3(-1.56877823632007,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.78952581384472,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.75767844978882,-2.33672706001804,-1.09455400000965) q[2];
u3(0.565932004829380,0.648490797598214,-0.127671742719904) q[10];
u3(1.33964496131680,3.36334648979326,-0.506452944043392) q[7];
u3(0.257290807916935,1.47279931455857,-1.63306113929971) q[11];
cx q[11],q[7];
u1(0.897533328966561) q[7];
u3(-3.55698909235659,0.0,0.0) q[11];
cx q[7],q[11];
u3(1.50125912455923,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.27281222481083,0.217071343416379,-0.610901135353581) q[7];
u3(1.29308486722415,-1.24322596213300,-5.00951621429963) q[11];
u3(0.197692990830785,0.473155688951480,-1.11467897188413) q[5];
u3(1.33191121581343,-0.520498745656142,-1.73806791503682) q[7];
cx q[7],q[5];
u1(3.42581869730884) q[5];
u3(-0.814854290511661,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.61237452733508,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.780621536648350,-2.33938419351992,3.73615978044273) q[5];
u3(1.20226495341067,1.44037556899077,4.31092693047274) q[7];
u3(2.26520476459930,0.333907359246947,-3.22351689066855) q[9];
u3(1.06226544828419,-2.41835862733954,2.27001384936263) q[10];
cx q[10],q[9];
u1(3.45857480609348) q[9];
u3(-1.30539074033618,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.32254980817808,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.03216573396908,2.02147363960678,-3.24475616837618) q[9];
u3(0.435270574092543,-1.57166518045007,-4.02821692684827) q[10];
u3(2.43476086207292,0.105977143331275,-2.14315046225563) q[1];
u3(1.83266745065925,-4.00278875469323,2.08671349899569) q[6];
cx q[6],q[1];
u1(2.45066905823694) q[1];
u3(0.316952556427036,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.89473695736885,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.97021870959418,-4.76830365702999,0.670028226257596) q[1];
u3(1.42453308423567,-3.95138762614442,1.68295984994384) q[6];
u3(1.15852051493675,-0.441349916245096,1.32666182372132) q[3];
u3(1.19565887911751,-0.713516051777410,-1.64564029169209) q[12];
cx q[12],q[3];
u1(-0.227919201724986) q[3];
u3(-1.71782413943161,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.07279116211526,0.0,0.0) q[12];
cx q[12],q[3];
u3(2.16769242190458,-3.34249877633322,1.90491636642827) q[3];
u3(1.52626099770146,-2.28969425665643,0.514847664618310) q[12];
u3(2.60278982797830,1.42826285713789,-1.53425003110845) q[11];
u3(2.45912161874537,1.49382377910020,-4.45537784278627) q[4];
cx q[4],q[11];
u1(1.49407688768306) q[11];
u3(-0.243423511044367,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.45301621762546,0.0,0.0) q[4];
cx q[4],q[11];
u3(0.747915437890346,-0.877955858826496,1.16744253158463) q[11];
u3(3.01163748459477,1.99133067250383,-2.58232945247313) q[4];
u3(1.03925317274113,1.75955007519038,-3.21941610530926) q[0];
u3(0.731447007546906,-2.72497409257817,3.25609218337396) q[8];
cx q[8],q[0];
u1(1.64132715936311) q[0];
u3(-0.303448945557131,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.98777993900134,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.21360350759060,2.58470561324594,-2.94418000419383) q[0];
u3(0.622879843287618,-3.35411034589042,-2.36434779160884) q[8];
u3(0.373091854653299,2.86406634135657,-3.41046129060628) q[5];
u3(1.19960317188652,-0.109303559802605,-1.81582424281652) q[10];
cx q[10],q[5];
u1(-0.139207983797637) q[5];
u3(-1.08881696791766,0.0,0.0) q[10];
cx q[5],q[10];
u3(1.45770499878501,0.0,0.0) q[10];
cx q[10],q[5];
u3(0.478485086349861,1.87096744515032,-2.40086690490970) q[5];
u3(0.400246172809550,0.0420234596213751,-3.97883869248800) q[10];
u3(1.26219696458299,-1.04538881459303,1.40693164356188) q[11];
u3(1.01119910818521,-1.10396235927065,-2.24764283160922) q[3];
cx q[3],q[11];
u1(0.720468549122991) q[11];
u3(-1.65199286439187,0.0,0.0) q[3];
cx q[11],q[3];
u3(2.81590491912776,0.0,0.0) q[3];
cx q[3],q[11];
u3(0.653935110878119,2.95354978436236,-0.735602578948365) q[11];
u3(2.61789663483994,-3.47123299702716,-1.90368920311610) q[3];
u3(1.01003832702510,0.575742359368229,1.92188125283272) q[7];
u3(1.28461107671254,-1.64785586282374,-1.31790867634939) q[6];
cx q[6],q[7];
u1(1.62528507012162) q[7];
u3(-2.49105608665249,0.0,0.0) q[6];
cx q[7],q[6];
u3(3.13294113781348,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.31157479110606,-4.50397577171311,0.320987143810514) q[7];
u3(1.47805652766893,0.874856834980383,5.37792806272425) q[6];
u3(1.11515744860936,1.08534619488174,-0.401956904985291) q[1];
u3(2.43705691013598,-1.06244105356984,-4.45282578373785) q[12];
cx q[12],q[1];
u1(2.71310390345845) q[1];
u3(-1.89675120071467,0.0,0.0) q[12];
cx q[1],q[12];
u3(0.578527614237058,0.0,0.0) q[12];
cx q[12],q[1];
u3(1.30092612788371,2.88962501171963,-0.996970038453745) q[1];
u3(1.13926043860210,-2.05909462067483,-1.80938833128583) q[12];
u3(1.95253797612320,-0.896028205200555,-2.00685673501495) q[4];
u3(1.60586899499319,1.37794987541187,-4.39069463352876) q[2];
cx q[2],q[4];
u1(2.09947531953725) q[4];
u3(-1.68177425593919,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.56460762097856,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.820929715343182,1.12874231032203,-0.694718189337847) q[4];
u3(1.34313999577242,-4.18797329045758,0.758651154589141) q[2];
u3(0.240594477418099,1.69118485855464,-3.00950731350681) q[8];
u3(1.57999228288889,-2.30914140518949,2.91261213057873) q[9];
cx q[9],q[8];
u1(-0.0706051875659202) q[8];
u3(0.526150527716457,0.0,0.0) q[9];
cx q[8],q[9];
u3(4.19274182462822,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.01688246042296,1.20567137436282,1.80502645256156) q[8];
u3(0.987925866058539,0.214569934976990,1.97749551646620) q[9];
u3(1.52546043305775,2.58162124201406,-1.42462660758098) q[3];
u3(2.07853248857478,0.0666068873062002,-3.01059258152130) q[4];
cx q[4],q[3];
u1(2.29997305306133) q[3];
u3(0.373033020352508,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.25146681959981,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.63026836196563,-2.57067708913153,0.111833800849570) q[3];
u3(1.76449137570457,3.08187902052000,1.97390308874993) q[4];
u3(1.06363566648406,1.53318544446629,-3.58075911633258) q[11];
u3(1.20485549869205,2.62384676784474,-2.78048741975524) q[2];
cx q[2],q[11];
u1(-0.120302126071642) q[11];
u3(-2.40173766711951,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.37423979222807,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.63137049719455,0.582326852944944,-0.935908830239304) q[11];
u3(0.949160129765430,-1.08821028327340,3.43555243426469) q[2];
u3(2.37764876002705,1.18873615465006,-3.44232388536492) q[1];
u3(1.81648922357838,3.26640386151435,-2.74290852664926) q[0];
cx q[0],q[1];
u1(2.25450123726497) q[1];
u3(-3.02707601755615,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.69907429243413,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.62524274944485,3.45308248091722,-2.09434406822768) q[1];
u3(1.45879165553800,3.29420986607673,-2.07595312779942) q[0];
u3(0.201714019019488,2.80256945329852,-2.66139328794190) q[8];
u3(1.24061387765054,1.70753218357275,-1.80104792442256) q[12];
cx q[12],q[8];
u1(3.36816861507787) q[8];
u3(-1.18377974357507,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.68792516168995,0.0,0.0) q[12];
cx q[12],q[8];
u3(1.54775059291293,-1.03872999072451,0.964881531683566) q[8];
u3(2.59139594171039,1.10045875528536,2.33258368197291) q[12];
u3(2.76872763957261,0.189236248280979,-2.49152418736823) q[10];
u3(2.28726933914971,3.48382376897331,-0.145263474009599) q[5];
cx q[5],q[10];
u1(2.08961578434065) q[10];
u3(-1.75236864950686,0.0,0.0) q[5];
cx q[10],q[5];
u3(3.67616761225319,0.0,0.0) q[5];
cx q[5],q[10];
u3(0.783535503865800,2.62585608803304,-0.511362301268180) q[10];
u3(1.14796030895321,0.442147196703906,-4.61551415847784) q[5];
u3(1.41463675882538,0.669997518860979,-2.70767355322040) q[9];
u3(1.14086478310856,-3.14984435647689,2.65254030244811) q[6];
cx q[6],q[9];
u1(1.05170359375235) q[9];
u3(-1.70581873128307,0.0,0.0) q[6];
cx q[9],q[6];
u3(-0.480392960085135,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.55983844840181,0.0890980761225260,-0.244697497460907) q[9];
u3(1.01685688425951,-2.70890531184408,3.41160688133866) q[6];
u3(2.23067207425182,0.146838180299717,-1.88250477577077) q[5];
u3(1.67813866457270,0.534122788342150,-4.17380655107202) q[6];
cx q[6],q[5];
u1(0.792998575327924) q[5];
u3(-0.499234742079566,0.0,0.0) q[6];
cx q[5],q[6];
u3(3.06415165850779,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.994227112518769,-3.93415061929210,1.16043438009455) q[5];
u3(0.968724570758339,0.550698805439407,-4.56265711862489) q[6];
u3(1.90417363913593,0.660425346568494,1.23272498859246) q[8];
u3(1.98446609855204,-1.64759715743374,-1.84967396169360) q[9];
cx q[9],q[8];
u1(0.00367903656250568) q[8];
u3(-1.38060099170002,0.0,0.0) q[9];
cx q[8],q[9];
u3(0.706998550770926,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.82923685665581,0.505066711959612,-1.90890257799084) q[8];
u3(1.79029244899559,-4.60486880154867,-1.49453807788711) q[9];
u3(1.91029560475173,1.13497805982993,0.821772301241489) q[7];
u3(1.60665056269812,-1.68335073880158,-1.54350334245401) q[1];
cx q[1],q[7];
u1(3.73824944476088) q[7];
u3(-3.57901164176494,0.0,0.0) q[1];
cx q[7],q[1];
u3(-1.02091614203450,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.938677482256052,-3.01480036262208,-0.599316948552924) q[7];
u3(2.53521462276574,-0.0645003582457832,-1.98754160517817) q[1];
u3(0.973427925080800,0.564504757662636,1.30158685044266) q[2];
u3(1.16238256092123,-1.35437643469389,-0.631614287506401) q[0];
cx q[0],q[2];
u1(0.0534135157578861) q[2];
u3(-1.54582869799826,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.70549094126188,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.20570475905529,-1.33474251381729,-0.682222698204828) q[2];
u3(2.27390715370640,-1.27906561107852,-1.01195715007633) q[0];
u3(2.26915281963807,0.782118958376092,0.228314886471841) q[11];
u3(0.706764925274102,-3.64705001451017,-0.589320316854415) q[10];
cx q[10],q[11];
u1(1.19390561779907) q[11];
u3(-0.683056902691259,0.0,0.0) q[10];
cx q[11],q[10];
u3(-0.282931013020730,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.61716846927792,-1.66061282791306,2.52204822640154) q[11];
u3(1.31813514703430,0.0468145308327006,-3.56880094865143) q[10];
u3(1.71407542572711,1.26676807197210,-3.90494220283392) q[12];
u3(2.02968134707225,-1.35443433508393,4.67869445682778) q[3];
cx q[3],q[12];
u1(-0.0253934438442009) q[12];
u3(-2.45878024882937,0.0,0.0) q[3];
cx q[12],q[3];
u3(1.47713616497523,0.0,0.0) q[3];
cx q[3],q[12];
u3(0.950772942519389,-1.12934136853125,1.93794988874913) q[12];
u3(2.41627551485860,0.881417611307111,0.843221471063438) q[3];
u3(0.413553718934132,2.50208326869978,-2.38442408056848) q[3];
u3(1.11127239767991,-0.365547319711963,-0.998283224502661) q[12];
cx q[12],q[3];
u1(0.336651233310223) q[3];
u3(-1.90706100017387,0.0,0.0) q[12];
cx q[3],q[12];
u3(1.37403850744645,0.0,0.0) q[12];
cx q[12],q[3];
u3(1.22747562479765,2.09587514520550,-0.149423545551443) q[3];
u3(1.70283539458438,-0.418154152506145,-2.70308384141443) q[12];
u3(0.570616221524674,-0.0739361574272802,-0.00137185711418746) q[11];
u3(1.36070386879414,-3.02866833549755,1.79832099266431) q[9];
cx q[9],q[11];
u1(1.57801880557079) q[11];
u3(-2.40408144890980,0.0,0.0) q[9];
cx q[11],q[9];
u3(3.81616442798839,0.0,0.0) q[9];
cx q[9],q[11];
u3(1.72677189711814,2.39442311409534,-2.09739456604626) q[11];
u3(0.794733674688616,-0.00556762537865696,4.07424212252444) q[9];
u3(1.83896680255628,-1.06719108763083,-0.530557683778645) q[2];
u3(0.602369738207709,-4.29767575852475,0.0260530083499826) q[0];
cx q[0],q[2];
u1(0.287653015702201) q[2];
u3(-1.33063229314497,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.61347269542091,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.53377549895481,3.03500451271655,-0.928222300291645) q[2];
u3(1.80541359063716,-1.96453917899340,-3.04020927406771) q[0];
u3(0.765071138534033,-2.05751603394453,0.370587467965406) q[10];
u3(2.12663052840485,-4.31501540951550,0.553246291750995) q[7];
cx q[7],q[10];
u1(1.74081786489380) q[10];
u3(-2.94143915708616,0.0,0.0) q[7];
cx q[10],q[7];
u3(3.15156592949334,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.389136314261778,-0.717775653985304,3.61808261270373) q[10];
u3(0.940586338848318,-4.82677501838364,-1.14262589459256) q[7];
u3(1.94290981228340,0.0482877933039103,0.712000344059009) q[8];
u3(0.571821384690610,-2.12679388807988,-1.91285324361491) q[5];
cx q[5],q[8];
u1(3.11291508297355) q[8];
u3(-1.35866090910521,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.61540511212112,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.58385956008944,-2.87765659717454,1.33995443747601) q[8];
u3(0.502738486071437,-0.0225806914362459,-0.504264355366343) q[5];
u3(1.33741072202910,2.33124877615975,-3.01567360922543) q[4];
u3(0.225589555240356,-2.72607659752947,2.50744230459682) q[6];
cx q[6],q[4];
u1(-0.0769627650992102) q[4];
u3(-1.55060120755207,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.90225986968527,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.15701085996299,3.97019765232626,-2.08901174619872) q[4];
u3(1.13374953682221,1.64766418457096,-3.26025715287530) q[6];
u3(1.75095487022931,1.70562847245480,-0.281103721987537) q[7];
u3(2.64240486797236,1.15034156158253,-3.48060310427409) q[6];
cx q[6],q[7];
u1(1.93365969744776) q[7];
u3(-2.97551342401289,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.563302516138775,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.34378256071356,0.696415123488884,-1.66370546138680) q[7];
u3(2.45530789804339,-0.492859683085479,0.680614666143502) q[6];
u3(0.783523867147622,-2.22217864595307,0.848405756977485) q[1];
u3(0.829624474297203,-2.38641791325510,0.916693584841088) q[5];
cx q[5],q[1];
u1(1.54434103069566) q[1];
u3(-3.45169375720951,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.17799998448509,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.731977101584726,-0.0320193735640648,2.62192156731186) q[1];
u3(2.49327272821207,2.03169877673858,-1.10213661670536) q[5];
u3(1.84046829287479,-1.49415240138976,-0.945829090757465) q[2];
u3(1.96027313819825,-1.95869522461649,0.339272084642329) q[3];
cx q[3],q[2];
u1(-0.703328828041320) q[2];
u3(1.10036295166855,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.58312425492514,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.29698214912059,-3.05349686745461,2.42203190763249) q[2];
u3(2.45278990284866,-0.414415719563351,0.276749270254285) q[3];
u3(2.36679500629007,2.05160287799715,0.189762209809494) q[8];
u3(1.88763841462923,-0.504024523361360,-2.43048620964404) q[11];
cx q[11],q[8];
u1(0.866477184024407) q[8];
u3(-0.319995288276772,0.0,0.0) q[11];
cx q[8],q[11];
u3(1.52898480859498,0.0,0.0) q[11];
cx q[11],q[8];
u3(1.60744448332494,1.72307514308171,-1.84842819462194) q[8];
u3(2.20065068390085,-4.25902892076176,-1.38431208063917) q[11];
u3(2.20390960323987,0.514642926936043,-2.35371159255081) q[0];
u3(2.23087700308338,-3.67920314416584,2.37324557772371) q[9];
cx q[9],q[0];
u1(0.172061558119426) q[0];
u3(-1.91368072370233,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.12387402721714,0.0,0.0) q[9];
cx q[9],q[0];
u3(0.942509713268636,-0.814703701664284,2.12866862435983) q[0];
u3(2.47359009508437,-1.45823165043194,-0.785552878433849) q[9];
u3(1.66484055475040,1.88554668531172,-2.98188643094054) q[4];
u3(2.28923667474684,1.88868982047931,-3.94618704078474) q[10];
cx q[10],q[4];
u1(0.663346901310725) q[4];
u3(-3.21336029599558,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.05068829853294,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.83959554121372,-1.72503726499257,1.41489321038355) q[4];
u3(0.513367581068066,4.27254552159612,1.21715222571719) q[10];
u3(1.75915241956564,2.53792599084513,-2.07598541903905) q[4];
u3(2.79716523586717,1.55232250654861,-0.0896399920822974) q[10];
cx q[10],q[4];
u1(2.72011823013592) q[4];
u3(-1.64290610003700,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.661761723550714,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.83646701362873,-0.306419508149018,-2.94900602555925) q[4];
u3(0.748593679022686,1.15000346053310,-2.06494621563868) q[10];
u3(0.383864957021950,3.11827857262951,-2.78409977046861) q[1];
u3(0.961217987965177,-0.645530941888551,-1.47460275486562) q[12];
cx q[12],q[1];
u1(0.153147564240057) q[1];
u3(-2.25293922759079,0.0,0.0) q[12];
cx q[1],q[12];
u3(1.31149347335654,0.0,0.0) q[12];
cx q[12],q[1];
u3(2.24552067318981,1.08775921722681,1.82726727765963) q[1];
u3(2.71441848880849,0.289460176376341,-3.90405797308944) q[12];
u3(1.85141227127402,2.85039563985463,-1.54797300287782) q[2];
u3(2.70023360604686,0.299851513586850,-1.75376339625186) q[8];
cx q[8],q[2];
u1(1.65983610778010) q[2];
u3(-2.43918575011530,0.0,0.0) q[8];
cx q[2],q[8];
u3(0.450575444965172,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.15202623171383,1.60432374896128,-2.38823179702837) q[2];
u3(1.67249554609013,-4.60914968608198,-1.18931330741198) q[8];
u3(1.62891809628075,2.11284104420396,-2.77096370733183) q[0];
u3(1.98962631943886,-2.35773355397985,3.12173199312677) q[11];
cx q[11],q[0];
u1(-0.231171266953556) q[0];
u3(-1.70081655276418,0.0,0.0) q[11];
cx q[0],q[11];
u3(1.96941949216300,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.59159498840786,0.249176207163629,1.08927514837969) q[0];
u3(0.834143158334138,-2.68701724104218,-3.21568163380257) q[11];
u3(1.73673003158423,0.336167139377396,2.45523669727551) q[7];
u3(0.785490991790974,-2.98382888696093,-2.76145571450172) q[3];
cx q[3],q[7];
u1(-0.256175060340896) q[7];
u3(-1.83218517809913,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.54333020519738,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.62729680895516,1.08812835240500,-1.65481834430839) q[7];
u3(0.869559012002601,-0.858467701973668,1.19470949444961) q[3];
u3(1.81772902513897,1.32032872936890,0.388829016637645) q[9];
u3(1.20605412035689,1.01053919859101,-2.36406826127502) q[6];
cx q[6],q[9];
u1(1.75733399757823) q[9];
u3(-2.53047456392063,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.970613076091753,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.71024505435904,-1.83523049044572,3.25310588267898) q[9];
u3(1.67131097618399,4.20053309234601,0.734653024761753) q[6];
u3(1.98421898368478,2.29996335947509,-1.70726429065638) q[8];
u3(1.65973198967096,1.08683450531556,-1.06133808206117) q[10];
cx q[10],q[8];
u1(1.61919059403253) q[8];
u3(-2.71747716094163,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.520850690207074,0.0,0.0) q[10];
cx q[10],q[8];
u3(1.44727081866489,0.780024104628202,0.398159203970349) q[8];
u3(1.31329929174356,1.96652460068941,1.72638573872888) q[10];
u3(0.255161003800498,-0.402409930846970,1.33241066292161) q[7];
u3(0.246557759291468,-2.51235001388156,0.371779460350421) q[11];
cx q[11],q[7];
u1(-1.10824613284033) q[7];
u3(0.283147124241281,0.0,0.0) q[11];
cx q[7],q[11];
u3(3.79226610502254,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.69624021925360,-0.572918924688244,2.41753043758163) q[7];
u3(2.37837046450512,-2.04664461278527,-3.99244856967497) q[11];
u3(2.33524327674056,0.196339968051246,1.21456589656104) q[5];
u3(1.39852388629482,-2.52555723225519,-2.31954117382615) q[9];
cx q[9],q[5];
u1(3.26681796721781) q[5];
u3(-1.44616961546299,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.65202554868865,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.54801763022478,1.93730647583840,-0.260893057669525) q[5];
u3(2.89321162020883,-2.69215850455280,-2.10810177125981) q[9];
u3(2.83045354904623,0.748436304665443,-1.50399719883335) q[4];
u3(2.94241476640477,3.02385594274657,-1.77393427240269) q[6];
cx q[6],q[4];
u1(3.65274033501405) q[4];
u3(-4.50309460566310,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.703333798542841,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.44904813122844,-2.92457600424703,-0.348017444235045) q[4];
u3(0.509747215914694,-3.46284487934783,-2.17865422790841) q[6];
u3(1.37223843868667,0.775636088342897,-2.67882126096773) q[1];
u3(0.277591120624377,2.69719091177187,-3.38148256290238) q[2];
cx q[2],q[1];
u1(0.0355425012936741) q[1];
u3(-1.59681268575295,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.21088418426670,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.08965614906248,0.0866236062264675,-0.103768867320585) q[1];
u3(1.86268159916940,3.21598479176596,-1.52724885732885) q[2];
u3(0.973681433953906,-2.05563446953581,3.38274937409545) q[0];
u3(2.01342406126076,1.97792977692498,-1.86095209206065) q[12];
cx q[12],q[0];
u1(2.27350372492081) q[0];
u3(0.270678269321485,0.0,0.0) q[12];
cx q[0],q[12];
u3(1.46486703354682,0.0,0.0) q[12];
cx q[12],q[0];
u3(0.630921963948531,2.76663263881504,0.0571934618591057) q[0];
u3(2.18518118202075,-1.55139497101863,2.21615235531882) q[12];
u3(1.76832692488305,-3.91876854705177,2.14391252587762) q[0];
u3(0.324085054431934,2.70048562320842,-0.465998340682953) q[5];
cx q[5],q[0];
u1(1.71734343329384) q[0];
u3(-2.94808641614788,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.735999643212608,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.64554886944868,-4.67702180609987,1.53143530743681) q[0];
u3(2.52342889929919,-4.17731078510756,-0.928450581989901) q[5];
u3(0.909035379979852,3.00696399306173,-1.80148557827586) q[2];
u3(1.38833486300168,1.26922052599361,-2.04865242449856) q[7];
cx q[7],q[2];
u1(0.661444710683501) q[2];
u3(-1.27934211369666,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.111526259480863,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.18559290547839,0.293405600799323,-1.07840542252147) q[2];
u3(2.36250958948849,3.78027417074248,2.19788309354114) q[7];
u3(2.57682729557479,-2.39759030786936,1.05874452296565) q[6];
u3(2.12122369120691,1.89620399731711,4.00984499361615) q[1];
cx q[1],q[6];
u1(1.19037541727914) q[6];
u3(-0.298994230383885,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.20882743313869,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.38685030383961,1.67090736572328,-2.04251892543612) q[6];
u3(1.82491761828702,-1.28901528990227,4.08547288180487) q[1];
u3(1.97621231481148,-2.26482233453605,0.313946647890744) q[12];
u3(2.05621460264064,-3.34756632645142,0.599934803647882) q[11];
cx q[11],q[12];
u1(1.68479226687223) q[12];
u3(0.361947354638069,0.0,0.0) q[11];
cx q[12],q[11];
u3(1.16421947809794,0.0,0.0) q[11];
cx q[11],q[12];
u3(2.03989599265922,-0.768997396055476,3.15917150706137) q[12];
u3(1.55062899549663,-0.110232798173979,3.13122506857735) q[11];
u3(2.25253011594442,4.28106151355845,-1.77524029256147) q[8];
u3(1.22309818425484,-1.18645222229910,2.14367256478475) q[3];
cx q[3],q[8];
u1(3.07386120703501) q[8];
u3(-2.10605028845336,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.608011577928936,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.73362189951845,-2.54687990943357,1.35462257719270) q[8];
u3(2.37829729270251,0.355702818540141,1.65535592234769) q[3];
u3(1.63358416186499,0.783502409299009,-3.60718842023361) q[10];
u3(1.99179467176996,3.01080579934767,-2.03352899777792) q[9];
cx q[9],q[10];
u1(3.40202788643450) q[10];
u3(-1.18314540728225,0.0,0.0) q[9];
cx q[10],q[9];
u3(1.66962169924375,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.49714579185654,-0.464822705766692,1.69187370130672) q[10];
u3(2.25789582013210,5.57981781813270,0.636355670879163) q[9];
u3(2.69624525481538,1.50232204348572,-1.27781983037862) q[7];
u3(2.09547180392861,-0.671523581701268,-5.26574657589035) q[5];
cx q[5],q[7];
u1(2.09493055592743) q[7];
u3(-3.23407824283149,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.592908603829362,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.32175753857310,4.22073304858050,-1.61367339837027) q[7];
u3(1.98366537114037,1.65457091986067,0.128506085188053) q[5];
u3(3.06631515892703,-1.49459678083638,-1.60526811995028) q[6];
u3(1.48288934112645,0.534230525943070,-5.37610643657977) q[3];
cx q[3],q[6];
u1(2.40185732907600) q[6];
u3(-2.84860443384820,0.0,0.0) q[3];
cx q[6],q[3];
u3(1.15741532086265,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.00644336925700,0.807052963126767,-1.48775966247739) q[6];
u3(2.42887449573799,-0.703421336582682,-2.33360988236281) q[3];
u3(1.44557526196559,-0.835923463775136,0.823531433116012) q[12];
u3(2.01221051995437,-1.85254828926856,-1.44892516832972) q[10];
cx q[10],q[12];
u1(3.71852818780852) q[12];
u3(-1.26133440934382,0.0,0.0) q[10];
cx q[12],q[10];
u3(2.15390679995381,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.00633005930986,2.18852208506118,-0.615844228732477) q[12];
u3(1.54055745294554,4.44658395246742,1.39914725338289) q[10];
u3(2.56724967992490,-2.77589164551613,-0.0885671673084485) q[0];
u3(2.04942926612390,0.430060238219014,1.13869737107767) q[4];
cx q[4],q[0];
u1(2.96980515517029) q[0];
u3(-1.29575130409545,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.277121205622416,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.52907664706007,-0.768099247261042,0.566373122637175) q[0];
u3(0.866460487582165,-0.488050512794823,2.15490160788505) q[4];
u3(2.55763272583337,3.62906563739338,-1.15783428281018) q[9];
u3(1.46899743494978,1.59792550982887,-0.470429025221276) q[1];
cx q[1],q[9];
u1(-1.06088630922821) q[9];
u3(0.587103539094985,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.24341232191816,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.10120656071406,-1.84402553027621,0.438528864968744) q[9];
u3(1.47957238084702,1.18457933165373,1.08836618406770) q[1];
u3(1.55491217286631,-4.27999664092886,1.17750103820839) q[2];
u3(0.937447934624887,-0.666016414265554,3.40380317837447) q[8];
cx q[8],q[2];
u1(1.41360780183508) q[2];
u3(-3.33369288511163,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.61451486511233,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.02481868301750,2.02428510811103,-3.58379984684953) q[2];
u3(0.463307730326767,-0.720936089552796,5.09529170925978) q[8];
u3(2.82473273666028,-0.216747828943423,0.410231100549431) q[2];
u3(0.432838067783909,-2.90264291731655,-1.83902285523536) q[7];
cx q[7],q[2];
u1(-0.0204099254716987) q[2];
u3(-1.73675702094001,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.995595363212355,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.801772966434711,-0.765444142743097,0.818339329177464) q[2];
u3(1.75790371461983,2.54406598660638,3.39475892775864) q[7];
u3(1.47710557803384,-0.660287884325166,0.399756186141442) q[12];
u3(1.65911527026643,-2.44427653382420,-1.39746507879837) q[1];
cx q[1],q[12];
u1(1.20894771981142) q[12];
u3(-0.773680265483134,0.0,0.0) q[1];
cx q[12],q[1];
u3(2.61004318527707,0.0,0.0) q[1];
cx q[1],q[12];
u3(2.24190904316569,-1.48781645802769,2.14172031829945) q[12];
u3(0.794623700278842,-1.92598125270803,-2.28157196180723) q[1];
u3(1.64030649271536,0.769505376490114,-2.66937698903466) q[8];
u3(1.22684381759406,-3.03086948530292,2.66302017633891) q[10];
cx q[10],q[8];
u1(2.69707126187628) q[8];
u3(-1.50784984085005,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.349037487148724,0.0,0.0) q[10];
cx q[10],q[8];
u3(2.81709901645123,-2.97073335428013,-0.210155006110009) q[8];
u3(0.748795270341119,1.30720238530510,1.30184512773537) q[10];
u3(1.90039558968505,1.08936517942872,-3.20118424426107) q[6];
u3(2.34495372065604,2.09192571004355,-2.83683848647687) q[3];
cx q[3],q[6];
u1(1.20988142020418) q[6];
u3(-0.149809709090947,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.22575422708346,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.85542967546667,3.87985407096649,-1.41220619858819) q[6];
u3(0.373717640720573,0.644367379841136,-1.59565400162498) q[3];
u3(2.23040599840492,2.54966616386100,-2.12433505133389) q[5];
u3(1.80547296745448,-3.38752096378239,2.53604818902705) q[11];
cx q[11],q[5];
u1(0.906942485255311) q[5];
u3(-3.39051695381043,0.0,0.0) q[11];
cx q[5],q[11];
u3(1.68534678427016,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.72021517127850,-1.60790388500351,1.65268322148100) q[5];
u3(0.188487187267978,-1.54878178729083,2.18113154023645) q[11];
u3(2.49316274121929,2.82210454933157,-1.41062844916964) q[0];
u3(2.02728602614969,-0.532192386456292,-5.46837156324657) q[4];
cx q[4],q[0];
u1(1.06254760671955) q[0];
u3(-3.60199564784598,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.13098847439504,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.687675148823248,1.65233734332049,-4.30018052780335) q[0];
u3(2.12376137276606,0.919270949307581,-4.19254586964858) q[4];
u3(1.92443492762554,-0.697939469986603,-2.12689362277900) q[3];
u3(1.65148402422702,-3.91858017872532,1.92969770788522) q[1];
cx q[1],q[3];
u1(2.83444998719952) q[3];
u3(-2.54022388847221,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.25872295116546,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.58650437745783,0.828428807007344,0.443063466300076) q[3];
u3(2.36610192077370,2.12776510311011,1.95132088871022) q[1];
u3(1.90636268583430,-0.186485192929164,2.02869851705897) q[4];
u3(1.52384334088264,-0.0211690752844320,-0.101979536337928) q[10];
cx q[10],q[4];
u1(1.67756837003024) q[4];
u3(-2.42205352572280,0.0,0.0) q[10];
cx q[4],q[10];
u3(0.322290560899590,0.0,0.0) q[10];
cx q[10],q[4];
u3(2.38938936615624,3.51413470639252,-0.817562846919622) q[4];
u3(1.02608554975162,-3.34810403881909,-0.878394388484781) q[10];
u3(0.837263650066828,1.93303537029988,-1.70898418859635) q[2];
u3(1.14230859761276,1.26868537950551,-1.30692619696821) q[12];
cx q[12],q[2];
u1(2.94621908846001) q[2];
u3(-1.82548170851655,0.0,0.0) q[12];
cx q[2],q[12];
u3(0.413029334079210,0.0,0.0) q[12];
cx q[12],q[2];
u3(1.30696371480777,-2.67057019409363,2.24997101928555) q[2];
u3(1.85686676496312,-2.67400456219284,-1.49116792947616) q[12];
u3(1.52096281986420,1.05989990380882,1.56776608405504) q[11];
u3(0.705009809183958,-0.369331736091521,-2.78505357758275) q[5];
cx q[5],q[11];
u1(2.28116924496381) q[11];
u3(0.189280680613785,0.0,0.0) q[5];
cx q[11],q[5];
u3(1.51603932788925,0.0,0.0) q[5];
cx q[5],q[11];
u3(2.62689608822492,-4.08430690010279,1.64841089313745) q[11];
u3(1.76574811090048,3.29972284032186,0.367797194755551) q[5];
u3(0.451140211541425,1.49921448535773,0.176473099006333) q[7];
u3(1.86448834971043,-0.327915881940495,-3.57277853548436) q[0];
cx q[0],q[7];
u1(-0.0520325710012914) q[7];
u3(-2.40390997404480,0.0,0.0) q[0];
cx q[7],q[0];
u3(0.959331021116783,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.75244239098239,2.20588112837300,-0.950986888527692) q[7];
u3(0.575655666088398,0.167444363669309,-5.91897893397377) q[0];
u3(0.961323615118267,-3.74182848499279,2.13519331745444) q[9];
u3(1.96750567437246,-2.46379942650126,3.24790750724277) q[6];
cx q[6],q[9];
u1(2.73929237955047) q[9];
u3(-1.97586385528422,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.56627712583282,0.0,0.0) q[6];
cx q[6],q[9];
u3(2.07050305112520,4.02742048711908,-1.51926115246315) q[9];
u3(2.07690903542637,3.42338379014859,-0.188989724363441) q[6];
u3(2.05761288162167,-0.0338007380040684,-2.61880163865494) q[5];
u3(2.01276880686275,0.0988741886604827,-4.66358612463524) q[3];
cx q[3],q[5];
u1(1.17821801968602) q[5];
u3(-3.26419537431812,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.68270683097111,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.43973799145937,-1.39313253081616,-0.733990004596927) q[5];
u3(0.900327330057154,1.04172348174427,4.85962631277870) q[3];
u3(1.52669355277390,-0.0767197467786404,-0.701761869580303) q[7];
u3(2.03200832541833,-5.01769246741726,0.980728373371716) q[10];
cx q[10],q[7];
u1(1.30621580070002) q[7];
u3(-3.17169889157163,0.0,0.0) q[10];
cx q[7],q[10];
u3(2.47283843444818,0.0,0.0) q[10];
cx q[10],q[7];
u3(2.52024324980745,-0.221207909136663,0.0213609668237459) q[7];
u3(2.04568415443740,-1.45265174315365,0.252283085119409) q[10];
u3(1.37340302702278,3.67163615330599,-0.576680008098170) q[11];
u3(2.23695659814151,3.30606672041011,-0.902240583294927) q[1];
cx q[1],q[11];
u1(0.971210151190229) q[11];
u3(-0.549708813855934,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.93316737832293,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.861182626118853,0.797643376214100,2.07936095252522) q[11];
u3(2.84422161778761,4.33504760703654,0.339595413180749) q[1];
u3(1.39278611189939,-2.15014046642638,-0.739919442074935) q[4];
u3(2.12948624745371,-4.20491285036118,-1.13231521290710) q[9];
cx q[9],q[4];
u1(2.92124853442990) q[4];
u3(-1.78692690664444,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.871220557768709,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.03307100972163,-0.713131330010890,3.14116843686858) q[4];
u3(0.375536271991239,-3.18468055086725,-0.201922994820759) q[9];
u3(2.03407336911817,-1.78208815355229,1.29405491112896) q[0];
u3(2.27707451826161,-3.15954388408103,0.0848775946790867) q[6];
cx q[6],q[0];
u1(3.43468787060966) q[0];
u3(-1.39548966133810,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.98802906726039,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.24744938575735,2.56693453179746,-3.58098629933617) q[0];
u3(1.06630600476260,-3.40848569904314,1.98454508561907) q[6];
u3(1.65209974753902,-1.83538203745082,-0.414507592807616) q[2];
u3(0.707563996236371,-3.48760494609069,0.904788170103622) q[12];
cx q[12],q[2];
u1(2.43761956098451) q[2];
u3(-1.92404373759531,0.0,0.0) q[12];
cx q[2],q[12];
u3(0.524445729194767,0.0,0.0) q[12];
cx q[12],q[2];
u3(1.13082608438755,-3.42282014425500,2.80180580136377) q[2];
u3(1.19350042116353,6.17829721203818,-0.0635274718424106) q[12];
u3(1.48220316981199,0.445763382703934,-3.29422091716090) q[5];
u3(0.305811516243401,-3.29308369363301,2.79561199521370) q[3];
cx q[3],q[5];
u1(1.51928645552202) q[5];
u3(-0.107658425830739,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.11618660762921,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.83651324505677,-0.584334576543265,-2.60191866658074) q[5];
u3(1.13593511711591,-4.29633664563307,1.56744632643863) q[3];
u3(2.73952208762843,-0.304090838035572,-0.749996380359570) q[11];
u3(1.31339245307178,0.237931777710381,-5.43789662568424) q[12];
cx q[12],q[11];
u1(1.65595834919687) q[11];
u3(-0.000639362691750245,0.0,0.0) q[12];
cx q[11],q[12];
u3(0.703679547028903,0.0,0.0) q[12];
cx q[12],q[11];
u3(1.47158457250897,1.13187326552721,-1.84120850691752) q[11];
u3(0.612515771383670,-1.89186359536138,-1.25216085942441) q[12];
u3(2.14829360928824,2.33512498209820,0.520182805187399) q[6];
u3(2.26813050446631,0.147620762556684,-2.08845361925624) q[8];
cx q[8],q[6];
u1(1.56027095126762) q[6];
u3(-3.40240744467418,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.18754839894833,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.86688547512199,0.460434423377493,-0.741653127267221) q[6];
u3(2.23239547703527,-0.743677792811976,1.27271804470414) q[8];
u3(1.11080202612359,2.74516478996948,-1.37310622285867) q[1];
u3(0.661707743429894,1.20856485917684,-3.27348389708403) q[0];
cx q[0],q[1];
u1(2.06499154274987) q[1];
u3(-2.69397934996988,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.42675123941203,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.02174013581022,-2.46087651039621,3.27905526557149) q[1];
u3(1.66709926470216,1.53678469011503,-1.71745300663044) q[0];
u3(0.615068719065415,0.367750118034934,-2.21475617468379) q[4];
u3(0.944861038955926,-3.85145256692185,2.42511730899677) q[7];
cx q[7],q[4];
u1(2.41116006534042) q[4];
u3(-1.84827430894202,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.142427304602262,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.16912908696026,-0.591931081339246,-2.38429402480052) q[4];
u3(1.90058448214741,-0.414927366469906,5.57349215109303) q[7];
u3(2.10320023791179,-4.36901393069237,1.45970700812115) q[10];
u3(0.481517024495782,1.59622604472978,0.0501316262107533) q[2];
cx q[2],q[10];
u1(2.94019707869342) q[10];
u3(-0.938319571704406,0.0,0.0) q[2];
cx q[10],q[2];
u3(1.76134711813047,0.0,0.0) q[2];
cx q[2],q[10];
u3(2.66057096852577,-2.07253021411571,0.346393947060722) q[10];
u3(0.391435196005403,-2.12186637282007,3.44879302293094) q[2];
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
