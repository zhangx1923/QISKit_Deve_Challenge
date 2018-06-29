OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.20877772938985,-2.64262148203699,1.65459937598413) q[9];
u3(0.922259846410419,0.927257969854409,-2.95123036323686) q[3];
cx q[3],q[9];
u1(1.40389287451312) q[9];
u3(-3.58990346685787,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.24556694886895,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.86150093650884,-2.40933640928055,3.12093868970662) q[9];
u3(1.30758484731973,-4.52078410324949,0.821889298573075) q[3];
u3(1.43282651599412,-0.0409757743334429,2.93743398685194) q[2];
u3(0.619198321582677,-0.748914152124556,-1.98930708729791) q[5];
cx q[5],q[2];
u1(1.39432823833672) q[2];
u3(-3.35588527121814,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.99306618279145,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.678785541306475,0.400265351258482,-0.107277781869699) q[2];
u3(1.08467039199761,-0.116279858762437,-2.30660344419544) q[5];
u3(2.07905888531849,1.89100933717250,-2.53007453606376) q[4];
u3(1.42995491416611,-2.72093534042688,2.78431470756855) q[7];
cx q[7],q[4];
u1(1.56584840117391) q[4];
u3(-1.16550446853857,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.418551178668853,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.26748545241010,2.36826295130826,1.50409794129924) q[4];
u3(1.96329795697372,-0.751917189529958,3.42903744679538) q[7];
u3(1.15613995932363,1.84924110762411,-0.0697555224964018) q[0];
u3(1.34494192135288,1.48871675966959,-1.03618522234815) q[1];
cx q[1],q[0];
u1(2.46140120427138) q[0];
u3(-1.76339395067532,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.110472790984944,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.10488275587230,2.94500633900061,-3.05501179955230) q[0];
u3(2.73858222864375,4.08131070506084,2.00905372027260) q[1];
u3(1.70388016133385,1.16164132782746,-2.89435808806370) q[8];
u3(1.44263077587623,-3.20857055480006,2.94358237016508) q[6];
cx q[6],q[8];
u1(1.94452000213792) q[8];
u3(-3.33438857126493,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.00525150082689,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.69962164257080,0.705648479356728,-4.60576701474903) q[8];
u3(0.806183159923383,-0.0632659297093490,0.217484261356184) q[6];
u3(2.92705427638880,0.307993644647889,-1.92736019208509) q[3];
u3(2.31951264686544,4.38931614713380,-0.0591503720784639) q[8];
cx q[8],q[3];
u1(0.682192028388333) q[3];
u3(-0.217759886162403,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.37766083219198,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.68736810026088,3.27951375391211,-2.34523864103580) q[3];
u3(2.34625295630624,-1.06530677403740,4.63564629527135) q[8];
u3(0.612923673170626,-0.581128641349869,0.936113379103187) q[0];
u3(0.175031412270712,-2.92778110439888,1.31276860751694) q[7];
cx q[7],q[0];
u1(1.35921768810196) q[0];
u3(-0.498430053328659,0.0,0.0) q[7];
cx q[0],q[7];
u3(2.92655913194777,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.68303954013347,-2.11373520754436,0.807999059842821) q[0];
u3(1.61403552491098,0.0777284325395324,-3.40350416971619) q[7];
u3(1.82128049890081,-0.390309068244900,2.53952882642940) q[9];
u3(2.14400083231008,-1.42958092300664,-1.45770690735239) q[4];
cx q[4],q[9];
u1(3.55108978612002) q[9];
u3(-3.94096271053690,0.0,0.0) q[4];
cx q[9],q[4];
u3(-1.15660867867199,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.21186152170687,1.48960496425644,-3.20913990432372) q[9];
u3(1.27897235682697,-0.279615866040913,-0.921730545439530) q[4];
u3(1.87766170637803,2.92432717943103,-1.64100006900441) q[6];
u3(1.80301464370193,1.89745733022105,-0.631138704059305) q[2];
cx q[2],q[6];
u1(1.59081500287016) q[6];
u3(-3.46021906474332,0.0,0.0) q[2];
cx q[6],q[2];
u3(2.20151437300416,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.71724218469469,-2.18348851553110,1.31899071659265) q[6];
u3(1.59652107702681,-0.408404606911280,-0.813315707455277) q[2];
u3(1.11629667080469,0.959642326578684,-4.05776621383186) q[5];
u3(1.92135302629029,4.25356486756684,-2.00570337617820) q[1];
cx q[1],q[5];
u1(1.54629014558977) q[5];
u3(-3.73250777399237,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.23672891311709,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.63718457167902,-2.35832471083611,0.134013986610489) q[5];
u3(0.434295698026114,-0.288813369660761,-0.154021879884976) q[1];
u3(2.25923404588845,1.61778177019599,0.931848207218812) q[7];
u3(1.39973654042080,-0.298022436800283,-3.11648431688333) q[1];
cx q[1],q[7];
u1(-0.250304474954839) q[7];
u3(-2.30196629241521,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.09656806968328,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.20422179099455,-4.29554022315079,1.83304612656561) q[7];
u3(0.856843551780055,-0.0136394104835509,-0.478273725125029) q[1];
u3(1.57683907252219,-0.282886304638987,2.06898417553091) q[8];
u3(1.79336940010744,-1.69069090042926,-0.723636929165231) q[6];
cx q[6],q[8];
u1(1.94966454794848) q[8];
u3(-2.60717427259535,0.0,0.0) q[6];
cx q[8],q[6];
u3(0.999163536442123,0.0,0.0) q[6];
cx q[6],q[8];
u3(2.05428158125989,-2.96941117437816,1.02151584812371) q[8];
u3(1.93762553500174,2.72122704446209,-1.78963127572153) q[6];
u3(1.27905357200638,0.830172964437379,-1.46160555636564) q[5];
u3(1.83964466361976,0.924513575073271,-4.99350740984240) q[9];
cx q[9],q[5];
u1(0.0482121422708193) q[5];
u3(-0.683408482713675,0.0,0.0) q[9];
cx q[5],q[9];
u3(2.99853802415147,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.47483649940900,1.55698074510850,-3.87637120289021) q[5];
u3(0.901051951989321,4.28217667767543,-0.399791289817649) q[9];
u3(2.37898359995804,-4.41480176345191,1.75359370730915) q[4];
u3(0.752720081084850,-1.87396258073646,3.46191270636815) q[3];
cx q[3],q[4];
u1(3.44500301647761) q[4];
u3(-4.34982454461283,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.162467165617613,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.486780726716321,4.24596645636862,-1.01683369934518) q[4];
u3(2.01601537741565,-0.884696916129442,3.62934164359263) q[3];
u3(1.11381220940050,1.67493490023500,-2.87373407579476) q[0];
u3(0.823331763286457,-2.95121824769402,3.06202228743294) q[2];
cx q[2],q[0];
u1(2.29371428099105) q[0];
u3(-1.73050874185740,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.24360375492374,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.79152760112170,0.120791439911067,0.550536333925239) q[0];
u3(1.57245889125611,-2.12454737091472,-0.782542176328294) q[2];
u3(0.655811999391107,-0.401135117699237,-0.562833902915625) q[0];
u3(1.22622606097408,-4.18009308568647,1.09684855099052) q[2];
cx q[2],q[0];
u1(-0.263116720728516) q[0];
u3(-1.90575827097951,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.803921019340685,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.16216585389523,-1.71253443069749,-1.40196940197567) q[0];
u3(1.22218095974632,2.23902563741703,-0.576993399424789) q[2];
u3(2.50763723700207,-0.432614941106902,1.27502660995850) q[4];
u3(1.84925577076059,-2.39647794801294,-1.83432998798555) q[5];
cx q[5],q[4];
u1(3.19279709869515) q[4];
u3(-1.83607843338747,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.54186562372568,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.80326979262916,-3.88190303280254,1.85424540175451) q[4];
u3(2.96161759931354,3.96328589051394,-0.136039461666792) q[5];
u3(2.38546643012436,-1.88952510304950,0.0356774515041611) q[3];
u3(2.11504433792245,-2.89994220721500,0.305438499495118) q[9];
cx q[9],q[3];
u1(0.991358367046725) q[3];
u3(-1.27178091602236,0.0,0.0) q[9];
cx q[3],q[9];
u3(-0.437176511862587,0.0,0.0) q[9];
cx q[9],q[3];
u3(2.37881332429811,1.01604209607370,-3.98721628248571) q[3];
u3(0.579480933073907,-5.41042092208191,0.125576031911245) q[9];
u3(1.15468723936821,0.938739938417488,1.39094130517896) q[1];
u3(1.15875203221428,-1.06254536231939,-2.24841497727469) q[8];
cx q[8],q[1];
u1(0.569531030309020) q[1];
u3(-1.21752087414805,0.0,0.0) q[8];
cx q[1],q[8];
u3(0.148893851168123,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.99946133228186,-3.52336464446056,2.12428303877370) q[1];
u3(0.802910029020350,3.31822839553679,1.44668306554209) q[8];
u3(1.92701494362858,-0.308300811847599,1.15034448532347) q[7];
u3(2.06945067830448,-0.447732869599616,-1.18615428392575) q[6];
cx q[6],q[7];
u1(0.969020726420808) q[7];
u3(0.0263451372828480,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.83351999557080,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.04122697560680,1.95310375550072,-2.99681750418010) q[7];
u3(0.319559928490952,-0.398906370025798,1.83509673399920) q[6];
u3(0.816375214592640,0.0650359540066332,-0.976384006858881) q[3];
u3(0.911407004265326,-2.76646177835017,1.18103921851004) q[0];
cx q[0],q[3];
u1(1.43199974414639) q[3];
u3(0.280489689868005,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.02138149005247,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.64754337159656,-0.0863842859431394,-3.12486213626696) q[3];
u3(2.04338547203356,-3.58066365005240,-0.784544524253626) q[0];
u3(1.88455357197578,2.54077201753758,-3.51110173381200) q[9];
u3(2.04642251858761,-3.62128825631341,2.61916470281052) q[2];
cx q[2],q[9];
u1(1.03310963537773) q[9];
u3(-1.59704235430390,0.0,0.0) q[2];
cx q[9],q[2];
u3(-0.635405959193816,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.703064493644319,-3.41409604130203,-0.515670228118525) q[9];
u3(1.34019961967793,-1.34434581594574,-2.08293503330691) q[2];
u3(1.01673601354437,-2.67814698155976,2.75624367395445) q[8];
u3(1.08735850640413,-2.96391272230317,2.14486507382618) q[7];
cx q[7],q[8];
u1(3.67592152520291) q[8];
u3(-1.21420989656541,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.64668263174888,0.0,0.0) q[7];
cx q[7],q[8];
u3(0.853045873365433,-3.19571686266733,2.15739381727366) q[8];
u3(2.29313108078203,2.43864321085459,-1.36834757548648) q[7];
u3(1.25648231834780,2.62516601897805,-0.261875068815117) q[1];
u3(1.95326261244947,0.0392148871533542,-4.04766779657109) q[4];
cx q[4],q[1];
u1(3.98929246448498) q[1];
u3(-3.30797437412953,0.0,0.0) q[4];
cx q[1],q[4];
u3(-0.516084193111565,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.976294369963260,-2.73102795351636,3.54001590650025) q[1];
u3(2.68842332625671,-1.42609437040686,-1.48909249173283) q[4];
u3(0.562054999067629,-2.11315289533874,2.87563314599064) q[5];
u3(0.504203660349532,1.33188056076514,-3.56415102150799) q[6];
cx q[6],q[5];
u1(1.85621654627917) q[5];
u3(-2.98433923763068,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.315241469045533,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.02271740244360,0.705718347825490,2.24336805667008) q[5];
u3(0.443428604111778,2.36207039599559,-2.74721567199934) q[6];
u3(1.10278835420355,-1.83587777009042,3.63400594751847) q[3];
u3(1.83606933527787,1.98019239312157,-2.17705885077439) q[7];
cx q[7],q[3];
u1(2.01933846066783) q[3];
u3(-2.73436958714389,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.02093756453526,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.43655589752484,-0.561414694118058,2.57166379899475) q[3];
u3(0.702179657668205,1.55437556041995,3.53922133266546) q[7];
u3(0.807426307600619,-1.63821537040722,-1.29979888946513) q[8];
u3(2.40549645930649,-2.38749995525042,2.63628857770285) q[5];
cx q[5],q[8];
u1(0.0638289906645519) q[8];
u3(-0.932727258146909,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.59311853173502,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.09725486142007,-0.227335063946863,0.528347377648260) q[8];
u3(0.191272577665972,-0.180296939930285,4.87338366048959) q[5];
u3(2.32506803141667,0.163560465805599,1.61734751223984) q[9];
u3(1.80492857320420,-2.22775101481886,-0.790259566895364) q[0];
cx q[0],q[9];
u1(2.46817501972303) q[9];
u3(-1.46905807766542,0.0,0.0) q[0];
cx q[9],q[0];
u3(3.21660714564841,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.03118746209004,-1.22769983246153,-0.256260379472782) q[9];
u3(1.90320944210906,3.90147081542452,2.12942638229389) q[0];
u3(2.30435769855404,4.04155456347813,-1.46027293320664) q[6];
u3(1.61217994387807,1.83651008670141,-0.653930431467156) q[4];
cx q[4],q[6];
u1(0.703346519838765) q[6];
u3(-1.27204866030800,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.17341854403774,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.09372183942965,-0.881688182190666,1.64046677111437) q[6];
u3(2.34453055935863,-2.71245624410474,-0.649666516009723) q[4];
u3(0.0998801232137059,-2.14593147478119,2.33040938138115) q[2];
u3(1.02148398763943,-3.51678059704618,1.13920122548939) q[1];
cx q[1],q[2];
u1(0.133986627556002) q[2];
u3(-1.72470641827027,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.508297403104500,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.965167476903550,3.26069370672895,-0.548298632046158) q[2];
u3(0.610829526280047,-1.49602271700315,1.10184766836774) q[1];
u3(2.15041759942115,-0.722207705952562,-0.953112314111919) q[7];
u3(1.70965623601566,-3.22801995836974,-0.230401052079128) q[9];
cx q[9],q[7];
u1(3.37511712254667) q[7];
u3(-0.897789820552680,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.39165683416867,0.0,0.0) q[9];
cx q[9],q[7];
u3(0.842062275313995,0.0192881081523705,-0.843587514750440) q[7];
u3(1.57575183421206,-2.13553418902326,-3.32851672207746) q[9];
u3(1.19654026813217,2.71134196906083,-3.42523780044875) q[4];
u3(1.68045925123453,-2.86232991718838,3.14387818000010) q[8];
cx q[8],q[4];
u1(2.24986219027973) q[4];
u3(-1.78873568717067,0.0,0.0) q[8];
cx q[4],q[8];
u3(2.89536339947066,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.70121593381051,0.885756464108377,-3.06592102372284) q[4];
u3(0.471565892746387,-2.23423413426871,-3.03618958525193) q[8];
u3(0.245511695310818,-0.490962114168306,0.261064211701437) q[5];
u3(0.273037054652420,-3.33036702778080,0.879722857365839) q[3];
cx q[3],q[5];
u1(0.935140259343216) q[5];
u3(-1.49224111975939,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.0229159363777909,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.77894653301224,2.70636564013700,-2.14074887375870) q[5];
u3(2.29649002070872,-4.73876046334846,0.786774949327923) q[3];
u3(2.09934704049331,0.928506484502341,-2.51302350772447) q[2];
u3(2.35079843799141,-3.46771737088564,2.58065025941720) q[6];
cx q[6],q[2];
u1(2.86051288797142) q[2];
u3(-2.54931084681133,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.30222596581067,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.31023877080388,-3.65467241585339,1.08417787755600) q[2];
u3(2.62777199556391,-1.64628542466133,-1.44853726596032) q[6];
u3(1.68703887380907,1.06413682027116,0.511326042935698) q[1];
u3(1.37621310306628,-1.78042157112220,-1.21973751190060) q[0];
cx q[0],q[1];
u1(2.77158576885290) q[1];
u3(-1.67835425014008,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.673068947125363,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.309627095314610,-0.0151540905081444,0.626124657818527) q[1];
u3(0.996006097393935,-2.36410838553774,-2.00211930849179) q[0];
u3(2.91227498723941,0.100285920492366,2.33166939650032) q[7];
u3(2.73939442068156,0.287306719492214,2.43051484207963) q[8];
cx q[8],q[7];
u1(2.45702064889472) q[7];
u3(-1.55152011041287,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.347890887829365,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.42884773937331,3.05039137560946,-1.12919013019369) q[7];
u3(1.01513385395014,2.27409192105701,-3.91423340111371) q[8];
u3(1.97197474668353,2.18286657827874,-3.95061795102653) q[1];
u3(0.876931733336515,1.85323599518785,-1.26314132229063) q[3];
cx q[3],q[1];
u1(1.60154161164018) q[1];
u3(-0.786936518506412,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.88238298642201,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.55989368844790,0.176650170351149,2.44287709539239) q[1];
u3(1.49592069727013,-1.16405463092851,4.88713909263546) q[3];
u3(1.65091574586644,1.11742060496630,-0.267137338107822) q[4];
u3(1.05399709323939,-0.501998698412002,-2.26604273750386) q[6];
cx q[6],q[4];
u1(3.29856989241431) q[4];
u3(-1.18362950660607,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.63266903038804,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.851976791516416,3.70712835946427,-0.0283456745809498) q[4];
u3(1.74531993466683,-1.94292875485842,-2.79910173409987) q[6];
u3(0.500718089721906,2.83045900678113,-2.57267663584982) q[5];
u3(0.656488878965772,-3.14845246310821,1.83783929746778) q[0];
cx q[0],q[5];
u1(0.0526628721677824) q[5];
u3(-1.01191649844549,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.47039908737479,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.279344201018399,-0.233240874126106,-0.413363403051342) q[5];
u3(1.36324882226579,2.58423488382807,-1.11756557866717) q[0];
u3(2.03034244460198,-1.45941788985558,-0.220304158623949) q[2];
u3(2.30352403880832,-3.06358335675341,-0.894392437727612) q[9];
cx q[9],q[2];
u1(2.23948415354641) q[2];
u3(-2.70219304952763,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.15494143068543,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.29409599391364,-1.42156451782424,2.41971569788884) q[2];
u3(1.82730221445062,1.35723194205229,-3.05113791875517) q[9];
u3(0.938090674909607,1.82296381568666,-2.57021846082855) q[1];
u3(0.972187822074187,-3.17782921721139,2.57212017862848) q[9];
cx q[9],q[1];
u1(2.05070701114984) q[1];
u3(0.705460652606061,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.61217596440080,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.07334623156765,3.88215737472315,-1.95066348678635) q[1];
u3(0.185079906989019,-1.81877683622603,3.82932503013231) q[9];
u3(0.953817697925301,0.491578383434489,0.977443599599506) q[6];
u3(1.40011540371785,-0.219020484801630,-3.70605316364902) q[7];
cx q[7],q[6];
u1(0.338976816464889) q[6];
u3(-1.85559778989194,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.12107881459367,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.89945425949327,4.44442612418851,-1.07765729846798) q[6];
u3(1.03320620724306,1.79938465260684,2.27261264268018) q[7];
u3(2.48885060816585,-1.29057011276787,4.38177481783546) q[0];
u3(0.603584709524095,-0.674545560411948,1.82579218532929) q[5];
cx q[5],q[0];
u1(4.42622990751541) q[0];
u3(-3.75293144469114,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.484745317125867,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.23992864045966,3.79551331656329,-2.08449840221610) q[0];
u3(0.491437947033578,0.603505456364257,-0.728224251674461) q[5];
u3(1.79393346129245,0.449328593021063,2.00480243918686) q[8];
u3(1.11644885645854,-2.71627251512564,-2.47650709095224) q[3];
cx q[3],q[8];
u1(1.32158903831369) q[8];
u3(-0.961601325753580,0.0,0.0) q[3];
cx q[8],q[3];
u3(-0.270972969596891,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.89961912883626,-1.24108983395582,2.14765583996706) q[8];
u3(1.32143239733112,2.10776020244916,-1.73292402561193) q[3];
u3(2.02104010954461,0.370943248028724,-0.957644200321451) q[2];
u3(1.63336208681200,-0.635935637885241,-3.41593071520296) q[4];
cx q[4],q[2];
u1(1.48709593898772) q[2];
u3(-2.74682424314078,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.937790195806707,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.685779130616661,-1.44942606306041,3.93123531164062) q[2];
u3(0.435798706099425,2.49234568427275,0.414080627816541) q[4];
u3(2.95436545094611,2.92613894923442,-2.17581456787468) q[9];
u3(1.55179411760910,-3.10937865854242,3.14245583887385) q[1];
cx q[1],q[9];
u1(0.0251059026211677) q[9];
u3(-0.850468141304386,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.13806703618027,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.39041583689455,-1.87146812917270,-0.238326179155963) q[9];
u3(2.20010889719094,4.80386535437583,0.669008414520047) q[1];
u3(2.54809290189562,-0.571300718437939,-0.942664136167457) q[8];
u3(0.551013760778779,-0.167533974724509,-4.13685898202082) q[4];
cx q[4],q[8];
u1(1.83601597438105) q[8];
u3(-2.55771553378537,0.0,0.0) q[4];
cx q[8],q[4];
u3(3.31763941004935,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.90179517170988,0.252916708451220,-3.15641754857204) q[8];
u3(2.99017553411470,-3.19150818446204,-2.33084040433296) q[4];
u3(2.93157464371224,2.30593328082321,-1.89045597580890) q[3];
u3(2.04800423268106,-0.427461819318332,-5.29888083470949) q[2];
cx q[2],q[3];
u1(3.39139166981232) q[3];
u3(-1.49757808327007,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.03412339038163,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.77169748109947,1.35730667511568,0.595367446286594) q[3];
u3(1.22479530834322,-2.46352769925844,1.30259684519492) q[2];
u3(0.730470586839162,0.453882516191745,0.697070274723661) q[6];
u3(0.792009475279778,-0.204998693897020,-1.74255033357536) q[0];
cx q[0],q[6];
u1(-0.683401934872990) q[6];
u3(1.41946544535966,0.0,0.0) q[0];
cx q[6],q[0];
u3(3.69681480103029,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.996198247901566,-1.79857535973116,4.17170544529682) q[6];
u3(2.35215611916604,-1.40474843650486,3.97217075755435) q[0];
u3(0.546684412347817,-1.91003998492941,2.17297876682800) q[7];
u3(0.327777233846190,1.70842977768836,-4.21224881072359) q[5];
cx q[5],q[7];
u1(2.44417839206857) q[7];
u3(-1.72450157063121,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.156106338325407,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.67442870253279,1.85320908324504,-2.72980889745487) q[7];
u3(1.72662033982504,0.154277286101328,2.49793236861004) q[5];
u3(2.46548440507341,2.42588391859550,-1.26873987466098) q[1];
u3(2.45993155743297,1.90270664552398,-3.09943775004535) q[6];
cx q[6],q[1];
u1(2.32673727982538) q[1];
u3(-2.88542447024655,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.35201009551990,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.870608970577298,-3.76642232592198,0.748340216894547) q[1];
u3(2.40239069106604,0.968300104145804,-2.95776186965251) q[6];
u3(2.36790915720645,1.70531419751258,0.405000796121922) q[7];
u3(1.18313685924363,0.385977739690198,-2.54010838097068) q[2];
cx q[2],q[7];
u1(3.29620679632554) q[7];
u3(-1.13443522764811,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.37416684857912,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.63461723388697,0.768122994004013,-0.733284717076616) q[7];
u3(0.676107085737970,3.22309729189813,2.86382463994949) q[2];
u3(0.765403489180194,-1.78911823476354,0.0578119629465110) q[3];
u3(0.958063041934684,-1.75518843803131,-0.280482439969204) q[5];
cx q[5],q[3];
u1(2.74184720840806) q[3];
u3(-2.98713806929623,0.0,0.0) q[5];
cx q[3],q[5];
u3(-1.15408678816763,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.96003561238396,-0.941883412893270,3.52455110557789) q[3];
u3(0.513441395860867,0.135621780987147,4.06187537098010) q[5];
u3(1.29456567068036,-0.830441097975304,2.80822996637528) q[8];
u3(0.478120726150588,-2.14285298641212,-1.67941892362728) q[4];
cx q[4],q[8];
u1(-0.0230763365188917) q[8];
u3(-1.34448414511222,0.0,0.0) q[4];
cx q[8],q[4];
u3(2.26367810943324,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.837241679725773,-3.43013213145111,-0.315793210261250) q[8];
u3(1.89080973822491,0.268574918541648,3.81844232926146) q[4];
u3(1.95153326556072,-0.562887870119426,1.38619021628284) q[0];
u3(1.64154410605057,-2.38274774948014,-1.80209033842147) q[9];
cx q[9],q[0];
u1(0.279417769871781) q[0];
u3(-1.46490283232628,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.73775243754579,0.0,0.0) q[9];
cx q[9],q[0];
u3(0.805410068305254,1.19155331804188,-1.78372408517140) q[0];
u3(1.56508440251521,-2.54801942584486,3.44887152021499) q[9];
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
