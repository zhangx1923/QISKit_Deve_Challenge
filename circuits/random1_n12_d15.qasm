OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.09874290606585,1.42737775518862,-3.49524570142763) q[1];
u3(1.33368502774978,-2.44683146249142,3.49363517971919) q[4];
cx q[4],q[1];
u1(1.47044850725492) q[1];
u3(-3.28435633394562,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.27789168341666,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.18616457909294,2.34594531714404,-2.69868822619772) q[1];
u3(1.91653646453273,0.985031148117917,-0.159020615530536) q[4];
u3(2.01744774125947,-2.78576514555218,1.05784719546975) q[8];
u3(2.69484876877076,-1.99665250864200,-0.557829376440781) q[6];
cx q[6],q[8];
u1(2.45483345582416) q[8];
u3(-1.62087411888461,0.0,0.0) q[6];
cx q[8],q[6];
u3(3.38871143395345,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.91661908997267,-0.734708688314001,1.36839118587433) q[8];
u3(2.92298142435044,-1.33422842957612,-4.64287177545876) q[6];
u3(1.57750856376160,0.665984891201546,1.72960117231132) q[0];
u3(2.14441593353543,-1.34931433168716,-0.982980545877026) q[11];
cx q[11],q[0];
u1(0.401969248713897) q[0];
u3(-0.920378992315581,0.0,0.0) q[11];
cx q[0],q[11];
u3(3.13956223808190,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.03483590659958,4.54344031795638,-0.907952015439367) q[0];
u3(1.59304520366624,3.44613981973954,-0.477386180492405) q[11];
u3(0.526340093531581,-0.603866143406947,-1.08941957356665) q[9];
u3(1.41404762802690,-3.71207282185411,1.80888681947022) q[10];
cx q[10],q[9];
u1(3.51478415293213) q[9];
u3(-1.49509152597924,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.28154516186426,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.46154659871447,0.315070425660803,1.32354452434689) q[9];
u3(1.39635558272562,-3.73111413246744,0.678918157436296) q[10];
u3(1.48232264401592,1.36619419564635,-0.580896186946253) q[2];
u3(1.53375027672900,0.383367177944826,-3.95311237891865) q[5];
cx q[5],q[2];
u1(0.168803793253681) q[2];
u3(-0.743711564334017,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.28280082260798,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.81251287878727,-1.08041563127037,-0.758682264527666) q[2];
u3(2.81562033531502,2.39696543952350,-2.84536962930780) q[5];
u3(1.29879425440855,0.966721018583196,0.756867493223734) q[7];
u3(0.847243320551757,-0.834554316861006,-1.95016424837425) q[3];
cx q[3],q[7];
u1(2.48432079994348) q[7];
u3(-2.10786030817619,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.490335892315913,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.28907607700260,-1.14680562839064,3.81617389833856) q[7];
u3(0.689681065732151,0.835136086342384,-2.94735411211655) q[3];
u3(1.51269872330417,1.00644374051606,-3.54143823454873) q[10];
u3(1.36947472471456,-2.23463132322598,3.17857744627277) q[0];
cx q[0],q[10];
u1(0.0506087785453013) q[10];
u3(-1.45986686566956,0.0,0.0) q[0];
cx q[10],q[0];
u3(0.963767123454378,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.276021790631274,-0.468158165802239,-0.815534353667571) q[10];
u3(2.05196325699352,0.368572109487600,0.713455560683261) q[0];
u3(1.11665619333475,3.06593794986463,-2.56295121770209) q[8];
u3(1.29975125131387,1.25983106898784,-1.11600122831073) q[2];
cx q[2],q[8];
u1(3.00702425783165) q[8];
u3(-1.94303418665658,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.598241966536174,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.26235224480020,-0.530604097055628,2.27371259554422) q[8];
u3(2.27247492967449,4.14647079233553,1.98238595602106) q[2];
u3(2.70838230834544,-0.551026134984728,2.02308374340950) q[11];
u3(1.69941263159284,-2.43878219617183,-1.07802517453760) q[1];
cx q[1],q[11];
u1(0.287282776127041) q[11];
u3(-1.48010778098304,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.37647576170251,0.0,0.0) q[1];
cx q[1],q[11];
u3(1.41709745600604,2.43378383673829,-3.50929414982722) q[11];
u3(0.968052572473033,1.26129577290785,-0.369882959532039) q[1];
u3(1.91707290026225,2.04592883640888,-0.533170213465936) q[7];
u3(2.93333194879447,3.99585239975433,0.106837325654380) q[6];
cx q[6],q[7];
u1(0.661343479281093) q[7];
u3(-1.56248149430376,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.94221056559366,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.14011113952062,1.38085132669166,-0.358881401325945) q[7];
u3(0.415124834394009,1.72252353362536,-3.33723173256346) q[6];
u3(1.96450503236216,2.96035522673773,-2.91072241562900) q[9];
u3(0.879109598633755,3.39693161306164,-1.49857677513359) q[3];
cx q[3],q[9];
u1(-0.264375388211187) q[9];
u3(-1.73834341059740,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.879713364889542,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.999781007625906,-1.49626417063483,3.74847479443173) q[9];
u3(0.996977253586107,-0.859601026424698,0.868842843083822) q[3];
u3(2.86282872828584,-3.29495737322125,0.760228436814328) q[4];
u3(2.36137055675130,1.11832259636110,2.35369241833021) q[5];
cx q[5],q[4];
u1(1.89494041532993) q[4];
u3(-2.07275192931347,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.0720881777721301,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.76362310108154,-1.84409419489199,4.01695085096405) q[4];
u3(0.658040284541566,0.274184439545206,0.603528654779168) q[5];
u3(2.83380012096048,-4.18002289524619,1.83232782498999) q[10];
u3(1.27554181145478,-0.0624602185992593,2.55349634611692) q[4];
cx q[4],q[10];
u1(3.06922890180331) q[10];
u3(-1.40681502853556,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.20069592139439,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.61241671457458,1.79238192713303,2.14006416486085) q[10];
u3(1.06310960765944,3.03712046784622,2.75289421971846) q[4];
u3(2.74694114259111,-2.35030514955831,1.94328279238014) q[9];
u3(2.00745149630776,1.14865436076235,3.04447142263651) q[5];
cx q[5],q[9];
u1(2.50462764700233) q[9];
u3(-2.60755119317342,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.77023944450265,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.88729781333734,-1.24912868351257,-0.506181609198881) q[9];
u3(1.60923975529423,-2.55246144714835,1.90342076332433) q[5];
u3(2.58694568345141,1.47539520939196,-2.23833117294905) q[2];
u3(2.48841545513163,0.211656828643002,-5.23295238657534) q[8];
cx q[8],q[2];
u1(3.44989647881877) q[2];
u3(-4.28427321848209,0.0,0.0) q[8];
cx q[2],q[8];
u3(-0.349750281178964,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.95782463442383,3.38221841034917,0.819564864218489) q[2];
u3(0.273463046022160,1.05080921927090,2.17734133758422) q[8];
u3(1.21923039716074,-1.37703293555990,-0.466448884675045) q[6];
u3(1.67392690188810,-3.33673200659151,-1.04575208034046) q[7];
cx q[7],q[6];
u1(4.21119644440966) q[6];
u3(-3.83475649769796,0.0,0.0) q[7];
cx q[6],q[7];
u3(-0.532912475929305,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.58013131039200,1.73109605638686,0.546441411292709) q[6];
u3(2.76568228382087,-3.01290094354080,0.829711734820015) q[7];
u3(1.22856634687679,-0.147889643467192,1.32763455078930) q[3];
u3(0.531166838671549,-2.09987133100652,-2.40022954999649) q[0];
cx q[0],q[3];
u1(2.93386876882176) q[3];
u3(-2.07682665820796,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.49897640114584,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.04777679269338,-0.204652355270505,1.46126511915384) q[3];
u3(2.05327157621539,-2.97710096955288,-2.10372744253257) q[0];
u3(2.09663857749551,1.02078609399247,-3.58270225694847) q[11];
u3(0.959306754752683,-0.977434231065543,5.18539248261612) q[1];
cx q[1],q[11];
u1(1.74838555537215) q[11];
u3(0.845925015999202,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.27276761306284,0.0,0.0) q[1];
cx q[1],q[11];
u3(2.62376678244235,-3.09778586315064,2.55868837124340) q[11];
u3(2.10115942076073,-1.64346058001529,1.40833679899235) q[1];
u3(2.51015210143157,-2.96376390708506,2.21521507223767) q[2];
u3(0.478585429771337,3.11595420342732,-1.20111354224657) q[3];
cx q[3],q[2];
u1(1.31181874593420) q[2];
u3(-0.618052844339512,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.76674951923774,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.54325659081078,0.252665364267085,1.92472082143697) q[2];
u3(2.47301814119653,0.435758597466944,-2.27303254516985) q[3];
u3(1.07726531993460,-2.28627663162315,1.99081856175132) q[1];
u3(0.754394783046910,1.63755658265116,-2.90503657508292) q[6];
cx q[6],q[1];
u1(1.43434109372732) q[1];
u3(-3.35999089499246,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.49061474736969,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.26897446324146,-0.652293093530074,-2.13431865651672) q[1];
u3(2.41038688809220,-2.11260661437763,2.68984741444512) q[6];
u3(2.10226996227820,-2.82067804157843,0.161460362667674) q[5];
u3(2.27682825935465,0.440556654438162,1.72806110379966) q[4];
cx q[4],q[5];
u1(1.12273174038565) q[5];
u3(-0.355656989372755,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.42379113146568,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.75501299291473,-2.71276389524960,0.548081762774770) q[5];
u3(1.12415357216732,-3.55019551810056,-2.60261570356116) q[4];
u3(2.15321795646723,1.38927603896719,1.39722421682670) q[8];
u3(0.621177219859396,-0.780563308402625,-3.28583350869319) q[10];
cx q[10],q[8];
u1(1.71138977103693) q[8];
u3(-2.16476407136347,0.0,0.0) q[10];
cx q[8],q[10];
u3(0.810152061323684,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.932291672874113,-2.99367269983384,3.02945337612503) q[8];
u3(1.48844876753879,-1.73496169863816,-2.17203626409618) q[10];
u3(1.39478428967424,1.23006025462025,-2.22334949127795) q[7];
u3(2.27696344574844,2.25768483789996,-3.72383511200279) q[9];
cx q[9],q[7];
u1(1.36345172512665) q[7];
u3(-0.601538998790793,0.0,0.0) q[9];
cx q[7],q[9];
u3(-0.00428383745419003,0.0,0.0) q[9];
cx q[9],q[7];
u3(2.35419279941934,2.31292261695907,-0.741415234396942) q[7];
u3(0.820838592453123,-0.556360173191977,5.34020170715901) q[9];
u3(1.75328095984662,-2.50000093600301,2.89486144802625) q[11];
u3(2.84689307073455,-2.68751907181960,1.25940841037466) q[0];
cx q[0],q[11];
u1(1.58301051686015) q[11];
u3(-2.48832719330198,0.0,0.0) q[0];
cx q[11],q[0];
u3(0.961072557465815,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.753547638853761,-2.37850360409437,3.40601587029121) q[11];
u3(2.39229405734639,-0.400821850450045,2.30390574351517) q[0];
u3(1.01583557044230,3.04328513011119,-2.59199380920309) q[9];
u3(0.427851865073227,1.69092378116409,-2.69164340675821) q[1];
cx q[1],q[9];
u1(2.67768961828606) q[9];
u3(-1.60321689937617,0.0,0.0) q[1];
cx q[9],q[1];
u3(0.797442083659478,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.50304136949359,2.20040514007313,-0.865869814892709) q[9];
u3(0.956543892690573,-1.58061680928333,2.21319070452934) q[1];
u3(2.30029410964270,1.05044015654384,-1.43019913048680) q[5];
u3(1.22731747081366,-4.88230981850556,1.11569660458834) q[11];
cx q[11],q[5];
u1(0.628139443532924) q[5];
u3(-1.39192056362471,0.0,0.0) q[11];
cx q[5],q[11];
u3(2.37705015427441,0.0,0.0) q[11];
cx q[11],q[5];
u3(0.0674860120061567,-4.91470418581165,1.20778972175565) q[5];
u3(1.66716822142056,3.17311731647658,-1.09318919395755) q[11];
u3(2.10898773567980,0.667783684641609,-0.617119264054195) q[10];
u3(1.67546692194669,-0.751750141359573,-3.37265630951727) q[4];
cx q[4],q[10];
u1(1.36171798021566) q[10];
u3(-3.30545614451769,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.29443435033248,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.772687867047909,0.104557497858650,-1.60645942317275) q[10];
u3(1.23896707578096,1.21622966794527,-3.56321874005726) q[4];
u3(1.54725170190318,0.232193047075548,2.81105713942598) q[2];
u3(1.87994844115277,-2.43190453909010,-1.66855953425578) q[3];
cx q[3],q[2];
u1(2.58972818001240) q[2];
u3(-1.67397186651045,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.250631204613386,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.47590518771577,2.50226107469769,-2.73419230047741) q[2];
u3(2.61554359010018,1.71479353418441,-3.75150702635896) q[3];
u3(2.31489370432146,-2.29162055752888,0.990533023641875) q[7];
u3(2.89469254194150,0.429394679850203,1.07988028465263) q[8];
cx q[8],q[7];
u1(1.28681361591469) q[7];
u3(-0.0743804003651332,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.48645915048257,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.33381110300912,0.0891684299806669,0.0379174014814211) q[7];
u3(1.57661624815093,-0.946677901878049,2.55734517131759) q[8];
u3(2.14039459014959,2.09951891392809,-3.50164640753088) q[0];
u3(2.26006494688845,2.96353886749083,-2.94454981326728) q[6];
cx q[6],q[0];
u1(2.41088532351447) q[0];
u3(-1.68601695468408,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.30278290720849,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.37504311345289,3.52987632579365,-2.33957003473656) q[0];
u3(2.14041757967802,-2.43010112844920,-3.09475607035373) q[6];
u3(2.41635289414246,4.17444750427654,-2.09653710062044) q[9];
u3(1.39605505981633,0.641356465055846,0.236827122248336) q[0];
cx q[0],q[9];
u1(3.45490831954700) q[9];
u3(-1.25431220804221,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.81133329710138,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.66948480220289,-2.39502649448950,2.37121147597101) q[9];
u3(1.99996148685159,4.21724427454073,0.240178722267601) q[0];
u3(0.542076218673517,-3.08445523929230,1.97560098780771) q[5];
u3(0.653810028141427,1.84276068047029,-3.21624812062973) q[10];
cx q[10],q[5];
u1(1.44801951604806) q[5];
u3(-3.02479851520679,0.0,0.0) q[10];
cx q[5],q[10];
u3(2.75387599984963,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.94965113400382,-3.05145258369242,1.20935874895664) q[5];
u3(2.26914808845981,0.560538076452699,4.07426522486323) q[10];
u3(1.82858802484117,0.646138308580723,0.724977731371650) q[6];
u3(0.805092602942029,-1.16385786242570,-1.81360919102936) q[4];
cx q[4],q[6];
u1(0.0468025031524253) q[6];
u3(-1.30715399014668,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.68323958528931,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.56704590818865,2.98716578369072,0.451169252622065) q[6];
u3(1.10725602546630,1.97998823365322,3.42869176454084) q[4];
u3(2.14893175248013,0.894228877962220,-2.36980652312733) q[7];
u3(2.70437831038213,3.93177621638177,-1.36675132262986) q[3];
cx q[3],q[7];
u1(2.63584458489327) q[7];
u3(-1.94777864351209,0.0,0.0) q[3];
cx q[7],q[3];
u3(0.937138194090729,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.69002939032209,0.0640453021417481,0.632904545665747) q[7];
u3(1.69713899514646,-1.84890930124329,-0.747729130479658) q[3];
u3(1.28804880054253,3.28116973783214,-0.755383163125520) q[2];
u3(1.40569074363779,1.24844077077134,-1.08766374826546) q[8];
cx q[8],q[2];
u1(-0.234030606574097) q[2];
u3(-1.43908006018140,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.90432034986856,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.52958553422571,-1.79908208975168,-1.38992917798796) q[2];
u3(2.00507604010144,-3.41646943257161,2.21547728778566) q[8];
u3(1.05385891977696,0.0637788444064310,-0.862474702310515) q[1];
u3(2.10562230384979,-3.69951030570982,1.31251703473860) q[11];
cx q[11],q[1];
u1(2.41805731828550) q[1];
u3(-1.95325289780101,0.0,0.0) q[11];
cx q[1],q[11];
u3(0.0569540508035224,0.0,0.0) q[11];
cx q[11],q[1];
u3(1.19602648408495,-1.18146214271738,0.707468601270561) q[1];
u3(1.73638734936221,-5.72631994629868,0.439411905757336) q[11];
u3(0.399937782577486,-2.66109237535083,1.57336555955917) q[11];
u3(1.08142305265086,1.56112751049273,-3.17730333173144) q[2];
cx q[2],q[11];
u1(1.60844865859381) q[11];
u3(-2.89209477794681,0.0,0.0) q[2];
cx q[11],q[2];
u3(0.670340752281123,0.0,0.0) q[2];
cx q[2],q[11];
u3(2.08395871653233,-1.87370128644375,0.755530637054200) q[11];
u3(1.06506373834074,-4.29093815866024,-0.314455446644006) q[2];
u3(1.61289238247498,0.984581378979863,-3.50858331369477) q[8];
u3(2.77700630285795,2.72424179642509,-2.91757491880055) q[5];
cx q[5],q[8];
u1(3.46570282650092) q[8];
u3(-1.46980421878037,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.07703303959095,0.0,0.0) q[5];
cx q[5],q[8];
u3(0.727861316068666,-1.28420637587997,-1.50963669824646) q[8];
u3(2.15829171029118,3.55313933302491,-0.859906346930063) q[5];
u3(1.26338349450142,1.12712644001783,-1.08711053890688) q[9];
u3(0.225572639723258,1.41935185761386,-3.71793736273452) q[10];
cx q[10],q[9];
u1(0.0374338040472086) q[9];
u3(-1.48509637366486,0.0,0.0) q[10];
cx q[9],q[10];
u3(2.42074265971922,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.70855565163176,-1.22166100508815,-1.69305349181932) q[9];
u3(2.11573303554273,4.23486091538555,-0.916608076570645) q[10];
u3(0.926727897062005,-2.46258003726408,1.94538673320980) q[1];
u3(0.963950332263112,0.836494484751050,-2.70329554913523) q[6];
cx q[6],q[1];
u1(3.45441738765694) q[1];
u3(-0.967407556622396,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.67488591124092,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.13628155202521,-0.399476682166797,2.64140153037897) q[1];
u3(0.947956747378833,-1.31859532805913,-1.06916281763915) q[6];
u3(1.29520172105701,0.966494336774829,0.899468450116833) q[0];
u3(1.46471498668321,-0.208725824554885,-2.58140134124873) q[4];
cx q[4],q[0];
u1(0.828598388186891) q[0];
u3(-0.224584657019188,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.58844909254834,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.806262916764612,4.20273446486646,-2.04716458228168) q[0];
u3(2.59479862332712,-0.979457271162695,-3.83876269932267) q[4];
u3(1.59050267870735,3.83066718924358,-1.78486789213790) q[7];
u3(2.73598138976839,1.72617811543763,-2.00787901769165) q[3];
cx q[3],q[7];
u1(2.01008389391229) q[7];
u3(-1.80021450894464,0.0,0.0) q[3];
cx q[7],q[3];
u3(-0.0385629810267478,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.861799643486412,0.780137319225103,-2.55476404414409) q[7];
u3(1.42656605568639,0.0590661086833417,-2.47259259735095) q[3];
u3(2.62084127232780,-2.36531379580531,3.56014931679948) q[0];
u3(0.654367568095044,-0.945321340456445,2.98472714852164) q[1];
cx q[1],q[0];
u1(1.33049406163774) q[0];
u3(-1.17561612973698,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.264048388559870,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.127832076716815,-4.38490225933213,0.969204767647247) q[0];
u3(2.10669997758989,0.658628784019230,-2.44381983529021) q[1];
u3(2.70030024725840,3.02030762622084,-2.68726872402244) q[3];
u3(1.79355250389092,2.81107299987902,-2.90691528388432) q[6];
cx q[6],q[3];
u1(0.120016982703599) q[3];
u3(-1.93255685927126,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.464784885904038,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.321619582307945,1.88478176337218,-0.120874508951443) q[3];
u3(0.962178985023604,3.08448429747812,1.83618832971910) q[6];
u3(1.69598121168560,-0.567701295105865,1.91552030225275) q[9];
u3(1.41246480143961,-2.49655103176212,-1.39424976456452) q[4];
cx q[4],q[9];
u1(2.99439906958241) q[9];
u3(-4.37577311967257,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.245224164908103,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.34183786076901,1.72147171275624,-3.42772843637298) q[9];
u3(2.28837508176839,-3.18140762755707,0.947153696800286) q[4];
u3(1.25142191033158,1.66169965526443,-2.61752791571491) q[11];
u3(2.01996418982026,-2.05935699624081,2.85547335773789) q[8];
cx q[8],q[11];
u1(3.30528083252150) q[11];
u3(-1.41984989522626,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.82650659528687,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.88608007238413,0.369390423665245,-4.71950345672109) q[11];
u3(1.90637432855939,-0.616361453952204,-2.85095221936903) q[8];
u3(2.47818325679443,2.92841645350898,-2.40621322956383) q[10];
u3(1.72783054111105,1.75951600386014,-2.12543252112709) q[5];
cx q[5],q[10];
u1(-0.376886864712265) q[10];
u3(-2.42822162651044,0.0,0.0) q[5];
cx q[10],q[5];
u3(1.35003334348027,0.0,0.0) q[5];
cx q[5],q[10];
u3(1.56473604948482,-2.78864743218724,-0.302001060923567) q[10];
u3(0.889724348909202,-2.62722914795080,-0.911716185767855) q[5];
u3(1.32043539983218,-0.668412509476932,1.60189068973100) q[2];
u3(0.374536546962065,-0.409083267012694,-0.944195797423314) q[7];
cx q[7],q[2];
u1(1.74251598486946) q[2];
u3(-3.13358948539531,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.777002042683306,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.662170267305632,2.45240126940850,0.459262747584868) q[2];
u3(0.965092632041931,4.37002626880070,-0.463856832633073) q[7];
u3(0.783521095463464,-1.07388008644358,3.96702686743407) q[9];
u3(2.04671988505379,1.67848554012866,1.45719717115733) q[7];
cx q[7],q[9];
u1(2.34998206780594) q[9];
u3(-3.04750923577209,0.0,0.0) q[7];
cx q[9],q[7];
u3(1.11057157245469,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.06269222825492,1.13023910547569,0.411978375047492) q[9];
u3(0.590519679207827,-2.09995502463115,-3.11420703333376) q[7];
u3(2.16705145982552,0.467324090803038,2.10164437840629) q[1];
u3(1.31123850558686,2.56047370235007,3.02683776626778) q[6];
cx q[6],q[1];
u1(1.54761398913107) q[1];
u3(-2.95186515055603,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.609165459241173,0.0,0.0) q[6];
cx q[6],q[1];
u3(0.721846709061089,1.28389056114218,-3.06912453836039) q[1];
u3(1.98381866385441,0.211914655135038,3.74966303209562) q[6];
u3(2.41055706473844,3.47733487227849,-2.56142227734693) q[5];
u3(0.748237592443360,0.0188910593775904,1.07634169820764) q[4];
cx q[4],q[5];
u1(0.239889989598619) q[5];
u3(-1.29407817049834,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.26647808886794,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.52247715023030,-2.10898915125655,0.867282045588065) q[5];
u3(0.653392835497479,-1.49907580512181,-0.750150302355890) q[4];
u3(1.63153939245985,1.21971125266199,-0.380837399665631) q[8];
u3(1.30889287229352,1.36938826133513,-4.52809157468579) q[10];
cx q[10],q[8];
u1(0.887428099698540) q[8];
u3(-1.15711416433047,0.0,0.0) q[10];
cx q[8],q[10];
u3(2.40024864327032,0.0,0.0) q[10];
cx q[10],q[8];
u3(2.17642201754655,1.04298188258378,-3.31129368932882) q[8];
u3(0.759043972373825,-0.0231994898587275,-0.570340851854358) q[10];
u3(0.888403945974739,0.979170216892071,-0.115231938367880) q[11];
u3(0.761406984006845,-0.619604647684100,-1.79855003978213) q[2];
cx q[2],q[11];
u1(3.47301495219801) q[11];
u3(-1.55707841899358,0.0,0.0) q[2];
cx q[11],q[2];
u3(2.14010144103054,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.44010273136361,3.43081017468438,-2.27200919090932) q[11];
u3(1.16833462380539,0.790206902815405,-2.53457042853188) q[2];
u3(0.641779471201071,2.92979594176905,-2.57206774691323) q[0];
u3(1.21590733341116,0.686324746113271,-1.30695168307277) q[3];
cx q[3],q[0];
u1(3.15664078293521) q[0];
u3(-1.32355384648656,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.16865021913186,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.27235818054156,1.14009031491118,-4.13516132948257) q[0];
u3(1.53021105844507,1.16764296743398,3.19252413455588) q[3];
u3(0.682144691171190,0.0509449527267538,0.797962835378163) q[11];
u3(1.14419596457622,-0.816404567353841,-2.05746391911436) q[7];
cx q[7],q[11];
u1(2.75791876352581) q[11];
u3(-1.75756307671140,0.0,0.0) q[7];
cx q[11],q[7];
u3(-0.00172029963845644,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.14451729159484,-1.26046129561424,4.88178228471389) q[11];
u3(2.55741653838156,-0.886966508783516,2.37625820601240) q[7];
u3(1.27157424067232,-0.821291511042798,-0.666065818827101) q[4];
u3(1.40851920998569,-2.42097986145168,0.178855149664437) q[8];
cx q[8],q[4];
u1(0.794280226340267) q[4];
u3(-1.37994852409097,0.0,0.0) q[8];
cx q[4],q[8];
u3(3.03428430593435,0.0,0.0) q[8];
cx q[8],q[4];
u3(2.08756683388611,-4.35924503956783,0.409047817440176) q[4];
u3(1.18331711893392,-2.14767618505161,-1.73078021433763) q[8];
u3(2.54712572301945,-3.71800509893540,2.34904913957866) q[3];
u3(0.968296628402109,3.43735160552484,-2.46352714214938) q[1];
cx q[1],q[3];
u1(2.74196649355577) q[3];
u3(-1.70057036784307,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.565666601341122,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.48753775192841,-2.20772881088355,3.08303054185921) q[3];
u3(1.92636877548261,-1.50372221500494,2.43003660342418) q[1];
u3(1.30210927621403,-1.11163109051210,-0.751176053977062) q[0];
u3(2.62485107937850,0.732072869081569,-4.95692891842248) q[2];
cx q[2],q[0];
u1(1.53772778390710) q[0];
u3(-2.95370102065415,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.823500914385830,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.13044440586053,-0.257063204376270,1.40650995760818) q[0];
u3(1.55342135462908,0.876413631149659,-2.21874719891001) q[2];
u3(0.672385521293736,1.37106951608102,0.291775178391497) q[6];
u3(1.46932562612251,-0.586239723490739,-3.58157565995670) q[5];
cx q[5],q[6];
u1(1.61396629506582) q[6];
u3(-2.74220556600380,0.0,0.0) q[5];
cx q[6],q[5];
u3(3.29676398867066,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.100936013146275,2.19287543766203,-0.184251842279292) q[6];
u3(1.59964879246739,-3.95314370662373,2.15793422678565) q[5];
u3(2.41832546160529,-0.500631492438176,-0.926176540747929) q[9];
u3(0.720617429802931,0.498223345523832,-4.83696489100618) q[10];
cx q[10],q[9];
u1(0.244087914403525) q[9];
u3(-1.00053046928153,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.56869964948101,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.60693424972961,-1.67522438508211,3.99919612582124) q[9];
u3(2.52211936122657,0.974609224034324,-4.13683632306032) q[10];
u3(1.42188481414300,0.966231709152809,-3.72597677653535) q[8];
u3(1.47786199308202,-1.32422195634660,4.76354966267376) q[2];
cx q[2],q[8];
u1(2.65249963781282) q[8];
u3(-1.66436737415161,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.19761986372739,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.63173446638081,-0.983408642827260,-0.342441101903637) q[8];
u3(2.13177240793142,-0.596957219139636,2.66408029727745) q[2];
u3(0.826443844325176,-2.70872387636727,3.10163183421368) q[6];
u3(0.300808642630273,-2.83001968130549,1.72095576340861) q[11];
cx q[11],q[6];
u1(1.14031389754363) q[6];
u3(-2.89775951458193,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.94100374145326,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.87378859515846,0.717520574747717,-2.90244090087689) q[6];
u3(0.0555844643069161,-0.564961204697975,-1.64475010181317) q[11];
u3(2.52048900450125,1.10082043040687,-2.42969624574626) q[3];
u3(2.68123589415855,3.48263659868296,-2.26065196429618) q[1];
cx q[1],q[3];
u1(1.37274399976940) q[3];
u3(-3.60668606316382,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.65833362281327,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.00557117158390,2.06589686357163,1.16882752654419) q[3];
u3(2.44982002562613,-2.00736398991488,1.34338915050603) q[1];
u3(0.819006814477285,1.01373065675587,0.920563612922578) q[9];
u3(0.870432043183150,-0.426976092368003,-3.19897490699331) q[5];
cx q[5],q[9];
u1(3.42073779839587) q[9];
u3(-4.50362453157005,0.0,0.0) q[5];
cx q[9],q[5];
u3(-0.487069515639930,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.77201163019881,-1.75064588574193,2.69560034469644) q[9];
u3(0.925219192418984,-1.17110518933893,0.138186563218192) q[5];
u3(2.30991290976228,-0.757430405054570,1.02806030692156) q[0];
u3(1.50177009964744,-2.53038164893023,0.282130539802958) q[4];
cx q[4],q[0];
u1(0.572223925201052) q[0];
u3(-1.33834763768301,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.79943973833537,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.50016941073107,2.86514943788200,-1.31336669862656) q[0];
u3(1.63232909185427,-0.0333300636722038,0.682415135696866) q[4];
u3(1.56202473590694,3.64448500046809,-2.21249923003426) q[7];
u3(2.08636875680440,1.66523078919174,-1.93696216198491) q[10];
cx q[10],q[7];
u1(2.07774898397247) q[7];
u3(0.0871938051007541,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.742502156816645,0.0,0.0) q[10];
cx q[10],q[7];
u3(1.50339815520930,1.97650453036074,-3.97174188427607) q[7];
u3(0.527719402297854,-1.05523865270345,3.95915986782448) q[10];
u3(2.52621721053434,0.0687741562532322,-1.81515694180588) q[10];
u3(1.30891459763591,1.27009217120089,-3.74784456472515) q[6];
cx q[6],q[10];
u1(0.747707427840062) q[10];
u3(-3.27512078470173,0.0,0.0) q[6];
cx q[10],q[6];
u3(1.84514847253845,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.13894329909828,0.371563521702471,-1.81045000819398) q[10];
u3(1.01487560594355,-1.45635944655244,0.992771125928429) q[6];
u3(1.82977116239672,2.56407538104016,-1.93489636769127) q[3];
u3(1.39832345753077,-2.86942753827839,2.70837134814842) q[8];
cx q[8],q[3];
u1(1.48471263189032) q[3];
u3(0.0952309627370704,0.0,0.0) q[8];
cx q[3],q[8];
u3(2.65036126715993,0.0,0.0) q[8];
cx q[8],q[3];
u3(1.72892941244476,0.320004249628958,0.307224352568055) q[3];
u3(1.12589585240872,0.324330467757059,0.439901604030045) q[8];
u3(1.21221078487476,0.606990353098295,-0.704448246402268) q[5];
u3(0.282740787033403,-4.15425765604500,1.53128882151350) q[1];
cx q[1],q[5];
u1(2.62276720231008) q[5];
u3(-1.97302627046855,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.986820890008013,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.97894763009721,3.26778461013250,-2.97591644228875) q[5];
u3(0.162669409897657,-4.10417270099359,-0.258018479693200) q[1];
u3(1.00376938490307,2.62491200240607,-2.25995375466288) q[4];
u3(0.385526559023266,-3.30924049681274,2.75041826320146) q[0];
cx q[0],q[4];
u1(-0.363633630487707) q[4];
u3(-1.89390452142779,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.56124424165900,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.44929527280452,1.01531531463177,1.43042037144059) q[4];
u3(1.30707924972853,2.95751021886348,-2.57346187891753) q[0];
u3(1.02688490799289,-1.26780536246492,0.496583861263631) q[11];
u3(0.898614404148207,-1.66474550006611,0.769784491544398) q[9];
cx q[9],q[11];
u1(2.93167154102513) q[11];
u3(-2.62235237721724,0.0,0.0) q[9];
cx q[11],q[9];
u3(0.903513469260581,0.0,0.0) q[9];
cx q[9],q[11];
u3(2.57906871810984,1.86337748201645,-2.53632167876698) q[11];
u3(1.61047099776654,1.22712801986070,-3.63896277874528) q[9];
u3(0.818170972867754,1.95348146403163,-2.29372203320740) q[7];
u3(0.757366719855118,2.13423254506589,-3.92911162150915) q[2];
cx q[2],q[7];
u1(3.49771635108516) q[7];
u3(-0.651457171264280,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.62251263908212,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.30265849835346,-2.46076331816472,-0.265765141646099) q[7];
u3(0.668028757090911,2.80934072914856,1.67116450712944) q[2];
u3(1.80607477157603,-0.167753413506556,2.41299026338032) q[4];
u3(1.92638290807171,-2.27851210160013,-1.72445124691538) q[3];
cx q[3],q[4];
u1(2.85334940302333) q[4];
u3(-1.84277924318813,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.46620615208211,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.20826011945179,2.73193958380057,-2.62162585093727) q[4];
u3(0.690520058372020,0.350531289694401,4.94541365394752) q[3];
u3(1.29687824865541,-2.46289820833110,0.273289433010788) q[7];
u3(1.27727791715631,-3.61285276010289,-0.176821035877767) q[11];
cx q[11],q[7];
u1(2.50873995365625) q[7];
u3(-2.77341186841901,0.0,0.0) q[11];
cx q[7],q[11];
u3(0.779685010983908,0.0,0.0) q[11];
cx q[11],q[7];
u3(2.41765146363731,-0.318035015455201,-1.13420829681726) q[7];
u3(1.63769237042230,-0.420267567086422,3.33831632768082) q[11];
u3(2.32556370391420,-1.29701853740945,-1.07269869724107) q[9];
u3(0.367066795387259,0.776499396836392,-5.29497413372325) q[10];
cx q[10],q[9];
u1(2.99855337981006) q[9];
u3(-2.39033112688815,0.0,0.0) q[10];
cx q[9],q[10];
u3(1.13786457169085,0.0,0.0) q[10];
cx q[10],q[9];
u3(0.358447962564295,2.51997138406267,-2.98941375702157) q[9];
u3(0.688764261111430,2.92178670071239,1.31186145903962) q[10];
u3(0.839308501252441,-1.61251690497569,2.56091231310807) q[1];
u3(0.175155106278763,-1.44355344748342,-0.175442759282957) q[5];
cx q[5],q[1];
u1(1.52328086073049) q[1];
u3(-4.20123192793336,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.78930737543507,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.43490892720803,1.95144651666372,0.746884505741326) q[1];
u3(0.873662015034433,-1.81354210936464,2.34840381543808) q[5];
u3(1.14151834190065,1.98236161768184,-1.30937933114188) q[6];
u3(0.871971146495782,-2.10709098523093,0.295965573765746) q[0];
cx q[0],q[6];
u1(0.541150855192393) q[6];
u3(0.0224964357225663,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.15285135124571,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.23926618403445,1.79696765016766,-1.45201435123806) q[6];
u3(1.90856233507640,-2.05948494524008,-3.68436607183360) q[0];
u3(1.97248041479437,-1.01961038017519,0.956812704116659) q[2];
u3(2.55198852052803,-0.767150168840203,-1.38909622710590) q[8];
cx q[8],q[2];
u1(-1.12434881817592) q[2];
u3(0.403721903342420,0.0,0.0) q[8];
cx q[2],q[8];
u3(3.43199734736107,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.34034059951890,2.65344887194343,-0.257516440472870) q[2];
u3(2.40082981924105,1.37268216062366,-2.59338273670395) q[8];
u3(1.04236935860389,-0.687121781806527,1.97771891372104) q[8];
u3(0.885201565296423,-0.750583693545350,-0.0543142297682041) q[5];
cx q[5],q[8];
u1(3.44622986807837) q[8];
u3(-4.18637145625772,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.737282076357521,0.0,0.0) q[5];
cx q[5],q[8];
u3(2.46671903595033,-0.188543269005112,4.51472495469377) q[8];
u3(1.57657645113832,2.57084523737921,1.23349775571446) q[5];
u3(0.227713842714437,2.13717172615157,-2.19869951921012) q[4];
u3(0.576378645816540,-0.835023145371476,-1.64185573103161) q[7];
cx q[7],q[4];
u1(1.04418667226219) q[4];
u3(-0.262626861302127,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.59567728101368,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.08945277909377,1.65263674316729,-0.299055671170849) q[4];
u3(1.74921537966005,-3.50047629007562,-0.119829002197249) q[7];
u3(1.47200434699587,2.14528210873478,0.183771742440193) q[3];
u3(0.732308101246018,-0.978133869772496,-1.81243654117934) q[9];
cx q[9],q[3];
u1(-0.336104351096776) q[3];
u3(1.19115783834220,0.0,0.0) q[9];
cx q[3],q[9];
u3(3.71127945588848,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.65436537183081,-1.89880907539796,2.39736796403123) q[3];
u3(0.360230280528061,-0.941332700379824,2.58374586041920) q[9];
u3(1.84631413855427,1.73519456586407,-0.906615943807803) q[10];
u3(2.36283485887265,-0.613098438991485,-3.39154624830178) q[11];
cx q[11],q[10];
u1(0.913429856327830) q[10];
u3(-3.32953389814525,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.68685540606767,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.71061545969891,3.99313918797888,-1.91512639167487) q[10];
u3(2.30639824600391,1.28015215165955,-3.26794227317445) q[11];
u3(2.36041609741622,-0.466683016631791,2.19205084220984) q[0];
u3(1.87979678469097,-2.63081310451110,-2.19405834231235) q[6];
cx q[6],q[0];
u1(2.75041942458251) q[0];
u3(-1.57720121641963,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.499831996533594,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.18129483984799,0.709703630270705,-1.90783714763712) q[0];
u3(2.44979685462680,2.56567347380575,1.64824957076224) q[6];
u3(2.31522246019595,0.301032328578016,0.573843768040552) q[2];
u3(1.57748472547264,-1.98911464846121,-1.00682297431964) q[1];
cx q[1],q[2];
u1(3.46052029833377) q[2];
u3(-1.97442797940604,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.07506662385995,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.167033160574892,0.441354666858814,-0.626121798870269) q[2];
u3(1.36981400787481,2.47411902811818,-2.15264867455744) q[1];
u3(0.182316406421507,-2.18476798508508,2.13223037474900) q[6];
u3(0.843717791729021,0.0313280976880214,-1.65322202851846) q[2];
cx q[2],q[6];
u1(1.64751821636120) q[6];
u3(-0.631900217464607,0.0,0.0) q[2];
cx q[6],q[2];
u3(-0.0477027877609881,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.399099146640934,1.19418745202061,1.70322850000145) q[6];
u3(2.15650443108522,3.87783731263888,2.38081931047320) q[2];
u3(0.966471980305762,-0.232481519464031,1.88277351617365) q[8];
u3(1.08436055820662,-1.55384029765710,-2.51473683005138) q[11];
cx q[11],q[8];
u1(1.06506196740458) q[8];
u3(-0.921137303538581,0.0,0.0) q[11];
cx q[8],q[11];
u3(-0.0480166366729218,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.43586647029261,2.32998280439628,-3.49937086345252) q[8];
u3(0.950343108840368,1.00468216519618,-5.24086794578142) q[11];
u3(0.276385521038569,0.942179535235549,-1.96313523998818) q[9];
u3(0.866278111767496,0.299282394451174,-2.06247147949302) q[7];
cx q[7],q[9];
u1(1.40104946300419) q[9];
u3(-3.44286834657549,0.0,0.0) q[7];
cx q[9],q[7];
u3(2.20799598320662,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.51732142918724,-0.884449066915603,1.29503316874871) q[9];
u3(2.03144276802850,1.94180263364801,-2.36788295482419) q[7];
u3(0.717954479804975,1.39567289478803,-1.93294787380234) q[0];
u3(0.0972774983563109,1.76316201380607,-4.11491603808639) q[4];
cx q[4],q[0];
u1(1.31284218124061) q[0];
u3(-0.619743083869543,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.16693369770769,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.31310631896797,-3.38513663575099,0.0331440677816428) q[0];
u3(1.18515494850673,-0.604455470589385,0.880906999894200) q[4];
u3(1.81470508802373,0.129356779753320,-1.53590345903522) q[1];
u3(0.782200807404236,0.364538596852779,-4.59048518209174) q[10];
cx q[10],q[1];
u1(-0.105841322796089) q[1];
u3(-2.25923107346777,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.10623021186585,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.75863405006718,2.14615195691939,0.0994203804641984) q[1];
u3(2.22700968374867,1.07496630327688,-2.92784109158716) q[10];
u3(0.644713367575123,-1.41746982373449,2.43363444754724) q[5];
u3(0.464131928504073,-2.36771558264423,1.10715809581363) q[3];
cx q[3],q[5];
u1(3.54310722887901) q[5];
u3(-1.59249659276566,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.82093737587462,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.87670424884120,1.28164556474265,-3.01189704672456) q[5];
u3(2.00328788757876,5.57012527120219,0.302967082046565) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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