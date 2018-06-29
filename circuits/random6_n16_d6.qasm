OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(1.12479538831618,0.663848267020873,-0.222511101742944) q[1];
u3(1.29022822707382,0.0990722205160666,-3.13440079392964) q[13];
cx q[13],q[1];
u1(2.95031240671244) q[1];
u3(-1.83818795178890,0.0,0.0) q[13];
cx q[1],q[13];
u3(1.16577357407157,0.0,0.0) q[13];
cx q[13],q[1];
u3(1.45136107175987,0.465736339610794,1.70740208809617) q[1];
u3(1.22694762602718,3.41132798369384,-2.52341335102224) q[13];
u3(2.16616462807070,0.127942443629400,2.55508786932002) q[7];
u3(2.57866683071983,-2.01302989000707,-0.594330746340228) q[4];
cx q[4],q[7];
u1(1.82002594182082) q[7];
u3(0.227826311077733,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.19657273477433,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.39748787661498,1.27252989594147,-1.93341048767896) q[7];
u3(1.40422701452581,1.97183597185690,-2.32500995394824) q[4];
u3(0.737641800336966,2.69493919522712,-1.26540556476402) q[11];
u3(1.71679911339492,0.175174319953628,-3.37677782888258) q[12];
cx q[12],q[11];
u1(2.43448700312133) q[11];
u3(-1.91990760003562,0.0,0.0) q[12];
cx q[11],q[12];
u3(0.274051388285716,0.0,0.0) q[12];
cx q[12],q[11];
u3(2.03139698273335,3.66859941747449,-2.18814527339350) q[11];
u3(0.227602413247144,-1.26379174880427,3.77615927641322) q[12];
u3(1.45149157299571,3.14393343917272,-2.44924981100193) q[14];
u3(2.21435310790479,1.13046654231642,-1.82408692583073) q[2];
cx q[2],q[14];
u1(0.689444701307913) q[14];
u3(-0.499572764723311,0.0,0.0) q[2];
cx q[14],q[2];
u3(2.90398305963579,0.0,0.0) q[2];
cx q[2],q[14];
u3(1.80350328481205,-2.23121659902585,-1.92368224244980) q[14];
u3(2.35133034025920,-2.19680167578444,-3.58170633732897) q[2];
u3(1.03437191249600,0.0984544336227751,-1.85702855980246) q[5];
u3(1.31959158631114,-5.40483685929228,0.802502354405759) q[6];
cx q[6],q[5];
u1(-0.361008413977526) q[5];
u3(-1.46968068894527,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.921787905897759,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.331122352432409,-0.358348242106617,0.248976213647096) q[5];
u3(1.27167129891570,3.24075835103115,-1.55874929314544) q[6];
u3(1.56567116586301,-1.44960112018989,-0.000595011088056596) q[15];
u3(0.235213857028996,-2.68051291209291,0.368106502625733) q[9];
cx q[9],q[15];
u1(2.04076273161399) q[15];
u3(0.425732027821290,0.0,0.0) q[9];
cx q[15],q[9];
u3(1.30791981798276,0.0,0.0) q[9];
cx q[9],q[15];
u3(2.41975751148030,1.15667310139081,-3.49275880905102) q[15];
u3(1.29499167594138,1.38781332082427,-4.54787228509545) q[9];
u3(2.31705954002355,0.166245023886097,2.66436781007998) q[0];
u3(2.33123461621958,-2.54838881005423,-1.98312411021549) q[3];
cx q[3],q[0];
u1(3.32550527394192) q[0];
u3(-0.825813152018084,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.13639096715608,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.960447306677710,0.893356952253452,-3.14462451078879) q[0];
u3(1.78267631549349,5.46057338603318,-0.0808246189354147) q[3];
u3(1.09551526814338,0.171319951142083,-1.55082279109935) q[10];
u3(2.07749292711411,0.541899973026534,-5.53114914363061) q[8];
cx q[8],q[10];
u1(0.266853464309048) q[10];
u3(-1.32908005438184,0.0,0.0) q[8];
cx q[10],q[8];
u3(2.21580684826620,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.34940492931100,1.59608646342923,-3.28415360562906) q[10];
u3(1.24680146325503,-1.50897175996656,3.19869890177118) q[8];
u3(1.73839078923011,2.08421879787136,-0.766065031379639) q[1];
u3(1.64927209081393,1.19583245423321,-1.51098253376400) q[8];
cx q[8],q[1];
u1(1.55775139199357) q[1];
u3(-3.53381692656702,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.33169769967153,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.625581761291672,-1.09568903779016,3.07593637211392) q[1];
u3(1.42719052541189,-1.08247292138289,0.816091628919721) q[8];
u3(0.626062436326545,2.53044283998215,-0.823432484404814) q[4];
u3(1.22100604655819,2.43723044534345,-1.42998325882132) q[10];
cx q[10],q[4];
u1(0.498518320419737) q[4];
u3(-1.20514777095443,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.28362272111774,0.0,0.0) q[10];
cx q[10],q[4];
u3(0.298162014882396,-1.75121876496544,2.08377774319111) q[4];
u3(2.70631314484272,-2.36473110542064,1.47473825184063) q[10];
u3(2.76565188047498,2.05958055103744,0.860334001840799) q[3];
u3(1.78127697647840,0.369220186697474,-3.71771945386428) q[7];
cx q[7],q[3];
u1(3.02158290354678) q[3];
u3(-2.39563318106851,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.58895542800088,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.03130418872245,-1.41450927419816,0.923435296877505) q[3];
u3(2.29597480745079,0.786590814768197,0.404258137358337) q[7];
u3(0.833832906231075,-1.52744694312990,1.47511571605652) q[13];
u3(0.567702855274126,-3.20567165792469,1.75530523460378) q[12];
cx q[12],q[13];
u1(1.66670814265804) q[13];
u3(-3.53438870812912,0.0,0.0) q[12];
cx q[13],q[12];
u3(2.53582982996833,0.0,0.0) q[12];
cx q[12],q[13];
u3(2.38330238851009,0.234973365858198,-3.72371329460763) q[13];
u3(0.402052113474576,5.22779092782521,0.633650189442163) q[12];
u3(1.49180380231042,-0.429045016695128,0.578834787013147) q[5];
u3(1.23830748517562,-1.88274685597602,-1.15719276821322) q[11];
cx q[11],q[5];
u1(2.27840459244294) q[5];
u3(-1.45106980058808,0.0,0.0) q[11];
cx q[5],q[11];
u3(0.207923765972871,0.0,0.0) q[11];
cx q[11],q[5];
u3(2.77614865148833,2.56446170327375,-2.34340238228423) q[5];
u3(1.72758079695089,-0.0613698746175353,3.14898219460412) q[11];
u3(2.08501962033142,2.44790589403377,-2.33307499538459) q[0];
u3(1.58976573996344,2.49634054350768,-3.68471348465450) q[6];
cx q[6],q[0];
u1(2.17564998310298) q[0];
u3(-1.50339611686670,0.0,0.0) q[6];
cx q[0],q[6];
u3(3.31268720626476,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.906414898363389,-3.95377180179031,2.06379865336902) q[0];
u3(0.262214398692686,-2.14546815753508,-1.36656950309462) q[6];
u3(0.794302576872923,-1.71009247002894,1.17740744411214) q[9];
u3(1.37190673067517,-2.50586668480172,0.318076720376688) q[15];
cx q[15],q[9];
u1(-1.17119306896376) q[9];
u3(0.379488673242836,0.0,0.0) q[15];
cx q[9],q[15];
u3(3.85820158400184,0.0,0.0) q[15];
cx q[15],q[9];
u3(1.16679591180719,1.93974012366287,-2.07493533113635) q[9];
u3(1.77041066711735,-3.11675027372468,-2.71288001239463) q[15];
u3(2.10808021613242,-2.73961303614122,0.810711091460612) q[14];
u3(2.57226920998256,-3.69832359259450,-2.04242876380736) q[2];
cx q[2],q[14];
u1(2.24897090561054) q[14];
u3(-1.92414460779638,0.0,0.0) q[2];
cx q[14],q[2];
u3(3.00029840870321,0.0,0.0) q[2];
cx q[2],q[14];
u3(1.35818273583168,0.796364047993178,-0.381479839434183) q[14];
u3(1.40144756773924,1.35224579048775,-2.89738713847394) q[2];
u3(2.13971799767637,1.63035050535426,-2.48382913017110) q[12];
u3(1.95097587384291,1.65133178239390,-3.68495989355893) q[9];
cx q[9],q[12];
u1(2.16331755274112) q[12];
u3(0.0309849638881194,0.0,0.0) q[9];
cx q[12],q[9];
u3(1.13653335794514,0.0,0.0) q[9];
cx q[9],q[12];
u3(1.55752008780636,0.567331601048634,-3.14035981292097) q[12];
u3(1.68952129130016,0.674964283673161,4.81686496693161) q[9];
u3(0.748086616314155,-1.70812286031402,-0.127129224716372) q[11];
u3(2.06850162673587,-4.41942772495999,0.354996427140807) q[8];
cx q[8],q[11];
u1(-0.137165900638757) q[11];
u3(-1.68228776745367,0.0,0.0) q[8];
cx q[11],q[8];
u3(0.769164741943225,0.0,0.0) q[8];
cx q[8],q[11];
u3(2.36991423642552,4.51567316040498,-1.70182757267811) q[11];
u3(1.97619256983133,-2.84213396323872,0.211510159973106) q[8];
u3(2.13956462198061,0.545839206337099,2.55364090508319) q[10];
u3(1.53894778893154,2.64918502310760,3.51984035110120) q[13];
cx q[13],q[10];
u1(2.24761754539786) q[10];
u3(0.290134851215013,0.0,0.0) q[13];
cx q[10],q[13];
u3(1.47478719777856,0.0,0.0) q[13];
cx q[13],q[10];
u3(1.57992783078724,-2.02678780146588,4.10637948542371) q[10];
u3(1.97450881629652,-1.29324803197146,-3.81778344722651) q[13];
u3(2.33729360396491,-1.94566460416547,3.80134888954607) q[5];
u3(0.598999248787391,2.42601076415977,0.414577898224122) q[15];
cx q[15],q[5];
u1(0.0874852881713943) q[5];
u3(-0.996483680617288,0.0,0.0) q[15];
cx q[5],q[15];
u3(2.66126329348466,0.0,0.0) q[15];
cx q[15],q[5];
u3(1.68434770417066,0.562378342705274,-2.84819090699075) q[5];
u3(1.67983325449446,3.77733644932172,-1.93512720478189) q[15];
u3(1.92501043313610,2.73202271843912,0.194500520768733) q[4];
u3(0.761344091553330,5.07617708771646,1.07985980572915) q[3];
cx q[3],q[4];
u1(0.409740677205363) q[4];
u3(-1.02415400862875,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.67203771861881,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.383219235504540,-4.10665200232635,1.53544625927718) q[4];
u3(2.29961803994650,2.70615035630427,-2.08004478406860) q[3];
u3(0.700206823863976,-2.01414671284137,-0.313319468140658) q[6];
u3(1.05698292517454,-2.51960326395130,-0.753484280041542) q[14];
cx q[14],q[6];
u1(1.34040307477074) q[6];
u3(-0.971408349846032,0.0,0.0) q[14];
cx q[6],q[14];
u3(2.87041264372441,0.0,0.0) q[14];
cx q[14],q[6];
u3(1.16380591564200,1.50737056636048,-2.14368206917063) q[6];
u3(0.711462861718326,-1.15858312283206,-4.71637824724978) q[14];
u3(1.77967700473447,0.150277981630079,-0.870089976372095) q[7];
u3(1.20545816266783,-5.15519159633506,1.12474502524888) q[2];
cx q[2],q[7];
u1(1.44313339813768) q[7];
u3(-3.53108436694792,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.84956153997796,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.32866372014137,-0.979141571235972,0.195134774614199) q[7];
u3(1.18985256178962,0.878846193682110,-2.19325228874999) q[2];
u3(1.12814981883572,0.955014777468565,0.825014563526521) q[1];
u3(1.55312967062430,-0.326598340151911,-3.80290441620696) q[0];
cx q[0],q[1];
u1(1.77723193160747) q[1];
u3(-0.307058614666764,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.68533130985686,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.50189377251006,0.365112618379004,-0.811551364477783) q[1];
u3(1.51289913413568,4.12639312246267,1.52125976164403) q[0];
u3(1.22387349502526,2.61152973229284,-1.48045818181684) q[4];
u3(1.25795982194574,1.11983658619304,-0.115980192549864) q[6];
cx q[6],q[4];
u1(1.39343540545412) q[4];
u3(-1.30229925150629,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.0603015059664045,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.90634389862462,2.48964244494987,-1.84054531470938) q[4];
u3(1.83304987885650,-3.53767545393729,-1.42229392715944) q[6];
u3(1.78844464582238,3.86995112668935,-1.07246811780970) q[7];
u3(1.72305908154730,1.34937452518961,-1.34770591355111) q[1];
cx q[1],q[7];
u1(3.33434382271604) q[7];
u3(-1.12906574951817,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.81021020362061,0.0,0.0) q[1];
cx q[1],q[7];
u3(0.866498537743195,0.626512919040507,-4.03901309313743) q[7];
u3(0.438141367615065,-4.14149858687226,1.50547765416653) q[1];
u3(2.08237603311811,-1.85310591728245,0.315748920142554) q[5];
u3(2.00229987170997,-2.12411540879279,-1.24996701050895) q[11];
cx q[11],q[5];
u1(1.69102281253919) q[5];
u3(0.391650593282901,0.0,0.0) q[11];
cx q[5],q[11];
u3(0.648251844248304,0.0,0.0) q[11];
cx q[11],q[5];
u3(2.05542686557449,-1.69649276520842,0.389682923470016) q[5];
u3(2.16406664069504,-4.37994914677113,0.844909769312690) q[11];
u3(2.26944352003522,2.50870739240217,0.444334155495429) q[0];
u3(2.76912942010507,3.42333302631396,-1.04137625216624) q[10];
cx q[10],q[0];
u1(1.59294575400988) q[0];
u3(0.197805980893319,0.0,0.0) q[10];
cx q[0],q[10];
u3(0.925022001841978,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.56256243798942,-1.97463855786206,-0.478281278269626) q[0];
u3(1.48880704407313,0.346024259051677,-4.21531186294753) q[10];
u3(1.05094322784787,-2.04317583870300,-0.593090148415676) q[8];
u3(2.16596606324660,-4.60883321672515,1.16024526666627) q[12];
cx q[12],q[8];
u1(2.31744000324670) q[8];
u3(-2.06965992516495,0.0,0.0) q[12];
cx q[8],q[12];
u3(1.32928998215642,0.0,0.0) q[12];
cx q[12],q[8];
u3(2.01630107670673,-2.99045831847799,0.481687821099545) q[8];
u3(0.796868665245926,1.13821166188777,-3.31161402859782) q[12];
u3(0.496442237383789,-2.49310980129741,1.45508891769055) q[2];
u3(0.441371793470865,0.757839050968821,-2.44407337434631) q[14];
cx q[14],q[2];
u1(3.03938428788516) q[2];
u3(-1.59624873084217,0.0,0.0) q[14];
cx q[2],q[14];
u3(2.38261472990590,0.0,0.0) q[14];
cx q[14],q[2];
u3(2.48616536643455,-2.88752296345321,-0.539516965086705) q[2];
u3(1.47442403583125,-0.360256536754415,3.64258080423654) q[14];
u3(2.44705705813079,3.28465164540355,-2.99407125485248) q[13];
u3(1.23591018697741,-0.576789612315341,2.23955352300802) q[3];
cx q[3],q[13];
u1(3.37603015814600) q[13];
u3(-1.25140794640972,0.0,0.0) q[3];
cx q[13],q[3];
u3(2.11946081931676,0.0,0.0) q[3];
cx q[3],q[13];
u3(1.16726284053084,-1.04733974003964,0.0199347177564444) q[13];
u3(0.625269670274822,-3.65214024464468,2.35380592321590) q[3];
u3(1.50042938351647,-0.391494484815101,-0.501274920845028) q[9];
u3(2.69927077676806,-4.29281173496553,1.92657461004945) q[15];
cx q[15],q[9];
u1(2.83660221826661) q[9];
u3(-2.23825820422521,0.0,0.0) q[15];
cx q[9],q[15];
u3(1.48330260550867,0.0,0.0) q[15];
cx q[15],q[9];
u3(2.19978004300596,-2.19498437073523,3.15958787826967) q[9];
u3(1.39771362534208,2.61668976518609,-3.04130738386727) q[15];
u3(1.02663255541330,1.06158560670009,0.291872918058226) q[10];
u3(1.79028725420066,-0.0822022081594977,-2.28509729602741) q[9];
cx q[9],q[10];
u1(1.81389266771615) q[10];
u3(-2.18299820496341,0.0,0.0) q[9];
cx q[10],q[9];
u3(3.42349189815110,0.0,0.0) q[9];
cx q[9],q[10];
u3(2.27991505996077,-2.30427162271019,3.06174060302890) q[10];
u3(0.971387645209892,-3.48256450761733,2.11602714649177) q[9];
u3(2.19972812966534,4.02191352462256,-1.15204812382742) q[8];
u3(1.74015892250121,1.77308896578114,-0.346313154983093) q[13];
cx q[13],q[8];
u1(1.14723174008558) q[8];
u3(-0.232316980713978,0.0,0.0) q[13];
cx q[8],q[13];
u3(1.82378329063346,0.0,0.0) q[13];
cx q[13],q[8];
u3(1.86803518646172,-0.458070477373987,1.09330045770091) q[8];
u3(1.30284299894392,2.60361846611137,-0.872171532571136) q[13];
u3(0.538663271909218,0.858654278494636,-0.759297160213282) q[6];
u3(0.750097702197243,-1.28750964519231,-0.116706516713607) q[7];
cx q[7],q[6];
u1(0.0307396054258684) q[6];
u3(-1.38536029873735,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.89369282617742,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.492393754439487,1.21258913497413,0.150540308327566) q[6];
u3(1.16212825052612,-3.74858342584850,-0.603281869641004) q[7];
u3(0.917280849895958,-0.0746057832815540,1.67020158182639) q[15];
u3(1.01012577396607,-0.887150175839510,-2.25345434096491) q[1];
cx q[1],q[15];
u1(2.73444768617791) q[15];
u3(-3.00398889649448,0.0,0.0) q[1];
cx q[15],q[1];
u3(1.17695683370841,0.0,0.0) q[1];
cx q[1],q[15];
u3(2.28269139982147,1.84576716292451,-0.864690334647504) q[15];
u3(2.56233562588256,-3.25661667795084,2.25750516840366) q[1];
u3(1.24237375850159,1.84542029604130,-2.40312864702088) q[12];
u3(2.28442913620508,-2.54595993553893,2.98990670025849) q[5];
cx q[5],q[12];
u1(1.75565878515606) q[12];
u3(-2.52634879079756,0.0,0.0) q[5];
cx q[12],q[5];
u3(0.411233498120760,0.0,0.0) q[5];
cx q[5],q[12];
u3(1.39009574288393,4.48580459999559,-1.48567992554641) q[12];
u3(1.24402028414003,-0.701631924776203,1.79841932231411) q[5];
u3(1.19228665884056,2.25313238819240,-1.79866053174546) q[14];
u3(0.325018604049667,-2.24875099599141,1.45060622111517) q[11];
cx q[11],q[14];
u1(2.42531595376855) q[14];
u3(-0.0265863096457857,0.0,0.0) q[11];
cx q[14],q[11];
u3(1.29116328821290,0.0,0.0) q[11];
cx q[11],q[14];
u3(1.50028205253028,1.02945815678839,-1.48941533599158) q[14];
u3(0.786657674746882,-0.480046699450232,-4.11696065622894) q[11];
u3(1.68066585315952,-0.128741520980581,1.28865912604844) q[0];
u3(1.58627801616212,-1.13422198304584,-0.801198858319466) q[2];
cx q[2],q[0];
u1(1.86309330660884) q[0];
u3(-0.549569495448891,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.00746983311945515,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.92143494365747,1.48918675592730,0.250208285260948) q[0];
u3(0.902397891431093,4.19625224280688,0.208329524605044) q[2];
u3(0.949492003035209,-1.97047397971573,2.11659609781073) q[4];
u3(0.414099400250729,0.114628583044308,-2.60404493301773) q[3];
cx q[3],q[4];
u1(-0.142208819277676) q[4];
u3(-1.13934860845227,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.72699147234519,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.41812401590956,-0.520689839845084,2.65057728479689) q[4];
u3(2.90781417618283,0.700987277891921,-5.30968986627700) q[3];
u3(2.28406246782113,3.86079874376414,-2.27983426297869) q[7];
u3(1.69516198069656,1.25972476434448,-1.54733696925584) q[14];
cx q[14],q[7];
u1(3.22774581621406) q[7];
u3(-0.870665049465359,0.0,0.0) q[14];
cx q[7],q[14];
u3(1.82123398252952,0.0,0.0) q[14];
cx q[14],q[7];
u3(1.61171768456144,-1.51747621087188,-0.527833216309371) q[7];
u3(1.12959162611707,-0.849027585780330,3.26346136307194) q[14];
u3(1.87926901211313,-0.544486418470470,-1.26192106896645) q[15];
u3(2.27714556008508,-3.83386679149799,2.10039976466467) q[6];
cx q[6],q[15];
u1(3.58093747904567) q[15];
u3(-1.26709225127857,0.0,0.0) q[6];
cx q[15],q[6];
u3(2.22139577201287,0.0,0.0) q[6];
cx q[6],q[15];
u3(2.79985752156145,-2.22233263293017,-0.602212898481643) q[15];
u3(2.25441857679330,-3.29044829042465,0.406868719064896) q[6];
u3(1.54425522746111,3.35893889087723,-1.41472586085055) q[12];
u3(2.15836782903034,1.57377631313741,-1.08861382712678) q[1];
cx q[1],q[12];
u1(3.07815591114946) q[12];
u3(-2.14151223347997,0.0,0.0) q[1];
cx q[12],q[1];
u3(0.625022056935673,0.0,0.0) q[1];
cx q[1],q[12];
u3(2.11752437935306,0.785985791260509,2.23662106695709) q[12];
u3(1.35480571691292,1.08742061434724,-0.803564809953186) q[1];
u3(2.24269061071714,0.970905324667833,-0.406625512385341) q[11];
u3(2.04896183335465,-0.273934486022404,-2.28481418066362) q[5];
cx q[5],q[11];
u1(1.50726547693381) q[11];
u3(-0.515034673344192,0.0,0.0) q[5];
cx q[11],q[5];
u3(-0.326136538159307,0.0,0.0) q[5];
cx q[5],q[11];
u3(0.940080097376503,-0.803625427192880,0.0946764033812837) q[11];
u3(2.75131936668304,-0.346983380588801,-5.64232296108733) q[5];
u3(0.785810011687586,0.913961674039532,-2.43631583573791) q[8];
u3(1.82411095093506,2.41864833447455,-3.05405500861747) q[4];
cx q[4],q[8];
u1(0.618913390591118) q[8];
u3(-1.27248577360173,0.0,0.0) q[4];
cx q[8],q[4];
u3(-0.319331564557374,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.760915687397671,-3.54524726678380,1.74165837416921) q[8];
u3(1.87453774700302,-0.244272091579101,0.673286812300235) q[4];
u3(2.43268846075265,0.651226586194257,-1.11605408164993) q[13];
u3(1.77186581749768,0.427067580996566,-4.04009029406685) q[10];
cx q[10],q[13];
u1(2.43852809444114) q[13];
u3(-2.85894007250013,0.0,0.0) q[10];
cx q[13],q[10];
u3(1.82107531190255,0.0,0.0) q[10];
cx q[10],q[13];
u3(1.75805561746517,-3.73636115092280,1.19132039683087) q[13];
u3(1.29526814723394,2.99562840001332,-1.44137858605762) q[10];
u3(1.48777708237592,1.31663342574605,-4.36682987465017) q[9];
u3(0.793030879733779,-1.19811849486384,3.98219396829854) q[2];
cx q[2],q[9];
u1(1.85027746611930) q[9];
u3(-2.67816569699733,0.0,0.0) q[2];
cx q[9],q[2];
u3(0.0685686555626557,0.0,0.0) q[2];
cx q[2],q[9];
u3(2.20114876744258,-2.83469833197927,2.39160225982786) q[9];
u3(2.02219450474351,-5.01827512502584,-0.194052536366020) q[2];
u3(2.20223484998875,-3.56399247914952,2.00366626122947) q[3];
u3(0.604739673523598,4.08767168122468,-1.99135792819916) q[0];
cx q[0],q[3];
u1(0.571212871451766) q[3];
u3(-1.39120401739060,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.16011909871474,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.396943658935716,1.23454863102441,0.271403205491667) q[3];
u3(2.14481714960556,0.0278246410434868,3.33250594430479) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
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
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
