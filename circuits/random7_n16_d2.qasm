OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0.351039391912656,2.26676439882378,-1.92031925474749) q[12];
u3(1.01038065885932,0.666884584617056,-2.85125380248817) q[6];
cx q[6],q[12];
u1(2.10369381008341) q[12];
u3(-2.92610939390933,0.0,0.0) q[6];
cx q[12],q[6];
u3(0.479474411373282,0.0,0.0) q[6];
cx q[6],q[12];
u3(1.61832541363586,-1.28244422084144,1.23876653710980) q[12];
u3(1.27529757484345,0.795395569515028,3.31567664534397) q[6];
u3(1.39044834001006,-0.122706241376330,1.78737908960277) q[15];
u3(1.92944203945267,-0.396929137815929,-2.10833859082311) q[1];
cx q[1],q[15];
u1(1.19514562846211) q[15];
u3(-0.208114647974164,0.0,0.0) q[1];
cx q[15],q[1];
u3(2.20712258268298,0.0,0.0) q[1];
cx q[1],q[15];
u3(0.751353709942866,-1.27841947089996,1.97335405199123) q[15];
u3(0.653240420225162,-1.45767783440697,-1.17184464831403) q[1];
u3(2.56888960370733,1.28136291331516,-4.17924807613181) q[13];
u3(1.95640456709557,4.06959025790956,-2.14978196347657) q[9];
cx q[9],q[13];
u1(3.34101039900309) q[13];
u3(-4.26405794805416,0.0,0.0) q[9];
cx q[13],q[9];
u3(-0.481097415970545,0.0,0.0) q[9];
cx q[9],q[13];
u3(1.32705083753431,0.250845015195957,-1.45111319946888) q[13];
u3(1.83752619613279,0.828064759660356,5.19672833375569) q[9];
u3(1.76287784823033,-0.841436919448591,0.871509201453286) q[14];
u3(1.58532938073912,-3.22695366045226,0.0154712830040318) q[11];
cx q[11],q[14];
u1(3.03491197259959) q[14];
u3(-2.03471156872383,0.0,0.0) q[11];
cx q[14],q[11];
u3(0.469324305844152,0.0,0.0) q[11];
cx q[11],q[14];
u3(0.982696133080063,1.19926111780803,0.511648203313007) q[14];
u3(0.669532969970675,-0.876627669276836,-2.87407361617765) q[11];
u3(1.91065390711238,1.84965183632437,-2.38345557693752) q[0];
u3(1.81442382169324,2.07195672180249,-3.73085707154328) q[8];
cx q[8],q[0];
u1(1.52822590992559) q[0];
u3(-0.434878880507366,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.05088930940841,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.211288841276755,-1.71449617474729,-0.382315761337389) q[0];
u3(1.45814081360803,0.322396708318035,-3.23246857191684) q[8];
u3(1.82393352253374,-0.559457755789998,3.38097352429658) q[7];
u3(1.11613618847029,1.58977527116480,1.62196503461204) q[4];
cx q[4],q[7];
u1(3.33148204705604) q[7];
u3(-4.24747455970008,0.0,0.0) q[4];
cx q[7],q[4];
u3(-0.555483578425158,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.74347199351355,-3.54377020123287,2.36323254254031) q[7];
u3(1.74726901160996,3.17595680400544,0.0454970153778615) q[4];
u3(1.84472208578619,0.535306260686987,-3.65468422949793) q[3];
u3(1.44475017516305,2.60256865326289,-2.94003886003373) q[10];
cx q[10],q[3];
u1(3.07264183914365) q[3];
u3(-2.22882986930464,0.0,0.0) q[10];
cx q[3],q[10];
u3(1.59407601214456,0.0,0.0) q[10];
cx q[10],q[3];
u3(1.21265566141633,0.625497519830243,1.41862782622345) q[3];
u3(1.88556052191797,2.30560643002169,1.02620657546369) q[10];
u3(0.813304005257128,1.89887671781585,-0.961118750714803) q[5];
u3(1.61715520491321,0.743361548117448,-3.33045326720325) q[2];
cx q[2],q[5];
u1(0.262700820143583) q[5];
u3(-0.934789554987432,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.55126932777236,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.68391134140129,-2.24216325096444,1.29457156590627) q[5];
u3(0.904845206702575,1.29598569796078,-1.72815154940394) q[2];
u3(1.32112455600922,-0.978512037486353,2.94605604495892) q[10];
u3(1.01461023843714,-2.10656499119396,-1.87513251806917) q[8];
cx q[8],q[10];
u1(3.39923719573225) q[10];
u3(-0.668251576795906,0.0,0.0) q[8];
cx q[10],q[8];
u3(1.21076240856792,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.59018682643950,2.22457158587215,1.66318912968059) q[10];
u3(1.30084110804719,1.39424941471149,4.30250862438694) q[8];
u3(1.32535443323391,0.713362110027540,-3.76797363207498) q[2];
u3(1.01437749585531,2.68046008953319,-2.80873769508047) q[9];
cx q[9],q[2];
u1(2.86600821064790) q[2];
u3(-1.99629777588817,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.60776206087227,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.64792899904047,0.450010576015209,1.88918923873353) q[2];
u3(1.05126534187773,1.72135462266401,1.46742594780592) q[9];
u3(0.658883945476250,0.529323158175406,0.638872482644224) q[13];
u3(1.88660574735583,-0.496693557557643,-2.92375782804567) q[15];
cx q[15],q[13];
u1(1.80911620317847) q[13];
u3(-2.23389229624655,0.0,0.0) q[15];
cx q[13],q[15];
u3(-0.128296263876011,0.0,0.0) q[15];
cx q[15],q[13];
u3(2.67723148460942,0.873268658239144,1.18371251506831) q[13];
u3(1.70321296358993,-3.64226425269977,-1.97037795034535) q[15];
u3(1.95828317722354,-0.0625642154877493,2.79640711083304) q[0];
u3(2.97345195543305,1.02034754923127,2.20163761997211) q[6];
cx q[6],q[0];
u1(1.89855136227543) q[0];
u3(0.165036343187748,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.604594789276928,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.35767957376099,2.98079611981343,-2.68672277218612) q[0];
u3(1.95905719103223,-2.81283850561474,-2.66600609502039) q[6];
u3(0.468602566836074,0.471836729873014,2.01567959399189) q[5];
u3(1.80116511351893,-2.67260163732667,-2.29430502510544) q[7];
cx q[7],q[5];
u1(3.12121906374124) q[5];
u3(-3.62228811021918,0.0,0.0) q[7];
cx q[5],q[7];
u3(-0.910804086281895,0.0,0.0) q[7];
cx q[7],q[5];
u3(2.09058916971855,-2.67248666002027,1.80657801947019) q[5];
u3(2.22491261292139,3.29619179855508,0.219889893212900) q[7];
u3(2.12453419747040,0.666783111412439,-2.96335131254907) q[12];
u3(2.61773486820013,0.943963210425033,-4.67907096224528) q[1];
cx q[1],q[12];
u1(-0.263522774669096) q[12];
u3(-1.53023841260868,0.0,0.0) q[1];
cx q[12],q[1];
u3(0.415159745612872,0.0,0.0) q[1];
cx q[1],q[12];
u3(1.62406180418657,-2.85165859061502,0.942494346501210) q[12];
u3(0.961074273855735,-2.13604094243595,-0.227396873189147) q[1];
u3(0.320387154867092,3.07219350345929,-2.13386575635223) q[3];
u3(1.17895830972916,-3.23738276082273,1.22967476330142) q[4];
cx q[4],q[3];
u1(2.58436181830557) q[3];
u3(-2.04619972944692,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.20311725846888,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.58728626453111,-1.05158408472918,1.48829373116237) q[3];
u3(1.90413427785788,3.16174805066099,2.22551971978358) q[4];
u3(2.01102503671665,4.01500838734468,-1.56716698700584) q[14];
u3(0.688252136795697,2.82066384059168,-0.388383127532177) q[11];
cx q[11],q[14];
u1(-0.0804676788357999) q[14];
u3(-1.88674746809225,0.0,0.0) q[11];
cx q[14],q[11];
u3(0.790389017693334,0.0,0.0) q[11];
cx q[11],q[14];
u3(1.70868812682210,0.120366286734398,-0.578162646930772) q[14];
u3(1.76744956073689,-2.46269447307262,-0.356870219488177) q[11];
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
