OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.89260362572727,-1.54026733530807,-0.270600987351745) q[9];
u3(1.70287681816246,-2.11090010387600,1.09022595521357) q[8];
cx q[8],q[9];
u1(1.75506740176182) q[9];
u3(-2.49127507392927,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.984041048809406,0.0,0.0) q[8];
cx q[8],q[9];
u3(0.988263786367734,2.32412449796715,-1.21453012080124) q[9];
u3(1.56178692932729,1.81221794838873,-1.66684401636221) q[8];
u3(2.24056200819068,-0.300596347609755,-0.0464453211333127) q[6];
u3(0.654766475145629,-5.09976903369665,-0.262700239590992) q[0];
cx q[0],q[6];
u1(1.16624270276947) q[6];
u3(-0.577221920875224,0.0,0.0) q[0];
cx q[6],q[0];
u3(3.01398258743634,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.06706931135636,1.04858910240351,0.862010903908139) q[6];
u3(1.72249895147984,-5.56884573974392,0.283316686553074) q[0];
u3(0.467686571868624,-1.28060864464785,1.98203473471169) q[3];
u3(0.749974412595138,-2.20545099523755,-0.349758467222955) q[1];
cx q[1],q[3];
u1(3.06383925243101) q[3];
u3(-0.713824001879964,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.86023835222632,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.68418300528982,-0.450983592635468,-0.464994086885541) q[3];
u3(1.75929381780460,1.49753327462534,3.31734528706278) q[1];
u3(1.37008276074872,0.234815854665842,2.23953477414428) q[7];
u3(0.409525842902665,-0.990495365008665,-1.79428636246938) q[4];
cx q[4],q[7];
u1(0.124049970227741) q[7];
u3(-1.35109464041050,0.0,0.0) q[4];
cx q[7],q[4];
u3(2.94860460139782,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.54200833555217,0.551714432793924,-0.330035479298128) q[7];
u3(1.78913865713421,4.57921448106272,-1.12646400942200) q[4];
u3(1.64267729029271,2.74788324813488,-2.81176832525062) q[11];
u3(1.80556636515240,1.08381004365755,-1.85925746198228) q[5];
cx q[5],q[11];
u1(1.67876388071980) q[11];
u3(-0.416081213200804,0.0,0.0) q[5];
cx q[11],q[5];
u3(0.0922674409298527,0.0,0.0) q[5];
cx q[5],q[11];
u3(1.43010019676851,2.24330778128281,-2.67029967785060) q[11];
u3(1.48438488225359,2.26533820915001,3.07599497359812) q[5];
u3(1.35989329587997,0.370968498329302,0.794362377962598) q[2];
u3(1.83568907908767,-0.569234680106089,-1.07937304284814) q[10];
cx q[10],q[2];
u1(1.63070864032850) q[2];
u3(-2.97125800647002,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.12440380632839,0.0,0.0) q[10];
cx q[10],q[2];
u3(1.72788741386858,0.0971989911220509,1.42354877605600) q[2];
u3(0.977296203296918,0.544889549753579,4.78632468295449) q[10];
u3(1.57023232268698,2.63155711282619,-1.43150552442152) q[0];
u3(0.764478723998986,1.18592164167606,-2.50386767378806) q[7];
cx q[7],q[0];
u1(3.46806128033746) q[0];
u3(-1.57508032177974,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.13697916168856,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.05036655985690,-0.904354863209427,0.379261027809882) q[0];
u3(1.56893564714437,1.57494543914713,-1.61927261486192) q[7];
u3(1.12442794850871,1.03601067302648,-0.444257212395593) q[6];
u3(0.721210672621128,0.0335137464200788,-2.98607262735761) q[10];
cx q[10],q[6];
u1(3.10696862673307) q[6];
u3(-1.35302439227464,0.0,0.0) q[10];
cx q[6],q[10];
u3(2.35106634017084,0.0,0.0) q[10];
cx q[10],q[6];
u3(2.84856009605901,1.78613879054594,0.550004495603802) q[6];
u3(1.18407158306782,1.33473889637310,-3.52242763064325) q[10];
u3(2.27565301131803,0.973896879371275,-3.25669677262039) q[5];
u3(1.49436511910008,2.52994211885095,-2.86995617711849) q[3];
cx q[3],q[5];
u1(3.15551977413141) q[5];
u3(-2.25203182572253,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.56027159641755,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.79673142399758,-0.577116745559039,-1.00082529796015) q[5];
u3(0.428541495221134,-2.76148969867834,-3.47175679103888) q[3];
u3(2.05685477546337,0.665447002532453,1.95769479972439) q[9];
u3(1.22593028717442,-1.05773820631720,-0.890629389627877) q[8];
cx q[8],q[9];
u1(1.50617413055827) q[9];
u3(-0.423993236092457,0.0,0.0) q[8];
cx q[9],q[8];
u3(-0.0664271661070157,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.35373974859897,-1.45788020659266,4.47998341688644) q[9];
u3(0.484478283719962,3.76715732881726,0.0785415121219639) q[8];
u3(1.96980071885821,2.40262123039035,-2.25768487128662) q[1];
u3(0.854939344642602,2.49561540015120,-2.85141503342991) q[4];
cx q[4],q[1];
u1(1.99742026090105) q[1];
u3(-2.48931375618051,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.804601973922449,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.0436668323764499,-0.0441641552808236,-0.564610456803710) q[1];
u3(2.47679298625837,-0.985135295400622,0.715666650654607) q[4];
u3(0.532330244825291,1.60953216739264,-1.11752695960334) q[2];
u3(1.13856024686162,1.22655808539413,-2.03234806084582) q[11];
cx q[11],q[2];
u1(0.0559412188990467) q[2];
u3(-0.577591036887700,0.0,0.0) q[11];
cx q[2],q[11];
u3(2.70001570261176,0.0,0.0) q[11];
cx q[11],q[2];
u3(1.04972234095198,0.457451521845167,-3.21085951669226) q[2];
u3(1.53490675039030,0.810746204876438,5.09919406900598) q[11];
u3(1.61844211701680,-0.836449337643914,-1.67530520028090) q[5];
u3(2.20181787178960,1.11573729353188,-4.62785976449337) q[0];
cx q[0],q[5];
u1(0.130291500629266) q[5];
u3(-1.49229850781410,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.14257327490046,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.78950501009992,-1.55071176956449,0.169006822905341) q[5];
u3(1.91455675770149,2.55871605534888,1.08991926160504) q[0];
u3(1.77660169679678,1.74807081626420,-3.21422753747049) q[8];
u3(2.35633110542199,2.01419745162194,-3.66654551573093) q[4];
cx q[4],q[8];
u1(0.351738922970630) q[8];
u3(-1.52298967987909,0.0,0.0) q[4];
cx q[8],q[4];
u3(1.16474117063301,0.0,0.0) q[4];
cx q[4],q[8];
u3(2.34508535744255,0.765155087856163,-3.67665912575486) q[8];
u3(0.770452987445331,0.173465718344238,2.28991336372084) q[4];
u3(0.635294005070782,-0.256392860465493,0.260352143701177) q[11];
u3(1.33536372405600,-0.331539857242735,-1.06793857362191) q[3];
cx q[3],q[11];
u1(3.03915676554326) q[11];
u3(-0.926092708929972,0.0,0.0) q[3];
cx q[11],q[3];
u3(2.02460677863749,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.14293840705386,1.25093633484829,-2.03193262990726) q[11];
u3(2.24139632837465,-2.05303962499825,2.09246331413979) q[3];
u3(0.421372874358158,0.484895373613514,-0.432923858071198) q[6];
u3(0.840670961237516,-2.24063692489708,1.51750408266221) q[7];
cx q[7],q[6];
u1(1.11274900780886) q[6];
u3(-1.30269254812494,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.636310348628302,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.39766812341992,1.02894728920221,-1.01886236172681) q[6];
u3(2.81329338954763,1.72820349950490,-3.21630396100899) q[7];
u3(2.06263311354043,0.114557459137622,2.21732867361472) q[2];
u3(1.33143110555370,-0.864687252369436,-1.24414340937778) q[1];
cx q[1],q[2];
u1(3.19869215516996) q[2];
u3(-2.06478084212760,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.479242545231447,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.33783823955391,3.30537606449328,-1.89664896628774) q[2];
u3(2.85783766831527,-1.66509585509842,-3.69404457822701) q[1];
u3(2.64096396834640,1.92047374756306,-3.83806181939956) q[10];
u3(0.519788234303231,1.77847980351247,-0.883734790313579) q[9];
cx q[9],q[10];
u1(2.20614463019907) q[10];
u3(0.415772579928851,0.0,0.0) q[9];
cx q[10],q[9];
u3(1.77908244066130,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.10447333210798,-2.33651173525160,3.34638392963305) q[10];
u3(1.63106704410167,4.72987368555152,-0.538881165648867) q[9];
u3(0.870475320196194,-3.44105048110650,2.06144913901549) q[10];
u3(1.52431780099630,-2.72675056750774,3.40385341293030) q[11];
cx q[11],q[10];
u1(-0.0286453452455580) q[10];
u3(-2.16990165551049,0.0,0.0) q[11];
cx q[10],q[11];
u3(1.44081583629766,0.0,0.0) q[11];
cx q[11],q[10];
u3(2.21513956371236,-0.509095213326653,0.883347411749101) q[10];
u3(1.14907907613731,0.591269637205509,2.08075007207357) q[11];
u3(0.911861031952773,1.07252304642842,-3.13155357453624) q[3];
u3(1.66215207178143,-2.86554539454306,3.31263149098859) q[8];
cx q[8],q[3];
u1(-0.178525656007772) q[3];
u3(-2.10262862480460,0.0,0.0) q[8];
cx q[3],q[8];
u3(1.48113732332251,0.0,0.0) q[8];
cx q[8],q[3];
u3(0.835296651866473,-3.10873203434239,0.773975699913042) q[3];
u3(1.99010972839717,4.52845383074992,-1.43571243931884) q[8];
u3(2.11505160806563,1.10005618639251,-3.87953199809555) q[0];
u3(1.22952248673372,2.40630668025410,-2.53572055357398) q[9];
cx q[9],q[0];
u1(3.66777040360316) q[0];
u3(-0.985616585070096,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.58866036124551,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.76348237104399,3.34396963498091,-1.26105121244797) q[0];
u3(1.91315168322060,2.68880650882238,-2.20814302495221) q[9];
u3(0.116665917620173,1.13439201074414,-0.238160709642970) q[7];
u3(0.391228488402337,1.42205993568581,-2.42951682261932) q[2];
cx q[2],q[7];
u1(1.70498619617514) q[7];
u3(-0.0317867077968217,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.670071305133183,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.14961612910982,-1.97688432462012,-1.61267867660407) q[7];
u3(1.75987912872533,-1.10597917356204,4.70754521501187) q[2];
u3(1.00400822857294,-2.17695201030816,4.03337904369764) q[4];
u3(1.90121126254269,1.62535483640215,-2.76579724613272) q[1];
cx q[1],q[4];
u1(-1.06064778928256) q[4];
u3(0.0646728766537485,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.39198386112526,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.61382090456738,2.07380586806531,-1.74325725878295) q[4];
u3(1.21649315101690,-5.13592685295970,-0.0900897973707622) q[1];
u3(2.14329475089127,1.27452497165688,0.560566388619389) q[5];
u3(0.836855188032387,-5.27293553525143,0.285005639434100) q[6];
cx q[6],q[5];
u1(3.36198943745254) q[5];
u3(-0.784242695083587,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.83041548159597,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.06382499390649,0.00350155001447827,1.24457656347141) q[5];
u3(1.00618085304042,-4.44284561811500,1.62591882760227) q[6];
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