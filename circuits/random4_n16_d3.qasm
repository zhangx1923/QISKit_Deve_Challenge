OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(1.75290393696278,-0.882103238853484,0.345100856015886) q[6];
u3(2.72005922284980,-1.42021003813556,-2.23479402725854) q[4];
cx q[4],q[6];
u1(0.445814136917439) q[6];
u3(-1.13201572638977,0.0,0.0) q[4];
cx q[6],q[4];
u3(3.14942122979873,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.93662400340563,-3.19051361204358,-1.02372910550955) q[6];
u3(1.25346387539497,-1.43062940326031,2.62550323194567) q[4];
u3(2.04931910019714,1.24361221367076,-3.44312538283746) q[3];
u3(1.22379235441626,2.75661907104154,-2.64783633156285) q[7];
cx q[7],q[3];
u1(3.58336560039420) q[3];
u3(-1.59947286581644,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.35068312900120,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.724990313925300,-4.92155226310631,0.510975211851778) q[3];
u3(0.750040266132474,-0.401759294108109,-4.25004151118304) q[7];
u3(1.04903325047326,1.08669031584745,-2.05560150752344) q[11];
u3(0.377246913068231,1.97196620792996,-3.15762233072950) q[2];
cx q[2],q[11];
u1(2.44132621750216) q[11];
u3(-3.14219103899175,0.0,0.0) q[2];
cx q[11],q[2];
u3(1.23569999287537,0.0,0.0) q[2];
cx q[2],q[11];
u3(1.84152839387540,-1.01464732993896,1.19754471877184) q[11];
u3(1.13867220781779,-4.10133405722210,-2.02422654374339) q[2];
u3(2.47008196419501,-3.36785321094974,1.19724253047401) q[8];
u3(2.28983565366751,1.32771498647038,2.96435593194467) q[10];
cx q[10],q[8];
u1(1.41309059255942) q[8];
u3(-0.759766569228332,0.0,0.0) q[10];
cx q[8],q[10];
u3(-0.288090774023804,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.514365498390685,-4.76375900945768,1.15301140118083) q[8];
u3(1.49296439206789,5.30718380786176,0.486120168281301) q[10];
u3(0.893385300818848,-0.839197631274056,2.37401280815535) q[9];
u3(1.12504310994033,-1.96471326545038,-1.82788023528533) q[5];
cx q[5],q[9];
u1(1.56513597919898) q[9];
u3(-0.199270864720872,0.0,0.0) q[5];
cx q[9],q[5];
u3(0.961672517751283,0.0,0.0) q[5];
cx q[5],q[9];
u3(0.224981552378677,0.113406079636705,-0.685767606944758) q[9];
u3(1.41777639692590,3.51371650389987,-2.38956226704410) q[5];
u3(1.04574767436603,1.30487570674090,0.961078555587535) q[15];
u3(0.613733513174028,-0.945468718060252,-2.88131358874621) q[14];
cx q[14],q[15];
u1(0.199472901911557) q[15];
u3(-2.00967348507551,0.0,0.0) q[14];
cx q[15],q[14];
u3(1.58346973283820,0.0,0.0) q[14];
cx q[14],q[15];
u3(0.590634008347315,0.335943556978676,-3.82782883449418) q[15];
u3(0.861288715844957,-5.42315901950925,0.0780351703164723) q[14];
u3(0.832384746243131,3.40394292914613,-1.50066910552123) q[0];
u3(1.88869987047304,1.31640083300865,-1.94904576231659) q[1];
cx q[1],q[0];
u1(1.69932107264334) q[0];
u3(0.529017922240513,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.944239808045231,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.24704145889120,-2.15915064251846,1.17653132485387) q[0];
u3(2.04366021234232,-0.433298744437478,0.163166622061841) q[1];
u3(2.61367725302788,-2.12873004473998,3.70309497029399) q[12];
u3(1.46756502291807,1.32503615776674,1.14768300746963) q[13];
cx q[13],q[12];
u1(1.90175835575971) q[12];
u3(-2.13053765899618,0.0,0.0) q[13];
cx q[12],q[13];
u3(3.54032545172659,0.0,0.0) q[13];
cx q[13],q[12];
u3(0.794025194214149,3.71436010902878,-1.44395868827565) q[12];
u3(0.588346633063472,-2.66101033436059,1.64257580280609) q[13];
u3(2.27904643869529,-3.07228531338720,3.17759915612400) q[3];
u3(1.44997456678144,2.77683741125136,-2.03928864332491) q[0];
cx q[0],q[3];
u1(2.06606940822763) q[3];
u3(-1.61533752970350,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.47028854579969,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.946458956100539,2.10662141989282,-0.861737369125915) q[3];
u3(1.35680729571243,0.0753686563700551,3.50237374941042) q[0];
u3(2.85591116913552,-4.58964027744730,1.53630556733361) q[11];
u3(1.43246963675973,1.61068689579110,-1.17854414406480) q[14];
cx q[14],q[11];
u1(2.37548390232634) q[11];
u3(0.184494209901997,0.0,0.0) q[14];
cx q[11],q[14];
u3(1.13325573826499,0.0,0.0) q[14];
cx q[14],q[11];
u3(1.07415556017185,-0.850502728640974,-1.29799854786211) q[11];
u3(0.851656720189472,4.30446594500752,-1.29241894181351) q[14];
u3(2.24190598534313,-1.55059712881134,0.199434636902056) q[5];
u3(2.09501103572379,-3.05057393386821,-1.10853280157402) q[12];
cx q[12],q[5];
u1(1.84851677976265) q[5];
u3(-3.05525941922821,0.0,0.0) q[12];
cx q[5],q[12];
u3(1.10727502592683,0.0,0.0) q[12];
cx q[12],q[5];
u3(2.17592992024356,0.450888021142008,-1.05530090394629) q[5];
u3(0.350299564977080,-1.08483869606398,3.62831895816853) q[12];
u3(0.491963064717145,-0.509107130750299,-0.871753744874029) q[4];
u3(1.46320906521275,-3.52627995922584,1.97871340097272) q[9];
cx q[9],q[4];
u1(-0.606009209756809) q[4];
u3(0.726561312195010,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.30685131531663,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.988437877848890,0.684918969045673,-2.92277420493315) q[4];
u3(1.69808912864718,2.85054004415876,0.460405722911760) q[9];
u3(1.55647829072978,1.64523946531615,-2.74680633656597) q[2];
u3(1.13147323790085,-2.23551193786285,2.09511288400517) q[6];
cx q[6],q[2];
u1(1.13299747132875) q[2];
u3(-0.246110066443344,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.45235000684001,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.03337130669689,-4.09831117309282,0.401085432972865) q[2];
u3(1.38248112565548,3.88358045481856,-0.150409518568064) q[6];
u3(2.21019159058598,2.56534182579848,0.0830690256767002) q[13];
u3(2.67892216922694,1.00890689515995,-2.92376504860304) q[7];
cx q[7],q[13];
u1(2.04816135425543) q[13];
u3(-2.66110020582620,0.0,0.0) q[7];
cx q[13],q[7];
u3(0.678185119039862,0.0,0.0) q[7];
cx q[7],q[13];
u3(1.76797139313525,-0.403811713332586,-2.16131699015541) q[13];
u3(0.323717238674093,-0.0603358160975820,-2.11191454826304) q[7];
u3(0.570439410314420,0.429488432236811,0.198700398972315) q[15];
u3(1.15070391473530,-0.816385128701868,-2.28288531439445) q[10];
cx q[10],q[15];
u1(-0.311490751194328) q[15];
u3(0.191439940325632,0.0,0.0) q[10];
cx q[15],q[10];
u3(4.09802369699882,0.0,0.0) q[10];
cx q[10],q[15];
u3(1.44103263463834,-1.05315753032540,-0.198246101005255) q[15];
u3(1.76111313473017,-2.18594904947257,1.99507944832075) q[10];
u3(2.73682688665349,2.24804420991158,-3.08391491078372) q[1];
u3(0.887830881757254,-1.81254691770265,3.29427440746490) q[8];
cx q[8],q[1];
u1(1.52012115783303) q[1];
u3(-3.27537241632155,0.0,0.0) q[8];
cx q[1],q[8];
u3(2.29458242164106,0.0,0.0) q[8];
cx q[8],q[1];
u3(0.999257882621242,2.25257196329409,-4.02635005883308) q[1];
u3(1.01061462769489,-2.56621713886750,0.689864030256385) q[8];
u3(1.94678264092503,0.886145541198760,-0.0787076699874628) q[1];
u3(1.98861240903755,-0.177822872039592,-4.38720945075078) q[4];
cx q[4],q[1];
u1(1.78247647636310) q[1];
u3(-1.97038442995713,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.73585690342996,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.15131795248206,1.15833455529680,-3.17409007941070) q[1];
u3(0.847749343874054,-1.38224006730083,-3.15702537659812) q[4];
u3(2.71849453651867,-1.54052372953129,1.69887497761022) q[12];
u3(2.08958331245504,-1.29833272979189,-1.21817700235915) q[2];
cx q[2],q[12];
u1(1.31383696225953) q[12];
u3(-0.750645745739924,0.0,0.0) q[2];
cx q[12],q[2];
u3(0.252927938002013,0.0,0.0) q[2];
cx q[2],q[12];
u3(1.21349569967380,-1.14288752266743,2.44543506889837) q[12];
u3(1.48743624256624,1.44777907022584,3.49346174824725) q[2];
u3(1.76911499058707,0.262570440297503,1.88089646585331) q[0];
u3(2.05447732262992,-1.88012977953699,-0.755739934002973) q[15];
cx q[15],q[0];
u1(2.52254975404165) q[0];
u3(-2.85624935077098,0.0,0.0) q[15];
cx q[0],q[15];
u3(1.77563998170712,0.0,0.0) q[15];
cx q[15],q[0];
u3(1.99018730287273,1.46170302501265,1.19441800518876) q[0];
u3(1.53852839759149,4.32508675410592,0.354431049204237) q[15];
u3(1.88020675506497,0.870158849483725,-1.13463726967203) q[3];
u3(1.27325443364490,-4.63876000586861,1.02052312605143) q[9];
cx q[9],q[3];
u1(2.19848115555306) q[3];
u3(-1.53283082649325,0.0,0.0) q[9];
cx q[3],q[9];
u3(0.710692308578361,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.369999300050418,2.19672778545984,-2.11394563185364) q[3];
u3(1.93868884750526,0.194510070087632,-1.93012983179508) q[9];
u3(1.15439869612523,-0.0683454880949406,1.73962928874840) q[5];
u3(1.66707555269264,-0.355730256166544,-2.27262583609156) q[10];
cx q[10],q[5];
u1(3.88012694554999) q[5];
u3(-4.36346977790042,0.0,0.0) q[10];
cx q[5],q[10];
u3(-0.478380947496329,0.0,0.0) q[10];
cx q[10],q[5];
u3(1.39963870344446,1.61885359284354,-0.563817825850408) q[5];
u3(2.47735801488375,-3.23214306479501,-1.55856130615883) q[10];
u3(1.40502424538627,-0.803804083954519,-0.464215714606326) q[6];
u3(1.34012542194987,-4.09080373614665,0.450233818687346) q[13];
cx q[13],q[6];
u1(0.488368777139828) q[6];
u3(-1.45269787406675,0.0,0.0) q[13];
cx q[6],q[13];
u3(2.56499704003922,0.0,0.0) q[13];
cx q[13],q[6];
u3(0.492030599271383,-2.12609344050777,-0.126818682210858) q[6];
u3(2.12979076347025,3.80465715307336,-1.56111846320645) q[13];
u3(1.32323323959063,3.57006517342066,-2.18481152305951) q[7];
u3(1.27402883080545,1.88579392326387,-1.64030272673345) q[14];
cx q[14],q[7];
u1(2.09090744494298) q[7];
u3(-2.96874208336069,0.0,0.0) q[14];
cx q[7],q[14];
u3(0.659783310358855,0.0,0.0) q[14];
cx q[14],q[7];
u3(1.18535302716164,-0.391955375561861,-0.329639890009145) q[7];
u3(1.46959206569324,0.995615959107251,-0.0959472917624905) q[14];
u3(1.26224872964987,0.245399847550004,2.15085540138963) q[11];
u3(1.40163014406284,-0.623979466982888,-1.57234510051576) q[8];
cx q[8],q[11];
u1(2.34111319219034) q[11];
u3(-1.67878312733804,0.0,0.0) q[8];
cx q[11],q[8];
u3(-0.0424449016348858,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.34365244510707,0.835495508794420,-4.11756256479261) q[11];
u3(0.855447985626650,-1.54642662192877,3.81493639905677) q[8];
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
