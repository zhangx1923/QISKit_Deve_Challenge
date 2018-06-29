OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(2.49847611776493,0.153143817394887,2.70285541874310) q[5];
u3(2.73463922694357,-2.56349507106501,-0.766498871725293) q[3];
cx q[3],q[5];
u1(1.90780438179713) q[5];
u3(-2.43375833753791,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.0136484159628429,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.57572173352222,-1.73042998881016,3.76069074430918) q[5];
u3(1.21711278363999,-4.35451290019379,-1.71634456946513) q[3];
u3(2.47216636130943,0.987994050347834,-0.741671856236862) q[7];
u3(1.88277277496909,0.799653698727103,-3.93827346257081) q[4];
cx q[4],q[7];
u1(-1.38415627101757) q[7];
u3(0.460969311722131,0.0,0.0) q[4];
cx q[7],q[4];
u3(3.58408368732398,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.19093504413637,1.16030053910804,1.67358037022297) q[7];
u3(1.32266529941499,1.95700517523735,-3.69837786714023) q[4];
u3(2.82632280851189,-1.44510697819362,1.53330023155080) q[8];
u3(2.45817741602415,-2.59028424580720,0.0498987635823058) q[1];
cx q[1],q[8];
u1(1.93085059466099) q[8];
u3(-3.01164393516188,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.42364822508194,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.71974257104997,-0.444143995367604,1.07388536082623) q[8];
u3(2.63764109881972,-1.70029219606708,2.65157356731706) q[1];
u3(0.942320323561449,2.51655119083444,-1.03695545961563) q[9];
u3(1.44304862777575,1.41939751367246,-2.49209148519523) q[6];
cx q[6],q[9];
u1(-0.105861423597838) q[9];
u3(-1.72927208806213,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.793294998152114,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.23895120340393,1.29148965955819,-3.06539378732096) q[9];
u3(1.42086587550357,4.21403651691138,-1.36935753619941) q[6];
u3(2.56667689966398,-1.27249856180436,-0.654059228384316) q[11];
u3(2.12309330163795,-2.57828188612410,1.00579679792588) q[10];
cx q[10],q[11];
u1(3.34219521086850) q[11];
u3(-1.15732960392671,0.0,0.0) q[10];
cx q[11],q[10];
u3(2.46251458836718,0.0,0.0) q[10];
cx q[10],q[11];
u3(1.51443293152716,-4.00771630265977,2.25892114645264) q[11];
u3(1.68289610786351,3.19118586636215,-1.36516615684572) q[10];
u3(1.20458829511715,-0.651226783655249,-1.18240212552330) q[0];
u3(2.00702908866601,-4.53961075340836,1.65006978439393) q[2];
cx q[2],q[0];
u1(2.48924056197005) q[0];
u3(-2.84925182873238,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.52940217952895,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.40655233475216,-2.31166401796390,1.69020142056322) q[0];
u3(0.605230326366357,-0.352832530114156,-2.12739024560057) q[2];
u3(1.39419417173781,1.44560464060098,-3.25781690945357) q[8];
u3(2.34426004716649,2.48227088658518,-2.89608927608110) q[4];
cx q[4],q[8];
u1(0.483547808418282) q[8];
u3(-1.45460326308186,0.0,0.0) q[4];
cx q[8],q[4];
u3(2.89526169638981,0.0,0.0) q[4];
cx q[4],q[8];
u3(0.283980579268267,-0.788410453744285,2.52818428493498) q[8];
u3(1.75774272698712,-4.91579660191582,0.976951009752976) q[4];
u3(2.47377792613478,-0.844404631185937,2.10431991204350) q[1];
u3(2.02964642329059,1.34926417458491,3.95666701325805) q[5];
cx q[5],q[1];
u1(1.10038489173352) q[1];
u3(-0.123124559184006,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.54376473332645,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.07563422775849,1.41776839392588,-2.34149764501248) q[1];
u3(1.06520827443713,1.56424076459174,1.15425604371136) q[5];
u3(2.28709741865532,-1.11273401848740,0.697748558260591) q[0];
u3(1.23341979203284,-2.84902513204737,0.0395249154141011) q[2];
cx q[2],q[0];
u1(0.958076993106707) q[0];
u3(-0.447085567243520,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.95406561337213,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.10443649216299,-0.257107674481726,0.709502907983782) q[0];
u3(0.825690195898236,-2.63992959180242,-1.22760764387613) q[2];
u3(1.28839218011350,1.70975532669438,-3.48616866535178) q[9];
u3(0.464668885427500,-2.33939337037754,3.15867580104004) q[11];
cx q[11],q[9];
u1(0.548495692466066) q[9];
u3(-1.57639488742708,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.34183904997438,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.42180557449915,-2.29054568061344,2.15930479758410) q[9];
u3(2.26529594675305,-1.20431255854396,1.90447791417913) q[11];
u3(1.33682791690128,0.0241654783002048,1.48254813221040) q[10];
u3(1.68963955778668,-0.109709813335406,-2.99534459345490) q[7];
cx q[7],q[10];
u1(1.49275693157913) q[10];
u3(-0.113282761088993,0.0,0.0) q[7];
cx q[10],q[7];
u3(2.43943052736361,0.0,0.0) q[7];
cx q[7],q[10];
u3(1.94715242094481,-0.989058392266937,2.53775784978639) q[10];
u3(0.877732962537359,-3.85457288347019,-2.04557607239444) q[7];
u3(0.744281538522154,0.190166326580046,-1.24339718446586) q[6];
u3(1.55097636834714,-3.52384192464390,1.42864266018274) q[3];
cx q[3],q[6];
u1(1.48695166296980) q[6];
u3(-3.20687717053036,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.76106214811550,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.14557353958513,0.140093717145342,0.470771043793432) q[6];
u3(2.49988208810275,1.87110568309127,0.335229244139345) q[3];
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
