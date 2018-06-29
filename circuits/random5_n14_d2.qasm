OPENQASM 2.0;
include "qelib1.inc";
qreg q[14];
creg c[14];
u3(1.64021641512571,1.80280045574633,-3.34041833251764) q[13];
u3(1.49153853699531,2.38640801195975,-2.84156246299798) q[2];
cx q[2],q[13];
u1(1.88717050327894) q[13];
u3(-1.46267409619287,0.0,0.0) q[2];
cx q[13],q[2];
u3(-0.185989775987152,0.0,0.0) q[2];
cx q[2],q[13];
u3(1.00963664265413,-3.23942666753999,1.55746388699336) q[13];
u3(1.62550983728322,-2.66740292008894,-2.75339692797210) q[2];
u3(2.37265878718899,1.37713364577817,0.854405138497454) q[1];
u3(2.08676272185697,-0.469018739935894,-3.55450560454621) q[11];
cx q[11],q[1];
u1(2.28108493072838) q[1];
u3(-1.72077197269294,0.0,0.0) q[11];
cx q[1],q[11];
u3(0.0839049167176675,0.0,0.0) q[11];
cx q[11],q[1];
u3(0.885762915555180,2.58597885031886,-1.99689787171078) q[1];
u3(1.07172764103709,-1.51994470523723,4.46833583977268) q[11];
u3(2.03197253587053,2.05461633428892,-4.00223506842209) q[4];
u3(2.52151068143997,-3.96340726924776,2.27277869265145) q[8];
cx q[8],q[4];
u1(1.61400279881757) q[4];
u3(0.0677681880079883,0.0,0.0) q[8];
cx q[4],q[8];
u3(0.712731537978718,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.11313092641469,1.37351193120464,-1.47903207853068) q[4];
u3(1.49906567463140,5.65156051022208,0.631366769946338) q[8];
u3(2.12734004590499,-1.19008269827238,0.419974275203232) q[7];
u3(1.68584069930741,-2.10528546462641,0.597328230528632) q[6];
cx q[6],q[7];
u1(2.00920778042450) q[7];
u3(-2.32363753232672,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.556553161272693,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.34847283516073,-1.27078733829643,4.00188586953003) q[7];
u3(0.821225882741847,3.29584937286599,0.647130956531642) q[6];
u3(2.23662990032779,-0.277619931979683,0.0657260135141046) q[10];
u3(1.67564590802791,-2.58256142066718,-1.61406834653738) q[3];
cx q[3],q[10];
u1(3.18353816071285) q[10];
u3(-2.49763019526238,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.45216961541938,0.0,0.0) q[3];
cx q[3],q[10];
u3(1.31426473204222,-0.0763928993068993,0.195977971368224) q[10];
u3(1.88273390302730,1.19320131218380,1.10928100227302) q[3];
u3(2.22120075224493,3.64084898212045,-1.35262528894216) q[12];
u3(1.33161049093453,1.25804917888973,-0.853417326975452) q[5];
cx q[5],q[12];
u1(1.63415442622085) q[12];
u3(-0.752266729051510,0.0,0.0) q[5];
cx q[12],q[5];
u3(-0.288787064615293,0.0,0.0) q[5];
cx q[5],q[12];
u3(2.84562375847645,0.196135135021661,-1.03238022010050) q[12];
u3(1.67447664902670,2.87682212504383,2.27284880011290) q[5];
u3(0.485044442261362,0.995588979776579,-2.31596592952766) q[0];
u3(1.21941591676152,2.52497168823467,-3.56249363599742) q[9];
cx q[9],q[0];
u1(3.29850155203826) q[0];
u3(-1.76491450054518,0.0,0.0) q[9];
cx q[0],q[9];
u3(2.64627661369193,0.0,0.0) q[9];
cx q[9],q[0];
u3(1.10874389364540,-3.86645684224067,1.77757403343206) q[0];
u3(1.90485427628643,-3.53376398688229,-2.63485816913236) q[9];
u3(0.791778196119410,1.83566548023456,-1.30695637950441) q[10];
u3(0.329107679331019,-2.71194405848120,1.25781628525934) q[4];
cx q[4],q[10];
u1(3.12914398985929) q[10];
u3(-1.19568635862783,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.69122528377599,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.79271892982822,0.0199703301792105,1.32582749586234) q[10];
u3(0.968844453667510,-1.94048546519038,1.49456273356191) q[4];
u3(2.91367486495641,1.47693253523950,1.09983500365879) q[8];
u3(0.949763476854933,-2.13649012945163,-3.01814236555249) q[11];
cx q[11],q[8];
u1(1.83071888842839) q[8];
u3(0.473741273732990,0.0,0.0) q[11];
cx q[8],q[11];
u3(0.928946668209605,0.0,0.0) q[11];
cx q[11],q[8];
u3(2.27562644842489,-3.05495204764419,2.38419478792586) q[8];
u3(1.36432601626010,0.491682996556899,-3.07827571478495) q[11];
u3(1.37031242748709,-0.205962671942480,-1.04511067113669) q[3];
u3(2.04314597641451,-5.37457528549042,0.821537619411036) q[1];
cx q[1],q[3];
u1(1.14719607720986) q[3];
u3(-0.0304653228526979,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.83741729960195,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.64443072437934,-1.76577601008076,4.48399701590235) q[3];
u3(2.82220541825046,3.19585217630087,-0.143046742979649) q[1];
u3(1.00123588786338,3.26939916813790,-2.16676572144489) q[0];
u3(1.56794466587051,1.48347607547936,-1.48126517838148) q[6];
cx q[6],q[0];
u1(2.99318963229074) q[0];
u3(-1.82500701827065,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.874957882731496,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.33990512694641,0.910327476026659,-1.70881438238837) q[0];
u3(2.20453166621589,-0.646971898900399,-3.41138131377524) q[6];
u3(1.72923180365883,1.91660039126506,-0.134328504392897) q[7];
u3(0.760269642719869,0.0469366205435813,-2.56247407098000) q[12];
cx q[12],q[7];
u1(2.37329489340709) q[7];
u3(-2.61760056086541,0.0,0.0) q[12];
cx q[7],q[12];
u3(1.20643252910064,0.0,0.0) q[12];
cx q[12],q[7];
u3(1.72811667359400,0.629260027109457,-3.57858687055466) q[7];
u3(1.34034602882997,-2.42069587317928,2.63850043433118) q[12];
u3(2.13970236572954,3.69850948113363,-2.23417792917179) q[9];
u3(2.10028146126726,1.78147086665056,-2.15662265861958) q[2];
cx q[2],q[9];
u1(-1.11868472091845) q[9];
u3(0.537444714964030,0.0,0.0) q[2];
cx q[9],q[2];
u3(3.36514060586372,0.0,0.0) q[2];
cx q[2],q[9];
u3(1.31874982481545,2.51499244482553,0.168522284601829) q[9];
u3(2.45892361311270,-3.85565611190282,-0.736072625032984) q[2];
u3(0.978255505602103,-0.867863681067902,1.93986120238068) q[5];
u3(0.921458227784909,-1.78783806791191,-1.18210134419013) q[13];
cx q[13],q[5];
u1(2.69116572917572) q[5];
u3(0.0959666429584889,0.0,0.0) q[13];
cx q[5],q[13];
u3(1.55084989946446,0.0,0.0) q[13];
cx q[13],q[5];
u3(1.75779802376127,-2.50799931945307,-0.414721154690767) q[5];
u3(2.25025018421353,3.44317243343990,-2.34063577535930) q[13];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13];
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
