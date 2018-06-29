OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[16];
u3(0.963824682238498,0.856863368040812,-2.30558674034923) q[0];
u3(0.972406122451954,-3.71000379159998,2.31154538141994) q[8];
cx q[8],q[0];
u1(2.99410929246618) q[0];
u3(-1.27344067151745,0.0,0.0) q[8];
cx q[0],q[8];
u3(0.768736614791710,0.0,0.0) q[8];
cx q[8],q[0];
u3(2.07139027129377,-0.211245629394546,2.66026391926937) q[0];
u3(1.04732116995869,-0.107820699200345,-2.56973624174633) q[8];
u3(0.940753102832162,1.92401572635175,-0.314478818231972) q[5];
u3(0.696106508331821,1.15025401661680,-2.55054286053996) q[3];
cx q[3],q[5];
u1(1.25359859935498) q[5];
u3(-0.509655661842459,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.90674345149522,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.69478670434888,-2.97607363942903,-0.869627391941533) q[5];
u3(1.70672305807934,-1.75456868788350,-4.17121044391077) q[3];
u3(1.46157023422651,-0.326346072983877,2.91989420383997) q[10];
u3(2.51388027981629,-2.11139323359050,-1.52573766031749) q[11];
cx q[11],q[10];
u1(2.88101890827038) q[10];
u3(-2.02710750073599,0.0,0.0) q[11];
cx q[10],q[11];
u3(0.347449551875239,0.0,0.0) q[11];
cx q[11],q[10];
u3(1.88349372702226,3.68235489558655,-0.384612946724836) q[10];
u3(2.28906349025366,1.88922851706358,3.66505080620004) q[11];
u3(2.47881920906999,-1.38800000789141,-1.45260254455604) q[15];
u3(0.616423937466734,-4.99534960837562,0.159361017676495) q[12];
cx q[12],q[15];
u1(2.61366347010049) q[15];
u3(-2.01080628021268,0.0,0.0) q[12];
cx q[15],q[12];
u3(0.372750681060659,0.0,0.0) q[12];
cx q[12],q[15];
u3(0.849010838629816,3.38298151914991,-1.92836101992048) q[15];
u3(0.823159577338045,-2.22826421740698,-0.316390585534988) q[12];
u3(1.65886402673064,-1.42976536405997,-1.34755198131187) q[13];
u3(0.344147339801438,-4.67452697119664,0.646229341076179) q[14];
cx q[14],q[13];
u1(3.31844551018894) q[13];
u3(-0.802859478142038,0.0,0.0) q[14];
cx q[13],q[14];
u3(1.71260091837769,0.0,0.0) q[14];
cx q[14],q[13];
u3(1.81705413165268,-2.13674366876712,1.24982389074139) q[13];
u3(2.45549070154808,-3.27764314280629,-1.82005315669546) q[14];
u3(2.18625080522203,-0.463513580595115,0.949624222011793) q[7];
u3(2.29796591416419,-1.95407250954862,-2.73030060901836) q[1];
cx q[1],q[7];
u1(0.950777803539236) q[7];
u3(-0.266781926958066,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.50407460349795,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.69390417687676,0.555863058852592,-2.15435288808499) q[7];
u3(0.877766028901340,-2.95738232963919,-1.21544256293726) q[1];
u3(0.855790996773116,2.31299720290798,-1.22609042394805) q[4];
u3(1.94984067977069,0.537049826764216,-3.43062079133869) q[9];
cx q[9],q[4];
u1(1.21796792383548) q[4];
u3(-0.484347859855434,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.04367419398328,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.936651522445520,2.94439845555162,-1.77622519100052) q[4];
u3(0.491518788166037,-3.96690672398206,-1.44858957035837) q[9];
u3(1.27930465572030,3.59316974039462,-0.827601180109635) q[2];
u3(1.37810426772800,2.10763880340055,-1.67493105212897) q[6];
cx q[6],q[2];
u1(3.04615724970780) q[2];
u3(-1.25276368404235,0.0,0.0) q[6];
cx q[2],q[6];
u3(2.30909502132881,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.594371324113961,-1.34646134348084,2.11735465384332) q[2];
u3(0.740464975814940,3.20867484607533,-0.920781366691141) q[6];
u3(0.962033029559693,1.17952247814169,-2.42318432122419) q[15];
u3(1.37842030679246,1.99310412793985,-4.24626213464301) q[0];
cx q[0],q[15];
u1(3.80194061356577) q[15];
u3(-4.36187757386480,0.0,0.0) q[0];
cx q[15],q[0];
u3(-0.476297669131126,0.0,0.0) q[0];
cx q[0],q[15];
u3(0.367876609844550,-1.77316983873631,-2.25377223043772) q[15];
u3(0.597291113247984,1.52969030866314,3.97983959739204) q[0];
u3(1.75588827342469,0.854139133838875,1.13996590361739) q[6];
u3(0.556839675304969,-4.66888389640302,-0.180167748781359) q[8];
cx q[8],q[6];
u1(3.85149192002591) q[6];
u3(-3.27435513471035,0.0,0.0) q[8];
cx q[6],q[8];
u3(-0.936876276161658,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.04237342783016,-0.703112235640355,-0.101313065844188) q[6];
u3(0.917060502215581,2.16898514100421,2.02602754477764) q[8];
u3(2.67939756361530,1.15379633551100,-2.41491711100018) q[5];
u3(2.20613435511620,1.79043172653755,-4.16477776425484) q[2];
cx q[2],q[5];
u1(3.09471927255522) q[5];
u3(-1.20073260674259,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.416855769836385,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.886343623928498,3.56565230990826,-0.734973802625466) q[5];
u3(0.163942133689129,0.684598452719117,2.06661196438302) q[2];
u3(1.26825659172666,0.410034947480346,-1.44199432343048) q[4];
u3(0.991625607960893,1.32251462135075,-3.89999522785676) q[14];
cx q[14],q[4];
u1(0.618316596543033) q[4];
u3(-1.25716223044124,0.0,0.0) q[14];
cx q[4],q[14];
u3(2.69224925442184,0.0,0.0) q[14];
cx q[14],q[4];
u3(0.886253100101938,1.12196885954883,-2.72669659988558) q[4];
u3(2.51591524951345,1.39078464560850,2.95616794940430) q[14];
u3(1.83262526589351,2.14137202410490,-0.701697692321635) q[9];
u3(1.84904835036870,0.505679828383526,-3.83252620796543) q[7];
cx q[7],q[9];
u1(1.01207616769705) q[9];
u3(-1.40836444607502,0.0,0.0) q[7];
cx q[9],q[7];
u3(-0.818121058674058,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.20912237801478,-1.06471057265969,-3.17085214785173) q[9];
u3(1.14872438216826,-0.853413203839563,-2.27187032217709) q[7];
u3(0.811568829000776,1.21809180279932,-0.223585288887958) q[12];
u3(1.40027953161325,-0.0292567509146173,-2.11255699044912) q[10];
cx q[10],q[12];
u1(0.702575487672396) q[12];
u3(-1.36929650846213,0.0,0.0) q[10];
cx q[12],q[10];
u3(3.23410967189436,0.0,0.0) q[10];
cx q[10],q[12];
u3(2.86163978111044,-2.76708101338183,1.77511506117977) q[12];
u3(1.43049838180529,2.12919405077212,-1.43008905586689) q[10];
u3(1.27218156383051,1.08682987155031,0.153104980126849) q[11];
u3(1.61249142401393,-0.223992076281575,-2.91089347790653) q[1];
cx q[1],q[11];
u1(0.900633944569090) q[11];
u3(-0.466394011265926,0.0,0.0) q[1];
cx q[11],q[1];
u3(2.08389837028984,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.903744419927205,3.28024784479122,-2.74143540722235) q[11];
u3(1.36602312394637,-4.28174147101213,-1.76583487746462) q[1];
u3(1.65323844047519,0.864006363964554,-1.98531340464944) q[3];
u3(2.55831711797034,-4.29463027086109,1.40416457407431) q[13];
cx q[13],q[3];
u1(1.87227725141680) q[3];
u3(-2.51268680754823,0.0,0.0) q[13];
cx q[3],q[13];
u3(0.0735201978402247,0.0,0.0) q[13];
cx q[13],q[3];
u3(2.74307307953467,-2.64724280254560,0.322064043840820) q[3];
u3(1.60458860366189,3.08264884493268,-1.60997026130888) q[13];
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
