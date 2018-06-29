OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(1.57212498955423,-3.52373451137722,2.50172796919468) q[0];
u3(0.794633550483217,3.32602133627410,-2.29529572576855) q[6];
cx q[6],q[0];
u1(1.66501599127241) q[0];
u3(-3.01120677264237,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.701893803676882,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.51032513047501,0.908713692947406,-3.83533821142074) q[0];
u3(2.29248896865640,1.34611788286001,4.21247565825765) q[6];
u3(1.45089177164815,0.700502556637745,-3.17011809301797) q[5];
u3(1.01256595941733,2.78754650997363,-2.88960173338977) q[1];
cx q[1],q[5];
u1(2.44866705996647) q[5];
u3(-1.50447143953390,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.317859317269974,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.38178277168788,-2.02727742546521,0.786274422134511) q[5];
u3(1.60546297523918,3.65204287204814,0.384631999595777) q[1];
u3(2.67789668321824,1.44913475043423,-4.08227536984997) q[2];
u3(1.21579192947363,3.72623864776838,-2.41075727998040) q[3];
cx q[3],q[2];
u1(1.94261519362385) q[2];
u3(-2.69020340362570,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.367499235425703,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.86039765320217,-0.600284126478637,-0.0493487187394533) q[2];
u3(1.66532025745064,1.06202956991120,2.80558744205776) q[3];
u3(1.32868694629882,1.80460310701282,-1.86732907006509) q[0];
u3(0.398565593317751,-2.71650434810003,3.07605961926265) q[2];
cx q[2],q[0];
u1(3.04113112420002) q[0];
u3(-1.68451089591020,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.517028732665197,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.17856228746570,-2.85525633542024,2.59092799640447) q[0];
u3(1.84590805720115,-0.463818018711844,-4.61168877623505) q[2];
u3(1.55130125779422,-1.55767595237960,4.61714354303092) q[6];
u3(1.86591597957974,1.63319659788288,1.88845847817757) q[4];
cx q[4],q[6];
u1(0.735248378243887) q[6];
u3(-0.233820770964150,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.84185134902987,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.90156599740314,2.09380613141512,-3.92353438595755) q[6];
u3(0.590636710669732,-0.564827810147495,-5.55046662284077) q[4];
u3(2.40099278226878,-1.27097339367151,3.50197792052331) q[5];
u3(1.32988275616746,1.98562334317225,1.43379904149077) q[3];
cx q[3],q[5];
u1(2.27264440364066) q[5];
u3(-2.39965238894937,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.53942195261462,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.47384894702290,-0.785285398946656,-1.04685365699745) q[5];
u3(1.59937336316268,-4.83965546760900,-0.184915561423553) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];