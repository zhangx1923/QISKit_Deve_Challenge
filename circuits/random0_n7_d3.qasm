OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0.797506663057536,3.87391090145146,-1.00971979024786) q[3];
u3(1.95876223553054,2.96475830117634,0.0911514022365247) q[0];
cx q[0],q[3];
u1(-0.0664074001704831) q[3];
u3(-1.97299439160014,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.863183858520349,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.779016785344574,-3.70461955103865,1.57391283748533) q[3];
u3(1.36291623739864,0.750456651783064,-5.17475599339327) q[0];
u3(1.54618110499155,-2.18327714885512,-0.244711891551814) q[6];
u3(2.34383015114563,-3.89047932755388,1.12062620331921) q[1];
cx q[1],q[6];
u1(1.60353908679827) q[6];
u3(0.412854506453627,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.22872295099228,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.69569276917487,-1.83497404115203,2.14007118051428) q[6];
u3(2.47313380426423,3.47016199426379,-0.436564674635593) q[1];
u3(0.125899296368553,-1.61994035006321,0.687759956389356) q[5];
u3(0.467588070691573,-2.96218614622620,1.34759965664882) q[2];
cx q[2],q[5];
u1(1.58485039404676) q[5];
u3(-0.962215558088438,0.0,0.0) q[2];
cx q[5],q[2];
u3(-0.332496773962076,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.78306779234510,-0.369441828533942,-2.21767800195767) q[5];
u3(2.94025873049128,-2.78484440975411,-1.42798803892121) q[2];
u3(1.67814461371242,1.78539441479523,-3.58169094964598) q[0];
u3(2.15556734756410,-3.40664449652161,2.76807335819739) q[3];
cx q[3],q[0];
u1(-0.174060343765198) q[0];
u3(-0.868539361223188,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.02468724076038,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.00065807966436,2.53760957220002,-3.04063820786802) q[0];
u3(2.75968703770699,0.268508506638633,-0.781123566052594) q[3];
u3(0.744039602829611,0.923473668494588,-0.227879231450978) q[1];
u3(1.71795243050694,-0.713686166712820,-2.84884811842510) q[6];
cx q[6],q[1];
u1(2.19175251421039) q[1];
u3(-1.79903450308433,0.0,0.0) q[6];
cx q[1],q[6];
u3(0.308124247456818,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.03489121632478,0.170158074777828,-0.698128959629423) q[1];
u3(1.52316514049999,-1.90201565890318,-0.860987680243487) q[6];
u3(1.43804417409407,1.05964362644170,1.56074704171016) q[4];
u3(0.640610154593015,-0.874353820810072,-2.58148256986951) q[5];
cx q[5],q[4];
u1(3.00372165759535) q[4];
u3(-2.22085550342651,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.593921764903535,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.82336802834623,-3.34848654744045,2.48187374097132) q[4];
u3(2.50816766975235,-2.65417861388033,-2.43534703556202) q[5];
u3(3.02489878870899,2.13657774640977,-0.948985913843609) q[5];
u3(2.67322880945764,5.15354725561641,-0.964818188938703) q[2];
cx q[2],q[5];
u1(0.429517063546964) q[5];
u3(-1.57601327397842,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.25675272466369,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.733531283644854,-2.06044421296783,4.17414044023390) q[5];
u3(0.887847947243759,3.61117219276150,2.35720988266423) q[2];
u3(2.52384135090594,-1.93702797394325,1.62187982958730) q[3];
u3(2.02005455072992,-1.74739665737810,0.335946808299302) q[6];
cx q[6],q[3];
u1(2.77691829951318) q[3];
u3(-1.90466521442159,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.955132509357556,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.589336109205236,-2.71079443708378,2.91116664634229) q[3];
u3(1.49909116506325,-3.95254531762638,0.406886264145606) q[6];
u3(2.17118277607383,-0.466569109129961,-1.41363224100165) q[4];
u3(0.757827559403793,-3.63920714251972,1.33776179293286) q[1];
cx q[1],q[4];
u1(0.578492284369635) q[4];
u3(-1.50898719567107,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.15055436138922,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.39201539106680,3.09990770743493,-2.68239048708415) q[4];
u3(1.41337531391318,1.20087067851707,3.19512534592287) q[1];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];