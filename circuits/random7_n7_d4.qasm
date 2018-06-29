OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0.673008695126294,-2.73486773529473,1.96445567754722) q[4];
u3(0.797876866366962,-0.806514259213050,-1.75409625643834) q[0];
cx q[0],q[4];
u1(1.74050009313816) q[4];
u3(-3.13519580153297,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.579442390330815,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.893546291360420,0.803114635483585,-1.57959682579639) q[4];
u3(0.918848895193359,4.06963472961261,0.930171339802445) q[0];
u3(2.59501138961404,-2.33843821240744,1.41260634057373) q[5];
u3(2.05052509671163,0.989173037751594,2.02050819628211) q[2];
cx q[2],q[5];
u1(0.302973485743034) q[5];
u3(-1.11747524820484,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.62558114970702,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.22952728702535,-2.15185900760372,2.69206712994547) q[5];
u3(1.34175605318519,1.69759536609737,-1.02057377450212) q[2];
u3(1.55628078665058,0.442550462405822,1.45095584840413) q[3];
u3(1.60002518185158,-1.24220753390986,-2.29579495573951) q[6];
cx q[6],q[3];
u1(2.16015464958052) q[3];
u3(-2.03568528573864,0.0,0.0) q[6];
cx q[3],q[6];
u3(0.157962542267827,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.47427314671931,-0.464976805780068,1.11512064466538) q[3];
u3(1.94499178148747,-3.69417456239190,2.07779384352449) q[6];
u3(0.787163818173458,1.92157421091911,-0.995392033904586) q[0];
u3(1.47627197704626,1.56166599226172,-0.459722898966557) q[3];
cx q[3],q[0];
u1(1.50707461948279) q[0];
u3(-0.233934765535356,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.11104426239624,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.18417628314187,0.351691054422919,0.178885944573828) q[0];
u3(2.25431575998131,-0.0337316392741323,-4.12056726756734) q[3];
u3(2.69869189695615,3.37897741414790,-2.82280871880972) q[4];
u3(1.50722156419797,-0.224524857988623,1.75121294937233) q[5];
cx q[5],q[4];
u1(1.61040518132102) q[4];
u3(-2.85035425331459,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.833688929723428,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.90649165557387,1.37067014111753,-4.80386201066780) q[4];
u3(2.82921323725759,0.332974775022963,1.25292870697616) q[5];
u3(1.90222872263497,0.704007028258391,-1.85759740222292) q[6];
u3(2.26175978142561,0.513341448102134,-3.55816922489418) q[1];
cx q[1],q[6];
u1(2.14818694807905) q[6];
u3(-2.64331298159985,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.18666203800310,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.35992531683193,-4.06323601671175,2.05947101994286) q[6];
u3(1.77713507802197,-3.00979237813982,0.962175014761583) q[1];
u3(1.39933036095087,-1.36019543206216,-0.0101999103987135) q[5];
u3(1.31952061743010,-2.53224151648493,-0.899150115144575) q[3];
cx q[3],q[5];
u1(0.915441964501528) q[5];
u3(-0.280222264595819,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.50680075853115,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.27634494307779,-2.33803503852302,2.06954636526229) q[5];
u3(0.907290033375270,-5.88607490865917,0.205846966924043) q[3];
u3(2.16491315191361,2.51680517918791,0.377101607697925) q[6];
u3(1.18625402179955,-0.253497426914403,-2.92362516941332) q[2];
cx q[2],q[6];
u1(4.05951366979655) q[6];
u3(-3.52929381660144,0.0,0.0) q[2];
cx q[6],q[2];
u3(-0.492523781757628,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.10309461204734,3.27975902370483,-2.37496723942752) q[6];
u3(0.913348531353313,1.88108826263661,3.15420467590325) q[2];
u3(1.13573447216533,-0.213231527689512,-0.940095886940566) q[1];
u3(2.23101812097974,-3.59931618912182,1.92073643580563) q[4];
cx q[4],q[1];
u1(1.31378090642571) q[1];
u3(-3.31656083724856,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.49446514410571,0.0,0.0) q[4];
cx q[4],q[1];
u3(0.705188414140459,4.33901498807381,-1.35608105756328) q[1];
u3(2.03469192748809,-0.270965559889300,0.317857786023980) q[4];
u3(2.10887258318950,1.20132787801052,1.23827220134892) q[3];
u3(1.80243139233623,-1.62353426842929,-1.50069058038248) q[4];
cx q[4],q[3];
u1(-0.0522072173261408) q[3];
u3(-2.15491693770943,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.50426361028070,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.47985656880359,0.945771314118790,-4.47771770208538) q[3];
u3(1.87380325267798,-2.83218025185006,-2.90837542776125) q[4];
u3(0.724913490996391,0.400755451241229,0.138165343013392) q[5];
u3(1.44079353766172,-0.415001613162731,-1.51070145589403) q[0];
cx q[0],q[5];
u1(1.29078953390708) q[5];
u3(-1.42405472263740,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.31464452937685,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.591001841479378,-2.08779041314409,2.57431863740536) q[5];
u3(1.33598954303164,-0.329902504027649,5.76857618943412) q[0];
u3(0.180242300865106,-3.37931168603702,2.12617875274056) q[1];
u3(1.01565741051514,0.0635846473869848,-1.86616529425182) q[6];
cx q[6],q[1];
u1(3.44479731255468) q[1];
u3(-1.18577479181190,0.0,0.0) q[6];
cx q[1],q[6];
u3(2.33112958625306,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.16179249961918,-3.93268401093875,1.64711693511605) q[1];
u3(1.83817655782377,-1.34021814043621,-2.12681635622389) q[6];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
