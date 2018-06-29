OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(0.986870582389085,-3.13056281315472,2.46873461912606) q[4];
u3(1.23549630626223,0.194230367421192,-1.77848800866712) q[3];
cx q[3],q[4];
u1(-1.17392434373175) q[4];
u3(0.703357475649660,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.54183624664937,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.955617236839955,0.813904783511160,2.39011658662463) q[4];
u3(1.45967490237465,-1.04561136323003,-2.34892449448483) q[3];
u3(1.59588610589937,0.541230133552970,-1.53850390266073) q[1];
u3(1.01146544519561,-3.86560205957195,1.64145049326524) q[9];
cx q[9],q[1];
u1(2.58734053350449) q[1];
u3(-3.05387095331197,0.0,0.0) q[9];
cx q[1],q[9];
u3(0.771830732966281,0.0,0.0) q[9];
cx q[9],q[1];
u3(0.927593527083934,3.50751665696499,-2.27136678697226) q[1];
u3(2.62628376014544,0.559422333298036,5.09323866865126) q[9];
u3(1.13307497746973,-1.87863713843051,0.617192398771257) q[8];
u3(1.05442953764237,-1.99464281866013,0.927746395215256) q[5];
cx q[5],q[8];
u1(4.32697467235991) q[8];
u3(-3.41838909379729,0.0,0.0) q[5];
cx q[8],q[5];
u3(-0.540973163961264,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.40574248173523,0.369532037904038,-2.29922741345739) q[8];
u3(0.918426369825054,-1.91364519623350,2.97241961909471) q[5];
u3(2.08789419661407,1.28110923323965,-2.74846001196795) q[2];
u3(1.45745344465653,2.43118536241223,-3.81905039664798) q[7];
cx q[7],q[2];
u1(0.708150544849525) q[2];
u3(-1.48782231080532,0.0,0.0) q[7];
cx q[2],q[7];
u3(3.01170703667033,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.20381447561444,-1.60163166184752,4.10760348475022) q[2];
u3(1.73663536515033,0.287368465404943,-0.157097121570135) q[7];
u3(0.682801917333681,0.907436747305625,-2.21968245838004) q[6];
u3(1.53038196980255,-3.59955684465784,2.39150821831760) q[0];
cx q[0],q[6];
u1(1.47640159969431) q[6];
u3(-2.95121840176075,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.09214683557619,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.34086663763039,0.584502796977862,1.57353174318562) q[6];
u3(1.14643399806924,1.16339351814030,-0.0806976033998226) q[0];
u3(0.415348797671132,-0.570628549453605,1.12480057966487) q[0];
u3(0.903964863072561,-2.31695312980275,0.914550356317632) q[9];
cx q[9],q[0];
u1(3.41698938463231) q[0];
u3(-1.78498527706336,0.0,0.0) q[9];
cx q[0],q[9];
u3(1.47031494421108,0.0,0.0) q[9];
cx q[9],q[0];
u3(2.43348769092641,1.01514095478208,-0.548362838962488) q[0];
u3(1.34155247456269,1.63966677265860,-2.93308280445079) q[9];
u3(1.53067145756282,1.37196683533765,-1.67770661303383) q[4];
u3(0.536076480682373,1.68350995977675,-3.26537108153288) q[1];
cx q[1],q[4];
u1(0.146484015794532) q[4];
u3(-1.72123906570491,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.584501402442403,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.98527712262611,-0.887865266425629,2.89659201967627) q[4];
u3(0.825207262539613,3.42422958563568,1.27058517246970) q[1];
u3(2.25372979869574,0.808704061676671,-2.74467175555229) q[3];
u3(2.12809894048226,0.522620142313610,-4.93873409677381) q[5];
cx q[5],q[3];
u1(3.72814843908576) q[3];
u3(-3.27648642871970,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.852589907877172,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.47155944633564,2.70015981714500,-3.02480831523487) q[3];
u3(2.27285199724592,-3.57272365285295,1.20876033266584) q[5];
u3(2.19671186652966,-2.56285484236101,-0.163427641638672) q[6];
u3(2.66292588160943,-3.17430819086422,-2.09230907737692) q[7];
cx q[7],q[6];
u1(1.38808248884889) q[6];
u3(-0.173614598424350,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.29654568241618,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.533563202722146,2.95508549048618,-0.417976806664474) q[6];
u3(2.71703918434851,0.437073418635064,3.25687907502474) q[7];
u3(1.97790726755949,3.11735386662184,-2.39377712068662) q[2];
u3(2.13486307327414,1.25658354148106,-2.31708927887535) q[8];
cx q[8],q[2];
u1(2.99890603152210) q[2];
u3(-2.39668018847304,0.0,0.0) q[8];
cx q[2],q[8];
u3(1.15078233216761,0.0,0.0) q[8];
cx q[8],q[2];
u3(1.34926442899648,0.815960947751437,-2.41972055423258) q[2];
u3(1.68093478472462,-5.27306886253920,-0.421954525174294) q[8];
u3(2.70902203937590,2.37153602787728,-0.331965622082693) q[4];
u3(2.23918023331455,4.14021446093689,-0.104473834158621) q[9];
cx q[9],q[4];
u1(2.86248660225002) q[4];
u3(-1.87444944786969,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.652915592517086,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.97807995437382,2.97128736072313,-0.325133970114072) q[4];
u3(1.97194499043390,-3.37777416503885,-0.0643604996241864) q[9];
u3(1.86021111147983,-0.909013961770035,-0.0379007804199462) q[3];
u3(2.19267514379855,-2.62462515936796,1.18910140815922) q[7];
cx q[7],q[3];
u1(0.551838979157004) q[3];
u3(-1.29186765677859,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.0341625100525125,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.77541873425454,-1.03298347406208,-0.388884167985866) q[3];
u3(2.21938231614817,0.112730367991538,5.49585329376067) q[7];
u3(1.57677158253428,-1.72219538867004,-0.213083795568256) q[0];
u3(2.95163479503216,-2.30723747962993,-0.708563010369729) q[8];
cx q[8],q[0];
u1(1.48880182587682) q[0];
u3(0.0821951878282143,0.0,0.0) q[8];
cx q[0],q[8];
u3(2.53517947938901,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.929556361856097,1.67057835879409,-2.40756090748441) q[0];
u3(2.44251144133723,0.335154784346247,-2.94581674816987) q[8];
u3(1.53263337656542,1.98077136632601,-3.53363537790184) q[1];
u3(0.699167404394747,-2.50300456395874,3.16383130365965) q[5];
cx q[5],q[1];
u1(1.57379924211545) q[1];
u3(-0.745820174378439,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.141649964039342,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.85511543677705,0.601808967441297,0.658410776227783) q[1];
u3(2.99059100656209,0.974389291429859,-0.427978898420467) q[5];
u3(1.86198560791320,-0.0838972413541477,2.60896794651395) q[2];
u3(2.22712848492806,-1.85374736947292,-2.16003051239200) q[6];
cx q[6],q[2];
u1(1.47911999189807) q[2];
u3(-3.41561654492896,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.10654849794870,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.57328185254874,-0.372104151543037,0.527944640311510) q[2];
u3(0.721739223171978,-4.95962878883868,0.164552393682218) q[6];
u3(1.28534075373562,2.44092382580338,-0.737476976838738) q[3];
u3(0.886803760787635,-0.0592925953501373,-3.49499875245725) q[9];
cx q[9],q[3];
u1(1.67353528692488) q[3];
u3(-0.247510676874550,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.79211116879643,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.741742580788716,-2.85902721466427,1.03655831153399) q[3];
u3(1.44291733109390,-1.02079825992799,-1.68829174906189) q[9];
u3(1.85586174063586,0.936528168641722,-2.11588323986382) q[6];
u3(1.96532422135742,-2.76033727968947,2.94866776364462) q[8];
cx q[8],q[6];
u1(0.735461093346343) q[6];
u3(-0.00678867243378489,0.0,0.0) q[8];
cx q[6],q[8];
u3(1.79771970777683,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.53735448744595,-3.36480749871244,2.10599587737958) q[6];
u3(1.41757934924277,-0.106564234608708,-2.40517926451353) q[8];
u3(2.17217027793090,1.96172077729210,-1.27723644812408) q[1];
u3(2.06279320643045,4.74316103319927,0.813196209653869) q[2];
cx q[2],q[1];
u1(2.72368471176762) q[1];
u3(-1.82477094253073,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.27770128106321,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.74871393157704,3.41487214695419,-2.36977545306708) q[1];
u3(1.64645605000912,-4.99973906154139,1.14090655578308) q[2];
u3(0.539815468245551,-0.960343321159215,0.819845832170245) q[7];
u3(0.575694679960301,-0.896252437435781,-0.802097969367768) q[4];
cx q[4],q[7];
u1(2.99934201398606) q[7];
u3(-4.25553500811773,0.0,0.0) q[4];
cx q[7],q[4];
u3(0.581920836778655,0.0,0.0) q[4];
cx q[4],q[7];
u3(3.02721129813261,-0.194039141042377,0.138492900437481) q[7];
u3(2.90081575951223,1.40963361698557,2.26685333625594) q[4];
u3(2.09959147018107,3.40108878982116,-2.28231164310535) q[5];
u3(1.69655126692273,1.40533087700863,-1.82259809971156) q[0];
cx q[0],q[5];
u1(-0.0218408751386678) q[5];
u3(-1.65465351734923,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.903828248401429,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.91535281923240,2.55204197784020,-2.36747736577485) q[5];
u3(1.09177797666607,-0.433220758940933,-5.68128297823335) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
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
