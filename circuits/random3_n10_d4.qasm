OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.08009497728011,1.49666121486241,1.48702407546769) q[4];
u3(1.08739521713887,-4.78801231862447,0.565378412779909) q[0];
cx q[0],q[4];
u1(-1.16826546511239) q[4];
u3(0.161653190117883,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.50498884064100,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.49279557111070,-0.864913112735263,0.179993546659448) q[4];
u3(1.56376185847429,0.388808125726627,1.50297411624880) q[0];
u3(2.72425705354806,0.610040929199716,-1.98047448667061) q[9];
u3(1.98087792436990,-3.82145792177757,2.37510173586826) q[6];
cx q[6],q[9];
u1(2.90570304502647) q[9];
u3(-2.29813690042643,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.63090665463907,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.793436294017087,-0.0624875478564900,0.329753259625587) q[9];
u3(0.257451540774973,0.587016936730157,2.70739581214124) q[6];
u3(1.89458970478938,-1.01640474159067,-1.04143040931001) q[8];
u3(2.23544267385516,1.83037455872654,-3.57555736907078) q[3];
cx q[3],q[8];
u1(1.07687165639777) q[8];
u3(-0.0322964319719383,0.0,0.0) q[3];
cx q[8],q[3];
u3(1.76919579653285,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.20388678608585,-1.48440822426850,1.85281424922786) q[8];
u3(1.02842323352184,-2.20687212177334,-2.27247728589026) q[3];
u3(1.21750947084489,-0.731869916199580,1.75882132911177) q[2];
u3(0.420530046341466,1.89951028665129,-3.37886827805933) q[5];
cx q[5],q[2];
u1(2.99446544236108) q[2];
u3(-2.52264201376982,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.905357813573637,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.687637007285121,-1.15365032267243,1.97669170245990) q[2];
u3(2.04907381963446,-0.635719167811015,-0.440536171979985) q[5];
u3(1.60908441717274,1.66254095390565,-2.86128168414465) q[7];
u3(0.389144694871454,-1.64869412991436,2.81042123788186) q[1];
cx q[1],q[7];
u1(0.801952431455170) q[7];
u3(-0.550375510013940,0.0,0.0) q[1];
cx q[7],q[1];
u3(2.79206857146997,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.64622675583167,0.454903102755165,-3.69261177425560) q[7];
u3(2.72423222525432,0.561022341526089,-0.361967671253775) q[1];
u3(0.741146857025770,-0.512853800134755,0.962155414760125) q[4];
u3(1.04675649514010,-0.955219045961196,-1.74384309786168) q[9];
cx q[9],q[4];
u1(1.74980399814211) q[4];
u3(-3.19112678250781,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.69045286256219,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.72196621180722,2.67938329676582,-0.605751333003792) q[4];
u3(0.513758287300152,3.96800181402471,1.24979630232035) q[9];
u3(1.89982615706112,-0.364230450837551,1.69333679498161) q[8];
u3(1.31601407555829,-1.81483924095062,-2.61642074185086) q[5];
cx q[5],q[8];
u1(1.42837855052112) q[8];
u3(0.486340851648436,0.0,0.0) q[5];
cx q[8],q[5];
u3(0.699956472991753,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.24590475650666,-0.562550741240563,-1.22146986213994) q[8];
u3(1.79521651919045,-1.60065414622264,-4.42614139892881) q[5];
u3(1.41685826611189,-0.812429585289001,3.30032920545571) q[1];
u3(0.958554591462822,-0.236102149863324,-0.957916716963970) q[3];
cx q[3],q[1];
u1(1.48176684444160) q[1];
u3(-3.59884692115250,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.27615593611204,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.43355390515084,0.295071320359230,-0.147776543971273) q[1];
u3(0.959115624214634,1.90133456102675,-1.36145032629821) q[3];
u3(1.73479854036126,0.766019380962806,-3.46263871271014) q[6];
u3(2.45187296087393,3.08182135106107,-2.34293240858391) q[0];
cx q[0],q[6];
u1(1.47060996806472) q[6];
u3(0.127635320682399,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.63365887114876,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.87583783049205,2.51800722297723,-3.40125614592786) q[6];
u3(0.837956892113778,2.44099176784350,0.852378233731981) q[0];
u3(1.96112715678164,0.563924299156485,1.82416530096008) q[2];
u3(2.23859503546948,-1.97147823682843,-1.12907318253839) q[7];
cx q[7],q[2];
u1(0.776934762847493) q[2];
u3(-3.67858507393794,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.74107100520961,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.54469469807892,1.11022662853840,0.210526897349172) q[2];
u3(2.64723554946170,-1.47149715784194,-3.05298553417212) q[7];
u3(1.35198968346772,0.379728509634085,0.629340434444916) q[2];
u3(1.85240958800850,-1.69368310432296,-1.38697971503203) q[0];
cx q[0],q[2];
u1(3.07108151603292) q[2];
u3(-2.39937317004716,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.876574640133879,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.994176477658310,0.674058500378188,-4.35761300157974) q[2];
u3(1.91178131339550,1.40063429108082,-1.03480884608402) q[0];
u3(1.32990855572456,1.45048571826166,-0.0653916479639263) q[1];
u3(1.60234151738763,0.865384100040704,-3.36187947333444) q[6];
cx q[6],q[1];
u1(4.40607379357855) q[1];
u3(-3.85183275013831,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.806582453751413,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.19351341109942,0.144819444916149,2.19056706559890) q[1];
u3(0.777347833633104,2.90370075913158,-3.10788617129254) q[6];
u3(1.36378110092630,0.660594553019660,0.608514477314682) q[3];
u3(1.89807010661244,-0.718937572298717,-1.19055856517889) q[5];
cx q[5],q[3];
u1(0.0450301159130737) q[3];
u3(-0.653634467742292,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.83158571115102,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.06472527119096,-1.87684740774519,-0.370938556761928) q[3];
u3(1.93654003275120,-1.97996427289096,2.69734232893984) q[5];
u3(1.29465242519973,1.21056000572261,1.33425117718119) q[4];
u3(1.35299029480616,-2.03563495377044,-0.647848225389187) q[9];
cx q[9],q[4];
u1(1.54693253437874) q[4];
u3(-0.649224650096764,0.0,0.0) q[9];
cx q[4],q[9];
u3(3.18233906442452,0.0,0.0) q[9];
cx q[9],q[4];
u3(2.08473840744688,-0.298962842909362,0.314776056003518) q[4];
u3(1.73715999486948,-2.77754040302632,-0.489487618734801) q[9];
u3(0.920921783849929,1.00500222057037,2.09493661179047) q[7];
u3(1.90802753582169,-1.58385838102464,-1.29586874441979) q[8];
cx q[8],q[7];
u1(2.18188684581238) q[7];
u3(-1.89473692127000,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.452055630355069,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.95210995894221,4.58700464374115,-1.67248249626361) q[7];
u3(1.67747121930263,4.77511328209076,-1.25838108910660) q[8];
u3(1.94507609116614,-0.562483874744315,1.90525779736841) q[8];
u3(1.97328250117273,-1.15946572706577,-1.38878806197900) q[9];
cx q[9],q[8];
u1(0.878330055864726) q[8];
u3(-1.43506970488011,0.0,0.0) q[9];
cx q[8],q[9];
u3(-0.406629397456221,0.0,0.0) q[9];
cx q[9],q[8];
u3(2.96751623330715,0.0440839239119913,-0.725249810951696) q[8];
u3(2.35439027536620,2.61149763353221,-1.69266357032728) q[9];
u3(0.827308474702376,1.21867346372283,-1.65774643556197) q[7];
u3(0.431455341127046,0.813587485305078,-3.00881124915951) q[1];
cx q[1],q[7];
u1(2.86306732816135) q[7];
u3(-1.83436142669844,0.0,0.0) q[1];
cx q[7],q[1];
u3(3.11852767901764,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.31810885158232,1.07270318538837,1.69692832479507) q[7];
u3(0.471709714194096,1.63949270708437,-2.39485150621938) q[1];
u3(1.35822801267456,1.26146700689508,-3.19405904237905) q[5];
u3(2.09880853708248,2.23887424239139,-3.56492601142051) q[0];
cx q[0],q[5];
u1(2.33701596306061) q[5];
u3(0.126996414019609,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.38556995466913,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.13048871311143,0.616947072843730,-0.713421322359368) q[5];
u3(2.37854827531217,-0.636136805233000,4.65166905113654) q[0];
u3(0.769870611183498,1.42507315792029,-1.57339965866004) q[6];
u3(0.489596276990151,-1.75121253160858,1.05322099857607) q[4];
cx q[4],q[6];
u1(2.37737211137770) q[6];
u3(0.597477102480371,0.0,0.0) q[4];
cx q[6],q[4];
u3(1.68549568571305,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.83343483619794,-0.556927240816567,-2.24799186899056) q[6];
u3(2.45084389565065,-0.669544805091599,-0.900496460098976) q[4];
u3(1.52173143557205,1.07546979437980,1.54014584053237) q[2];
u3(1.28639627395642,-0.534933562705948,-3.02891951622378) q[3];
cx q[3],q[2];
u1(2.82968032058647) q[2];
u3(-1.63124624314790,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.646196536666942,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.28753670159992,-0.309574878407150,-1.62509398347923) q[2];
u3(1.68731315609366,-0.0299364395255603,-0.174481098402179) q[3];
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
