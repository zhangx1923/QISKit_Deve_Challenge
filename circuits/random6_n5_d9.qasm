OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.761763460858733,-2.05849259790356,2.85913159420657) q[2];
u3(1.29396719239909,0.857518037016975,-1.43129312196075) q[3];
cx q[3],q[2];
u1(0.806801545961210) q[2];
u3(-0.0386831649153263,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.69013439121359,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.41832799382373,0.205604774561027,-2.97776144453283) q[2];
u3(1.74719627354786,1.12585594980419,-3.77961585131973) q[3];
u3(1.71069590972749,-0.883374628513109,1.74391635897824) q[0];
u3(1.67899838128836,-1.57758411292184,-2.32421461436501) q[1];
cx q[1],q[0];
u1(3.01918615423715) q[0];
u3(-1.35634746524056,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.679315498342108,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.17140711377974,-0.876960472366532,1.78223337718704) q[0];
u3(2.39491873161540,-2.49979291129588,-2.19006094020638) q[1];
u3(2.03433519619600,-2.29059469445853,-0.822105212969917) q[1];
u3(1.27359704355594,-4.06696966704323,-0.0410188660214619) q[0];
cx q[0],q[1];
u1(0.601026488280613) q[1];
u3(-1.04481449501812,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.44511715021857,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.66728884336816,-2.54606387610085,2.48811300641647) q[1];
u3(2.19603832009364,2.40495513558591,0.0676124133536906) q[0];
u3(1.10237013201134,2.50841257076275,-2.60979989857009) q[4];
u3(0.843933269022330,2.52690469393836,-2.74634249672784) q[2];
cx q[2],q[4];
u1(3.03740048604235) q[4];
u3(-1.70688959427369,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.18019913838396,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.64385833445940,1.96697762383887,-0.380574028151607) q[4];
u3(1.12587594496776,1.10001463043848,1.08574611697238) q[2];
u3(2.37325254964588,0.127314041422319,-0.205846581032812) q[3];
u3(1.29383290573413,-3.30264972096643,-0.904481715767965) q[0];
cx q[0],q[3];
u1(3.00676486305663) q[3];
u3(-2.14416210554825,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.35657980186202,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.01316799699511,1.06921075708443,1.62756966748012) q[3];
u3(2.23588677813457,3.95532016405362,1.18777853052798) q[0];
u3(1.94128902893236,3.11032990576021,-2.31360181026095) q[4];
u3(1.55688995851006,2.93430511883179,-2.99023285275027) q[2];
cx q[2],q[4];
u1(2.10399899918789) q[4];
u3(0.111778716329925,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.56324343881255,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.93676126400915,1.13566634988447,-3.34508486862280) q[4];
u3(1.15392497998763,0.298983730125210,0.235863379264642) q[2];
u3(1.84377345532256,-1.53607312728953,0.905193709139715) q[3];
u3(1.85178271407598,-3.53497135308144,0.377797741077035) q[0];
cx q[0],q[3];
u1(2.10040354869463) q[3];
u3(0.472331869579843,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.60485199961606,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.279282728803058,0.656756058466307,-4.74262051601067) q[3];
u3(2.13174042672141,4.20893387193331,-0.994804981248454) q[0];
u3(1.42156065980970,1.53598623778779,-2.41518855426113) q[1];
u3(0.782870751697692,-2.59259056813177,2.19787195322509) q[2];
cx q[2],q[1];
u1(2.43450684213333) q[1];
u3(-2.56476241455705,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.34907139756235,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.987876885662435,-2.83826217511274,-0.440225688938818) q[1];
u3(1.64710310058156,0.148660465375910,-3.13278697108775) q[2];
u3(2.47892526677797,2.70983735381367,-3.12913981812506) q[3];
u3(1.38236045717267,-2.66447042431064,3.02141354961269) q[1];
cx q[1],q[3];
u1(4.46619094996336) q[3];
u3(-3.38199139709863,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.314911636818512,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.31161614510080,-0.324761767739347,-0.477852853593172) q[3];
u3(0.892636802202179,1.20211891214879,3.82976169269343) q[1];
u3(1.12047359038124,0.615108545456874,2.37862643925359) q[2];
u3(1.81490672845149,-2.66711120525734,-2.32692926093752) q[4];
cx q[4],q[2];
u1(1.00054681361974) q[2];
u3(-3.19977869265302,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.56272628071706,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.930653944666721,-0.814565453031934,3.56891322334904) q[2];
u3(1.49667169475256,0.359671937142648,4.11391726474512) q[4];
u3(1.81848804399288,0.960553861440444,-2.18212186555519) q[3];
u3(2.89682167536353,-3.54729521860182,2.65059036285711) q[1];
cx q[1],q[3];
u1(1.26725948020740) q[3];
u3(-3.91693821218958,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.95377129248223,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.696326075588727,-1.56239808155616,1.95142332040864) q[3];
u3(1.21754160478076,2.93451432231389,3.18006651555118) q[1];
u3(1.64169906478551,-0.332631413044463,1.51215956668626) q[2];
u3(1.33638979532168,-2.10043953420315,-2.26526989492245) q[0];
cx q[0],q[2];
u1(1.94870752883807) q[2];
u3(-0.319099986466215,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.61421925923657,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.27885523405137,3.77011681990888,-1.84704848827416) q[2];
u3(1.95817518520953,-0.772726366978227,-2.25634316277390) q[0];
u3(1.61120729058229,1.81208797980585,-0.352414998551315) q[2];
u3(2.87852504060219,0.332304881246918,-2.77166074560028) q[1];
cx q[1],q[2];
u1(2.67956934221887) q[2];
u3(-2.08904885410238,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.53746889826097,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.59855756375394,-0.586351184383380,1.61967373357667) q[2];
u3(2.21055663851444,0.580138402168416,1.28570713969427) q[1];
u3(0.856814684698005,-2.23942741433464,1.10077083299867) q[0];
u3(0.503925836641989,-0.867090658302349,-0.976734178159501) q[4];
cx q[4],q[0];
u1(0.580137264808171) q[0];
u3(-3.22985769208872,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.90679625066726,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.59679319085786,-0.532686927690876,1.63542111881061) q[0];
u3(1.75381332124394,-0.125647497554529,2.95758893110530) q[4];
u3(1.39166880353987,-1.10089659426421,1.22701993785887) q[3];
u3(0.891708308052995,-3.07252950198669,-0.287738909103911) q[2];
cx q[2],q[3];
u1(1.38506952568911) q[3];
u3(-0.439996657940859,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.50476559539835,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.48537828103992,3.71677281355000,-1.01907970457579) q[3];
u3(0.720687517799197,-2.98663899764346,-1.90012744894626) q[2];
u3(1.31617755377097,1.46287710830011,-1.00084951174124) q[1];
u3(1.11615123524187,0.464924241554319,-2.92363216422247) q[4];
cx q[4],q[1];
u1(1.51710501944164) q[1];
u3(-0.451529871175583,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.07084549438858,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.18270666798816,-0.707211534332876,3.03294748150528) q[1];
u3(0.273181629384541,-0.0920938019201217,-0.237821040216966) q[4];
u3(1.77087185525120,2.62936853681379,-2.27709108666026) q[3];
u3(1.12253229106885,2.72261547438778,-2.65195375480588) q[1];
cx q[1],q[3];
u1(3.13273746011070) q[3];
u3(-2.31407673626334,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.41169007702404,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.50908215972067,-1.82922028913561,0.640614118350932) q[3];
u3(0.949667193543385,2.33063774417270,1.43910997057583) q[1];
u3(1.92082261250783,3.26670479935155,-0.896944064701568) q[0];
u3(2.06744513039523,0.945653185164139,-1.09224275922096) q[2];
cx q[2],q[0];
u1(3.12940283377257) q[0];
u3(-2.46577579713618,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.869224695632789,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.794714983773319,2.39296136450521,-2.95154572349114) q[0];
u3(0.714701624102886,1.89931729918136,2.68833242063489) q[2];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
