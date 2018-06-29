OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.29451054936094,0.149121433813559,1.74989217990295) q[4];
u3(1.74149457630082,-1.63799947426859,-0.466608051529630) q[2];
cx q[2],q[4];
u1(1.61718112707660) q[4];
u3(-0.801719455063438,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.80998062904411,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.85323516898946,-1.68527285884300,1.95079373906399) q[4];
u3(0.711841499839246,-2.03587323707239,3.67320522882763) q[2];
u3(2.12586474604011,1.62697686319425,-4.07514460490557) q[3];
u3(1.88568782949484,3.49195685010908,-2.60119130214582) q[0];
cx q[0],q[3];
u1(2.65891502015510) q[3];
u3(-1.56429204290330,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.0977993858942630,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.22766511892414,1.61077954284438,-4.19360416672788) q[3];
u3(0.362470782350635,1.52381913251169,3.81480894393416) q[0];
u3(1.57146858687991,3.82121255701226,-1.37483988594000) q[4];
u3(0.834697486568767,2.44966202594561,-3.07449392875917) q[0];
cx q[0],q[4];
u1(0.675493842289072) q[4];
u3(-1.00141655256652,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.13764453700225,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.70402900805116,-0.541264571857718,0.725152633272918) q[4];
u3(0.541152687909630,-2.64204248308388,-2.50081869671747) q[0];
u3(1.53933413799322,1.23340145933429,-1.02250112191471) q[2];
u3(1.03089095137058,1.55830082789992,-4.64797783649757) q[1];
cx q[1],q[2];
u1(1.98934128534081) q[2];
u3(-2.53725272575449,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.29835445716007,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.21747877601220,0.469968055177782,0.661058363421056) q[2];
u3(2.69156383650306,-1.40372485642207,1.25012680983798) q[1];
u3(0.900149339498227,-1.41713780389704,1.36812540407449) q[4];
u3(0.938723273637491,2.10317174152930,-2.72470497875146) q[2];
cx q[2],q[4];
u1(-0.426785027274204) q[4];
u3(1.16216897376681,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.18540184739921,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.429287790558688,-1.62232452002040,3.63548929225599) q[4];
u3(1.70712008571215,0.831984839494269,-4.34668329062843) q[2];
u3(1.74097272500248,0.517135537414956,0.974226877720797) q[1];
u3(1.19416663334575,-2.92248619684171,-0.592992689213521) q[0];
cx q[0],q[1];
u1(1.34768960925131) q[1];
u3(-3.49776030685546,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.98134656277773,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.497288998107964,-0.270271310074822,-4.35378613131780) q[1];
u3(2.67437519557062,2.85640347338253,0.00124630452395080) q[0];
u3(1.59148879729087,1.50275690335490,-0.823299019540435) q[4];
u3(0.192910969687336,-0.968403118019065,-1.35324103603795) q[3];
cx q[3],q[4];
u1(0.113191586523218) q[4];
u3(-1.25439104486477,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.75259102704351,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.746013661638898,-0.139140260790366,-0.366418456592407) q[4];
u3(1.54719445277819,-2.89096034311430,-1.97413742573241) q[3];
u3(1.06692224501253,-1.29764365702385,-0.121678526987786) q[1];
u3(1.16375250479978,-2.20869590824098,1.09576800013836) q[0];
cx q[0],q[1];
u1(1.00871956774497) q[1];
u3(-1.55673390079254,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.78123001888364,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.18622893249102,-2.32083822812162,0.162304858433294) q[1];
u3(2.28702764840134,0.208421480651697,0.999370021505662) q[0];
u3(0.739071383077811,-2.66969504222056,2.14173264037630) q[1];
u3(0.752266377378686,1.62983734781660,-3.40601370112635) q[2];
cx q[2],q[1];
u1(-1.08275573909987) q[1];
u3(-0.0785167755018497,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.76877080064592,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.291883886799778,-0.413094501422880,0.597745004686059) q[1];
u3(0.965417673484953,0.221586462072476,5.38161795718747) q[2];
u3(0.839104892045351,-0.242468402106445,-2.22613605459747) q[4];
u3(1.79192731878157,1.04767308687193,-4.49946809340453) q[3];
cx q[3],q[4];
u1(1.27354513286273) q[4];
u3(-0.578750567851743,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.69573878092620,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.636543280568809,3.08361131835207,-1.63750444235426) q[4];
u3(1.72694693734902,1.54326041374699,0.0260385846795429) q[3];
u3(0.639698046942883,-1.61715872019237,2.42191582285273) q[0];
u3(0.262666331454216,2.04820882440146,-4.15458834438971) q[1];
cx q[1],q[0];
u1(-0.357290889985071) q[0];
u3(-2.13372393742482,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.17000087185466,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.37836074886618,3.16238038443365,-0.589148199146292) q[0];
u3(2.59034634295261,-3.32868542279368,-1.67918107689387) q[1];
u3(2.36232764814377,-0.432307448733340,2.32422977532780) q[3];
u3(1.94578871130365,-1.09463091538823,-1.53821239198189) q[2];
cx q[2],q[3];
u1(2.10651662539915) q[3];
u3(-0.196882276737161,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.42655872915429,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.43305458332958,3.41553663690288,-0.933634197168086) q[3];
u3(0.770476229209407,2.31058082277469,-2.02903842279172) q[2];
u3(0.379698426385298,-2.52067311789436,3.33594347327327) q[3];
u3(1.15691768768380,-2.68058047750031,1.05309910492834) q[1];
cx q[1],q[3];
u1(1.34679617123601) q[3];
u3(-0.0450374250888970,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.46120343332775,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.71888707571339,-1.83863167483066,4.42026231409553) q[3];
u3(1.60038111871715,-0.633820891319275,-3.19577976684149) q[1];
u3(0.575226134020774,-0.983415794238435,2.04842024814823) q[0];
u3(0.825040181727713,-1.77905395551973,-0.289886071354189) q[4];
cx q[4],q[0];
u1(1.55919490987086) q[0];
u3(-2.34473202935785,0.0,0.0) q[4];
cx q[0],q[4];
u3(3.58438342130924,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.15773746393595,-0.105935007520900,-0.282688980649054) q[0];
u3(1.83037985626236,1.89305886749562,1.94031442969168) q[4];
u3(1.75514084530248,3.60827340184451,-1.66549769829029) q[2];
u3(2.64811156002416,2.29317503708185,-1.17381091077020) q[0];
cx q[0],q[2];
u1(2.29014617075397) q[2];
u3(-2.03081310228069,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.96808755062115,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.21417238957409,1.99271379300401,-2.89601002001917) q[2];
u3(1.78430915934769,-1.09823372335897,-3.35719226138426) q[0];
u3(1.77712819991303,2.58606067564555,-2.70617511197013) q[4];
u3(1.04444156224049,3.18848302274370,-2.14487964849013) q[3];
cx q[3],q[4];
u1(1.53469124414415) q[4];
u3(0.201322989577449,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.06739927122750,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.40002688899956,-0.407858083512663,-1.44577355719135) q[4];
u3(1.46618614698179,-1.77690370039563,3.71902650712115) q[3];
u3(0.407340679902954,-0.405622596493026,0.675808640956160) q[2];
u3(0.776677023768214,-2.23482811282154,0.975409153964081) q[3];
cx q[3],q[2];
u1(2.63301685905091) q[2];
u3(0.159738871515250,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.28958206081938,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.13982774978052,-1.24993515329037,1.30519627533682) q[2];
u3(1.14289814926022,-2.46530003904135,1.41540304163325) q[3];
u3(1.07387146196607,1.21787600861458,1.26481986645564) q[4];
u3(1.14719616453969,-0.220615106692872,-3.16405690370929) q[1];
cx q[1],q[4];
u1(2.49576183303433) q[4];
u3(-1.98692481871290,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.32078324811790,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.81897986107453,2.27916730897373,-3.20453785927638) q[4];
u3(1.64535777963093,-1.67973504635791,2.32200963489404) q[1];
u3(0.683177068471669,0.981858105229771,-3.28787000321943) q[3];
u3(1.78456995581658,3.15133850221442,-3.11242608145244) q[1];
cx q[1],q[3];
u1(1.67651946107788) q[3];
u3(-0.219477168129250,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.08396621688730,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.16472705149422,2.19068300747305,-0.835292707075470) q[3];
u3(2.36169865305020,0.178709677103230,3.93665789593989) q[1];
u3(2.09745864842717,2.79154058345216,-0.983809084828411) q[4];
u3(2.11407979673461,1.79062301973650,-1.00924610095830) q[2];
cx q[2],q[4];
u1(3.72995358396079) q[4];
u3(-3.36185547763190,0.0,0.0) q[2];
cx q[4],q[2];
u3(-1.00670751294543,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.700044854911246,-2.46848793654519,-0.924533539610335) q[4];
u3(1.76576922269803,1.13555047108745,0.367152326665715) q[2];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];