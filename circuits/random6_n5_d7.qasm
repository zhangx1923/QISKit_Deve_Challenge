OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.59006784732561,0.568945649737647,-1.44113670409802) q[1];
u3(2.07510858888924,-4.68317354055339,1.04311213748244) q[2];
cx q[2],q[1];
u1(2.09873680666400) q[1];
u3(0.510623545128596,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.38491054761100,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.35712026625613,-2.59755605144521,0.0641998806088766) q[1];
u3(1.42322235863357,3.86157624796307,-2.11584032021911) q[2];
u3(0.226146119132026,1.12777845728019,-0.347136697737411) q[4];
u3(1.64366109424621,-0.654141496108892,-3.78085126649335) q[3];
cx q[3],q[4];
u1(2.19775478707782) q[4];
u3(-1.77095219858665,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.409917773511713,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.38475307341693,1.28816314776206,-0.852602612873967) q[4];
u3(0.533330138731120,-1.16424486883316,2.94210963588937) q[3];
u3(1.77607521849433,1.86467009569999,0.463775850392467) q[4];
u3(1.89949404364026,0.836801123781062,-2.19924946176312) q[1];
cx q[1],q[4];
u1(0.636501911328095) q[4];
u3(-1.04990834405161,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.73499525041972,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.647194859317384,1.91349808516845,-0.400016221488622) q[4];
u3(1.74164384113932,-2.34680353926640,-0.854205652339543) q[1];
u3(0.609010591156570,0.0423639702725001,1.85063029653005) q[0];
u3(1.32447793260157,-1.95580088367056,-0.599666585595796) q[2];
cx q[2],q[0];
u1(1.67432686033107) q[0];
u3(-2.34844730011966,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.270734035120623,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.854148431845868,-2.25236467704516,3.91111334087305) q[0];
u3(0.763399502266570,-0.462324851992386,-0.961147158005484) q[2];
u3(1.09460014445724,2.00908681685550,-1.03786289912237) q[2];
u3(1.26994102454507,0.406428204845632,-3.58613105399487) q[4];
cx q[4],q[2];
u1(3.89660933002706) q[2];
u3(-3.71030734106590,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.179915269017599,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.24131510443942,-0.999538201017042,0.698386782961993) q[2];
u3(0.859846343700963,4.89617802420227,-0.775433411702991) q[4];
u3(2.03383550383608,2.86409203076041,-1.25127640437870) q[1];
u3(1.80968400611852,1.33971779560665,-1.80407127234713) q[0];
cx q[0],q[1];
u1(1.31297491564300) q[1];
u3(-1.01631857656790,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.81368232560696,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.913258059111180,3.24286110128278,-1.68414014683672) q[1];
u3(0.316246716243504,-0.522354337567636,3.71613126514040) q[0];
u3(1.10729081761762,-3.60842489621864,1.59831570276837) q[1];
u3(1.84999009678935,-0.938550414389456,0.508693869356447) q[0];
cx q[0],q[1];
u1(2.96427701268464) q[1];
u3(-2.19474946936360,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.851938317455266,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.36405700686614,0.860884520012289,-3.22422548401687) q[1];
u3(2.11935887545525,1.21362550059611,-2.87846106559706) q[0];
u3(2.10413373382811,2.94723326080390,-0.884925606215858) q[3];
u3(2.57290699137420,1.11215000091867,-1.95009605512698) q[2];
cx q[2],q[3];
u1(1.10358767935005) q[3];
u3(-3.80475239042519,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.83859200324813,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.87971330315884,2.96875725811143,0.0698223751451217) q[3];
u3(2.95257795692577,-0.503460067982259,-2.61126246268242) q[2];
u3(1.45912662440889,1.00462130812163,0.947182672410151) q[0];
u3(0.374064991037373,-5.10186501810129,0.840542197771586) q[1];
cx q[1],q[0];
u1(1.47577334396006) q[0];
u3(-0.323826583045643,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.95673426815707,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.77840972759383,2.28703081629287,-0.0362413260046659) q[0];
u3(0.548719184977544,-4.30611700542168,1.57315012024909) q[1];
u3(0.742479606453709,0.921146436482376,-0.489612271879042) q[4];
u3(0.628307003574479,-1.13266331592872,-0.903183971303808) q[2];
cx q[2],q[4];
u1(3.26028182545222) q[4];
u3(-1.83505230288857,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.497187467482966,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.72133045728940,2.05214714438944,0.314371203149302) q[4];
u3(0.766730097988269,1.58241105925576,-3.71187640641331) q[2];
u3(1.85732560931865,0.818468528230385,-1.33246853183519) q[0];
u3(1.78316810564317,-4.13915126555581,1.28103864548923) q[4];
cx q[4],q[0];
u1(0.0140350399792126) q[0];
u3(-1.43261615326615,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.07046002102128,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.557978405747574,-0.801732021210172,2.07041400971245) q[0];
u3(0.717812756733027,-1.05917845232009,-2.46627742353340) q[4];
u3(1.10225513962469,-0.345529307300768,-1.94300642433714) q[2];
u3(1.14023372959318,0.562249667747749,-4.41947827766541) q[3];
cx q[3],q[2];
u1(1.05745486793736) q[2];
u3(-1.62601355161666,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.245933194254273,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.973677120642250,1.82413730627442,-1.58153677957239) q[2];
u3(1.42655936966418,5.61327976780154,0.0819295179846882) q[3];
u3(1.89210022049810,0.701352206449400,1.46486214433222) q[2];
u3(2.15355569210116,-1.42350138331261,-0.0191397866623614) q[3];
cx q[3],q[2];
u1(-0.359893056130663) q[2];
u3(-2.03307946648268,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.954083316910426,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.73631303649033,2.88165943079550,-1.16828385126197) q[2];
u3(2.05616347648106,-0.774069955587029,-5.12707023083622) q[3];
u3(1.89728888238779,0.326376465709060,1.63594581738516) q[0];
u3(1.56589768659866,-2.67921690986946,-1.89761854669388) q[1];
cx q[1],q[0];
u1(2.66178374654377) q[0];
u3(-2.98315604608653,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.99399920927062,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.496982908874568,-1.67567405608859,-0.277985453392254) q[0];
u3(1.27180176116520,-2.04021214162638,-3.49413984451879) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
