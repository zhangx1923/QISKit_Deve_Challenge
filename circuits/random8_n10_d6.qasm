OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.20227708546244,0.848684432898546,2.10921374751824) q[7];
u3(1.95146896824673,-1.84704613701596,-2.73272880526251) q[8];
cx q[8],q[7];
u1(1.83029105595184) q[7];
u3(-0.119967109450674,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.68148772102747,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.759613165878186,-4.40978795221924,0.731686574187251) q[7];
u3(2.70219768373102,0.105028579357126,-2.98925420653111) q[8];
u3(0.981967349741184,3.40128112998169,-1.50344517197324) q[5];
u3(1.29556303355298,2.14238356810945,-2.61061335256794) q[0];
cx q[0],q[5];
u1(1.21915044866726) q[5];
u3(-0.203809332811339,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.49422951905726,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.46340861351089,1.96295300206961,-1.40280741318057) q[5];
u3(0.592162098311473,-2.54669917961887,-0.527395754673680) q[0];
u3(0.519250465180159,0.260874124633271,-1.48011451378643) q[1];
u3(1.35845507031466,2.09790531311587,-4.04876582262049) q[6];
cx q[6],q[1];
u1(1.59813749273249) q[1];
u3(-3.42038850807581,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.94555705667443,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.73431875600212,-1.51693152117323,-3.04796409985634) q[1];
u3(1.76176809748893,-0.766050231909369,2.02172204258713) q[6];
u3(2.86009299698711,2.11137371503903,-0.157928778214571) q[4];
u3(1.23665457106594,0.468282659994984,-3.03926176542692) q[2];
cx q[2],q[4];
u1(1.47428557703951) q[4];
u3(0.112164828557910,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.593531163768263,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.45874800986555,-2.00625261302563,-1.76339348292046) q[4];
u3(0.741297557131792,4.99772396214428,0.823472210799225) q[2];
u3(1.00093389400961,-0.00452419442324870,1.53294660592117) q[3];
u3(1.08079522069213,-2.04358871161899,-1.33615973126349) q[9];
cx q[9],q[3];
u1(1.76870619026415) q[3];
u3(-0.534358166690785,0.0,0.0) q[9];
cx q[3],q[9];
u3(2.20660673248942,0.0,0.0) q[9];
cx q[9],q[3];
u3(0.0977907453701838,-2.04178507995158,1.94384045617053) q[3];
u3(0.707285810342425,0.183725513182004,-5.89796679546700) q[9];
u3(0.833799727678861,-2.18901290366446,0.194855626944849) q[2];
u3(1.84537222153234,-3.38698274352228,-0.943215839204615) q[8];
cx q[8],q[2];
u1(4.18277960404343) q[2];
u3(-3.41350007359265,0.0,0.0) q[8];
cx q[2],q[8];
u3(-0.646667512604655,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.04123006895360,-2.25621371457846,0.988065602379643) q[2];
u3(1.35293528499274,4.11841599814517,0.270241808189684) q[8];
u3(1.52983513649209,0.270830982519046,0.799191986159702) q[7];
u3(1.25117546168345,-2.26193903167177,-1.11552535733450) q[1];
cx q[1],q[7];
u1(3.37693624643244) q[7];
u3(-2.34152935768656,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.59443750087859,0.0,0.0) q[1];
cx q[1],q[7];
u3(1.17898585219357,1.51010543729371,-1.34094502878688) q[7];
u3(0.610666206160730,0.716000283166405,2.72222460495724) q[1];
u3(1.69781166128196,2.03859003992600,-2.59257585067721) q[5];
u3(0.736087368040727,-3.09877892717923,2.72493214272575) q[4];
cx q[4],q[5];
u1(3.16874038043994) q[5];
u3(-1.80735093453215,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.873485756900789,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.00339855813527,-1.77507517894726,-0.383446802265388) q[5];
u3(2.34134721343633,5.55362033052992,0.692914475533028) q[4];
u3(2.06096674265815,-1.27357471714592,4.11743207398042) q[0];
u3(0.581798174439975,-1.43747842700191,2.73984257604578) q[3];
cx q[3],q[0];
u1(3.01066181773782) q[0];
u3(-1.74433663767967,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.947108151650918,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.34259863520049,-0.276651446240324,3.03291457674471) q[0];
u3(1.65128997491034,-3.73371520381405,-0.665493595983045) q[3];
u3(0.825015935978611,-3.46733761166880,2.68583769057157) q[9];
u3(1.95701842244128,3.58702413839071,-2.62515519725008) q[6];
cx q[6],q[9];
u1(1.75608108816384) q[9];
u3(-2.80232818669896,0.0,0.0) q[6];
cx q[9],q[6];
u3(0.106509689881486,0.0,0.0) q[6];
cx q[6],q[9];
u3(2.86722110355050,-1.99276726427387,2.81680628343680) q[9];
u3(2.00346835300787,2.20102894142773,-3.37972356414026) q[6];
u3(1.26994495775049,-0.682065370646048,-0.868233899814754) q[5];
u3(2.86852116105783,0.975704161238796,-4.40992884013881) q[0];
cx q[0],q[5];
u1(1.98892866106395) q[5];
u3(-2.99066318864413,0.0,0.0) q[0];
cx q[5],q[0];
u3(0.700250849264339,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.14453137148794,-1.36190969023558,2.81433038060031) q[5];
u3(2.50472140965216,-2.43418535337300,-0.756308521515528) q[0];
u3(1.05004444375304,-2.05489535791045,0.724739262229353) q[2];
u3(0.927616627175563,-2.21129863873661,0.101953726287366) q[6];
cx q[6],q[2];
u1(2.20965472975347) q[2];
u3(0.0969492824808114,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.39533659333122,0.0,0.0) q[6];
cx q[6],q[2];
u3(0.543904592943153,-2.30380849660387,2.27202217055288) q[2];
u3(2.33021100221598,-2.33440538467970,-2.32126448530705) q[6];
u3(1.43323807365744,-1.40557417569710,0.679333728275024) q[4];
u3(0.964045081615222,-1.51819750440276,-0.814343227713760) q[9];
cx q[9],q[4];
u1(2.79211534173590) q[4];
u3(-4.54517186515608,0.0,0.0) q[9];
cx q[4],q[9];
u3(0.475164159683301,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.12170378463074,-0.613450406900012,-2.30866802787090) q[4];
u3(0.711228316958299,1.10516037424473,-3.66833805971234) q[9];
u3(2.52747661932344,1.59884227804861,-3.63391974408914) q[8];
u3(0.670920146241357,-1.89785575708816,3.63071889093608) q[7];
cx q[7],q[8];
u1(1.60937710068943) q[8];
u3(-2.97278490111555,0.0,0.0) q[7];
cx q[8],q[7];
u3(0.811748868512337,0.0,0.0) q[7];
cx q[7],q[8];
u3(1.12870219354695,0.843283621644324,-1.49027114338982) q[8];
u3(2.57423975848246,0.0848110136277600,2.17304960537648) q[7];
u3(0.785960486774683,2.93698839977088,-0.980341253647005) q[3];
u3(1.50980431817653,1.51542002158337,-1.33437808545123) q[1];
cx q[1],q[3];
u1(0.515819794096954) q[3];
u3(-0.307827920788584,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.94932326856729,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.913078919207189,2.68584915174427,-1.80379806672277) q[3];
u3(0.442288887747288,-0.170331801856395,-3.12983027839236) q[1];
u3(2.41279228595037,-0.405132651399480,-1.72040631971322) q[7];
u3(1.35559697614286,-3.53745938515824,1.56273022629964) q[4];
cx q[4],q[7];
u1(2.36088894345672) q[7];
u3(-2.78629432684291,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.24902044984572,0.0,0.0) q[4];
cx q[4],q[7];
u3(0.911935654566517,-1.90865179275756,1.56418250191854) q[7];
u3(0.632386207412158,3.54922063813419,-1.99531670421058) q[4];
u3(2.81324451290897,-1.29039245003271,1.94796956361039) q[3];
u3(2.43296292992874,-2.95963230110764,-0.782130788869575) q[0];
cx q[0],q[3];
u1(3.14050598326223) q[3];
u3(-1.69882314500959,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.496907951337220,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.57533086464320,-4.63061165360268,1.30898733702665) q[3];
u3(1.65408537521056,-4.69509004177712,-0.615587036016185) q[0];
u3(0.704270811332503,-1.24270949457857,-0.0911588026087276) q[2];
u3(0.833856218299919,-2.87609182296833,1.01722407980353) q[1];
cx q[1],q[2];
u1(1.25241672301509) q[2];
u3(-3.20495391312466,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.20013981172169,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.38568424583939,-1.15481703787254,0.745456839566416) q[2];
u3(0.512245776868621,-1.51540030858659,-4.41744955181579) q[1];
u3(2.94865072613086,0.176901619538881,0.581273567967626) q[8];
u3(1.28735374160152,-2.79160044783144,-2.38754528619187) q[6];
cx q[6],q[8];
u1(2.19413690214982) q[8];
u3(-2.89928247378577,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.84884335187084,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.324378292748411,1.16694630511408,-1.05396470710303) q[8];
u3(0.782862429996796,2.84200302621896,-0.618781378849374) q[6];
u3(1.28901226414760,1.63171122993189,-3.44128351889084) q[9];
u3(1.94619324114981,-2.13503244821869,3.83555827496783) q[5];
cx q[5],q[9];
u1(3.18468989868155) q[9];
u3(-0.726788487011316,0.0,0.0) q[5];
cx q[9],q[5];
u3(1.55136299529960,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.53506224008438,1.22569649632241,0.306988431342414) q[9];
u3(1.80481647707328,-2.59487539207925,0.493208100031282) q[5];
u3(1.70304872230210,0.728378115382350,-1.86997169855400) q[1];
u3(0.930729365094234,1.72785350436211,-4.38614122218653) q[3];
cx q[3],q[1];
u1(0.541123500625719) q[1];
u3(-0.194446209545624,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.31371858974569,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.82743334717551,1.72481443316424,-0.851684567549081) q[1];
u3(1.41936603637253,-4.73458903228436,-0.817833082387900) q[3];
u3(2.15996330508020,0.659349895663544,1.27751422500919) q[9];
u3(1.00933224499713,-2.29930253475488,-2.51099532540312) q[4];
cx q[4],q[9];
u1(3.26811445637423) q[9];
u3(-2.47781683172294,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.897873947159663,0.0,0.0) q[4];
cx q[4],q[9];
u3(2.74609777552428,-2.42600921682847,2.29754189797482) q[9];
u3(2.44573729428963,4.55426739219881,-1.22939475421810) q[4];
u3(1.42210045830647,0.412155704626728,-1.99088633371842) q[2];
u3(1.39777317755126,0.548232014224663,-3.86151779020214) q[0];
cx q[0],q[2];
u1(-0.670251079078452) q[2];
u3(1.29959664994003,0.0,0.0) q[0];
cx q[2],q[0];
u3(4.03672744294578,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.868448486199315,-0.660501453900860,-2.34048381500914) q[2];
u3(2.14471889165467,2.73479413584852,-0.149745936898849) q[0];
u3(2.00741095451908,1.61213326688835,-4.61018296103417) q[5];
u3(0.710978963103688,-1.65953881150748,2.99945341333296) q[6];
cx q[6],q[5];
u1(1.78759343969256) q[5];
u3(-2.57820810922542,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.915278317477931,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.12240866954940,3.51389152851776,-2.50594074840230) q[5];
u3(0.0718429163548519,1.56325823608013,0.370676082842487) q[6];
u3(1.56463800047460,1.56834555246671,0.447388994005974) q[7];
u3(2.58176989705650,1.04434137655569,-1.35914721489270) q[8];
cx q[8],q[7];
u1(-0.316672789230823) q[7];
u3(0.611781621549329,0.0,0.0) q[8];
cx q[7],q[8];
u3(4.23536793450574,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.17685857495672,1.87854005981042,1.19097103495816) q[7];
u3(2.28478239939829,-4.67312864292237,-0.883743265947760) q[8];
u3(2.07265911127407,-2.25360493775969,-0.0421400657053472) q[1];
u3(1.48449866525361,-4.73456085284530,-1.51535267443806) q[6];
cx q[6],q[1];
u1(2.68759660111490) q[1];
u3(-2.92334539843380,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.62091729549218,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.81014477797130,-0.766504569512646,0.730739259885415) q[1];
u3(0.264097453610506,0.00598557283734119,-0.597381966262601) q[6];
u3(0.549237206187006,-0.506613905246525,-0.182108356293659) q[4];
u3(1.20976499935700,-2.94696825143923,1.65434383564632) q[5];
cx q[5],q[4];
u1(0.520563022179359) q[4];
u3(-1.06133510515034,0.0,0.0) q[5];
cx q[4],q[5];
u3(3.01770814022375,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.09094278208808,1.19754578287962,1.66846533405873) q[4];
u3(2.80123382503432,0.148657308180559,-4.42400096781015) q[5];
u3(1.58755612523724,3.65774669594617,-0.935757532635159) q[7];
u3(1.56893590170960,2.13652676142151,-0.757384665311290) q[2];
cx q[2],q[7];
u1(0.171640733109125) q[7];
u3(-1.81496072993759,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.41187773778225,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.452512801670834,-1.23600321926292,3.33984180636400) q[7];
u3(1.13583340049044,-5.66119484012482,0.186371000094222) q[2];
u3(2.70266597048720,1.37805130162799,0.875365383050622) q[0];
u3(1.04184429364722,-2.32633646836645,-2.47641220005798) q[8];
cx q[8],q[0];
u1(-0.0344641413662845) q[0];
u3(-2.13762689007170,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.13731919134476,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.24657162618314,-0.719241817406971,2.20945725772263) q[0];
u3(1.78591091523696,0.463332940049244,2.96308874533896) q[8];
u3(2.13043245032699,0.926880752307060,-4.04793054930464) q[9];
u3(1.40776087248327,5.37079800065742,-0.663358854304160) q[3];
cx q[3],q[9];
u1(3.42028999212221) q[9];
u3(-1.10314830264900,0.0,0.0) q[3];
cx q[9],q[3];
u3(2.31968842132133,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.274907550368830,-0.145621033846538,3.03681709429867) q[9];
u3(2.04330002368174,-2.94851466486141,1.79103724531641) q[3];
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
