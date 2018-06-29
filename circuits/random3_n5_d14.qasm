OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(2.07358060594756,2.76150449372882,-2.83201140834807) q[2];
u3(1.75529094445871,-3.46818760989254,1.86025014909393) q[0];
cx q[0],q[2];
u1(1.70115989738535) q[2];
u3(-3.51792469057610,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.31087974881329,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.95776048569977,3.32192194550508,-2.89709576122981) q[2];
u3(0.697342485644378,-2.92968589524454,3.09437256964760) q[0];
u3(1.20344642601794,-0.116409118217658,1.45293523217841) q[1];
u3(1.64248865176784,-1.05377738457278,-2.80323478106213) q[3];
cx q[3],q[1];
u1(1.81304390224753) q[1];
u3(-3.06887809303977,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.99907737377714,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.02889097551721,2.87717963468234,-0.221246631853933) q[1];
u3(1.94011596940606,-3.64357470422918,-2.05532988114199) q[3];
u3(1.14989691696026,1.67404915173255,-2.69077284461745) q[1];
u3(0.816188638744678,-3.31319856470732,2.72034817963548) q[3];
cx q[3],q[1];
u1(0.445949931791363) q[1];
u3(-1.58602437989831,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.14692588223445,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.71360731942969,-1.28343118448927,-0.479249230059774) q[1];
u3(2.35991422144443,-5.61876173832239,-0.434064392481052) q[3];
u3(1.73701858324212,0.875543489302492,-2.06608308886006) q[4];
u3(1.00014968469129,1.32028146189030,-4.55254164915989) q[0];
cx q[0],q[4];
u1(2.80395883160557) q[4];
u3(-1.69557173599734,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.28493810846766,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.83632335895953,-4.35671975535135,1.65310126555554) q[4];
u3(1.41588026395746,-2.73618507546512,1.70722823542122) q[0];
u3(2.16592563339466,-1.19372072111794,1.21475139853195) q[1];
u3(1.44135453112557,-1.84122955778372,-0.446575587588288) q[0];
cx q[0],q[1];
u1(0.624518051647327) q[1];
u3(-3.37935535564133,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.93408850289735,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.34272032801431,1.13450113939047,-4.38042831852700) q[1];
u3(0.387885552540009,-0.845995337440985,0.135033290163622) q[0];
u3(2.72640032446709,0.457959217853523,-1.58935304380277) q[4];
u3(2.43076558356475,5.11178448226800,0.242060915036593) q[3];
cx q[3],q[4];
u1(2.14574901005594) q[4];
u3(-1.75775467510004,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.06348480681176,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.620681916306794,-2.91414395491163,1.52996909610277) q[4];
u3(1.28172450779004,-3.09956197202902,-1.86509952973646) q[3];
u3(0.899350599167070,2.29359115543426,-0.805239651559413) q[1];
u3(1.29656963468804,1.65151914016671,-1.11125978094921) q[2];
cx q[2],q[1];
u1(1.38821967490523) q[1];
u3(-0.771242955647163,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.45175276339709,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.63194507060999,1.98559937577431,-0.420359162993651) q[1];
u3(2.16124146087671,0.108189871279208,4.64557111500968) q[2];
u3(0.391796526223371,0.124520565252006,0.839774188162484) q[0];
u3(0.959561312872230,0.899061142182044,-2.75728174484414) q[3];
cx q[3],q[0];
u1(0.690265208885407) q[0];
u3(-1.44231754610843,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.0594240048906489,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.403619528322979,2.45698087627853,-2.94129221037568) q[0];
u3(1.64700978860684,1.28049503199524,2.63289425887429) q[3];
u3(1.96371760699395,4.04902617908296,-1.11194037159616) q[4];
u3(2.10104853082455,1.68556964259262,-0.842450555066325) q[0];
cx q[0],q[4];
u1(1.59650118130156) q[4];
u3(0.191070999042714,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.01186614897564,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.38123022249922,-0.0666392805078468,-2.20294908771557) q[4];
u3(0.369627701137895,-3.31250695993147,0.894885155172824) q[0];
u3(1.77405586886135,-0.416990388689469,2.67192501378838) q[1];
u3(1.54591238222284,-0.625771120881985,-1.56572426096065) q[3];
cx q[3],q[1];
u1(0.321539431191913) q[1];
u3(-1.14267278749621,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.21633358115397,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.03344456461269,-1.68026689550377,3.63599680768303) q[1];
u3(2.37294843972114,-1.60626001668937,-2.34887155942879) q[3];
u3(2.16631067204185,1.44743141665457,-0.195574627996765) q[4];
u3(0.731502499696338,-0.332403957619410,-2.70736221476755) q[2];
cx q[2],q[4];
u1(1.73492253857468) q[4];
u3(-2.16687942210507,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.206532208588780,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.71463542986320,4.34920415361298,-1.50570334009549) q[4];
u3(2.28789256872700,-6.14519476117704,-0.00779388268609083) q[2];
u3(1.62281670037201,-1.57926159523109,-1.30828745867218) q[3];
u3(0.427535325332211,-4.48356811738415,1.07895351945276) q[0];
cx q[0],q[3];
u1(3.38399033783558) q[3];
u3(-3.88149505427813,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.221556064644457,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.28443999336286,0.256100346021565,-3.59652296290045) q[3];
u3(1.51211363930395,0.412218788654326,2.31882778429128) q[0];
u3(0.843011164275295,2.51013191856631,-0.702436244076252) q[2];
u3(1.16080372619590,-0.0621031764868816,-2.73328685684808) q[3];
cx q[3],q[2];
u1(3.21724028379497) q[2];
u3(-2.24766437452152,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.02231893374492,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.54808256759525,1.02368641510881,1.18237576304494) q[2];
u3(2.19991701393122,1.96370839125504,-0.0750685151411516) q[3];
u3(1.15091033385386,-0.528782031577280,-0.529870077367709) q[0];
u3(1.27495266324350,-2.97212856321914,-0.240375262213272) q[1];
cx q[1],q[0];
u1(2.57386791105244) q[0];
u3(-2.11150624223595,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.312561369339654,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.69825052801775,2.77697727433592,-1.78352517517899) q[0];
u3(1.83997307381447,1.30324473969489,-4.27390948364069) q[1];
u3(1.54805730455577,1.92629385241952,-4.15465548252706) q[2];
u3(2.66967796501669,3.49720653229279,-2.74895241343971) q[3];
cx q[3],q[2];
u1(0.0137393235965224) q[2];
u3(-1.13081210225856,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.18136959932407,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.59064678249926,2.46784930173964,-0.377233816789886) q[2];
u3(1.05708197004755,0.183764993345686,-3.57091349603642) q[3];
u3(1.54118805142322,0.697445696104572,1.68065613536374) q[1];
u3(0.988409135762518,-1.74786017809630,-2.03085839930929) q[4];
cx q[4],q[1];
u1(3.01500678612260) q[1];
u3(-1.83374908183347,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.17871727580094,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.31416769128930,0.626958610694419,1.70357380978759) q[1];
u3(1.79429812376768,-1.33545297914753,4.87502300030377) q[4];
u3(1.41947275826804,-0.0190558280333480,-1.70614224838940) q[1];
u3(1.22667237836129,-3.55346775308127,1.82143264605586) q[2];
cx q[2],q[1];
u1(0.955190593268241) q[1];
u3(-3.38916454381059,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.01347380787294,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.31893422163571,4.62363334535328,-1.11511488430655) q[1];
u3(1.64238213343701,-2.84606889554183,2.25032594574007) q[2];
u3(0.952230751479619,0.529444874283653,-1.03268295388736) q[0];
u3(0.746285676494418,-1.49075618660672,0.579496788367721) q[4];
cx q[4],q[0];
u1(3.03273947500573) q[0];
u3(-1.58043005113484,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.934371595932496,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.488308821291939,-0.306404850162162,-0.181310520574726) q[0];
u3(1.82503977236363,1.97890922636259,0.407357638307726) q[4];
u3(2.95781022952644,1.77779113202812,-2.21062170102464) q[4];
u3(1.53680933779069,2.39597630016176,-3.14817351516130) q[2];
cx q[2],q[4];
u1(3.11318505643624) q[4];
u3(-1.61714633120901,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.741422715709961,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.57369554533467,1.39921058916174,-3.64501446192969) q[4];
u3(2.59566600623750,1.22668996380779,4.96965534355688) q[2];
u3(1.86763037910597,-4.31675345591855,1.53837456037547) q[3];
u3(0.507900710087922,-1.17959985987741,2.83126926238346) q[0];
cx q[0],q[3];
u1(1.05906182598719) q[3];
u3(-0.0586640999052932,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.39026026245814,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.94516585562691,-2.66882804030925,2.22399480091499) q[3];
u3(1.95396834322333,-4.01275825504472,-0.496813581234322) q[0];
u3(1.22446787044322,-1.11708527234968,1.73260012053608) q[4];
u3(0.187195637528102,1.33490452751337,-2.68758203857296) q[2];
cx q[2],q[4];
u1(1.21994549264723) q[4];
u3(-0.570820023705194,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.0560609461943000,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.15894488992518,0.767563773808949,-4.91600068038894) q[4];
u3(2.00863552212012,1.98617146480221,1.32338641136730) q[2];
u3(0.448405562948718,2.09173405923234,-2.43875839328142) q[0];
u3(0.402136996369842,-4.21613453334387,1.67343606174789) q[3];
cx q[3],q[0];
u1(0.452085503998845) q[0];
u3(-1.62725825655452,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.90864865543179,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.79584199676152,1.60144811167114,-1.79099609419403) q[0];
u3(1.65480722022292,5.72125557373792,0.182788969741424) q[3];
u3(1.17206045060626,-1.62605348090514,2.46834498072732) q[3];
u3(0.678451506440524,-2.20771114839448,1.30907588923279) q[4];
cx q[4],q[3];
u1(-0.0309529624454163) q[3];
u3(-1.63730951595087,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.835705944310917,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.39547512071178,-1.84321432556536,-1.05554574916970) q[3];
u3(0.742962077829541,-2.47519333181856,-0.912002991688486) q[4];
u3(1.87840496780036,2.24194738867809,-3.76810077471997) q[0];
u3(1.25202813086344,3.58583668165547,-2.41923048130782) q[2];
cx q[2],q[0];
u1(-0.230902082509685) q[0];
u3(-2.34361959223975,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.21033559216256,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.787344792234166,-1.95200531841712,1.87556019975194) q[0];
u3(0.777445627778936,0.630296608828869,1.41749806742338) q[2];
u3(1.95305670297495,3.33826396802631,-2.69231721681988) q[1];
u3(0.100319419791895,-2.29586892734054,3.74331892542372) q[0];
cx q[0],q[1];
u1(1.15843933171928) q[1];
u3(-0.478926786575838,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.53232833477674,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.897923161695381,1.98168891842125,-4.01725008181818) q[1];
u3(2.90353511101042,0.379260645575474,1.54282985088743) q[0];
u3(1.25779922407141,-0.424901090170599,2.31455515503253) q[2];
u3(0.888092220055796,-2.06029078048921,-1.78610875019017) q[4];
cx q[4],q[2];
u1(-0.459646980614116) q[2];
u3(1.20021622782270,0.0,0.0) q[4];
cx q[2],q[4];
u3(3.81867506911646,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.354232320102802,0.516486777669684,-0.389628205449326) q[2];
u3(0.645186263353638,-3.42596961276824,2.73677848375907) q[4];
u3(1.02502000423083,1.08935429622607,-0.368911920327431) q[2];
u3(0.729240415405390,-1.03993419697607,-2.35114108695457) q[0];
cx q[0],q[2];
u1(3.06539381400880) q[2];
u3(-2.50660750996926,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.16217377746817,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.64993765207346,1.47182269923789,2.53127345657271) q[2];
u3(1.83639295032957,-1.80049153310890,-0.295901237502845) q[0];
u3(1.89331130375763,0.126638441579571,2.06342074665454) q[3];
u3(1.74602278463195,-3.10794536732343,-2.74063928138542) q[1];
cx q[1],q[3];
u1(-0.105564249597276) q[3];
u3(-1.66503511434379,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.990934118439850,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.22465733959296,-2.53727126556806,2.51660917103683) q[3];
u3(1.97372303909645,1.63278346412072,-0.953886379915530) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];