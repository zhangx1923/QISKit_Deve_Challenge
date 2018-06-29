OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.11315788613728,4.08037031687743,-1.27737423704144) q[6];
u3(0.954702787472467,0.940049957913789,-0.639655603497493) q[5];
cx q[5],q[6];
u1(1.58611795502217) q[6];
u3(-3.22430759553895,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.75939386226172,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.72550375408428,0.820463006015724,2.21919235264145) q[6];
u3(2.39694323142966,1.65850363245810,1.53912584853053) q[5];
u3(1.16636560981009,2.23443837184863,-0.639799189449779) q[4];
u3(1.30774095574041,1.03228385747599,-0.680225162593625) q[3];
cx q[3],q[4];
u1(3.70524078810795) q[4];
u3(-4.23192519834225,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.110520971877965,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.63481511465006,-2.46812031317326,3.21196338329230) q[4];
u3(2.61526964279526,-1.06754778899309,-4.29865981866055) q[3];
u3(1.04816882727951,0.605642682298069,-1.25480182765137) q[0];
u3(1.32878871385752,-0.684953597042015,-0.160619336279697) q[7];
cx q[7],q[0];
u1(2.74787338369405) q[0];
u3(-1.94985333667210,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.70143499842576,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.55956560834616,-0.859732423493600,0.978804684227678) q[0];
u3(0.967389537013145,-2.35247706427502,1.42754794001370) q[7];
u3(0.521499114608616,1.92047637345061,-0.846894116091978) q[2];
u3(1.25759893655775,0.0631957479298666,-2.21290605544969) q[1];
cx q[1],q[2];
u1(2.37380436834571) q[2];
u3(-2.95247843332403,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.17365795613645,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.07641556728758,0.389906047877832,-3.47972107748864) q[2];
u3(2.34939787800068,-2.08523337577506,3.51612232709588) q[1];
u3(1.72280741462977,0.187246947786468,1.08235694090833) q[7];
u3(1.56133106671830,-2.70355057220809,-1.85982533405320) q[4];
cx q[4],q[7];
u1(0.489992445037229) q[7];
u3(-1.63980963287277,0.0,0.0) q[4];
cx q[7],q[4];
u3(-0.226964193321291,0.0,0.0) q[4];
cx q[4],q[7];
u3(1.32163648442620,-1.34180080796992,-0.631158099167971) q[7];
u3(0.811290904659027,1.94457323095571,4.15504875172668) q[4];
u3(1.46222352686171,2.17310969039243,-3.23695408410908) q[1];
u3(2.26344527234346,-3.15073658732644,2.83014088126207) q[3];
cx q[3],q[1];
u1(-0.0597456279103741) q[1];
u3(-2.03877699100437,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.56799726886761,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.71617671219672,2.98898217246321,-2.92160028017500) q[1];
u3(2.32296711917095,-1.75069012884607,-3.51745964576961) q[3];
u3(1.63418162223240,2.33105821280114,-2.83959302079261) q[6];
u3(0.689086357575282,2.91098999405990,-2.93328142283232) q[0];
cx q[0],q[6];
u1(2.25467056155505) q[6];
u3(-1.82937242499880,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.98951242228247,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.04901307339525,1.91304636445713,1.12465379689842) q[6];
u3(0.206219544458968,-2.36718273171789,1.86078310447293) q[0];
u3(1.97835674810402,0.760464994267119,-3.42957705204267) q[5];
u3(1.59263392979201,2.44154248557982,-3.79184695576607) q[2];
cx q[2],q[5];
u1(2.76954985528427) q[5];
u3(-2.14438332533991,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.524965906130898,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.699792996691529,-1.07061978863809,4.26851881212294) q[5];
u3(1.14671605832492,-4.48716491268484,0.183746509035260) q[2];
u3(1.75465585373259,1.34815429529810,0.808124261973534) q[1];
u3(0.722711146155469,-0.292964995204912,-3.30708158390103) q[5];
cx q[5],q[1];
u1(0.250942144364048) q[1];
u3(-1.31690700448929,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.61909175530247,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.79017271280355,3.67458890817981,-2.16215422498426) q[1];
u3(1.22385835927159,0.458647204407690,4.16093480043519) q[5];
u3(2.85216431307512,-0.589213136681822,2.33712722687518) q[6];
u3(2.04009277049579,-2.30458763071061,-1.84329036746274) q[3];
cx q[3],q[6];
u1(-0.369170508769481) q[6];
u3(1.36400544469283,0.0,0.0) q[3];
cx q[6],q[3];
u3(3.37707824395425,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.718826038326125,0.287653886033905,-1.97091173161196) q[6];
u3(2.54540860034033,3.70554339808573,1.90550912041240) q[3];
u3(2.18886325704370,1.58407973307357,1.40566756145294) q[0];
u3(0.571980169590041,-0.405055850626300,-2.94948156276594) q[4];
cx q[4],q[0];
u1(3.05545431373538) q[0];
u3(-1.80344966543470,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.593204446611107,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.78109253230092,-1.98340624827819,3.15133239328175) q[0];
u3(1.04604357451565,4.54449321229555,1.30749007212096) q[4];
u3(0.387941571618221,2.51072780098218,-1.79834034670459) q[2];
u3(0.430098687649644,0.597672777353556,-2.47305185306145) q[7];
cx q[7],q[2];
u1(3.45414074600579) q[2];
u3(-1.60926911151919,0.0,0.0) q[7];
cx q[2],q[7];
u3(2.19739487407811,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.09586810946650,-2.04314400069695,-0.283451299907183) q[2];
u3(1.14249070483655,2.92702509575642,-1.11713848067843) q[7];
u3(2.15090204256054,0.881934922438585,1.93151894182972) q[7];
u3(1.77637489073223,-0.896015835184664,-0.908388058564329) q[2];
cx q[2],q[7];
u1(1.72080366878394) q[7];
u3(-2.45569065535096,0.0,0.0) q[2];
cx q[7],q[2];
u3(3.50675501360177,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.33117593429582,2.06416439169968,-3.36372280534399) q[7];
u3(0.382975862250862,-3.94469179130017,-1.72735051211252) q[2];
u3(0.298652381628917,0.640224817142693,-1.34180089213362) q[1];
u3(0.673663114010879,-1.18395952044498,-0.433999607687753) q[4];
cx q[4],q[1];
u1(1.84328056382523) q[1];
u3(-3.13524660970186,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.674556304823162,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.06562892046750,-2.00104945068297,-2.16050378660866) q[1];
u3(1.09425340895875,-1.46018679581565,1.66806086201654) q[4];
u3(0.810055332663621,-1.54081343037942,0.480390896153037) q[5];
u3(0.805289634472295,-1.86503978805715,-0.308362321737373) q[6];
cx q[6],q[5];
u1(0.113302114222211) q[5];
u3(-1.45398372944037,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.56923250192213,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.64525688958629,4.16450645581016,-1.95164841647840) q[5];
u3(1.07258138101949,1.19953389398343,-2.80825369913496) q[6];
u3(1.81371962448890,-1.34665124643813,1.02580160631816) q[3];
u3(2.53915452999246,-1.06723015423133,-1.30439621358224) q[0];
cx q[0],q[3];
u1(1.84069389295309) q[3];
u3(-2.92023739439483,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.0779130162735375,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.75107712990288,3.14141496989825,-3.11920332352136) q[3];
u3(0.464105718177902,-1.19870366708020,-3.56698195896585) q[0];
u3(1.23669726937617,3.24402614072058,-1.58848983158712) q[1];
u3(0.791032662322240,1.07412760221163,-2.51322859683993) q[5];
cx q[5],q[1];
u1(3.09175053810663) q[1];
u3(-2.00735594257144,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.458614210990714,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.61062351742145,-0.187200047664955,0.905738795123618) q[1];
u3(1.12417049522345,-2.75494928268388,3.02638748760695) q[5];
u3(1.98171089870364,0.0443680480122406,0.738368583487585) q[2];
u3(0.397682774672997,-2.04051761851361,-2.25781419004150) q[6];
cx q[6],q[2];
u1(-1.29303352972557) q[2];
u3(0.879963739431227,0.0,0.0) q[6];
cx q[2],q[6];
u3(3.83758505436717,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.23960853172276,-0.796042054788433,-0.819005244120973) q[2];
u3(1.30632185106419,2.73621943445966,1.94565831094896) q[6];
u3(0.816504589208377,1.59863732336206,-0.0942061648216259) q[4];
u3(1.48114252681960,0.559266963481844,-2.02923723563458) q[7];
cx q[7],q[4];
u1(3.55452553823679) q[4];
u3(-4.17575712363804,0.0,0.0) q[7];
cx q[4],q[7];
u3(-0.861386418826864,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.77401194581274,0.650743317596566,-3.58359709733182) q[4];
u3(1.37801162492887,-0.541566286897446,5.60999332350236) q[7];
u3(1.13019478475926,-0.241153300572488,2.77483031681632) q[0];
u3(1.41104040088196,-2.22466083883258,-1.76945097598905) q[3];
cx q[3],q[0];
u1(1.03646094325132) q[0];
u3(-0.499607104743433,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.89056812005688,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.99265383341101,0.507931447758099,-0.807087478962185) q[0];
u3(2.12409508991927,5.51986852201393,0.616420258144929) q[3];
u3(1.94020588015709,1.70521579250859,-0.509481664083218) q[7];
u3(1.62753789975807,0.223455357039740,-2.24853406219500) q[2];
cx q[2],q[7];
u1(-0.701195730143793) q[7];
u3(-1.72417465050202,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.870711775221429,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.14967370262691,-0.280492390479453,2.67312674947985) q[7];
u3(0.744525797322325,-1.77478401814385,-0.760222412231993) q[2];
u3(1.47049254143829,0.629604942736281,-3.34781278597813) q[4];
u3(1.22086173199523,-2.67009631285729,2.36112418825695) q[6];
cx q[6],q[4];
u1(2.88948784049781) q[4];
u3(-2.27480898882423,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.13680010339981,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.41062859153721,-2.39059795749741,2.59369779042867) q[4];
u3(2.30472518795158,-3.63735770529620,-0.0308696665259902) q[6];
u3(1.96938251729313,3.12679721591305,-1.07879326957291) q[5];
u3(1.09433350008503,0.216660528986364,-0.558534892226506) q[0];
cx q[0],q[5];
u1(3.66772657157524) q[5];
u3(-1.50634901514710,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.19549116335530,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.977993481618097,-1.09910660118319,0.0943983158787222) q[5];
u3(2.42785184350693,-1.65030794308936,-2.12959144705162) q[0];
u3(1.64171164810802,3.73381981105344,-1.17287644697074) q[1];
u3(1.92239095504822,1.89397888033601,-1.37204294012717) q[3];
cx q[3],q[1];
u1(1.68737651842342) q[1];
u3(-2.14881828588812,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.75258793715334,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.277494355091821,-1.28671974645927,1.88337017649334) q[1];
u3(1.76953411736672,-0.283004923799374,-4.84494084750756) q[3];
u3(0.809648244238913,1.52893291259681,-0.835645463085020) q[4];
u3(1.06915777448097,0.419463127695955,-3.30746759869741) q[3];
cx q[3],q[4];
u1(-0.117930495117049) q[4];
u3(-1.32724000155191,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.61844695745158,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.66883834157107,-1.11886595814607,0.934370676220927) q[4];
u3(1.38544330760695,-5.76077215867042,0.223648073373379) q[3];
u3(1.75205446560619,2.12428160238746,-2.34889765826771) q[2];
u3(1.26376613315397,-2.41796334660721,2.81471647746823) q[0];
cx q[0],q[2];
u1(2.31695135031341) q[2];
u3(-1.58711357385666,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.187962704940586,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.21295804506081,3.85895766633109,-1.54642361049315) q[2];
u3(0.967762138141623,-3.37508136945028,2.85157758302951) q[0];
u3(0.645920704169678,-2.25062658998029,0.475811524864761) q[6];
u3(1.94329326126167,-3.69383372570273,0.453545489638529) q[1];
cx q[1],q[6];
u1(1.09105522238797) q[6];
u3(-0.0623942772796344,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.86214812000070,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.50033021787522,0.509854381160811,2.57354450464475) q[6];
u3(0.673765633316101,4.56349327155721,0.748676232233383) q[1];
u3(1.87480856807043,2.49128740765927,-0.924942192560219) q[7];
u3(2.86864219944694,0.162632316312408,-4.39001048258859) q[5];
cx q[5],q[7];
u1(2.11542345042798) q[7];
u3(-2.98353385223227,0.0,0.0) q[5];
cx q[7],q[5];
u3(0.869156398875644,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.69570579004792,0.456351953592577,3.10641012601749) q[7];
u3(2.56046299659983,5.76475005175258,0.0627450112184080) q[5];
u3(0.552159709457521,-1.76677898466436,2.05141379547655) q[7];
u3(0.710263170001706,-3.09667919778746,1.49113601752261) q[6];
cx q[6],q[7];
u1(1.68388098375535) q[7];
u3(-2.54990864348410,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.270998766512230,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.53240724304555,0.416179661983092,-1.88153314754886) q[7];
u3(2.09018415072023,2.43013185103381,2.56686254755123) q[6];
u3(2.08134083039251,0.572269257909883,-3.20297976012569) q[5];
u3(2.58493824951601,4.80196900855269,-0.731574716838484) q[4];
cx q[4],q[5];
u1(1.56947615858165) q[5];
u3(-2.31326206221190,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.699943507025539,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.05808114560910,0.920704195660029,0.128578262537050) q[5];
u3(1.92664136596301,1.02459902680610,2.47854115558301) q[4];
u3(0.338922243732025,1.56484549487559,-1.65914442210701) q[2];
u3(0.827685263768075,-3.34601449868281,1.15939846085834) q[1];
cx q[1],q[2];
u1(2.83078924737416) q[2];
u3(-2.30651912716982,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.16062996915570,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.21655657706897,0.593115378848169,-0.237740894729487) q[2];
u3(1.84354220978336,2.36475711135919,0.328539717171579) q[1];
u3(1.81259483620013,-0.671969329155874,1.84407493852086) q[0];
u3(2.05606302097376,-1.55960883038028,-0.243494988264961) q[3];
cx q[3],q[0];
u1(1.67725146511027) q[0];
u3(-2.24645414045226,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.287523277553189,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.76504568888362,-2.43759538598013,-0.469841173524605) q[0];
u3(2.24052041446359,0.0604414209495459,0.856230882970217) q[3];
u3(1.40125165904929,-0.437692741981313,1.59487121064108) q[6];
u3(1.60973126313951,-1.21421862246560,-1.69150383627694) q[1];
cx q[1],q[6];
u1(1.95088338132983) q[6];
u3(-0.618006327277945,0.0,0.0) q[1];
cx q[6],q[1];
u3(3.15388311567407,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.94594702016136,-3.20013797688377,0.696166876027775) q[6];
u3(0.913697952012617,-0.361156804824067,3.87752872121895) q[1];
u3(1.64353526192160,3.73270373968594,-0.757003356028870) q[3];
u3(1.44312373913227,2.94968452998392,0.329526793721235) q[0];
cx q[0],q[3];
u1(1.89502132605863) q[3];
u3(-2.91270086028385,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.953324980664981,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.20884566463078,-0.903599026682519,1.27712881105255) q[3];
u3(1.46572958306551,2.71511277896794,1.71584859835625) q[0];
u3(1.53465057412561,-0.0253759223582388,-0.895122486009804) q[5];
u3(0.641928284298738,-3.46507336860705,0.444037804360902) q[2];
cx q[2],q[5];
u1(1.71635497508774) q[5];
u3(-2.36739098729045,0.0,0.0) q[2];
cx q[5],q[2];
u3(3.08911617156207,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.31953147684738,-1.90893479034482,3.47280522730059) q[5];
u3(1.67913910524996,-1.48684869724745,-1.96295666878262) q[2];
u3(1.44789771181015,0.842113765857206,-3.90325304828394) q[7];
u3(2.75889570859150,2.64858341158483,-2.75273018345692) q[4];
cx q[4],q[7];
u1(-0.810836406677257) q[7];
u3(0.303384698106256,0.0,0.0) q[4];
cx q[7],q[4];
u3(4.20325040266672,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.05309351071266,1.36962802255920,-0.429756559478158) q[7];
u3(0.995411328035113,1.52708828042214,3.17106052903273) q[4];
u3(2.24341503491142,3.24656681043205,-1.72883394796895) q[5];
u3(1.14910826428139,1.36793447009036,-0.972512413358904) q[6];
cx q[6],q[5];
u1(-0.920240334320052) q[5];
u3(0.0877452073152443,0.0,0.0) q[6];
cx q[5],q[6];
u3(3.40633481827861,0.0,0.0) q[6];
cx q[6],q[5];
u3(2.06511235127487,1.12343246650729,2.44860375802885) q[5];
u3(1.22089847982777,1.19467561439666,-0.316448413557593) q[6];
u3(1.68504611453995,-0.604442140799463,0.841046501859864) q[2];
u3(1.25127134291665,-1.51941446923060,-2.32191007948347) q[4];
cx q[4],q[2];
u1(1.69859222144452) q[2];
u3(-0.278360098857945,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.49620053053984,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.69409261610437,-2.92771363349804,1.34209884733839) q[2];
u3(0.985214107476794,0.296868313492461,5.09355726093780) q[4];
u3(1.52725252299988,-1.44480623651691,0.731841320859072) q[3];
u3(1.63823579710967,-3.05050994784799,-0.165807480945180) q[0];
cx q[0],q[3];
u1(-0.712408851032055) q[3];
u3(1.23176080960806,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.61139377858147,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.32241594022431,2.20358468388164,0.638813138015665) q[3];
u3(1.69000741496456,-0.433899957210033,0.0541731848745768) q[0];
u3(0.877983897279841,-0.0717156614411261,1.84231350417583) q[7];
u3(1.17464102661729,-1.31775378460617,-0.984456717368698) q[1];
cx q[1],q[7];
u1(2.90449379731919) q[7];
u3(-1.87574722563044,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.37431136129291,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.26313632989795,-0.978408931707944,-3.24128175384188) q[7];
u3(1.35725473816575,-2.22285840890430,3.43707126781343) q[1];
u3(2.62111517009932,-4.55463936489999,1.46398740679292) q[0];
u3(0.997359398376366,3.45723773285947,-1.87440087197451) q[4];
cx q[4],q[0];
u1(2.99125465094865) q[0];
u3(-1.41572084892804,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.19964591985418,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.69728592265606,-2.81249762881144,0.0189145506569228) q[0];
u3(2.81182748669803,3.52275846714785,-2.06553096021901) q[4];
u3(0.616094452818310,2.24090543677667,-3.53648506491592) q[5];
u3(0.690020689769036,1.07318124245618,-2.83733561397291) q[6];
cx q[6],q[5];
u1(1.59952196444544) q[5];
u3(-0.248517996833509,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.03750718141234,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.350812798853490,4.77343831424427,-0.326177580856278) q[5];
u3(2.37673633991211,2.45999677015982,3.11203460330072) q[6];
u3(0.103844437364760,1.23040761740755,-3.01311197278478) q[1];
u3(1.49237037091311,-3.29357744196024,2.67430879726259) q[3];
cx q[3],q[1];
u1(0.438429217992448) q[1];
u3(-1.53468121640409,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.33491733378786,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.498544812171138,-0.0236107396911658,-0.621548139528362) q[1];
u3(0.884891559854943,1.84952090052077,0.846926739655168) q[3];
u3(2.27109672188227,-0.501351392665398,2.57981504041054) q[7];
u3(2.45047378725188,-1.20448879031987,-1.49319008741720) q[2];
cx q[2],q[7];
u1(1.21872167340702) q[7];
u3(-0.630925463761722,0.0,0.0) q[2];
cx q[7],q[2];
u3(-0.114969807950990,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.23086072615901,0.597208461415206,1.98109267134172) q[7];
u3(1.48968356906351,-2.25785652355394,-0.630644109573342) q[2];
u3(1.44015693888776,1.75980323790948,-3.43462137396130) q[6];
u3(0.823031573395561,-3.08338734445101,3.11973628442104) q[0];
cx q[0],q[6];
u1(2.83523749966715) q[6];
u3(-1.85785621197935,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.32233725100556,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.21874944569978,1.49995089445509,-0.965539146496146) q[6];
u3(1.11843568717679,-1.10197031433619,-1.67934166922395) q[0];
u3(0.854605742195660,1.49217190268045,-3.05430756457646) q[5];
u3(1.71323326677270,-3.52862282838089,2.63252075375735) q[4];
cx q[4],q[5];
u1(1.43466255679762) q[5];
u3(-2.68117114006054,0.0,0.0) q[4];
cx q[5],q[4];
u3(-0.116117746895439,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.88177602314348,-1.43878058725115,-1.55141588021599) q[5];
u3(0.289619521791367,1.76195470383658,4.30028383153865) q[4];
u3(0.969244673311439,2.02638455207667,-0.478354076434255) q[7];
u3(1.22475628165011,-0.247730173062098,-4.30864171342045) q[3];
cx q[3],q[7];
u1(0.777669527265298) q[7];
u3(-3.56913559606383,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.43706477556775,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.43078863675778,0.968219106206628,2.92517048555820) q[7];
u3(2.42004768427995,-1.22105314623998,-0.738758946833411) q[3];
u3(1.91861259489832,-1.14764742541254,0.625185225417996) q[2];
u3(1.97761455486218,-1.83828724863352,-1.27722253350096) q[1];
cx q[1],q[2];
u1(2.44476824999988) q[2];
u3(0.276080104504533,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.45996546993029,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.58402950575975,1.40099559111822,-1.40827322111917) q[2];
u3(1.46439585545200,1.59444530512176,-0.0154765726187495) q[1];
u3(1.63083195991763,-0.668264755734686,0.180646543482910) q[5];
u3(2.05102674260155,-1.23078663401082,-1.75568008007590) q[6];
cx q[6],q[5];
u1(0.179020604673410) q[5];
u3(-2.36731551232517,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.18399016375917,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.915481549146197,-2.36053038271822,3.40558179967351) q[5];
u3(1.05918514473300,2.54300729626012,-0.801111729882406) q[6];
u3(0.664349435748412,2.21270630408526,0.332803399474596) q[2];
u3(1.65515423004574,0.566812563773575,-3.99414714534288) q[4];
cx q[4],q[2];
u1(0.0715551897688362) q[2];
u3(-1.37714451459962,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.82583474085711,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.17501320945565,-4.16918843107132,0.798769815374849) q[2];
u3(2.05947390516072,-4.10679295696023,-0.661533545716787) q[4];
u3(2.27938502567555,0.616921865638489,1.06620442673733) q[0];
u3(0.811053696512704,-3.98046292831527,-0.920671693902491) q[3];
cx q[3],q[0];
u1(1.09780848889175) q[0];
u3(-0.180492686222524,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.48718917782254,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.44398962186662,0.326728086384702,-2.13195593346685) q[0];
u3(0.964909698255623,-3.53359215544108,-1.86247402337867) q[3];
u3(2.52839738966990,-2.37099443237585,0.609091585777295) q[1];
u3(2.66965628920596,1.17515950707288,1.80989290609247) q[7];
cx q[7],q[1];
u1(1.20332207296002) q[1];
u3(-0.274456017723909,0.0,0.0) q[7];
cx q[1],q[7];
u3(2.78256497494189,0.0,0.0) q[7];
cx q[7],q[1];
u3(0.804676130754381,-2.79963864538994,2.59334494371450) q[1];
u3(2.68102024601501,2.62547601163404,0.803625050432255) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
