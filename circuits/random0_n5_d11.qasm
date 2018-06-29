OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.68049707592225,-0.756201447370905,0.725924891366394) q[3];
u3(2.35590820804382,-0.753346543936281,-1.84561489565328) q[1];
cx q[1],q[3];
u1(0.567553883615995) q[3];
u3(-1.46155884452322,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.0451638229655718,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.59777634731739,-1.12572145314970,3.26533699647565) q[3];
u3(0.687544301666434,2.13072176814042,2.62289444799619) q[1];
u3(2.21629782251311,1.11662536458005,-2.72713208365690) q[4];
u3(2.46034731758808,0.938513860514959,-4.40399609846126) q[0];
cx q[0],q[4];
u1(1.88834158220931) q[4];
u3(-2.64267236632990,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.257488309428834,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.189749356122133,4.22326896783684,-1.68937375278518) q[4];
u3(1.72851239933309,-0.0612544161273723,1.56814567060825) q[0];
u3(2.69203900529110,-0.548475082419758,3.13219129629184) q[3];
u3(2.53979387350694,-1.81319506281275,-0.427772821353652) q[1];
cx q[1],q[3];
u1(-0.312421339961407) q[3];
u3(-2.48147132338214,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.42808259650650,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.990422962101040,-0.286762303270663,-2.22569471685334) q[3];
u3(1.40341598186806,5.28609499336068,0.650365941164303) q[1];
u3(2.19844874264905,-0.147849321897644,-0.488193036769851) q[0];
u3(0.909140731481330,-0.0809048477573198,-4.44635329955240) q[4];
cx q[4],q[0];
u1(3.19167197619937) q[0];
u3(-1.39470914994574,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.37307181981766,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.636632135180527,-0.380571248536259,-0.518945550640765) q[0];
u3(2.05098102773931,5.81065647279571,-0.401024270287372) q[4];
u3(1.19890554933259,-1.14180056627141,1.44305427749294) q[1];
u3(0.0473498071807912,-1.14965187100339,0.437900452984252) q[3];
cx q[3],q[1];
u1(3.25771275982215) q[1];
u3(-1.55118059133755,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.537376980865067,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.52712342601061,1.22078159943916,0.283559155738003) q[1];
u3(2.16618091795274,3.77707659178584,1.10349092230689) q[3];
u3(2.35631946797039,-1.06401864270097,1.77976839105002) q[0];
u3(1.82071238910975,-1.65231835099353,-1.23802160899421) q[2];
cx q[2],q[0];
u1(2.31204001901936) q[0];
u3(0.341022829785911,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.56772234217238,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.45450335313248,0.898748043680937,-1.29528761344302) q[0];
u3(1.24371327397776,-0.0930286584287439,1.17427438239662) q[2];
u3(1.70427664816075,-2.03177019000281,0.132630797206967) q[1];
u3(1.52598140657495,-2.20408722925184,-0.0586567306660009) q[2];
cx q[2],q[1];
u1(-0.437516317274790) q[1];
u3(0.615295195380994,0.0,0.0) q[2];
cx q[1],q[2];
u3(3.91109886189767,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.14893389157523,-0.596031144519076,-2.19423519564401) q[1];
u3(1.96388860506278,-0.553274106939185,1.15044903723176) q[2];
u3(1.40686013057705,-0.0501696252094328,2.22181882116185) q[3];
u3(1.66357015736696,-2.55211060069004,-1.96715907094596) q[0];
cx q[0],q[3];
u1(1.82646701801377) q[3];
u3(0.0393899489434599,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.619872076372762,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.79173881897885,1.72406629426662,-4.18619145287319) q[3];
u3(0.178678343979038,-1.49443122401415,-3.75228159927670) q[0];
u3(1.10384166297031,1.77312611806740,-4.33431598862131) q[4];
u3(1.15118553550491,2.30613764106424,-2.73528921284150) q[3];
cx q[3],q[4];
u1(1.12967213908910) q[4];
u3(-2.87019905623369,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.71364400439318,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.99767065913094,-0.525394668013177,0.935466963185228) q[4];
u3(1.51333364177883,4.63917876152132,-1.16132630416742) q[3];
u3(1.22225536807462,0.629350426102889,1.68266215583271) q[0];
u3(1.49121768960930,-2.03541836122446,-1.43032907056244) q[1];
cx q[1],q[0];
u1(-0.0241006373707999) q[0];
u3(-2.03256605018027,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.56635604512248,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.42806945323314,1.23521950139520,-4.39123214517921) q[0];
u3(1.47748207498994,2.23690384968936,3.85177326179916) q[1];
u3(1.01453957891002,-2.48416530223334,-0.164782661803171) q[2];
u3(0.975198266060535,-2.78961553205775,0.395460021418789) q[3];
cx q[3],q[2];
u1(2.46327390360363) q[2];
u3(-3.10531972735615,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.30070444084568,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.98012105996664,1.04259331010161,-1.87121009675193) q[2];
u3(0.270171110003300,-1.53057847866140,-0.111220452406067) q[3];
u3(2.20810345152912,2.62675231805576,-1.82480010847154) q[1];
u3(2.35646334395769,1.54239272723284,-0.884339169521802) q[0];
cx q[0],q[1];
u1(-0.588828009317155) q[1];
u3(-2.04780301119224,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.40975980610333,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.45811562626934,-0.276968784551016,-0.733543826491652) q[1];
u3(1.45370646762183,1.56034799821556,-4.26423877897365) q[0];
u3(1.60487817483902,-0.730198670568102,-1.33657537474464) q[0];
u3(2.87741561823037,2.19693868299026,-3.73641406061254) q[1];
cx q[1],q[0];
u1(4.11715631312383) q[0];
u3(-3.68033051363972,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.193850276065126,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.11462425284413,3.99532153242245,-2.04879845282833) q[0];
u3(0.246115920532955,1.21902576050076,4.11627746593102) q[1];
u3(0.789359619721882,2.30635454100929,-0.930966587865570) q[2];
u3(0.197886446568985,-0.0398463890378368,-1.40698896438092) q[3];
cx q[3],q[2];
u1(2.46386894005356) q[2];
u3(-1.79739579879255,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.99658618186064,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.01372369243622,2.70454388766196,-1.82028818643377) q[2];
u3(2.36160952390764,0.207113368296607,0.303800052849766) q[3];
u3(1.75641161684438,3.25464236009330,-0.997936762922185) q[0];
u3(0.808257037339725,1.04973399054981,-1.10408222164714) q[2];
cx q[2],q[0];
u1(0.356101819673751) q[0];
u3(-0.757849834736328,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.12184359152405,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.803635177545233,1.17729864550632,-3.17657943189565) q[0];
u3(1.47452837513638,2.83268845020944,-3.03633824844839) q[2];
u3(2.35139988242637,1.39786442615126,-0.314677723391576) q[3];
u3(1.81904438355291,-0.459902511896783,-4.25351889940883) q[1];
cx q[1],q[3];
u1(1.64129847470895) q[3];
u3(-2.16008144565908,0.0,0.0) q[1];
cx q[3],q[1];
u3(3.61691226528383,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.68690140631882,-0.935850318260342,-1.70572366946765) q[3];
u3(0.940487593767524,2.16161873163717,3.61370399300665) q[1];
u3(0.361736901600982,1.91963894145915,-1.46815456146147) q[2];
u3(0.465570046263706,0.274489208771251,-2.40131243745954) q[4];
cx q[4],q[2];
u1(-0.0824834157842844) q[2];
u3(-1.49402202072852,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.344530063669213,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.38503123109674,-1.32381481241784,-1.60944838996861) q[2];
u3(1.74505002063426,1.85096452767950,-0.998676174305934) q[4];
u3(1.92984600409394,-1.70826683307164,0.596641195536923) q[3];
u3(2.34345085783259,-3.69051939037990,-0.287842258506787) q[0];
cx q[0],q[3];
u1(1.68023981063449) q[3];
u3(0.430628840852480,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.12754746010790,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.647415532310235,-0.796904365184635,-1.72239935603296) q[3];
u3(0.707785143135414,3.99792671511920,-1.98964637934356) q[0];
u3(1.01301007692280,3.35559524944890,-1.64745144495954) q[3];
u3(1.06380691930656,1.58703035793613,-0.995573231869760) q[4];
cx q[4],q[3];
u1(2.52208211744881) q[3];
u3(0.352470690035745,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.79146211945757,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.22385873262144,-2.69421277011063,2.46688580853098) q[3];
u3(0.672287373018175,1.69146664371140,4.12370733709786) q[4];
u3(0.351990269564649,-0.303858353094289,-1.41452831326960) q[2];
u3(1.56189919474216,1.68900122662663,-3.33960506108474) q[0];
cx q[0],q[2];
u1(3.26131577102761) q[2];
u3(-1.32007640619903,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.69991791256474,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.604214389390053,-1.54817366689514,-0.196327095673665) q[2];
u3(0.418085716867297,-0.959562866845742,-5.30481637858708) q[0];
u3(0.766167940467445,2.20493637573240,-1.40243139913450) q[3];
u3(0.860416376371459,0.339913059846281,-1.82439101012565) q[2];
cx q[2],q[3];
u1(0.443549627932197) q[3];
u3(-1.36433669411057,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.04454945692688,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.41795280807115,1.25178532960997,0.197563727876358) q[3];
u3(1.91929105058781,1.72177411834285,-3.76275543841564) q[2];
u3(1.79129556786077,0.603924485917101,1.25071755513190) q[1];
u3(0.248275201771572,-1.44895922689884,-3.51632038345116) q[0];
cx q[0],q[1];
u1(1.08441623435637) q[1];
u3(-0.427773723553032,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.63495248712538,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.722658553916936,3.25045617994434,-1.00876860424411) q[1];
u3(0.836277682007887,0.416132928948714,-3.33495230807290) q[0];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
