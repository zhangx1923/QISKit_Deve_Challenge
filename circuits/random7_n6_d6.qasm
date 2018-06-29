OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.56291280366094,0.347291573478041,0.258306674862317) q[2];
u3(2.49839355671112,-0.535278856657502,-3.98640040343303) q[1];
cx q[1],q[2];
u1(-0.170475031500005) q[2];
u3(0.639697324287470,0.0,0.0) q[1];
cx q[2],q[1];
u3(4.03562500811738,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.90404509904624,-0.188365865701045,1.42829627152526) q[2];
u3(2.84491418943231,1.59028763831783,1.38287844126871) q[1];
u3(0.704000843061147,-1.42329928458280,2.13814857068918) q[4];
u3(0.390594164483226,-0.405096076788089,-2.05927288268994) q[3];
cx q[3],q[4];
u1(1.83910211363507) q[4];
u3(-0.0149708402073097,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.946969611332453,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.04971823593004,0.0990122560774060,-2.00805827239910) q[4];
u3(1.58216390946561,0.782248513522372,-3.75803210194535) q[3];
u3(1.43688003202904,0.972404121122051,-2.81593789327986) q[5];
u3(2.53502328799020,2.88670203098048,-3.06549290338611) q[0];
cx q[0],q[5];
u1(3.12127741850362) q[5];
u3(-1.81426893974815,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.66738873835402,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.34358500124867,1.75932747925689,-3.57506797119942) q[5];
u3(1.37038940503343,1.82459455638194,-3.63635940625609) q[0];
u3(0.800300281744879,-0.697201958098042,0.466707591760426) q[0];
u3(1.00646766431963,-3.13081802783361,-0.204266739464449) q[5];
cx q[5],q[0];
u1(4.29463458099054) q[0];
u3(-3.20126517632973,0.0,0.0) q[5];
cx q[0],q[5];
u3(-0.256057633351104,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.86490406127279,3.53286294153506,-0.629097399034559) q[0];
u3(1.57165151431700,4.51394324224788,-1.15826087015125) q[5];
u3(1.45133227388789,1.40605130957973,0.134507796937742) q[4];
u3(0.991483032732528,0.602956949373956,-4.56170726682881) q[3];
cx q[3],q[4];
u1(2.56453161848441) q[4];
u3(-2.01648740674878,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.911940360143654,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.47384967803441,-1.58276304537426,1.03537443991239) q[4];
u3(1.21774109201334,2.35804612202443,3.16714611711493) q[3];
u3(2.46941856522707,2.69875792935393,-3.47314441777824) q[2];
u3(1.23595436449901,0.524193641025632,1.26838564476425) q[1];
cx q[1],q[2];
u1(0.467545483091557) q[2];
u3(-0.177642031939364,0.0,0.0) q[1];
cx q[2],q[1];
u3(4.31132053858978,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.65071521741348,2.25387185246160,0.100683486726232) q[2];
u3(1.80504129929823,-2.93125029421562,1.41773456014391) q[1];
u3(2.73729506579249,-2.38669804521694,3.88645119092441) q[1];
u3(0.997128456983195,-0.344778160940970,2.07457138405645) q[5];
cx q[5],q[1];
u1(1.16458130566661) q[1];
u3(-3.60095212087372,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.66493324701473,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.92181948914662,-1.62010768688943,2.79271856869702) q[1];
u3(2.38406351253798,2.11575738532746,1.90914688524056) q[5];
u3(1.70629703285333,0.728357264058260,1.86816348580962) q[3];
u3(1.38202854313105,-0.930839239918138,-0.411815951709478) q[4];
cx q[4],q[3];
u1(1.60752655849629) q[3];
u3(-2.88992468988311,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.926967916843940,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.03721601145041,1.63905298439913,-0.180134411027881) q[3];
u3(2.34601261238752,4.03063561386546,-0.444740922240371) q[4];
u3(2.57225160901153,3.78323097323948,-2.40681542305025) q[0];
u3(1.21480464867648,-1.04248696165025,2.69675753169461) q[2];
cx q[2],q[0];
u1(-0.376878490356953) q[0];
u3(1.00362686195524,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.13721120823286,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.95780237188339,-0.816919260867255,-2.15350637619810) q[0];
u3(1.07798783475622,1.05589813821927,0.0803923441077011) q[2];
u3(1.59799112719045,3.39015669130057,-1.30725038984502) q[1];
u3(2.54998595870629,2.42123827376317,-0.129581780311801) q[4];
cx q[4],q[1];
u1(1.39301718536973) q[1];
u3(-3.39413862403990,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.91231963049853,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.90610370956635,-1.94386996928740,0.429257337729620) q[1];
u3(1.98791334295833,5.07824183341661,0.119571502604459) q[4];
u3(2.35010850689716,0.181809795910071,-1.51451580448030) q[3];
u3(1.61406168726916,-4.09657626960404,1.41517822483471) q[2];
cx q[2],q[3];
u1(2.95451464983893) q[3];
u3(-2.89060510450473,0.0,0.0) q[2];
cx q[3],q[2];
u3(-1.48307729189855,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.688337594124845,2.69998653281684,-1.48207511680629) q[3];
u3(2.07687814379330,-1.71408074479056,2.83193955989695) q[2];
u3(1.32875270055093,1.93346809718544,0.810152357917084) q[5];
u3(1.92679313866059,0.396209722606527,-2.36397408399710) q[0];
cx q[0],q[5];
u1(-0.223255392608951) q[5];
u3(0.909377200045137,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.74514091457645,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.57648908380680,0.0376137967523340,-3.22297454783436) q[5];
u3(1.09981444193108,0.455853387363470,-2.63617207999350) q[0];
u3(1.43154992114963,1.01611678692803,-3.04929502514763) q[2];
u3(0.355085502748218,2.84685706975950,-3.23932887289746) q[5];
cx q[5],q[2];
u1(1.67194809708578) q[2];
u3(-3.18334296617265,0.0,0.0) q[5];
cx q[2],q[5];
u3(0.831887274301363,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.42678497022749,3.62598942516592,-0.409433364185984) q[2];
u3(1.58731252713411,-5.67750483034865,-0.561320045240286) q[5];
u3(0.565502348228958,-0.751291609536084,2.48478838053401) q[4];
u3(1.40525063054832,-1.22202929404745,-2.85491327519658) q[1];
cx q[1],q[4];
u1(1.76954044003442) q[4];
u3(0.0622858120926901,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.534510476734459,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.85590703513011,3.21139921845280,-2.41497918340926) q[4];
u3(0.830445876234489,2.19956270845245,-2.27899274737006) q[1];
u3(0.851153520917638,2.70660898642195,-0.290803465987780) q[0];
u3(1.43877108352633,1.41813289478731,-1.55100131113227) q[3];
cx q[3],q[0];
u1(0.412214661216437) q[0];
u3(-0.846793587014123,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.19867317090848,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.51189295984990,-0.714378516208145,-2.66547530850027) q[0];
u3(1.97883540930748,-1.37212260259681,4.09239864295403) q[3];
u3(0.915246765514898,-0.779604010221170,1.03137395268982) q[3];
u3(0.993859058836318,-3.49145918768776,2.42037713511151) q[5];
cx q[5],q[3];
u1(1.51512415297471) q[3];
u3(-2.99056822188436,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.283917874854419,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.60117857141882,0.970776271907563,-2.84412846697145) q[3];
u3(0.489719427381861,-1.11457882419437,-0.921936741638961) q[5];
u3(1.50379417956441,3.53033368748275,-1.42945300058778) q[4];
u3(2.69217196441803,0.896461140096772,-2.64145397668801) q[1];
cx q[1],q[4];
u1(0.562194249782057) q[4];
u3(-1.30709375697808,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.04460791928545,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.706283429467441,1.48260887971936,-1.93830703888290) q[4];
u3(1.39419269214528,-1.73569513483526,-2.19247718621130) q[1];
u3(1.23008510013598,2.39355895116048,0.382125712297879) q[2];
u3(1.85736074711233,0.375506767987785,-2.61130982535378) q[0];
cx q[0],q[2];
u1(1.02229579674072) q[2];
u3(-0.187491755452389,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.40244010727670,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.75968340945393,-2.02651843187304,1.25845516153405) q[2];
u3(1.13011738883458,3.47770417881052,0.940380404093478) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
