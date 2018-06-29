OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(1.03068596057096,1.12261871441205,-2.89642515152095) q[0];
u3(2.27905761334234,-3.21994995855881,2.82742071264069) q[7];
cx q[7],q[0];
u1(1.88427401657598) q[0];
u3(-3.13334709484702,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.749059819976300,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.18831579481671,-0.882572974385406,0.127047834916076) q[0];
u3(2.86648619451285,2.34361648514197,-3.85415037761325) q[7];
u3(1.47613492212378,3.32690330469108,-2.19348652508738) q[6];
u3(1.42330470393497,2.29313045481067,-2.40799514671250) q[5];
cx q[5],q[6];
u1(2.87836047925728) q[6];
u3(-2.07197927960222,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.41354241704635,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.913786331895841,4.50641547611540,-0.824489048117204) q[6];
u3(2.28061303265360,0.284040330374045,5.69439922601176) q[5];
u3(2.06333759561082,-0.363840633121954,0.546176811831292) q[8];
u3(1.82104933291760,-2.57874392785446,-1.80836101281722) q[9];
cx q[9],q[8];
u1(3.59606601375247) q[8];
u3(-4.18301629667964,0.0,0.0) q[9];
cx q[8],q[9];
u3(-0.00866756607560371,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.47879564301847,4.04625018686663,-0.697069254381516) q[8];
u3(0.320334250742322,1.98229585637882,-2.68820229031522) q[9];
u3(1.06678627813684,3.38464218129107,-1.79631878697327) q[3];
u3(0.729110149445112,1.53014444691835,-2.08413592670011) q[4];
cx q[4],q[3];
u1(1.43717305983199) q[3];
u3(-3.65205231911764,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.44808833454564,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.02262601095310,-0.164554595290817,0.715797229218708) q[3];
u3(1.50171920631294,-1.35860452553255,-2.49806016034242) q[4];
u3(0.657770474359129,0.477961208828117,-2.17100214374151) q[1];
u3(1.73316816813915,-2.96721963076636,2.37531059133179) q[2];
cx q[2],q[1];
u1(3.11096778005136) q[1];
u3(-1.96354840335580,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.822399639227105,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.37391622830170,0.441146407734503,0.381811828835173) q[1];
u3(1.46704180009590,1.76204002954753,-2.72681079584395) q[2];
u3(2.71249274928500,-0.253368753534666,2.44247245962705) q[1];
u3(2.57686227711834,-1.28282055184392,1.11377229027569) q[4];
cx q[4],q[1];
u1(0.579125589796275) q[1];
u3(-1.44031723357749,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.15393601972219,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.09654547248869,-2.51676011657241,-0.334576937952231) q[1];
u3(1.09446243835722,4.32091834293051,1.09578105137085) q[4];
u3(2.56406586820750,-0.400949383982688,2.72995146717768) q[6];
u3(2.83535047654112,-2.01203500632588,1.07358159626668) q[8];
cx q[8],q[6];
u1(2.93005201286131) q[6];
u3(-1.53674479829820,0.0,0.0) q[8];
cx q[6],q[8];
u3(0.994590806844312,0.0,0.0) q[8];
cx q[8],q[6];
u3(2.40037623563526,2.20616749519244,-4.06436787258616) q[6];
u3(1.84037956644286,4.57608705599278,-1.50123062648281) q[8];
u3(1.51504382317929,3.27575933588738,-1.23030221932368) q[0];
u3(2.70250027311745,2.80460037023052,-1.58281068344734) q[3];
cx q[3],q[0];
u1(1.28001637951321) q[0];
u3(-0.526996624512607,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.99956933810353,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.05198929757682,1.05878687622280,1.96425887470685) q[0];
u3(0.236430157162876,0.823214542375811,4.35325525879830) q[3];
u3(0.914196499694437,3.05709808860111,-1.93688453881217) q[7];
u3(0.541269388260191,1.87216290402185,-3.06410704235717) q[2];
cx q[2],q[7];
u1(1.12386593017460) q[7];
u3(-0.0780499896154641,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.14677119393169,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.24476418185646,1.58833805804679,-3.32807526396869) q[7];
u3(1.95116261697933,-0.106902645013318,2.87356688624701) q[2];
u3(2.59104532820943,-3.11403355902629,3.04099446957723) q[9];
u3(1.62413497246737,2.59097002715256,-1.45928555538655) q[5];
cx q[5],q[9];
u1(3.97064140161541) q[9];
u3(-4.44399232332981,0.0,0.0) q[5];
cx q[9],q[5];
u3(-0.988922199045124,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.26449446828522,3.11427477949795,-2.85041164749054) q[9];
u3(2.58965036702050,3.63744971985010,-0.803029545550496) q[5];
u3(1.56006954810967,0.832735713148610,2.20219570219750) q[6];
u3(0.925574244198539,-1.63950587209835,-0.828671227813450) q[5];
cx q[5],q[6];
u1(1.88632749582693) q[6];
u3(-3.17598288787090,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.51528287930153,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.35630398453404,-3.92695099349631,1.72839280330073) q[6];
u3(0.510841576071156,2.42845062501983,-0.761699240681035) q[5];
u3(1.14977870329233,3.00164488059955,-2.36534852297502) q[9];
u3(1.64427254919880,1.09150994021984,-1.95651945895425) q[0];
cx q[0],q[9];
u1(-0.134015421975215) q[9];
u3(-2.59319333752593,0.0,0.0) q[0];
cx q[9],q[0];
u3(1.52439315619662,0.0,0.0) q[0];
cx q[0],q[9];
u3(1.31874258852036,2.13844114773520,-1.07217311432071) q[9];
u3(2.15976008882819,0.358210486871416,-3.33211447196678) q[0];
u3(0.418558168781645,1.80381834948547,-1.68936359913915) q[8];
u3(0.958601955451186,-3.23362627665325,1.27973893708086) q[2];
cx q[2],q[8];
u1(1.93446852934044) q[8];
u3(-3.20385214885665,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.47631163054757,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.63552346226168,2.78414162048231,-2.38486547513466) q[8];
u3(1.99029594201635,4.04387190937504,-0.600087707734672) q[2];
u3(2.14156609468306,0.889178509650340,-0.403366495422227) q[3];
u3(0.0923085065681890,-4.16343431004392,-0.533483354034777) q[4];
cx q[4],q[3];
u1(3.74536614983155) q[3];
u3(-4.27411208631271,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.686417028986096,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.85137050687784,1.89194251486544,0.283131334367086) q[3];
u3(1.50524914886159,-0.589229803319996,1.52418143188222) q[4];
u3(1.01891355776238,2.18864784247051,-0.497835766778997) q[1];
u3(1.24926279068895,0.540806424159240,-4.08280225987030) q[7];
cx q[7],q[1];
u1(2.82548716762345) q[1];
u3(-2.26738563850606,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.0989493585682830,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.33663081852495,-0.543906052218618,3.18293592398148) q[1];
u3(2.04302754119757,2.15927910131658,-2.50668633636958) q[7];
u3(1.68428413018096,-0.928219983380991,-2.17478569895877) q[1];
u3(1.42275426363371,-4.44573550984172,1.53944262777779) q[6];
cx q[6],q[1];
u1(0.687649735122432) q[1];
u3(-1.39511009525755,0.0,0.0) q[6];
cx q[1],q[6];
u3(3.10352542547224,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.73251730857833,-4.29661407296596,0.262122541976209) q[1];
u3(2.62530779180082,3.22496285433987,0.169137938996121) q[6];
u3(1.81376732549218,2.82684394153442,-0.627204829840184) q[9];
u3(1.38351436048712,1.22484965392828,-1.17643229385727) q[8];
cx q[8],q[9];
u1(1.52823406350814) q[9];
u3(-3.12689796801708,0.0,0.0) q[8];
cx q[9],q[8];
u3(2.77149581053910,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.68230780860552,2.68553801191282,-0.990657694033055) q[9];
u3(1.22316309667595,-1.70997889694678,-4.38920625227099) q[8];
u3(1.43209095401450,-0.278255372763918,-1.75590091605100) q[7];
u3(1.24654998411192,-4.26560699114035,1.68821913598477) q[0];
cx q[0],q[7];
u1(0.706733777625031) q[7];
u3(-3.39459258296892,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.88832526370782,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.23828225048571,0.0893078359551649,-0.454844801038483) q[7];
u3(1.90064076997700,0.562947253449511,2.55928293302803) q[0];
u3(2.76340727939776,2.62067404931279,-2.67670000809939) q[5];
u3(1.72552401639869,-2.53895139497755,3.16071983070461) q[4];
cx q[4],q[5];
u1(2.91786807209178) q[5];
u3(-2.14054126463386,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.72262866284642,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.04698061349678,-1.86193725819686,3.83611331654136) q[5];
u3(2.70237472199716,-0.537005834344931,-1.93338099771390) q[4];
u3(1.68775565897124,0.193742382649553,1.88524201785476) q[3];
u3(1.32540188967237,-2.81100469181780,-2.28768075806550) q[2];
cx q[2],q[3];
u1(-0.0284159658166103) q[3];
u3(-0.594297520694083,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.45837989082301,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.44291614818444,-1.76783195527772,0.394549020105712) q[3];
u3(1.57212475841998,5.73968360892214,-0.434522452428442) q[2];
u3(0.662760524819155,-0.843910621279413,0.132938372990555) q[1];
u3(0.824405853371783,-3.16820962635186,0.728479622296197) q[9];
cx q[9],q[1];
u1(0.153452415493752) q[1];
u3(-0.743826161336711,0.0,0.0) q[9];
cx q[1],q[9];
u3(1.29790243038346,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.09164167089062,3.49031538034649,-2.33260594880375) q[1];
u3(2.97796434989779,3.44155321672320,0.174055823514794) q[9];
u3(0.413799421964422,2.82634118799478,-2.45032169317764) q[2];
u3(0.263873272235841,1.00613754422146,-1.97621533460593) q[5];
cx q[5],q[2];
u1(0.689577271252583) q[2];
u3(-0.313167627300217,0.0,0.0) q[5];
cx q[2],q[5];
u3(2.09723081071794,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.28996313593294,-0.334435438729307,1.11670277040025) q[2];
u3(1.84653431748413,-1.27339009411474,3.91099827677543) q[5];
u3(1.58228326601714,-2.69377507402037,-0.0202899070918419) q[0];
u3(2.29094258987688,-3.56274571291963,0.0153264626239702) q[3];
cx q[3],q[0];
u1(1.77745058896596) q[0];
u3(-3.23009950477440,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.24211607902042,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.71217306031986,-2.15113457546350,2.17635784299623) q[0];
u3(1.56702086906192,3.99402876219173,-1.36172130487076) q[3];
u3(1.87621465272392,0.0334282594861195,2.08844488252384) q[7];
u3(1.43486324091762,-0.769599922518032,-2.36026146968196) q[8];
cx q[8],q[7];
u1(1.96378936101133) q[7];
u3(-0.0370348780154670,0.0,0.0) q[8];
cx q[7],q[8];
u3(1.14142983070217,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.41026811238384,1.09133550865822,0.405037967383483) q[7];
u3(2.55774655378484,-0.0121251564892175,0.559932182584145) q[8];
u3(0.243854788504962,-1.45418553802679,-1.55785086182683) q[6];
u3(1.73218152963693,-3.15154974762824,0.250130883217305) q[4];
cx q[4],q[6];
u1(0.135731385455743) q[6];
u3(-1.44977387535613,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.57285235237310,0.0,0.0) q[4];
cx q[4],q[6];
u3(1.49922843114737,-1.01829653769939,5.04187807312932) q[6];
u3(2.20176943448090,-1.20449740256955,-2.98981596212319) q[4];
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
