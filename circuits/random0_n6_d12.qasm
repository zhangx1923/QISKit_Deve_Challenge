OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(2.93284789922020,1.55467663932720,1.33691713081566) q[2];
u3(1.28076665538053,-2.62317518479312,-1.69030308260787) q[4];
cx q[4],q[2];
u1(0.0574136901303806) q[2];
u3(-2.38005901254708,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.09735081023558,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.539572105910565,2.25035914600327,-3.04225965954633) q[2];
u3(1.79020094791040,-1.89143363322726,-2.40116754046976) q[4];
u3(1.73164403939838,-2.38295334854234,0.683749279682151) q[3];
u3(1.51498236244666,-2.49607440242390,0.483169667506401) q[1];
cx q[1],q[3];
u1(2.00019028278140) q[3];
u3(-2.84985603030164,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.775215207047808,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.50339667904831,-2.53641259905597,0.482765645863918) q[3];
u3(1.87153123256646,1.27707257292476,-0.984317938733402) q[1];
u3(1.81178462706527,-2.05870630750539,-0.780452755251127) q[5];
u3(1.91240060271789,-3.81830657482913,-0.0354402628583235) q[0];
cx q[0],q[5];
u1(1.81501649074593) q[5];
u3(-3.13062906341209,0.0,0.0) q[0];
cx q[5],q[0];
u3(1.44941057406421,0.0,0.0) q[0];
cx q[0],q[5];
u3(2.06069866369286,2.50708136817892,-1.95614355948759) q[5];
u3(2.24722156174141,1.98752544645073,-4.05024281920728) q[0];
u3(0.902987606369737,-3.67564049108359,2.52706358083060) q[0];
u3(1.41473848400571,2.91801538453334,-2.76729398544883) q[4];
cx q[4],q[0];
u1(2.44985130667966) q[0];
u3(0.00814292075036893,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.30481700737079,0.0,0.0) q[4];
cx q[4],q[0];
u3(0.859436172915610,-1.84035482046038,3.70255729095705) q[0];
u3(2.24488423842366,-0.364637668906760,-5.18540174856914) q[4];
u3(2.95583862938549,1.31571321194042,1.55707412306897) q[3];
u3(1.36555889017957,0.331311671889358,-5.35115417067480) q[2];
cx q[2],q[3];
u1(2.45328901839187) q[3];
u3(-1.64989270457535,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.311193965349398,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.790036753282717,-2.82959755793328,2.13670458894115) q[3];
u3(2.07615716732289,-2.08317345983050,3.89985252396859) q[2];
u3(0.894857714110646,-0.540848292851028,2.39972024742515) q[5];
u3(1.74792749562656,-2.23171616813580,-1.22260791612345) q[1];
cx q[1],q[5];
u1(-0.572021775814161) q[5];
u3(0.692280261238880,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.65088380006295,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.71712585836479,-3.81934720047707,1.09757465883945) q[5];
u3(1.56245805031663,-1.85070426736270,-4.31665719907699) q[1];
u3(1.06170326523038,1.75990100218177,-1.93027683586719) q[0];
u3(0.887506263530055,1.43402242296269,-2.80339042636088) q[1];
cx q[1],q[0];
u1(-0.117462779986945) q[0];
u3(-1.05458346940262,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.29795912448895,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.14203109598314,-1.74204579396608,3.40907456816957) q[0];
u3(2.35182569356695,-0.325467593040177,2.82542029080118) q[1];
u3(1.62261187769460,0.0546586866621807,1.01835352249545) q[3];
u3(1.53305326906884,-2.34846816685472,-1.56547544151916) q[2];
cx q[2],q[3];
u1(1.31637707861517) q[3];
u3(-2.64505445846071,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.129071389093308,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.787773279813053,-0.360378733309183,-0.157479737378784) q[3];
u3(1.75777968281923,0.779645098672157,-0.297193658307259) q[2];
u3(0.462958191973042,-2.69798717326515,3.52651872680405) q[4];
u3(1.36354644796558,1.38825772227312,-2.20184662947188) q[5];
cx q[5],q[4];
u1(3.47247573297321) q[4];
u3(-1.21162493763963,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.69200256242953,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.894319399088769,-0.474741193727963,1.14330284375790) q[4];
u3(2.37626455157013,0.856097399714030,4.27812963943640) q[5];
u3(2.80052997473470,-0.377542779587075,2.13216267021868) q[5];
u3(1.50587677996792,-2.22937930702090,-0.655013792717810) q[2];
cx q[2],q[5];
u1(0.0984462771612300) q[5];
u3(-1.20870916770900,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.43441082918467,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.77558542020082,-1.43490830682448,4.73921145200075) q[5];
u3(2.05622355647458,3.78765859521264,1.15175878408115) q[2];
u3(1.33365988698247,-2.01165351610919,-0.224148421296253) q[1];
u3(1.43646478157810,-4.10923811634385,-0.182709740389419) q[0];
cx q[0],q[1];
u1(1.45450971709103) q[1];
u3(-2.99540721410737,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.32051740946608,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.51443680151741,2.04084436690511,-2.71774699571911) q[1];
u3(1.03499376176615,1.53819351206937,-2.19023017693280) q[0];
u3(2.37858732522158,-1.77363095241631,1.72571218759376) q[3];
u3(2.51903555904667,1.01127631770682,3.23624316648064) q[4];
cx q[4],q[3];
u1(2.20099595963124) q[3];
u3(-3.26224046922598,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.21390340720193,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.04293197764397,-2.34556354770546,1.86223863609751) q[3];
u3(2.35103673728638,0.751581733085309,1.76513428931427) q[4];
u3(1.82236032387473,-1.82086115701232,-0.789881075596190) q[0];
u3(0.824189509041164,-4.43182322858669,-0.653876049721115) q[5];
cx q[5],q[0];
u1(1.03365965109734) q[0];
u3(-0.262865270293622,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.67688567264943,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.03540459591128,2.02094186441340,0.0337469576077758) q[0];
u3(3.04133004486218,-1.99656602925147,-0.458181952149559) q[5];
u3(1.33117297444283,1.23374019984898,-2.07224292543970) q[3];
u3(1.67968710095610,-2.08872591770345,2.94883534140709) q[1];
cx q[1],q[3];
u1(3.15769157017495) q[3];
u3(-1.41986505484785,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.31447121031175,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.810965576683143,-1.65014801342831,3.04312321407940) q[3];
u3(1.99567794354183,0.200566942846087,-2.82221438362802) q[1];
u3(1.39974775029665,0.0440058076254822,0.949711821952408) q[2];
u3(1.43642488919866,-1.19564719341807,-1.91593312571008) q[4];
cx q[4],q[2];
u1(3.11437518092299) q[2];
u3(-0.885812947451540,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.13230636079053,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.816852175556279,-2.07411857041824,0.963457849331902) q[2];
u3(1.71467791560808,-0.585549709013957,-4.53840399153655) q[4];
u3(1.43332414450429,-3.98993137688276,0.889364893529043) q[3];
u3(2.49694238495534,0.617002414981143,4.10767803721698) q[4];
cx q[4],q[3];
u1(0.0599747175908587) q[3];
u3(-0.662320031366957,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.76686262641089,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.392097764672795,2.51189002236511,-3.24771471039988) q[3];
u3(2.46926958257836,0.387840789426523,2.95650785545570) q[4];
u3(0.328305322633267,1.56737074331050,-1.43664623322218) q[1];
u3(0.555774163224211,-3.20868732750548,1.85773667565143) q[2];
cx q[2],q[1];
u1(2.30948310323470) q[1];
u3(0.0778466241472715,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.53882289374335,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.31904619279211,-1.52351781795672,2.91368044324117) q[1];
u3(2.21666121593965,-3.18194206568992,-1.74490994582249) q[2];
u3(0.805464494001188,0.441648764035022,2.04681315963142) q[0];
u3(0.991561974063672,-0.504273440166231,-3.05747978251375) q[5];
cx q[5],q[0];
u1(2.30597683159460) q[0];
u3(-1.77596298980261,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.175554160811533,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.37245180295733,-3.65700568743464,1.35449128779092) q[0];
u3(1.18774244457489,-4.55899024435412,1.56753228033784) q[5];
u3(2.44029160538307,-0.789944876921806,1.73587165101477) q[2];
u3(2.39052447077350,-3.08725449962947,-2.56949708162351) q[4];
cx q[4],q[2];
u1(0.887621898998349) q[2];
u3(-1.39978621318227,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.354019383035230,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.48854117878003,-3.23153716988892,0.733940238903506) q[2];
u3(2.13629857712521,-0.829049494765646,-2.31410500065882) q[4];
u3(1.75963223341583,-2.82745678709656,0.617270197689859) q[5];
u3(2.78332802747575,-3.73460016082916,-1.49000455871568) q[1];
cx q[1],q[5];
u1(-0.00732522448607997) q[5];
u3(-1.70897408602777,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.672095587884215,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.87471748699099,-1.82831084998753,1.79993932631905) q[5];
u3(2.12060270699430,-1.55263945134728,3.53623713589743) q[1];
u3(2.96485434050922,3.09794363465746,-2.84495167901077) q[3];
u3(1.29881245291485,-0.466465913588707,1.89481573285707) q[0];
cx q[0],q[3];
u1(1.61854741826475) q[3];
u3(-2.96807863610899,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.510902415199066,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.09847419271858,-2.32799127959670,-1.32702671270523) q[3];
u3(1.49524526880956,3.83119163313303,1.78772033350103) q[0];
u3(2.28111454663061,2.52175871730846,-0.305701632903266) q[3];
u3(2.49493539414595,1.11285042850748,-4.57849403803650) q[1];
cx q[1],q[3];
u1(1.93224616221392) q[3];
u3(-2.30318525721252,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.195700403881262,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.08109350321517,1.32748886089504,-3.95741960778363) q[3];
u3(1.56276134310061,-0.0783032921967808,1.13104982480848) q[1];
u3(2.17476898256416,1.40902703419799,-2.40061644793419) q[0];
u3(1.15390349787005,-1.86952647752737,2.58471866473851) q[4];
cx q[4],q[0];
u1(2.81889718484361) q[0];
u3(-1.45085687162033,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.406095800882443,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.94007821775777,0.270377731818673,0.973301105244227) q[0];
u3(0.857420422707993,2.73750485525098,2.93479164981121) q[4];
u3(1.22804010942742,0.0994628335283831,1.58780791390534) q[2];
u3(2.28176153127029,-0.730487387703051,-1.92065952149441) q[5];
cx q[5],q[2];
u1(1.93506791409168) q[2];
u3(-3.00047127320977,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.29146448825429,0.0,0.0) q[5];
cx q[5],q[2];
u3(2.34445681604375,-4.59535414579831,1.06962764746223) q[2];
u3(1.96748481110736,-2.52428324890890,-0.899926811925029) q[5];
u3(1.27478381806909,2.58012398442777,-1.33463263551132) q[2];
u3(0.906890444301778,1.25425937855750,-0.894998521208127) q[5];
cx q[5],q[2];
u1(-0.751085681680823) q[2];
u3(0.426173631368821,0.0,0.0) q[5];
cx q[2],q[5];
u3(4.33039395182827,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.01232777831552,2.67834170803491,-3.42901926493100) q[2];
u3(3.00499199794679,2.18366790364805,-0.678432968392796) q[5];
u3(1.32115778996224,2.07479417452030,-3.60197797398879) q[4];
u3(1.91045624844393,2.72419871883132,-3.05921077579371) q[1];
cx q[1],q[4];
u1(1.74240632075766) q[4];
u3(-3.40844429406881,0.0,0.0) q[1];
cx q[4],q[1];
u3(2.29271096257334,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.08292948316623,-0.595428844611068,-1.95156376007273) q[4];
u3(1.62739218070428,0.618237387967246,-4.73551322804905) q[1];
u3(2.30600343533304,-1.69535112819782,4.56220054255993) q[3];
u3(0.506168728288336,2.16669452148073,-0.719273455333316) q[0];
cx q[0],q[3];
u1(1.44129058240191) q[3];
u3(-0.143561982392010,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.44048938032727,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.00781043991611,2.65706495617005,-0.0918756525751576) q[3];
u3(1.93562580282332,2.43251176594289,-2.30511033882327) q[0];
u3(1.69878978277460,-1.05704998481178,1.88591924015891) q[5];
u3(1.86867544817063,-2.20600477877140,-2.00694402195769) q[0];
cx q[0],q[5];
u1(0.317393471156736) q[5];
u3(-1.37185379543022,0.0,0.0) q[0];
cx q[5],q[0];
u3(3.11599029038517,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.705004987152978,-1.57156756149042,1.41825065690837) q[5];
u3(2.77677671081094,0.602622742549629,-0.258498316756519) q[0];
u3(0.803538599072256,0.0845347172510168,1.43988822174732) q[1];
u3(1.35644755309838,-0.331411614628705,-2.45611043439369) q[2];
cx q[2],q[1];
u1(0.829682065906819) q[1];
u3(-1.15782460776437,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.71362863676645,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.545634637321690,-1.97913519720721,1.61943015577459) q[1];
u3(0.987856348072462,-0.110437930472718,1.69101222814926) q[2];
u3(1.88223502662673,-2.01906111586017,-0.106138955257684) q[4];
u3(2.60873913750636,-3.13071775946406,-1.46076637246523) q[3];
cx q[3],q[4];
u1(1.95759616862377) q[4];
u3(0.357847775444051,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.828991372758454,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.08853601918843,0.741188466804126,1.14276379568020) q[4];
u3(0.674820358024266,-2.06877881808422,2.82760069148007) q[3];
u3(1.72888057999695,-0.754391914545306,2.83061544061079) q[2];
u3(1.41870079101265,-1.36873040328536,-1.77053618614040) q[3];
cx q[3],q[2];
u1(3.09123848601429) q[2];
u3(-0.676070850580520,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.02686275544961,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.19299773542257,-1.78880833000679,1.24699256432423) q[2];
u3(1.64830289927423,-1.12965924309638,-1.42682362764457) q[3];
u3(1.73717102000521,2.32087388892058,-3.44210063317208) q[0];
u3(0.503520101330572,-0.597520267534049,1.48390120260128) q[4];
cx q[4],q[0];
u1(3.31208556148044) q[0];
u3(-1.16499654442255,0.0,0.0) q[4];
cx q[0],q[4];
u3(2.45718388815230,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.65506566697891,-1.68011925870866,3.19184921359751) q[0];
u3(2.75064810623025,-1.31772219101990,-0.0884865304222084) q[4];
u3(0.625605438963294,2.11961814604430,-1.62365358150420) q[5];
u3(0.731571686054474,0.682411017067116,-3.33666029835219) q[1];
cx q[1],q[5];
u1(1.52626195637255) q[5];
u3(0.0839599401759721,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.906336908318755,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.394089695211839,-0.931828455496081,1.71472720184082) q[5];
u3(0.439839573479203,0.849392616242796,2.29597326477546) q[1];
u3(1.14692540672044,1.54498623443145,-3.77679302499957) q[3];
u3(0.587891778426076,2.69696888773843,-2.39066181707206) q[1];
cx q[1],q[3];
u1(1.69367615212473) q[3];
u3(-2.29697438589555,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.97547484870073,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.70883550581260,2.29497662363739,1.34781511141279) q[3];
u3(2.98914277426148,2.83079070190098,-2.46042584868912) q[1];
u3(0.548290089513538,1.30033508783392,-1.94302405590726) q[0];
u3(0.849679129069152,2.11420743356253,-3.56281063176271) q[5];
cx q[5],q[0];
u1(1.76726156242298) q[0];
u3(-3.81200458427676,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.27616372816870,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.77796671397840,1.65974071923590,-1.66349378486598) q[0];
u3(0.827054325053024,-0.874744779400598,-2.26022602117634) q[5];
u3(2.28219942104327,3.57918462735864,-1.76045218424649) q[4];
u3(1.28682555825796,0.754973071632544,-0.919104709069738) q[2];
cx q[2],q[4];
u1(3.23959668694487) q[4];
u3(-0.873821206566268,0.0,0.0) q[2];
cx q[4],q[2];
u3(2.02094682596376,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.66756218368485,0.404534804683572,2.62024826702209) q[4];
u3(1.31214163326559,-1.16764056657165,-4.58072684330020) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];