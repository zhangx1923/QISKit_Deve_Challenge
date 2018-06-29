OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.80936005348599,-0.336252953630898,0.763405639918589) q[7];
u3(1.55689969715668,-2.38570256130666,-0.965491197548602) q[3];
cx q[3],q[7];
u1(1.38627631330121) q[7];
u3(-0.116970647668569,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.57814253478221,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.63485765856116,-0.892520597145178,1.86597754231265) q[7];
u3(1.56241653383608,-3.97952317669204,-0.228010828003463) q[3];
u3(1.11080282823896,1.73536431325269,-2.17054983465641) q[1];
u3(0.805854428081098,1.35551657496736,-3.05437881035920) q[6];
cx q[6],q[1];
u1(1.10402044288854) q[1];
u3(-1.51397895924267,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.337437758119301,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.62745541584397,0.740551416364092,-2.81555466650946) q[1];
u3(1.80349669720851,1.56211310454194,-3.21982385883494) q[6];
u3(1.50463786530599,1.36306866523440,-0.0693957181244954) q[4];
u3(2.30747823296547,-0.874181905965889,-2.36004184117145) q[5];
cx q[5],q[4];
u1(-0.269843238252224) q[4];
u3(-2.10938850014503,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.16365846274557,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.67014500586677,-1.95214787375308,1.47701602130501) q[4];
u3(1.38801753641459,3.66593880889593,-1.90659368283847) q[5];
u3(0.696084059791931,-1.48205272407421,0.376396344664372) q[0];
u3(1.43536194290106,-1.75898458726970,-0.423282511490754) q[2];
cx q[2],q[0];
u1(1.01584866627592) q[0];
u3(-1.50099498995212,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.57633930026724,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.85545290854250,-2.88412523551563,2.30366639773529) q[0];
u3(1.37494339801785,2.73920563076860,1.46941179244969) q[2];
u3(0.986792530956029,-2.85356267704999,2.01606390941545) q[2];
u3(0.839397029283636,0.597633345713859,-2.15803106237429) q[4];
cx q[4],q[2];
u1(3.16477962886146) q[2];
u3(-1.20544649801924,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.422681641261997,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.37063069227940,-2.33062561138948,0.440784453031870) q[2];
u3(0.734783266495234,3.90877441002989,-1.55934123328483) q[4];
u3(0.697737927711428,-1.77543205174974,1.44944974814936) q[6];
u3(0.464362976582438,0.975567586245852,-2.01015143889095) q[5];
cx q[5],q[6];
u1(1.38221398709596) q[6];
u3(-0.304430751145528,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.81980460748481,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.424872844472006,-2.27931255143861,3.47150161851409) q[6];
u3(2.04442742206460,1.01105326033006,-3.55227856911730) q[5];
u3(2.60122887448291,-1.58081393008237,1.93149708859281) q[1];
u3(2.08441736137624,1.20659152993968,2.51077795020457) q[7];
cx q[7],q[1];
u1(1.71014191961806) q[1];
u3(0.418688521115473,0.0,0.0) q[7];
cx q[1],q[7];
u3(1.08454338701081,0.0,0.0) q[7];
cx q[7],q[1];
u3(1.36508636350073,-0.468522604628933,1.20608933187606) q[1];
u3(0.797923836357232,2.07185528469061,-1.93737242973468) q[7];
u3(2.52227477069571,1.77777984543509,-1.92594037829372) q[3];
u3(1.53661774493479,2.55845339249273,-3.23876858493328) q[0];
cx q[0],q[3];
u1(1.71510540616815) q[3];
u3(-2.44106072258900,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.0383866272934175,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.671454321119596,-3.05068667933138,-0.584452890995435) q[3];
u3(0.863826506289188,-0.736197708571036,-5.17033855122427) q[0];
u3(2.29247539838744,1.13308645396214,-3.44746886387863) q[0];
u3(1.29897804010614,-2.46427154541869,2.47951991975592) q[7];
cx q[7],q[0];
u1(3.16738600849833) q[0];
u3(-2.33491455030932,0.0,0.0) q[7];
cx q[0],q[7];
u3(0.975136949462294,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.57586830296690,0.723853471258579,-1.30496400222383) q[0];
u3(2.49290658270884,-2.18622600237231,3.53956014184827) q[7];
u3(2.74222962072932,4.22991492371394,-1.95382390085333) q[4];
u3(1.35005705488285,1.49485650948847,0.780030436264805) q[3];
cx q[3],q[4];
u1(0.699083366841851) q[4];
u3(-1.29964935528305,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.13113996566682,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.79980082838424,0.507521566569562,1.10245597962619) q[4];
u3(1.81295914197708,0.234363906216788,4.82098634230812) q[3];
u3(2.19496252991929,1.20769677231193,-0.0400968027421016) q[2];
u3(2.04074842227847,-0.631711402197551,-2.25903794358773) q[6];
cx q[6],q[2];
u1(1.42676036431394) q[2];
u3(-0.349097698647683,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.99886624648666,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.01656354173567,4.10186442725093,-2.17575551684132) q[2];
u3(0.410836016538977,-2.02018194984468,-1.24527399873053) q[6];
u3(2.09609958545087,-2.87273644051461,0.680794338854813) q[5];
u3(1.26080809344786,-3.34757291470968,0.339736818434770) q[1];
cx q[1],q[5];
u1(-0.0113088417665141) q[5];
u3(-1.00832680329582,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.91168049444047,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.74102668400628,-0.956454338647397,2.27492978564468) q[5];
u3(0.989433601343982,-1.61241761257725,1.53724755713420) q[1];
u3(1.76357413994897,0.230452673687155,1.07094799871806) q[2];
u3(1.66128329550661,-1.15083713563059,-1.59516666385527) q[7];
cx q[7],q[2];
u1(2.19677601554843) q[2];
u3(-2.84184377464175,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.67156858746094,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.831679566306212,0.437576128324888,2.74786188461428) q[2];
u3(1.27966298327624,-1.42216171856593,-1.85966290834365) q[7];
u3(0.366785467234374,-3.31428637072358,2.93246459695877) q[6];
u3(0.547909964488790,-3.43780618279158,1.16586371665051) q[4];
cx q[4],q[6];
u1(1.67092063629028) q[6];
u3(-0.117954775278287,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.20402718899038,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.21330710452058,1.55472378762267,-2.97870892344365) q[6];
u3(1.21870720370798,3.98461811096440,0.287472195999824) q[4];
u3(2.82327049002639,-2.65651609087788,1.00285583539409) q[0];
u3(2.28687935562610,1.16907149619099,3.30124654277510) q[3];
cx q[3],q[0];
u1(-0.745280567207157) q[0];
u3(-1.96858616150944,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.75511077024770,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.28202827451511,1.18776626491149,2.14413513558805) q[0];
u3(2.81462808298789,0.387414173377455,2.19095124599410) q[3];
u3(1.71736981464830,3.75142372326567,-0.737542280778213) q[5];
u3(2.03552122248892,2.30651381105730,-1.59289002383308) q[1];
cx q[1],q[5];
u1(2.00672982620796) q[5];
u3(-2.35684759338705,0.0,0.0) q[1];
cx q[5],q[1];
u3(3.27950734532134,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.82583974314518,-2.28578159541726,0.433305908544960) q[5];
u3(1.76015236303930,-3.71441853026385,2.53864622167308) q[1];
u3(1.45731265393144,2.50548594860235,-1.90456714960659) q[1];
u3(0.0904352637033468,-3.15129173610233,2.36874402152803) q[3];
cx q[3],q[1];
u1(2.95185888114744) q[1];
u3(-2.67600832350331,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.07069186914076,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.562948143949983,-0.838936387179969,3.12397266925590) q[1];
u3(2.59537292373294,-2.80448960317363,3.26208726117080) q[3];
u3(2.81908006515410,0.0263055118291711,-1.71166131056567) q[2];
u3(1.44715089800489,1.54579592270272,-3.97904726766383) q[5];
cx q[5],q[2];
u1(1.08384072270715) q[2];
u3(-0.629552147765543,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.70552080388803,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.69261213518937,3.77546536018010,-2.45080659153805) q[2];
u3(1.76152208844663,1.71480865440807,0.614173987407534) q[5];
u3(0.629585166024216,2.34221297397340,-1.73306846849099) q[4];
u3(1.51916522779072,1.51824537723120,-0.877175056980443) q[6];
cx q[6],q[4];
u1(4.04477242819387) q[4];
u3(-3.67470418966032,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.412996032774364,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.758268374136378,2.17575747995560,-2.22181177341689) q[4];
u3(1.94339962161099,0.0635807349606330,2.12991278065325) q[6];
u3(1.63279600893728,1.88636685160800,-3.42967055158282) q[0];
u3(0.686338725545947,2.05727789415253,-2.09805573171209) q[7];
cx q[7],q[0];
u1(3.43244779994458) q[0];
u3(-0.637889492685143,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.80304193448172,0.0,0.0) q[7];
cx q[7],q[0];
u3(2.42205039608898,3.99032994159003,-1.24326821289307) q[0];
u3(2.38101822141490,-1.90389739288800,3.66875875047278) q[7];
u3(1.90247625137697,-1.97840374934420,-0.776996203036632) q[7];
u3(2.15081781639783,-2.24473190315071,0.223879659674743) q[3];
cx q[3],q[7];
u1(1.69788310065054) q[7];
u3(0.0350248761396035,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.01988263208415,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.49369498749342,0.723760440392157,2.78429547295582) q[7];
u3(1.00316328851626,4.98790881646514,-0.903288757285628) q[3];
u3(0.995401294135919,2.17059586631589,0.948886065944478) q[6];
u3(0.648208964293690,-0.115681364453144,-2.92874982040376) q[2];
cx q[2],q[6];
u1(-0.464939441707799) q[6];
u3(-1.70943823523082,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.10984160360724,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.39152649300906,-3.96311757649977,0.467212991482683) q[6];
u3(1.56557052569823,3.47380119019787,-1.38890081883325) q[2];
u3(2.63438038433560,-4.24879302946939,1.86710788860911) q[5];
u3(0.333805229532501,2.77393602012798,-0.894410316219859) q[4];
cx q[4],q[5];
u1(-0.662122334192271) q[5];
u3(-1.77418951571662,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.46483931957887,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.66090843753377,2.06282928298093,-2.02256784948330) q[5];
u3(2.63122049439136,1.74480408168040,-2.10284002175511) q[4];
u3(1.00413103320380,0.166729197299844,1.21153730341286) q[0];
u3(1.41567037613831,-0.221990034586675,-2.21817251342835) q[1];
cx q[1],q[0];
u1(0.647368666235841) q[0];
u3(-1.42603461958047,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.23154118311814,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.49119508423086,3.05376062952459,-1.66718182494202) q[0];
u3(0.656138933702605,-0.595004227783564,-2.60698151905688) q[1];
u3(0.703352239626993,-2.35021218859410,2.10438093710210) q[7];
u3(0.664511747142650,0.763060157837277,-2.28046766516888) q[6];
cx q[6],q[7];
u1(2.41223559834252) q[7];
u3(-1.42669171417456,0.0,0.0) q[6];
cx q[7],q[6];
u3(3.58422813123847,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.77833936716613,1.21497047186007,-3.14841054138745) q[7];
u3(0.464351480557858,4.08494981425941,1.54806203102061) q[6];
u3(1.49299938952115,-1.14630986033008,0.264031375404562) q[1];
u3(1.56450083628994,-2.93694649422080,0.162657531235801) q[4];
cx q[4],q[1];
u1(2.88455795811000) q[1];
u3(-1.77543347662444,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.846932235024801,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.39107485042555,2.75536927343117,-0.706929933782043) q[1];
u3(1.02983093420920,3.09275829352397,-1.68039088745634) q[4];
u3(1.57428322449882,3.80393819278321,-0.700315656547428) q[3];
u3(1.50495055483952,3.16705834385315,0.390631164159579) q[5];
cx q[5],q[3];
u1(0.973242622940583) q[3];
u3(-1.60152222738765,0.0,0.0) q[5];
cx q[3],q[5];
u3(-0.356242777163240,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.62911589935391,-0.421945084317366,-0.364278014764977) q[3];
u3(1.83611527878757,0.203131245614371,3.53256639726801) q[5];
u3(0.994868231155019,2.06687883905194,-0.906708805991631) q[0];
u3(1.25666433762618,1.55651304347590,-0.609570137547652) q[2];
cx q[2],q[0];
u1(1.59801598408567) q[0];
u3(-3.02645453363206,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.12037601791516,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.87493335177356,1.16848729042948,-4.29521287337825) q[0];
u3(2.33334737336562,1.32987160724553,-3.73551843756592) q[2];
u3(0.596012156284960,0.00197065102230021,1.17281718899976) q[7];
u3(0.818245616146165,-1.62178427730516,-0.262586680445612) q[5];
cx q[5],q[7];
u1(1.60825521198382) q[7];
u3(-2.34299138574397,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.11089894446343,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.568161011382638,-0.838324842396755,4.29287716507976) q[7];
u3(1.57876860312352,-0.0665324595522183,1.38585959296577) q[5];
u3(1.83706391007095,-0.515190113466795,0.988820266952558) q[6];
u3(1.11169145748351,-2.14055103172009,-1.74348183467740) q[3];
cx q[3],q[6];
u1(0.333841757748445) q[6];
u3(-1.13890836477980,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.21022407044481,0.0,0.0) q[3];
cx q[3],q[6];
u3(0.447018441841844,-0.832128095527063,-1.63624116915195) q[6];
u3(1.65476727059704,5.09347640070591,-0.800016576758561) q[3];
u3(1.93347652487325,-0.579173493331554,1.37700488574479) q[1];
u3(2.18085593116759,-1.61682174699769,-0.684931318836344) q[2];
cx q[2],q[1];
u1(-0.111554466948548) q[1];
u3(-1.50176711668175,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.46060661944592,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.06493948569912,-3.34984616439786,-0.0691356714229514) q[1];
u3(0.497313150338996,0.349357217984047,-3.21104342094685) q[2];
u3(2.83626181394572,1.49711259075917,-2.23777165371993) q[4];
u3(1.89248216283481,1.56468138545298,-3.03023904383674) q[0];
cx q[0],q[4];
u1(-0.192741736119418) q[4];
u3(-1.72830312881666,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.902172882015558,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.02482095500429,1.43398433193910,0.531182543213717) q[4];
u3(0.918280379401999,-0.0265667720711025,5.47941145394615) q[0];
u3(1.89935915083243,-0.401498219656834,2.21616017407467) q[2];
u3(2.48154996600638,-2.10112402736101,-0.725699034460357) q[4];
cx q[4],q[2];
u1(1.92793834895252) q[2];
u3(-2.10809919968380,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.833932975851274,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.30909719645205,2.55055449860119,-2.91257688146415) q[2];
u3(1.96130643990166,1.48149266584843,4.12147076425051) q[4];
u3(0.579458221695613,2.33914902945871,-2.47459481122244) q[6];
u3(0.641192296321870,0.602022008439262,-1.07233604684554) q[7];
cx q[7],q[6];
u1(3.18199370984522) q[6];
u3(-1.44410191311052,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.29831386568450,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.12251861003589,2.89205605063393,0.694835410338198) q[6];
u3(1.95904660935989,-0.123936548012971,1.00641093531500) q[7];
u3(1.25894505324267,-1.78083211879878,0.144390368972910) q[1];
u3(1.00796022103761,-4.28102034996468,-0.404289060113824) q[5];
cx q[5],q[1];
u1(3.02204230745664) q[1];
u3(-2.13178508415417,0.0,0.0) q[5];
cx q[1],q[5];
u3(1.49942956186934,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.70564326747279,2.57760232023088,-0.663366467075783) q[1];
u3(2.02539234068869,4.08269127305828,0.181926320914655) q[5];
u3(3.03153404585272,2.36075921409806,-1.60168191062640) q[3];
u3(1.63755416800687,1.02420118138563,-1.79549749193523) q[0];
cx q[0],q[3];
u1(2.71047819217040) q[3];
u3(-2.81922811059663,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.16441939491118,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.849932602101333,-2.80450449097362,1.36756021480950) q[3];
u3(2.05084872556464,-3.04839010016928,-0.937778636045957) q[0];
u3(1.30160924600539,2.17800468874110,-0.168492170172122) q[6];
u3(2.25467355755791,0.661914049962393,-3.02564317768906) q[7];
cx q[7],q[6];
u1(0.659300827502531) q[6];
u3(-0.301710880876503,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.37756524104485,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.84519660176107,-2.89158030419491,-1.14850945803973) q[6];
u3(1.96429587532202,-4.21134814209858,-1.16784932216110) q[7];
u3(0.502034716451674,0.459238929772681,-2.11265143847491) q[0];
u3(1.49630729844902,-3.46321691870478,1.97317098623079) q[1];
cx q[1],q[0];
u1(3.11620062209918) q[0];
u3(-1.85475261662720,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.558041334538817,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.83875380621973,3.15766331856686,0.111469767100630) q[0];
u3(2.64919640840531,0.433974320486927,1.30510087199692) q[1];
u3(2.04366556642291,-3.03342929834573,0.707707845083974) q[5];
u3(2.46819702516595,-2.98616492068147,-2.18992328008606) q[2];
cx q[2],q[5];
u1(2.01563289176890) q[5];
u3(-2.87737193366551,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.688898966064039,0.0,0.0) q[2];
cx q[2],q[5];
u3(2.22338405763684,-3.55985363386150,2.56447601261339) q[5];
u3(1.21323659131214,-1.20782593728160,3.20568453888290) q[2];
u3(1.75201579953742,-2.40302920049618,0.505289072586612) q[3];
u3(1.84828106785914,-3.37394769138064,0.807408335674689) q[4];
cx q[4],q[3];
u1(2.11882870286823) q[3];
u3(0.197628358114327,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.55239635238983,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.37546010458079,-3.80489363618282,1.17967053030466) q[3];
u3(1.58348030632159,-0.657309564776074,5.17617656165293) q[4];
u3(1.97191424695613,-1.00831121379136,3.32354290361874) q[7];
u3(0.886057009106311,1.67454385470741,2.12090979089859) q[6];
cx q[6],q[7];
u1(3.47761809199746) q[7];
u3(-0.821057025830129,0.0,0.0) q[6];
cx q[7],q[6];
u3(1.92005248091472,0.0,0.0) q[6];
cx q[6],q[7];
u3(0.759737771926181,-0.154549354540149,-1.12548236560804) q[7];
u3(1.02560068588142,1.92422222992585,1.25361625113236) q[6];
u3(1.27359985572353,-3.31277220922031,2.90108325039909) q[4];
u3(0.715232634667721,-0.403689019116072,-1.63882864828986) q[3];
cx q[3],q[4];
u1(3.69469264072109) q[4];
u3(-3.17276316590821,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.720770451677959,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.83906850763069,3.67890066555806,-2.54817371583161) q[4];
u3(0.302392124446663,-0.497104052298044,-3.11934852912320) q[3];
u3(2.24561086131229,0.245831199186760,-1.84464087471636) q[2];
u3(2.55436219367708,-0.0941130112906547,-4.88138129296593) q[5];
cx q[5],q[2];
u1(2.31493361787179) q[2];
u3(-0.0187175648449984,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.04698011780816,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.469941423986979,-0.140958011570995,-3.47136098835848) q[2];
u3(1.90041310430941,4.58020873520867,0.387322327402769) q[5];
u3(1.51450177779536,1.97992797337187,-0.000163651824496069) q[1];
u3(2.22453511276139,1.03445808909118,-2.29531366336660) q[0];
cx q[0],q[1];
u1(-0.479212141184898) q[1];
u3(1.17876798884942,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.89403523399535,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.77289110526451,1.15111545645726,-2.29352766406820) q[1];
u3(2.63784100362205,-2.65796967756842,-2.19318028506826) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];