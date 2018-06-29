OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.03271761986032,-0.388685515225535,0.324637454513847) q[2];
u3(0.860129590327875,-2.30281445453229,-0.691339630148175) q[1];
cx q[1],q[2];
u1(3.08847194477879) q[2];
u3(-4.21432456465104,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.139984694797425,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.50653527832718,0.729927713668813,-4.65863038725427) q[2];
u3(2.52328587692002,2.67958499891868,-3.51973785896757) q[1];
u3(2.82452654221988,0.242624287816202,0.329888990240871) q[3];
u3(1.21803358236194,0.219392995822401,-5.28486163487126) q[4];
cx q[4],q[3];
u1(1.11370539940017) q[3];
u3(-3.51168427291259,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.48991536113432,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.32166869953155,-0.408439245435059,0.100168334110868) q[3];
u3(1.93445965390313,-1.74338595285617,-2.39305201953928) q[4];
u3(0.593018321884281,0.551849077263205,-0.0813227737652349) q[0];
u3(0.852861235132505,-0.722869782212280,-1.38090748994976) q[1];
cx q[1],q[0];
u1(2.01011792982342) q[0];
u3(-2.32284180911500,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.511856375905680,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.65695364551087,0.218378561244656,3.61797956805451) q[0];
u3(2.67470680291016,2.88524229616506,-1.76699332198517) q[1];
u3(0.255510755538725,-1.98267293586840,2.32494431497155) q[2];
u3(1.43143662385902,1.02886018883098,-1.61781535378015) q[3];
cx q[3],q[2];
u1(2.39388153233272) q[2];
u3(-1.79078091184996,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.93994991767763,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.26921521524974,-0.231254730840165,-1.45531789566896) q[2];
u3(1.69695855280124,2.57281727856543,-1.11254454950048) q[3];
u3(0.843378430043220,2.20142436271308,-1.72672418963314) q[2];
u3(0.644429744062567,1.08279456473143,-1.46872932034942) q[1];
cx q[1],q[2];
u1(2.84802758840697) q[2];
u3(-2.26089865984516,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.40712177781404,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.856972648424866,1.05627113742093,1.08123968393012) q[2];
u3(0.231867802449651,0.406954092897255,-1.46578588671400) q[1];
u3(1.40081604959522,-0.256241097631829,-1.37440033086878) q[4];
u3(1.77217331601339,0.271565700921713,-4.45490634328169) q[0];
cx q[0],q[4];
u1(2.67698567425751) q[4];
u3(-1.27870768818365,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.523363120165306,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.40971189158607,-1.17101655777993,3.53073998315920) q[4];
u3(2.02162254048708,2.97349528097024,0.186990050336506) q[0];
u3(1.21463130274157,-1.13498152464544,-0.636241132329653) q[3];
u3(1.20217968110076,-2.97583089268967,0.702964519448337) q[2];
cx q[2],q[3];
u1(3.38319734393409) q[3];
u3(-0.565837282293133,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.59501091431958,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.98869409688121,-3.18900040587936,-0.268621633048468) q[3];
u3(0.726307894023586,1.50058361321896,-3.56393471775422) q[2];
u3(0.172433340436333,2.54355192955795,-3.22077721651079) q[1];
u3(1.53046250103812,-2.86290912528425,2.88955014009799) q[4];
cx q[4],q[1];
u1(1.74390174475306) q[1];
u3(-0.00666291692495569,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.25171180414947,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.14691531948069,-0.821488328315468,0.272698792453177) q[1];
u3(1.96512743393253,-0.277340037337977,-1.82940725114894) q[4];
u3(2.67103813568182,-1.78287282299666,0.379713100945493) q[4];
u3(2.39773452491427,-2.13288855697194,0.0948216979934426) q[0];
cx q[0],q[4];
u1(2.69716862523328) q[4];
u3(-2.13580240702783,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.33773359830560,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.707648863638216,1.61496032438587,-0.644542547088145) q[4];
u3(2.25764456010535,-0.288257475980079,1.57947569550160) q[0];
u3(0.238789449398980,1.20338595188197,-0.857683608377393) q[2];
u3(0.183104189853499,-4.10141897556426,1.45510920899996) q[3];
cx q[3],q[2];
u1(3.00768435322566) q[2];
u3(-0.630610385688019,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.20382637166731,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.63976697934595,-0.221196178217682,3.08816673470420) q[2];
u3(1.41515183605734,1.80994507871549,2.03336975753137) q[3];
u3(2.04850659575599,1.84672856150207,-3.18644872006871) q[1];
u3(0.852684565147655,-1.78352935947050,2.69055243247764) q[2];
cx q[2],q[1];
u1(3.22910360934118) q[1];
u3(-0.966859436541708,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.53116790717255,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.42669330433706,-2.31593087561065,3.32615756357350) q[1];
u3(1.29255811992309,3.45569914733677,-0.603267789317222) q[2];
u3(2.34985882257803,1.35017966306096,-2.59647805900409) q[4];
u3(0.772921801822376,-2.66150277023882,2.15360769251972) q[3];
cx q[3],q[4];
u1(2.60442564768655) q[4];
u3(-1.80141457346866,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.248687917676440,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.51949231471890,2.31658722988397,0.321006013827970) q[4];
u3(1.11951198157406,-0.653928667783832,3.25470694856350) q[3];
u3(2.62346090925151,0.0602092762330543,0.186096411883683) q[4];
u3(0.709832517935885,-2.96458138120570,-1.33088069260109) q[2];
cx q[2],q[4];
u1(3.18598181489452) q[4];
u3(-2.02349695576066,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.59385937052614,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.877580795356633,-0.269803917144473,3.58269893122987) q[4];
u3(0.467134774854320,-4.96292564011318,0.0973175315050088) q[2];
u3(2.15462676516106,1.62421449912586,0.0388501048218893) q[3];
u3(2.53225674800924,0.543365312064487,-4.06170997327994) q[1];
cx q[1],q[3];
u1(2.47172642447817) q[3];
u3(-2.04833825231915,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.292279753478362,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.21881736768090,-2.84866878726296,2.73390637643406) q[3];
u3(1.13418207818212,3.64526823730883,2.25074667635640) q[1];
u3(2.16772343988535,0.568245038728524,-1.36385263922149) q[4];
u3(1.74203695908860,0.532706283971440,-4.48153554867479) q[0];
cx q[0],q[4];
u1(0.536626099338901) q[4];
u3(-0.830836989873533,0.0,0.0) q[0];
cx q[4],q[0];
u3(4.28125681920532,0.0,0.0) q[0];
cx q[0],q[4];
u3(3.01021414928511,2.76098778645715,-1.60551147560738) q[4];
u3(2.08881184203945,-3.85014957159075,2.21365897212075) q[0];
u3(0.642424127767886,2.00326875517771,-4.20795229400522) q[1];
u3(1.74943546029904,2.82162521481643,-2.75982914748581) q[3];
cx q[3],q[1];
u1(-0.117925344352212) q[1];
u3(-0.978085714434351,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.74157052166532,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.84984933082754,2.44809246495670,-2.35591995054685) q[1];
u3(1.46490683921286,1.50265049195448,0.309831779812055) q[3];
u3(0.787868979834942,3.16652619712440,-2.76630464197613) q[1];
u3(0.415929857160230,1.37750039668766,-2.30701135330688) q[2];
cx q[2],q[1];
u1(-0.163116910683146) q[1];
u3(-1.79389965740532,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.00065681188403,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.87051612803175,0.244372912760733,1.19647467684579) q[1];
u3(1.88201032570761,-1.60444589760787,-0.0889303155631950) q[2];
u3(2.05403478964268,0.0353123291976005,-1.59863613381538) q[3];
u3(2.41086195361009,0.697454578025495,-4.55676285761950) q[4];
cx q[4],q[3];
u1(1.76833002083669) q[3];
u3(0.138102485372203,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.920405329155884,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.252021122633140,-2.75754104467956,2.02810455925159) q[3];
u3(1.66150635065219,0.612581630395447,-0.171898031589109) q[4];
u3(1.55422244321113,1.98089188302227,0.296614532014100) q[4];
u3(2.66341095750538,1.12862122845979,-1.46506079334075) q[2];
cx q[2],q[4];
u1(1.37769842722293) q[4];
u3(-0.635569751933340,0.0,0.0) q[2];
cx q[4],q[2];
u3(-0.242474614772421,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.88497247308991,1.60575756182048,-2.43230888404807) q[4];
u3(1.22259188940081,-2.46579040422873,-3.08213282509490) q[2];
u3(0.885598326570219,-1.57110251837467,0.262276397622473) q[3];
u3(0.774158658976638,-1.40572910184521,0.187088518520213) q[1];
cx q[1],q[3];
u1(-0.201539112110612) q[3];
u3(0.557243850040579,0.0,0.0) q[1];
cx q[3],q[1];
u3(4.04266090424476,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.48106031516915,-3.53925655017042,2.47340319098435) q[3];
u3(1.75001028999418,4.17963363167418,-1.49844203483087) q[1];
u3(0.859852813317000,-1.07814303936851,2.00683817484457) q[3];
u3(0.279506620135715,1.24858858840154,-2.65152978043206) q[1];
cx q[1],q[3];
u1(2.96331301997583) q[3];
u3(-1.77898910708680,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.93590774855513,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.55403491188182,1.52724233750665,-1.61563532075480) q[3];
u3(2.64143131439748,0.876120866347368,4.83942352063337) q[1];
u3(1.87074152825746,-2.16370472469278,-0.462328105271537) q[2];
u3(1.59663693367264,-2.53223409349446,1.08614310165835) q[0];
cx q[0],q[2];
u1(2.59819442584074) q[2];
u3(-2.95675715604170,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.822905561825826,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.75823523950227,-0.798546364790638,-1.77790469548826) q[2];
u3(0.376741395822595,-1.90180787908072,-0.0999859817807768) q[0];
u3(1.25009077193177,0.598814463491414,0.796937419327606) q[0];
u3(1.58437051468093,-0.981152125838668,-1.37078746614584) q[3];
cx q[3],q[0];
u1(1.43502513443373) q[0];
u3(-1.02728002246669,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.09727376337282,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.39198447738150,-3.22213551780989,1.84702028847042) q[0];
u3(2.13017110101281,-0.996135076714689,3.37926323781106) q[3];
u3(1.39629872433248,1.88076878435914,-1.48919334937889) q[2];
u3(0.365606098091520,1.31087421939133,-1.76591491030365) q[1];
cx q[1],q[2];
u1(0.781204940331888) q[2];
u3(-1.52518772977353,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.07697420349260,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.39252707016303,-3.31928411442396,2.27458065972470) q[2];
u3(1.81377767254926,-0.481797081492530,0.363721440890833) q[1];
u3(1.57051141932655,3.23053327005248,-1.21964328300828) q[0];
u3(2.19798082579836,2.33348395386755,-0.672825225281405) q[2];
cx q[2],q[0];
u1(-0.486668789283022) q[0];
u3(-2.18655709001049,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.68637411135091,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.35879058965364,0.0632077575311090,0.449289595039447) q[0];
u3(2.57768271618202,0.887965833829919,-0.402491163407657) q[2];
u3(2.97372772450583,1.37619589303760,-0.855431207029479) q[3];
u3(2.40656002032975,4.10646525444979,-0.835496432755893) q[1];
cx q[1],q[3];
u1(0.899725631615912) q[3];
u3(-3.55346956466105,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.91603030546999,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.70332532378490,-2.41416018498048,1.53211582642127) q[3];
u3(0.635604387975539,-0.351260471042461,2.35595930485186) q[1];
u3(1.72354686665348,0.532311597476793,1.10214756626135) q[3];
u3(1.27303245096653,-0.735654549298867,-1.94359683987733) q[2];
cx q[2],q[3];
u1(1.58895317812464) q[3];
u3(-2.49514035310377,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.466189105283081,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.40718996502630,-2.17307195730373,0.712256142344379) q[3];
u3(2.59229630443764,-1.58276820716721,-1.80762931803676) q[2];
u3(1.22296926687191,-2.04504474659616,-0.953398229319149) q[4];
u3(1.56134968366026,-3.51486945834301,0.0525421337683112) q[1];
cx q[1],q[4];
u1(2.12225261964659) q[4];
u3(0.207779075621661,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.11998403365943,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.39332115627030,-0.360423780638908,0.812167080386400) q[4];
u3(2.20781705323254,-2.32375316454743,-0.587702250677704) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
