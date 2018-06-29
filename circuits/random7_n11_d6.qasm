OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.25195970044903,-0.187600135689948,0.904965474444069) q[4];
u3(1.30155163150340,-2.05654282661907,-1.55459686498563) q[5];
cx q[5],q[4];
u1(1.30649839082621) q[4];
u3(-0.388008798583250,0.0,0.0) q[5];
cx q[4],q[5];
u3(-0.00461974074516447,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.17883714153569,4.28567588788696,-0.741748765970754) q[4];
u3(0.941542049079355,4.06200947441714,1.82362754965732) q[5];
u3(0.557556214026886,1.60021080073984,0.0242325554255649) q[10];
u3(1.66489838164577,-0.126148030291763,-4.35538723174632) q[8];
cx q[8],q[10];
u1(2.26828359724849) q[10];
u3(0.290790529459019,0.0,0.0) q[8];
cx q[10],q[8];
u3(1.58480045824480,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.32481866055423,0.778173692694814,-1.44114737943372) q[10];
u3(1.85723110619627,-0.0567414496616991,-0.102254721770023) q[8];
u3(0.845622976315475,-0.330463386990643,1.58108169648041) q[1];
u3(1.04555398249430,-1.60147617934640,-2.05024422698396) q[9];
cx q[9],q[1];
u1(1.28212753071970) q[1];
u3(-0.948794898355808,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.96164537354467,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.45251800868656,1.72362065483116,0.628844839676670) q[1];
u3(0.702309334749738,0.527489395409755,-5.59911681087242) q[9];
u3(2.46141755451870,1.00282547553811,-0.584016106999983) q[6];
u3(1.75481014056061,0.0607440328198394,-4.14086991001022) q[2];
cx q[2],q[6];
u1(1.84312945158668) q[6];
u3(0.429772676449037,0.0,0.0) q[2];
cx q[6],q[2];
u3(0.822587204153600,0.0,0.0) q[2];
cx q[2],q[6];
u3(2.16470833989840,-2.02987805854316,1.21000816662382) q[6];
u3(1.26876991605362,-4.45440905192661,0.454006894771242) q[2];
u3(1.61421777550478,-1.15098053940644,-1.38372032754100) q[3];
u3(2.25525796189622,1.95891433901478,-4.08419177123954) q[7];
cx q[7],q[3];
u1(3.63772865322619) q[3];
u3(-3.20052843115041,0.0,0.0) q[7];
cx q[3],q[7];
u3(-1.06838013790920,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.50248760788322,-2.21063077482180,1.16417189957199) q[3];
u3(2.90917472091781,-1.86962448845464,-2.37910860180040) q[7];
u3(2.58102464234895,-1.76699358409514,-0.167094765606946) q[8];
u3(1.49823687233567,-5.26241138070241,-0.992232738419529) q[0];
cx q[0],q[8];
u1(1.48980474996512) q[8];
u3(-0.547063149041339,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.92589893614937,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.129448389723584,1.74057103066016,0.484447851177158) q[8];
u3(1.77076903795595,2.57005804909927,2.20299271671277) q[0];
u3(0.439557945153735,0.723054273513704,-1.12852185454675) q[9];
u3(0.873892872413269,-3.25928060285558,1.28758107314804) q[10];
cx q[10],q[9];
u1(-0.224929546225006) q[9];
u3(0.986304094037737,0.0,0.0) q[10];
cx q[9],q[10];
u3(3.72627844607374,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.20286306893182,-1.39148617310516,-2.64635474462383) q[9];
u3(1.63027469885404,1.99075458824658,-0.839795415211200) q[10];
u3(1.83468418095029,-1.14878070401550,0.426034055908518) q[3];
u3(1.69177142023470,-2.09849369288821,-0.267619880481450) q[4];
cx q[4],q[3];
u1(0.625155262023706) q[3];
u3(-1.08368758879679,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.10530078519094,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.749734360367884,-0.681416868004293,-0.0282646425253288) q[3];
u3(1.54020916522037,-3.36033573245477,2.79292039565060) q[4];
u3(1.23018756768533,-1.45253406229486,-0.162003526601225) q[6];
u3(1.15560642042499,-3.16185061160172,-0.891060750280528) q[7];
cx q[7],q[6];
u1(0.213895586954352) q[6];
u3(-0.920139538977157,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.40786271108898,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.313112119869183,-0.276080591439203,0.811832233835185) q[6];
u3(1.77113262796334,5.32229196599611,-0.660221520258974) q[7];
u3(1.69129823007815,2.08964135535827,-3.34999200791720) q[2];
u3(1.15033954326438,2.20197387101732,-2.41728433907062) q[1];
cx q[1],q[2];
u1(0.0548310270921353) q[2];
u3(0.323257341401005,0.0,0.0) q[1];
cx q[2],q[1];
u3(4.04505019592522,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.10717094642778,1.56299838729934,-0.235523410354022) q[2];
u3(1.29006092867249,2.68893527965797,-3.50200417162447) q[1];
u3(0.515747209956559,1.98693963065592,-3.42108971815958) q[8];
u3(1.32446201157385,3.32415214862397,-2.37750288161082) q[4];
cx q[4],q[8];
u1(2.30236760329526) q[8];
u3(-1.79511493127729,0.0,0.0) q[4];
cx q[8],q[4];
u3(3.37811121115343,0.0,0.0) q[4];
cx q[4],q[8];
u3(1.13303197616761,0.987953691667681,-1.12927083761028) q[8];
u3(1.93390453366023,3.21348548512847,2.12708837613157) q[4];
u3(0.737028123117788,2.66600102145343,-3.26061099760643) q[0];
u3(1.63656512671631,-2.81163620515762,2.78618011956167) q[2];
cx q[2],q[0];
u1(3.36461245776327) q[0];
u3(-1.48802691853787,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.01976255011451,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.60146617969352,1.39427767093522,0.493641014089377) q[0];
u3(1.95262526044803,3.90617859152977,0.523837869256782) q[2];
u3(0.486907297944672,-2.55538504787223,0.255745633641166) q[7];
u3(1.18895919951449,-2.56237068067447,-1.03046267972630) q[3];
cx q[3],q[7];
u1(-0.881198558439089) q[7];
u3(0.527464254390862,0.0,0.0) q[3];
cx q[7],q[3];
u3(3.38948982031312,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.69567517808116,-3.36452117540233,1.69142048468260) q[7];
u3(2.15365745755461,0.435148824301786,0.151424843760596) q[3];
u3(2.18048579202484,2.43922476788965,-1.27981754490940) q[9];
u3(2.36963416668411,0.681731070112824,-2.54302116820126) q[1];
cx q[1],q[9];
u1(0.161068261147976) q[9];
u3(-0.834819765594686,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.47262957678264,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.48561528070534,2.17202361631428,-2.02257075033621) q[9];
u3(1.66183533134590,2.94685424115537,2.51084329568149) q[1];
u3(0.734725578273977,2.26636364963359,-3.83967650395119) q[5];
u3(1.36550108453335,2.50153032247827,-3.04581082575160) q[6];
cx q[6],q[5];
u1(0.800709686344892) q[5];
u3(-0.236632451166980,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.07439614221953,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.957796922156838,-0.186067455775615,0.284101220457051) q[5];
u3(2.78231893587702,-1.64543335689626,-0.461389887587736) q[6];
u3(1.98034981775765,2.46980056649494,-2.58129111197694) q[2];
u3(1.85829466148836,1.44943854444635,-1.23376919080179) q[0];
cx q[0],q[2];
u1(1.39696521509038) q[2];
u3(-0.323972356917999,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.44937314402680,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.27147559367919,4.23729134837908,-1.33287087800625) q[2];
u3(0.940274921017214,-0.769844077586184,2.20423281566443) q[0];
u3(2.34885464468681,1.55775431827306,-0.978239988699414) q[9];
u3(1.93287139328347,1.20876856335276,-3.44877402720188) q[6];
cx q[6],q[9];
u1(2.04955021487593) q[9];
u3(-2.49321567289346,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.17282414793565,0.0,0.0) q[6];
cx q[6],q[9];
u3(0.669358968703525,-2.94069866455835,0.136903861063244) q[9];
u3(0.435175162970814,-0.868911685068451,0.171483058823106) q[6];
u3(1.57667236597347,0.932352454277281,-3.80689978721287) q[1];
u3(0.895827346904664,5.17624566804077,-1.08880950182920) q[8];
cx q[8],q[1];
u1(-1.37039048432966) q[1];
u3(0.342469538096679,0.0,0.0) q[8];
cx q[1],q[8];
u3(3.67018919249683,0.0,0.0) q[8];
cx q[8],q[1];
u3(1.66325508119966,-1.34786830393438,0.342179314005079) q[1];
u3(1.98880545838452,0.490819567286132,-3.32542588704184) q[8];
u3(0.303189727480528,-1.30089267929996,-1.02613549392266) q[3];
u3(1.25802768689600,-3.21969067530879,-0.348648509123016) q[4];
cx q[4],q[3];
u1(2.70461277038706) q[3];
u3(-1.89426129319050,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.902806410913077,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.23625550814593,-1.34916265554776,4.01121623188619) q[3];
u3(2.33184636356379,2.12746830336723,-2.05324479820374) q[4];
u3(1.89651652071730,0.761478115651213,1.02158333713481) q[10];
u3(0.794182986664062,-5.85953576339538,-0.0923232335699629) q[7];
cx q[7],q[10];
u1(2.17190807822238) q[10];
u3(-1.85792996573393,0.0,0.0) q[7];
cx q[10],q[7];
u3(0.652514432396783,0.0,0.0) q[7];
cx q[7],q[10];
u3(0.955008597632702,4.07096978231637,-0.0690176324005369) q[10];
u3(0.547056408806257,1.60306876317268,-1.47424783411730) q[7];
u3(0.872003680944969,-1.76895950422541,1.50913859863397) q[3];
u3(0.615613955626012,0.879313679287316,-2.80768514442746) q[5];
cx q[5],q[3];
u1(0.149296017332629) q[3];
u3(-0.673784179711486,0.0,0.0) q[5];
cx q[3],q[5];
u3(2.35748950668944,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.48601905219766,3.04497411860041,-3.03882606397568) q[3];
u3(1.20286701502815,0.00387481620154784,-4.14316842268583) q[5];
u3(2.45670365829912,2.10934581772322,-2.78155345685367) q[4];
u3(0.991839754503236,-3.06220090779846,3.12111675774366) q[7];
cx q[7],q[4];
u1(2.39702995529761) q[4];
u3(-1.85930810554393,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.0992290735177259,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.54979320295379,-0.927024020211850,-0.924476654617011) q[4];
u3(1.84483486480219,3.83382661900464,-1.32152215427292) q[7];
u3(2.42464539253235,2.49841582082232,-3.17118719598726) q[6];
u3(1.48476797623724,1.48688579548781,-1.62659437384852) q[8];
cx q[8],q[6];
u1(0.227320388165266) q[6];
u3(-1.21666692276462,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.38412722199274,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.775237373024900,-3.37163136493812,1.07317321285455) q[6];
u3(2.08190042637490,-3.83119942652059,-2.35471479493207) q[8];
u3(1.92842710377767,-0.739122876593114,1.04101687750526) q[10];
u3(2.03924734133577,-1.57716368023043,-1.44749634975173) q[1];
cx q[1],q[10];
u1(0.678928504087261) q[10];
u3(-0.212053809163297,0.0,0.0) q[1];
cx q[10],q[1];
u3(1.98837254928605,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.73259686810823,3.83319661325354,-2.34716991772517) q[10];
u3(1.27407263078550,-1.18912298475961,-1.93385053526289) q[1];
u3(2.48864202022856,0.193729981297576,-1.47010392583386) q[2];
u3(1.62741351043311,-3.26856068249426,1.12029247908178) q[0];
cx q[0],q[2];
u1(2.90947871967132) q[2];
u3(-2.65525755751691,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.46797873283704,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.66026770787468,3.48379004611840,-2.26437504532269) q[2];
u3(0.547354472290652,-0.471606808338062,0.734382018566921) q[0];
u3(1.22959018803858,-0.199950758563703,2.22088237786061) q[10];
u3(2.13750357387377,-1.55432898219161,-1.08378761996433) q[6];
cx q[6],q[10];
u1(2.93497831671543) q[10];
u3(-0.947417971342986,0.0,0.0) q[6];
cx q[10],q[6];
u3(2.04911399970177,0.0,0.0) q[6];
cx q[6],q[10];
u3(2.47241589635647,0.536308751598088,0.666231842660544) q[10];
u3(1.08550839054273,4.14755353995339,-1.27042288273310) q[6];
u3(2.92528097545090,2.43536514663778,-0.738479704267557) q[8];
u3(2.28080600838466,0.937591894780823,-4.17801684086176) q[1];
cx q[1],q[8];
u1(-0.565025494768093) q[8];
u3(-1.98401752522376,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.43734778888747,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.978816760154507,0.0953071957663656,2.43347767612035) q[8];
u3(0.702662032176928,-4.60602207553599,1.05452772373078) q[1];
u3(1.50454099157838,0.917444005281888,-1.34516468350693) q[0];
u3(2.02576728078597,-4.88166353141691,1.25949317957307) q[4];
cx q[4],q[0];
u1(2.75129870317739) q[0];
u3(-1.85301646845964,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.04113596711036,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.38488152508195,-1.06818050126837,-1.77006443732239) q[0];
u3(1.39594440761616,0.555944749333401,-4.72008609527017) q[4];
u3(2.08025668689936,-1.17666340102367,1.00466591345151) q[3];
u3(1.63160221949066,-1.77224231251079,-2.04107595235597) q[5];
cx q[5],q[3];
u1(0.969685145351178) q[3];
u3(-0.284522723801279,0.0,0.0) q[5];
cx q[3],q[5];
u3(1.79221953020001,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.18461370917823,-1.56678754048838,2.67436616407298) q[3];
u3(2.37435045955616,4.86218387592750,-1.37971292026900) q[5];
u3(1.90457075241736,0.749685200344293,-2.72612577858804) q[7];
u3(2.87610804278791,2.74011375712589,-3.08295556047658) q[2];
cx q[2],q[7];
u1(1.88063553612695) q[7];
u3(-1.72571529316087,0.0,0.0) q[2];
cx q[7],q[2];
u3(0.590946190741492,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.51608447060321,-2.78860134256660,-1.70169794172062) q[7];
u3(1.49887627959409,1.53397266778346,2.85482300971967) q[2];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10];
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
measure q[10] -> c[10];
