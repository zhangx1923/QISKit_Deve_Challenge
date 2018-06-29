OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(1.68248099577100,2.12923605734030,0.437597377441467) q[0];
u3(2.59687458915872,0.0354926094936188,-2.86894375751280) q[1];
cx q[1],q[0];
u1(1.31290883781486) q[0];
u3(-0.599423889228592,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.80941796015716,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.97173326233102,-2.91060629502931,-0.449868600493826) q[0];
u3(1.63444732012519,1.09337049127174,1.58764018249206) q[1];
u3(2.27890279345304,-0.183672614780074,1.87893101854402) q[2];
u3(2.82702037541241,1.57596726217924,3.03323228822871) q[3];
cx q[3],q[2];
u1(2.37733543345152) q[2];
u3(-1.95323770218798,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.651732198099570,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.79875299587496,-0.160471138851035,4.21553214009897) q[2];
u3(1.24522739937619,-0.783413475611852,1.31383915556830) q[3];
u3(0.797289303850212,2.58221909656871,-2.92287689781812) q[2];
u3(0.904027101798975,2.57632702615864,-3.34421876611546) q[3];
cx q[3],q[2];
u1(2.68734048992965) q[2];
u3(-1.62994204539286,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.0847587547706903,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.611502994511496,1.64030998929426,-0.657371939865162) q[2];
u3(1.29214971304402,-2.64390610072668,-1.64351894989724) q[3];
u3(1.64606623808864,0.346656274185138,2.29112983218311) q[0];
u3(2.09830835870835,-2.11653220120081,-1.47962029656272) q[1];
cx q[1],q[0];
u1(2.30120249324990) q[0];
u3(-1.60566739092567,0.0,0.0) q[1];
cx q[0],q[1];
u3(3.14091693842051,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.93557703465232,2.96594843767391,-2.50333994824236) q[0];
u3(0.549160132764429,4.97387540180769,-0.988380698644068) q[1];
u3(1.47666490490658,1.66571042033290,-2.33898943739075) q[2];
u3(1.15583860839278,-2.26081670035368,2.46877289994256) q[1];
cx q[1],q[2];
u1(0.0754602016702335) q[2];
u3(-1.44235687497302,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.50336647953989,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.183576825794851,-0.979413103967131,2.17413406460323) q[2];
u3(1.18637595003884,-1.61724947734730,-2.63264939142360) q[1];
u3(2.00883830192212,1.67560074931595,-0.461812095925542) q[0];
u3(0.676746772364842,-0.618240808742083,-2.15415327671094) q[3];
cx q[3],q[0];
u1(3.75500637164034) q[0];
u3(-1.36541635539980,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.19426201866136,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.17875851941375,1.77286486539623,-3.76392822753466) q[0];
u3(2.04281373070890,0.561213713467009,-1.07259534059655) q[3];
u3(0.660668220946107,2.03504859161768,-0.424796921241582) q[2];
u3(1.58551221806045,-1.01873013508316,-3.38535968768775) q[0];
cx q[0],q[2];
u1(2.88017505195720) q[2];
u3(-1.99282404501963,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.45763727675238,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.00262598435361,1.62900027794006,-0.363027925936669) q[2];
u3(2.26117124267822,-4.06096702693679,-0.671083538361057) q[0];
u3(1.35236486901572,-0.847225599141883,-1.10945805640155) q[3];
u3(1.16371923559425,-2.69760056726428,0.352019209005171) q[1];
cx q[1],q[3];
u1(2.52799294539526) q[3];
u3(-2.98394650024324,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.50938661896689,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.81640980163483,-0.726677284679636,-1.88868725918069) q[3];
u3(1.94542633395178,2.11876193085216,-2.55247830377695) q[1];
u3(2.01068060156652,1.88829438368275,-0.495814765199972) q[3];
u3(2.06851054110009,-0.974544955531154,-5.17867286891963) q[0];
cx q[0],q[3];
u1(3.40466633953095) q[3];
u3(-4.34645721448272,0.0,0.0) q[0];
cx q[3],q[0];
u3(-0.0693661332803492,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.55878489684603,2.39249338710020,-3.77510464728992) q[3];
u3(2.04588855297998,3.36572229276917,1.43710572491411) q[0];
u3(1.65551705907818,-0.643291434937107,-0.682200263225761) q[2];
u3(2.70190765967794,2.23752804930173,-3.51270529361830) q[1];
cx q[1],q[2];
u1(2.95677228178102) q[2];
u3(-1.66270349696830,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.520542515276598,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.21551437383662,-3.53408591709582,1.27502561513932) q[2];
u3(1.22218994185643,1.46374754539414,-2.51970247121750) q[1];
u3(0.665580662643005,1.25625135784440,0.0613967812069779) q[3];
u3(1.41585893284933,0.439837158508032,-2.18131017461080) q[1];
cx q[1],q[3];
u1(2.18277037646768) q[3];
u3(0.341807020201476,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.62588719930591,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.76057031933577,3.97252766167349,0.0142296174749172) q[3];
u3(2.37448586510092,1.16306485824855,4.36873883541443) q[1];
u3(1.49201704145654,1.74729748157904,-2.72835232068053) q[2];
u3(1.30730010272472,-2.82368360761900,2.70408908668716) q[0];
cx q[0],q[2];
u1(-0.175716312593776) q[2];
u3(-1.68555138674820,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.48003992501139,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.80612490123951,1.25416573114320,-1.49632915226590) q[2];
u3(2.10489085336520,-5.27962841771330,0.487434836036007) q[0];
u3(1.51800552742043,0.564795862043914,0.898832652963622) q[0];
u3(1.21220294476405,-0.838708759067371,-2.57251966419702) q[2];
cx q[2],q[0];
u1(2.94730255287159) q[0];
u3(-1.61188375886814,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.840142340161241,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.59770329050309,-3.53632171806711,2.07195174451230) q[0];
u3(1.80723227392345,-1.39132902972470,-1.04457571649068) q[2];
u3(2.38938578024091,-0.846701162566758,0.261974590657316) q[1];
u3(2.19418895795748,-0.743451804364821,-1.48402169249521) q[3];
cx q[3],q[1];
u1(0.851788820362367) q[1];
u3(0.138470381951025,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.83841224559017,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.722289662030008,1.97550096936664,-1.33656581801805) q[1];
u3(0.121253659888336,2.12305664375364,0.703571150872705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
