OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(2.56530277691536,1.80293551223712,-3.36943847613654) q[2];
u3(1.33153412307874,2.17668866554628,-1.01130778845096) q[4];
cx q[4],q[2];
u1(2.89834175148333) q[2];
u3(-2.17329043367969,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.42990749778068,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.54540050837468,-1.35081097047048,0.285213950035046) q[2];
u3(0.883389834361324,0.230388939194364,-3.67728720041429) q[4];
u3(1.89317095753771,-0.827141944135645,0.0609006573018465) q[8];
u3(1.74482011330570,-2.49704917902274,-0.521243653825418) q[7];
cx q[7],q[8];
u1(0.0646013434917816) q[8];
u3(-0.900011300493891,0.0,0.0) q[7];
cx q[8],q[7];
u3(1.99206365614915,0.0,0.0) q[7];
cx q[7],q[8];
u3(2.36750413419749,-2.80641352230255,-0.834817335593478) q[8];
u3(2.49994342778670,-0.959835091532069,4.26194298634485) q[7];
u3(1.59741545798928,-0.768026949706647,-0.420867794298920) q[9];
u3(1.93007980916880,-3.29692028199409,0.578681258346240) q[1];
cx q[1],q[9];
u1(1.72519858879325) q[9];
u3(-3.03135285123562,0.0,0.0) q[1];
cx q[9],q[1];
u3(1.07178693813716,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.88227786763228,0.453063636594180,-2.12640847947099) q[9];
u3(1.37731431033625,2.17837666377080,-1.10703580226230) q[1];
u3(1.84646097533961,-1.20606854544742,0.529505236320970) q[3];
u3(1.70743544092417,-2.55206305750642,0.229138012104596) q[5];
cx q[5],q[3];
u1(3.08131110229886) q[3];
u3(-2.69730980166987,0.0,0.0) q[5];
cx q[3],q[5];
u3(0.902531641781961,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.64528947463668,0.993233001857753,-1.71770511483370) q[3];
u3(3.01270948050337,-1.87093300402809,4.34762684189358) q[5];
u3(2.16776134808037,0.274898488738712,-3.36082440149718) q[0];
u3(1.43714188868325,-2.91011602580902,3.15501855806698) q[6];
cx q[6],q[0];
u1(1.56073062045948) q[0];
u3(-3.34821841102502,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.49778637075402,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.35414949173058,0.475801357144002,1.63303666304874) q[0];
u3(0.934690006284744,0.321215942819063,2.90084468391652) q[6];
u3(1.38686519995937,2.42039022775889,-2.24014426575604) q[0];
u3(1.18736662531517,-3.60006620009489,2.44621576225129) q[3];
cx q[3],q[0];
u1(-0.162298183201621) q[0];
u3(-0.958540920187807,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.51052516952818,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.28784559743975,0.0775691965989460,-2.32397247071395) q[0];
u3(2.09133225081060,1.47264812693162,-4.02664980195249) q[3];
u3(1.72349751790430,2.85631459108894,-1.59281407406372) q[5];
u3(2.37334949397749,1.09838315873458,-1.57584735240109) q[6];
cx q[6],q[5];
u1(1.27580393835212) q[5];
u3(-3.62655959067474,0.0,0.0) q[6];
cx q[5],q[6];
u3(2.11294866396166,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.87689576883902,-1.00806848266090,3.48645000623971) q[5];
u3(2.14237774009457,-3.05295113003147,1.73333414082683) q[6];
u3(2.42676763436100,0.0815420206088688,-0.154717860056119) q[8];
u3(0.961092648458632,-3.61377316607286,-1.64684649239726) q[2];
cx q[2],q[8];
u1(2.81083257681466) q[8];
u3(-1.69932704548948,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.07894357008444,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.83629474632613,-1.30066403448529,2.81366758456250) q[8];
u3(1.66808341336258,0.822290714926868,0.306965817398312) q[2];
u3(0.490707033400566,-2.06076817664415,2.98198152831297) q[1];
u3(0.115055065469059,-2.75770996738705,1.60317412670622) q[9];
cx q[9],q[1];
u1(1.31608140123015) q[1];
u3(-0.236826953848215,0.0,0.0) q[9];
cx q[1],q[9];
u3(2.22066220682106,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.31994481491218,0.829813023463122,2.43070841143224) q[1];
u3(2.61086019261829,-1.77399617001747,3.98576424465658) q[9];
u3(1.65261059045919,0.0739965383226250,1.77692964743290) q[4];
u3(1.66105423979711,-1.54886464809027,-1.71491210192751) q[7];
cx q[7],q[4];
u1(1.47394800312481) q[4];
u3(-1.99670773224351,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.98439059634776,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.47052582319650,-1.79232747259955,3.75683667117879) q[4];
u3(2.12373154346668,1.11860810032056,2.81085613648064) q[7];
u3(2.63207039040488,0.622228494409320,-1.43832899200117) q[4];
u3(2.39946012660224,1.34455240782923,-3.25396069194982) q[0];
cx q[0],q[4];
u1(0.971723000032597) q[4];
u3(-0.159050326528059,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.83724118887867,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.54342236393194,2.00707853369530,-2.70708665828824) q[4];
u3(0.619975162457271,2.91995473961014,2.22310677904304) q[0];
u3(1.48837468292359,1.67196124201676,1.01912315095716) q[9];
u3(2.78168268384410,0.291503639359825,-2.51107948530628) q[7];
cx q[7],q[9];
u1(2.68992659963766) q[9];
u3(-1.87004806022745,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.0958744739200268,0.0,0.0) q[7];
cx q[7],q[9];
u3(2.03363421427640,2.00138677049093,-1.51411381771628) q[9];
u3(1.31533219206999,3.82355814017685,-1.12991736309037) q[7];
u3(1.89467321013229,-2.59654070619212,-0.106209427040553) q[8];
u3(2.79888394332822,-1.95291645587032,-0.348888548744606) q[1];
cx q[1],q[8];
u1(-1.08217202221398) q[8];
u3(0.901462152779287,0.0,0.0) q[1];
cx q[8],q[1];
u3(3.61085389730489,0.0,0.0) q[1];
cx q[1],q[8];
u3(0.539144143283461,-0.424054051369692,-1.25448297214725) q[8];
u3(1.32290053356485,-0.501529175499329,-3.37491301050013) q[1];
u3(1.36544715236796,2.47137904155121,0.278574995963728) q[6];
u3(1.94306497384564,0.602299231496458,-2.52162958860624) q[5];
cx q[5],q[6];
u1(1.53680975350713) q[6];
u3(-0.975854262827939,0.0,0.0) q[5];
cx q[6],q[5];
u3(-0.693664347182618,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.22726021346181,-1.34962636453899,4.54761773778632) q[6];
u3(2.72846490287909,-2.97973165461775,-2.90220915654991) q[5];
u3(1.44722542896164,1.85488845879843,-0.599047766930819) q[2];
u3(1.84346187950734,0.964018050161816,-3.27733235512418) q[3];
cx q[3],q[2];
u1(-1.17229598272428) q[2];
u3(0.653928137914771,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.92629760741108,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.02538641936410,-2.02757084069319,1.15632604934378) q[2];
u3(1.24573335982953,2.46018265606623,-0.962994952581277) q[3];
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
