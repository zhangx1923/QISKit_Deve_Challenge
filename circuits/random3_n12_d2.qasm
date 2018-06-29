OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.91750408370196,1.84979007464316,1.06102541367196) q[7];
u3(2.20999936819605,-0.338133481668511,-2.96419992164572) q[3];
cx q[3],q[7];
u1(1.31859920610074) q[7];
u3(-0.624337955050463,0.0,0.0) q[3];
cx q[7],q[3];
u3(-0.328234922273120,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.30969872521627,-1.15525841894166,-2.74702997124696) q[7];
u3(1.62501099115655,-4.46932615458401,0.126905620877144) q[3];
u3(2.05664330437211,-2.34369444215828,-0.776493667455210) q[8];
u3(1.19579245718534,-5.02528275238032,-0.788725408566116) q[9];
cx q[9],q[8];
u1(0.902323758359136) q[8];
u3(-3.25112473089640,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.34619078270831,0.0,0.0) q[9];
cx q[9],q[8];
u3(1.68899153003081,3.27435563430347,-2.12678184021395) q[8];
u3(1.22402424027765,1.59113897018294,0.0106253932448160) q[9];
u3(2.41867432268957,2.16750937349182,-2.02159300055516) q[2];
u3(2.11182393783377,1.24868746147065,-3.09850072911114) q[1];
cx q[1],q[2];
u1(0.576788732463042) q[2];
u3(-0.290270722630615,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.36830245417762,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.783513311801953,1.22383198255260,-4.06857760721212) q[2];
u3(1.41453170487905,-1.10717757850922,-4.95638436896323) q[1];
u3(2.42081944345656,0.115386266183847,-1.59441902985412) q[5];
u3(1.62667551189828,-3.72146556709935,1.18943978338152) q[4];
cx q[4],q[5];
u1(0.924401427172461) q[5];
u3(-0.875256919435714,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.505459103818746,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.86936066580068,2.13898075537127,1.84515128299490) q[5];
u3(1.43432564237614,4.18780537715322,-0.898387707637959) q[4];
u3(2.36031410104801,1.58125782241428,1.33573257570818) q[0];
u3(2.07551268116003,0.735593828875201,-2.70915174290471) q[6];
cx q[6],q[0];
u1(1.69465380840885) q[0];
u3(-2.19876982421320,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.269251661888658,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.929400132888838,2.21168198867583,-2.03250352311087) q[0];
u3(1.78849146917410,1.91865060439110,-0.0665197878384778) q[6];
u3(2.58276230315413,0.650893114396754,-2.94536677101269) q[10];
u3(1.69941423686125,2.88004541588005,-2.28830749399704) q[11];
cx q[11],q[10];
u1(2.26287548818205) q[10];
u3(-1.55708628058566,0.0,0.0) q[11];
cx q[10],q[11];
u3(0.672025093596350,0.0,0.0) q[11];
cx q[11],q[10];
u3(0.948848562808773,1.79398276448407,-0.973986380622910) q[10];
u3(2.05282962205411,1.08952172674139,-3.21030964486172) q[11];
u3(1.10872958114329,1.31998127306836,0.798430693889824) q[5];
u3(0.900510242961624,-1.30765498997550,-2.11458818377003) q[7];
cx q[7],q[5];
u1(1.46999867972614) q[5];
u3(-3.13517263583637,0.0,0.0) q[7];
cx q[5],q[7];
u3(2.19299485078000,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.557884901252637,-3.94350646816091,0.848215166151247) q[5];
u3(2.66273673005887,-2.02873701007300,1.71949298751963) q[7];
u3(1.02797857894723,-0.809732712135396,0.496143238946048) q[6];
u3(1.77264686390235,-3.19012422604944,-0.110911407340753) q[11];
cx q[11],q[6];
u1(0.454457375239869) q[6];
u3(-1.57162558515944,0.0,0.0) q[11];
cx q[6],q[11];
u3(2.49288349132597,0.0,0.0) q[11];
cx q[11],q[6];
u3(1.57346491628130,-1.82241580749147,0.912638431454331) q[6];
u3(2.76580621671596,-0.149006859085765,3.97100655525329) q[11];
u3(2.36310543816753,-2.85450489723428,2.68088087308892) q[2];
u3(0.521292578936748,0.690572367633615,1.12725060539129) q[1];
cx q[1],q[2];
u1(1.74587187024696) q[2];
u3(-2.44890323998015,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.981438298120273,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.608485042451005,0.258269401923666,2.33227242402228) q[2];
u3(0.377629462568020,1.59119502578516,-0.626038702233299) q[1];
u3(1.37596409739516,0.880904799818409,0.314545266286002) q[0];
u3(1.51539688534163,-0.651899963886622,-2.05405988213790) q[8];
cx q[8],q[0];
u1(2.40502836843370) q[0];
u3(0.166578446496913,0.0,0.0) q[8];
cx q[0],q[8];
u3(1.14983856336995,0.0,0.0) q[8];
cx q[8],q[0];
u3(0.202578417460982,0.221329846737757,2.00646571310555) q[0];
u3(2.73962560042158,0.160479752625814,3.44203809298959) q[8];
u3(1.70959922678403,-1.33364989349978,1.04048965145501) q[10];
u3(2.23987033226316,-1.48241056813348,-1.47048524789637) q[4];
cx q[4],q[10];
u1(0.293023028080928) q[10];
u3(-0.600152375593716,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.27220335956266,0.0,0.0) q[4];
cx q[4],q[10];
u3(0.785959315979418,3.59416950860345,-0.696197362159528) q[10];
u3(2.81444838716878,0.00295960457690247,-1.20403659158381) q[4];
u3(1.54951662395898,-2.43746137254683,0.503970117441745) q[9];
u3(2.23726534104206,-3.35884441200505,-1.15475774929847) q[3];
cx q[3],q[9];
u1(3.03926389879466) q[9];
u3(-1.81237570034604,0.0,0.0) q[3];
cx q[9],q[3];
u3(0.970383600391385,0.0,0.0) q[3];
cx q[3],q[9];
u3(0.579517120471161,2.08528164032145,-3.51402132301882) q[9];
u3(2.31919874471889,-0.362614360463940,3.13584296901682) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
measure q[11] -> c[11];