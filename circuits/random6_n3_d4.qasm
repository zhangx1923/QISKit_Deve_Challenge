OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(0.847166948368851,-3.37914595097797,2.08153500261116) q[0];
u3(1.04918260060753,0.128162103273667,-2.37323532099297) q[2];
cx q[2],q[0];
u1(1.13062183095735) q[0];
u3(0.238644231623963,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.08919588322040,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.879455195098819,-3.27916812906407,0.531729260739112) q[0];
u3(2.18231045621672,-2.07957959633492,-1.30113766328370) q[2];
u3(1.96663946480919,0.478358830340301,1.95875122945960) q[1];
u3(1.88869453883826,-2.04179087034419,-2.29305744774054) q[0];
cx q[0],q[1];
u1(2.06028024757490) q[1];
u3(-1.48309889422183,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.639837270238888,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.25784106439118,1.76203527854633,-2.63432090646038) q[1];
u3(0.724720189531658,-3.44660990212728,-1.63147177777679) q[0];
u3(2.18422973030011,-0.546614026383131,2.56933937801060) q[2];
u3(2.58953864569951,-0.166876185624836,2.28742373038119) q[0];
cx q[0],q[2];
u1(2.96546865767203) q[2];
u3(-2.20734458724113,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.26791110979865,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.41944759968484,-2.62012164594821,0.387952244544762) q[2];
u3(1.46882388608328,0.0616641333380059,2.97800315975907) q[0];
u3(1.47607577712486,1.84103209660813,-2.24976752485147) q[1];
u3(0.322111814749448,3.04219354765278,-3.14827983627058) q[0];
cx q[0],q[1];
u1(0.893039745248836) q[1];
u3(-1.55736874029990,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.430097157184121,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.25507043438305,0.222920930654651,-2.39671828292322) q[1];
u3(1.10055962848082,0.692463457924502,0.582869261035468) q[0];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];