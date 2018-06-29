OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
u3(0.747102156801104,0.584613958085837,-1.22253387473309) q[1];
u3(0.668162381238627,-1.82352642627045,0.0306107752394493) q[0];
cx q[0],q[1];
u1(1.72372833352419) q[1];
u3(-2.45727539012275,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.247108890700103,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.05102307278178,-0.0166629351893108,3.20440183155790) q[1];
u3(2.78941117768354,-1.12448459853860,1.93350001113820) q[0];
u3(1.37608802620765,-0.0548685303649437,1.56677653965775) q[2];
u3(1.33577035241995,-1.58800088497693,-1.75344484300363) q[3];
cx q[3],q[2];
u1(2.77897481881972) q[2];
u3(-1.65423005204312,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.34927671063623,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.47271006568175,0.954969361767831,-1.51651887580168) q[2];
u3(0.270890881940546,0.291894033938540,1.96338731571535) q[3];
u3(1.56928409548886,1.98129572599946,-3.73318338933839) q[2];
u3(2.93118154111282,-1.18928931583275,4.28863106085688) q[1];
cx q[1],q[2];
u1(1.36349762659129) q[2];
u3(-0.0249119549981276,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.58250862448472,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.943596391851834,2.42123418244966,-2.04141535382533) q[2];
u3(1.33291689632174,-1.88688839296762,3.18193196003823) q[1];
u3(0.749974695008427,1.92938765766246,-2.85723658417536) q[0];
u3(1.41831944555048,-3.92760228190192,2.31185440811956) q[3];
cx q[3],q[0];
u1(1.66407453586774) q[0];
u3(-3.01195324724583,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.991958754263139,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.09611753499689,-3.69587170343462,2.57163181603648) q[0];
u3(0.574654035642404,-1.51397882686900,2.37706945717607) q[3];
u3(2.19152050939046,-0.262451775983512,-2.43227388827840) q[0];
u3(2.49092705566277,5.12077610470716,0.697994617420350) q[3];
cx q[3],q[0];
u1(-0.189767033370461) q[0];
u3(-2.27371223624365,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.47943236301316,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.64713862869013,1.93651771792615,-2.67153756249630) q[0];
u3(1.83065216273745,-0.0596108462066112,5.33084802190340) q[3];
u3(1.55288576278359,1.27104675871894,-2.88135216696451) q[2];
u3(0.620568088924636,2.82060890134054,-3.02968710986012) q[1];
cx q[1],q[2];
u1(3.07298902046188) q[2];
u3(-1.58697855258489,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.377788327333476,0.0,0.0) q[1];
cx q[1],q[2];
u3(0.719563395236532,0.710571478407837,1.76156423630110) q[2];
u3(1.90716650970008,-0.354925195934128,-1.03742266004664) q[1];
u3(2.33237098123303,1.36655954093039,-3.52561137274603) q[3];
u3(1.45327728745235,-2.32550571554486,3.40966065043248) q[1];
cx q[1],q[3];
u1(1.45940427443400) q[3];
u3(-0.140685344333897,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.586184638558290,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.67555223385703,0.216801914160587,3.17474435106763) q[3];
u3(1.85888986578102,1.11257485490255,-1.42930666504080) q[1];
u3(2.10467731039093,1.79175870284588,-2.68185442446766) q[0];
u3(1.04665100779845,-2.19242275549420,2.51038741799733) q[2];
cx q[2],q[0];
u1(-0.831076037991654) q[0];
u3(0.475738681191283,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.17324224357070,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.51970658872099,-0.120780819213442,-0.774635642867115) q[0];
u3(0.902242905901655,-0.367821937322209,-4.12570693104906) q[2];
u3(1.97621992570488,-0.240977664737330,2.47435861804503) q[3];
u3(2.40632181009553,-3.42262319880996,-1.59217447874361) q[2];
cx q[2],q[3];
u1(1.20005115584572) q[3];
u3(-3.52869705416212,0.0,0.0) q[2];
cx q[3],q[2];
u3(1.78736820150511,0.0,0.0) q[2];
cx q[2],q[3];
u3(0.174638702339727,1.03064927881540,0.619508755821978) q[3];
u3(0.558476584702839,-0.228845976097502,1.49008465831273) q[2];
u3(1.21186197693879,-1.48487636032012,1.60158538246500) q[1];
u3(0.327155606685056,-1.23991845064094,-1.74205929939572) q[0];
cx q[0],q[1];
u1(0.381151660345806) q[1];
u3(-1.56701522427552,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.185081185256143,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.50525302995234,1.39037971921224,2.34490442849804) q[1];
u3(1.64725323978953,0.769782250964087,3.89098085531940) q[0];
u3(2.02662481196251,2.78841424402947,-2.39092009122783) q[3];
u3(1.16065981957179,2.18889167811050,-2.58284254485020) q[2];
cx q[2],q[3];
u1(0.805914634970504) q[3];
u3(-1.35344475196143,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.0930582669157880,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.34476889143589,0.146230616552549,-3.24208263576287) q[3];
u3(1.23088825911275,-1.58597880724456,-0.526132286317870) q[2];
u3(1.96240094194758,-0.638157876328866,1.70439871400482) q[0];
u3(2.30636200273060,-1.17939835492647,-0.648418678090382) q[1];
cx q[1],q[0];
u1(3.14213948848729) q[0];
u3(-1.75594621829421,0.0,0.0) q[1];
cx q[0],q[1];
u3(-0.173340366192920,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.714481524259657,0.311955339012368,0.492301494012215) q[0];
u3(1.23708021887153,-1.88207629358731,4.25849305390397) q[1];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
