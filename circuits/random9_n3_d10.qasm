OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
u3(0.563274756663646,-3.29684576842153,2.60159282670212) q[2];
u3(1.17409367692332,-2.90009234113118,2.70843932492500) q[1];
cx q[1],q[2];
u1(0.548133152379027) q[2];
u3(-1.38482678397028,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.17020042032997,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.61233416495651,-1.46315997729936,-0.0782889627743186) q[2];
u3(1.46333492360740,2.47646252125564,2.69538794154071) q[1];
u3(1.09085770836735,3.07153589096136,-1.64805812434302) q[1];
u3(0.664431478790669,2.04560322115843,-2.93810420511957) q[0];
cx q[0],q[1];
u1(0.816965864072915) q[1];
u3(-3.44672180193224,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.35356534473208,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.571134104853188,-0.151925840826247,3.16237867231022) q[1];
u3(1.45810172047080,2.54640306570901,0.592323751267296) q[0];
u3(1.74808828818397,-1.99021144930986,1.05054578294968) q[1];
u3(2.09585484694686,-3.84138366928067,0.289222067149761) q[2];
cx q[2],q[1];
u1(2.06375815411575) q[1];
u3(0.786853671260523,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.61273977117496,0.0,0.0) q[2];
cx q[2],q[1];
u3(2.50199254604067,4.24873287025688,-0.494955904344685) q[1];
u3(1.18406151137858,0.0673588392074933,4.83245287719758) q[2];
u3(2.13192558376433,-0.0602388870518190,-0.327157107411538) q[2];
u3(0.828303587588124,-0.419081538880816,-5.15390953665396) q[0];
cx q[0],q[2];
u1(1.86579093691847) q[2];
u3(-2.48838719144779,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.157141375430701,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.06815848398267,-2.69260713006815,0.484165103415670) q[2];
u3(0.734016000763162,-0.142032041307177,2.54329859807588) q[0];
u3(2.02469975513580,2.20055878536357,-3.37358733164005) q[0];
u3(2.20326887336243,2.41265176662373,-3.50249203668964) q[2];
cx q[2],q[0];
u1(0.764168012405195) q[0];
u3(-0.0995333535041147,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.08860332092229,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.05561788988817,-0.165125742496966,2.94935600765146) q[0];
u3(2.22002623776649,0.824297532285130,-3.43846216064247) q[2];
u3(1.28498600011361,0.482854794261833,0.373206037497241) q[2];
u3(1.31107478333172,-1.25953045681649,-1.88325766596362) q[0];
cx q[0],q[2];
u1(1.23373798110443) q[2];
u3(-1.53509098940411,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.28971134780130,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.48431569567809,-3.15623983004759,1.44072520441993) q[2];
u3(2.51988728754106,2.44802014653485,1.25854792524082) q[0];
u3(2.61677436920301,-2.86788209081619,2.77619647928024) q[1];
u3(1.77615758017659,-0.903795883471160,2.50021169081453) q[0];
cx q[0],q[1];
u1(1.96214357233219) q[1];
u3(0.188111736062498,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.04122751955856,0.0,0.0) q[0];
cx q[0],q[1];
u3(2.40185015792988,-1.36857804510991,-0.762906355197290) q[1];
u3(1.93840425006856,-2.13246081670015,-2.37950756023778) q[0];
u3(1.38572816490843,-2.07011728915373,0.734277005276130) q[0];
u3(1.56313972959906,-3.67354475741118,-0.118291205583024) q[2];
cx q[2],q[0];
u1(0.903551256816666) q[0];
u3(-3.25633681322312,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.96362755451100,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.43390894788021,-0.440956141354924,3.18816220199639) q[0];
u3(1.08801047674827,-0.118911258344782,-5.52621388169586) q[2];
u3(1.52445979994744,-1.35465938248367,0.420932730707481) q[0];
u3(1.63457834456717,-1.70933174456051,-0.605830751937972) q[2];
cx q[2],q[0];
u1(1.99431143709993) q[0];
u3(0.0904508804169610,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.48978902999033,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.127527264808186,-0.918690813707721,3.61804501417856) q[0];
u3(1.12234004361385,-4.09770837189848,2.06425275497515) q[2];
u3(1.73348167099163,-1.00927173221122,-1.02414746617204) q[0];
u3(0.502300013118336,-4.34960768701310,1.01109665134536) q[2];
cx q[2],q[0];
u1(2.46897419488552) q[0];
u3(-2.88536906846760,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.02592133711850,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.51308387298282,0.248185275556815,-0.981133459469749) q[0];
u3(1.39323790347281,-2.82139356969036,-2.72590468498617) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
