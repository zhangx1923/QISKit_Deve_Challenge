OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(1.35417129861041,0.244751042524190,1.61131750985704) q[4];
u3(1.59024914231534,-1.41713309831401,-2.65392462071897) q[1];
cx q[1],q[4];
u1(1.15205070328699) q[4];
u3(-3.23630579709420,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.93758406872071,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.10868304720744,-0.0109257584559797,4.31068628964772) q[4];
u3(1.92887463717576,0.368239730994008,-5.07733096459409) q[1];
u3(1.53935865285301,0.250681454377337,2.87678658311859) q[2];
u3(2.64253046729072,-3.03801237921710,-2.57254459877368) q[0];
cx q[0],q[2];
u1(1.05480052714601) q[2];
u3(-0.370220655635004,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.31137662810187,0.0,0.0) q[0];
cx q[0],q[2];
u3(0.884036038044735,2.37470906661718,-1.83860414272599) q[2];
u3(2.82577190766497,1.52179506866485,2.53885784227164) q[0];
u3(0.659880652075232,0.550621952681961,-1.22200293125061) q[3];
u3(0.339461039330643,-0.465651588382337,-0.673872831629843) q[4];
cx q[4],q[3];
u1(3.31630084699080) q[3];
u3(-2.62866853485711,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.10182440976491,0.0,0.0) q[4];
cx q[4],q[3];
u3(0.914391426701813,1.11044619489532,1.38892666210375) q[3];
u3(0.951369755026316,0.571513052782640,5.23752903222690) q[4];
u3(1.43912937889574,1.97447775027808,0.0251339788089988) q[1];
u3(0.837481259696175,0.788915608321124,-2.62487202693546) q[0];
cx q[0],q[1];
u1(1.42179027031884) q[1];
u3(-0.475967240425738,0.0,0.0) q[0];
cx q[1],q[0];
u3(2.14797123782507,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.61246546098514,-0.111595437365812,-1.66741474430976) q[1];
u3(1.08055103632712,0.674703501870742,-4.38079574150656) q[0];
u3(2.03490076734165,1.16197260800381,0.365784905428470) q[3];
u3(1.11238011178278,-0.452547058263443,-2.24838249730922) q[2];
cx q[2],q[3];
u1(1.84717536027516) q[3];
u3(0.200909548400796,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.754256275879184,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.32109363301120,-1.93493314064626,2.77853863494537) q[3];
u3(1.40949872519109,-0.705984706678991,5.18031236326628) q[2];
u3(1.31965144177677,0.279308187984834,2.17640985258004) q[0];
u3(1.26886568489548,-1.67048825160405,-1.72460695270039) q[1];
cx q[1],q[0];
u1(1.40425976021858) q[0];
u3(-2.99783695940898,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.702701581224422,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.00825041687115,2.12658519425322,-1.17568738897122) q[0];
u3(2.01317032371022,-5.38865777548922,-0.135764611776834) q[1];
u3(2.28348899564696,-1.41798871355189,-1.21352966720457) q[0];
u3(0.388794668377522,-2.67865129925377,-2.40191624530941) q[2];
cx q[2],q[0];
u1(-0.537765690382880) q[0];
u3(-1.77675101908383,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.875505255589619,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.578348466970688,1.01784466327965,-0.837364871160181) q[0];
u3(1.17072133324729,-4.50369215890248,0.777219244599163) q[2];
u3(0.147408283134750,-0.327376083195193,0.307416028410817) q[3];
u3(0.891095814206054,1.10098973343536,-1.57015239030569) q[4];
cx q[4],q[3];
u1(0.440211888533533) q[3];
u3(-1.05854273833216,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.08219197987158,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.27355525573813,-1.71405435645069,4.38477607544077) q[3];
u3(2.00961999630831,-2.39012740284731,-1.35040071400690) q[4];
u3(1.04753559011601,1.36836998720298,0.800404152577830) q[2];
u3(0.675748206600987,-0.433941596870115,-3.72351733783308) q[3];
cx q[3],q[2];
u1(-0.662680900642444) q[2];
u3(-1.79356230483659,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.997963417620631,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.683184365064943,0.742738242624286,-1.54163498999808) q[2];
u3(2.09585015882105,-0.736271554629738,0.486635759125743) q[3];
u3(0.718443918767387,2.96780294543401,-2.43196493784658) q[4];
u3(0.204107822174649,1.71657341194148,-2.84773472980062) q[0];
cx q[0],q[4];
u1(1.58833706144279) q[4];
u3(-0.601308546405931,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.199936646082755,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.33859269398580,-0.659808291255123,0.0545315275164539) q[4];
u3(2.52664950048306,-3.62597297299664,1.23046732380289) q[0];
u3(1.75129483839204,2.09622960196430,-2.76401197294182) q[3];
u3(0.990421838269183,-3.09655775106808,2.53033061814716) q[0];
cx q[0],q[3];
u1(2.46142392184447) q[3];
u3(-1.69318798635409,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.27283870870133,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.47924864193625,0.484121254059602,0.617545401355568) q[3];
u3(2.59483855351442,1.24221313871620,-0.732226597893710) q[0];
u3(0.898362766106827,-1.81220711835837,1.19114410727323) q[1];
u3(0.259648549313088,-2.39748374269539,1.76614750266781) q[2];
cx q[2],q[1];
u1(0.481215248906464) q[1];
u3(-1.33207981159528,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.02649684704975,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.87871074587489,-3.47959230304034,1.44723995171063) q[1];
u3(2.51466844612010,-1.44510863112807,0.217321512959066) q[2];
u3(1.58512599642358,2.70476425277694,-1.71846828134174) q[2];
u3(1.79721254937858,1.44367956949351,-2.38583986064036) q[1];
cx q[1],q[2];
u1(3.32851116772749) q[2];
u3(-1.33465337047619,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.25171067851434,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.70063613356309,1.28829759163372,-2.07672138963593) q[2];
u3(2.04911839034805,-1.81911617676350,-0.0129573807340176) q[1];
u3(1.15545485103934,1.10484766205898,-2.44653379521238) q[3];
u3(1.64054783554925,1.61939699888522,-4.27495544804595) q[0];
cx q[0],q[3];
u1(2.50818401917373) q[3];
u3(-2.82840261863197,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.73282862579012,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.82738681382804,-0.820383285214146,0.969123295852720) q[3];
u3(1.70720285747600,-1.92022849020278,1.37652910884769) q[0];
u3(1.84078414783867,2.08261976672904,-3.32440354327533) q[3];
u3(0.390151198533174,2.65122895026415,-1.35532174364297) q[0];
cx q[0],q[3];
u1(-0.0481951849672182) q[3];
u3(-2.42100246290617,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.43588649050539,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.29632527384311,1.83215209958526,-0.826742979277730) q[3];
u3(0.931203968125220,1.13262253746354,1.44344072630188) q[0];
u3(1.78370296974258,2.77059389956143,-1.18056067881866) q[4];
u3(0.960312812235081,1.09585003828961,-0.865959571118773) q[1];
cx q[1],q[4];
u1(-1.24338930435430) q[4];
u3(0.148915622758470,0.0,0.0) q[1];
cx q[4],q[1];
u3(3.58731757926037,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.28950797817261,2.20376767488619,-3.13089226422554) q[4];
u3(2.26508892567388,-0.539141379636114,-0.978344982948747) q[1];
u3(1.99295841701936,3.17043318120049,-2.53775134160390) q[1];
u3(1.37097075357958,2.49432224254835,-2.07674169873515) q[4];
cx q[4],q[1];
u1(0.383489523793104) q[1];
u3(-1.05982851074254,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.51173108208497,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.12401805080973,0.739059665990116,2.18214523694717) q[1];
u3(1.80340437639222,4.71987566411289,-0.0790187760781911) q[4];
u3(1.35192837166401,1.47603717993348,-3.46561294266158) q[2];
u3(0.412713476667404,2.71969828173718,-3.55144667631429) q[3];
cx q[3],q[2];
u1(0.573806373317758) q[2];
u3(-1.47290559140554,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.45670909190195,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.83890317923215,2.72803590494388,-2.11614193365470) q[2];
u3(1.40971046621115,0.330830044859580,-1.41896264880123) q[3];
u3(1.84074473276204,-1.75732001320467,4.20731758899873) q[2];
u3(0.963932521483348,-1.65280427036003,3.34984921350943) q[0];
cx q[0],q[2];
u1(3.28137146635163) q[2];
u3(-1.06995078618227,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.84225957727586,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.05505996198421,-0.896219366926963,1.87900617784717) q[2];
u3(2.20761963305132,0.222554166123950,3.64402163100431) q[0];
u3(0.739143610300719,-0.208010781322307,-0.799133491168611) q[1];
u3(1.21337314799260,-3.20879230394513,1.01210749974609) q[3];
cx q[3],q[1];
u1(0.680056609277409) q[1];
u3(-1.34948418781338,0.0,0.0) q[3];
cx q[1],q[3];
u3(-0.322631581072824,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.52105157979415,-1.19063548271188,-1.71759333282080) q[1];
u3(1.33817747916056,1.63980551382792,-2.56692410568200) q[3];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
