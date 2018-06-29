OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(2.32757736488546,1.63394782149497,-1.16522874282597) q[4];
u3(1.66774204260583,-0.443919609918639,-2.83444898784032) q[0];
cx q[0],q[4];
u1(0.00422293826500009) q[4];
u3(-1.31348081378512,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.35771621517564,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.25104207334467,2.77939043661608,-1.79948163878535) q[4];
u3(0.670765053837922,-4.07315063896383,-0.617444095224179) q[0];
u3(0.214688229611384,2.25964703516611,-1.58989074408278) q[8];
u3(0.200275420735201,2.53161360344281,-3.71514090553866) q[2];
cx q[2],q[8];
u1(2.44019958389637) q[8];
u3(-1.75985486296349,0.0,0.0) q[2];
cx q[8],q[2];
u3(3.56764752385156,0.0,0.0) q[2];
cx q[2],q[8];
u3(0.816808766052050,2.13291483343420,-2.05354340936026) q[8];
u3(1.21246408707475,1.00816327002838,-2.58731929359912) q[2];
u3(1.20848167765072,-1.07822932735302,-1.04397207190741) q[5];
u3(1.82243441244876,1.27274229043879,-4.90721487564759) q[3];
cx q[3],q[5];
u1(-0.433417170209619) q[5];
u3(-2.08192818367962,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.87447990040755,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.10972745153549,1.32552417940031,-1.95697393319033) q[5];
u3(1.50365718516280,-3.68042300737398,0.845598289689970) q[3];
u3(0.922805180169657,2.32352319896353,-1.56646641488648) q[9];
u3(1.59065720698527,0.979983404510861,-1.48044291271709) q[1];
cx q[1],q[9];
u1(-0.624096473048289) q[9];
u3(0.310446494833550,0.0,0.0) q[1];
cx q[9],q[1];
u3(3.94991755370466,0.0,0.0) q[1];
cx q[1],q[9];
u3(2.48202275968244,-1.92784427137093,-1.21111930245164) q[9];
u3(2.25713493493418,-2.00487283673930,2.55516700406821) q[1];
u3(1.64160737442256,0.207970284733870,1.48051348753879) q[10];
u3(1.62161527514544,-0.791222318062503,-1.22447276322317) q[6];
cx q[6],q[10];
u1(3.35758293527258) q[10];
u3(-4.16868243222856,0.0,0.0) q[6];
cx q[10],q[6];
u3(-0.347844427180219,0.0,0.0) q[6];
cx q[6],q[10];
u3(1.34395828152285,1.41955640473145,0.969935937498747) q[10];
u3(1.47733561729530,3.40983813027621,-1.93185118308746) q[6];
u3(0.990442885127101,-1.64881001398614,2.14188034800231) q[9];
u3(0.239240099039155,-1.21629022730058,-0.998812955651729) q[8];
cx q[8],q[9];
u1(1.36554456818942) q[9];
u3(-2.42896274831966,0.0,0.0) q[8];
cx q[9],q[8];
u3(0.116800512997101,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.17940385491265,-2.56341103600146,3.64747819205446) q[9];
u3(1.20972594612528,0.458426872011377,1.88866085110541) q[8];
u3(1.02750533332856,-2.11215901839568,0.650350108535624) q[1];
u3(0.564955713774224,-2.46174862334313,0.970897664665903) q[3];
cx q[3],q[1];
u1(1.09575600911933) q[1];
u3(-0.685057290184266,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.60897240430024,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.69194768735997,1.27851473920683,0.704104726812075) q[1];
u3(2.64483021284802,-2.09368965602436,-2.07628333779459) q[3];
u3(2.45384165493386,-1.35303555252171,2.50537814033997) q[6];
u3(2.35707651572318,-2.38663466726232,-0.226842688239426) q[0];
cx q[0],q[6];
u1(1.53981998985778) q[6];
u3(-0.341996655903383,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.09592051887243,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.29685195925127,3.05825195048775,0.992289513175232) q[6];
u3(1.41567778471861,2.82669465633438,1.58676589141421) q[0];
u3(0.732510925371167,-0.166514884989081,0.0762643896878186) q[7];
u3(2.22644421280616,-3.53916877131357,-0.525031571783329) q[10];
cx q[10],q[7];
u1(3.01433923098953) q[7];
u3(-1.66745330020580,0.0,0.0) q[10];
cx q[7],q[10];
u3(0.793072303718015,0.0,0.0) q[10];
cx q[10],q[7];
u3(2.37286387541063,1.33492482138254,-4.78665039973540) q[7];
u3(1.26840857202373,-5.31738036610029,-0.177037179690426) q[10];
u3(1.98296395034532,0.287614380919250,1.62779816139383) q[5];
u3(1.46526465936480,-1.16269047095107,-1.05585991678234) q[4];
cx q[4],q[5];
u1(3.09405058146730) q[5];
u3(-0.577567785337841,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.90468969700225,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.53365554474733,3.36802664448754,-1.41296874177746) q[5];
u3(0.560895249973500,-0.981797335312516,-1.76065931227114) q[4];
u3(1.88474150556793,-1.91346150294414,-0.290614601757436) q[10];
u3(2.68288838162575,-1.84295897302388,-0.224063196632832) q[4];
cx q[4],q[10];
u1(2.64349127104110) q[10];
u3(-1.98149873922299,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.298516571118066,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.74215624892219,-2.50449370402334,1.56965248177918) q[10];
u3(0.542307733807086,1.78948266440540,0.469867593002992) q[4];
u3(1.76438654011888,0.518484628858519,-3.33802826061689) q[6];
u3(1.30628534054775,2.70937462488520,-2.79976488205782) q[0];
cx q[0],q[6];
u1(3.22291056871026) q[6];
u3(-1.47791462953065,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.79164729383010,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.49413960883184,-1.10743513235769,0.0731507028053620) q[6];
u3(1.68433452768145,3.26447249765664,-0.712329033727395) q[0];
u3(2.14274192834274,-0.0331672134625174,2.69789563738069) q[9];
u3(2.92256826704493,-1.52767979514362,-1.05733399900752) q[7];
cx q[7],q[9];
u1(2.27122171723579) q[9];
u3(-1.83511656978031,0.0,0.0) q[7];
cx q[9],q[7];
u3(0.418830676508724,0.0,0.0) q[7];
cx q[7],q[9];
u3(1.38655716674113,-3.53706380544651,2.53307885670796) q[9];
u3(1.16313501034760,4.35023828910755,1.32569815158361) q[7];
u3(1.74558963634698,2.48963887653326,-2.13275349388731) q[1];
u3(0.887056245233295,1.77667364984283,-2.60204364788942) q[3];
cx q[3],q[1];
u1(3.52351705473358) q[1];
u3(-0.715099101080231,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.85876616456147,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.61967645256828,-1.83447865325231,-1.06528754740775) q[1];
u3(1.27266293334784,-0.429006000398081,-1.81809939551684) q[3];
u3(2.57676416893708,0.717252584888796,-1.05331621933222) q[8];
u3(1.55640261270540,0.0440050971302470,-3.07476380454515) q[2];
cx q[2],q[8];
u1(0.268647547559933) q[8];
u3(-1.06278320919388,0.0,0.0) q[2];
cx q[8],q[2];
u3(1.72358764776638,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.51964047692194,-3.65293240278202,0.811741816211001) q[8];
u3(1.42138470698058,2.45055047478380,-0.181719962323567) q[2];
u3(0.929854215206651,-3.54091831719817,2.44595345389697) q[4];
u3(1.79473908604892,2.73312285332707,-3.30106618776602) q[9];
cx q[9],q[4];
u1(1.23346011906285) q[4];
u3(-0.207199291211945,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.48009393414784,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.62954894288741,-1.54997538702168,1.01484992040099) q[4];
u3(1.75851609847167,-2.70474287931879,-0.859544399732125) q[9];
u3(1.96062614564732,2.93453513240770,-1.66330570092540) q[1];
u3(1.71399214950841,1.05809374380918,-2.73043408918815) q[2];
cx q[2],q[1];
u1(1.65275224197766) q[1];
u3(-0.290397575735243,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.07203263101769,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.782236746931291,-4.36764783633467,1.56195438178975) q[1];
u3(2.94954324250838,0.325955157410849,-0.0130406188184099) q[2];
u3(1.74474166050581,0.856485277375331,0.873298801525445) q[3];
u3(0.695792708057860,-0.533553607681728,-2.64161258349835) q[7];
cx q[7],q[3];
u1(1.50481740154734) q[3];
u3(-0.796254660451966,0.0,0.0) q[7];
cx q[3],q[7];
u3(2.48960051558526,0.0,0.0) q[7];
cx q[7],q[3];
u3(0.904727258316314,1.14679706843101,1.11647748462262) q[3];
u3(2.61323306810521,-5.27051863436012,-0.942191200287410) q[7];
u3(0.913744211897129,1.84161759783095,-1.64179661884823) q[6];
u3(0.870752832177949,-0.746447731154240,0.0185453283264522) q[8];
cx q[8],q[6];
u1(0.946532902194360) q[6];
u3(-0.359426319522273,0.0,0.0) q[8];
cx q[6],q[8];
u3(1.60622011116662,0.0,0.0) q[8];
cx q[8],q[6];
u3(1.71904292421244,2.95699932222862,-2.58698177328230) q[6];
u3(1.45966633102970,-2.30726555362965,-3.03045955918213) q[8];
u3(0.727689734756531,-2.36511824314475,3.04798581240467) q[10];
u3(0.952473506832636,-3.18061637130515,1.11649366854852) q[0];
cx q[0],q[10];
u1(3.30721052887751) q[10];
u3(-0.758533126227509,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.71250845472425,0.0,0.0) q[0];
cx q[0],q[10];
u3(0.885668994136885,0.0977684179614233,-2.94197840610394) q[10];
u3(1.81736860498587,-0.432679814926832,-2.34695046661182) q[0];
u3(0.879365008668941,1.73517977417868,-0.162514651461827) q[1];
u3(1.30003554467888,0.196340183684281,-4.34485452203070) q[10];
cx q[10],q[1];
u1(0.743504279853455) q[1];
u3(-0.139405146569714,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.83343962393134,0.0,0.0) q[10];
cx q[10],q[1];
u3(1.43978903506653,0.0349756641763546,4.37087838022717) q[1];
u3(2.49595227446587,-1.32456440116886,1.94116619938489) q[10];
u3(1.87685287161536,-0.910262943503071,0.652610520830053) q[8];
u3(0.998638729237785,-4.17439201625901,-0.680981047527006) q[2];
cx q[2],q[8];
u1(2.09916769486366) q[8];
u3(-2.67851628513846,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.653218190606166,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.31243735725984,3.24689732284100,-0.907628108335811) q[8];
u3(0.659713081677832,1.84242890092922,1.22480634240420) q[2];
u3(2.22648616450260,2.59768107441446,-2.54859259906185) q[3];
u3(1.66241285294363,1.85306972699556,-2.23567275673489) q[7];
cx q[7],q[3];
u1(0.896100541510071) q[3];
u3(-1.41048677612070,0.0,0.0) q[7];
cx q[3],q[7];
u3(-0.0717410100368343,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.18618375974036,-2.96683032686734,1.97060031162703) q[3];
u3(1.69483675223946,-2.42619754377083,-0.557580586426028) q[7];
u3(1.60927498144376,0.300671536068581,-2.26854553133669) q[6];
u3(1.68806136215804,2.99721087973172,-3.14507136398605) q[0];
cx q[0],q[6];
u1(3.81107459081017) q[6];
u3(-3.55090929382887,0.0,0.0) q[0];
cx q[6],q[0];
u3(-0.887571838816410,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.48614231481517,-2.55903559584713,3.69732689614487) q[6];
u3(1.50814991821251,3.25380635530719,-0.868943967353367) q[0];
u3(0.853097469357040,0.00897853024899820,-2.50603049876590) q[9];
u3(2.69222167048738,-2.37426554682411,3.38106240676745) q[5];
cx q[5],q[9];
u1(3.03148407816560) q[9];
u3(-1.95817090075625,0.0,0.0) q[5];
cx q[9],q[5];
u3(0.872452835836601,0.0,0.0) q[5];
cx q[5],q[9];
u3(1.65872171942531,-0.0603906997079894,1.12375774324681) q[9];
u3(1.11262871647153,-0.521771707365700,-3.67031766907006) q[5];
u3(1.92347230988180,3.60852586292827,-0.993781322139935) q[2];
u3(2.50341722390081,2.64223751962530,-1.18231574870086) q[0];
cx q[0],q[2];
u1(3.91249337906600) q[2];
u3(-4.35706118539820,0.0,0.0) q[0];
cx q[2],q[0];
u3(-0.629888641732011,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.88513700902323,-2.20494100109638,1.74680748531496) q[2];
u3(0.742845371177522,-1.77352384586821,-4.16363877712507) q[0];
u3(1.31646475514521,0.838265071664873,0.812477562832897) q[5];
u3(1.18485651276303,-0.0288138401481521,-2.19201775724609) q[1];
cx q[1],q[5];
u1(-0.266085887942010) q[5];
u3(-2.22332536837457,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.64189152275356,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.34572668160891,-2.10266175381008,2.15829513182297) q[5];
u3(2.65618801255783,-0.613969304585789,-1.10651341921409) q[1];
u3(1.92259946466194,-1.02136293376604,1.73611536448607) q[10];
u3(2.17287781386119,-1.98219899592096,-2.25216121330161) q[3];
cx q[3],q[10];
u1(1.73339827733531) q[10];
u3(-3.01467089760167,0.0,0.0) q[3];
cx q[10],q[3];
u3(1.15526773583699,0.0,0.0) q[3];
cx q[3],q[10];
u3(2.60382107751132,-0.307669411781022,0.00180324847894830) q[10];
u3(2.81471899725418,-3.07326679994066,2.86711527511130) q[3];
u3(1.90935200952293,-0.773749621609843,1.27100892276729) q[7];
u3(1.92328080118647,-1.84559849843876,-0.402202749628457) q[4];
cx q[4],q[7];
u1(0.820589024795445) q[7];
u3(-0.474122608338719,0.0,0.0) q[4];
cx q[7],q[4];
u3(1.82087251148555,0.0,0.0) q[4];
cx q[4],q[7];
u3(2.01591170373036,-3.37355287585091,2.57833394745846) q[7];
u3(1.39962508496318,0.305890329719287,-0.947530610662957) q[4];
u3(1.55830537593302,3.88439502812269,-1.29683084355955) q[9];
u3(0.398005685414626,1.12812151551954,-0.495200684535946) q[8];
cx q[8],q[9];
u1(1.31380333533672) q[9];
u3(-1.00786392129463,0.0,0.0) q[8];
cx q[9],q[8];
u3(-0.0478910709738005,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.94689277620214,0.408523252569311,-2.01598523838078) q[9];
u3(2.00829879205087,-0.665636886872649,-3.89301021413636) q[8];
u3(2.53117884937314,1.67070399936384,-2.36709245764494) q[5];
u3(1.82272564241046,-2.16469961414361,2.14112523381382) q[1];
cx q[1],q[5];
u1(3.49064301799063) q[5];
u3(-3.23293557222778,0.0,0.0) q[1];
cx q[5],q[1];
u3(-1.01339458727839,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.678914735594723,1.93484983404439,-1.12002516188034) q[5];
u3(1.47450307115039,-2.78887738117909,0.150469900845073) q[1];
u3(2.20240526070426,-1.56404465887607,0.652665324133489) q[2];
u3(1.74449018806363,-3.55168577370509,0.569668160246236) q[0];
cx q[0],q[2];
u1(2.32305319469170) q[2];
u3(-1.91903934682434,0.0,0.0) q[0];
cx q[2],q[0];
u3(3.50723288627713,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.20443353687480,1.74384366415398,-0.554572458309217) q[2];
u3(1.48179286586567,2.04809458500027,2.90660110584451) q[0];
u3(1.81628868676784,2.99652240947097,-0.848381456297151) q[8];
u3(0.691002730648217,0.859427042400391,-0.542769720165645) q[6];
cx q[6],q[8];
u1(3.03102982371433) q[8];
u3(-2.24882610722152,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.47886575065322,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.88997169256944,1.29185393738493,-0.643665633351347) q[8];
u3(1.17604157245026,2.54261004120436,3.20368653627539) q[6];
u3(0.995711089929313,-1.67745792552097,1.07222918801823) q[9];
u3(1.82194652578918,-2.23722304904908,0.884244018187657) q[3];
cx q[3],q[9];
u1(2.60445352278171) q[9];
u3(-1.71679999547647,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.01170846218136,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.32042764834852,-1.23738460429061,0.558715300883515) q[9];
u3(1.12539163054767,-0.821434990681140,-1.53136631417815) q[3];
u3(2.40565555461281,2.69469873347739,-1.41435845496703) q[4];
u3(1.79326128153024,0.625836785521243,-2.08684546179712) q[10];
cx q[10],q[4];
u1(0.630523678794164) q[4];
u3(-1.28418647138434,0.0,0.0) q[10];
cx q[4],q[10];
u3(2.83670317210190,0.0,0.0) q[10];
cx q[10],q[4];
u3(1.63546551093967,0.832744148041828,-2.07674448378189) q[4];
u3(1.20147405004261,5.23999411667316,-0.0884426180270972) q[10];
u3(1.05157527145950,1.55422743109097,-3.32545502961057) q[5];
u3(1.20378569337958,-1.95043886008660,3.78050844247955) q[1];
cx q[1],q[5];
u1(1.53263648626940) q[5];
u3(-3.24544322654116,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.09950245087662,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.92702572626296,1.42944065413009,-0.0117542491424544) q[5];
u3(2.27907782037888,1.96695568868163,-2.10633571248335) q[1];
u3(1.73627557757301,0.0481141691207413,-0.932732625493183) q[7];
u3(1.64048120428679,-3.47309656065856,0.827754722526745) q[3];
cx q[3],q[7];
u1(1.33978768119255) q[7];
u3(-3.56553866169297,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.34793605200855,0.0,0.0) q[3];
cx q[3],q[7];
u3(2.28140168883812,1.54734327339037,-3.73604185208694) q[7];
u3(1.80078462009546,-0.578717300762113,4.43583405476451) q[3];
u3(1.22893591888268,0.642557671146006,-0.866264081675737) q[6];
u3(1.01407376137991,0.198921075773907,-3.01510741786165) q[10];
cx q[10],q[6];
u1(1.70592157652952) q[6];
u3(-2.82446248279908,0.0,0.0) q[10];
cx q[6],q[10];
u3(0.594331919284634,0.0,0.0) q[10];
cx q[10],q[6];
u3(1.22812583600872,-0.229590372212330,3.51706824665437) q[6];
u3(1.61019521990214,-1.05758192420707,1.39473161531652) q[10];
u3(1.84766008419483,2.98228837266449,-1.96348374238018) q[2];
u3(2.20215138929896,2.26062558767633,-2.11852822741581) q[4];
cx q[4],q[2];
u1(1.54486844595331) q[2];
u3(-2.92436578907586,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.272137935165687,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.51907600597694,-3.11016786100601,0.721793913782265) q[2];
u3(1.49234120795757,6.01851417884102,-0.183862045059409) q[4];
u3(1.39789620181868,1.14135582458993,-3.48482898901686) q[8];
u3(2.01415304327958,2.95504399697943,-2.47882898357606) q[0];
cx q[0],q[8];
u1(3.04246312927794) q[8];
u3(-1.33504632025877,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.36490617557556,0.0,0.0) q[0];
cx q[0],q[8];
u3(2.64457099037946,-2.03360060299134,0.829787133826506) q[8];
u3(1.81422839872139,-2.10600556121165,-1.11473297944118) q[0];
u3(1.08781758693188,-0.991365763718028,2.08923758493427) q[9];
u3(1.70315128292144,-1.35172083492397,-2.18709932891903) q[5];
cx q[5],q[9];
u1(1.61833215064189) q[9];
u3(-1.11951040051760,0.0,0.0) q[5];
cx q[9],q[5];
u3(-0.664384245667351,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.09253331823941,2.32055310423529,-1.37491749822681) q[9];
u3(2.19847826722829,1.36486827284100,-0.974861978596222) q[5];
u3(2.39259740338605,-3.73645550712401,2.06123576643281) q[10];
u3(0.340845398523844,-1.48667578823424,2.93901491857121) q[1];
cx q[1],q[10];
u1(-0.185490312721109) q[10];
u3(-1.80013955757506,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.995791741377439,0.0,0.0) q[1];
cx q[1],q[10];
u3(1.38921684679736,3.10343946448458,-0.0303795653354235) q[10];
u3(1.21673281520768,-0.765107463376074,3.70548287117458) q[1];
u3(0.578687321475422,1.92856724126087,-2.38073775768957) q[6];
u3(0.506876342880215,2.06039954560873,-3.77782256096491) q[8];
cx q[8],q[6];
u1(1.68202336408830) q[6];
u3(-0.250590318030941,0.0,0.0) q[8];
cx q[6],q[8];
u3(2.60349965096427,0.0,0.0) q[8];
cx q[8],q[6];
u3(0.896634062691421,1.41883356689377,-0.0925749263288637) q[6];
u3(1.93803098457647,0.00685544489918222,2.92898497049183) q[8];
u3(1.12497391659071,2.93913340594369,-1.62764985256851) q[2];
u3(1.25352439408393,1.09455120627326,-2.55087810535289) q[7];
cx q[7],q[2];
u1(1.55748292295107) q[2];
u3(0.533382249362153,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.922482587374473,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.28327045844342,-1.58079939198707,1.18662819115176) q[2];
u3(2.21589798260162,-2.38324109076956,-3.79001963747301) q[7];
u3(1.35745355752837,3.13351675029180,-2.02654055626848) q[3];
u3(2.49049721301528,0.352526182289063,-1.88195661383907) q[0];
cx q[0],q[3];
u1(1.04004165826198) q[3];
u3(-0.403789143926970,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.20347259348725,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.25231556254922,0.192911236301314,-1.99368498239853) q[3];
u3(0.592608649798990,2.38336944150678,-2.96428144614572) q[0];
u3(1.89613386503395,0.0516058036758308,1.25322481435666) q[0];
u3(2.08866992721892,-1.20187702631198,-1.95001494423285) q[10];
cx q[10],q[0];
u1(1.27511733274775) q[0];
u3(-0.905306918021512,0.0,0.0) q[10];
cx q[0],q[10];
u3(3.24069978661190,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.12472305962393,4.20234607366443,-0.620559465649740) q[0];
u3(1.87926144666111,1.82757783081602,-1.24761608781351) q[10];
u3(0.796326048227910,2.55630287290681,-2.95537988274873) q[2];
u3(0.826573961434924,2.10373596072902,-3.67739188476589) q[9];
cx q[9],q[2];
u1(4.47030412311767) q[2];
u3(-2.68012956749592,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.528286299286391,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.28062649618782,0.140603760991167,-2.32525912216924) q[2];
u3(1.64166808624919,-1.81006767904779,-1.65948696184635) q[9];
u3(2.55980683634710,2.08579880484735,0.668687150098313) q[4];
u3(1.20674230857473,-0.368212106947974,-3.01865324927372) q[3];
cx q[3],q[4];
u1(1.05186675786184) q[4];
u3(-1.40319715700828,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.54502470597313,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.35068362936099,2.05128779209808,-0.331274707970629) q[4];
u3(1.68049850836386,-0.353701333465220,0.977714292573913) q[3];
u3(0.618938794357048,0.233055961751656,0.191453406961408) q[6];
u3(0.659508292385104,-1.87688804541219,1.44114821716015) q[5];
cx q[5],q[6];
u1(0.832044529916854) q[6];
u3(-0.564921804001147,0.0,0.0) q[5];
cx q[6],q[5];
u3(2.12796887855597,0.0,0.0) q[5];
cx q[5],q[6];
u3(0.882639427772471,2.16647598385763,-3.48649884741533) q[6];
u3(1.05777788863076,-0.192841886658244,-2.70533815678597) q[5];
u3(1.43240021473012,3.53547151304326,-1.70199311498711) q[8];
u3(0.640117900296088,1.67975232078040,-0.663165247193141) q[1];
cx q[1],q[8];
u1(-0.258317803368453) q[8];
u3(0.345335887908867,0.0,0.0) q[1];
cx q[8],q[1];
u3(3.95674125796737,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.96425511774427,-2.88671585347089,2.41377790108001) q[8];
u3(1.58712531577188,-0.570444085373537,-0.446838511854574) q[1];
u3(0.114637948049242,-0.132808824015262,0.426193041453881) q[4];
u3(0.714449189129452,0.445908758626421,-1.63303711591916) q[5];
cx q[5],q[4];
u1(1.52167455847365) q[4];
u3(0.0916486526500018,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.17259029717653,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.60489359847057,-0.250760795817325,2.23810349166140) q[4];
u3(1.85833683424297,3.20176136610903,2.66588123720324) q[5];
u3(1.94267000721211,3.41364808130459,-2.23315826576638) q[8];
u3(1.23507576172684,1.03565244375701,-2.02200920111548) q[3];
cx q[3],q[8];
u1(-0.366356300474144) q[8];
u3(0.948324971240952,0.0,0.0) q[3];
cx q[8],q[3];
u3(4.02056925717257,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.33246040498392,2.44819152658157,-2.09807475343623) q[8];
u3(2.15789802870199,-2.46550221815928,-3.26403965330397) q[3];
u3(0.967060162834763,2.45868917420607,-3.14590810150161) q[10];
u3(1.83081383707965,-2.87880068234584,3.07864762631700) q[1];
cx q[1],q[10];
u1(2.99026992196988) q[10];
u3(-2.79940795853359,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.681285629291962,0.0,0.0) q[1];
cx q[1],q[10];
u3(0.859260633297380,-3.51034357364849,-0.199199011100767) q[10];
u3(2.23566114861659,0.369309548432843,5.56730274466516) q[1];
u3(0.539599146916626,1.48938380230157,-1.86556689140756) q[0];
u3(1.86686233750583,1.97183969188304,-4.19344476466457) q[7];
cx q[7],q[0];
u1(3.18973541106026) q[0];
u3(-0.989373230490857,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.81174013738009,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.18087742067404,3.65024977166401,-1.73302497143473) q[0];
u3(1.61815985615398,0.450167797371406,-2.87423887826981) q[7];
u3(1.38274039247567,2.02519077861341,-2.22920181465133) q[9];
u3(0.894456699030249,2.68787428799949,-3.13176838891378) q[6];
cx q[6],q[9];
u1(3.15818852818299) q[9];
u3(-1.00566146844623,0.0,0.0) q[6];
cx q[9],q[6];
u3(2.36400251573581,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.46007288687035,-1.44885364479578,-0.729122288405631) q[9];
u3(0.769447195173669,0.306066286725878,-1.17758019147302) q[6];
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