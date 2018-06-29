OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(1.50789219470150,3.09293578660134,-0.455618704121902) q[3];
u3(2.73035257527671,0.149819691655204,-2.82891021081488) q[0];
cx q[0],q[3];
u1(-0.0604490603530723) q[3];
u3(-1.94128548115097,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.514375391906346,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.613489350872380,1.13533197903161,-2.09448024192857) q[3];
u3(1.05750237226182,-2.60709794368052,-0.690490114681800) q[0];
u3(1.14376524578284,2.86447433893212,-2.65446722747076) q[7];
u3(0.578456846821180,2.94876898093635,-3.13753849442274) q[2];
cx q[2],q[7];
u1(2.52564370711071) q[7];
u3(0.192317986648847,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.57483425226544,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.448736142435748,-2.50016716322395,2.51477544983445) q[7];
u3(2.54284497952218,-0.673465899620364,-2.86928762773325) q[2];
u3(1.35748466459897,0.748812594999519,-3.68183623648628) q[6];
u3(2.26405435691198,-1.78234490569067,4.07993543916740) q[1];
cx q[1],q[6];
u1(-0.0331913649956728) q[6];
u3(-0.780879756363768,0.0,0.0) q[1];
cx q[6],q[1];
u3(1.76322123720875,0.0,0.0) q[1];
cx q[1],q[6];
u3(1.92809817932141,-2.79435549095884,3.35786868302774) q[6];
u3(0.935360813453582,1.87832630616850,-0.987439365215229) q[1];
u3(1.61714416881662,0.0972089093107000,0.593224676616404) q[5];
u3(0.752096680870024,-0.0269302781837224,-5.02174449887905) q[4];
cx q[4],q[5];
u1(1.41539141252513) q[5];
u3(-3.75164601703625,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.35774679281029,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.16196563288161,4.08114959793491,-0.224908502526885) q[5];
u3(0.811518528037585,-0.124801637185960,5.22396759281856) q[4];
u3(2.23497057427551,-1.21301463490136,2.85247508942296) q[5];
u3(2.86724451203925,-1.21580278452356,0.926282112398996) q[6];
cx q[6],q[5];
u1(1.90127013992123) q[5];
u3(-0.0457724357854341,0.0,0.0) q[6];
cx q[5],q[6];
u3(0.465235424821901,0.0,0.0) q[6];
cx q[6],q[5];
u3(0.976813509949707,1.57363728420627,-0.240875118529610) q[5];
u3(2.18836150863253,0.458476991247392,1.32922015913859) q[6];
u3(1.47025455969126,0.511402241932898,0.916650955542797) q[0];
u3(1.63991658673843,-0.451801817694260,-3.05772519156460) q[1];
cx q[1],q[0];
u1(3.00233101702365) q[0];
u3(-2.11091412460993,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.45063737675732,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.90502887620419,-0.764276404117708,2.14673255303738) q[0];
u3(1.85813343949655,-0.226636716750533,0.501554037553809) q[1];
u3(1.27835215924014,3.58503383633433,-0.812358203683038) q[4];
u3(0.524138016009329,1.96616678957148,-1.67526410870752) q[2];
cx q[2],q[4];
u1(2.32764322815498) q[4];
u3(-1.78759840764735,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.37509179090593,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.83142969206388,-1.40318450366536,-1.37855478422715) q[4];
u3(1.86639686905291,0.993370405909821,-3.65545347994042) q[2];
u3(1.45915012360474,1.97933620482100,-1.94344206249918) q[3];
u3(0.435087374671055,2.42201403955069,-3.22945361313463) q[7];
cx q[7],q[3];
u1(1.71944958522819) q[3];
u3(-0.101972099291388,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.03992288776971,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.73403302742658,-2.37832902879041,0.927465482465395) q[3];
u3(1.12387412877847,2.79136576203244,-2.96605843134013) q[7];
u3(1.65114979266035,-0.947403179155978,-1.19625694683580) q[6];
u3(1.37140923411869,-4.07391927983748,0.907130095853092) q[2];
cx q[2],q[6];
u1(0.952428775835944) q[6];
u3(-0.780979008848499,0.0,0.0) q[2];
cx q[6],q[2];
u3(1.87244610121923,0.0,0.0) q[2];
cx q[2],q[6];
u3(0.693520668072923,1.40028482538807,-1.36614102150576) q[6];
u3(0.684102142758862,1.19411016256579,-2.86954775055117) q[2];
u3(0.728591543909007,2.24663854785808,-1.61915665841421) q[3];
u3(0.593698902746587,-2.80905547208499,1.79108767823298) q[4];
cx q[4],q[3];
u1(1.62019265738639) q[3];
u3(-0.848293583555370,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.41083223366148,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.86547929220753,-0.982355090296224,-0.643826632072749) q[3];
u3(0.516709840496290,-1.74781086068189,3.18331031517177) q[4];
u3(1.79981022961506,0.000119552092463726,-2.01505478680501) q[1];
u3(1.29929707340698,-3.44681030783948,2.07003073668903) q[5];
cx q[5],q[1];
u1(1.32425228702829) q[1];
u3(-0.523818731078063,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.19694646762713,0.0,0.0) q[5];
cx q[5],q[1];
u3(2.94216045237776,3.44339695144215,-2.56591207725395) q[1];
u3(1.52021174950933,-0.306905110437919,5.15872953305136) q[5];
u3(1.57344733142092,0.375153130047897,2.60105160923574) q[0];
u3(2.13769354046206,-2.81967454564404,-2.06249275837048) q[7];
cx q[7],q[0];
u1(1.09046701919713) q[0];
u3(-3.23197822262268,0.0,0.0) q[7];
cx q[0],q[7];
u3(1.66550457073072,0.0,0.0) q[7];
cx q[7],q[0];
u3(1.68981204684502,-1.67036294838108,3.69027319200449) q[0];
u3(1.61182549944928,-1.51881149081551,4.24842656340827) q[7];
u3(0.793789834852790,2.05831565281192,0.938594738217660) q[2];
u3(1.17655159639571,-0.0318698635726682,-3.33262898162406) q[7];
cx q[7],q[2];
u1(3.38542129307151) q[2];
u3(-0.620774801796776,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.89872370770274,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.53976699144142,-2.13992951629859,-0.322764268736101) q[2];
u3(1.65246086271690,0.251923164332345,-3.16233686728038) q[7];
u3(0.815699691407799,0.281639238166151,-0.916298154340717) q[1];
u3(0.906228460328951,-3.45256193980181,1.36943412046467) q[5];
cx q[5],q[1];
u1(2.21065550485101) q[1];
u3(-3.16321757152822,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.913736563132597,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.987329918290216,-2.13305112662322,2.72334161917783) q[1];
u3(2.34817524196446,-4.54709837017373,0.464699514721388) q[5];
u3(2.44378118951685,-0.173285834980640,2.73085606363774) q[6];
u3(2.63177657442830,0.977165172772392,2.15069951533553) q[0];
cx q[0],q[6];
u1(0.988085257563013) q[6];
u3(-3.11804129890916,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.98093697067943,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.664027988945204,-1.04836997252602,-2.69004446580604) q[6];
u3(1.91482382607916,-1.92592941261284,-0.804785221062549) q[0];
u3(1.69210610082847,-2.09082280038649,-0.385247394660584) q[3];
u3(1.73905240428324,-3.64596025513916,1.00226028971341) q[4];
cx q[4],q[3];
u1(1.98117479115569) q[3];
u3(-2.58105565842057,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.771350950900020,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.41565783607832,0.0504273219012725,-2.92862461721923) q[3];
u3(2.23502439711482,2.24629397499933,-0.00328494230807119) q[4];
u3(1.14931601133090,0.420707696329127,0.536242406944887) q[4];
u3(1.89986357154150,-0.0687626315083456,-2.09481595230415) q[6];
cx q[6],q[4];
u1(1.60992046508387) q[4];
u3(-2.87506869915400,0.0,0.0) q[6];
cx q[4],q[6];
u3(0.998104322602593,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.256266791338718,1.50578780670281,-0.493410538705009) q[4];
u3(2.35812713356340,0.473101379359086,0.556860342228559) q[6];
u3(0.687066709496210,-1.44050279674029,0.763301377105744) q[7];
u3(0.771571322246815,-2.96641819994256,1.38421214389409) q[5];
cx q[5],q[7];
u1(1.55304273062592) q[7];
u3(-3.85503191350842,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.16254178267591,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.98285950271105,0.282938364269137,-2.27729144808479) q[7];
u3(2.28285835877933,-6.15253522315714,0.0322581739495580) q[5];
u3(2.34450503757767,1.80525096143550,-4.01145611753044) q[2];
u3(1.54311419103988,-1.96716270607530,3.93254365680206) q[0];
cx q[0],q[2];
u1(3.12000708594998) q[2];
u3(-1.58771332728732,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.826504230852324,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.70946075202963,-2.95825996958897,0.0211284800774367) q[2];
u3(2.07273296201674,-1.07295704203740,-2.81069661693208) q[0];
u3(1.26859093363165,-1.52799026820205,-0.418360997742643) q[3];
u3(1.37111768673015,-1.77135728274348,-0.600461322433183) q[1];
cx q[1],q[3];
u1(3.25464903624126) q[3];
u3(-3.70498218723183,0.0,0.0) q[1];
cx q[3],q[1];
u3(-0.876192703003989,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.23218257314468,-4.52605657324855,0.439384360486311) q[3];
u3(0.576579040919946,3.62538029584558,-0.100033695245976) q[1];
u3(2.21418862690467,3.18609882404699,-2.41090984837916) q[5];
u3(0.291910943382253,-1.12694229416101,2.81062170009428) q[1];
cx q[1],q[5];
u1(1.22398074740615) q[5];
u3(-0.255437534451741,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.43859571720782,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.42627192024315,0.984610664507605,-4.57735906419518) q[5];
u3(1.72638953599058,4.31661549967559,1.59944464212125) q[1];
u3(1.76013963745224,0.143213432552449,-2.41212553092811) q[0];
u3(2.03056482573902,1.93901918569627,-4.00650357428288) q[3];
cx q[3],q[0];
u1(3.69696135210846) q[0];
u3(-4.17046436571724,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.108149121733156,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.853458281882412,0.711388650739253,1.80487690254706) q[0];
u3(1.43645994674879,-1.28175216756518,-2.02911946008307) q[3];
u3(0.851242108260269,1.77457644151616,-1.32502809349062) q[2];
u3(0.638438785141696,-0.308214907924077,-2.11824044859502) q[4];
cx q[4],q[2];
u1(3.25727334512944) q[2];
u3(-0.823083879614414,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.85783617705887,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.340978348237602,-2.05410812971555,-0.395235837672987) q[2];
u3(1.78392684701157,-3.79847840117181,-1.99684260751480) q[4];
u3(0.553581917122993,-0.601879175820439,0.640559749434308) q[7];
u3(0.544286003229022,-3.24141240763732,0.310839752001654) q[6];
cx q[6],q[7];
u1(0.340390276725941) q[7];
u3(-1.59095692127631,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.51187913063702,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.97729590382432,1.19089740676362,0.193595541680804) q[7];
u3(2.00138087682109,3.31656438009749,-0.632979091083290) q[6];
u3(2.25379903271427,-0.475836612460467,-1.15668189064524) q[6];
u3(0.752054663998403,-0.525355159787353,-4.12605681466711) q[7];
cx q[7],q[6];
u1(2.09336544533319) q[6];
u3(0.227122163478082,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.86094574213214,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.12877612585257,-2.25399644049924,2.23702487527250) q[6];
u3(0.770678753751050,1.30008752072789,-4.23736355961819) q[7];
u3(2.55136952184857,3.04736994012563,-1.50128792459384) q[1];
u3(1.46754737673623,1.98052241400769,-3.01167374651471) q[5];
cx q[5],q[1];
u1(1.65295837633545) q[1];
u3(-3.06205506549055,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.700775842981002,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.59156531352509,1.02938348823927,0.135203129179961) q[1];
u3(1.39032731746896,-3.33768485914718,2.32197132249334) q[5];
u3(0.434208946755898,-1.10268711764818,-0.460416338985357) q[3];
u3(1.77608701627680,0.701811551699302,-4.51913692255282) q[4];
cx q[4],q[3];
u1(4.06761297556275) q[3];
u3(-3.21020899662992,0.0,0.0) q[4];
cx q[3],q[4];
u3(-0.680210139577985,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.55627694301769,2.13958136598811,1.07568196204218) q[3];
u3(1.59435879137601,1.47360643848988,0.577944225772759) q[4];
u3(1.38693954089345,0.451584100140251,-0.911914633576627) q[2];
u3(0.443931931679150,0.00408756510539088,-3.17349853197587) q[0];
cx q[0],q[2];
u1(0.839627474240579) q[2];
u3(-0.446323241610374,0.0,0.0) q[0];
cx q[2],q[0];
u3(3.07789864039258,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.16821795232689,-3.58004481761860,0.320417751009415) q[2];
u3(2.94876609226679,0.781823313753329,-3.90329714449405) q[0];
u3(2.48637133137692,3.31321433337389,-2.83672514936984) q[1];
u3(1.28626244182963,2.91027500357220,-1.79844796423161) q[2];
cx q[2],q[1];
u1(1.71899188527003) q[1];
u3(-0.0807079154816235,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.53964496630522,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.92658361130087,1.10836591447572,0.182201250727858) q[1];
u3(1.96477092736658,-0.557738654723319,-3.23855235162374) q[2];
u3(2.52521578186509,0.906528083209964,-0.978791652188212) q[6];
u3(2.12483061873941,-4.20732276130472,1.35880088628533) q[5];
cx q[5],q[6];
u1(0.844107134659206) q[6];
u3(-1.67383643866343,0.0,0.0) q[5];
cx q[6],q[5];
u3(-0.393582305173417,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.88852314638437,1.69661202285798,-3.66927285700338) q[6];
u3(2.07808702030141,-1.48871780222859,0.461328359644372) q[5];
u3(2.46608065855922,1.01796590560822,0.221264101172869) q[4];
u3(1.88541816974492,-0.479180921765235,-3.39303941192426) q[7];
cx q[7],q[4];
u1(0.266023531733118) q[4];
u3(-1.02797373619070,0.0,0.0) q[7];
cx q[4],q[7];
u3(2.36092489993589,0.0,0.0) q[7];
cx q[7],q[4];
u3(2.60641438535312,-1.61442666688935,1.39568622601989) q[4];
u3(0.818510839778130,-1.99986782147883,0.0812197432420732) q[7];
u3(0.684580788041992,-1.98104919296682,-1.00950094022588) q[0];
u3(1.46961398220014,-2.97809053362709,0.0296018997463710) q[3];
cx q[3],q[0];
u1(1.04681724372668) q[0];
u3(-0.715008604299013,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.94319163423318,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.532117075156978,0.129090421375115,-0.771858392710290) q[0];
u3(1.37505207080733,-2.17003839554623,-3.77241065404739) q[3];
u3(2.31426453903243,1.57063999166925,0.132189025682189) q[0];
u3(1.52938686244923,0.0204149825334901,-4.01905520900954) q[3];
cx q[3],q[0];
u1(1.10639899929101) q[0];
u3(-1.41881462434619,0.0,0.0) q[3];
cx q[0],q[3];
u3(-0.635616723947376,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.36023665798102,1.48496825374253,-1.27514615353477) q[0];
u3(2.82936745217164,-4.35565171578047,-0.139236260312144) q[3];
u3(0.620534428304369,-0.886371337206940,1.48093275082358) q[1];
u3(0.771848505791184,0.450647710809231,-1.64768138804669) q[6];
cx q[6],q[1];
u1(1.60309563195241) q[1];
u3(-0.604070627703394,0.0,0.0) q[6];
cx q[1],q[6];
u3(-0.292809323498434,0.0,0.0) q[6];
cx q[6],q[1];
u3(2.56287893439815,-1.01089042107257,3.79159067914729) q[1];
u3(1.40220768939980,-1.53241105254555,-1.37561710213457) q[6];
u3(1.65715735437880,-0.454851601893218,1.43218623526811) q[2];
u3(1.42929109889572,-0.623314745517630,-0.911498866863683) q[4];
cx q[4],q[2];
u1(2.04738879284941) q[2];
u3(-2.50887294220681,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.0139171473805895,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.26408771879512,4.11426013676763,-0.357822948109682) q[2];
u3(1.23811235857330,1.36303612072946,1.18374299929625) q[4];
u3(0.938058918919248,-2.50691545056335,0.703362277595741) q[7];
u3(1.75963479049453,-3.90617495296255,0.760476077471544) q[5];
cx q[5],q[7];
u1(1.48931322156887) q[7];
u3(-3.19603371724260,0.0,0.0) q[5];
cx q[7],q[5];
u3(2.32877999014449,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.79564869248744,2.18146226051254,-0.361523210313063) q[7];
u3(0.397056462416752,2.00596268788608,-1.32462454654796) q[5];
u3(2.07425539951778,1.77032126471315,-0.0177939516963156) q[5];
u3(2.09604577112003,0.568016536308256,-1.89880611018565) q[3];
cx q[3],q[5];
u1(2.09672533345728) q[5];
u3(-2.67903378174957,0.0,0.0) q[3];
cx q[5],q[3];
u3(0.724046866915064,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.68924587639808,-2.35464032170579,1.34666115041636) q[5];
u3(2.48311272700380,1.03816596782154,4.23215020523068) q[3];
u3(2.72873858393316,0.893070709796607,-3.22599602632663) q[6];
u3(1.56624197593434,2.73286796227005,-2.53341327005768) q[7];
cx q[7],q[6];
u1(3.23135771691090) q[6];
u3(-1.49197118645502,0.0,0.0) q[7];
cx q[6],q[7];
u3(2.20590955052382,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.07782060977663,4.72607235422400,-0.969405282776557) q[6];
u3(2.03351605618656,-1.17105847587487,-2.13904449561400) q[7];
u3(2.53336746530137,0.650080542698223,0.417251088523305) q[0];
u3(1.12433325704265,-3.49793007666538,-0.891194706695627) q[2];
cx q[2],q[0];
u1(0.204705706442889) q[0];
u3(-0.604307452919526,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.12565410710267,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.33429688635611,0.438789106236538,-3.94538049026618) q[0];
u3(1.33450947386186,2.01400235551624,-1.31832853558183) q[2];
u3(2.79337274856088,2.47783881512089,-1.19444406565368) q[1];
u3(1.97812480903914,1.69343838244460,-2.56842604103414) q[4];
cx q[4],q[1];
u1(-1.10802815815433) q[1];
u3(0.207379873933610,0.0,0.0) q[4];
cx q[1],q[4];
u3(3.60905369113133,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.39261434912790,0.957594928129638,0.398081401198928) q[1];
u3(2.95105134444502,-4.54808857654441,-1.40996434692750) q[4];
u3(1.36485809585962,3.67227038455096,-1.73900260572576) q[7];
u3(1.12386289152642,1.65876074465805,-2.58266389801844) q[1];
cx q[1],q[7];
u1(2.61263059004748) q[7];
u3(-1.83594335070454,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.933840527455893,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.15951808100424,-1.40065464513901,0.436738812233902) q[7];
u3(2.12911170362688,-3.89849686561443,-0.131379425139988) q[1];
u3(2.55089494651639,-2.90433405398053,2.87854363889290) q[4];
u3(1.10350871709632,-0.713684394236940,2.75257112314834) q[0];
cx q[0],q[4];
u1(1.39064032969647) q[4];
u3(-0.148829185416577,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.38271611462014,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.11686874278014,-0.403103110899499,-2.29526680972741) q[4];
u3(1.92343170629594,-0.390468574829117,-0.281367913418428) q[0];
u3(1.97509654591004,-1.20489880526661,-0.264635503236225) q[6];
u3(1.82817052358487,-3.98571073503757,-0.0374736385083925) q[5];
cx q[5],q[6];
u1(2.28830316225694) q[6];
u3(-1.77642350417770,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.281782906623035,0.0,0.0) q[5];
cx q[5],q[6];
u3(1.92480469684325,-2.96410456488245,-0.100880474565761) q[6];
u3(0.983518648435477,0.419444204037021,-4.12840740849918) q[5];
u3(0.933899154187786,-1.38365300059959,1.69880190723250) q[3];
u3(0.461487417122752,-2.50402596260205,0.549519186509188) q[2];
cx q[2],q[3];
u1(1.38227993922574) q[3];
u3(-3.19931477058486,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.37488180605497,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.04758881952362,-2.21467258739912,3.21237887451362) q[3];
u3(1.55417203108100,-4.62967802135242,0.546817448113284) q[2];
u3(1.74759863703226,-0.429463650414601,1.24320510354889) q[0];
u3(1.77306882990459,-1.68496849781930,-2.06843150943854) q[6];
cx q[6],q[0];
u1(1.55479832572200) q[0];
u3(-2.81954782743963,0.0,0.0) q[6];
cx q[0],q[6];
u3(0.135603458985984,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.89133583771610,-4.27726426131466,1.41498126823341) q[0];
u3(2.00148892966410,4.76478859513723,0.534603794200336) q[6];
u3(0.941899964607241,-0.988695979448164,-1.46329722151377) q[2];
u3(0.727737187396993,-4.60282165574959,1.19845230544478) q[5];
cx q[5],q[2];
u1(1.05496858034645) q[2];
u3(-3.35000192660813,0.0,0.0) q[5];
cx q[2],q[5];
u3(1.72115403966663,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.08852253728410,2.53524055641668,-1.16865027994501) q[2];
u3(1.14663774644864,-0.131610409190393,-2.95746035800096) q[5];
u3(0.607299401620817,-1.22820870854022,-0.737150214401039) q[1];
u3(2.16053364542502,-4.50533426098443,-0.218075990418582) q[4];
cx q[4],q[1];
u1(0.938323730811224) q[1];
u3(-1.47053789656939,0.0,0.0) q[4];
cx q[1],q[4];
u3(2.41972595368508,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.34032732588421,1.39374779217531,2.15636076173283) q[1];
u3(1.59259708964136,5.62815967098880,-0.143357348919306) q[4];
u3(0.595521310016939,-1.18856221753420,2.24183341918957) q[3];
u3(0.574734491918215,0.544407176127995,-2.35353249577425) q[7];
cx q[7],q[3];
u1(2.51247194288703) q[3];
u3(-1.64767512614592,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.526229179714778,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.23608711308567,-3.96896123485046,1.26227215264594) q[3];
u3(1.01089899829704,1.63104223120210,-2.37364828383642) q[7];
u3(1.00580462545830,2.32476770825419,-2.25425314141064) q[1];
u3(1.57483070291277,-2.90552456832628,2.85454197091534) q[0];
cx q[0],q[1];
u1(0.909299240587443) q[1];
u3(-3.06710472398783,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.98456116431538,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.69748865532596,-3.31887254824717,-0.247656877625528) q[1];
u3(0.909022869694996,-0.194975697991165,-1.29470691253484) q[0];
u3(1.53384961451087,0.453923431822398,2.54126082915105) q[5];
u3(1.66708284090568,-1.53172462230045,-1.78121485644284) q[2];
cx q[2],q[5];
u1(1.07690272200891) q[5];
u3(-3.25797889173673,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.80867738885559,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.60196089856738,1.64857884279551,-1.11266721574949) q[5];
u3(1.71181170927444,-2.13425081218923,-3.24773719180844) q[2];
u3(1.55985010276312,1.78784957265817,-3.79484979189138) q[6];
u3(2.63990589051477,3.34115731966710,-2.65314272199886) q[3];
cx q[3],q[6];
u1(3.68749210459411) q[6];
u3(-3.18638624725113,0.0,0.0) q[3];
cx q[6],q[3];
u3(-0.956928407599300,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.71954819576825,1.37795458264857,0.135068157433538) q[6];
u3(0.630482738190273,-3.74419304331956,-2.17492269820524) q[3];
u3(0.871822616762041,1.71681003629749,-0.898494055016404) q[4];
u3(1.15354760311866,0.0320434221209622,-3.46943885605427) q[7];
cx q[7],q[4];
u1(3.60209167328533) q[4];
u3(-0.951798404156332,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.81243117965460,0.0,0.0) q[7];
cx q[7],q[4];
u3(3.05064189616590,-3.56875622925064,1.32369726271166) q[4];
u3(1.13294343615421,-2.83917094864992,-3.03072448154397) q[7];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];