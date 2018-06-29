OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(2.49645259808434,-0.528317388413675,-0.403164262458361) q[0];
u3(0.799921815731443,-5.28771578689449,0.969159464275232) q[4];
cx q[4],q[0];
u1(3.55636772809183) q[0];
u3(-3.84062640845521,0.0,0.0) q[4];
cx q[0],q[4];
u3(-1.07863467920883,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.87924218788787,1.51210742093185,1.85866036186639) q[0];
u3(0.576142140262267,3.03027719010220,0.398284627315435) q[4];
u3(1.31825094281565,0.455206475195347,0.527929649345268) q[1];
u3(1.19828610994257,-1.49682561810768,-1.20748145605314) q[5];
cx q[5],q[1];
u1(-0.297107364433347) q[1];
u3(-1.66477449429159,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.535107671291670,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.99567140519700,-1.06792686201440,0.433173248706743) q[1];
u3(1.76165897166958,0.658812976039733,2.98676600032394) q[5];
u3(1.62160907851616,1.23882705707599,0.790749194447457) q[2];
u3(1.49056895446658,-0.831816475768218,-2.61484861673135) q[3];
cx q[3],q[2];
u1(0.899532552080878) q[2];
u3(-3.39125432559650,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.69928369785971,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.681735409777325,-1.38856395872311,1.13318659355765) q[2];
u3(1.79898243599759,-1.28254389230871,0.929339577769274) q[3];
u3(0.882772309027609,1.23027117147019,-0.954757744558374) q[1];
u3(0.840667935032030,-1.17952723256977,-0.231475106010921) q[3];
cx q[3],q[1];
u1(2.96349016126856) q[1];
u3(-4.50292187255488,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.0331452177697868,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.28761404111486,1.22789513820098,-1.16475316300369) q[1];
u3(1.90089587096833,2.64484665407707,3.48086908677342) q[3];
u3(0.574625666905654,-0.795302470608364,1.12413945195337) q[2];
u3(0.428869811034214,-1.21541662872536,-0.241020812926736) q[0];
cx q[0],q[2];
u1(0.0258942208502002) q[2];
u3(-1.87053178137102,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.832923091464920,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.10086093371881,0.0385223488992216,-2.01999177084941) q[2];
u3(1.80706806870021,-2.96097333383866,1.26572393784898) q[0];
u3(2.13113634673334,2.51000212873829,-0.616356464478823) q[4];
u3(2.62886204293034,2.26575250859315,-1.73324449651484) q[5];
cx q[5],q[4];
u1(-1.10411577833256) q[4];
u3(0.411592857333734,0.0,0.0) q[5];
cx q[4],q[5];
u3(4.00640106730418,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.776613812742481,0.559219613511806,-0.349697045644848) q[4];
u3(1.30130012497004,-2.43005564020307,-2.63772435948334) q[5];
u3(1.19769416445568,-1.81412098859888,-0.350975244924607) q[4];
u3(1.56672567703365,-1.90213276402016,-0.378322652107466) q[2];
cx q[2],q[4];
u1(0.425502179664968) q[4];
u3(-0.883022229743947,0.0,0.0) q[2];
cx q[4],q[2];
u3(1.35619447949266,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.00343206671813,-1.53291259162420,2.14240234794383) q[4];
u3(1.14504477956129,5.69381847030006,-0.0837852711173728) q[2];
u3(0.470827737689088,1.82921777718979,-2.42252107599839) q[5];
u3(0.239381408198224,-2.71846066315168,2.02993585431242) q[3];
cx q[3],q[5];
u1(3.17733348445960) q[5];
u3(-2.55133448098100,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.80006599423122,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.09688858947104,-1.26853466170870,-2.47243130022077) q[5];
u3(1.03221438444895,1.27695834414658,4.03202199103051) q[3];
u3(2.31414700718305,-0.741984185267597,-0.972489390351300) q[0];
u3(1.00028939141020,-1.59691271716130,-3.25109181587766) q[1];
cx q[1],q[0];
u1(2.09367555731123) q[0];
u3(0.385145038742565,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.66191793881032,0.0,0.0) q[1];
cx q[1],q[0];
u3(0.407227083452494,1.26055325073245,-2.98162413598167) q[0];
u3(1.07028305317273,2.31949314274382,0.845106255441191) q[1];
u3(1.54442013123406,1.91331913951280,-3.16373554640920) q[1];
u3(1.35225710529766,-2.97373541638837,2.90260116327456) q[4];
cx q[4],q[1];
u1(2.28720003488167) q[1];
u3(-1.39163631469335,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.372583738641241,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.63027341267955,-1.35297104978865,3.05491570878498) q[1];
u3(2.25618544579814,1.72289139898806,1.12627051891679) q[4];
u3(1.96985120951709,-2.20816615278383,0.714939770086658) q[2];
u3(0.948365961146132,-3.80499680584780,0.163871032300662) q[3];
cx q[3],q[2];
u1(-0.0302810897594663) q[2];
u3(-2.13057407407520,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.810144965992286,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.52639267692293,-4.44722893974543,1.42005470199122) q[2];
u3(0.529953754275196,0.498340923771195,4.24732602142251) q[3];
u3(1.77988926459199,1.47526401781670,-0.469555332209198) q[0];
u3(0.688178869340503,-0.560615734444704,-2.06815706367224) q[5];
cx q[5],q[0];
u1(1.15363526928923) q[0];
u3(-0.610243055653035,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.10659952162370,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.71617347854119,2.00656778962690,-4.18134278205543) q[0];
u3(0.739437683289678,2.03014106123935,-3.66854978003635) q[5];
u3(2.73570968153100,1.99594978673018,-0.303905947900103) q[2];
u3(2.64965591233576,0.743589041167607,-3.67254215270025) q[1];
cx q[1],q[2];
u1(1.83675968388452) q[2];
u3(-2.22350478750285,0.0,0.0) q[1];
cx q[2],q[1];
u3(-0.202897186928517,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.34674386161244,2.45650184537479,-0.583055464360142) q[2];
u3(1.39398096384310,-0.477556938499294,1.07271469403050) q[1];
u3(1.62845421473416,1.02870322100928,-3.97366646420360) q[5];
u3(1.86801726836634,-1.37810553530712,4.17226795620008) q[4];
cx q[4],q[5];
u1(1.40029321120565) q[5];
u3(-0.0863637481309365,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.81365317187244,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.44377125360333,1.12727472999081,-1.01179634286284) q[5];
u3(1.70339949514532,1.13377223446127,-1.71234111118522) q[4];
u3(0.152429973283758,-1.06157956328143,0.886731748136387) q[3];
u3(0.250572843384918,-1.74899299466792,-0.918526535540715) q[0];
cx q[0],q[3];
u1(2.34444360528567) q[3];
u3(-3.29152085769203,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.923729832121750,0.0,0.0) q[0];
cx q[0],q[3];
u3(0.993934313592129,2.57705768360363,0.201407583479615) q[3];
u3(2.30731853503797,5.32917454676272,0.269297125875435) q[0];
u3(1.26958926581660,1.46606517907300,-3.34459274641423) q[1];
u3(0.922831392548194,-2.34828727649964,2.62552576508072) q[5];
cx q[5],q[1];
u1(-0.973089162285136) q[1];
u3(-0.366208318778521,0.0,0.0) q[5];
cx q[1],q[5];
u3(3.18374174831669,0.0,0.0) q[5];
cx q[5],q[1];
u3(0.520327171162276,-1.56111408567503,0.187944212280836) q[1];
u3(1.03638948311474,-2.58121569564994,0.207883253362076) q[5];
u3(1.65239473393117,0.167248355699331,1.39989539483200) q[2];
u3(1.95490210577269,-1.75050448158052,-0.516789411106572) q[3];
cx q[3],q[2];
u1(1.71603499244463) q[2];
u3(-0.931970920901068,0.0,0.0) q[3];
cx q[2],q[3];
u3(-0.739587247204510,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.964008821667500,-1.72329696746938,0.522198637208463) q[2];
u3(2.74057471914374,-0.735009148457717,2.02731107648074) q[3];
u3(0.749254218516797,1.75609774737077,-3.45382034578379) q[4];
u3(1.66831232882711,-2.31964000311305,3.43150556663901) q[0];
cx q[0],q[4];
u1(0.468135011882677) q[4];
u3(-0.796407068319657,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.92165212146864,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.657787519096605,-3.40595117142519,-0.245349418775734) q[4];
u3(1.34350568855018,0.669511842804599,5.59355612782026) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
