OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
creg c[11];
u3(1.34464341115587,1.78323393904218,-1.36117341829324) q[2];
u3(0.282487836544945,1.68115448912328,-3.56541006164106) q[10];
cx q[10],q[2];
u1(1.69000002208993) q[2];
u3(-1.97757678204616,0.0,0.0) q[10];
cx q[2],q[10];
u3(0.525660017903217,0.0,0.0) q[10];
cx q[10],q[2];
u3(2.03240961272428,-0.513695478591670,3.02960040535219) q[2];
u3(1.78856854494987,4.00375659252574,-1.83420629961496) q[10];
u3(0.691559080741435,-1.79255837189903,2.30514593349884) q[8];
u3(1.18286732199129,1.58284578895010,-3.20004755332816) q[0];
cx q[0],q[8];
u1(0.910291519164182) q[8];
u3(-1.23518366062885,0.0,0.0) q[0];
cx q[8],q[0];
u3(2.39891714692927,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.153364460780913,-2.70439314900102,1.02136012249317) q[8];
u3(0.634291596344652,-1.30156263067468,-0.568790563099226) q[0];
u3(1.52070235467814,1.83864808805069,-0.114536185712345) q[6];
u3(2.33524300828037,0.416407368067205,-2.39321208935678) q[9];
cx q[9],q[6];
u1(0.192719744130203) q[6];
u3(-0.861897016736503,0.0,0.0) q[9];
cx q[6],q[9];
u3(1.76182445354829,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.659842220946075,0.735691582541923,2.15849781114296) q[6];
u3(1.86397045921926,3.83431862632936,-0.994072124724897) q[9];
u3(2.80966329791593,-3.02749959749065,0.0730738501602666) q[1];
u3(2.20386644375388,-1.96267238732251,-1.22287423267143) q[3];
cx q[3],q[1];
u1(0.203420738354114) q[1];
u3(-1.14768412972604,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.75210466446858,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.25844830755139,4.25913442914745,-0.966821159505005) q[1];
u3(2.24565204051299,3.81649790721736,-0.0935350744099057) q[3];
u3(1.39427992050259,-0.965345357395720,0.884629904352174) q[4];
u3(2.31216943111306,-1.99964109334222,-2.17113646170078) q[7];
cx q[7],q[4];
u1(-0.228281368342751) q[4];
u3(-1.76372959966654,0.0,0.0) q[7];
cx q[4],q[7];
u3(0.819047066252381,0.0,0.0) q[7];
cx q[7],q[4];
u3(0.743239988537156,-1.71926657004465,2.45273545835904) q[4];
u3(1.90089413866480,3.40310992535154,0.946197425933989) q[7];
u3(0.239325273099953,-1.23391023634671,2.25256938605974) q[7];
u3(0.831789726130070,-3.06962245324263,1.63907678284894) q[8];
cx q[8],q[7];
u1(0.734141152096826) q[7];
u3(-0.405489916487856,0.0,0.0) q[8];
cx q[7],q[8];
u3(3.06138801321304,0.0,0.0) q[8];
cx q[8],q[7];
u3(1.13043934104335,0.465284897140143,-1.43038984016738) q[7];
u3(1.72538841376440,0.703055861640534,-1.91699373237334) q[8];
u3(2.21175938430920,1.14894190329793,1.75097359592917) q[10];
u3(0.419923083655713,-0.947519748557764,-3.34460728721244) q[0];
cx q[0],q[10];
u1(0.0568255335992605) q[10];
u3(-0.848794842517157,0.0,0.0) q[0];
cx q[10],q[0];
u3(1.93321167941531,0.0,0.0) q[0];
cx q[0],q[10];
u3(1.74307651479103,-2.18227252673261,3.61265121012275) q[10];
u3(1.43104087599716,0.797579235534316,2.63613490551186) q[0];
u3(2.80776956630938,0.531156388894000,-3.39562857150732) q[6];
u3(1.94175835196642,2.64387346620518,-2.08178152143163) q[2];
cx q[2],q[6];
u1(-0.759378742450513) q[6];
u3(0.401371513161220,0.0,0.0) q[2];
cx q[6],q[2];
u3(4.10252990528640,0.0,0.0) q[2];
cx q[2],q[6];
u3(1.59516868195606,-1.44144274285767,2.42384554775043) q[6];
u3(2.42524155175742,-0.617582314034305,4.57558602409810) q[2];
u3(1.48957122330973,-2.25894279348113,0.973322565046127) q[3];
u3(1.51457879928431,-3.24943080969206,-0.00375823951168730) q[9];
cx q[9],q[3];
u1(2.24020478619678) q[3];
u3(-2.99163485723831,0.0,0.0) q[9];
cx q[3],q[9];
u3(1.17717679079252,0.0,0.0) q[9];
cx q[9],q[3];
u3(1.70742525208346,-2.54419509663414,-1.06665896029084) q[3];
u3(1.52314411025252,-0.857573620068440,-3.24165430237818) q[9];
u3(2.13042929302722,2.60306578764196,-3.33437234447395) q[4];
u3(0.392163309425672,-1.62928276334513,2.54233035240360) q[1];
cx q[1],q[4];
u1(4.30265537930175) q[4];
u3(-3.74813909484819,0.0,0.0) q[1];
cx q[4],q[1];
u3(-0.256320387457671,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.57998857835142,1.83600233734778,-0.418516714414680) q[4];
u3(0.979203072820909,-2.79684295573324,-2.74221336378958) q[1];
u3(1.61504510680124,-0.216542397950770,0.348673705546186) q[7];
u3(1.87806135700954,-0.373843672271971,-1.67357835922022) q[0];
cx q[0],q[7];
u1(0.525544754771839) q[7];
u3(-1.25596101558632,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.58123660582480,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.52282549666021,3.06438901474864,-1.07182878306714) q[7];
u3(2.14963579890711,0.978191507312043,-1.51857370825377) q[0];
u3(1.81899921989727,2.48415416992193,-0.631652733062135) q[6];
u3(2.08813406101281,2.16968342436075,-1.08090276743684) q[10];
cx q[10],q[6];
u1(-0.0331772302816316) q[6];
u3(-1.78163175940845,0.0,0.0) q[10];
cx q[6],q[10];
u3(0.932131924628652,0.0,0.0) q[10];
cx q[10],q[6];
u3(0.596965006352620,0.793448988249732,-3.50528805288873) q[6];
u3(1.26649264747415,2.20745767340506,-1.73939418637738) q[10];
u3(0.441081477370774,2.93403572755879,-0.910881981562117) q[8];
u3(1.63261090866113,2.35154441114339,-1.94884362139840) q[2];
cx q[2],q[8];
u1(3.40706955059924) q[8];
u3(-0.803207657697751,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.05477740602744,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.23373802982547,-2.03871773522217,0.466686283800285) q[8];
u3(0.713923347413878,1.01318078172593,-2.91014621508356) q[2];
u3(0.833162183204250,1.71675188397983,-0.658674502545804) q[5];
u3(1.81880669429765,0.606031841853499,-3.53363498256902) q[3];
cx q[3],q[5];
u1(3.44152349705652) q[5];
u3(-1.29290949914074,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.27304733280819,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.29278829142620,-0.504510067941023,-1.38372820567083) q[5];
u3(2.91381782195901,4.26033927201379,0.937886076348223) q[3];
u3(2.18289169640132,-1.30513224716396,2.60169994647576) q[4];
u3(2.89119756213740,-3.56781429033799,-2.07051278415952) q[9];
cx q[9],q[4];
u1(1.42736205632383) q[4];
u3(-0.537359850247921,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.42890393188074,0.0,0.0) q[9];
cx q[9],q[4];
u3(0.227095166666314,2.08171350988075,-1.75898172220084) q[4];
u3(1.24778672798361,0.383747451647850,-5.29770181638987) q[9];
u3(1.17001528972468,1.12476405903290,-3.46925929600617) q[0];
u3(0.727589946994657,2.49954210746963,-2.49111792047187) q[5];
cx q[5],q[0];
u1(2.87641365312121) q[0];
u3(-2.19319858225700,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.75894870064508,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.59878208969776,-0.945298460928294,-2.60020043481888) q[0];
u3(2.60363860545459,2.69234980910308,-1.87044099876928) q[5];
u3(1.53044389722246,1.64235341107348,-1.76284521222319) q[6];
u3(0.228004700428406,-1.38851005135498,0.911157957125675) q[1];
cx q[1],q[6];
u1(0.0122104674439336) q[6];
u3(-1.80299051176202,0.0,0.0) q[1];
cx q[6],q[1];
u3(0.707318285979836,0.0,0.0) q[1];
cx q[1],q[6];
u3(2.26331732478611,-0.437043714870555,3.95552130107464) q[6];
u3(1.71265265826513,-1.83483147722785,-3.72679619229787) q[1];
u3(2.12193659973080,-0.892888877675679,-0.885007972570592) q[2];
u3(1.84830056839478,-3.23116091763075,-0.410771890485927) q[7];
cx q[7],q[2];
u1(-0.145665476351154) q[2];
u3(-1.65658614293280,0.0,0.0) q[7];
cx q[2],q[7];
u3(0.794703566768614,0.0,0.0) q[7];
cx q[7],q[2];
u3(2.77539540005491,-1.70344632836720,3.71022255319312) q[2];
u3(2.65539786222326,1.04992275260414,1.93557128947868) q[7];
u3(1.66532716462705,-0.207363056416202,-1.34450014273106) q[8];
u3(1.96707879857543,-5.28179506893544,0.848903695972777) q[3];
cx q[3],q[8];
u1(0.809386683933932) q[8];
u3(-1.25813327616384,0.0,0.0) q[3];
cx q[8],q[3];
u3(-0.412376896339225,0.0,0.0) q[3];
cx q[3],q[8];
u3(1.36017593726798,-0.792292808516188,1.03595908883022) q[8];
u3(2.07042574029040,0.446438768693827,2.04681934621689) q[3];
u3(1.79458248737844,0.840093828938845,1.47857973135433) q[9];
u3(1.45517376715306,-1.54153661800745,-1.98523375501738) q[4];
cx q[4],q[9];
u1(1.72840996970465) q[9];
u3(-0.158731668440149,0.0,0.0) q[4];
cx q[9],q[4];
u3(0.882402862666307,0.0,0.0) q[4];
cx q[4],q[9];
u3(1.07851000039902,3.19059249657135,-2.30224926940788) q[9];
u3(0.860794209346005,3.57924084659455,-0.992114959862633) q[4];
u3(2.24082751492463,3.70167453013179,-1.97841517701990) q[0];
u3(0.265177820227100,-1.36385108000717,2.40695549156514) q[1];
cx q[1],q[0];
u1(1.78018035647572) q[0];
u3(-2.53499207486987,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.520784188101367,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.55927408780719,-1.93330959729811,1.06428699517298) q[0];
u3(1.01499195396553,1.90073376107285,2.77138610272380) q[1];
u3(1.84004889915351,1.65731420516982,-2.46743848971684) q[3];
u3(1.08566796226264,1.80549243861699,-3.18113402047722) q[2];
cx q[2],q[3];
u1(1.57029695239666) q[3];
u3(-0.792278517978824,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.0975321836857854,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.74497721913418,0.289948425488075,-0.254200837600914) q[3];
u3(1.16238830078037,3.30005832201326,1.42955561857609) q[2];
u3(0.954105870255198,-0.759620223818515,0.735240038313765) q[5];
u3(0.298125291943746,-1.70423431619782,0.982846972151616) q[4];
cx q[4],q[5];
u1(1.48123719254324) q[5];
u3(-3.03811460766856,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.01048772596442,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.91581815163444,-2.49940647070513,-0.241930404255345) q[5];
u3(2.49808920927674,3.02956428642830,2.13404628992949) q[4];
u3(1.57850132180255,2.66196278706848,0.194697326771921) q[9];
u3(2.09248508914971,0.0466208163735444,-3.46476400513297) q[10];
cx q[10],q[9];
u1(1.75256252874633) q[9];
u3(-2.63501477053572,0.0,0.0) q[10];
cx q[9],q[10];
u3(0.279326349439971,0.0,0.0) q[10];
cx q[10],q[9];
u3(1.71832420894731,0.815632611688910,0.0975435175806454) q[9];
u3(2.23043970553969,3.29399861570625,0.702190005416072) q[10];
u3(1.96307906280260,1.63337871692612,1.35096823107769) q[8];
u3(0.289616342880516,-0.999345962412913,-2.95466429428076) q[6];
cx q[6],q[8];
u1(0.258353485656899) q[8];
u3(-1.25353757634674,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.51924149564507,0.0,0.0) q[6];
cx q[6],q[8];
u3(0.911581001162394,3.46311049634402,-2.61602783041445) q[8];
u3(1.77909386457518,2.55754563644045,-0.689382957718865) q[6];
u3(1.16576334008042,-0.652781021682083,1.06990671951939) q[9];
u3(1.35083348920717,-1.74985075361026,-2.09632045562408) q[10];
cx q[10],q[9];
u1(1.20848331544989) q[9];
u3(-1.00252041794695,0.0,0.0) q[10];
cx q[9],q[10];
u3(3.64354867306552,0.0,0.0) q[10];
cx q[10],q[9];
u3(2.57532440937863,-2.97883810733502,0.0292763017507471) q[9];
u3(2.03804731774215,-3.87708772193772,0.601803516357843) q[10];
u3(1.86968609965926,0.471888766904415,-3.49125991202799) q[4];
u3(0.919621077911687,-3.08616499304569,2.33872579717239) q[5];
cx q[5],q[4];
u1(2.72461404411200) q[4];
u3(-3.04200435020725,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.11683400243792,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.77691724214332,-1.42004178800585,1.09788896859641) q[4];
u3(2.11184767785022,0.972235381321526,-2.09490140685046) q[5];
u3(1.19184826212020,-1.27030410352473,1.19287954748270) q[2];
u3(1.22229276072583,-1.47667689403250,-1.75129067943997) q[3];
cx q[3],q[2];
u1(2.35014169134998) q[2];
u3(-2.87844232885843,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.32358771085169,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.56821235681213,0.210001405035529,1.22431001201239) q[2];
u3(1.41719881276064,0.872229876689157,3.43053507949582) q[3];
u3(1.48438740190844,-1.25770574253943,2.67479658061431) q[1];
u3(1.19312628129718,-1.47430455988842,-1.58758995547021) q[0];
cx q[0],q[1];
u1(2.01431530547063) q[1];
u3(-3.16013115844264,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.55687060916357,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.43806481905929,-1.71037271256744,3.80611232129134) q[1];
u3(2.61213472727582,-1.02054272820227,2.95478968175264) q[0];
u3(1.40478136803958,3.39777405713379,-1.80419286493018) q[7];
u3(1.03390526820289,2.90192974967348,-2.72917872215753) q[8];
cx q[8],q[7];
u1(2.14710657358495) q[7];
u3(-1.68203495062528,0.0,0.0) q[8];
cx q[7],q[8];
u3(0.767855585214369,0.0,0.0) q[8];
cx q[8],q[7];
u3(0.783112977529975,-2.22460908761535,-0.746348149706710) q[7];
u3(2.29561655288374,-3.82216086637622,2.23711152433102) q[8];
u3(0.345732785293022,2.19697181191146,-2.77618967703789) q[0];
u3(0.762708510024015,-3.22157910844677,2.51383353977674) q[5];
cx q[5],q[0];
u1(0.712554066383621) q[0];
u3(-1.78298915743060,0.0,0.0) q[5];
cx q[0],q[5];
u3(2.77956123705056,0.0,0.0) q[5];
cx q[5],q[0];
u3(2.25394926256196,2.39059253867169,0.0651890069696031) q[0];
u3(0.799014707637653,-0.581160512286192,-3.57905931414350) q[5];
u3(0.966170673057622,1.73826868936476,-3.08412072295301) q[8];
u3(0.839037086368265,-2.57584794184230,2.71992373419042) q[2];
cx q[2],q[8];
u1(3.27440797935650) q[8];
u3(-1.79633724796842,0.0,0.0) q[2];
cx q[8],q[2];
u3(2.44748748354530,0.0,0.0) q[2];
cx q[2],q[8];
u3(2.68304137268380,2.23140227459734,-1.19702867692061) q[8];
u3(2.52629044951055,1.58150909963548,-1.66209418925186) q[2];
u3(0.649840820258747,3.14061567465431,-0.843769808323449) q[1];
u3(1.57366604143753,1.54492560105936,-0.765460815773472) q[6];
cx q[6],q[1];
u1(0.362767399691863) q[1];
u3(-1.14558835364318,0.0,0.0) q[6];
cx q[1],q[6];
u3(1.66881504527679,0.0,0.0) q[6];
cx q[6],q[1];
u3(1.80041333810788,0.0997518426464142,0.951966397931657) q[1];
u3(1.12213287024087,-4.37549313345168,1.15384159592342) q[6];
u3(1.60566038293074,-0.353802973950857,0.0232737319688635) q[10];
u3(1.16403406603844,-3.45220463362375,-0.0431759558070652) q[4];
cx q[4],q[10];
u1(1.66925021292606) q[10];
u3(-3.07811362024560,0.0,0.0) q[4];
cx q[10],q[4];
u3(0.791830187445924,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.86583976974435,2.60184430396873,-2.29903143167363) q[10];
u3(1.20971873353374,4.92607087670033,-0.0963210454806318) q[4];
u3(2.04498037368556,1.74562379434660,-2.33779690197661) q[7];
u3(0.469976668721719,-2.54534586923808,1.88287352891775) q[3];
cx q[3],q[7];
u1(1.26529197482053) q[7];
u3(-0.700733994008744,0.0,0.0) q[3];
cx q[7],q[3];
u3(2.57793652059316,0.0,0.0) q[3];
cx q[3],q[7];
u3(1.59798965573591,0.991595872269813,1.37851157926106) q[7];
u3(2.18679437814872,-4.06330300590279,2.03102149233513) q[3];
u3(1.82259430222461,-4.58829452578495,1.49870244204854) q[8];
u3(0.453125953607684,1.93635724905023,-0.777500338821337) q[6];
cx q[6],q[8];
u1(2.17170083207024) q[8];
u3(-2.64635331508801,0.0,0.0) q[6];
cx q[8],q[6];
u3(1.25775665454654,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.31905364144535,-0.111610278584095,0.489771463880327) q[8];
u3(0.829513972180101,1.98944291291553,2.96417275319808) q[6];
u3(1.86090469697949,2.22102386658529,-2.25570695305871) q[4];
u3(1.55148345755366,2.92549724361962,-3.35613780232480) q[9];
cx q[9],q[4];
u1(1.20192564280549) q[4];
u3(-0.139353500311427,0.0,0.0) q[9];
cx q[4],q[9];
u3(2.63344191709114,0.0,0.0) q[9];
cx q[9],q[4];
u3(1.66242372418053,-0.164681945681266,2.06022191439296) q[4];
u3(2.56787176152718,-1.16705072942591,-4.77287554985199) q[9];
u3(0.990661500306805,3.30815754981952,-1.64406555282340) q[10];
u3(1.46583967597550,2.92320217193686,-2.36570228624889) q[1];
cx q[1],q[10];
u1(1.89537584447884) q[10];
u3(-2.75621616308429,0.0,0.0) q[1];
cx q[10],q[1];
u3(0.607233690305980,0.0,0.0) q[1];
cx q[1],q[10];
u3(2.24785963388310,-4.73474719584568,1.25352713653436) q[10];
u3(0.688068718238371,3.12497546766515,1.59138312748735) q[1];
u3(2.65504290604094,-0.855417708584370,-1.50965850627103) q[5];
u3(1.26482344140341,-5.67622239816743,0.178284246818219) q[0];
cx q[0],q[5];
u1(1.32380614108783) q[5];
u3(-0.210817266896221,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.26272030752414,0.0,0.0) q[0];
cx q[0],q[5];
u3(1.35870783426387,-1.09139988390313,-2.37863531355666) q[5];
u3(0.730940380364103,-1.45269730272794,-2.92064643535871) q[0];
u3(2.68623493250132,0.959865792240535,-2.75272590391349) q[2];
u3(2.70275803703243,4.07468232517730,-1.53509552675629) q[3];
cx q[3],q[2];
u1(0.475215270805551) q[2];
u3(-1.70240381983917,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.08463457776171,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.453263588220247,1.38800560096545,-3.45343183955911) q[2];
u3(1.24105976581979,4.95167685745066,-1.30052453947372) q[3];
u3(1.13737749401503,1.86657928567013,-1.31539854578601) q[10];
u3(0.157936600597474,-0.441298676201119,-0.412396539496044) q[2];
cx q[2],q[10];
u1(0.736729151769310) q[10];
u3(-1.35017391355646,0.0,0.0) q[2];
cx q[10],q[2];
u3(3.22588850306721,0.0,0.0) q[2];
cx q[2],q[10];
u3(2.68409805903045,1.67638778462300,-3.78665915571965) q[10];
u3(1.47395262445033,-4.72202466200850,-0.436060501001654) q[2];
u3(0.974946228641038,3.12444581744450,-1.99123770442466) q[4];
u3(1.41817101753823,0.971902473334948,-2.34881906772444) q[0];
cx q[0],q[4];
u1(1.44580992811483) q[4];
u3(-2.93873648240548,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.17481254929324,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.50154033402919,-2.45008550377831,2.78000792671581) q[4];
u3(1.66613779895113,-4.36523088640010,-1.57878355040946) q[0];
u3(0.681917083300393,2.02346542787437,-0.804194728420721) q[7];
u3(2.15277311974529,-0.0437639734764965,-3.30155608222683) q[6];
cx q[6],q[7];
u1(1.69589827005083) q[7];
u3(-2.66592503065822,0.0,0.0) q[6];
cx q[7],q[6];
u3(0.785682740782307,0.0,0.0) q[6];
cx q[6],q[7];
u3(2.78801586523299,0.733721828657836,-2.28896535762872) q[7];
u3(0.970115233756988,-0.761754850773739,0.443177335332902) q[6];
u3(1.04173345520823,2.53027229226776,-1.32580254469301) q[5];
u3(0.668832996614878,1.08481824718257,-2.30713441204340) q[8];
cx q[8],q[5];
u1(2.02451233150099) q[5];
u3(0.475085631518814,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.40088290161068,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.55876435121287,1.25433445276725,-1.77872361227357) q[5];
u3(2.16460665580139,-4.95861514538616,-0.437278231052405) q[8];
u3(0.558165889301953,0.690445294126014,-2.03699831549418) q[1];
u3(1.44808128745874,-2.42348138843519,3.16024875513046) q[3];
cx q[3],q[1];
u1(0.881279142846002) q[1];
u3(0.0307853461223295,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.49626882706344,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.411213681148505,2.69132613902409,-0.356636475194388) q[1];
u3(2.24835198849795,3.94278746815697,0.697227248694363) q[3];
u3(1.94105361006742,2.88996081595119,-2.06715344865743) q[7];
u3(2.61957501969721,1.94901066734224,-0.123839616335645) q[8];
cx q[8],q[7];
u1(1.28210423019367) q[7];
u3(-0.414458654455028,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.83170882770250,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.37819838870704,1.91114450185135,-1.92275076463881) q[7];
u3(2.11055002103636,0.411615185958804,-3.51324051215546) q[8];
u3(1.39036235220233,2.89674632879349,-1.35451594497424) q[5];
u3(1.56981024821080,1.07669742987424,-2.79830780878158) q[6];
cx q[6],q[5];
u1(2.24078094072985) q[5];
u3(-2.88998667327038,0.0,0.0) q[6];
cx q[5],q[6];
u3(1.16595815067217,0.0,0.0) q[6];
cx q[6],q[5];
u3(1.84632601297798,1.66133698899652,-3.04112855427896) q[5];
u3(2.20084116716410,-1.98171896343372,4.22013736609327) q[6];
u3(1.49213435247075,1.18399822008833,-1.14780189064662) q[0];
u3(1.63574201604782,1.75041824165671,-4.00028959127886) q[3];
cx q[3],q[0];
u1(0.447057045364694) q[0];
u3(-1.69510359387251,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.51553195174167,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.56684607848989,0.683058392663686,-0.570859885681807) q[0];
u3(1.79549670705437,4.11538958840654,-1.07644814557277) q[3];
u3(2.28154545333978,-0.919754836607187,-0.621748825236916) q[2];
u3(0.663601209937765,-0.550331363177524,-4.70442960471682) q[9];
cx q[9],q[2];
u1(0.597665385158607) q[2];
u3(-3.32494698349042,0.0,0.0) q[9];
cx q[2],q[9];
u3(1.74628828594029,0.0,0.0) q[9];
cx q[9],q[2];
u3(1.36924493609483,-1.09669345196139,-1.91373688366937) q[2];
u3(0.737827330479512,-1.69828564708532,3.88360222704401) q[9];
u3(1.75798073626846,-2.49600632019331,0.181306238116065) q[10];
u3(2.59786171955614,-2.58605018973372,-1.20067892217897) q[4];
cx q[4],q[10];
u1(2.58547601844588) q[10];
u3(-2.86228309641441,0.0,0.0) q[4];
cx q[10],q[4];
u3(1.25161687931248,0.0,0.0) q[4];
cx q[4],q[10];
u3(2.43626198729999,-0.825620704909857,-0.718375409812906) q[10];
u3(1.89536510429219,3.08303123809694,0.447003108638669) q[4];
u3(2.80962563328690,1.72137016154783,-3.37889176604833) q[8];
u3(1.03718892025793,-1.68220080941623,2.87345561428400) q[5];
cx q[5],q[8];
u1(0.607594561712173) q[8];
u3(-1.23984506063146,0.0,0.0) q[5];
cx q[8],q[5];
u3(2.38758576824449,0.0,0.0) q[5];
cx q[5],q[8];
u3(1.68547540601048,0.714998266350755,-0.671027300253454) q[8];
u3(1.39834099923025,-4.44910232620261,-0.558832438631471) q[5];
u3(2.15855084444482,-1.05458930599961,-0.844614204612762) q[6];
u3(0.526098820008088,-2.50673594923038,-1.07387745691793) q[9];
cx q[9],q[6];
u1(3.53465327893593) q[6];
u3(-4.12410964708561,0.0,0.0) q[9];
cx q[6],q[9];
u3(-0.847119880237590,0.0,0.0) q[9];
cx q[9],q[6];
u3(2.85516411782203,-0.458302348635867,0.316780622903250) q[6];
u3(1.25444800972359,-0.0142113087726725,-2.59338061709767) q[9];
u3(1.59118100338164,-0.180247343047753,2.18814107720693) q[10];
u3(1.17984362106870,-1.02364662838720,-1.42197684255284) q[4];
cx q[4],q[10];
u1(3.13680805214401) q[10];
u3(-1.56075629663798,0.0,0.0) q[4];
cx q[10],q[4];
u3(2.58795601068174,0.0,0.0) q[4];
cx q[4],q[10];
u3(1.78805153467308,3.16426414246337,-0.689549917885309) q[10];
u3(1.57753343118466,1.72515275240681,-3.32186595438945) q[4];
u3(1.94338137172718,-0.227814214022327,1.96442234326817) q[7];
u3(1.79025853185760,-2.94977576726824,-2.19379129132539) q[3];
cx q[3],q[7];
u1(2.19221565368567) q[7];
u3(0.0700535045436799,0.0,0.0) q[3];
cx q[7],q[3];
u3(1.51619558107143,0.0,0.0) q[3];
cx q[3],q[7];
u3(0.950875023875663,0.861400255058631,-1.73448310591549) q[7];
u3(1.74591020404879,-1.45141725682731,4.68808911555823) q[3];
u3(1.10145591193388,-1.19603090581701,0.0199037383471894) q[1];
u3(1.60015050654212,-2.73642133302709,0.656089067870630) q[0];
cx q[0],q[1];
u1(2.35329313437797) q[1];
u3(-1.59499784939990,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.57715347678534,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.14822324760603,-0.798241566181517,-1.73304990799273) q[1];
u3(1.62344330640410,-1.45242579976936,0.899540757833113) q[0];
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
