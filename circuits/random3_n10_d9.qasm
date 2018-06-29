OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg c[10];
u3(0.191202955943525,2.75769968401348,-2.91666323200197) q[8];
u3(0.467112354552648,1.02588350923956,-1.78023242402434) q[1];
cx q[1],q[8];
u1(2.01675472414360) q[8];
u3(-0.0400563649130525,0.0,0.0) q[1];
cx q[8],q[1];
u3(0.599381025836557,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.41102279607360,-1.87452302913230,-0.862250475379004) q[8];
u3(1.52580163255154,-1.44034328540205,1.53205672901435) q[1];
u3(1.75309782223538,0.928958404256352,-2.77556251247611) q[6];
u3(1.00426616494729,-3.30855072502318,2.43815996176128) q[0];
cx q[0],q[6];
u1(3.49700113599251) q[6];
u3(-1.63915719602596,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.41342185544734,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.31484193255939,2.85955580344339,-0.244706949434546) q[6];
u3(0.686817939927257,-0.886816398562658,5.32044276215044) q[0];
u3(1.08188629429996,1.60021032203604,-0.854757415553761) q[2];
u3(0.404494614774007,-2.35913756500199,0.513490569883615) q[4];
cx q[4],q[2];
u1(1.35131858130387) q[2];
u3(-3.21824132822356,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.79493201824654,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.66508699233666,0.299805586150192,2.25665820035815) q[2];
u3(1.49385199715652,-1.42750276069062,-3.03615233431947) q[4];
u3(1.13597319668211,2.49627708115357,-2.34408636263697) q[5];
u3(1.23315296821340,1.33097259045945,-1.72810959570814) q[7];
cx q[7],q[5];
u1(-0.0592734615461761) q[5];
u3(-1.11042790354544,0.0,0.0) q[7];
cx q[5],q[7];
u3(1.81765424718347,0.0,0.0) q[7];
cx q[7],q[5];
u3(1.18155185738933,0.291351361422102,4.35962801952159) q[5];
u3(0.367018056405224,3.66660383370233,-0.700242919591834) q[7];
u3(2.82215608694090,1.34498878771954,-2.54095304433354) q[9];
u3(1.48820469866390,-3.25021661753093,2.48595908467610) q[3];
cx q[3],q[9];
u1(2.67613480152146) q[9];
u3(-2.45525881198656,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.31097635072439,0.0,0.0) q[3];
cx q[3],q[9];
u3(1.73149133906620,1.95245909252003,-1.62695702398025) q[9];
u3(0.156296932586039,-1.65443722713800,-3.61292861518930) q[3];
u3(2.04489559820197,-1.29525608477838,-0.238730281866362) q[4];
u3(1.43737136719750,-4.04149683664946,-0.945969917125809) q[3];
cx q[3],q[4];
u1(3.64960903259898) q[4];
u3(-4.37016475752912,0.0,0.0) q[3];
cx q[4],q[3];
u3(-0.222991869155932,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.66697091044196,1.39811065253806,-4.35394883991441) q[4];
u3(0.437602308065308,1.10449805762195,5.02797087194936) q[3];
u3(1.39890957460576,-0.939261341157771,0.850133619008232) q[9];
u3(1.73022622817557,-2.75240335230938,-0.123130082810555) q[1];
cx q[1],q[9];
u1(1.40717945678149) q[9];
u3(-3.29268160820976,0.0,0.0) q[1];
cx q[9],q[1];
u3(2.61218473722769,0.0,0.0) q[1];
cx q[1],q[9];
u3(1.78583273279887,0.332227233648429,-0.786784177633619) q[9];
u3(2.92450551909442,-3.48531404887420,-1.64880439518670) q[1];
u3(2.41259269168085,0.621075805316860,1.53099523679363) q[5];
u3(0.765604666334541,-1.86575116586563,-2.25738810539898) q[7];
cx q[7],q[5];
u1(2.01020568787656) q[5];
u3(-3.42334848988032,0.0,0.0) q[7];
cx q[5],q[7];
u3(-0.154875957139225,0.0,0.0) q[7];
cx q[7],q[5];
u3(0.339644523472557,-2.60858042497758,-0.886100582948753) q[5];
u3(1.70582506286811,-1.56702753610630,-0.955947504512503) q[7];
u3(2.56964292323586,0.291102164990546,1.85561678996482) q[8];
u3(1.50994765961679,-3.04360441135790,-2.82826792239929) q[0];
cx q[0],q[8];
u1(2.90596768876312) q[8];
u3(-1.55758201999736,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.93312090577866,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.92439686724480,-1.82526206870223,3.34822148864656) q[8];
u3(0.583247119790877,1.48351176604601,4.51625100461539) q[0];
u3(1.15729176565095,2.66252411759515,-1.59438039891680) q[2];
u3(1.45600532124037,0.806627472355812,-2.66322888686137) q[6];
cx q[6],q[2];
u1(1.45681077507299) q[2];
u3(-0.884096356046446,0.0,0.0) q[6];
cx q[2],q[6];
u3(-0.435822786928392,0.0,0.0) q[6];
cx q[6],q[2];
u3(1.20542162336851,-0.182709182732126,1.64759724424999) q[2];
u3(2.41653658186876,0.457930462441042,-2.71780907262211) q[6];
u3(2.41076489132670,2.16366003787360,0.0139577108281477) q[6];
u3(2.76737163118514,2.56531902812078,-2.07725483436350) q[7];
cx q[7],q[6];
u1(2.69226051891137) q[6];
u3(-1.73359312170044,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.681005261507359,0.0,0.0) q[7];
cx q[7],q[6];
u3(0.707116489162242,3.30021601231991,-1.30756137020281) q[6];
u3(0.928721717538556,1.26757882629838,1.63807657073410) q[7];
u3(2.12750663586388,-3.59982286063437,2.61945678255191) q[2];
u3(2.27220167793126,2.57303280681294,-3.64884255424742) q[3];
cx q[3],q[2];
u1(1.25110726445530) q[2];
u3(-3.41738263144294,0.0,0.0) q[3];
cx q[2],q[3];
u3(2.16050556368131,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.45405984138335,-1.90728614420255,0.340452339921777) q[2];
u3(1.53060338452532,-4.22583589769031,0.0954225182333595) q[3];
u3(0.980486663280962,-1.20667749325618,0.892405094325809) q[0];
u3(0.901084141398893,-1.19748307222538,-1.10102567658170) q[4];
cx q[4],q[0];
u1(2.25956790693276) q[0];
u3(0.169955282523308,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.07395840860268,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.13280444268969,2.63647066803125,-0.607150083784713) q[0];
u3(2.10438219651062,-0.0473586882389813,0.971696897901704) q[4];
u3(2.47492831894272,-3.24214874672451,0.190544225220544) q[8];
u3(3.03036544074863,0.683989386765447,1.64165827660797) q[1];
cx q[1],q[8];
u1(0.531954775465179) q[8];
u3(-1.50507046553022,0.0,0.0) q[1];
cx q[8],q[1];
u3(1.99999897341330,0.0,0.0) q[1];
cx q[1],q[8];
u3(1.16326516106788,-3.39929405158580,2.40461400128315) q[8];
u3(1.66555752699156,-1.94957695213210,-2.41195085988027) q[1];
u3(1.01336057683079,0.358503992601523,-2.12790947467697) q[9];
u3(1.47146673916331,-3.59942752344823,2.15224718427068) q[5];
cx q[5],q[9];
u1(3.63988663772432) q[9];
u3(-4.42160839099068,0.0,0.0) q[5];
cx q[9],q[5];
u3(-0.209177461492553,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.68122439134331,-0.437250014450117,3.59858752101568) q[9];
u3(1.10414825621319,0.340102637149336,1.53449263159577) q[5];
u3(2.18620324318805,-0.805615563741284,1.28834317533988) q[6];
u3(2.02549332206175,-1.64153348571658,-0.773669937221666) q[0];
cx q[0],q[6];
u1(3.12243301332792) q[6];
u3(-1.02273431163505,0.0,0.0) q[0];
cx q[6],q[0];
u3(1.41200307871147,0.0,0.0) q[0];
cx q[0],q[6];
u3(0.472505782659424,-0.292078218655170,-1.93147733110113) q[6];
u3(1.64101525357209,0.705681523471343,0.737881713630072) q[0];
u3(1.13806049143367,0.276441643793454,-2.18941246556175) q[1];
u3(1.70205063260461,2.57600533455385,-3.58966085533061) q[3];
cx q[3],q[1];
u1(-0.528207524859253) q[1];
u3(-2.12485066706190,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.27044263993121,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.756880627075491,1.18102032933843,0.287460937064180) q[1];
u3(1.49159583038404,-0.140976683243823,5.11810399577463) q[3];
u3(0.112184020616503,-0.160456470375766,0.605556385294947) q[4];
u3(0.416418607229590,2.56240624796578,-3.67682004913480) q[7];
cx q[7],q[4];
u1(2.20026712008585) q[4];
u3(-0.206205514200095,0.0,0.0) q[7];
cx q[4],q[7];
u3(1.31806298372030,0.0,0.0) q[7];
cx q[7],q[4];
u3(1.85888951892982,-0.446751587749188,3.05888506541748) q[4];
u3(0.874236914928858,-1.41848059229608,1.08859063108547) q[7];
u3(2.07827154446012,2.21750116999993,-2.40320854351289) q[2];
u3(1.00926929896514,-3.15754333966780,2.74464546139012) q[9];
cx q[9],q[2];
u1(3.06451237575604) q[2];
u3(-2.17384217024073,0.0,0.0) q[9];
cx q[2],q[9];
u3(0.458369152690152,0.0,0.0) q[9];
cx q[9],q[2];
u3(2.16094333977987,1.81180606733935,-0.697584075656914) q[2];
u3(0.957526714356500,-4.60880732366568,-0.500673868429790) q[9];
u3(1.98090285062571,3.11988932804719,-1.81289859791611) q[5];
u3(1.92030443248235,1.92669902867600,-0.648033147567758) q[8];
cx q[8],q[5];
u1(3.59384473876694) q[5];
u3(-0.992297709567166,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.55765423550510,0.0,0.0) q[8];
cx q[8],q[5];
u3(2.87373120209847,-1.01639694557950,1.34859125301099) q[5];
u3(2.15893496432145,2.97846009102299,-2.90339852846800) q[8];
u3(0.645734465944524,2.03178406014579,-2.81971652371371) q[0];
u3(0.891304053109965,0.330340255049022,-2.18171330304228) q[6];
cx q[6],q[0];
u1(3.12348098900753) q[0];
u3(-1.50011058239349,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.33429144897890,0.0,0.0) q[6];
cx q[6],q[0];
u3(2.45901004956863,-2.05367468538206,2.86779809223468) q[0];
u3(0.540563747177090,-0.261841681870911,-0.984006430122198) q[6];
u3(2.34360206221020,-1.36119967086008,2.84796381421854) q[8];
u3(2.85734189309038,2.12077757296163,4.07286351239804) q[3];
cx q[3],q[8];
u1(1.71975322047479) q[8];
u3(-2.93736183205245,0.0,0.0) q[3];
cx q[8],q[3];
u3(0.765397607614302,0.0,0.0) q[3];
cx q[3],q[8];
u3(0.476112149654893,-3.12841526970796,-0.0850835723989964) q[8];
u3(1.98129213723316,2.39994020946794,-0.760491325059826) q[3];
u3(2.27656205216566,1.27465271567157,-2.77951277258492) q[2];
u3(0.946741542734853,-2.29469728827157,1.92767518220559) q[7];
cx q[7],q[2];
u1(3.00699573724405) q[2];
u3(-2.24319945179904,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.55731280058883,0.0,0.0) q[7];
cx q[7],q[2];
u3(1.23048792286858,-0.431945009622562,0.256029737401356) q[2];
u3(2.94021923212917,-3.70810740335426,2.31146791597366) q[7];
u3(1.06943058211937,3.20131742205613,-0.434307937561700) q[5];
u3(1.46842119417597,0.716707823935703,-1.14095246014354) q[9];
cx q[9],q[5];
u1(3.08695432977332) q[5];
u3(-2.27568964729856,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.12482430052914,0.0,0.0) q[9];
cx q[9],q[5];
u3(2.21949102676925,0.722233204636336,-0.657804294588265) q[5];
u3(0.598667693739339,1.42307787424326,-3.74163487446281) q[9];
u3(0.425991753522825,-2.57779505777373,2.04757522086822) q[1];
u3(1.01500992476453,-3.22238407932456,2.86883025891747) q[4];
cx q[4],q[1];
u1(-0.00524047542450501) q[1];
u3(-2.08942953020153,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.846338632714109,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.27696554351689,-0.933234483842750,2.87188193297070) q[1];
u3(1.84364105179214,0.407317938113782,0.209112960881157) q[4];
u3(1.75420476959161,1.78845519016616,-0.194726035113310) q[6];
u3(1.12961913924968,0.464078689090671,-4.33328046997724) q[9];
cx q[9],q[6];
u1(0.737217926378051) q[6];
u3(-0.338446944187587,0.0,0.0) q[9];
cx q[6],q[9];
u3(1.94373124584810,0.0,0.0) q[9];
cx q[9],q[6];
u3(0.730147401398190,-3.06361872036450,0.937600306426939) q[6];
u3(1.42295189739143,-2.70910679960859,1.37690593972997) q[9];
u3(0.574001660600300,1.43634920342619,-1.62275851381184) q[7];
u3(0.912453786442504,-0.955741924636423,-1.50298528084790) q[2];
cx q[2],q[7];
u1(1.09634869078025) q[7];
u3(-0.104394603234547,0.0,0.0) q[2];
cx q[7],q[2];
u3(2.64932909174065,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.26608311150586,0.989702932717559,-0.462099967453312) q[7];
u3(1.06157936189453,-1.17386238764618,1.09179331781886) q[2];
u3(1.68550499166072,1.64956199766191,1.08107535463637) q[1];
u3(2.49795363394344,0.624345455316445,-3.11926105647031) q[3];
cx q[3],q[1];
u1(0.572482705604895) q[1];
u3(-1.30830233651211,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.00596865435864236,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.40129667069531,1.66414224341523,-1.66035015919249) q[1];
u3(0.0848707961572166,-3.07568119949916,2.28609224181883) q[3];
u3(1.40052617999459,-1.86483055197212,-0.668306303114715) q[4];
u3(1.98704319434498,-2.42441751656647,-0.613058652049116) q[0];
cx q[0],q[4];
u1(3.45033698037525) q[4];
u3(-0.622931076818772,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.79925261613225,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.53002309422759,1.21901422783239,0.755446147507128) q[4];
u3(0.589334445428429,-1.44865243684839,3.15031674007142) q[0];
u3(1.88924894425510,-2.18784824551448,-0.526059204559500) q[5];
u3(0.697000131268159,-4.30221052275713,-0.0169094327106905) q[8];
cx q[8],q[5];
u1(1.12885975506878) q[5];
u3(0.0650736166955950,0.0,0.0) q[8];
cx q[5],q[8];
u3(1.73065177230264,0.0,0.0) q[8];
cx q[8],q[5];
u3(1.57701624541955,-1.85712470582656,2.35145812634013) q[5];
u3(2.09065266602765,-3.82767161806249,2.42664401174662) q[8];
u3(1.39150398241642,-3.85292188585715,0.940329310752695) q[4];
u3(2.47658539909075,-0.0837687241546645,2.13026396591068) q[1];
cx q[1],q[4];
u1(0.0351188510020546) q[4];
u3(-1.98647532314261,0.0,0.0) q[1];
cx q[4],q[1];
u3(0.831896962144287,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.573538237708627,-1.64162550970098,0.956684199854931) q[4];
u3(1.97783591386093,2.39274472262344,0.491976998937193) q[1];
u3(2.70892970612971,-1.02927910755988,-1.81931898170408) q[5];
u3(1.17047027725585,0.635560178288178,-5.36995235561769) q[9];
cx q[9],q[5];
u1(4.28500608403137) q[5];
u3(-3.65543131196106,0.0,0.0) q[9];
cx q[5],q[9];
u3(-0.897669947319593,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.25433997869597,-1.48888264491217,4.07654646109683) q[5];
u3(1.23258627995383,-6.12564637242767,-0.000718560231190768) q[9];
u3(1.96926364860025,-1.88018768735333,0.757905533532674) q[3];
u3(1.80702936226626,-3.47886921321809,-0.508066803510245) q[6];
cx q[6],q[3];
u1(0.791486538440751) q[3];
u3(-0.111803959147874,0.0,0.0) q[6];
cx q[3],q[6];
u3(2.69665104832150,0.0,0.0) q[6];
cx q[6],q[3];
u3(0.359633147779792,1.26361874359574,-1.53336747388877) q[3];
u3(2.09586852691194,0.951323291591431,2.55186280042355) q[6];
u3(1.57371008543945,0.809591148128320,-3.95020046064664) q[8];
u3(1.96013446378900,2.93313702368747,-2.43066149824143) q[0];
cx q[0],q[8];
u1(0.943789822522835) q[8];
u3(0.115326309881792,0.0,0.0) q[0];
cx q[8],q[0];
u3(1.91365829158877,0.0,0.0) q[0];
cx q[0],q[8];
u3(1.09843097959338,1.11067145918521,-0.980023142893194) q[8];
u3(1.70612829381246,-1.69837085897873,1.82289766055616) q[0];
u3(0.981763576546990,2.90149005260177,-3.04253248685341) q[2];
u3(1.03209750074618,1.21426319403430,-1.57941428969086) q[7];
cx q[7],q[2];
u1(1.88067373236148) q[2];
u3(-2.53092494117417,0.0,0.0) q[7];
cx q[2],q[7];
u3(3.36420298156137,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.861380184297799,2.92583153553109,-3.21974496048525) q[2];
u3(3.07813185361371,3.76944003557552,0.171450625198875) q[7];
u3(3.01181773761422,-2.43109929504379,0.923144663517059) q[4];
u3(2.84207843616297,-3.43609091511375,-1.94620421699185) q[5];
cx q[5],q[4];
u1(1.11691419012825) q[4];
u3(-0.149377165649917,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.48704874581882,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.74498578170074,-2.37629693876040,-0.342930291382793) q[4];
u3(1.34074367151570,2.19486316949611,-0.861757494849935) q[5];
u3(1.55566844161892,1.33054896349027,-3.01868879791491) q[2];
u3(2.40180790554847,-2.01560120848748,2.92288532797334) q[9];
cx q[9],q[2];
u1(1.53013370475628) q[2];
u3(-0.785817607010350,0.0,0.0) q[9];
cx q[2],q[9];
u3(-0.391400950574200,0.0,0.0) q[9];
cx q[9],q[2];
u3(0.946072952856584,1.00915564839897,-1.99494849068346) q[2];
u3(2.69101947514117,-1.98896043698558,-1.43402669742560) q[9];
u3(2.78632000356246,1.63762460026050,-4.53459543083661) q[3];
u3(1.19020389277514,-2.03554377813626,3.52038560610232) q[1];
cx q[1],q[3];
u1(2.34652053675008) q[3];
u3(0.402522996432758,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.24159090150174,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.03818251516965,-0.144429163627743,3.03318144813513) q[3];
u3(0.755990006340924,-4.35243729038257,0.858954689783084) q[1];
u3(2.01232595754447,1.14072745652963,-0.694619079321949) q[8];
u3(2.12668224395280,-0.444956332341175,-3.67998667730380) q[6];
cx q[6],q[8];
u1(0.179442819332933) q[8];
u3(-1.32329971642704,0.0,0.0) q[6];
cx q[8],q[6];
u3(2.57714510022975,0.0,0.0) q[6];
cx q[6],q[8];
u3(1.84286506097758,0.646638953934645,1.39413127867382) q[8];
u3(0.404512707699966,-1.61346519982498,-3.23958426745853) q[6];
u3(1.87553733950161,0.310866190760684,0.774630290441196) q[7];
u3(2.03252775296546,-0.558782322525977,-1.53289427445896) q[0];
cx q[0],q[7];
u1(0.243190172803853) q[7];
u3(-0.616473138855893,0.0,0.0) q[0];
cx q[7],q[0];
u3(1.27399369890564,0.0,0.0) q[0];
cx q[0],q[7];
u3(0.609766678893049,2.01947424189035,-1.16398010578760) q[7];
u3(1.76832670505903,-3.56261025814629,-1.03992722551883) q[0];
u3(1.33153593693057,-0.770642484009263,0.652816006642793) q[9];
u3(1.85029523177944,-1.46530861807799,-1.43101759953143) q[8];
cx q[8],q[9];
u1(-0.285718263973818) q[9];
u3(-2.48370847806731,0.0,0.0) q[8];
cx q[9],q[8];
u3(1.84414079612926,0.0,0.0) q[8];
cx q[8],q[9];
u3(1.57415641755301,-2.55172143728955,3.09758486017117) q[9];
u3(1.96773059138994,1.87321336708925,-2.72472409285435) q[8];
u3(0.424902135412357,0.263839852080440,0.0898855815590808) q[2];
u3(0.690323300677785,-0.995160455660766,-1.47218368383680) q[1];
cx q[1],q[2];
u1(1.08260655112419) q[2];
u3(-1.53441928279856,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.27710072760208,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.17311379230440,0.232776132524098,0.994834813986510) q[2];
u3(1.57869132039948,-0.166195831455737,-2.45220178583546) q[1];
u3(2.48104338347374,-1.33167446475313,1.70553556315073) q[5];
u3(2.68161199739639,2.77046399578786,3.03737370158834) q[4];
cx q[4],q[5];
u1(2.70797928996602) q[5];
u3(-1.74007938592895,0.0,0.0) q[4];
cx q[5],q[4];
u3(0.112803461079485,0.0,0.0) q[4];
cx q[4],q[5];
u3(2.89059923460668,-2.47647697362116,0.0284924901082502) q[5];
u3(2.46112648265867,1.96957702898439,-2.85067182175693) q[4];
u3(1.23445301794078,0.784960760144344,1.37233098930459) q[6];
u3(1.20853922874388,-0.491160361407960,-3.01625177907896) q[7];
cx q[7],q[6];
u1(4.44057961290254) q[6];
u3(-3.26678642772806,0.0,0.0) q[7];
cx q[6],q[7];
u3(-0.387889504629154,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.65245420465736,-3.55997475865908,0.164167205270796) q[6];
u3(1.65530213395875,-0.603775363810847,3.74238677806312) q[7];
u3(0.869149685024043,0.655046101267635,-0.934607530300972) q[3];
u3(0.130914409546620,-4.12376083161196,1.65870488978715) q[0];
cx q[0],q[3];
u1(1.08563953207170) q[3];
u3(-3.08442736423505,0.0,0.0) q[0];
cx q[3],q[0];
u3(1.66002622845813,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.73361782133372,1.69662466132887,1.22795263548095) q[3];
u3(1.42110052848122,-3.95133165438952,0.378633128863902) q[0];
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
