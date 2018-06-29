OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
u3(2.78044021679119,1.87356315259913,0.168471010968469) q[6];
u3(1.50236730819628,-0.959096938105201,-2.60335412542689) q[4];
cx q[4],q[6];
u1(-0.688934390982197) q[6];
u3(-1.66916356547755,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.953665661108094,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.11287385617145,-1.09976552964049,4.34723219244874) q[6];
u3(1.28257575898590,-2.86114639377771,-3.41041519692266) q[4];
u3(2.22896772649366,1.62724972704930,1.40861486763701) q[3];
u3(0.213514665826181,-2.93029736710519,-0.906467228698948) q[7];
cx q[7],q[3];
u1(1.46872870128004) q[3];
u3(-2.44245034816934,0.0,0.0) q[7];
cx q[3],q[7];
u3(0.519033107719403,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.33490011064150,2.06519740744663,-3.97049473272335) q[3];
u3(0.944024806049953,0.354756482044631,-3.67947118807621) q[7];
u3(2.02839485536539,-2.10866642619179,0.553234172194199) q[0];
u3(2.91046526736969,-3.34712742849731,-1.94884524500503) q[2];
cx q[2],q[0];
u1(0.885300944624110) q[0];
u3(-1.41743342053581,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.15095908244984,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.64187391098556,2.00253532431680,0.423296082420042) q[0];
u3(0.970465692025623,4.32507812425896,-0.152031852054319) q[2];
u3(1.67890493757843,3.27164085951368,-0.823596600525387) q[1];
u3(0.920719135110853,2.20648298130685,-2.13960798779038) q[5];
cx q[5],q[1];
u1(1.68756611102192) q[1];
u3(-2.85199963536900,0.0,0.0) q[5];
cx q[1],q[5];
u3(0.978765608364123,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.93083648332046,-1.09915826217582,-0.981478270767305) q[1];
u3(0.133611523609736,0.744048434004327,3.65823937737818) q[5];
u3(2.16688222633266,2.75540251909683,-3.26702523232725) q[5];
u3(1.03604229870922,-0.146470969060963,1.27032152209625) q[1];
cx q[1],q[5];
u1(-0.154330813670435) q[5];
u3(-1.88943954246110,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.665434158137107,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.66877633916046,4.18737821058136,-1.83833067053411) q[5];
u3(0.633202785300870,-1.64954684298446,3.91133554022566) q[1];
u3(1.11403527420100,0.612540818524380,1.71406861948091) q[7];
u3(1.33722760048282,-1.76502125014168,-0.867266636519728) q[2];
cx q[2],q[7];
u1(0.774026063912134) q[7];
u3(-3.51638720996254,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.21492321382980,0.0,0.0) q[2];
cx q[2],q[7];
u3(1.10224150957430,4.06705795035766,-1.88519874949008) q[7];
u3(2.38662131891368,-5.21528221093930,1.03642399927083) q[2];
u3(1.51971608824596,1.20224001962896,-2.92051920869200) q[0];
u3(0.770114132878428,1.91475632821798,-2.01668755220572) q[4];
cx q[4],q[0];
u1(4.34347966278819) q[0];
u3(-3.80910999026630,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.936453796687176,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.56338494160104,2.06961266748316,-2.68663580653713) q[0];
u3(2.69670244490397,-2.92599274771212,2.90873736281890) q[4];
u3(1.89440920219329,-0.136806534890517,1.48383191316238) q[6];
u3(1.99075823104743,-2.34167931203459,-2.30253980383794) q[3];
cx q[3],q[6];
u1(-0.654071640459954) q[6];
u3(1.34139503455485,0.0,0.0) q[3];
cx q[6],q[3];
u3(3.64404524984247,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.86069521401003,1.71182295323380,-4.20596469639952) q[6];
u3(0.758898164384018,3.55438039324902,0.986097702740122) q[3];
u3(0.127238274274214,1.66419801721973,-1.62437236898275) q[6];
u3(0.0635232685066098,0.977347036666594,-2.75943559620544) q[0];
cx q[0],q[6];
u1(1.36778873709174) q[6];
u3(0.236051728451778,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.19635119367994,0.0,0.0) q[0];
cx q[0],q[6];
u3(1.01424508530914,-0.435896692608693,-0.476311379187825) q[6];
u3(1.67361716569666,1.07473763001365,-3.98965674435314) q[0];
u3(2.96764028527654,-0.300263561608103,2.40894143486402) q[2];
u3(2.52146611454640,-0.611581816175889,0.0272337246607385) q[7];
cx q[7],q[2];
u1(-0.0423078878345315) q[2];
u3(-2.48211545491440,0.0,0.0) q[7];
cx q[2],q[7];
u3(1.11996590288687,0.0,0.0) q[7];
cx q[7],q[2];
u3(0.799154093368076,-0.534933575187695,3.71039144919855) q[2];
u3(2.55063091498251,2.86339766130874,1.76289005426823) q[7];
u3(1.32409789579909,-0.748569517494488,1.56703121035629) q[4];
u3(1.21832213977461,-1.06069101658922,-1.40309834603316) q[5];
cx q[5],q[4];
u1(1.58761870242511) q[4];
u3(-2.53227353402692,0.0,0.0) q[5];
cx q[4],q[5];
u3(0.243636834563237,0.0,0.0) q[5];
cx q[5],q[4];
u3(2.58282050262432,1.09011263939383,-2.35323134048395) q[4];
u3(0.366920412752013,0.828706434417007,-5.09787257081628) q[5];
u3(1.87998943649597,1.09170045172827,-1.93801041248265) q[1];
u3(1.90206298478444,1.76487574885779,-4.32366091152833) q[3];
cx q[3],q[1];
u1(2.11031177344579) q[1];
u3(-1.31428818390235,0.0,0.0) q[3];
cx q[1],q[3];
u3(0.719879525637108,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.41192056398017,-3.62833964600145,1.43453376477795) q[1];
u3(1.24673130636332,-0.736029111071107,-0.840819863434449) q[3];
u3(1.07361360298038,2.06909341864091,-0.492879341734347) q[7];
u3(0.590656027418282,1.33992070667817,-3.71751964641903) q[1];
cx q[1],q[7];
u1(-0.158461306450419) q[7];
u3(-0.805422536805138,0.0,0.0) q[1];
cx q[7],q[1];
u3(1.98802459200093,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.45858703143840,-0.242915451966130,-0.675057091472807) q[7];
u3(1.56048839243887,-1.64503947920621,1.61003163313107) q[1];
u3(1.25457541922629,1.48776176652650,-3.70508016392796) q[0];
u3(0.742762954531295,2.05135160166317,-1.79542664664838) q[4];
cx q[4],q[0];
u1(0.637623256607047) q[0];
u3(-3.33345398892531,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.87867430273697,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.32140789669912,-1.94075653108577,0.883973972013677) q[0];
u3(1.27234079080062,-3.50993901897015,1.17264508644784) q[4];
u3(1.59887223952879,0.0135151638139908,2.74080814156191) q[3];
u3(0.607202484205987,-0.550514039584337,-1.84193359281977) q[2];
cx q[2],q[3];
u1(-0.0179159629486316) q[3];
u3(-1.71861573911082,0.0,0.0) q[2];
cx q[3],q[2];
u3(0.562170136719235,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.46727364168741,3.55565337846671,-1.54341481131541) q[3];
u3(2.50741930596123,-3.78249643592069,0.00117719936644556) q[2];
u3(1.91273784018565,-0.609211704845440,2.09626921527921) q[6];
u3(1.79196295728922,-2.52199939473991,-2.36657059196544) q[5];
cx q[5],q[6];
u1(-0.195494970073097) q[6];
u3(-1.90825296965299,0.0,0.0) q[5];
cx q[6],q[5];
u3(1.09206122555226,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.04305437274783,-2.19958129893373,1.68887401420548) q[6];
u3(2.35270422537690,-2.33619723582730,-0.471688669968640) q[5];
u3(1.66613413918897,2.06529025379959,-3.47248463352665) q[5];
u3(1.39568021297040,3.50387115571427,-2.74085417374041) q[4];
cx q[4],q[5];
u1(0.347609292081324) q[5];
u3(-1.63761734734339,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.38899960554599,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.28267197024709,-0.849536868479544,2.30329011013433) q[5];
u3(0.0876987833350168,0.00252948463333103,3.53925928636212) q[4];
u3(1.73350881429579,0.473880898250160,-2.20576647212090) q[7];
u3(2.76816121572348,2.49933732103524,-3.16065947453826) q[1];
cx q[1],q[7];
u1(2.36341200086520) q[7];
u3(-1.66863978116990,0.0,0.0) q[1];
cx q[7],q[1];
u3(0.296056240654009,0.0,0.0) q[1];
cx q[1],q[7];
u3(2.09957100297374,2.74816174949054,-0.746624116793026) q[7];
u3(0.397240229063553,-0.182849263533222,-2.55489965688828) q[1];
u3(2.25764259000667,-0.221555754841711,-0.938065576184697) q[3];
u3(0.741949078704602,0.346311797226125,-5.33280121641867) q[6];
cx q[6],q[3];
u1(2.21229124736302) q[3];
u3(0.232655750146348,0.0,0.0) q[6];
cx q[3],q[6];
u3(1.47216980186143,0.0,0.0) q[6];
cx q[6],q[3];
u3(2.36397243207651,-1.22092587228526,-0.587094596546120) q[3];
u3(0.952028887869615,0.220120745605319,-0.383091627824859) q[6];
u3(1.10736295968759,-0.401449664025150,-0.834645725230026) q[0];
u3(1.72231120603503,-4.78304800806503,1.11208991602981) q[2];
cx q[2],q[0];
u1(-0.472384020137636) q[0];
u3(1.10386284289143,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.73288658209794,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.45721563190493,-2.61258046694955,1.24451123876870) q[0];
u3(2.48804870264932,-3.39093024658853,-0.841867395191649) q[2];
u3(2.60814731023964,3.09128381871665,-2.73426298504511) q[3];
u3(1.96038324878290,1.74156262453010,-1.65519803273976) q[1];
cx q[1],q[3];
u1(1.63107247159487) q[3];
u3(-2.97488754745066,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.904052510151516,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.76802430839801,0.496966393345220,0.640051643362868) q[3];
u3(0.760604258021561,3.18376499587269,-0.119744366377004) q[1];
u3(2.54817251010275,3.59142501262728,-2.29274050509993) q[2];
u3(0.884915134658708,-0.475244733095195,1.64932979038162) q[4];
cx q[4],q[2];
u1(3.06484970483218) q[2];
u3(-2.66875281516942,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.931188268878133,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.38365366608328,0.772268696846410,2.90377206501882) q[2];
u3(1.55622016195014,5.59282876545596,0.0524024492852160) q[4];
u3(1.56483512145742,-0.307019685155140,1.03991039880323) q[7];
u3(1.82223927460288,-1.29976287465430,-1.96668319096137) q[5];
cx q[5],q[7];
u1(0.709849442500808) q[7];
u3(-3.15943264321893,0.0,0.0) q[5];
cx q[7],q[5];
u3(1.10326233293480,0.0,0.0) q[5];
cx q[5],q[7];
u3(1.49272645878651,4.17410074534215,-1.18094945900067) q[7];
u3(1.82982622505893,-3.11320829302309,-2.82927813206025) q[5];
u3(1.94784586467867,-2.15745435101274,0.0670204070524862) q[0];
u3(2.27673602773580,-3.57202532791180,-1.29226959887360) q[6];
cx q[6],q[0];
u1(2.34576243570561) q[0];
u3(0.316235425984131,0.0,0.0) q[6];
cx q[0],q[6];
u3(1.15980940547519,0.0,0.0) q[6];
cx q[6],q[0];
u3(1.85916108337036,-1.92425801967325,4.07511660564184) q[0];
u3(1.46966100214950,-2.46555604494737,-0.404039660539006) q[6];
u3(1.17355968351969,0.733833563096087,-1.25320135423972) q[7];
u3(0.664387282558568,-1.87103471094070,0.232284617033347) q[2];
cx q[2],q[7];
u1(2.71891533607071) q[7];
u3(-2.87230156244486,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.20909645907764,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.09356849804215,-0.603903730793847,-2.24703070225594) q[7];
u3(1.27807082145137,-2.13130072572301,1.26298066712502) q[2];
u3(1.02692764237289,-1.20541921064866,-0.620715237519488) q[6];
u3(2.45600836671056,0.688785009645126,-4.53698467835534) q[5];
cx q[5],q[6];
u1(0.145857646860773) q[6];
u3(-1.93751016175897,0.0,0.0) q[5];
cx q[6],q[5];
u3(0.898139079298256,0.0,0.0) q[5];
cx q[5],q[6];
u3(2.03549559345461,3.75658742554086,-2.33055554892698) q[6];
u3(2.59856801770563,2.70107134538948,3.13095168976183) q[5];
u3(0.256456627628947,2.98398987623745,-1.94369782742837) q[0];
u3(1.68665603799739,2.05305652292686,-2.10013875846780) q[1];
cx q[1],q[0];
u1(2.83577219137018) q[0];
u3(-1.92954005823444,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.780014898791232,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.28936601530110,-0.346362332959840,0.139335494326479) q[0];
u3(0.888265360263794,-3.77161900890993,2.09885878675065) q[1];
u3(2.14640928984976,-1.36168829367061,3.88165634230241) q[4];
u3(1.38999424767111,1.51757981772145,2.48813783735411) q[3];
cx q[3],q[4];
u1(-0.127373539593602) q[4];
u3(-1.71643687607934,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.16964214502227,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.51929195957179,-0.986189116312539,0.841121027520063) q[4];
u3(2.10816633470811,-0.314146872863639,5.63939947735491) q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
