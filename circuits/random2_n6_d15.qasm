OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(0.646755592381333,3.18468898084356,-0.0507551798456214) q[4];
u3(1.93576387695561,3.54643696246037,0.828473226460575) q[5];
cx q[5],q[4];
u1(0.628449260818930) q[4];
u3(-1.30189829749412,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.71360859832827,0.0,0.0) q[5];
cx q[5],q[4];
u3(3.06723421713177,2.10943409458883,-2.69367494096860) q[4];
u3(1.71428465255978,4.96101267457260,0.794650387076114) q[5];
u3(1.70047378073116,0.352134291681904,2.36249854389070) q[1];
u3(3.02293317975692,-2.85074188620926,-2.14638456514493) q[0];
cx q[0],q[1];
u1(1.65172815112849) q[1];
u3(-0.341648334957163,0.0,0.0) q[0];
cx q[1],q[0];
u3(-0.191871846539827,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.00630848688685,-3.37402657346400,2.59049440684048) q[1];
u3(2.71142792755043,-1.21298893907741,2.98783804214021) q[0];
u3(2.10872251799518,2.34860795369417,0.323387777426556) q[3];
u3(1.81562201728608,0.426042854775644,-1.82352221581677) q[2];
cx q[2],q[3];
u1(3.39193102139404) q[3];
u3(-0.917396512835627,0.0,0.0) q[2];
cx q[3],q[2];
u3(2.00700738851355,0.0,0.0) q[2];
cx q[2],q[3];
u3(1.13821254148223,-0.690256029231764,0.401585280391204) q[3];
u3(1.41340130827250,4.44287893881098,1.64268526424023) q[2];
u3(2.56590973403342,-0.144631304246559,-1.04560661098320) q[5];
u3(1.79899235223524,1.74574841063319,-3.96608275531970) q[3];
cx q[3],q[5];
u1(0.537391006190830) q[5];
u3(-1.48371907534598,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.358062641237442,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.00168831901228,-1.99874362625460,1.06170815924694) q[5];
u3(0.981486290605868,-0.802606257437594,-4.44844639806668) q[3];
u3(0.705052920925648,2.39718038191840,-3.39300407283621) q[4];
u3(1.30912609197587,2.56440511210161,-2.96957082316589) q[2];
cx q[2],q[4];
u1(1.53303594906498) q[4];
u3(-2.58672302122227,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.40533604835858,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.08020818244187,0.172440645830447,-0.794197148228589) q[4];
u3(2.78251624928657,-4.70426106786491,-0.539613394751161) q[2];
u3(2.35536460988825,0.860985431754852,0.501988554355766) q[0];
u3(1.03283186599389,-1.95091252759435,-1.80250308100195) q[1];
cx q[1],q[0];
u1(2.33097463137255) q[0];
u3(0.0814382127174171,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.959856043544116,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.17720576241639,0.368070945028153,0.703629186402225) q[0];
u3(0.987950758911383,4.75929757079576,0.719202850865839) q[1];
u3(1.79345817897020,1.23815729096421,-0.857045964549030) q[4];
u3(1.28639416020267,-4.60054610021482,0.800686330567063) q[5];
cx q[5],q[4];
u1(3.49778309945230) q[4];
u3(-1.35888319575992,0.0,0.0) q[5];
cx q[4],q[5];
u3(2.52040783902295,0.0,0.0) q[5];
cx q[5],q[4];
u3(1.49464338240860,0.697496821759018,1.58318372494683) q[4];
u3(0.674981752654153,0.839330840199586,-0.140403367699547) q[5];
u3(2.25956075967447,-2.86393658478698,0.598678801021003) q[0];
u3(2.51602218712914,1.41744401751051,2.06331979571457) q[2];
cx q[2],q[0];
u1(2.63182999151871) q[0];
u3(-2.28874897035411,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.41274641074165,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.19231123354111,-1.82982701088374,-1.79632056041447) q[0];
u3(0.316929542034524,2.77064998785756,-1.98379852497801) q[2];
u3(1.92978873060984,2.40371901000386,-2.18117551868538) q[3];
u3(0.921740349769189,-3.12295541657056,2.14572088514654) q[1];
cx q[1],q[3];
u1(0.338120507725213) q[3];
u3(-1.23862280838840,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.60753514614005,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.18004418344512,2.72997656328901,1.31810486699472) q[3];
u3(1.82567906086736,-0.0541388567894949,3.53497887263154) q[1];
u3(1.36931623986807,0.506023229388821,0.840654991753047) q[0];
u3(2.04249925437288,-0.747781830923023,-2.56016356757810) q[5];
cx q[5],q[0];
u1(-0.269890764923623) q[0];
u3(-2.40478556775547,0.0,0.0) q[5];
cx q[0],q[5];
u3(1.62068381254316,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.51289643353918,-1.70939851778395,2.05348267352406) q[0];
u3(1.42406770598540,-2.31391274502790,-1.04497997996581) q[5];
u3(2.05908779659850,2.48980274497524,-3.53223814774175) q[2];
u3(2.43061433227454,2.31254858648092,-3.42227721228779) q[3];
cx q[3],q[2];
u1(1.41469217718927) q[2];
u3(-0.228338528482115,0.0,0.0) q[3];
cx q[2],q[3];
u3(0.974476496379430,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.57737677739537,-4.01018384523625,0.804006179791975) q[2];
u3(1.46311817836870,2.00816784586229,2.62150703691430) q[3];
u3(1.52209704005177,-0.405806922655305,1.07147484538734) q[1];
u3(1.99832334349295,-2.14804086505294,-1.65285379357158) q[4];
cx q[4],q[1];
u1(1.54165809894244) q[1];
u3(-2.52242552034945,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.159076155087406,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.87160405055700,3.45091039471986,-0.928861653063929) q[1];
u3(2.32431009541549,-1.18584313449652,-3.98752663327547) q[4];
u3(2.32102679047832,2.19588843784901,-2.21382875318350) q[1];
u3(2.15166957432010,-0.619702794880677,-5.52135047511112) q[5];
cx q[5],q[1];
u1(3.50681348440004) q[1];
u3(-1.38182805143968,0.0,0.0) q[5];
cx q[1],q[5];
u3(2.42994137371824,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.60445847607144,-1.64136667182609,3.32982306469199) q[1];
u3(1.97413170434490,-0.634499737927804,-4.71598493021146) q[5];
u3(1.91875184045122,-0.967993124454983,0.0156872271645267) q[2];
u3(2.51150840399861,-2.25115081786539,1.12909193581594) q[0];
cx q[0],q[2];
u1(0.671256473194764) q[2];
u3(-0.262624907969056,0.0,0.0) q[0];
cx q[2],q[0];
u3(2.19354892620631,0.0,0.0) q[0];
cx q[0],q[2];
u3(2.68079998404834,2.25510900645149,-2.71018556177070) q[2];
u3(1.57723717481443,4.02310487839606,1.09369443120312) q[0];
u3(2.52012140636149,1.36840014526740,1.57969324052174) q[3];
u3(1.12714176076895,0.183363885612926,-5.83871421115575) q[4];
cx q[4],q[3];
u1(-0.683292796524074) q[3];
u3(0.320095150217613,0.0,0.0) q[4];
cx q[3],q[4];
u3(3.26129202061221,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.49234755553690,-2.82790569054771,2.15651189653571) q[3];
u3(2.05473594089378,4.16038207461003,-1.75142420936317) q[4];
u3(0.996324871310747,-2.90018894204934,3.07028499747467) q[0];
u3(0.351924090227412,-0.0997438197115852,-2.64556000142881) q[1];
cx q[1],q[0];
u1(2.65427677454645) q[0];
u3(-0.0362596682494807,0.0,0.0) q[1];
cx q[0],q[1];
u3(1.11437926204590,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.76027378082155,2.10300077742960,-1.66964216422797) q[0];
u3(2.63752906252956,-3.02032830813483,1.37452701685509) q[1];
u3(1.99998705846017,-3.59915659820158,2.55928773646876) q[3];
u3(0.121250964640024,2.95849334805948,-1.13707918246181) q[5];
cx q[5],q[3];
u1(0.656467366596840) q[3];
u3(-1.28410073237166,0.0,0.0) q[5];
cx q[3],q[5];
u3(3.18777315760958,0.0,0.0) q[5];
cx q[5],q[3];
u3(1.74894477539908,1.13928606634542,-2.35904285389179) q[3];
u3(1.01169392230756,1.02250359873766,-3.14675841083552) q[5];
u3(2.70171216109386,-0.441264125552778,1.62911211296202) q[4];
u3(2.51391712847825,-1.24647172725330,0.635002184497000) q[2];
cx q[2],q[4];
u1(2.88784924605463) q[4];
u3(-2.18312476636234,0.0,0.0) q[2];
cx q[4],q[2];
u3(-0.0254813863272101,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.549849651421865,-0.618811952586448,4.87998286704982) q[4];
u3(0.153365090620267,0.387349655282605,5.67778844756325) q[2];
u3(1.30601862939489,1.60954873215594,-2.80364430088269) q[4];
u3(0.353817656633818,1.32392019170761,-2.77364666545663) q[3];
cx q[3],q[4];
u1(1.18624726589892) q[4];
u3(-3.23096501193472,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.65422842037954,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.702931353385772,-1.36121113829698,1.62606407945971) q[4];
u3(1.80678524091057,2.36735920902442,2.70070934638101) q[3];
u3(2.57426043516453,-1.38161263521354,3.47680354567744) q[0];
u3(1.91214164386977,0.903139890555287,1.91647999639373) q[5];
cx q[5],q[0];
u1(1.77752291618389) q[0];
u3(-2.67173808844383,0.0,0.0) q[5];
cx q[0],q[5];
u3(3.05288842105757,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.96306777001643,2.58811187454012,-2.34031064634007) q[0];
u3(2.48145997901519,-1.20527845181487,4.58908465995420) q[5];
u3(1.72080211831336,0.277556900704515,1.47366601754637) q[2];
u3(1.27775899095306,-0.951720357235534,-2.74522428557952) q[1];
cx q[1],q[2];
u1(1.62268711308159) q[2];
u3(-3.28763294391034,0.0,0.0) q[1];
cx q[2],q[1];
u3(2.56270629766323,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.82754112397425,-0.625875945993265,-0.598433497274713) q[2];
u3(2.27182443576308,-5.12689618419556,-0.694177780400161) q[1];
u3(1.65346401756714,1.07742622666950,-3.37181116605883) q[1];
u3(2.10434944883380,2.49490391603968,-3.15938968011983) q[3];
cx q[3],q[1];
u1(3.37438957703480) q[1];
u3(-1.30466143948420,0.0,0.0) q[3];
cx q[1],q[3];
u3(2.46909132482120,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.56914954550314,0.396101391055616,-1.02603575403376) q[1];
u3(2.25633046630668,-4.04593176953274,-0.622795950947565) q[3];
u3(1.86481261129141,-3.40375888348730,2.43875918785698) q[4];
u3(0.799663676900966,3.35906822625181,-2.41309208481550) q[0];
cx q[0],q[4];
u1(2.07805409180165) q[4];
u3(0.0285984145143761,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.55161170217937,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.82050013952827,4.56975836332158,-1.41769073749562) q[4];
u3(0.765479665344573,-2.62383157940214,-1.96244092727190) q[0];
u3(1.17264085891152,-0.308166730013293,1.76752760197083) q[5];
u3(0.786550040578833,-1.72346332341945,-2.22161229559965) q[2];
cx q[2],q[5];
u1(1.39644638714743) q[5];
u3(-0.0919686879582817,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.33098380238402,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.861651636149598,1.49234077139398,0.236887703833072) q[5];
u3(2.22820393155092,-1.99877916001837,3.22550177072732) q[2];
u3(1.40100807472714,1.24233511306537,-1.91692059242694) q[0];
u3(1.18715188989238,-2.01228211567033,2.77068076426262) q[2];
cx q[2],q[0];
u1(2.42747351723214) q[0];
u3(-1.59146652014925,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.241160713728490,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.819578178696473,0.558118785660317,-1.02091783349827) q[0];
u3(2.91974556814125,-3.99629692581579,1.79726874413959) q[2];
u3(1.53820624401971,0.745524575135172,0.719127250177302) q[5];
u3(1.04543950565875,-1.21245374874602,-2.63063944640953) q[4];
cx q[4],q[5];
u1(0.0825046278993735) q[5];
u3(-0.760548451119027,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.90406538596885,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.24943765832351,1.70025199031858,-1.42756954400521) q[5];
u3(0.893489131987865,-1.14043149000488,1.87990519490561) q[4];
u3(0.200697818607550,-3.11983696391649,2.51027467518459) q[3];
u3(0.401556576411898,0.150363481179937,-1.58859670057198) q[1];
cx q[1],q[3];
u1(3.04130182820729) q[3];
u3(-2.04832440968645,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.501174889606981,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.737309548870124,2.28075341614442,0.559408905006354) q[3];
u3(1.70102497452095,-1.64713233768241,-3.17949544078995) q[1];
u3(0.891012354872804,-0.0851295544706900,-0.876107165200267) q[0];
u3(2.06083086673096,0.481943525259744,-5.39929539084409) q[5];
cx q[5],q[0];
u1(1.89166936068910) q[0];
u3(-2.75936354497941,0.0,0.0) q[5];
cx q[0],q[5];
u3(0.596704717218601,0.0,0.0) q[5];
cx q[5],q[0];
u3(1.46289150642675,-0.0111936645604811,1.12447287857047) q[0];
u3(1.03129094519786,-0.754767914797862,4.19808459113446) q[5];
u3(2.68282716703439,-0.418672308436264,-0.902879725239952) q[2];
u3(1.14599579249397,-2.38716271876984,-1.79803457170349) q[1];
cx q[1],q[2];
u1(1.99825670719590) q[2];
u3(-2.51220873491736,0.0,0.0) q[1];
cx q[2],q[1];
u3(3.19292617530436,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.46299921575833,-1.06271838671784,1.11463586679515) q[2];
u3(0.365054478453603,-0.991596335204888,-2.23762334012186) q[1];
u3(1.28862455189653,0.545443127678173,-2.75718544462475) q[4];
u3(1.06989153065362,2.77874420085177,-3.45620545430856) q[3];
cx q[3],q[4];
u1(1.62416703521169) q[4];
u3(-2.31879366633468,0.0,0.0) q[3];
cx q[4],q[3];
u3(3.69191126843891,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.15347025953265,1.54501042554109,1.61196258033929) q[4];
u3(2.41745146459345,1.84955308313570,-3.42130225260050) q[3];
u3(1.22609048896251,3.39716197914382,-0.967413892795735) q[4];
u3(0.600150118961459,1.29977996648963,-1.31674524936147) q[0];
cx q[0],q[4];
u1(1.79240024957522) q[4];
u3(-3.20095215514475,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.380920554271436,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.365442876240777,0.881286800275274,-1.04234214164820) q[4];
u3(2.51484168151989,-0.264214173469221,-4.52886541617634) q[0];
u3(1.45007635711969,-1.92424103375550,-1.15826780634443) q[5];
u3(2.21341321256341,-2.76439434851475,-0.0749228223924236) q[1];
cx q[1],q[5];
u1(2.68933255456661) q[5];
u3(-2.42396276116987,0.0,0.0) q[1];
cx q[5],q[1];
u3(1.13973710215605,0.0,0.0) q[1];
cx q[1],q[5];
u3(0.941067781186026,-0.128802798736868,-0.887091696835441) q[5];
u3(0.307311072467413,2.68953882009747,2.56026022895764) q[1];
u3(0.714455209192435,3.47507322059080,-2.35617743879523) q[2];
u3(0.822540149432277,-2.77814612433185,0.967083707576450) q[3];
cx q[3],q[2];
u1(1.38139966837769) q[2];
u3(-0.540433319843974,0.0,0.0) q[3];
cx q[2],q[3];
u3(3.19516291954275,0.0,0.0) q[3];
cx q[3],q[2];
u3(0.798365693280659,0.157454638688928,3.36479674255586) q[2];
u3(1.81121198613385,2.15467440806452,3.78364219274745) q[3];
u3(1.20937710952904,0.733474688345537,-1.00339154954221) q[4];
u3(0.385814671655893,-0.938571760283180,-1.23016567640315) q[0];
cx q[0],q[4];
u1(1.64055327079752) q[4];
u3(-0.686116687585410,0.0,0.0) q[0];
cx q[4],q[0];
u3(2.93974440696404,0.0,0.0) q[0];
cx q[0],q[4];
u3(0.847120899956113,-2.72204054727188,2.09047133778130) q[4];
u3(1.29111537719151,4.60841322329819,-1.33239398516360) q[0];
u3(0.938910042629025,-1.14288676084014,0.163726410479853) q[2];
u3(0.804869317781866,-1.93102359629475,0.201593725049244) q[5];
cx q[5],q[2];
u1(1.42295437897006) q[2];
u3(-0.972212660335292,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.52015453003266,0.0,0.0) q[5];
cx q[5],q[2];
u3(0.323476277219774,2.25642668845342,-3.10304978078846) q[2];
u3(2.13347149503598,2.43685243897634,-3.44296943787137) q[5];
u3(1.65672759525942,0.123748716479048,1.24869269218603) q[3];
u3(1.85502729236073,-1.55930444492211,-0.374386978981189) q[1];
cx q[1],q[3];
u1(1.46504594617670) q[3];
u3(-2.17527294136695,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.736526004552379,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.49148587569158,1.44543409137531,-2.11624347841919) q[3];
u3(0.368199663305245,2.86660299756683,2.03292613525101) q[1];
u3(1.16529551244208,1.36959329444373,-2.98935915269958) q[3];
u3(2.15665455517036,2.67998192488160,-2.83244127664191) q[1];
cx q[1],q[3];
u1(0.729363202864307) q[3];
u3(-0.0986365016564936,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.85296080462192,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.07723736845691,-0.989978888711640,0.137356406866346) q[3];
u3(2.33210038442462,3.08117674454254,-3.07041427036105) q[1];
u3(0.738585864501151,1.62779004568535,-2.86088972109567) q[4];
u3(1.49335513225632,-2.45273472279263,3.42865614852922) q[5];
cx q[5],q[4];
u1(1.10356567900203) q[4];
u3(-3.07153001153926,0.0,0.0) q[5];
cx q[4],q[5];
u3(1.69322676841720,0.0,0.0) q[5];
cx q[5],q[4];
u3(0.493231656743357,-0.914062507494682,1.85126890502987) q[4];
u3(0.105313592111788,2.38800744396795,-2.58213643655853) q[5];
u3(2.85171606224244,-1.72040231065113,1.95439958006599) q[0];
u3(1.90319943897710,-1.55193682936452,0.173811848906548) q[2];
cx q[2],q[0];
u1(1.08732890734063) q[0];
u3(-0.437142722608625,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.26067060232469,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.0698677785586271,-0.263726679773731,-2.27534531092290) q[0];
u3(1.83346730756598,-1.10334114762256,-1.74165647202742) q[2];
u3(2.37385116413278,-1.87283456863386,0.544931312466793) q[4];
u3(2.84121512738586,-2.53519267265831,-1.66794103291290) q[0];
cx q[0],q[4];
u1(0.227275291816344) q[4];
u3(-0.963172935729840,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.68410257614051,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.17672509533363,-1.15026975274183,2.27361146251524) q[4];
u3(2.54460247902970,-4.36393085816053,1.49043166649394) q[0];
u3(1.06590746607271,1.07273303371874,-1.81855746682526) q[2];
u3(1.95930347980992,-4.51655861742412,1.44839556515063) q[1];
cx q[1],q[2];
u1(0.644340769640754) q[2];
u3(-3.54449137998040,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.73279819386608,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.50259019803538,-1.48595842024828,3.12343827534483) q[2];
u3(1.34399934321571,-4.61015264843528,0.419295386169153) q[1];
u3(2.83385427016562,-1.38440795051200,-1.69035381792843) q[5];
u3(1.14928235450267,-2.26191440191653,-1.86476865108463) q[3];
cx q[3],q[5];
u1(3.21595838685382) q[5];
u3(-4.33948792925344,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.357338880192229,0.0,0.0) q[3];
cx q[3],q[5];
u3(0.957494126291886,1.53016429246901,-3.03115914177873) q[5];
u3(1.02394417457460,-3.48351402310404,2.39914412191931) q[3];
u3(1.76298481869370,2.87042289201764,0.193238590735726) q[5];
u3(2.07749091380852,0.0399549413547904,-4.80782834665535) q[2];
cx q[2],q[5];
u1(3.33711855584699) q[5];
u3(-0.883457852397823,0.0,0.0) q[2];
cx q[5],q[2];
u3(1.72375330021148,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.892934784274240,1.35994771872761,-1.79571115812421) q[5];
u3(1.54293443608058,-0.0669088620583087,3.32728976700636) q[2];
u3(1.95619512374171,2.34974693265682,-2.48839325454407) q[3];
u3(0.779622867660970,2.17994393925855,-2.31065980069916) q[4];
cx q[4],q[3];
u1(1.57701340453088) q[3];
u3(-3.75171835622383,0.0,0.0) q[4];
cx q[3],q[4];
u3(1.98168197263003,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.68603055683645,4.46988701158958,0.0522775025062301) q[3];
u3(2.19133553461333,-1.25360677190248,1.58806303365317) q[4];
u3(0.323876736327216,-1.20979652579164,0.878261558024303) q[1];
u3(1.12255463604538,-2.73946829144084,-0.0119454930147553) q[0];
cx q[0],q[1];
u1(2.09952298159328) q[1];
u3(-3.14681355003624,0.0,0.0) q[0];
cx q[1],q[0];
u3(0.498077834007161,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.35648335390411,0.667990395913412,-1.60101033934953) q[1];
u3(1.71763028951794,0.916296723719266,-1.40057629105026) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
