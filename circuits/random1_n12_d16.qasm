OPENQASM 2.0;
include "qelib1.inc";
qreg q[12];
creg c[12];
u3(1.56710526059430,-1.54917799788562,0.655128967349044) q[10];
u3(1.62935078531105,-2.01028836852328,-0.816454572567026) q[9];
cx q[9],q[10];
u1(2.24458455318564) q[10];
u3(-1.66462600851137,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.0920185332280987,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.43856114751571,-1.74338407528739,-0.689207590472663) q[10];
u3(2.50176987064066,0.259727528770590,-1.02660084178781) q[9];
u3(1.39745029526274,1.33901385545559,1.75660483207388) q[5];
u3(1.47217660283809,-1.61496989206819,-0.325015741257937) q[1];
cx q[1],q[5];
u1(3.49431734738193) q[5];
u3(-1.46680359969535,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.41814351006635,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.74366524781042,0.398896064566522,2.38508936088144) q[5];
u3(1.13320416951773,1.21152109343980,2.42528256936326) q[1];
u3(0.817618236469471,1.37608780382415,-4.04002610762452) q[4];
u3(1.96315055331619,2.41563561264447,-2.92198501629755) q[6];
cx q[6],q[4];
u1(1.43887329860274) q[4];
u3(-3.40425666208526,0.0,0.0) q[6];
cx q[4],q[6];
u3(2.33955128049757,0.0,0.0) q[6];
cx q[6],q[4];
u3(0.960531329305930,2.24352892176812,-1.50910872468392) q[4];
u3(1.96649490149698,1.48144275869927,3.65063516746724) q[6];
u3(1.61872413847852,-1.63006196203277,-0.864615665472085) q[0];
u3(1.23236614938218,-4.24929636177111,0.321206102189344) q[8];
cx q[8],q[0];
u1(4.25510660366763) q[0];
u3(-4.00331191167296,0.0,0.0) q[8];
cx q[0],q[8];
u3(-0.804663595446804,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.93043453793583,-2.65888666747625,0.0848218278233988) q[0];
u3(2.24104271172004,-0.496165229418968,-1.41405476897724) q[8];
u3(1.51222359786778,1.08628602689113,-3.65280773187610) q[3];
u3(2.76328119129565,-2.00475044931790,3.74659448878297) q[11];
cx q[11],q[3];
u1(1.62460553201539) q[3];
u3(-2.71975077907120,0.0,0.0) q[11];
cx q[3],q[11];
u3(3.43440635692605,0.0,0.0) q[11];
cx q[11],q[3];
u3(0.694101082134221,3.38632546939852,-0.907014660709320) q[3];
u3(1.61194073055667,-1.88500272667205,1.92875640865520) q[11];
u3(1.21230486929101,1.79033215172658,-0.575144049677168) q[7];
u3(1.11879847214776,0.391447454516268,-4.05631226896941) q[2];
cx q[2],q[7];
u1(4.07258882394361) q[7];
u3(-3.18168372604615,0.0,0.0) q[2];
cx q[7],q[2];
u3(-0.582763444169669,0.0,0.0) q[2];
cx q[2],q[7];
u3(0.725291067163767,-0.537987563730349,-0.117402684077309) q[7];
u3(1.41920201835406,0.390695241310734,4.66270252378904) q[2];
u3(0.372898599904677,-0.148203661733752,-1.74345207879740) q[10];
u3(1.30540665895732,-4.21570211861438,1.64591351173777) q[9];
cx q[9],q[10];
u1(1.81165393032254) q[10];
u3(0.431225885954810,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.814698590960490,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.00902685749340,0.887626113555319,-3.46232339028270) q[10];
u3(0.671354900422512,1.35107542572854,4.31310595648246) q[9];
u3(1.02211588785023,2.29079923632428,-0.204981116243847) q[0];
u3(1.13374436885420,1.78640579189379,-1.45871920736330) q[6];
cx q[6],q[0];
u1(0.0787933181387632) q[0];
u3(-1.45981098745545,0.0,0.0) q[6];
cx q[0],q[6];
u3(2.43698146402591,0.0,0.0) q[6];
cx q[6],q[0];
u3(0.729195333741783,4.01014063760363,-1.41497161414374) q[0];
u3(1.47445421281754,4.09232335516937,-0.449267710793765) q[6];
u3(0.782212202679433,1.03074943491574,-0.899221566015706) q[1];
u3(0.377192518592368,-0.517741358112192,-1.08763734449715) q[4];
cx q[4],q[1];
u1(0.0617324780502946) q[1];
u3(-2.45948696111227,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.43994439014945,0.0,0.0) q[4];
cx q[4],q[1];
u3(2.19028043256046,-2.53457483611106,1.15194119194678) q[1];
u3(1.29567346789958,-5.66209936429517,-0.112489783334934) q[4];
u3(3.04800268996481,-1.24501967312780,1.88447265791333) q[11];
u3(1.51834758520042,0.608198192644516,2.99448028246692) q[8];
cx q[8],q[11];
u1(3.18762276337643) q[11];
u3(-1.41014386752737,0.0,0.0) q[8];
cx q[11],q[8];
u3(2.67843865557780,0.0,0.0) q[8];
cx q[8],q[11];
u3(1.46110017154835,0.0848371242290522,-2.04718472138903) q[11];
u3(2.12924163529042,-0.752629414862580,2.60376944302364) q[8];
u3(2.54747074507232,2.22500289053008,-3.42813316594234) q[3];
u3(0.393191279792165,4.03002225620741,-1.73587713845469) q[7];
cx q[7],q[3];
u1(0.00394317415849588) q[3];
u3(0.561591171390212,0.0,0.0) q[7];
cx q[3],q[7];
u3(3.54009708193780,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.66329330866620,-1.99910037661101,2.81492752312075) q[3];
u3(1.38193625567898,0.950374880720654,2.14773445751485) q[7];
u3(1.48444801987036,2.30837349395485,-0.533168331011296) q[2];
u3(1.78138441938343,0.0479698840308971,-4.25057276626415) q[5];
cx q[5],q[2];
u1(2.37621401408964) q[2];
u3(-1.87945275309757,0.0,0.0) q[5];
cx q[2],q[5];
u3(3.12264704651357,0.0,0.0) q[5];
cx q[5],q[2];
u3(1.87974384746159,-0.182180313882822,-3.01783181397815) q[2];
u3(2.30973468558939,-5.96051883180159,0.163055554166747) q[5];
u3(2.21800016146972,-3.09387669372828,2.20457620748911) q[4];
u3(0.925341413809159,3.42552769795300,-1.30539075109041) q[1];
cx q[1],q[4];
u1(2.74807863415134) q[4];
u3(-3.03199892819877,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.30578112651771,0.0,0.0) q[1];
cx q[1],q[4];
u3(0.703726152050842,-1.02245277412471,0.999069744813233) q[4];
u3(2.30591228606275,3.04310292386039,1.58861783119824) q[1];
u3(1.39744033900098,0.217367875101662,1.35223189218288) q[8];
u3(0.570644583018113,-1.94289570267991,-2.47543794246056) q[9];
cx q[9],q[8];
u1(3.38285917375208) q[8];
u3(-0.923358593287530,0.0,0.0) q[9];
cx q[8],q[9];
u3(1.83786582532608,0.0,0.0) q[9];
cx q[9],q[8];
u3(0.844848623116837,3.30368428591663,0.0626921292439580) q[8];
u3(2.22949445093617,0.522490898763657,4.28229503948561) q[9];
u3(3.02776847547123,-0.435624838886343,-0.150983985602947) q[5];
u3(1.07197063358614,-0.436299063198455,-3.52829691297045) q[3];
cx q[3],q[5];
u1(-0.260763951879929) q[5];
u3(-2.15287651639080,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.24770699604420,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.35646522098612,0.984328199364963,1.09617221369603) q[5];
u3(1.61089365399544,2.88343024975489,0.618283626428187) q[3];
u3(1.70282229703784,-2.35799572778965,-0.175827181931786) q[10];
u3(1.39238839650024,-3.92363287581842,1.25425340449005) q[2];
cx q[2],q[10];
u1(1.68202192949191) q[10];
u3(-2.45941130337601,0.0,0.0) q[2];
cx q[10],q[2];
u3(0.289392318751382,0.0,0.0) q[2];
cx q[2],q[10];
u3(0.218256786607824,2.80312726657592,-1.81352645046051) q[10];
u3(1.15774650896509,-1.21277992518197,0.300472720618522) q[2];
u3(2.65516984250320,1.98943986999487,0.986704989683370) q[11];
u3(1.15188693749417,-2.05870138384504,-2.26756424594590) q[0];
cx q[0],q[11];
u1(0.321426080892355) q[11];
u3(-1.05602729597082,0.0,0.0) q[0];
cx q[11],q[0];
u3(1.66682799640503,0.0,0.0) q[0];
cx q[0],q[11];
u3(0.855551817374629,0.450526958441206,-0.871306587457588) q[11];
u3(0.637262460897440,-1.15696179870953,-2.22673496847286) q[0];
u3(1.71449253977891,-0.0280671548394386,1.27268646811570) q[7];
u3(1.45348340633131,-0.324487724705985,-1.18750575579790) q[6];
cx q[6],q[7];
u1(0.856455151044955) q[7];
u3(-3.34670387839885,0.0,0.0) q[6];
cx q[7],q[6];
u3(2.13024613633402,0.0,0.0) q[6];
cx q[6],q[7];
u3(1.71194907405037,0.0863896504284506,1.67141379054226) q[7];
u3(2.31232305732910,-0.941628613595709,-0.490745639113246) q[6];
u3(0.614851520976713,2.20309427337310,-3.45703310594779) q[6];
u3(1.36356745440657,-3.72063503789569,2.14892164589300) q[4];
cx q[4],q[6];
u1(1.33389821852273) q[6];
u3(0.112714530272987,0.0,0.0) q[4];
cx q[6],q[4];
u3(0.504570032657237,0.0,0.0) q[4];
cx q[4],q[6];
u3(0.557960122069804,1.43980299914244,0.347882214951378) q[6];
u3(1.26925378470258,0.362746741269759,-5.68175557585282) q[4];
u3(2.48837385465018,1.86714165766181,-0.510252791928393) q[8];
u3(0.833442138011699,0.0662348352544389,-2.49329751387992) q[10];
cx q[10],q[8];
u1(1.67392001013369) q[8];
u3(-0.294282786440620,0.0,0.0) q[10];
cx q[8],q[10];
u3(2.62207922514016,0.0,0.0) q[10];
cx q[10],q[8];
u3(2.84420848206569,1.48495989598054,0.896342001752365) q[8];
u3(0.867071796432438,0.499564160496312,4.49356742779521) q[10];
u3(2.47760327603416,3.17571447017897,-2.01341488134822) q[1];
u3(1.44838606721364,2.68327367547830,-3.00815104136408) q[3];
cx q[3],q[1];
u1(3.18577973667527) q[1];
u3(-1.09999262389158,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.95255742891452,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.896511252436232,-2.73857041650172,2.49087550199727) q[1];
u3(1.72505157660740,-2.87344510217644,1.90355217380072) q[3];
u3(1.31321343531515,2.77102042869385,-0.239496305499762) q[5];
u3(1.16030441894816,0.808386188361222,-1.13003484628761) q[2];
cx q[2],q[5];
u1(1.21534690415207) q[5];
u3(-0.982975083132576,0.0,0.0) q[2];
cx q[5],q[2];
u3(0.501857351045218,0.0,0.0) q[2];
cx q[2],q[5];
u3(0.811414481319877,-0.393677067366495,-0.949962804553090) q[5];
u3(1.40458112436043,2.29039634910450,2.40218293552569) q[2];
u3(1.07960222412112,1.46106016156200,-0.774787717755969) q[7];
u3(0.328569275914065,-3.12473479296987,1.42560495854070) q[9];
cx q[9],q[7];
u1(0.224620962940242) q[7];
u3(-1.42500179068784,0.0,0.0) q[9];
cx q[7],q[9];
u3(2.17356802897228,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.03901831329033,3.31696149524510,-1.09514671469624) q[7];
u3(2.24103818847065,1.22509317129287,-2.57819557346888) q[9];
u3(1.25436442957607,-1.90059604371615,-0.510753860031838) q[0];
u3(2.08491022467865,-3.76036745206168,-0.621190154664270) q[11];
cx q[11],q[0];
u1(1.83028086594137) q[0];
u3(-2.69058894372054,0.0,0.0) q[11];
cx q[0],q[11];
u3(0.763440755456916,0.0,0.0) q[11];
cx q[11],q[0];
u3(1.31700406713039,2.70086902866053,0.241192064563322) q[0];
u3(2.01561182380445,3.15877626598017,0.748790412696808) q[11];
u3(2.45684642307010,1.72192024361470,0.653203573787580) q[1];
u3(1.77993998695946,0.217478969561617,-2.86507447436648) q[10];
cx q[10],q[1];
u1(2.44112523542680) q[1];
u3(-1.80831207104484,0.0,0.0) q[10];
cx q[1],q[10];
u3(0.0938518170161691,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.561829933737239,-1.79348432613089,0.821657808068793) q[1];
u3(2.42156623875358,-3.12360989463133,-2.35396565460717) q[10];
u3(0.894579133633724,2.82104713060041,-3.43283558597908) q[6];
u3(1.27685986553247,2.76270083344871,-3.05606721566838) q[0];
cx q[0],q[6];
u1(-0.0264622734120081) q[6];
u3(-0.956600011869451,0.0,0.0) q[0];
cx q[6],q[0];
u3(2.79681204573166,0.0,0.0) q[0];
cx q[0],q[6];
u3(2.02647563241696,-3.81466938490653,1.76263582855467) q[6];
u3(2.45986676461451,1.32093987314469,0.613898694203152) q[0];
u3(1.44549131441659,0.536110486685769,1.15470721280660) q[7];
u3(1.53490778511788,-1.21537250784157,-2.44373757274542) q[5];
cx q[5],q[7];
u1(1.26117560891500) q[7];
u3(-0.998023400109159,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.03008656100212,0.0,0.0) q[5];
cx q[5],q[7];
u3(0.824580934015100,1.44899720476053,-0.742645793214577) q[7];
u3(0.809759836364284,-1.40009736987461,4.45433925811705) q[5];
u3(1.33599869386875,0.308334634665866,2.67457702983801) q[4];
u3(1.69989797380202,-2.78219137676687,-2.36662524107212) q[2];
cx q[2],q[4];
u1(1.90016226872744) q[4];
u3(-1.65880845631737,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.90301976040774,0.0,0.0) q[2];
cx q[2],q[4];
u3(0.643938219532014,0.117555108626085,-0.929486351344427) q[4];
u3(2.63655009960805,-2.42454919229191,2.66111005618667) q[2];
u3(1.39100337638029,-2.01869492465702,3.03279148072938) q[11];
u3(0.593394299271776,1.49488492878522,-0.940384420083168) q[3];
cx q[3],q[11];
u1(2.62822633889206) q[11];
u3(-1.70007433743449,0.0,0.0) q[3];
cx q[11],q[3];
u3(1.24417157103313,0.0,0.0) q[3];
cx q[3],q[11];
u3(1.59737839578995,2.06975232650533,-3.04916466734918) q[11];
u3(0.530898409356436,3.05793143326723,-1.48777374109149) q[3];
u3(2.46181526114088,-2.23726173405810,0.248602608458089) q[9];
u3(1.74031151397159,-2.66891243087868,0.519694479808750) q[8];
cx q[8],q[9];
u1(0.773778374736512) q[9];
u3(-1.03750503287452,0.0,0.0) q[8];
cx q[9],q[8];
u3(-0.0342767558049815,0.0,0.0) q[8];
cx q[8],q[9];
u3(2.10855860126742,0.102134664260327,-4.06271678729033) q[9];
u3(2.21932359027761,-2.74841436008297,0.0352631145932929) q[8];
u3(2.27632309770014,2.44917184189134,-2.52760853673816) q[4];
u3(1.58189192812858,2.99011012697860,-3.17630287509953) q[8];
cx q[8],q[4];
u1(1.64893410188273) q[4];
u3(-2.48097685690864,0.0,0.0) q[8];
cx q[4],q[8];
u3(3.53947878121356,0.0,0.0) q[8];
cx q[8],q[4];
u3(1.28877892222046,-2.82426517421606,1.20807397585024) q[4];
u3(1.65169181773374,2.36792424157283,0.194220574240292) q[8];
u3(1.66387785582337,1.32551291647992,-1.67313453993940) q[2];
u3(1.44652634134536,1.56854408644399,-4.49540606196921) q[3];
cx q[3],q[2];
u1(1.74085293362807) q[2];
u3(0.122094264705530,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.14676881176030,0.0,0.0) q[3];
cx q[3],q[2];
u3(2.75417960333868,-2.71793081949670,2.23708757285659) q[2];
u3(2.10166776821736,-2.56326076658643,0.714146730533596) q[3];
u3(0.633098450025652,2.45543067502850,-1.84683599792823) q[6];
u3(0.182214225746401,2.10912636674936,-2.74885464297403) q[7];
cx q[7],q[6];
u1(2.80195345647289) q[6];
u3(-1.84551536066844,0.0,0.0) q[7];
cx q[6],q[7];
u3(0.673893793412701,0.0,0.0) q[7];
cx q[7],q[6];
u3(1.32527088849669,1.37345842657545,-3.15108934117951) q[6];
u3(1.64403129327536,3.37398515012107,-2.69173314984088) q[7];
u3(1.45643063424158,2.55948519357326,-2.32812316351607) q[1];
u3(0.432334902246624,1.60378099146171,-2.31722233094917) q[9];
cx q[9],q[1];
u1(-1.16709778282785) q[1];
u3(0.437211269482752,0.0,0.0) q[9];
cx q[1],q[9];
u3(3.31680881434109,0.0,0.0) q[9];
cx q[9],q[1];
u3(2.42249954290137,-0.755513310505903,2.65871359779209) q[1];
u3(2.49337336937434,1.46351059147888,3.52086894161630) q[9];
u3(2.39066076127951,0.154620309249912,1.97803884570125) q[5];
u3(2.77524235933891,-3.25880743255836,-1.68293335846655) q[11];
cx q[11],q[5];
u1(0.0934375945189294) q[5];
u3(-1.85403406682832,0.0,0.0) q[11];
cx q[5],q[11];
u3(2.27288570485368,0.0,0.0) q[11];
cx q[11],q[5];
u3(1.65186702922156,-1.45009070390284,1.29023972425514) q[5];
u3(1.43056081931096,-2.71283666575318,1.32319216770977) q[11];
u3(0.379149365482490,1.72856728009779,-2.08712502669335) q[0];
u3(0.988826154780736,-0.796717948940820,-1.74268619779669) q[10];
cx q[10],q[0];
u1(-1.05234107231545) q[0];
u3(0.902419110984252,0.0,0.0) q[10];
cx q[0],q[10];
u3(3.82749444432839,0.0,0.0) q[10];
cx q[10],q[0];
u3(1.34386940835700,-1.16331648595729,2.83366301709351) q[0];
u3(1.36134533833590,4.23460612536413,-1.95662324409846) q[10];
u3(1.89423818423857,0.470791984802942,-3.58721613430828) q[9];
u3(2.03349588109561,-2.94431061458453,3.11943357690632) q[3];
cx q[3],q[9];
u1(0.953578965002053) q[9];
u3(-3.01417251753670,0.0,0.0) q[3];
cx q[9],q[3];
u3(1.92860871893204,0.0,0.0) q[3];
cx q[3],q[9];
u3(2.42591418684784,-0.622941158279353,1.09285702547688) q[9];
u3(0.730758467667708,3.87295246515994,-1.36566311406892) q[3];
u3(0.890369601423157,2.12449197639575,-2.93455024277128) q[6];
u3(1.81973871352586,2.87533325392535,-3.19950112989597) q[4];
cx q[4],q[6];
u1(1.66158725202750) q[6];
u3(-0.150449025669776,0.0,0.0) q[4];
cx q[6],q[4];
u3(2.58426850143046,0.0,0.0) q[4];
cx q[4],q[6];
u3(2.36672404950642,2.54322490290540,-2.41204725069207) q[6];
u3(1.69999192127312,-0.0452895775795441,-0.326482019021012) q[4];
u3(1.82638824638073,0.100452495266747,1.66220534823557) q[11];
u3(1.34591681403768,-1.55509365288734,-2.16800933293322) q[7];
cx q[7],q[11];
u1(2.91585867103025) q[11];
u3(-1.31759501434605,0.0,0.0) q[7];
cx q[11],q[7];
u3(1.89733173777624,0.0,0.0) q[7];
cx q[7],q[11];
u3(1.15776491751082,-2.79621280258093,-1.51780411073589) q[11];
u3(1.39636664544001,4.93318544792093,-0.302949825756473) q[7];
u3(0.695491535560640,0.0465822030949362,-1.02103628867402) q[0];
u3(0.533771945296910,-2.54760598946828,0.338327503314267) q[10];
cx q[10],q[0];
u1(1.98179169922665) q[0];
u3(-2.70725196711042,0.0,0.0) q[10];
cx q[0],q[10];
u3(0.731781901705530,0.0,0.0) q[10];
cx q[10],q[0];
u3(2.07274697160469,-1.50918442966375,1.20872580203029) q[0];
u3(0.953190260406908,1.24195801101686,4.15113012770962) q[10];
u3(2.31081372337043,-2.02511372537227,3.49522501531149) q[5];
u3(0.648614683651655,-1.08904017680021,2.62632565409377) q[1];
cx q[1],q[5];
u1(1.47940731800188) q[5];
u3(-0.765854252726029,0.0,0.0) q[1];
cx q[5],q[1];
u3(-0.376793328180628,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.63865216534869,-1.58109679997673,3.27753790918700) q[5];
u3(2.34711681030489,-4.07317755716720,-0.562857905953613) q[1];
u3(1.69777745448332,1.54649893012012,-3.41319180786780) q[8];
u3(2.39283917679774,-2.81181691065013,3.18732791154460) q[2];
cx q[2],q[8];
u1(0.0983330286035999) q[8];
u3(-1.30843384570676,0.0,0.0) q[2];
cx q[8],q[2];
u3(0.302084538622558,0.0,0.0) q[2];
cx q[2],q[8];
u3(1.92207628218845,-1.79159949110179,-0.0619962618834581) q[8];
u3(1.24736978698240,-3.45289208642510,-1.90449846654590) q[2];
u3(1.19306926415512,0.599158807065151,0.690497110132090) q[8];
u3(1.76866434122480,-0.315734267319526,-2.90163813658242) q[10];
cx q[10],q[8];
u1(3.62413188823612) q[8];
u3(-0.929402708612360,0.0,0.0) q[10];
cx q[8],q[10];
u3(1.83186501931024,0.0,0.0) q[10];
cx q[10],q[8];
u3(0.228164559287038,-0.875816191762240,-1.16465643032112) q[8];
u3(2.19410669989994,0.611127184420915,4.84215083803895) q[10];
u3(0.830987817754976,-0.191699661712548,-0.732296515875852) q[11];
u3(1.56073475329110,-3.63758940159002,0.946572796947676) q[1];
cx q[1],q[11];
u1(0.594883234896979) q[11];
u3(-0.262010800950045,0.0,0.0) q[1];
cx q[11],q[1];
u3(1.82901557253837,0.0,0.0) q[1];
cx q[1],q[11];
u3(0.913131467035099,4.69065046946779,-1.33773178988522) q[11];
u3(0.930696950297370,-5.23184965537609,-0.182291321835333) q[1];
u3(0.482461043054760,-3.04149289742682,2.84064358678822) q[3];
u3(1.14774227615109,-4.33071575433181,1.52638863420101) q[0];
cx q[0],q[3];
u1(1.73485192102510) q[3];
u3(-2.29592973366127,0.0,0.0) q[0];
cx q[3],q[0];
u3(3.12035380415096,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.69137736328069,-0.721943534053983,-0.513088941372821) q[3];
u3(1.24779582222701,-2.22583532475412,-2.96978712728476) q[0];
u3(1.99858659204316,1.13728563905232,-2.62485530324252) q[7];
u3(1.93701218148656,-2.53317459795753,2.29333501512419) q[5];
cx q[5],q[7];
u1(-0.923856229173119) q[7];
u3(0.343094925428769,0.0,0.0) q[5];
cx q[7],q[5];
u3(3.84608732916311,0.0,0.0) q[5];
cx q[5],q[7];
u3(2.14143348304987,1.91273231969170,-2.94850942847272) q[7];
u3(1.60690236297893,-5.44544976536056,-0.665296933746264) q[5];
u3(1.91302442773734,2.68258831766647,-3.24569095013692) q[9];
u3(0.222821916863586,2.98360371716050,-1.82724593256053) q[6];
cx q[6],q[9];
u1(2.97722048749220) q[9];
u3(-2.36957950390743,0.0,0.0) q[6];
cx q[9],q[6];
u3(1.19942675040572,0.0,0.0) q[6];
cx q[6],q[9];
u3(1.25163459997435,-0.549041749793141,0.484573446925306) q[9];
u3(1.94156489818907,-3.81125557129779,2.31887392062761) q[6];
u3(0.236209001770744,2.66245377666293,-3.57111893839945) q[4];
u3(0.619930237551206,-3.20431433222217,1.64796830895254) q[2];
cx q[2],q[4];
u1(1.32396599850755) q[4];
u3(-0.640204805205755,0.0,0.0) q[2];
cx q[4],q[2];
u3(3.13474291539634,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.58689447307119,-2.51823051050591,2.05829652722120) q[4];
u3(1.80691876837164,3.44611992891948,1.28301569646180) q[2];
u3(1.28753781668257,-1.44622463016974,2.16869003658676) q[2];
u3(0.569525106067444,1.44204507186688,-2.62906422060681) q[8];
cx q[8],q[2];
u1(2.57352551413899) q[2];
u3(-2.07432975617326,0.0,0.0) q[8];
cx q[2],q[8];
u3(-0.187318772676152,0.0,0.0) q[8];
cx q[8],q[2];
u3(0.415731813976075,0.277682492993428,2.96020959433070) q[2];
u3(2.32856332696670,-0.473146572260705,-2.09859550607106) q[8];
u3(1.58868783782271,2.75446430112357,-2.62935712787717) q[7];
u3(1.04548952857840,3.06715427250109,-3.14512944628464) q[0];
cx q[0],q[7];
u1(1.31397928612142) q[7];
u3(-3.42965811128725,0.0,0.0) q[0];
cx q[7],q[0];
u3(2.30807637814759,0.0,0.0) q[0];
cx q[0],q[7];
u3(1.17012653571605,3.42289187274269,0.124958480576449) q[7];
u3(1.35709387548976,-2.29382916845492,-2.90357974325142) q[0];
u3(1.61151533352984,-2.60162442469176,-0.436283964431007) q[9];
u3(1.80517513650308,-3.32251398522887,-1.19959908956248) q[5];
cx q[5],q[9];
u1(0.717419937225072) q[9];
u3(-1.07344837738294,0.0,0.0) q[5];
cx q[9],q[5];
u3(-0.0274930642474667,0.0,0.0) q[5];
cx q[5],q[9];
u3(2.30056366232296,-2.32530046714229,2.94317726045862) q[9];
u3(0.870661193870542,-0.956986042651999,1.26817905176927) q[5];
u3(0.711531295308972,3.46193672060951,-1.66279086395431) q[4];
u3(1.31949133498281,1.61814459159528,-3.03113097934595) q[11];
cx q[11],q[4];
u1(2.97658443784956) q[4];
u3(-2.21173155861385,0.0,0.0) q[11];
cx q[4],q[11];
u3(0.559626860895540,0.0,0.0) q[11];
cx q[11],q[4];
u3(1.46762660390450,-1.96104065496650,3.32444669335103) q[4];
u3(1.95818255414311,-4.31305518668667,-1.96580290622869) q[11];
u3(2.97328888357960,0.487621859533777,-3.07961886560480) q[6];
u3(2.68575758712269,3.48599468617629,-0.466745671737370) q[3];
cx q[3],q[6];
u1(3.66723454210806) q[6];
u3(-1.55287045943804,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.41356604294547,0.0,0.0) q[3];
cx q[3],q[6];
u3(2.01256275957106,1.74344339008551,-0.519832387193043) q[6];
u3(0.984529956273987,-2.63182750353721,2.38483277160753) q[3];
u3(0.282389095717857,2.43141385882457,-3.38043732353524) q[1];
u3(0.503877072163199,0.735188353745720,-2.07262816999568) q[10];
cx q[10],q[1];
u1(3.11517593698951) q[1];
u3(-0.898577419144668,0.0,0.0) q[10];
cx q[1],q[10];
u3(1.66102028064005,0.0,0.0) q[10];
cx q[10],q[1];
u3(2.40597946746762,0.0623412041131625,0.809742253074070) q[1];
u3(0.524388561749852,-1.68306119593634,-4.43686118544957) q[10];
u3(1.40140249479353,0.758296900953181,-1.26129880779928) q[8];
u3(1.10191056078106,-4.13802913435900,1.32811444303244) q[0];
cx q[0],q[8];
u1(1.62048245638130) q[8];
u3(-2.14984442153376,0.0,0.0) q[0];
cx q[8],q[0];
u3(-0.0373589430149912,0.0,0.0) q[0];
cx q[0],q[8];
u3(0.924749463529384,2.86734465192810,-2.83655422846703) q[8];
u3(0.564021124934739,-0.488723970816544,1.84759228450873) q[0];
u3(2.20535174901544,0.233889905762628,2.34033792809776) q[9];
u3(2.75264130199275,-2.25905632169760,-0.939049252033115) q[2];
cx q[2],q[9];
u1(-1.06784741615460) q[9];
u3(0.314464122553582,0.0,0.0) q[2];
cx q[9],q[2];
u3(3.01388285123470,0.0,0.0) q[2];
cx q[2],q[9];
u3(0.858629422337729,2.27459473076679,-2.66617203001891) q[9];
u3(1.49992977435676,-3.78752710412279,-2.23341333399801) q[2];
u3(1.03748398991812,1.20008789565477,-0.330750714929468) q[6];
u3(0.845442434000373,-0.344174435553156,-3.25013315883410) q[11];
cx q[11],q[6];
u1(2.66937485289080) q[6];
u3(-2.23594245660029,0.0,0.0) q[11];
cx q[6],q[11];
u3(1.48596128082775,0.0,0.0) q[11];
cx q[11],q[6];
u3(0.529677139433289,-1.48567242919654,-1.74476089002443) q[6];
u3(0.778311765468438,-0.492589797506832,3.79356810987294) q[11];
u3(2.83429865618493,-1.12160608384520,-1.10558747471636) q[1];
u3(0.853102013076600,-1.21581259587508,-3.79584993306705) q[10];
cx q[10],q[1];
u1(2.12921476117283) q[1];
u3(-2.42028234608163,0.0,0.0) q[10];
cx q[1],q[10];
u3(3.04010606223015,0.0,0.0) q[10];
cx q[10],q[1];
u3(0.928728247881901,0.725517375830576,1.98037786872780) q[1];
u3(2.66773156546363,-3.39255068928998,1.12619862467639) q[10];
u3(0.909427696898499,0.0778428764480066,2.40030061642909) q[3];
u3(2.08439406015809,-2.42755553424091,-1.74680691468369) q[7];
cx q[7],q[3];
u1(-0.250106856596861) q[3];
u3(-2.17055626958283,0.0,0.0) q[7];
cx q[3],q[7];
u3(1.08915121386221,0.0,0.0) q[7];
cx q[7],q[3];
u3(1.24437668649581,1.42312834492949,-3.28048323127070) q[3];
u3(3.12728094770357,4.77556029305420,0.0295831454522744) q[7];
u3(2.22870572622643,2.95132634825735,-2.13546613344664) q[5];
u3(1.75341655566393,2.89941748567801,-3.34832222100996) q[4];
cx q[4],q[5];
u1(2.00208390371515) q[5];
u3(0.126470384954316,0.0,0.0) q[4];
cx q[5],q[4];
u3(1.58669293395357,0.0,0.0) q[4];
cx q[4],q[5];
u3(0.482477405553544,2.51373824011771,-2.11191170534946) q[5];
u3(1.47966392079446,-3.69663933048448,1.69408834549288) q[4];
u3(1.84545058157090,1.71988290774074,-2.78423588606669) q[7];
u3(1.54363880026065,1.96294951781749,-2.78942447673826) q[9];
cx q[9],q[7];
u1(2.37991376927739) q[7];
u3(-2.62402284553040,0.0,0.0) q[9];
cx q[7],q[9];
u3(1.08038510310015,0.0,0.0) q[9];
cx q[9],q[7];
u3(1.13265336438668,-2.89487313344463,0.422917648538079) q[7];
u3(1.59527815181521,0.0164441375019886,-1.55219410525835) q[9];
u3(1.06834985956850,-0.783983442500792,1.31482411367484) q[1];
u3(0.320291548902996,0.0948870445409319,-1.72018517022911) q[0];
cx q[0],q[1];
u1(1.60191672487623) q[1];
u3(-0.389580687057327,0.0,0.0) q[0];
cx q[1],q[0];
u3(3.16490871947636,0.0,0.0) q[0];
cx q[0],q[1];
u3(0.976427214200267,-2.57626367050817,3.37971537909545) q[1];
u3(1.58577906390655,0.257849385547571,5.23117810149559) q[0];
u3(2.60995377053472,0.0745089860115837,2.08363954771399) q[11];
u3(2.62674629287296,-2.71824037394273,-1.55411275434043) q[10];
cx q[10],q[11];
u1(0.865144697452175) q[11];
u3(-1.40595266619935,0.0,0.0) q[10];
cx q[11],q[10];
u3(-0.193879689018449,0.0,0.0) q[10];
cx q[10],q[11];
u3(2.45229437661345,3.26471096711061,-2.44938000042006) q[11];
u3(1.41466961158549,-0.0515096789955451,-2.26032417941746) q[10];
u3(1.10923841466375,-3.07534653088181,1.05870930473068) q[4];
u3(2.04943689500817,0.0610235320202377,5.47799829855292) q[6];
cx q[6],q[4];
u1(2.61295125147016) q[4];
u3(-2.15774205021004,0.0,0.0) q[6];
cx q[4],q[6];
u3(1.42695155764101,0.0,0.0) q[6];
cx q[6],q[4];
u3(1.29999029411451,-0.352636125045893,-0.0666439617134308) q[4];
u3(1.22098029670097,1.11968088675919,4.93291305762713) q[6];
u3(2.79353895317965,-2.55961325572389,3.56305597410038) q[5];
u3(0.591674964302616,0.966325852237580,0.778078059171047) q[3];
cx q[3],q[5];
u1(-0.0524614157553001) q[5];
u3(-1.25885724893649,0.0,0.0) q[3];
cx q[5],q[3];
u3(2.08924043689408,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.69612403535893,2.03966557479931,-0.408532905511228) q[5];
u3(1.53944364693268,-0.161918810239373,-2.58670346428494) q[3];
u3(2.54793358627351,-1.59167699966129,1.87932406529214) q[2];
u3(2.82344126613530,-1.90467889866530,0.118567166560618) q[8];
cx q[8],q[2];
u1(0.0891194355103062) q[2];
u3(-1.05321374169224,0.0,0.0) q[8];
cx q[2],q[8];
u3(2.34883685625744,0.0,0.0) q[8];
cx q[8],q[2];
u3(2.63977270836240,-4.30862105998874,1.61489994834413) q[2];
u3(1.22335065879556,0.621846231700438,-4.72568676203996) q[8];
u3(2.97589274862715,-0.445839397887211,2.96512925249460) q[10];
u3(2.35956997442337,0.574248276828541,2.91595493318231) q[6];
cx q[6],q[10];
u1(3.61717503994352) q[10];
u3(-1.50685310668523,0.0,0.0) q[6];
cx q[10],q[6];
u3(2.59206651638286,0.0,0.0) q[6];
cx q[6],q[10];
u3(0.267842811928068,0.113563528337436,-1.59084149364570) q[10];
u3(2.11860012010861,-2.06792732984231,2.22899607315396) q[6];
u3(2.54322367260401,-3.04593851284654,0.171454731118687) q[8];
u3(2.28969685503770,-4.22136222185034,-1.32522661437450) q[3];
cx q[3],q[8];
u1(0.827111106639658) q[8];
u3(-1.40256418361892,0.0,0.0) q[3];
cx q[8],q[3];
u3(3.08760924688872,0.0,0.0) q[3];
cx q[3],q[8];
u3(2.14884469852002,0.195397701540851,-2.30451840407747) q[8];
u3(1.78725526425154,-0.905396386001813,-1.42705881016317) q[3];
u3(1.55695462017894,0.256752196168667,1.98743404195722) q[7];
u3(0.823287665602619,-0.939568020144062,-1.11896284422023) q[2];
cx q[2],q[7];
u1(2.93901239313123) q[7];
u3(-1.49680267727816,0.0,0.0) q[2];
cx q[7],q[2];
u3(1.12579872536744,0.0,0.0) q[2];
cx q[2],q[7];
u3(2.06533367464194,-2.19534000520970,2.29824254829163) q[7];
u3(1.70195018733134,-1.35295204472469,0.0751508546666256) q[2];
u3(1.63685745298261,0.124534618562106,0.391645379097977) q[9];
u3(0.806485065614398,-1.87562779769433,-1.55771489495294) q[11];
cx q[11],q[9];
u1(0.921950075760310) q[9];
u3(-0.129433651544391,0.0,0.0) q[11];
cx q[9],q[11];
u3(2.08061758669832,0.0,0.0) q[11];
cx q[11],q[9];
u3(1.14505497929824,2.73787766736509,-3.17802567520937) q[9];
u3(0.535402830479156,0.722578674500659,-1.91508467156238) q[11];
u3(2.02761001521815,-0.553497783836115,2.67149753207643) q[0];
u3(2.21383126967467,-0.504750747459388,0.228771193529210) q[4];
cx q[4],q[0];
u1(2.88883236082880) q[0];
u3(-1.86035106160217,0.0,0.0) q[4];
cx q[0],q[4];
u3(0.718714344912748,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.84048638471621,-1.32675184767156,-0.156367614359624) q[0];
u3(0.611516110516570,-1.43264424766276,4.40268117496965) q[4];
u3(1.30690617347397,0.890535169112536,-2.64905234393047) q[5];
u3(2.39320073516048,-3.76669282208127,2.35120176270329) q[1];
cx q[1],q[5];
u1(0.568159603114144) q[5];
u3(-1.30464325604782,0.0,0.0) q[1];
cx q[5],q[1];
u3(0.0776985788577713,0.0,0.0) q[1];
cx q[1],q[5];
u3(2.37191300552803,-0.216579408503288,-0.460979357283717) q[5];
u3(0.870981871005001,0.783227249721585,-1.93329593407863) q[1];
u3(0.717033495965178,0.727013845639429,0.342956343456034) q[7];
u3(1.62158230826145,-1.70630836966217,-2.11230931925473) q[11];
cx q[11],q[7];
u1(3.35951341633507) q[7];
u3(-3.92102337693108,0.0,0.0) q[11];
cx q[7],q[11];
u3(-0.871257309879273,0.0,0.0) q[11];
cx q[11],q[7];
u3(1.80411181464980,-2.34909683035261,2.12466473213648) q[7];
u3(1.93506658576711,3.15544154880496,-1.41155601840584) q[11];
u3(2.62040907964999,-1.52390485652143,4.23119912195842) q[3];
u3(0.633812942944924,-1.70067056886671,3.50806411915522) q[1];
cx q[1],q[3];
u1(3.48711777116693) q[3];
u3(-0.849284034856644,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.51903519033931,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.32335850034966,-4.24255386670901,0.305258478360305) q[3];
u3(0.579630709503820,4.84284703219247,0.437696636154120) q[1];
u3(0.904043279011213,2.24554950868573,-0.953107295413560) q[10];
u3(1.53688262888031,-0.238806352181304,-3.32731803365991) q[8];
cx q[8],q[10];
u1(-0.249689783266699) q[10];
u3(-1.89639284045999,0.0,0.0) q[8];
cx q[10],q[8];
u3(0.596410239882890,0.0,0.0) q[8];
cx q[8],q[10];
u3(1.11001105274284,2.52150305487937,-0.968629849322704) q[10];
u3(2.26848031336134,-4.20496918740209,0.602711298272495) q[8];
u3(1.47843139320990,-2.03446647195620,-0.124257389029125) q[0];
u3(1.65486873124838,-3.99399764411522,-0.914525343841592) q[2];
cx q[2],q[0];
u1(2.18521360291120) q[0];
u3(-3.09125547488127,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.64532088008411,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.62645743975134,-0.00463828990151205,-1.64896048436208) q[0];
u3(1.39977947022168,-4.96776889282140,-0.702124302047894) q[2];
u3(1.35138589583527,1.54039528782263,-0.100017856307178) q[4];
u3(0.867592944308674,-1.01716281159178,-1.99654675622090) q[6];
cx q[6],q[4];
u1(1.64781563024997) q[4];
u3(-0.736114159625719,0.0,0.0) q[6];
cx q[4],q[6];
u3(-0.382200275708040,0.0,0.0) q[6];
cx q[6],q[4];
u3(2.68328136831434,-2.03545785355178,0.0717841808669898) q[4];
u3(1.36620195387562,-1.88884361989887,1.04805073174166) q[6];
u3(2.09478427088000,4.01653658126373,-0.965169631037447) q[5];
u3(2.06220290927815,1.64483578115891,-2.43234805073986) q[9];
cx q[9],q[5];
u1(-0.296732430310848) q[5];
u3(-2.08857541874115,0.0,0.0) q[9];
cx q[5],q[9];
u3(1.39660430849046,0.0,0.0) q[9];
cx q[9],q[5];
u3(1.91869706196207,2.33052156989707,-1.54260868170353) q[5];
u3(1.40423274220118,0.607788681448959,-1.05852281633653) q[9];
u3(1.93181118033493,3.74786383848933,-0.938872648462927) q[10];
u3(2.13507617248574,3.32723939522287,0.498411445356374) q[9];
cx q[9],q[10];
u1(2.40039175762573) q[10];
u3(-1.85585471990903,0.0,0.0) q[9];
cx q[10],q[9];
u3(0.528250716232434,0.0,0.0) q[9];
cx q[9],q[10];
u3(1.95679948071975,3.05047565914918,-0.470325264359143) q[10];
u3(1.20481459278284,1.19897964406393,-1.20945653969254) q[9];
u3(0.903391420501111,1.72484241580039,0.0588494158509476) q[5];
u3(1.20398997138966,-0.342419484444175,-3.73655607080768) q[3];
cx q[3],q[5];
u1(1.77367106653773) q[5];
u3(-2.57388510331349,0.0,0.0) q[3];
cx q[5],q[3];
u3(-0.00137952492464777,0.0,0.0) q[3];
cx q[3],q[5];
u3(1.34784650113322,-0.180307170681345,-0.261974956417681) q[5];
u3(0.453880128718256,2.45636270146001,3.03898076936856) q[3];
u3(0.726905656863307,1.71468645388195,-0.157109915580981) q[1];
u3(1.12720431802835,-0.0271705603602264,-3.15317228739030) q[7];
cx q[7],q[1];
u1(3.18487316944016) q[1];
u3(-1.82958189510441,0.0,0.0) q[7];
cx q[1],q[7];
u3(0.450018780556823,0.0,0.0) q[7];
cx q[7],q[1];
u3(2.02014944098869,1.48578213311041,-1.81298697964176) q[1];
u3(0.226708929495256,0.540906077395830,-1.91548769281499) q[7];
u3(1.67819705337091,-0.150477376713563,1.54619535603082) q[2];
u3(1.56469003164475,-2.40789173893132,-0.785421656086708) q[6];
cx q[6],q[2];
u1(2.45797434837936) q[2];
u3(-1.90896583510157,0.0,0.0) q[6];
cx q[2],q[6];
u3(1.51817323165850,0.0,0.0) q[6];
cx q[6],q[2];
u3(2.05019775672222,0.544920441889348,1.00698083003822) q[2];
u3(1.85914096755868,-0.282380707701166,3.07574820893210) q[6];
u3(2.09455866402232,-0.736972249132543,-0.208407301065180) q[0];
u3(0.698864654948834,-2.24021238667342,-2.02751770365527) q[8];
cx q[8],q[0];
u1(3.24284874734170) q[0];
u3(-4.03022776734473,0.0,0.0) q[8];
cx q[0],q[8];
u3(-0.653696511967268,0.0,0.0) q[8];
cx q[8],q[0];
u3(1.14269103277125,1.54418168104331,-1.43520628188753) q[0];
u3(1.54989672475261,-2.01425806056434,1.79231194632785) q[8];
u3(2.11488551254331,1.60821640135396,0.269628723253859) q[11];
u3(2.54060991767448,0.850865392822470,-3.82757999229858) q[4];
cx q[4],q[11];
u1(1.24524508594116) q[11];
u3(-0.111994402147218,0.0,0.0) q[4];
cx q[11],q[4];
u3(2.17047544042833,0.0,0.0) q[4];
cx q[4],q[11];
u3(0.877528926692109,-2.34503310575549,1.57628890768448) q[11];
u3(0.923179222960185,-3.38885806841315,-0.631082147849202) q[4];
u3(2.15023527120317,-0.193112313452750,-1.06303657933602) q[3];
u3(1.37334820809420,-3.19899408109124,0.800498329664521) q[0];
cx q[0],q[3];
u1(1.68918507631281) q[3];
u3(-0.200732707330399,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.583830916024836,0.0,0.0) q[0];
cx q[0],q[3];
u3(2.39141807981642,-3.22312773055157,-0.610349927668937) q[3];
u3(2.67999335907912,-2.41882072349828,-3.36207724452006) q[0];
u3(1.30326432879154,2.78862233532721,-0.974747513382124) q[1];
u3(0.764010692124340,0.630870675081039,-2.15363561657035) q[9];
cx q[9],q[1];
u1(1.36904575703016) q[1];
u3(-0.621952685173978,0.0,0.0) q[9];
cx q[1],q[9];
u3(-0.0260080963406240,0.0,0.0) q[9];
cx q[9],q[1];
u3(1.54014813333424,3.65156797581883,-1.55297990861625) q[1];
u3(2.74556149864641,1.37147948892881,-3.78679516811925) q[9];
u3(1.96594241052502,0.431369351773549,1.67562618422991) q[10];
u3(0.974715334677636,-2.27470023742786,-3.00285554006225) q[8];
cx q[8],q[10];
u1(2.43824117568102) q[10];
u3(0.164750147267271,0.0,0.0) q[8];
cx q[10],q[8];
u3(1.86658388628020,0.0,0.0) q[8];
cx q[8],q[10];
u3(2.80260682244489,2.43182496258961,-2.28388287538572) q[10];
u3(1.32481492618439,-0.399703218372577,0.817145969987815) q[8];
u3(1.76230506649936,-0.957825128219038,0.762254697180717) q[6];
u3(1.77468801062907,-1.85051051235425,-1.36315874728199) q[7];
cx q[7],q[6];
u1(0.254195125142720) q[6];
u3(-1.58924964575394,0.0,0.0) q[7];
cx q[6],q[7];
u3(1.92531880239952,0.0,0.0) q[7];
cx q[7],q[6];
u3(2.50266949787708,-3.09599397152494,-0.690181643682769) q[6];
u3(2.12695267928308,2.15183625756312,2.60827616682282) q[7];
u3(2.17412629712968,3.42074991669554,-1.08700695968657) q[5];
u3(1.51967169364330,2.15817123284841,-0.666349776276609) q[2];
cx q[2],q[5];
u1(0.0660034352682544) q[5];
u3(-0.795687881239835,0.0,0.0) q[2];
cx q[5],q[2];
u3(2.99189567900917,0.0,0.0) q[2];
cx q[2],q[5];
u3(1.11003914903552,-4.46234798902467,0.631661884652620) q[5];
u3(1.55550201881997,-4.53337463142639,1.29115271916006) q[2];
u3(0.618143190141925,1.85631668461256,-3.41169252498564) q[4];
u3(1.73726966366856,-2.12996788704076,3.12334153894633) q[11];
cx q[11],q[4];
u1(1.64836931630408) q[4];
u3(-0.280695024191060,0.0,0.0) q[11];
cx q[4],q[11];
u3(2.14403017540266,0.0,0.0) q[11];
cx q[11],q[4];
u3(2.70002628430492,1.97971771236212,-3.75239046599865) q[4];
u3(2.46876980662993,1.30705139669635,-4.55565208711160) q[11];
u3(0.514359394540534,-2.70909707059925,2.69046945794254) q[6];
u3(0.575545345628943,0.243117984583427,-2.68109524633423) q[3];
cx q[3],q[6];
u1(1.51011644101325) q[6];
u3(-1.96942224760720,0.0,0.0) q[3];
cx q[6],q[3];
u3(2.30741105200456,0.0,0.0) q[3];
cx q[3],q[6];
u3(1.35765558993958,-2.59764125138097,1.90744211205942) q[6];
u3(2.51033180870175,-0.140285208121831,2.39662830175820) q[3];
u3(1.12704642380519,0.0668265316737502,2.68268710599705) q[9];
u3(1.50620199183724,-1.16064462851457,-1.53070699246870) q[0];
cx q[0],q[9];
u1(3.22819419532565) q[9];
u3(-1.83257653981336,0.0,0.0) q[0];
cx q[9],q[0];
u3(0.941806708203996,0.0,0.0) q[0];
cx q[0],q[9];
u3(2.20852325407601,0.750306507782271,-5.03308346556026) q[9];
u3(0.648526473442627,5.45013207240818,0.591425365119014) q[0];
u3(0.958400315348985,-1.77919982844066,2.08949898817297) q[2];
u3(0.590454615396257,-3.60960290522563,2.62574257027913) q[10];
cx q[10],q[2];
u1(-0.460804673648786) q[2];
u3(-1.59273029633499,0.0,0.0) q[10];
cx q[2],q[10];
u3(1.71594291672132,0.0,0.0) q[10];
cx q[10],q[2];
u3(0.408764077433555,1.05599619410072,-5.17698897842970) q[2];
u3(0.965785897524480,-0.565103697310767,0.389787274400154) q[10];
u3(2.05341385015681,1.07686693195404,-1.90335160037489) q[4];
u3(0.956730040426094,1.42317376276184,-4.35138891336201) q[11];
cx q[11],q[4];
u1(0.661461872107339) q[4];
u3(-0.208132351477347,0.0,0.0) q[11];
cx q[4],q[11];
u3(1.73044972330905,0.0,0.0) q[11];
cx q[11],q[4];
u3(2.19951351090855,-0.266973025542770,-2.74620140157161) q[4];
u3(0.539567763114748,1.54677331676577,2.14818154379422) q[11];
u3(1.92170833441660,-2.20461667229459,-0.250936620716525) q[1];
u3(1.54685340714955,1.17724191592825,4.63621834847275) q[5];
cx q[5],q[1];
u1(4.54440008132384) q[1];
u3(-3.39881457604615,0.0,0.0) q[5];
cx q[1],q[5];
u3(-0.103119126504157,0.0,0.0) q[5];
cx q[5],q[1];
u3(1.72042837601298,0.681257014480459,-1.14088237221047) q[1];
u3(0.0856387644706791,-2.83012691609317,-0.474179641398861) q[5];
u3(2.21597887347835,-1.87495788164319,1.19001380137953) q[7];
u3(2.12380756760398,-2.59989876846774,0.273463115375409) q[8];
cx q[8],q[7];
u1(1.34590466105938) q[7];
u3(-0.677920748231085,0.0,0.0) q[8];
cx q[7],q[8];
u3(2.99041851773131,0.0,0.0) q[8];
cx q[8],q[7];
u3(2.45996922101446,-2.25819911481818,1.74842000745500) q[7];
u3(1.60254721598961,1.74779581344265,3.16982714170582) q[8];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11];
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
measure q[11] -> c[11];
