OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.823343835259078,-2.33546698722900,2.74923367229172) q[4];
u3(0.703458770762332,1.81727136872818,-3.01695713319604) q[3];
cx q[3],q[4];
u1(1.63065751319745) q[4];
u3(-2.99528405326592,0.0,0.0) q[3];
cx q[4],q[3];
u3(2.05416122764484,0.0,0.0) q[3];
cx q[3],q[4];
u3(0.298654816967334,-0.786621851122030,3.62034127431320) q[4];
u3(1.86230633952382,-0.241049707305972,-3.55770489452646) q[3];
u3(1.07814731887179,-1.04340414696324,1.18974747005127) q[1];
u3(1.71915154246244,-3.08634467873100,-0.0116390078696555) q[2];
cx q[2],q[1];
u1(-0.274800034329596) q[1];
u3(-1.99600183325031,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.52412811709607,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.83223195027033,0.325754397759549,-2.12324420784374) q[1];
u3(2.32242139197389,-0.0597350166922492,-4.39898894303820) q[2];
u3(0.798810139685665,0.470020374752968,1.62350364617889) q[3];
u3(1.19828534535090,-0.928917478878814,-0.993140989197317) q[1];
cx q[1],q[3];
u1(0.326982915099806) q[3];
u3(-1.36703122380480,0.0,0.0) q[1];
cx q[3],q[1];
u3(2.67467639806253,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.03232273046559,-1.01727086532014,5.13220028635439) q[3];
u3(1.39038021291087,-0.724410043814480,2.40797472739557) q[1];
u3(1.83520544187627,1.71490038423717,-2.80926662482699) q[4];
u3(2.33638829700488,1.73462511838562,-3.60616857444361) q[0];
cx q[0],q[4];
u1(1.86264306644431) q[4];
u3(0.232455299151509,0.0,0.0) q[0];
cx q[4],q[0];
u3(0.540810146994074,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.77380361019560,-1.94726542014218,0.819824962969829) q[4];
u3(1.59747853841669,-4.02846767876953,1.04962657757232) q[0];
u3(0.705706343828235,-1.07745262098126,1.39833010928061) q[1];
u3(1.26944930063390,-2.47337241161319,0.154667074596812) q[4];
cx q[4],q[1];
u1(1.58210174124408) q[1];
u3(-2.34485328869932,0.0,0.0) q[4];
cx q[1],q[4];
u3(0.200253178982521,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.33692237817335,-3.55642413523500,0.913942617243390) q[1];
u3(0.698457453336268,-0.910575912957971,-4.87370703696113) q[4];
u3(1.54561983511444,0.953800380271762,-2.58027956501291) q[3];
u3(1.86360862882349,1.55180660983535,-4.37856407929393) q[0];
cx q[0],q[3];
u1(1.78583175878887) q[3];
u3(-3.23845287386949,0.0,0.0) q[0];
cx q[3],q[0];
u3(0.695286480058047,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.11946294457994,-1.82730165072837,0.695235144740202) q[3];
u3(1.93904061206071,4.72075470556805,0.124411261196037) q[0];
u3(0.714160415758974,-0.708588250642423,-1.49653676376997) q[4];
u3(1.47626476645522,0.805851976732573,-4.18412343528965) q[2];
cx q[2],q[4];
u1(2.92102422398275) q[4];
u3(-1.86849965410919,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.563542811800197,0.0,0.0) q[2];
cx q[2],q[4];
u3(1.98517788572649,-2.79607297629197,0.231502740033389) q[4];
u3(1.99008841758588,2.55989470017758,3.16920782702026) q[2];
u3(1.22202774760474,0.398123825445680,-1.49208951877006) q[0];
u3(1.53085908132749,1.52459611614229,-4.34904368395033) q[3];
cx q[3],q[0];
u1(2.58685018694100) q[0];
u3(-1.65560869571475,0.0,0.0) q[3];
cx q[0],q[3];
u3(3.47133937824652,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.84511702172550,-1.27635786137860,1.63494234733224) q[0];
u3(1.68829905053890,2.63359552478237,1.00960259948985) q[3];
u3(0.839742325003050,2.08702646873460,-0.873857686027142) q[0];
u3(1.48703492780545,0.646661157638821,-2.23029631604726) q[4];
cx q[4],q[0];
u1(3.97549432263935) q[0];
u3(-4.39725617246699,0.0,0.0) q[4];
cx q[0],q[4];
u3(-0.520715045742563,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.53184432558190,3.02680225766436,-2.10272035196614) q[0];
u3(1.50241259907223,-0.395897144588060,2.95350798037044) q[4];
u3(1.99255154306290,2.43696774850263,-0.359711646643553) q[2];
u3(1.97159021555477,-0.593006870311875,-4.84909945741016) q[1];
cx q[1],q[2];
u1(-0.590137756001764) q[2];
u3(1.24780496749399,0.0,0.0) q[1];
cx q[2],q[1];
u3(4.05120873136505,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.53114366215936,2.41551370625909,-2.78292569152194) q[2];
u3(0.839763640922742,4.36162908695583,1.09401521483021) q[1];
u3(1.03470476772441,-2.87246049214488,3.37204613539121) q[1];
u3(1.41901492942658,0.504030011033410,-0.835356943525533) q[2];
cx q[2],q[1];
u1(2.01574082737162) q[1];
u3(-2.61525249016252,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.393912312922962,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.19164201277658,-0.0620454685881837,1.22202645228852) q[1];
u3(2.32092437294111,4.35080844896887,-1.51365737445745) q[2];
u3(0.158465152390288,0.995733964671950,0.769105820766418) q[0];
u3(1.37419268360343,-0.660747263284797,-3.45330006904890) q[4];
cx q[4],q[0];
u1(2.55761583430378) q[0];
u3(-2.01707458269768,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.48488575652429,0.0,0.0) q[4];
cx q[4],q[0];
u3(2.21849510369599,2.95801912521900,-1.14917699209887) q[0];
u3(0.670606707140257,0.0945522667919274,-6.17979096297919) q[4];
u3(1.90448004762355,-3.66244474006602,1.33602607566703) q[2];
u3(1.82734648160594,0.0996803169704241,3.12378520268243) q[1];
cx q[1],q[2];
u1(2.00970371145098) q[2];
u3(-2.72467465725541,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.950759411052366,0.0,0.0) q[1];
cx q[1],q[2];
u3(1.04949280213005,4.30322154386785,-1.79364747304855) q[2];
u3(2.22050560075194,-1.31670480576040,-4.46230796570091) q[1];
u3(2.05618114958786,0.867222591002942,-2.83346413540977) q[0];
u3(2.17887054200845,-2.50943742644071,3.29382405709252) q[3];
cx q[3],q[0];
u1(0.250031260434891) q[0];
u3(-1.59982663413108,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.22393639228563,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.63349988720722,-0.304385283041717,-1.31412630801769) q[0];
u3(2.00738958852722,1.75898286702687,1.11862121134528) q[3];
u3(2.56240129961484,3.91432185976701,-1.70119713030731) q[3];
u3(0.691974623779454,0.125513314131700,1.27080844310675) q[0];
cx q[0],q[3];
u1(0.635526507470762) q[3];
u3(-1.44267502281086,0.0,0.0) q[0];
cx q[3],q[0];
u3(2.17376333753239,0.0,0.0) q[0];
cx q[0],q[3];
u3(1.33611551737357,4.71555639579279,-0.776918318771729) q[3];
u3(1.19697042147596,1.80555234369391,-1.18295790720112) q[0];
u3(1.20971757694340,0.826859987851851,-1.08263925378839) q[2];
u3(0.611677686113490,-1.16416421317938,-1.30921988014678) q[1];
cx q[1],q[2];
u1(2.99855617426550) q[2];
u3(-1.94939933553770,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.64701402280198,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.46509358672032,2.49958519960660,0.633819861060413) q[2];
u3(2.28435407759102,-0.627811917162105,1.20839179733812) q[1];
u3(1.23795977957913,-2.12920851664592,0.816988565629855) q[2];
u3(1.81739280265019,-3.66938770860985,0.253128267596760) q[4];
cx q[4],q[2];
u1(3.38782511708454) q[2];
u3(-4.51300980599132,0.0,0.0) q[4];
cx q[2],q[4];
u3(-0.110302750223613,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.69527605673668,-0.555536038911666,-1.47154938839549) q[2];
u3(2.19830677789145,-1.69534208141234,3.44489776980079) q[4];
u3(2.06239138732491,2.88662246178017,-3.32977779196062) q[0];
u3(1.25788157348533,3.03781898090207,-3.07514187387396) q[3];
cx q[3],q[0];
u1(1.24319895553150) q[0];
u3(-0.127902249786588,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.63945395735036,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.280707319015118,-3.05223529029872,0.804763775490519) q[0];
u3(3.01593531153613,-3.27803394143918,-0.0305001660250896) q[3];
u3(1.06050259685048,2.67252423211686,-1.20826659601073) q[1];
u3(2.20536909926707,0.225271131192865,-3.64793553269344) q[2];
cx q[2],q[1];
u1(0.712232583483461) q[1];
u3(-1.18545312560425,0.0,0.0) q[2];
cx q[1],q[2];
u3(2.72292256269233,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.45733852590861,0.860302056602093,1.59815421063478) q[1];
u3(1.51149593726877,-1.32977017520523,-0.837106348810269) q[2];
u3(2.72487083004509,2.12142884661047,0.302441673273749) q[4];
u3(1.66740253523599,0.00113734652450037,-3.36772518783723) q[3];
cx q[3],q[4];
u1(1.46843368000058) q[4];
u3(-2.55934643133502,0.0,0.0) q[3];
cx q[4],q[3];
u3(0.329825434507202,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.98525783774613,2.84155736626950,-2.99651015262115) q[4];
u3(0.660448739728409,1.18308977349937,0.995755734996300) q[3];
u3(0.654530234238212,1.41872246062029,-0.820039901874957) q[1];
u3(1.86013717508767,-1.20922995244435,-3.76026856957779) q[3];
cx q[3],q[1];
u1(-0.984841751562784) q[1];
u3(0.798873334507632,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.66649527130429,0.0,0.0) q[3];
cx q[3],q[1];
u3(2.37054013295085,3.96503665135572,-0.348604497191805) q[1];
u3(2.22058976821592,3.42449211982982,2.56694028725116) q[3];
u3(0.750165899998202,-3.49457505244637,2.74722999470196) q[4];
u3(1.58173529425156,-3.63758874531849,2.64153038618083) q[2];
cx q[2],q[4];
u1(1.22421694964816) q[4];
u3(-0.0359059484378643,0.0,0.0) q[2];
cx q[4],q[2];
u3(0.686765564666630,0.0,0.0) q[2];
cx q[2],q[4];
u3(2.01668498156372,-0.0636156755486122,-0.951838524816806) q[4];
u3(1.30158167704882,4.29469870420882,-1.70895420664968) q[2];
u3(1.51584031069122,-1.33817750429893,0.457274400795746) q[2];
u3(1.39774676271632,-3.17842838220873,0.544796914366328) q[4];
cx q[4],q[2];
u1(0.664440967342931) q[2];
u3(-0.260255290715250,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.18897017071038,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.05031340666382,2.73494332116463,-1.72890892806478) q[2];
u3(1.98721291017480,-0.219154232445066,0.358235979577961) q[4];
u3(1.07913548422452,1.81846068154847,-0.113703957384621) q[0];
u3(1.19204005195438,0.504214368154669,-3.94430715818585) q[1];
cx q[1],q[0];
u1(0.0565710975302238) q[0];
u3(-1.37562510344467,0.0,0.0) q[1];
cx q[0],q[1];
u3(2.75923846802508,0.0,0.0) q[1];
cx q[1],q[0];
u3(2.37377330500395,0.470143441291703,-2.45493456965203) q[0];
u3(0.812492510092056,-3.74816499559349,-1.38128589113409) q[1];
u3(0.790634852188342,0.0542881012797021,-0.871638547303002) q[3];
u3(0.533589496287711,-2.29966513185297,0.713747710932651) q[1];
cx q[1],q[3];
u1(3.21887520672508) q[3];
u3(-2.30877870173829,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.51272653385168,0.0,0.0) q[1];
cx q[1],q[3];
u3(1.96304981479128,1.89845477341739,-1.45006268033151) q[3];
u3(0.710578513099152,1.75245461534983,1.54816058711740) q[1];
u3(1.26226019621602,1.74474316102510,-0.0827024527436336) q[2];
u3(2.57847074046539,-0.683903876978612,-4.07763747732016) q[4];
cx q[4],q[2];
u1(2.16414208838061) q[2];
u3(-1.69525829321218,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.0165982647378575,0.0,0.0) q[4];
cx q[4],q[2];
u3(0.930631730537947,-1.78057351914728,2.91620485517253) q[2];
u3(1.50338537710319,0.183149023985162,3.07845885475947) q[4];
u3(1.15260832168125,0.541715704328953,-1.21576250306872) q[0];
u3(2.40307606213791,1.52136361931715,-4.42140645284148) q[3];
cx q[3],q[0];
u1(0.158551697092453) q[0];
u3(-0.906232627342124,0.0,0.0) q[3];
cx q[0],q[3];
u3(2.41701849411497,0.0,0.0) q[3];
cx q[3],q[0];
u3(2.10685083897652,2.48539687860010,-3.42345982414654) q[0];
u3(1.09133194090282,-1.38352137851459,-1.62247333556775) q[3];
u3(1.80253243293310,2.34582596380201,-0.409443307913056) q[1];
u3(2.29207304736523,0.569269805137232,-2.57108088994002) q[2];
cx q[2],q[1];
u1(3.68558460523308) q[1];
u3(-3.47295332318164,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.449172427189520,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.53140181230023,0.889752236264910,-4.87010072844077) q[1];
u3(1.25598695211385,3.66120486176209,-0.676198386610074) q[2];
u3(1.68178477513308,3.05657559118885,-1.62557846903498) q[0];
u3(2.63569935001932,2.69958505100127,-0.219066068339792) q[4];
cx q[4],q[0];
u1(0.693355429333034) q[0];
u3(-0.325473132223479,0.0,0.0) q[4];
cx q[0],q[4];
u3(1.63904397415844,0.0,0.0) q[4];
cx q[4],q[0];
u3(1.72311230291464,2.53824369365843,-3.32101807684175) q[0];
u3(1.81054765827359,-0.345913690117605,1.34183780359593) q[4];
u3(0.763250296219413,1.28742947610673,-2.71297175265320) q[1];
u3(0.702385061144366,2.07888298050613,-3.57760085757856) q[2];
cx q[2],q[1];
u1(1.46719635171736) q[1];
u3(-0.00414920679188402,0.0,0.0) q[2];
cx q[1],q[2];
u3(0.487395173316103,0.0,0.0) q[2];
cx q[2],q[1];
u3(0.939171775317436,-2.30606110992074,2.69325903660581) q[1];
u3(1.04889523102749,0.352616940342289,2.55360767182970) q[2];
u3(0.753969571721889,1.67296182005137,-1.28756621845530) q[4];
u3(0.446351809875012,0.480739751890771,-3.05219891391108) q[0];
cx q[0],q[4];
u1(2.66542748404353) q[4];
u3(-1.66720613236426,0.0,0.0) q[0];
cx q[4],q[0];
u3(3.40075934836623,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.47234973731999,-2.07983647616604,-0.466263664070436) q[4];
u3(2.27160927924850,1.00605631245443,-1.59209303947432) q[0];
u3(1.39757739708858,-0.160227710233255,-0.356012427846129) q[3];
u3(1.56155333114175,-2.83978054713996,0.652690229914410) q[2];
cx q[2],q[3];
u1(4.11209827398991) q[3];
u3(-3.66338791753569,0.0,0.0) q[2];
cx q[3],q[2];
u3(-0.0877497446669289,0.0,0.0) q[2];
cx q[2],q[3];
u3(2.48826852942104,-2.88315024716489,-0.736212926921656) q[3];
u3(1.93757238926719,-0.549951524622580,-5.01981579774529) q[2];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
