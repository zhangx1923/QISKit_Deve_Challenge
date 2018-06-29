OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
u3(0.983149482210510,-2.91485327525368,2.57102783655643) q[1];
u3(0.119351650183050,-0.779001960050122,-1.79115544378359) q[3];
cx q[3],q[1];
u1(3.01142960308845) q[1];
u3(-2.36066612695593,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.17302493043931,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.22066063311997,1.50109484840639,-3.38344145466667) q[1];
u3(2.01529619943508,-4.68669183602393,-1.19145235746597) q[3];
u3(2.14260617828951,2.20753782310657,-0.0298446169155533) q[2];
u3(2.19255087024594,-0.707898932263724,-4.61530596100684) q[4];
cx q[4],q[2];
u1(1.61276945292684) q[2];
u3(-2.81468812344867,0.0,0.0) q[4];
cx q[2],q[4];
u3(0.698882907122108,0.0,0.0) q[4];
cx q[4],q[2];
u3(2.03853862301445,-1.03550210652728,0.971953235628279) q[2];
u3(0.459319329830280,3.36405939990436,-1.50324755446827) q[4];
u3(1.60452675733054,0.0783569758611526,-2.00471903982437) q[1];
u3(0.970487503725244,-4.04072452657999,2.09449855500211) q[0];
cx q[0],q[1];
u1(2.98918881762036) q[1];
u3(-2.55095977347641,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.26669602983586,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.13077442677057,1.98791418717704,2.05015287970911) q[1];
u3(2.19577252246947,2.54581954251950,0.346841439701515) q[0];
u3(0.931940623080685,-0.0905230846514568,0.382070637920208) q[4];
u3(0.845705349739625,-2.62447196637286,2.40413577316705) q[3];
cx q[3],q[4];
u1(3.18894038590906) q[4];
u3(-0.734080016782100,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.86596328910624,0.0,0.0) q[3];
cx q[3],q[4];
u3(2.07203964412107,2.01766189628454,-1.93822740294461) q[4];
u3(1.65155253464779,-3.84267924279942,-2.26005299419356) q[3];
u3(1.97062957502914,-1.03022462316189,3.82750272810130) q[1];
u3(1.50336076181696,0.999639702483212,1.59430515115625) q[3];
cx q[3],q[1];
u1(2.51230841656677) q[1];
u3(-1.56469963332892,0.0,0.0) q[3];
cx q[1],q[3];
u3(3.42558881728054,0.0,0.0) q[3];
cx q[3],q[1];
u3(0.312295314212147,-1.45276142400705,2.53730756468495) q[1];
u3(1.62882619739435,0.399368395598044,4.15446413583602) q[3];
u3(0.960090251440496,1.53591712516237,-0.977920669641869) q[2];
u3(0.0702003437183778,-0.136437696820608,-0.658805190944056) q[0];
cx q[0],q[2];
u1(2.73529322908377) q[2];
u3(-2.18542190530925,0.0,0.0) q[0];
cx q[2],q[0];
u3(0.0640476506267542,0.0,0.0) q[0];
cx q[0],q[2];
u3(1.92854750143552,-2.20568819879045,3.73778715338955) q[2];
u3(1.80583463931316,-0.913892085445317,1.27977807797978) q[0];
u3(2.30124578615086,0.470170143009390,0.850775946828110) q[3];
u3(2.08174109124857,-1.05343944938831,-1.44433941463152) q[4];
cx q[4],q[3];
u1(3.52574098554905) q[3];
u3(-1.52369909182712,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.48969807690230,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.18663343132455,2.61594303937723,-3.43021489680380) q[3];
u3(1.13661987432599,0.352282411515701,5.25228420227273) q[4];
u3(1.20590585648540,2.17631337843292,-0.606952333466849) q[2];
u3(1.46909294197710,0.327935937536341,-3.87126729863193) q[1];
cx q[1],q[2];
u1(1.57825483362429) q[2];
u3(-2.43835985086464,0.0,0.0) q[1];
cx q[2],q[1];
u3(0.0745320173004924,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.64350792489181,-3.29795007128181,-0.307611748896134) q[2];
u3(1.12691821858434,2.96570196930955,0.0775256947407192) q[1];
u3(1.90723452713104,-0.392167569550463,-0.455281357773098) q[0];
u3(0.915650001622166,-3.46870137824774,0.241700791859274) q[3];
cx q[3],q[0];
u1(3.49726906169625) q[0];
u3(-0.956574642578300,0.0,0.0) q[3];
cx q[0],q[3];
u3(1.77256805900342,0.0,0.0) q[3];
cx q[3],q[0];
u3(0.835270631139108,-2.73131105511450,3.31535633899481) q[0];
u3(1.87051083907764,-5.40606466188993,0.400035022899825) q[3];
u3(2.62860022617165,1.99292314647085,-3.48915678022656) q[1];
u3(1.21451480956377,3.08527667634168,-2.43646506199145) q[4];
cx q[4],q[1];
u1(3.23894405212853) q[1];
u3(-0.657650868750890,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.35507771065482,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.82097571139773,-0.740531855803809,1.17291763773531) q[1];
u3(0.906786204448733,1.88608770998073,-1.06203071880881) q[4];
u3(2.86153613879791,1.97706358828859,-0.791520370760945) q[1];
u3(2.62199883230014,4.56320034664420,-0.840253843650881) q[2];
cx q[2],q[1];
u1(4.33941885787883) q[1];
u3(-3.71919398122826,0.0,0.0) q[2];
cx q[1],q[2];
u3(-0.900083740072720,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.95164468187945,-0.455858601432892,-2.09064938330621) q[1];
u3(1.56615363350527,3.05892826195576,0.499876116964280) q[2];
u3(1.04578756197826,0.349638261113950,1.57681056529828) q[4];
u3(1.41212962140403,-0.851805811279783,-3.01166410887443) q[0];
cx q[0],q[4];
u1(3.70874915150844) q[4];
u3(-1.21940718147669,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.73729995405252,0.0,0.0) q[0];
cx q[0],q[4];
u3(2.79075098875217,-0.900431818448763,0.568611711230498) q[4];
u3(0.993893794414546,-2.81986474809605,1.04483006472719) q[0];
u3(2.23375619411017,-0.467532707970806,0.430547044304480) q[4];
u3(2.12156944004412,-2.83107619818253,-0.500857928485184) q[1];
cx q[1],q[4];
u1(0.960837508073456) q[4];
u3(0.0627075429101140,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.64479306089286,0.0,0.0) q[1];
cx q[1],q[4];
u3(2.22201089694510,2.27800943162803,-1.17477929288098) q[4];
u3(0.482333504242260,-0.901812306885656,-4.51296620913352) q[1];
u3(2.23136423830271,2.85005294043274,-1.97653808821982) q[0];
u3(1.76034059412006,1.13371001976487,-1.74421025944254) q[2];
cx q[2],q[0];
u1(3.45405387988283) q[0];
u3(-1.09022987706060,0.0,0.0) q[2];
cx q[0],q[2];
u3(1.77065195943058,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.10439397788741,1.57311099466202,-0.140739632221341) q[0];
u3(1.62359885287055,-1.82082549057158,-4.06901196978946) q[2];
u3(1.34700912797478,-3.62093000255696,0.885101466374499) q[0];
u3(1.81382563641453,-0.0198968552519381,3.40019508348741) q[1];
cx q[1],q[0];
u1(3.04140305940675) q[0];
u3(-2.31013694582427,0.0,0.0) q[1];
cx q[0],q[1];
u3(0.522028814969639,0.0,0.0) q[1];
cx q[1],q[0];
u3(1.46379689101211,0.885617325140635,-0.430551217996779) q[0];
u3(2.52110279327782,1.50412286695173,-2.15176099012989) q[1];
u3(0.513211286795876,1.36384707532675,-2.07925304216257) q[2];
u3(1.53618902384643,-4.45088934641377,1.82812132928773) q[3];
cx q[3],q[2];
u1(2.28208913080223) q[2];
u3(0.249117731819310,0.0,0.0) q[3];
cx q[2],q[3];
u3(1.36791867501947,0.0,0.0) q[3];
cx q[3],q[2];
u3(1.21557294071419,-2.70760834456451,-0.0215455158346369) q[2];
u3(0.512329891513651,-2.58366066300927,1.06147164639373) q[3];
u3(1.71326243200349,0.0341452954687442,-2.22006139227923) q[4];
u3(1.10894355263670,-3.44659252938247,2.00112391127654) q[3];
cx q[3],q[4];
u1(0.686891745016872) q[4];
u3(-3.15743644748572,0.0,0.0) q[3];
cx q[4],q[3];
u3(1.84769473152725,0.0,0.0) q[3];
cx q[3],q[4];
u3(1.95915520195862,-3.56971441553900,2.22705837454628) q[4];
u3(2.62832116447695,-1.71030314234283,2.15023251737697) q[3];
u3(3.12533497331033,-1.42221041399961,-0.672657619029136) q[1];
u3(1.02418446826085,-4.08142444892926,-1.27673037585792) q[2];
cx q[2],q[1];
u1(3.42174196256046) q[1];
u3(-0.704211358060739,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.93993534350475,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.97654713635315,1.30594458582652,-1.26262554349375) q[1];
u3(2.27024708929369,-2.73481353841452,-1.10170890056607) q[2];
u3(1.57298919608633,-0.537652612908643,2.39582741204667) q[1];
u3(1.26954691075136,-0.567021867972418,-1.53907648652089) q[0];
cx q[0],q[1];
u1(0.0154316548970121) q[1];
u3(-2.58337124387982,0.0,0.0) q[0];
cx q[1],q[0];
u3(1.14256712314989,0.0,0.0) q[0];
cx q[0],q[1];
u3(1.77998186826898,0.974017835301651,-3.91770294401918) q[1];
u3(1.72307263921180,1.42895195588610,1.84216112201377) q[0];
u3(2.50116521403148,2.33697502742108,-3.64094640934704) q[3];
u3(0.372643487351610,3.39057970250650,-1.17668596817748) q[4];
cx q[4],q[3];
u1(1.68260226873580) q[3];
u3(0.445362875304283,0.0,0.0) q[4];
cx q[3],q[4];
u3(0.740937852391233,0.0,0.0) q[4];
cx q[4],q[3];
u3(2.59533112024933,1.64987987785214,-2.64807046123516) q[3];
u3(0.872222437365689,-3.24031285543287,-2.00895114081855) q[4];
u3(0.760683806130042,-0.756264714169072,-0.303806105267735) q[2];
u3(1.97240416888640,-4.98569143337191,1.28744907050783) q[4];
cx q[4],q[2];
u1(0.283320984906488) q[2];
u3(-1.20400392736904,0.0,0.0) q[4];
cx q[2],q[4];
u3(2.16819985023413,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.65421125745080,-0.643160311723952,0.0708587839620041) q[2];
u3(2.43629247196450,-0.559928998203852,2.45901470082441) q[4];
u3(0.376676744546434,-0.782122091109345,1.13231873296622) q[3];
u3(0.558719084469743,-2.38777022647586,-0.268850156430576) q[1];
cx q[1],q[3];
u1(3.16430453605712) q[3];
u3(-1.74109239387713,0.0,0.0) q[1];
cx q[3],q[1];
u3(1.06351259054125,0.0,0.0) q[1];
cx q[1],q[3];
u3(2.20739817461327,-2.15368680600006,3.42623300526614) q[3];
u3(1.41641386495895,-3.97272609027277,-0.649112083334375) q[1];
u3(1.92230490108804,-1.06092550755851,1.45128729870716) q[0];
u3(1.97419687420663,-2.19505765677256,-2.12186307620882) q[2];
cx q[2],q[0];
u1(0.180752879337267) q[0];
u3(-0.950624559210570,0.0,0.0) q[2];
cx q[0],q[2];
u3(2.78780495770557,0.0,0.0) q[2];
cx q[2],q[0];
u3(0.764580290502516,3.82658271774415,-1.69401419410749) q[0];
u3(2.20321604242763,0.659861046784189,-1.42804500503856) q[2];
u3(1.27230076311620,1.26465520164813,-0.389355222096128) q[1];
u3(0.433234201882873,-0.479486446121475,-1.36181176031636) q[4];
cx q[4],q[1];
u1(2.48990856052166) q[1];
u3(-3.09623754722038,0.0,0.0) q[4];
cx q[1],q[4];
u3(1.55135515229414,0.0,0.0) q[4];
cx q[4],q[1];
u3(1.97622436942285,1.24018836571956,-2.88082645078032) q[1];
u3(1.27720944320035,1.27192074130929,-1.60542098157680) q[4];
u3(1.07673582663435,-0.737001301643942,1.51716567305488) q[1];
u3(0.880122124312216,-2.62740744443141,-0.601076479574137) q[2];
cx q[2],q[1];
u1(1.82120889930735) q[1];
u3(-2.61690006968020,0.0,0.0) q[2];
cx q[1],q[2];
u3(1.20323120681692,0.0,0.0) q[2];
cx q[2],q[1];
u3(1.47676953166924,-3.14197649063066,2.28803937958274) q[1];
u3(1.74737072960547,3.54632277571721,1.53706914741132) q[2];
u3(1.66478138283989,-1.12130986606220,-0.281389023430054) q[0];
u3(2.35386894422221,-1.90841688100921,1.19594576646424) q[3];
cx q[3],q[0];
u1(2.40326553301667) q[0];
u3(-1.69359501196799,0.0,0.0) q[3];
cx q[0],q[3];
u3(0.474039801589573,0.0,0.0) q[3];
cx q[3],q[0];
u3(1.76802673185910,-4.59071809652466,0.898186537661036) q[0];
u3(0.805912310456949,-2.59744557199625,2.27969788375589) q[3];
u3(1.47665730571730,1.12060708233847,-2.05758949447219) q[4];
u3(2.11321012444914,-4.74390216324179,1.22661323848422) q[0];
cx q[0],q[4];
u1(2.09588781395184) q[4];
u3(-2.56955155165077,0.0,0.0) q[0];
cx q[4],q[0];
u3(1.26860670473640,0.0,0.0) q[0];
cx q[0],q[4];
u3(1.39789797817904,2.30816042189718,0.989046431006010) q[4];
u3(1.36405566255902,0.506617718633538,2.12644716352749) q[0];
u3(2.39262415484394,1.25427471962920,0.937604984202793) q[2];
u3(1.57823413974479,0.459148369141283,-3.08396202433082) q[1];
cx q[1],q[2];
u1(2.71543177905680) q[2];
u3(-1.97533710766913,0.0,0.0) q[1];
cx q[2],q[1];
u3(1.22098871467418,0.0,0.0) q[1];
cx q[1],q[2];
u3(2.29985000916248,3.41579283558251,-0.729951649205621) q[2];
u3(1.79266811518837,3.92348701778297,0.124369325468332) q[1];
barrier q[0],q[1],q[2],q[3],q[4];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
