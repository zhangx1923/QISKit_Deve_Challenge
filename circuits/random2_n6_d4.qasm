OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[6];
u3(1.02013098385628,-0.393320243787515,2.17987952568103) q[2];
u3(1.00474380186496,-1.77219302635811,-2.08592750304601) q[4];
cx q[4],q[2];
u1(1.49321088299148) q[2];
u3(-2.89148935046288,0.0,0.0) q[4];
cx q[2],q[4];
u3(1.98845187793252,0.0,0.0) q[4];
cx q[4],q[2];
u3(1.10124154096607,-0.670427334733272,-1.16630591081945) q[2];
u3(1.02183379106099,-3.12506309402546,-0.618271421673124) q[4];
u3(1.90480482956290,-2.14046549241495,0.901435096081941) q[5];
u3(2.84761173815709,-2.43517389812689,0.269081759383683) q[0];
cx q[0],q[5];
u1(0.978763628131851) q[5];
u3(-0.551616499348217,0.0,0.0) q[0];
cx q[5],q[0];
u3(2.14671082103109,0.0,0.0) q[0];
cx q[0],q[5];
u3(0.907730727604461,0.272260062069381,-1.17998591543301) q[5];
u3(2.41329162411909,-2.60643487924524,1.23353855035241) q[0];
u3(0.516614590684074,1.82648132978185,-2.11252339539938) q[1];
u3(0.686021867211321,0.248268904107873,-1.05789651815158) q[3];
cx q[3],q[1];
u1(2.25534841680346) q[1];
u3(-2.96741761768326,0.0,0.0) q[3];
cx q[1],q[3];
u3(1.69214552996813,0.0,0.0) q[3];
cx q[3],q[1];
u3(1.91139654747801,-2.88908793150133,-0.374029116770294) q[1];
u3(1.14206901488703,0.654297940569700,-2.34910278641809) q[3];
u3(1.18536592500342,-0.208449571681425,-1.29337513583616) q[5];
u3(1.71749435967943,-3.57356336648048,2.47723045750367) q[1];
cx q[1],q[5];
u1(1.49206508956540) q[5];
u3(-3.60052767775186,0.0,0.0) q[1];
cx q[5],q[1];
u3(2.35160431130266,0.0,0.0) q[1];
cx q[1],q[5];
u3(1.46704616658277,0.880004150422045,2.62323136344918) q[5];
u3(1.32130914947293,-0.230084472940264,3.20806769786751) q[1];
u3(2.67133984317548,0.105920123799514,-2.28614266639468) q[3];
u3(2.18814661774048,4.56341666176179,0.288438349481631) q[4];
cx q[4],q[3];
u1(0.0798382199885992) q[3];
u3(-1.46655037259405,0.0,0.0) q[4];
cx q[3],q[4];
u3(2.78852473732678,0.0,0.0) q[4];
cx q[4],q[3];
u3(1.41898082301251,-1.12062053581306,-1.46916281547314) q[3];
u3(1.66461777637165,-0.826640138324021,-4.90745085948292) q[4];
u3(1.46360707914953,1.74666618058230,-3.31812376396650) q[0];
u3(1.87210057568600,2.12472876888221,-3.62764757162977) q[2];
cx q[2],q[0];
u1(0.921409734957452) q[0];
u3(-0.368742930813271,0.0,0.0) q[2];
cx q[0],q[2];
u3(3.15014783254648,0.0,0.0) q[2];
cx q[2],q[0];
u3(1.69671751593113,1.39263595067204,-3.66508228608237) q[0];
u3(2.18953680923474,-1.55761571888213,-3.27120386698966) q[2];
u3(1.22866511246230,1.24796838740909,-0.0501533253887732) q[0];
u3(1.36444166851335,0.227447111485870,-3.57463694246167) q[2];
cx q[2],q[0];
u1(3.14448681576324) q[0];
u3(-1.82399153215489,0.0,0.0) q[2];
cx q[0],q[2];
u3(0.856433824735067,0.0,0.0) q[2];
cx q[2],q[0];
u3(2.83377871967564,0.177793250678632,-4.04102711556531) q[0];
u3(1.58274613940860,-2.16036946143871,1.64682102391750) q[2];
u3(1.74050331682927,0.512459798525940,-0.445778802684106) q[5];
u3(1.44822985250982,-4.71824175218940,1.33916456264783) q[4];
cx q[4],q[5];
u1(3.67034888515527) q[5];
u3(-1.44195303792614,0.0,0.0) q[4];
cx q[5],q[4];
u3(2.25255405947788,0.0,0.0) q[4];
cx q[4],q[5];
u3(1.55546870632022,2.87517661772717,-0.233865544267406) q[5];
u3(0.129091730039324,0.389587856209819,-3.73240289073329) q[4];
u3(2.13291399761810,-0.509328892641801,1.10833947891442) q[3];
u3(2.21813982563044,-2.00817354188458,-2.27627574137509) q[1];
cx q[1],q[3];
u1(2.56674752758348) q[3];
u3(-2.15620266933666,0.0,0.0) q[1];
cx q[3],q[1];
u3(0.277810416413361,0.0,0.0) q[1];
cx q[1],q[3];
u3(0.869590997584930,0.102305400891768,-0.195958356597473) q[3];
u3(2.41611939746099,-2.43621938723098,0.629662981828665) q[1];
u3(1.18334917823155,0.965412870477379,-3.04201677547752) q[4];
u3(1.87209782481783,-3.21458572120547,2.75611751019513) q[1];
cx q[1],q[4];
u1(-0.407579722813587) q[4];
u3(-2.10854519434817,0.0,0.0) q[1];
cx q[4],q[1];
u3(1.21446057348982,0.0,0.0) q[1];
cx q[1],q[4];
u3(1.32255811037585,1.84273198247264,-1.70753355452961) q[4];
u3(0.776507981495984,1.11720605750475,-2.21987265049422) q[1];
u3(0.137244915152202,-2.06335258138727,2.36975236489651) q[5];
u3(0.228761255734349,-2.74162166631558,1.55270635032719) q[3];
cx q[3],q[5];
u1(0.671480799981783) q[5];
u3(-0.216449493121657,0.0,0.0) q[3];
cx q[5],q[3];
u3(1.87679111082676,0.0,0.0) q[3];
cx q[3],q[5];
u3(2.75451761202216,3.38715741457073,-1.23868131350846) q[5];
u3(1.23693310394336,-0.513849403257303,-2.52209422190924) q[3];
u3(1.43300060442205,1.52276499498014,-3.22585677729000) q[2];
u3(1.77931025967746,2.90164244985054,-3.28015307012536) q[0];
cx q[0],q[2];
u1(1.57565879474913) q[2];
u3(-3.03217070164505,0.0,0.0) q[0];
cx q[2],q[0];
u3(1.12578415986631,0.0,0.0) q[0];
cx q[0],q[2];
u3(3.05583026225875,-1.04128475841282,5.02126929936342) q[2];
u3(2.14021766424955,1.81847342390697,-0.400347584528048) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
