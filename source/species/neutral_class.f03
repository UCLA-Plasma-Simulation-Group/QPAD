module neutral_class

use parallel_module
use param
use sysutil_module
use options_class
use fdist2d_class
use part2d_class
use part2d_comm
use field_class
use field_e_class
use field_b_class
use field_psi_class
use field_src_class
use hdf5io_class

implicit none

private

!-----------------------------------------------------------------------------------------
! Ionization Parameters for ADK model

! Data taken from:

! 2010 CODATA Internationally recommended values of the Fundamental Physical Constants

! http://physics.nist.gov/cuu/Constants/
!
! NIST Atomic Spectra Database Ionization Energies Data
! Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2014). NIST Atomic Spectra
! Database (ver. 5.2), [Online]. Available: http://physics.nist.gov/asd [2014, December].

! National Institute of Standards and Technology, Gaithersburg, MD
!
! http://physics.nist.gov/PhysRefData/ASD/ionEnergy.html

double precision, parameter, dimension(3,1) :: H_param = &
    reshape( (/ 8.522542995398661d19, 342.53947239007687d0, 1.0005337056631487d0 /), &
                                                                     (/ 3, 1 /) )
double precision, parameter, dimension(3,3) :: Li_param = &
    reshape( (/ 3.460272990838495d21, 85.51998980232813d0, 2.1770706138013733d0, &
                3.6365138642921554d20, 4493.713340713575d0, 0.6964625952167312d0, &
                2.0659396971422902d22, 9256.32561931876d0, 0.9999745128918196d0 /), &
                                                                     (/ 3, 3 /) )
double precision, parameter, dimension(3,11) :: Na_param = &
    reshape( (/ 4.1003766070018457d21, 79.58018178025945d0, 2.2542264609655485d0, &
                6.646980062859022d21, 2221.166444870671d0, 1.1456178592204318d0, &
                1.4450031005839113d24, 4140.270846854678d0, 1.6151368370485608d0, &
                1.1938031389836218d26, 6722.159287368323d0, 1.966696637087852d0, &
                2.363515645181666d27, 11121.962149645387d0, 2.1353951423226762d0, &
                6.533560672234857d28, 15439.488740548195d0, 2.3727930406073865d0, &
                1.2810219583986103d30, 20565.362378303853d0, 2.5763114425119764d0, &
                7.256788653171509d30, 29333.019535393767d0, 2.630951747911784d0, &
                1.2248965512782968d32, 35469.501316660346d0, 2.8341894069142315d0, &
                3.7691814856879473d24, 383083.54115896183d0, 0.9273100377897379d0, &
                1.3468138246659564d25, 457288.7579154525d0, 0.998535765107404d0 /), &
                                                                     (/ 3, 11 /) )
double precision, parameter, dimension(3,19) :: K_param = &
    reshape( (/ 7.216074209379994d21, 61.77481632066079d0, 2.5408885716705307d0, &
                9.406761452174264d22, 1214.8511771634942d0, 1.6236448747182823d0, &
                6.779724188614074d25, 2117.481312873944d0, 2.2701233926960898d0, &
                1.811619844825418d28, 3247.768300305309d0, 2.7807769758524734d0, &
                1.0119877391508012d30, 5133.580528079887d0, 3.0570711651756124d0, &
                9.049414272968163d31, 6769.504077462043d0, 3.4396463834481468d0, &
                5.0242638226357555d33, 8706.961891869527d0, 3.7627586361816903d0, &
                2.578271439282718d34, 13165.223159585394d0, 3.74238337717886d0, &
                7.967483576694055d35, 15924.655928532302d0, 4.007280316055681d0, &
                3.196309283214214d30, 77214.0625914301d0, 2.28713302745643d0, &
                2.3563335866942047d31, 91884.31431199581d0, 2.4121515077162674d0, &
                1.4984572622004427d32, 108273.28574947864d0, 2.524170311816701d0, &
                5.8647348952001106d32, 130518.49747304033d0, 2.5873086510013654d0, &
                3.0398926235765757d33, 150600.54746568995d0, 2.6832845094542193d0, &
                1.4304663067720226d34, 172552.5863001384d0, 2.7713790342861593d0, &
                3.4929050463473155d34, 205618.1565293413d0, 2.7944506601793733d0, &
                1.85797716283929d35, 227299.68713051578d0, 2.8991091229466286d0, &
                9.046983127409002d25, 2.1387084800213743d6, 0.9555610969933448d0, &
                1.9843676239339244d26, 2.3674584145525303d6, 0.9954563352900228d0 /), &
                                                                     (/ 3, 19 /) )
double precision, parameter, dimension(3,37) :: Rb_param = &
    reshape( (/ 8.142038491103064d21, 58.31683597165995d0, 2.6095364358731663d0, &
                2.5535438019000036d23, 973.8041035058765d0, 1.824374032909187d0, &
                2.7000601553750055d26, 1679.526944428912d0, 2.5327128978496343d0, &
                9.685790844349714d28, 2576.2205018430946d0, 3.0842758054317407d0, &
                1.0748558554029407d31, 3864.21631810689d0, 3.459973470196984d0, &
                1.1603503988608595d33, 5155.954467164836d0, 3.861433025510535d0, &
                7.592811710363505d34, 6695.067652968157d0, 4.198712465440701d0, &
                3.228036166231273d35, 10453.804888994271d0, 4.121312432452735d0, &
                1.2294619485882202d37, 12628.28840362558d0, 4.409748725556149d0, &
                4.594628560441955d34, 31512.226042057184d0, 3.4315574796699924d0, &
                5.099862812557592d35, 37844.42220681386d0, 3.5860783832270737d0, &
                3.66051737982471d36, 45883.04950656357d0, 3.691876466328285d0, &
                2.4532542930217856d37, 54647.11820807218d0, 3.795166356205491d0, &
                1.6954926367367167d38, 63691.70382403349d0, 3.9070049315032813d0, &
                5.741143312291227d38, 76830.35737138857d0, 3.938894193953593d0, &
                3.1772643710662465d39, 88109.210912276d0, 4.033023518460716d0, &
                1.547282384799606d40, 100644.25410629017d0, 4.115663943965691d0, &
                6.970575981163561d40, 114246.76199511335d0, 4.192469483287085d0, &
                3.2341579464837933d41, 128139.88477866953d0, 4.275232097540516d0, &
                7.882013953979101d40, 171375.41143787434d0, 4.03999269039215d0, &
                4.383609296869856d41, 186065.5873110981d0, 4.148886949420316d0, &
                2.0604060864046853d42, 202832.3654692189d0, 4.2411470329012415d0, &
                7.118253469778463d42, 223834.59618026362d0, 4.302347413111842d0, &
                3.050179866504092d43, 242444.82069622068d0, 4.387529805709947d0, &
                1.6607075650654932d44, 257754.38432565072d0, 4.4986245825104865d0, &
                1.7436005002448962d44, 299172.63990271214d0, 4.441460377264627d0, &
                7.906799692579144d44, 318173.5516199869d0, 4.535945678561607d0, &
                3.491233570113727d36, 1.1980646266135697d6, 2.690178412171373d0, &
                8.395498820952659d36, 1.283768228078111d6, 2.7349535459621346d0, &
                1.8409179480873614d37, 1.379410655833714d6, 2.772299305036374d0, &
                4.105265551141743d37, 1.4754721916179487d6, 2.811542488265903d0, &
                6.475886707043335d37, 1.6096041792672188d6, 2.822020648909448d0, &
                1.3258054862939667d38, 1.7203230250054656d6, 2.855020456297321d0, &
                2.0212930745786477d38, 1.8686117460094132d6, 2.863865259228529d0, &
                4.89932668666366d38, 1.9644384793285115d6, 2.911752011774463d0, &
                2.947497999154891d27, 1.6918556303257223d7, 0.9628968654748764d0, &
                4.4985357459233346d27, 1.784120018639596d7, 0.9820279890950028d0 /), &
                                                                     (/ 3, 37 /) )
double precision, parameter, dimension(3,55) :: Cs_param = &
    reshape( (/ 1.0067145666103756d22, 52.487440693139064d0, 2.7385019836913114d0, &
                7.850381786372777d23, 761.2255911152253d0, 2.0660197638382733d0, &
                1.2352289567990247d27, 1306.4304282001383d0, 2.841273587760102d0, &
                8.285899003417227d29, 1926.1050102158758d0, 3.5000370081744165d0, &
                1.3943473132068207d32, 2862.591109127988d0, 3.9290850597544393d0, &
                1.6090811519096703d34, 3923.686937109447d0, 4.324790783549418d0, &
                1.2197188038907911d36, 5155.954467164836d0, 4.671671863095624d0, &
                7.772313542644456d36, 7890.401181613261d0, 4.624790779520253d0, &
                3.3985065640224123d38, 9616.423940737586d0, 4.924077392169907d0, &
                4.8014562020006024d36, 21279.60493046487d0, 4.051203833568829d0, &
                1.3222551053571997d38, 24294.678658917313d0, 4.31624572540669d0, &
                1.6574025779488334d39, 28803.02083574891d0, 4.479631003881846d0, &
                1.9601994628561338d40, 33560.16146953232d0, 4.641372183771167d0, &
                2.3511108954259957d41, 38371.42266737041d0, 4.809982392624659d0, &
                1.4586939574159844d42, 45111.915987085304d0, 4.898074435696151d0, &
                1.3185543170683242d43, 51000.251990244644d0, 5.039189024375879d0, &
                1.0773662920890716d44, 57332.69476301885d0, 5.1711232098416176d0, &
                8.026118769310927d44, 64123.5110697902d0, 5.294812832444929d0, &
                6.132316537917955d45, 70939.4319707329d0, 5.424516464841734d0, &
                5.978655302384787d44, 99641.15986441393d0, 5.038556747870522d0, &
                4.5918747463503677d45, 107758.9238855643d0, 5.177095189338379d0, &
                2.7596043753895266d46, 117405.54836165802d0, 5.288917363302168d0, &
                1.8217938026575976d47, 126509.84906429503d0, 5.4131166262677d0, &
                5.5986254360520656d47, 141146.37923668668d0, 5.45214435668305d0, &
                3.043832153610278d48, 151964.45007995746d0, 5.557557219247648d0, &
                3.2582850917570695d48, 176802.88895687615d0, 5.484248396924973d0, &
                1.5867266730381928d49, 189395.76723911107d0, 5.580967644950393d0, &
                3.730087044273447d43, 433902.22052782355d0, 4.176983985293046d0, &
                1.0956515613976363d44, 467015.92197585426d0, 4.232029490341937d0, &
                2.962740161056216d44, 503077.5298274912d0, 4.279899983800397d0, &
                7.314590621543037d44, 542663.0206076808d0, 4.3198701538914825d0, &
                1.8934246830379147d45, 581882.5146795525d0, 4.365221097852297d0, &
                4.566171043429717d45, 624310.0815640229d0, 4.404596034427387d0, &
                9.012250277206341d45, 674813.9638445469d0, 4.425839936597809d0, &
                2.1956241197801815d46, 719340.5988659592d0, 4.467715823989638d0, &
                4.956319009281406d46, 767770.6575843783d0, 4.503108536433492d0, &
                1.2087646202679079d47, 814213.8503829134d0, 4.546320683994775d0, &
                4.2393229796198117d46, 949310.504573465d0, 4.412072205833181d0, &
                1.186145565511963d47, 994239.5071187183d0, 4.469534289314101d0, &
                3.0198512028692754d47, 1.0442356373731942d6, 4.518781954155857d0, &
                8.138731102071833d47, 1.091704816186824d6, 4.5735452197837745d0, &
                9.024016205173133d47, 1.1893000683619515d6, 4.548831754866617d0, &
                2.199734638466941d48, 1.24464668345491d6, 4.595460589784465d0, &
                2.4640844209070167d48, 1.3494631911312118d6, 4.57333408005987d0, &
                6.188697514903444d48, 1.4053401408214383d6, 4.6234319550448095d0, &
                2.0790590117445198d39, 4.877709249302955d6, 2.7966620649773297d0, &
                3.683305380261343d39, 5.094504263667927d6, 2.8233726484803303d0, &
                6.047476196981079d39, 5.337996286620714d6, 2.84442355319902d0, &
                1.0471081885558927d40, 5.568011076425189d6, 2.8697133664308336d0, &
                6.273346384614805d39, 6.175820914079348d6, 2.8146484527263547d0, &
                9.995132845792308d39, 6.454298659538151d6, 2.834157108994108d0, &
                1.3448971695667585d40, 6.804266629699228d6, 2.841129803690224d0, &
                2.4075380464334907d40, 7.0459254669033475d6, 2.8697175674767808d0, &
                1.6534330222238202d28, 5.8504976845801085d7, 0.9470580668607365d0, &
                2.1839164977677438d28, 6.072400966498238d7, 0.9586580582564923d0 /), &
                                                                     (/ 3, 55 /) )
double precision, parameter, dimension(3,2) :: He_param = &
    reshape( (/ 7.2207661763501d18, 832.809878216992d0, 0.48776427204592254d0, &
                2.7226733893691d21, 2742.1316798375965d0, 1.0000920088118899d0 /), &
                                                                     (/ 3, 2 /) )
double precision, parameter, dimension(3,18) :: Ar_param = &
    reshape( (/ 4.583155091154178d19, 427.3616288942586d0, 0.858307468844768d0, &
                2.347069522730368d23, 992.0665907024354d0, 1.8069357289301329d0, &
                1.931637188507353d26, 1775.9423286421918d0, 2.4675897958480384d0, &
                2.3004409750295724d28, 3141.43462246603d0, 2.8229627269683206d0, &
                3.468592152886133d30, 4422.602894394151d0, 3.263766733070865d0, &
                2.965675081118506d32, 5958.158993507227d0, 3.6326551045811497d0, &
                2.1249901067576396d33, 9479.178143644607d0, 3.629746649292696d0, &
                8.991061313075681d34, 11737.03857874398d0, 3.92742350307885d0, &
                7.611312705906995d29, 59342.85460567607d0, 2.229751529947494d0, &
                6.602179339057654d30, 71781.6325368281d0, 2.3680482459080237d0, &
                4.815600903407338d31, 85819.59318777094d0, 2.490706079578982d0, &
                2.055113377686043d32, 105186.65111959413d0, 2.558310025830369d0, &
                1.1945429165015238d33, 122591.38628416909d0, 2.663021333332675d0, &
                6.212435859768618d33, 141745.69057711778d0, 2.7584389081812857d0, &
                1.6060852787486324d34, 170916.68217921766d0, 2.783373249974532d0, &
                9.327897444360857d34, 190110.7378894357d0, 2.894937691775023d0, &
                6.6878740700780965d25, 1.8068758775178685d6, 0.9536895901987263d0, &
                1.5247916822996254d26, 2.0115330349032478d6, 0.9959340925415332d0 /), &
                                                                     (/ 3, 18 /) )
double precision, parameter, dimension(3,7) :: N_param = &
    reshape( (/ 6.441568477283914d19, 378.4956025896395d0, 0.9350660865589433d0, &
                1.4701035287257558d23, 1100.1258100528623d0, 1.711847679118431d0, &
                4.963404971976876d25, 2232.374603536868d0, 2.2130314611564854d0, &
                1.429510866064162d27, 4658.081063546513d0, 2.3525381223431663d0, &
                1.3012849307544791d29, 6615.849778853221d0, 2.7281285044063783d0, &
                2.1409091975134166d23, 88606.4475724123d0, 0.8838466900060211d0, &
                1.4212916247536375d24, 117682.27130021283d0, 0.9994495038042972d0 /), &
                                                                     (/ 3, 7 /) )
double precision, parameter, dimension(3,8) :: O_param = &
    reshape( (/ 8.471015773053202d19, 343.2810702373722d0, 0.9990920676510977d0, &
                4.659951179950614d22, 1421.770921017062d0, 1.4896379815554006d0, &
                1.375943467286961d25, 2781.3610873184325d0, 1.985966133301126d0, &
                1.4410401568574448d27, 4652.670876621719d0, 2.353837077475246d0, &
                2.1812887948673846d28, 8303.41171218499d0, 2.4562134926799617d0, &
                1.0859737168339012d30, 11088.09514588445d0, 2.766300924086769d0, &
                5.135265231597041d23, 137319.4144041063d0, 0.8991975448040144d0, &
                2.7649244731504574d24, 175715.8649323391d0, 0.999259077262769d0 /), &
                                                                     (/ 3, 8 /) )
double precision, parameter, dimension(3,6) :: C_param = &
    reshape( (/ 1.8809310466441196d20, 258.1083104918942d0, 1.1984440550739555d0, &
                5.509352426405289d23, 822.523017116882d0, 1.987881643259267d0, &
                4.572825176241475d25, 2263.6763526780164d0, 2.198152910400111d0, &
                9.835564484116034d27, 3537.9525592618193d0, 2.674447703635764d0, &
                7.536123777733801d22, 53034.28979013096d0, 0.8628040191262976d0, &
                6.587972799911127d23, 74090.46603762524d0, 0.9996157827449552d0 /), &
                                                                     (/ 3, 6 /) )
double precision, parameter, dimension(3,10) :: Ne_param = &
    reshape( (/ 1.2380498982127483d19, 684.0495717312216d0, 0.5886207199036224d0, &
                1.6861743227006676d22, 1790.87087951275d0, 1.3052851455522987d0, &
                4.012462791033353d24, 3450.2505445538854d0, 1.7789909006599105d0, &
                1.4251452589514637d26, 6544.999998027671d0, 1.9932260907792414d0, &
                6.670617838240815d27, 9689.667572133209d0, 2.282840577038655d0, &
                1.9391656057538128d29, 13557.847041347119d0, 2.522116807457813d0, &
                1.3881927947789378d30, 20383.797283570835d0, 2.586898525643846d0, &
                3.1075134088626144d31, 25254.466329809107d0, 2.8167466774677434d0, &
                2.1020658703689154d24, 282468.0585535283d0, 0.9200040384519763d0, &
                8.391062111505151d24, 343429.9456106982d0, 0.9988031600227598d0 /), &
                                                                     (/ 3, 10 /) )
double precision, parameter, dimension(3,54) :: Xe_param = &
    reshape( (/ 1.3775354155861185d20, 288.57589581574035d0, 1.1181793520578327d0, &
                1.550226815196774d24, 656.1909430408086d0, 2.2215830587335996d0, &
                2.2768338286515254d27, 1181.8699332253861d0, 2.971739666983438d0, &
                1.0222240210041819d30, 1872.6040313895708d0, 3.542491127683178d0, &
                2.1623381482112875d32, 2721.1636112083115d0, 4.013040352645607d0, &
                2.6991402989153737d34, 3721.3055273127543d0, 4.419620648906867d0, &
                2.4622862093269497d35, 5988.533565050575d0, 4.395610921298916d0, &
                1.5012400890751945d37, 7452.459544515593d0, 4.7328801021744935d0, &
                5.376912395698707d35, 16474.290444557057d0, 3.9509630414793238d0, &
                1.3066584018868181d37, 19611.207359780554d0, 4.190565038153051d0, &
                1.8497011055940987d38, 23674.858558706495d0, 4.362240665366397d0, &
                2.6681724805371856d39, 27815.543371271946d0, 4.54372243455411d0, &
                3.577640165134184d40, 32176.348524525496d0, 4.721112847931664d0, &
                2.7089693415351963d41, 38007.71401995982d0, 4.828456151144133d0, &
                2.6559087850025507d42, 43392.8782290532d0, 4.974953378510279d0, &
                2.1899377984971629d43, 49406.56731263479d0, 5.1034375243621675d0, &
                1.8601406162201934d44, 55468.870845425416d0, 5.239482352335317d0, &
                1.9313682862298425d45, 60908.78844944052d0, 5.403664925396409d0, &
                1.5426588652718157d44, 87869.02232374942d0, 4.982156226229307d0, &
                1.1590871169485781d45, 95909.52818579777d0, 5.115878165169938d0, &
                8.030370382706514d45, 104435.5470090333d0, 5.241935268399854d0, &
                5.382165292841963d46, 113200.23085498222d0, 5.365848876529456d0, &
                1.8217938026575976d47, 126509.84906429503d0, 5.4131166262677d0, &
                1.077849728298342d48, 136393.60209955153d0, 5.526234224434658d0, &
                1.145516874883936d48, 159811.20754403196d0, 5.44842538130848d0, &
                6.014132302364118d48, 171384.410225511d0, 5.5518758215736685d0, &
                1.759085923766389d43, 394064.2066635542d0, 4.1549467350722304d0, &
                5.305157185638429d43, 425345.21856075496d0, 4.211470223441139d0, &
                1.4911016576131055d44, 459078.08137395576d0, 4.2620127128282785d0, &
                3.746338544125226d44, 496648.9293453104d0, 4.3025833188760005d0, &
                1.0200735803192588d45, 533001.5075147987d0, 4.351821613265934d0, &
                2.4271623447664046d45, 574235.1060028574d0, 4.3889333565493684d0, &
                4.9683984737540496d45, 621542.8884528814d0, 4.412604825238724d0, &
                1.2484414815313788d46, 663477.9943518823d0, 4.456566977558425d0, &
                2.894901889446211d46, 709203.4692145847d0, 4.493643943279613d0, &
                7.180428844235488d46, 753475.3199613485d0, 4.53769337085137d0, &
                2.448645505910615d46, 882711.0267638017d0, 4.398978788274832d0, &
                7.090658854982765d46, 925001.57509592d0, 4.459072393688684d0, &
                1.860229435466672d47, 972223.3666808198d0, 4.510512643085584d0, &
                5.064124828052308d47, 1.0180486819614146d6, 4.5657012421984975d0, &
                5.9565042899701805d47, 1.1084286245427015d6, 4.545372106494971d0, &
                1.4905062982964533d48, 1.1608077808083333d6, 4.593864401786412d0, &
                1.660914855730416d48, 1.2615303421739575d6, 4.570386159811437d0, &
                4.2520741495048395d48, 1.3148814271135165d6, 4.621771988602169d0, &
                1.5790651376647514d39, 4.579524976873413d6, 2.7930489801003286d0, &
                2.8386327726136714d39, 4.786413465857762d6, 2.8206493017311245d0, &
                4.700763570275373d39, 5.0203491847052295d6, 2.8421056006525185d0, &
                8.200335703674566d39, 5.242021070951296d6, 2.8677440883669703d0, &
                5.20578174350226d39, 5.804149130601667d6, 2.8165062116744424d0, &
                8.391595879280377d39, 6.070112038854565d6, 2.8366647300643155d0, &
                1.135074333575459d40, 6.406097404576862d6, 2.8437495112068905d0, &
                2.0553597017159955d40, 6.637543043722892d6, 2.8730250674268727d0, &
                1.542287349895057d28, 5.520489369472207d7, 0.9483460319455781d0, &
                2.0492457824596215d28, 5.733209089671276d7, 0.9602460194755702d0 /), &
                                                                     (/ 3, 54 /) )

double precision, parameter, dimension(3,18) :: Yb_param = &
    reshape( (/ 1.998729582135397d21, 1.068604450949123d02, 1.949272871592323d0, & 
                6.266630953160455d25, 2.903961394160274d02, 3.226884022732426d0, &
                1.633490350054849d28, 8.567480112611549d02, 3.420697880666071d0, & 
                7.033126033817732d29, 1.967626338656452d03, 3.467518893419359d0, &
                1.813176768210167d31, 3.630105579722960d03, 3.553208464047432d0, &
                9.495491157219170d31, 6.730018983991275d03, 3.447675593140554d0, &
                5.358792557003665d33, 8.646538707567524d03, 3.773144244929749d0, &
                2.432849794882799d35, 1.071677592250407d04, 4.078349104741582d0, &
                2.980985396748150d36, 1.421820933251101d04, 4.199339269018307d0, &
                9.090128804380945d37, 1.677531096121317d04, 4.467184823115996d0, &
                1.110784555706743d39, 2.064348074619395d04, 4.612017181201895d0, &
                6.552418197839887d39, 2.604043512816300d04, 4.666120192229148d0, &
                4.124326514324229d40, 3.183976123298108d04, 4.740380764911407d0, &
                6.913068158587253d41, 3.567904659544683d04, 4.951743930475629d0, &
                9.931519433973369d42, 3.984565888597911d04, 5.146362971528592d0, &
                5.430823951587641d43, 4.666779146289682d04, 5.219682021168320d0, &
                3.021561896285397d44, 5.384015187193020d04, 5.300873756949119d0, &
                1.794945766837406d45, 6.113351998916520d04, 5.394893554476560d0 /), &
                                                                        (/ 3, 18 /) )

public :: neutral

type neutral

  private

  class(fdist2d), pointer :: pf => null()
  class(part2d), pointer :: part => null()
  class(part2d_buf), pointer :: part_add => null()

  class(field_rho), allocatable :: q
  class(field_jay), allocatable :: cu, amu
  class(field_djdxi), allocatable :: dcu

  ! max ionization level
  integer :: multi_max
  ! index of neutral gas residue in multi_ion array
  integer :: idx_neut
  ! index of total ion in multi_ion array
  integer :: idx_ion
  ! wp is plasma frequency in SI unit
  real :: wp, dt
  ! pusher type
  integer :: push_type
  
  double precision, dimension(:,:), pointer :: rate_param

  ! in-cell ionization level. There are multi_max + 2 levels, the multi_max+1 level 
  ! is for the neutral, the first to multi_max levels are for ions, and the last
  ! level is for the total ion.
  real, dimension(:,:,:), allocatable :: multi_ion
  real, dimension(:,:), allocatable :: ion_old

  ! the ion charge density for diagnostic 
  class(field_rho), allocatable :: rho_ion, rho_ion_add

  contains

  procedure :: alloc  => alloc_neutral
  procedure :: new    => init_neutral
  procedure :: renew  => renew_neutral
  procedure :: update => update_neutral
  procedure :: del    => end_neutral
  procedure :: qdp    => qdeposit_neutral
  procedure :: amjdp  => amjdeposit_neutral
  procedure :: push   => push_neutral
  procedure :: psend  => psend_neutral
  procedure :: precv  => precv_neutral
  procedure :: wr     => writehdf5_neutral
  procedure :: wrq    => writeq_neutral
  procedure :: wr_ion => write_ion_neutral
  procedure :: cbq    => cbq_neutral
  procedure :: get_multi_max
  procedure :: ion_deposit => ion_deposit_neutral

end type neutral

! This extended type only store the position for neutral deposition
type, extends( part2d ) :: part2d_buf

   contains

   procedure :: new      => init_part2d_buf
   procedure :: end      => end_part2d_buf
   procedure :: pipesend => pipesend_part2d_buf
   procedure :: piperecv => piperecv_part2d_buf

end type part2d_buf

save

character(len=10) :: cls_name = 'neutral'
integer, parameter :: cls_level = 2

contains

subroutine alloc_neutral( this )

   implicit none

   class(neutral), intent(inout) :: this

   if ( .not. associated( this%part ) ) then
      allocate( part2d :: this%part )
      allocate( part2d_buf :: this%part_add )
   endif

end subroutine alloc_neutral

subroutine init_neutral( this, opts, pf, max_mode, elem, max_e, qbm, wp, s, &
  push_type, smth_type, smth_ord )
! element is atomic number. e.g.: For Li, element = 3
! max_e is the maximum number of electrons that the programmer allow the atom
! to lose due to the ionization. It should be less or equal to element
   
  implicit none

  class(neutral), intent(inout) :: this
  type(options), intent(in) :: opts
  class(fdist2d), intent(inout), target :: pf
  integer, intent(in) :: max_mode, elem, max_e, push_type
  real, intent(in) :: qbm, wp, s
  integer, intent(in), optional :: smth_type, smth_ord

  ! local data
  character(len=18), save :: sname = 'init_neutral'
  integer :: i, nrp, n_theta
  real :: pi2_ntheta

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%pf => pf
  this%wp = wp
  this%dt = opts%get_dxi()
  this%push_type = push_type

  allocate( this%q, this%cu, this%amu, this%dcu, this%rho_ion, this%rho_ion_add )

  if ( present(smth_type) .and. present(smth_ord) ) then
    call this%q%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%rho_ion%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%rho_ion_add%new( opts, max_mode, p_ps_linear, smth_type, smth_ord, has_2d=.false. )
    call this%cu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%dcu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
    call this%amu%new( opts, max_mode, p_ps_linear, smth_type, smth_ord )
  else
    call this%q%new( opts, max_mode, p_ps_linear )
    call this%rho_ion%new( opts, max_mode, p_ps_linear )
    call this%rho_ion_add%new( opts, max_mode, p_ps_linear, has_2d=.false. )
    call this%cu%new( opts, max_mode, p_ps_linear )
    call this%dcu%new( opts, max_mode, p_ps_linear )
    call this%amu%new( opts, max_mode, p_ps_linear )
  endif

  call this%part%new( opts, pf, qbm, this%dt, s, if_empty=.true. )
  call this%part_add%new( opts, pf, qbm, this%dt, s )

  this%q  = 0.0
  this%cu = 0.0
  this%rho_ion = 0.0
  this%rho_ion_add = 0.0

  ! call this%q%ag()

  ! TODO: WHY NEED THIS?
  if ( id_stage() == 0 ) then
    call this%q%copy_slice( 1, p_copy_1to2 )      
  endif

  this%multi_max = max_e
  this%idx_neut  = this%multi_max + 1
  this%idx_ion   = this%multi_max + 2

  select case ( elem )
    
  case( p_neut_H )

    this%multi_max = min( this%multi_max, size( H_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = H_param( :, 1:this%multi_max )

  case( p_neut_He )

    this%multi_max = min( this%multi_max, size( He_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = He_param( :, 1:this%multi_max )

  case( p_neut_Li )

    this%multi_max = min( this%multi_max, size( Li_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Li_param( :, 1:this%multi_max )

  case( p_neut_C )

    this%multi_max = min( this%multi_max, size( C_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = C_param( :, 1:this%multi_max )

  case( p_neut_N )

    this%multi_max = min( this%multi_max, size( N_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = N_param( :, 1:this%multi_max )

  case( p_neut_O )

    this%multi_max = min( this%multi_max, size( O_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = O_param( :, 1:this%multi_max )

  case( p_neut_Ne )

    this%multi_max = min( this%multi_max, size( Ne_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Ne_param( :, 1:this%multi_max )

  case( p_neut_Na )

    this%multi_max = min( this%multi_max, size( Na_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Na_param( :, 1:this%multi_max )

  case( p_neut_Ar )

    this%multi_max = min( this%multi_max, size( Ar_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Ar_param( :, 1:this%multi_max )

  case( p_neut_K )

    this%multi_max = min( this%multi_max, size( K_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = K_param( :, 1:this%multi_max )

  case( p_neut_Rb )

    this%multi_max = min( this%multi_max, size( Rb_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Rb_param( :, 1:this%multi_max )

  case( p_neut_Xe )

    this%multi_max = min( this%multi_max, size( Xe_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Xe_param( :, 1:this%multi_max )

  case( p_neut_Cs )

    this%multi_max = min( this%multi_max, size( Cs_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Cs_param( :, 1:this%multi_max )

  case( p_neut_Yb )

    this%multi_max = min( this%multi_max, size( Yb_param, 2 ) )   
    allocate( this%rate_param( 3, this%multi_max ) )
    this%rate_param = Yb_param( :, 1:this%multi_max )

  case default
     
    call write_err( 'Invalid neutral gas species!' )

  end select

  ! initialize the ionization level array
  nrp = opts%get_ndp(1)
  n_theta = pf%num_theta
  allocate( this%multi_ion( nrp, n_theta, this%multi_max + 2 ) )
  allocate( this%ion_old( nrp, n_theta ) )
  this%multi_ion = 0.0
  this%multi_ion( :, :, this%idx_neut ) = 1.0

  ! normalized longitudinal density
  ! this%den_lon = this%pf%get_den_lon(0.0)
  
  ! normalized charge (coordinate of the particle array) 
  ! this%qm = this%multi_max * this%density * this%pf%den * this%pf%qm / &
  !   ( abs(this%pf%qm) * real(this%pf%ppc1) * real(this%pf%ppc2) * real(n_theta) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_neutral

subroutine update_neutral( this, e, psi, s )

  implicit none
  class(neutral), intent(inout) :: this
  class(field_e), intent(in) :: e
  class(field_psi), intent(in) :: psi
  real, intent(in) :: s
  ! local data
  character(len=18), save :: sname = 'update_neutral'
  integer :: i

  call write_dbg(cls_name, sname, cls_level, 'starts')

  this%ion_old = this%multi_ion( :, :, this%idx_ion )

  call ionize_neutral( this%multi_ion, this%rate_param, e, this%wp, this%dt, this%pf%ppc )

  ! here multi_ion(:,:,idx_ion) and ion_old should be discrete.
  call add_particles( this%pf, this%part, this%part_add, this%multi_ion, this%ion_old, s )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine update_neutral

subroutine ionize_neutral( multi_ion, adk_coef, e, wp, dt, ppc )

  implicit none

  real, intent(inout), dimension(:,:,:) :: multi_ion
  double precision, dimension(:,:), intent(in) :: adk_coef
  class(field_e), intent(in) :: e
  real, intent(in) :: wp, dt
  integer, intent(in), dimension(2) :: ppc

  ! local data
  character(len=18), save :: sname = 'ionize_neutral'
  real, dimension(:,:), pointer :: e_re => null(), e_im => null()
  integer :: i, j, k, m, nrp, n_theta, noff, ppc_tot, max_mode, multi_max, idx_neut, idx_ion
  real :: eij, cons, inj, dens_temp, theta, pi2_ntheta
  real :: w_ion(20)
  real, dimension(:), allocatable :: e1, e2, e3
  complex(kind=DB) :: phase, phase_inc
  logical :: shoot

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp        = e%rf_re(0)%get_ndp(1)
  n_theta    = size( multi_ion, 2 )
  noff       = e%rf_re(0)%get_noff(1)
  max_mode   = e%get_max_mode()
  ppc_tot    = ppc(1) * ppc(2)
  pi2_ntheta = 2.0 * pi / real( n_theta )
  multi_max  = size( multi_ion, 3 ) - 2
  idx_neut   = multi_max + 1
  idx_ion    = multi_max + 2

  allocate( e1(nrp), e2(nrp), e3(nrp) )

  do k = 1, n_theta

    theta = pi2_ntheta * ( k - 1 )
    phase_inc = cmplx( cos(theta), sin(theta) )

    ! calculate the in-cell values of electric fields
    ! m = 0 mode
    e_re => e%rf_re(0)%get_f1()
    phase = cmplx( 1.0, 0.0 )
    do j = 1, nrp
      e1(j) = 0.5 * ( e_re(1,j) + e_re(1,j+1) )
      e2(j) = 0.5 * ( e_re(2,j) + e_re(2,j+1) )
      e3(j) = 0.5 * ( e_re(3,j) + e_re(3,j+1) )
    enddo

    ! m > 0 mode
    do m = 1, max_mode

      e_re => e%rf_re(m)%get_f1()
      e_im => e%rf_re(m)%get_f1()
      phase = phase * phase_inc

      do j = 1, nrp
        e1(j) = e1(j) + ( e_re(1,j) + e_re(1,j+1) ) * real(phase) - ( e_im(1,j) + e_im(1,j+1) ) * aimag(phase)
        e2(j) = e2(j) + ( e_re(2,j) + e_re(2,j+1) ) * real(phase) - ( e_im(2,j) + e_im(2,j+1) ) * aimag(phase)
        e3(j) = e3(j) + ( e_re(3,j) + e_re(3,j+1) ) * real(phase) - ( e_im(3,j) + e_im(3,j+1) ) * aimag(phase)
      enddo

    enddo

    do j = 1, nrp

      ! total electric field in unit of GV/m
      eij = sqrt( e1(j)**2 + e2(j)**2 + e3(j)**2 ) * wp * 1.708d-12

      if ( eij > 1.0d-6 ) then

        if ( multi_ion( j, k, idx_ion ) < real(multi_max) ) then

          shoot = .false.

          ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
          ! w_ion is in normalized unit now
          do i = 1, multi_max
            w_ion(i) = adk_coef(1,i) * eij ** (-adk_coef(3,i)) * exp( -adk_coef(2,i) / eij ) / wp
          enddo

          ! Use the wrong 2nd order Runge Kutta method in OSIRIS
          ! TODO: check with Yujian
          cons = multi_ion(j, k, idx_neut) * w_ion(1) * dt * ( 1.0 + 0.5 * w_ion(1) * dt )

          ! overshoot
          if( cons > multi_ion(j, k, idx_neut) ) then
             shoot = .true.
             cons = multi_ion(j, k, idx_neut)
          endif
          ! subtract from 0.
          multi_ion(j, k, idx_neut) = multi_ion(j, k, idx_neut) - cons

          ! For the levels in the middle
          do i = 1, multi_max - 1

            ! remember how much charge will be added to i-th level
            inj = cons

            ! overshoot in previous level (we go to time centered densities)
            ! TODO: is dens_temp correct?
            if (shoot) then
              dens_temp = inj * 0.5
              shoot = .false.
            else
              dens_temp = 0.0
            endif

            ! ionizing m-1. level (index m) adding to m. level
            ! 2nd order Runge-Kutta (wrong)
            cons = ( multi_ion(j,k,i) + dens_temp ) * w_ion(i+1) * dt * ( 1.0 + 0.5 * w_ion(i+1) * dt )

            ! overshoot in current level (temp might come from previous overshoot)
            if ( cons > multi_ion(j,k,i) + dens_temp ) then
              shoot = .true.
              cons = multi_ion(j,k,i) + dens_temp
            endif

            ! update density for i-th level
            ! Physically, multi_ion(j,k,i) - cons + inj will never exceed 1.0
            ! min is only used to round off the tiny error introduced by a double number (~10^-16)
            multi_ion(j,k,i) = min( multi_ion(j,k,i) - cons + inj, 1.0 )

          enddo ! all levels in the middle

          ! update the last level
          multi_ion(j, k, multi_max) = min( multi_ion(j, k, multi_max) + cons, 1.0 )

          ! End use Runge Kutta method

          ! calculate the total ion charge density by summing over the charge density from all 
          ! types of ions 
          multi_ion(j, k, idx_ion) = 0.0
          do i = 1, multi_max
            multi_ion(j, k, idx_ion) = multi_ion(j, k, idx_ion) + real(i) * multi_ion(j,k,i)
          enddo

          ! We want the total ion level that corresponds to release an electron at 'half' step (like OSIRIS and qpic2.0)
          ! Now the multi_ion(:,:,idx_ion) is discrete
          multi_ion(j, k, idx_ion) = real(multi_max) / real(ppc_tot) * &
            int( multi_ion(j, k, idx_ion) * real(ppc_tot) / real(multi_max) + 0.5 )

        endif 

      endif ! eij > 1e-6

    enddo
  enddo

  deallocate( e1, e2, e3 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine ionize_neutral

subroutine add_particles( prof, part, part_add, multi_ion, ion_old, s )

  implicit none

  class(fdist2d), intent(inout) :: prof
  class(part2d), intent(inout) :: part
  class(part2d_buf), intent(inout) :: part_add
  real, intent(in), dimension(:,:,:) :: multi_ion
  real, intent(in), dimension(:,:) :: ion_old
  real, intent(in) :: s

  ! local
  character(len=18), save :: sname = 'add_particles'
  integer :: nrp, noff, ppc_tot, i, j, k, pp1, pp2, ppc_add, multi_max, idx_ion, n_theta
  real :: rn, dr, theta, dtheta, den_lon, den_perp, coef
  real, dimension(2) :: x_tmp

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  nrp          = prof%nrp
  n_theta      = prof%num_theta
  noff         = prof%noff
  dr           = prof%dr
  dtheta       = 2.0 * pi / real( n_theta )
  ppc_tot      = product( prof%ppc )
  multi_max    = size( multi_ion, 3 ) - 2
  idx_ion      = multi_max + 2
  part_add%npp = 0

  call prof%get_den_lon( s, prof%prof_pars_lon, den_lon )
  coef = real(multi_max) * sign(1.0, prof%qm) / ( real(ppc_tot) * real(n_theta) )

  do k = 1, n_theta
    do j = 1, nrp

      ! calculate the # of particles to inject based on the difference in ion density between 
      ! the previous time step and the current time step.
      ! here multi_ion(:,:,idx_ion) and ion_old should both be discrete.
      ppc_add = int( ( multi_ion(j, k, idx_ion) - ion_old(j, k) ) / multi_max * real(ppc_tot) + 0.5 )

      pp1 = part%npp
      pp2 = part_add%npp
      do i = 1, ppc_add

        ! randomly generate injection position
        ! r     = rand() + j + noff - 1
        ! theta = ( rand() + k - 1.5 ) * dtheta
        rn    = j + noff - 1 + ( real(i) - 0.5 ) / ppc_add
        theta = ( k - 1 ) * dtheta
        x_tmp(1) = rn * dr * cos(theta)
        x_tmp(2) = rn * dr * sin(theta)
        
        call prof%get_den_perp( x_tmp, s, prof%prof_pars_perp, prof%prof_pars_lon, den_perp )

        if ( den_lon * den_perp * prof%density < prof%den_min ) cycle

        pp1 = pp1 + 1
        pp2 = pp2 + 1

        part%x(1,pp1)   = x_tmp(1)
        part%x(2,pp1)   = x_tmp(2)
        part%q(pp1)     = rn * den_lon * den_perp * prof%density * coef
        part%p(1,pp1)   = 0.0
        part%p(2,pp1)   = 0.0
        part%p(3,pp1)   = 0.0
        part%gamma(pp1) = 1.0
        part%psi(pp1)   = 1.0

        ! store the position of ionized particles for ion charge deposition
        part_add%x(1,pp2) = part%x(1,pp1)
        part_add%x(2,pp2) = part%x(2,pp1)
        part_add%q(pp2)   = -part%q(pp1) ! note the sign

      enddo
      part%npp = part%npp + ppc_add
      part_add%npp = part_add%npp + ppc_add

    enddo
  enddo

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine add_particles

subroutine renew_neutral( this, s )

  implicit none

  class(neutral), intent(inout) :: this
  real, intent(in) :: s

  ! local data
  character(len=18), save :: sname = 'renew_neutral'
  integer :: i
        
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%part%renew( this%pf, s, if_empty=.true. )

  this%q = 0.0
  this%cu = 0.0
  this%rho_ion = 0.0

  this%multi_ion = 0.0
  this%multi_ion( :, :, this%idx_neut ) = 1.0

  ! TODO: WHY NEED THIS
  ! if ( id_stage() == 0 ) then
  !   call this%q%copy_slice( 1, p_copy_1to2 )       
  ! endif

  ! normalized longitudinal density
  ! this%den_lon = this%pf%get_den_lon(s)
  
  ! normalized charge (coordinate of the particle array) 
  ! this%qm = this%multi_max * this%density * this%pf%den * this%pf%qm &
  !   / ( abs(this%pf%qm) * real(this%pf%ppc1) * real(this%pf%ppc2) * real(this%pf%num_theta) )

  this%part%npp = 0
  this%part_add%npp = 0

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine renew_neutral

subroutine qdeposit_neutral( this, q_tot )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_rho), intent(inout) :: q_tot
! local data
  character(len=18), save :: sname = 'qdeposit_neutral'
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%q = 0.0

  call this%part%qdeposit( this%q )
  call this%q%acopy_gc_f1( dir=p_mpi_forward )
  call this%q%smooth_f1()
  call this%q%copy_gc_f1()
  
  call add_f1( this%q, q_tot )
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine qdeposit_neutral

subroutine ion_deposit_neutral( this, q_tot )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_rho), intent(inout) :: q_tot
  ! local data
  character(len=18), save :: sname = 'ion_deposit_neutral'
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%rho_ion_add = 0.0

  call this%part_add%qdeposit( this%rho_ion_add )
  call this%rho_ion_add%acopy_gc_f1( dir=p_mpi_forward )
  call this%rho_ion_add%smooth_f1()
  call this%rho_ion_add%copy_gc_f1()
  
  call add_f1( this%rho_ion_add, this%rho_ion )
  call add_f1( this%rho_ion, q_tot )

  ! clear up the buffer
  this%part_add%npp = 0
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine ion_deposit_neutral

subroutine amjdeposit_neutral( this, e, b, cu, amu, dcu )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_jay), intent(inout) :: cu, amu
  class(field_djdxi), intent(inout) :: dcu
  class(field_e), intent(in) :: e
  class(field_b), intent(in) :: b
  ! local data
  character(len=18), save :: sname = 'amjdeposit_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%cu = 0.0
  this%dcu = 0.0
  this%amu = 0.0
  select case ( this%push_type )
    case ( p_push2_robust )
      call this%part%amjdeposit_robust( e, b, this%cu, this%amu, this%dcu )
    case ( p_push2_clamp )
      call this%part%amjdeposit_clamp( e, b, this%cu, this%amu, this%dcu )
    case ( p_push2_robust_subcyc )
      call this%part%amjdeposit_robust_subcyc( e, b, this%cu, this%amu, this%dcu )
  end select
  
  call this%cu%acopy_gc_f1( dir=p_mpi_forward )
  call this%dcu%acopy_gc_f1( dir=p_mpi_forward )
  call this%amu%acopy_gc_f1( dir=p_mpi_forward )
  call this%cu%smooth_f1()
  call this%dcu%smooth_f1()
  call this%amu%smooth_f1()
  call this%cu%copy_gc_f1()
  call this%dcu%copy_gc_f1()
  call this%amu%copy_gc_f1()
  call add_f1( this%cu, cu )
  call add_f1( this%dcu, dcu )
  call add_f1( this%amu, amu )

  call write_dbg( cls_name, sname, cls_level, 'ends' )
  
end subroutine amjdeposit_neutral

subroutine push_neutral( this, e, b )
  
  implicit none
  
  class(neutral), intent(inout) :: this
  class(field_e), intent(in) :: e
  class(field_b), intent(in) :: b
! local data
  character(len=18), save :: sname = 'push_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  select case ( this%push_type )
    case ( p_push2_robust )
      call this%part%push_robust( e, b )
    case ( p_push2_clamp )
      call this%part%push_clamp( e, b )
    case ( p_push2_robust_subcyc )
      call this%part%push_robust_subcyc( e, b )
  end select

  call this%part%update_bound()
  call move_part2d_comm( this%part )
  
  call write_dbg( cls_name, sname, cls_level, 'ends' )
     
end subroutine push_neutral

subroutine end_neutral( this )

  implicit none
  class(neutral), intent(inout) :: this

end subroutine end_neutral

subroutine psend_neutral( this, tag, id )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in), dimension(4) :: tag
  integer, intent(inout), dimension(4) :: id
  ! local data
  character(len=18), save :: sname = 'psend_neutral'
  integer :: i, idproc_des, ierr, count
           
  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  call this%part%pipesend( tag(1), id(1) )
  call this%part_add%pipesend( tag(2), id(2) )
  call this%rho_ion%pipe_send( tag(3), id(3), 'forward' )

  if ( id_stage() == num_stages() - 1 ) then
    id(4) = MPI_REQUEST_NULL
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  idproc_des = id_proc() + num_procs_loc()
  count = size( this%multi_ion )

  call mpi_isend( this%multi_ion, count, p_dtype_real, idproc_des, tag(4), &
    comm_world(), id(4), ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_ISEND failed.' )
  endif

  call stop_tprof( 'pipeline' )
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine psend_neutral

subroutine precv_neutral( this, tag )

  implicit none

  class(neutral), intent(inout) :: this
  integer, intent(in), dimension(4) :: tag
  ! local data
  character(len=18), save :: sname = 'precv_neutral'
  integer, dimension(MPI_STATUS_SIZE) :: stat
  integer :: i, idproc_src, ierr, count

  call write_dbg( cls_name, sname, cls_level, 'starts' )
  call start_tprof( 'pipeline' )

  call this%part%piperecv( tag(1) )
  call this%part_add%piperecv( tag(2) )
  call this%rho_ion%pipe_recv( tag(3), 'forward', 'replace' )

  if ( id_stage() == 0 ) then
    call stop_tprof( 'pipeline' )
    call write_dbg( cls_name, sname, cls_level, 'ends' )
    return
  endif

  idproc_src = id_proc() - num_procs_loc()
  count = size( this%multi_ion)

  call mpi_recv( this%multi_ion, count, p_dtype_real, idproc_src, tag(4), &
    comm_world(), stat, ierr )
  ! check for error
  if ( ierr /= 0 ) then
    call write_err( 'MPI_RECV failed.' )
  endif
           
  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine precv_neutral

subroutine writehdf5_neutral( this, file )

  implicit none

  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in) :: file
  ! local data
  character(len=18), save :: sname = 'writehdf5_neutral'

  call write_dbg(cls_name, sname, cls_level, 'starts')

  call this%part%wr(file)

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine writehdf5_neutral

subroutine writeq_neutral( this, files, rtag, stag, id )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in), dimension(:) :: files
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id
  ! local data
  character(len=18), save :: sname = 'writeq_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%q%write_hdf5( files, 1, rtag, stag, id ) 

  call write_dbg( cls_name, sname, cls_level, 'ends' )  

end subroutine writeq_neutral

subroutine write_ion_neutral( this, files, rtag, stag, id )

  implicit none
  
  class(neutral), intent(inout) :: this
  class(hdf5file), intent(in), dimension(:) :: files
  integer, intent(in) :: rtag, stag
  integer, intent(inout) :: id
  ! local data
  character(len=18), save :: sname = 'write_ion_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call this%rho_ion%write_hdf5( files, 1, rtag, stag, id ) 

  call write_dbg( cls_name, sname, cls_level, 'ends' )  

end subroutine write_ion_neutral
   
subroutine cbq_neutral( this, pos )

  implicit none
  
  class(neutral), intent(inout) :: this
  integer, intent(in) :: pos
  ! local data
  character(len=18), save :: sname = 'cbq_neutral'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  call add_f1( this%cu, this%q, (/3/), (/1/) )
  call this%q%copy_slice( pos, p_copy_1to2 )

  call this%rho_ion%copy_slice( pos, p_copy_1to2 )

  call write_dbg( cls_name, sname, cls_level, 'ends' )

end subroutine cbq_neutral

function get_multi_max(this)

  implicit none

  class(neutral), intent(in) :: this
  integer :: get_multi_max

  get_multi_max = this%multi_max

end function get_multi_max

subroutine init_part2d_buf( this, opts, pf, qbm, dt, s, if_empty )

  implicit none

  class(part2d_buf), intent(inout) :: this
  type(options), intent(in) :: opts
  class(fdist2d), intent(inout) :: pf
  real, intent(in) :: qbm, dt, s
  logical, intent(in), optional :: if_empty

  ! local data
  character(len=18), save :: sname = 'init_part2d_buf'

  call write_dbg( cls_name, sname, cls_level, 'starts' )

  this%qbm = qbm
  this%dt  = dt
  this%part_dim = 3

  ! this is max number of particles to be ionized
  this%npmax = product( pf%ppc ) * pf%num_theta * opts%get_ndp(1)
  this%npp   = 0

  this%dr   = opts%get_dr()
  this%edge = opts%get_nd(1) * this%dr

  ! only allocate position and charge, no need to initialize with profile
  allocate( this%x( 2, this%npmax ), this%q( this%npmax ) )

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine init_part2d_buf

subroutine end_part2d_buf(this)

   implicit none

   class(part2d_buf), intent(inout) :: this
   character(len=18), save :: sname = 'end_part2d_buf'

   call write_dbg(cls_name, sname, cls_level, 'starts')

   deallocate( this%x, this%q )

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine end_part2d_buf

subroutine pipesend_part2d_buf(this, tag, id)

  implicit none

  class(part2d_buf), intent(inout) :: this
  integer, intent(in) :: tag
  integer, intent(inout) :: id

  ! local data
  character(len=18), save :: sname = 'pipesend_part2d_buf'
  integer :: des, ierr, i
  real, dimension(:,:), allocatable, save :: sbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(sbuf) ) then
    allocate( sbuf( this%part_dim, this%npmax ) )
  endif

  des = id_proc() + num_procs_loc()

  if ( des >= num_procs() ) then
    id = MPI_REQUEST_NULL
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! to be implemented if using tile
  ! call this%pcb()

  do i = 1, this%npp
    sbuf(1:2,i) = this%x(:,i)
    sbuf(3,i)   = this%q(i)
  enddo

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_ISEND(sbuf, int(this%npp*this%part_dim), p_dtype_real, des, tag, &
    comm_world(), id, ierr)

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI_ISEND failed')
  endif

  call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine pipesend_part2d_buf

subroutine piperecv_part2d_buf(this, tag)

  implicit none

  class(part2d_buf), intent(inout) :: this
  integer, intent(in) :: tag
  ! local data
  character(len=18), save :: sname = 'piperecv_part2d_buf'
  integer, dimension(MPI_STATUS_SIZE) :: istat
  integer :: nps, des, ierr, i
  real, dimension(:,:), allocatable, save :: rbuf

  call write_dbg(cls_name, sname, cls_level, 'starts')

  if ( .not. allocated(rbuf) ) then
    allocate( rbuf( this%part_dim, this%npmax ) )
  endif

  des = id_proc() - num_procs_loc()

  if (des < 0) then
    call write_dbg(cls_name, sname, cls_level, 'ends')
    return
  endif

  ! NOTE: npp*xdim might be larger than MAX_INT32
  call MPI_RECV(rbuf, int(this%npmax*this%part_dim), p_dtype_real, &
    des, tag, comm_world(), istat, ierr)

  call MPI_GET_COUNT(istat, p_dtype_real, nps, ierr)

  this%npp = nps/this%part_dim

  do i = 1, this%npp
    this%x(:,i)   = rbuf(1:2,i)
    this%q(i)     = rbuf(3,i)
  enddo

  ! to be implemented if using tile
  ! call this%pcp(fd)

  ! check for errors
  if (ierr /= 0) then
    call write_err('MPI failed')
  endif

   call write_dbg(cls_name, sname, cls_level, 'ends')

end subroutine piperecv_part2d_buf

end module neutral_class