from Crypto_tools import *
from itertools import permutations


P_bits = 444
Q_bits = 666
R_bits = 333
e = 0x1337
N = [8321077117329356263954581766837194016859681833859374146551469738742553789565498761528408178000096341991081753628879035591190841107228873036782248755852096597317053559269854941020999105514186022112075838112491499884564745335454966665835001848999256218403570429047541524272647861099813598603292295695775504244505874838019009730562620216, 9663141503982563384103774905603769762205667685102275298721284964403449121449261483138514307090449027807047697811539118959328065885920230514670112839967221129701708335087378871176539521374006686377418843364889059913595942583737991465545688834167085579154677350865488342245644093471665857007588133415608450554035129609049971856915687905, 7080633525505006454857949889380886258474613936169915325357991912983821798902257837234148311383635716165386646093418183743215120431715933036921480432793786600194625124412063608565640368381660643929081066605712749838630092722581230189543696229548387666046034403406721477818265752443487173947232032487026509033018565048660685068628813900, 6348260112940191945095264085450804431350547836448100928568733296493334845533262312663784701113046746633754397388581143523779371919889446537883618910341310362947454041409541741124231406605654693961025815677260091930522280737378333554694234586235628097018111636491016261274584391201625445101930372416766825601778759258114504767453644116, 5453076441876067965962987075376616480678826248967242473452690159966124023105342358461519279607336831463834252487811766366887983597918775466645199024629945952612311092154451713992362747056489054444283302530500475646346968235522866557033556131387030809018889094871163097069588840212411234716020934787724757538389436231926479459435139729, 8539092301764573132139384894241535432591998166686651428176862041680365196821019488767353225937458669267710968146783359781905669306801237140282934737328995064410005908343087032652676207615356474601039341695241288937409773110450995503958663731870314493254162942656333393836208952363635746218345484752702086447693272390570444229475023063, 6519174659211289290465989985638494640591837577268694359892571134942820094011179335155100770258746122812579164176239810640710311061280175847537113068566174335469669480601670136633042879486860594176408555229541384726863069317134887988109513500264430505092783429319336134768414472691544453196941818593187717904852926051114502346647888426, 6961711680362025924587083752271982856615461409839316941574792747717174299272141804413488210456939607971950020060573061131683587521057101549266612028332060247289306929784809695126285828514889617483817766102136513569391338432827451306976412448973958662677599051104723480081915666230391091496248279567627159875132838397868850500912784991, 5464727156582411007360377345208743900616053705663005668786499961992377236151787056734059201141143278419642045996194426551642956434137382420547801325331080754615348522599387497424727626522195285281557448617564720973185851757117201639221059630112719699600469361046245232476531622290326229811709079901764559526124144866364301870192468062]
C = [800378461059400239726680783421062702546581299113618553895453491207714321944554499622887232532612118204284779120928524046451494597619154079853122057618867592408424421335915888671560524092660578952242621890439766919785431411789789232309134048322721650012432166587969915464252995890054635969469155870141839815222805619769926841873928532, 8685468246369062574820183134847029157229023858170863526469628501966638181721681547114662091162797572149161013458000532984909663639626346493828947027439012131912176125653717020650233650230608573276523862941298063827867869000918623143520066067119099918633173584693454642685133071154989133688921507896223776765538029556643440655490373815, 4635611296372235589362291842711945807825964919968727011796279830725567747087132786100965922682161492876463568645940638975728831156672106717718242995621775763972524561170035488180440169190421072680121175269363490806991700969253196228364471655905426168859215651384013432550900570720981720919343406680667172355882395426739478340888619526, 742099161415136628218807400531862454374875770332166710320711769923774839345990297388615457974047093967700994503615622499139058644697776431451063778211507061619987487339659448529973693084099004007650792690143707738139278995489030445825999666762580518420739517354868021396496747289677543653758101688499365196709600349831855157276274803, 4273006447766599851029197343910625305964779588947130545729882009677080459892767902139074897266046998378948193319329776150920821074305998905666368175737032487336505440974393061632257356697260043478850055918798196360043557336723402834256592984312957607731622934368169520289756716653736422872196593367101920308166527352079412344590114695, 1337615323422531101514598853478737615483725265103339849469231329210692205474781484172946466402355800190175297435923447189337380569240535257395618353403524999296122247241720498058900285176383928871155893061316083476812641934026019749574784965397208046520128690081232283708190247657782663414534596381371506499151016470118707739609022847, 6153442337491399463730666360630451953223069596225489087488144486112604050064388945200405772892415725900338512969888610774357818442533883395713661512978158729825590405567077341407532331228257558692149544940358814513930330757334104407683768945074261746314376336567331685071086776886141587264981001140113933152084084711732651014186070962, 3604640305526232611907645024550062043296875756707299105733076926046725900367109300556792007038226350415916164581084668832286944703572581980333135146219609212135443758296278947833224676930513359604019639508959092586550129343133801670578311053639151186604019506764426302411065588655525183584977810609190948096510237056027272167828516985, 2085633143403792178757363870459578017278494765962824809917819295576034993395441252397370078671500834434293581806342828037037877739374959864038090394888141183851040814572501861997433316552606005140700519616064369570511427453278726929921212633252656552135711332324923683242632912449218888750151556281762432994484741408165895663710356342]
X = 2**R_bits
m = len(N)

PR = PolynomialRing(ZZ, names=[str('x%d' % i) for i in range(1, 1 + m)])

h = 3
u = 1
variables = PR.gens()

gg = []
monomials = [variables[0]**0]
for i in range(m):
    gg.append(N[i] - variables[i])
    monomials.append(variables[i])

print('monomials:', monomials)

B = Matrix(ZZ, len(gg), len(monomials))
for ii in range(len(gg)):
    for jj in range(len(monomials)):
        if monomials[jj] in gg[ii].monomials():
            B[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj]([X] * m)

B = B.LLL()
print('-' * 32)

new_pol = []
for i in range(len(gg)):
    tmp_pol = 0
    for j in range(len(monomials)):
        tmp_pol += monomials[j](variables) * B[i, j] / monomials[j]([X] * m)
    new_pol.append(tmp_pol)

if len(new_pol) > 0:
    Ideal = ideal(new_pol[:m-1])
    GB = Ideal.groebner_basis()
    function_variables = var([str('y%d' % i) for i in range(1, 1 + m)])
    res = solve([pol(function_variables) for pol in GB], function_variables)

    print('got %d basis' % len(GB))
    print('solved result:')
    print(res)
    for tmp_res in res:
        PRRR.< x, y> = PolynomialRing(QQ)
        q = abs(PRRR(res[0][0](x, y)).coefficients()[0].denominator())
        p = N[-1] // q
        flag_tail = long_to_bytes(pow(C[-1], inverse(e, (p - 1) * (q - 1)), p * q))[-5:]
        flag = b''
        for i in range(m-1):
            then_res = PRRR(res[0][i](x, y))
            q = abs(then_res.coefficients()[0].numerator())
            flag += long_to_bytes(pow(C[i], inverse(e, (p-1)*(q-1)), p * q))[-5:]
        print(flag + flag_tail)

