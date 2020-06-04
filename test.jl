ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
using PyPlot
using PyCall
@pyimport matplotlib.pyplot as pyplt
using VoronoiDelaunay
sp = pyimport("scipy.spatial")

function draw_an_i_cell(i_cell::Array{Float64})
    x_pos = i_cell[1]
    y_pos = i_cell[2]
    radius= i_cell[3]
    circle = pyplt.Circle((x_pos, y_pos), radius, alpha=0.1 , color=:green, lw=0.3)
    return circle
end

function main()
    cluster::Array{Float64} = [4.029966192209636 -24.235622245175954 4.555829122323626; 50.347088866538336 -13.611809643785415 4.492776500152409; 24.208947111778123 -18.88937379034569 3.578333192307571;
        28.38317356883336 3.297257180872464 3.6430265861049587; -15.057990818610167 -20.078342111767103 4.518997719847969; -24.373116057970712 -7.880064285757995 3.5662842812237963; -4.887588900847609 -28.113390474784218 4.904353145360522;
        22.217305221345026 -24.66106336202428; 4.6756110739497405;21.455506389107608 -31.267355860717373 3.9149677413197197; -12.326044056970936 -10.540403475618195 3.5578022595436796;
        -7.734072830254266 -54.16214437215862 4.1272978223737065; -34.70109716839603 10.554126332327321 3.893402364925022; 23.629898028882963 -10.139149138195437 4.386622671242524;
        -4.187816010547823 23.500327549958286 4.98605853809991; 23.05158393864919 -17.438722060260577 4.834849002834547; 1.4392390336935907 -46.19987122119058 4.874328802055369;
        -18.76573500773088 -41.73726441435922 3.600736738336242; -2.528665536776914 -26.308941165296556 4.494434405841675; -0.5106103364410826 -15.283004050286982 4.5363653761654845;
        31.98570033739687 -7.699503996106661 3.92229172916434; 15.227061236262522 -22.437995769452026 4.618723996676776; -1.8226283473727012 4.629162929799875 4.992054827092701;
        43.257521598458894 -2.8334842340654554 4.311297118916599; 8.880472901697216 -22.72418171642875 4.8831168928756234; 6.410208292693223 -45.001640755842445 4.997584539862925;
        10.664196591174465 -4.4255659518269725 4.182850960161763; 1.5653930033452406 -10.431025602994 4.876481219836166; -19.389365107308514 -1.6110258780831654 4.977289869418323;
        0.11535362974912752 -2.026694341359384 4.812512833817315; -8.798996184194124 -33.79523925578274 4.7820001449477765; 22.998558644509316 17.066833338950985 3.6744358334197837;
         28.434134778677656 -38.15168080828297 4.782271679048689; -25.292255591668933 -29.85348566756057 4.91888546653187; 5.37271363457914 -13.684039627603537 4.566088138274815;
         39.84511030442959 -9.205575466124177 4.182944003620574; 21.838130181354703 -14.609748192931063 4.384579755420322; 7.374247595061265 -36.082038952568894 4.381233173125821;
          19.162252191288232 -44.28950597572376 4.618562935329358; -17.748739266982167 3.335305958840007 4.770561774522136; -7.385415651168388 -28.553208502723496 4.342012377946507;
           32.84816425217077 -7.418762784048873 3.9276141650096217; 14.255094497693875 -12.911644711370547 4.4181752246880155; 21.888852418265017 3.569342703527516 4.12430262868386;
            -9.814164089798652 15.476802418123457 4.094436845043389; -18.916245953632405 -38.80305666988572 3.8786559767199424; -18.42349405683296 28.853045602398605 4.29862905785202;
             -21.886243390702024 4.973782585436631 3.634650236609859; 17.594305162343062 -16.60919903291428 3.6205616888373378; -29.13060987491825 -9.21265643272673 3.6308587673170765;
             15.45029222146308 -19.365696838276076 3.9557197523471026; -21.3264023699128 -41.77506255685126 3.6483919145201575; 0.9605460158736332 1.7338592547354523 3.6894867860911145;
              6.552013513835905 -34.35158139670638 3.6566312537581362; 8.161372228410439 12.720315546914165 4.534043037807636; -11.507827587327009 -14.061942313470656 4.070884103682171;
              -3.6605217133743557 5.368790437033393 4.683135814939198; 4.1940679265090415 -14.580632903533228 4.179523216485913; 21.104577356379828 -21.75138152436704 4.47932104873039;
               -14.418496181216389 -20.38384950013092 4.717340483358165; -1.2975057904107496 -21.451841355686014 4.272470385654995; 28.527456558113233 -31.001323385876365 4.500742374760389;
                29.90866196135948 -21.540009190241037 4.5667522180445514; 17.145444751598102 -5.363211670939717 3.798085205132899; 0.022395052022904656 0.8202360237863595 4.809751178095939;
                 41.71008075008411 -12.180744035438217 4.3504819353699276; 5.5113502742160865 -7.0181422162862885 4.04859676308962; 21.620345240016302 5.550860310713909 4.375973964664479;
                  28.46658402752934 2.699623273390648 4.235163566433475; 2.72602489133173 -10.656005709134297 4.894583665771015; -15.497060767650597 -2.570122990125118 3.8747339394689844;
                  -7.914201830155655 -7.3273944595503835 4.49027142092459; 24.903411311101245 -23.21419207010225 4.162682269016984; 32.87910460175724 8.87991340785395 4.956138955052321;
                  -1.0598108517060456 -42.75957533384135 4.689700733341145; 4.807867050608603 -48.607869451345174 4.01910546047683; 37.78247582035044 -26.40514099443926 4.388986434195438;
                   34.563825084986426 6.218338017832526 4.751668852661795; 48.530530121120044 -11.84174125796154 4.621396482891163; 26.90657861252012 0.9381371412430773 4.1049613098062165;
                    28.585654815936035 14.30272702647767 4.236698867221782; -8.256635583627004 -20.832910038131434 4.048551368011697; 31.63532049494433 -1.6189453220604013 4.115321815732534;
                     1.331662746998606 -12.268780703907401 3.866679736418493; 6.503098544003463 0.5479991553301152 3.94551155571875; 7.607625977046018 -24.65841008749908 4.259093851330552;
                      1.3077481441670862 -41.297584589626474 4.499159144422204; 38.62529027319201 -26.9068429846219 3.7657931993141536; 27.801200588479194 -8.348389927431114 3.9299891302772614;
                      -12.26759980758789 -38.46763967662847 4.221756616138911; -21.033210409720848 19.1599509517179 4.122341517983498; -3.7755888275654987 10.091275208142651 3.70361745669897;
                      -18.213898531493385 -24.743849366563197 3.572529395273718; -27.024057660499896 -12.363839281241392 3.8206788395086964; -17.393485675553034 -29.8696339812102 4.268715427666669;
                      3.3037645757911993 -7.334009792595194 3.7618150217234056; 4.177902450105164 -8.974023274603445 4.879836162519129; 15.390303470414432 -9.717827971379332 4.906666278590004;
                       -5.535741246289203 -27.20951309986036 4.8495161455181774; 26.91320635719302 8.728000724303913 3.771856409798692; 1.2103478185273269 -32.45789317386747 3.5897129663826277;
                        23.871255609750474 23.582361130924323 4.810909077454327; 2.588454136973457 -34.84242234544557 3.681909223575819; 10.524212840728186 -6.212353303983394 3.5899364537107554;
                        5.418603388982019 -3.005477984839357 3.989911270811795; 27.785196616155517 -17.96020022086642 3.6269744542041833; -15.812590819922823 -2.2685067338808933 3.797717448248255;
                         13.675852295754519 -17.699524778700994 4.567108431491967; 11.59766391582035 -18.636523805401236 3.5355339059327373; -3.1972696735344335 24.583972386352595 4.657582316212511;
                          12.014038627210455 -33.31128079473813 4.9224505392309; -12.744184201909764 -54.02750983511864 3.7708299146136675; -1.3745755009676142 0.32220544314905536 4.826791529156798;
                          13.25999360071263 -33.6605122962531 4.846464901948001; 1.8334046020924388 8.617080249711538 4.839497479999835; -8.343006048876925 -28.65216190031076 4.729066155289801;
                           5.944326359502884 10.432799053073369 3.5585841159472285; -16.773174451078216 -15.772994792651605 4.960446871933323; 6.255300386107505 20.38968745012692 4.6301889109468615;
                            9.209093473014685 -23.416036629247376 4.972215563555721; 11.827923541308769 -32.91943123233946 3.5355339059327373; 18.978145541242117 -43.151073179229385 4.42005716675757;
                             7.316965119164588 -10.582102393953729 4.532478422444256; -14.297588213701973 -15.133028947857566 3.640609026523554; 6.160063936296565 -5.540190331452347 4.690881938274542;
                              14.451256761315214 -41.309828798609615 4.61672720047058; -1.7923992762358592 12.817155023079149 4.514360212303638; 38.628809595723126 -23.957886099274454 4.347552121002803;
                              35.59245980751547 22.155454136011198 4.770691337508814; 2.6929702148754378 -19.164032777260083 4.8777246048825145;2.346505212503763 -15.04571868217956 4.690171682801096;
                               15.43303051663189 -4.937843318249347 4.766316498800751; 15.994595379828365 -22.665380970283135 4.331104462773446; 32.49160119963068 -30.287114930699023 4.569862737683683;
                               5.1283382020911965 -44.922470055218355 4.558616170871972; -10.219892297841437 -8.013539216251367 4.181315152946079; -5.529390074141899 -22.578195469813828 4.760518550949003;
                                -3.6363393010270895 -29.78348961638505 4.5923822525254705; 12.634236429279783 0.347679556771148 4.432975649960139; 18.687828088749885 -4.697269296746138 4.520358549916619;
                                -29.045832122938425 -27.82827364772563 4.52485740537866; 7.158499097037244 -33.6996639239853 4.806462698379949; 28.430065777996678 -31.89776494236193 4.4248541644111725;
                                28.850389870779658 -18.457041149400368 4.577698318962922; 18.489604982076745 -5.799769089294766 4.643436096248969; 48.31370380339915 -11.49335199805605 4.436792091119345;
                                 6.333202967763874 -11.521490062899222 4.496312851338547; 25.992560255178347 -9.28561171254019 4.798203324696759; 42.63627445424992 -14.847118202599688 4.619638274910993;
                                 -11.620679549929205 -19.13167699871712 4.4179561658327176; 22.21104406287517 -18.801651156960197 4.52900717368697; 15.547567303716697 -17.80476121600689 4.607894690025552;
                                  -19.296909268226443 25.472816017934854 4.303966951102513; 27.958377920343388 1.9536243514952056 4.430630379774139; 38.379088292329236 -6.838349645178887 4.292745380486154;
                                   40.248303211615166 -17.43174377654395 4.390755272737553; 23.9428240488954 -18.20486045563731 4.207683397216301;31.962476387582313 -6.942308719387139 4.045822028775442;
                                   -0.5256481715233867 -22.474928978175836 4.114719858298652; 2.732718376964927 -26.38746687757149 4.131857495936751; 16.603190462881646 -7.000087693309547 4.433125919711209;
                                    10.285650374858388 -54.819103310122216 4.442498452500032; -3.0429760268424486 -53.14306614901409 4.381928262687933; 18.733939697419338 -1.3288328609309432 4.009436438428586;
                                     -8.205901216830375 -20.06970387152325 4.21937627933251; 8.643272391476478 12.962000507126035 4.0422037694528; 21.735178482593646 -3.2238885827598796 3.83872912047904;
                                     -3.2148723490301143 -2.8777568825051656 3.9264923319945733; 17.681652168115917 -0.9220862774924355 3.9845514936955415; -9.741599609212216 -47.77500260949072 3.849618995913424;
                                      10.857972749442562 -19.358032176519618 4.027215588856022; 38.91460804921429 -6.675001398874645 3.911699770153965; -22.633174125161396 14.297618313677741 3.9998615089238236;
                                       -13.534003145584743 -32.89952594925196 3.929760115837831; -33.51503108378137 22.445583629560648 4.0364001603583235; 6.183252194946779 -13.579336493427931 4.0527726799522;
                                        -7.57327676375003 -7.737155777982903 3.911737597148557; 9.894774208253791 -16.4809264750813 3.9494578531146805; -26.82018649390868 -36.93301573048284 3.7707109055944636;
                                         28.176744839218266 -30.134084007442326 3.9722173161309913; 7.919342342076702 3.2787543618282955 3.8610553833579146; 18.778082349805242 -7.471495419050877 3.7522920949239578;
                                         1.0374453604649254 -2.5919099044716534 3.8217939782927117; 29.425390530179154 -7.543293293772885 3.833274145983436; -16.198833254893533 2.2436429363880963 3.945005230952917;
                                         -24.582637988661663 -10.676639282932776 3.8196726436033885; -20.574153870452907 -44.81653964125077 3.6990312474460154; 0.04878870881266742 9.255831834913108 3.9612574546962214;
                                          -14.983633943715583 -15.347371714503357 3.7376073907747958; 37.34502657460135 -33.51339281207573 3.7003223864846206; -22.389183242319962 6.735551059174586 3.7416089760440374;
                                          -26.22326979007686 -44.18157378240172 3.6542752750052117; 26.134087607018888 12.71944877491571 3.8949372175115418; 31.88838384272797 5.377288184813411 3.6011576874021447;
                                           -0.9233066658447696 -4.220167111219 3.626203852666109; 12.46173164975903 -30.48572031348784 3.6185336878602152; -33.14969632337744 -8.828784765806128 3.6506548676334423;
                                           -16.552038670498025 -17.812424419894874 3.6023516944605287; -19.545338278608494 -25.454378549351592 3.6276951909383293; 20.558565026852392 -18.118926516292404 3.6821872488611023;
                                            -0.6430527538617317 -32.27885030496037 3.581644398082727; 0.9940998314310713 -34.80473155108466 3.6915311951117906; 9.57977408517352 -32.04717657431207 3.5487829125737953;
                                            30.455893112412586 -15.604553844416793 3.6997584824705556; -19.302305996122683 -37.98569792159171 3.6532799332198462; 23.43031025112344 -18.57381235091496 3.6127919924090275;
                                            28.588059762778723 3.9887905241398456 3.6444817672331906; 11.254462340953992 -7.206775200806682 3.5355339059327373; 8.843894195724378 12.38154552668474 3.5363270620900065;
                                            -10.848176625715599 -10.261102884427224 3.5355339059327373; -23.66979121585819 -8.62079066774822 3.5478172599657958; 11.400737723704719 -18.993924400059687 3.5355339059327373]
    plt.clf()
    ax=pyplt.gca()
    ax.set_xlim((-50, 50))
    ax.set_ylim((-50, 50))
    ax.set_aspect("equal")

    for k::Int64 in (1:1:size(cluster,2))
        cell_to_draw = cluster[:,k]
        ax.add_artist(draw_an_i_cell(cell_to_draw))
    end

    x_coords = cluster[1, :] #only CoM coords
    y_coords = cluster[2, :] #only CoM coords
    points = hcat(x_coords, y_coords)
    triang = dela(points)
    plt.triplot(points[:,1], points[:,2], triang.simplices, color = :blue)
    plt.plot(points[:,1], points[:,2], lw= 0, marker="o", linestyle="")
    fig = gcf()
    display(fig)
end
