import numpy as np
from main_sim import *
import matplotlib.pyplot as plt

def shearFlowGraph1(craft):
    #verification model shear flow, Sy = 1, Sz = 0
    q1Ver = np.array([ 1.34282076,  1.34241032,  1.34172636,  1.34076907,
        1.33953868,  1.3380355 ,  1.33625991,  1.33421235,  1.33189335,
        1.32930349,  1.32644341,  1.32331383,  1.31991556,  1.31624943,
        1.31231638,  1.30811739,  1.30365353,  1.2989259 ,  1.29393572,
        1.28868422,  1.28317274,  1.27740266,  1.27137543,  1.26509257,
        1.25855567,  1.25176636,  1.24472636,  1.23743744,  1.22990143,
        1.22212024,  1.21409581,  1.20583018,  1.19732542,  1.18858366,
        1.17960713,  1.17039806,  1.16095879,  1.15129168,  1.14139918,
        1.13128376,  1.12094799,  0.88590517,  0.87513653,  0.8641555 ,
        0.85296484,  0.84156737,  0.82996596,  0.81816353,  0.80616305,
        0.79396754,  0.78158007,  0.76900376,  0.75624178,  0.74329733,
        0.73017369,  0.71687414,  0.70340205,  0.68976079,  0.67595381,
        0.66198458,  0.64785662,  0.63357348,  0.61913877,  0.6045561 ,
        0.58982916,  0.57496166,  0.55995733,  0.54481995,  0.52955333,
        0.51416132,  0.4986478 ,  0.48301665,  0.46727183,  0.4514173 ,
        0.43545704,  0.41939508,  0.40323546,  0.38698224,  0.37063952,
        0.35421141,  0.33770205,  0.32111559, -0.05183435, -0.06856245,
       -0.08535507, -0.10220798, -0.11911693, -0.13607767, -0.15308593,
       -0.17013742, -0.18722787, -0.20435295, -0.22150836, -0.23868979,
       -0.2558929 , -0.27311337, -0.29034686, -0.30758903, -0.32483554])
    q2Ver = np.array([3.61411753, 3.61367894, 3.61294795, 3.61192457,
       3.6106088 , 3.60900062, 3.60710006, 3.6049071 , 3.60242174,
       3.59964399, 3.59657384, 3.5932113 , 3.58955636, 3.58560903,
       3.58136931, 3.57683719, 3.57201267, 3.56689576, 3.56148645,
       3.55578475, 3.54979066, 3.54350417, 3.53692528, 3.530054  ,
       3.52289033, 3.51543426, 3.5076858 , 3.49964494, 3.49131168,
       3.48268603, 3.47376799, 3.46455755, 3.45505472, 3.44525949,
       3.43517186, 3.42479185, 3.41411943, 3.40315462, 3.39189742,
       3.38034782, 3.36850583, 3.35637144, 3.34394466, 3.33122548,
       3.31821391, 3.30490994, 3.29131358, 3.27742482, 3.26324367,
       3.24877012, 3.23400418, 3.21894585, 3.20359511, 3.18795199,
       3.17201647, 3.15578855, 3.13926824, 3.12245553, 3.10535043,
       3.08795294, 3.07026305, 3.05228076, 3.03400608, 3.015439  ,
       2.99657953, 2.97742767, 2.95798341, 2.93824675, 2.91821771,
       2.89789626, 2.87728242, 2.85637619, 2.83517756, 2.81368653,
       2.79190311, 2.7698273 , 2.74745909, 2.72479849, 2.70184549,
       2.6786001 , 2.65506231, 2.63123212, 2.60710955, 2.58269457,
       2.55798721, 2.53298744, 2.50769528, 2.48211073, 2.45623378,
       2.43006444, 2.40360271, 2.37684857, 2.34980205, 2.32246312,
       2.29483181, 2.2669081 , 2.23869199, 2.21018349, 2.18138259])
    q3Ver = np.array([ 1.85654705,  1.81296295,  1.76982134,  1.7271222 ,  1.68486554,
        1.64305136,  1.60167965,  1.56075043,  1.52026368,  1.48021941,
        1.10791236,  1.06875305,  1.03003621,  0.99176186,  0.95392998,
        0.91654058,  0.87959365,  0.84308921,  0.80702724,  0.77140775,
        0.73623074,  0.70149621,  0.66720415,  0.63335458,  0.59994748,
        0.56698286,  0.26224733,  0.23016767,  0.19853048,  0.16733577,
        0.13658354,  0.10627379,  0.07640652,  0.04698172,  0.01799941,
       -0.01054043, -0.03863779, -0.06629267, -0.09350508, -0.12027501,
       -0.14660245, -0.17248742, -0.19792992, -0.43465145, -0.45920899,
       -0.48332405, -0.50699663, -0.53022673, -0.55301435, -0.5753595 ,
       -0.59726216, -0.61872235, -0.63974006, -0.6603153 , -0.68044805,
       -0.70013833, -0.71938612, -0.73819144, -0.75655429, -0.92570431,
       -0.9431822 , -0.9602176 , -0.97681053, -0.99296099, -1.00866896,
       -1.02393445, -1.03875747, -1.05313801, -1.06707607, -1.08057165,
       -1.09362476, -1.10623538, -1.11840353, -1.1301292 , -1.14141239,
       -1.2429909 , -1.25338914, -1.2633449 , -1.27285818, -1.28192898,
       -1.2905573 , -1.29874315, -1.30648652, -1.3137874 , -1.32064582,
       -1.32706175, -1.3330352 , -1.33856618, -1.34365468, -1.3483007 ,
       -1.35250424, -1.38651124, -1.38982982, -1.39270593, -1.39513956,
       -1.39713071, -1.39867938, -1.39978558, -1.4004493 , -1.40067054])
    q4Ver = np.array([-1.40067054, -1.4004493 , -1.39978558, -1.39867938, -1.39713071,
       -1.39513956, -1.39270593, -1.38982982, -1.38651124, -1.35250424,
       -1.3483007 , -1.34365468, -1.33856618, -1.3330352 , -1.32706175,
       -1.32064582, -1.3137874 , -1.30648652, -1.29874315, -1.2905573 ,
       -1.28192898, -1.27285818, -1.2633449 , -1.25338914, -1.2429909 ,
       -1.14141239, -1.1301292 , -1.11840353, -1.10623538, -1.09362476,
       -1.08057165, -1.06707607, -1.05313801, -1.03875747, -1.02393445,
       -1.00866896, -0.99296099, -0.97681053, -0.9602176 , -0.9431822 ,
       -0.92570431, -0.75655429, -0.73819144, -0.71938612, -0.70013833,
       -0.68044805, -0.6603153 , -0.63974006, -0.61872235, -0.59726216,
       -0.5753595 , -0.55301435, -0.53022673, -0.50699663, -0.48332405,
       -0.45920899, -0.43465145, -0.19792992, -0.17248742, -0.14660245,
       -0.12027501, -0.09350508, -0.06629267, -0.03863779, -0.01054043,
        0.01799941,  0.04698172,  0.07640652,  0.10627379,  0.13658354,
        0.16733577,  0.19853048,  0.23016767,  0.26224733,  0.56698286,
        0.59994748,  0.63335458,  0.66720415,  0.70149621,  0.73623074,
        0.77140775,  0.80702724,  0.84308921,  0.87959365,  0.91654058,
        0.95392998,  0.99176186,  1.03003621,  1.06875305,  1.10791236,
        1.48021941,  1.52026368,  1.56075043,  1.60167965,  1.64305136,
        1.68486554,  1.7271222 ,  1.76982134,  1.81296295,  1.85654705])
    q5Ver = np.array([3.61426373, 3.61411753, 3.61367894, 3.61294795, 3.61192457,
       3.6106088 , 3.60900062, 3.60710006, 3.6049071 , 3.60242174,
       3.59964399, 3.59657384, 3.5932113 , 3.58955636, 3.58560903,
       3.58136931, 3.57683719, 3.57201267, 3.56689576, 3.56148645,
       3.55578475, 3.54979066, 3.54350417, 3.53692528, 3.530054  ,
       3.52289033, 3.51543426, 3.5076858 , 3.49964494, 3.49131168,
       3.48268603, 3.47376799, 3.46455755, 3.45505472, 3.44525949,
       3.43517186, 3.42479185, 3.41411943, 3.40315462, 3.39189742,
       3.38034782, 3.36850583, 3.35637144, 3.34394466, 3.33122548,
       3.31821391, 3.30490994, 3.29131358, 3.27742482, 3.26324367,
       3.24877012, 3.23400418, 3.21894585, 3.20359511, 3.18795199,
       3.17201647, 3.15578855, 3.13926824, 3.12245553, 3.10535043,
       3.08795294, 3.07026305, 3.05228076, 3.03400608, 3.015439  ,
       2.99657953, 2.97742767, 2.95798341, 2.93824675, 2.91821771,
       2.89789626, 2.87728242, 2.85637619, 2.83517756, 2.81368653,
       2.79190311, 2.7698273 , 2.74745909, 2.72479849, 2.70184549,
       2.6786001 , 2.65506231, 2.63123212, 2.60710955, 2.58269457,
       2.55798721, 2.53298744, 2.50769528, 2.48211073, 2.45623378,
       2.43006444, 2.40360271, 2.37684857, 2.34980205, 2.32246312,
       2.29483181, 2.2669081 , 2.23869199, 2.21018349, 2.18138259])
    q6Ver = np.array([-0.32483554, -0.30758903, -0.29034686, -0.27311337, -0.2558929 ,
       -0.23868979, -0.22150836, -0.20435295, -0.18722787, -0.17013742,
       -0.15308593, -0.13607767, -0.11911693, -0.10220798, -0.08535507,
       -0.06856245, -0.05183435,  0.32111559,  0.33770205,  0.35421141,
        0.37063952,  0.38698224,  0.40323546,  0.41939508,  0.43545704,
        0.4514173 ,  0.46727183,  0.48301665,  0.4986478 ,  0.51416132,
        0.52955333,  0.54481995,  0.55995733,  0.57496166,  0.58982916,
        0.6045561 ,  0.61913877,  0.63357348,  0.64785662,  0.66198458,
        0.67595381,  0.68976079,  0.70340205,  0.71687414,  0.73017369,
        0.74329733,  0.75624178,  0.76900376,  0.78158007,  0.79396754,
        0.80616305,  0.81816353,  0.82996596,  0.84156737,  0.85296484,
        0.8641555 ,  0.87513653,  0.88590517,  1.12094799,  1.13128376,
        1.14139918,  1.15129168,  1.16095879,  1.17039806,  1.17960713,
        1.18858366,  1.19732542,  1.20583018,  1.21409581,  1.22212024,
        1.22990143,  1.23743744,  1.24472636,  1.25176636,  1.25855567,
        1.26509257,  1.27137543,  1.27740266,  1.28317274,  1.28868422,
        1.29393572,  1.2989259 ,  1.30365353,  1.30811739,  1.31231638,
        1.31624943,  1.31991556,  1.32331383,  1.32644341,  1.32930349,
        1.33189335,  1.33421235,  1.33625991,  1.3380355 ,  1.33953868,
        1.34076907,  1.34172636,  1.34241032,  1.34282076,  1.34295758])

    Sz = 0
    Sy = 1
    discret = Discretization(99,198,99,198)
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4 = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, discret.n1,discret.n2,discret.n3,discret.n4)

    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],-q5Ver))

    relError1 = ((q1-q1VerM)/q1VerM) * 100
    relError2 = ((q2-q2VerM)/q2VerM) * 100
    relError3 = ((q3-q3VerM)/q3VerM) * 100
    relError4 = ((q4-q4VerM)/q4VerM) * 100

    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=1,Sz=0)", fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=1,Sz=0)", fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1,label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2, label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3, label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4, label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Relative error [%]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Relative error [%]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.legend()
    plt.show()
    return

def shearFlowGraph2(craft):
    #verification model shear flow, Sy = 0, Sz = 1
    q1Ver = np.array([ -0.07218748, -0.07835894, -0.08452878, -0.09069619,
       -0.09686037, -0.1030205 , -0.10917577, -0.11532538, -0.12146854,
       -0.12760442, -0.13373224, -0.13985119, -0.14596048, -0.15205931,
       -0.15814689, -0.16422244, -0.17028516, -0.17633427, -0.182369  ,
       -0.18838857, -0.19439221, -0.20037915, -0.20634862, -0.21229988,
       -0.21823215, -0.2241447 , -0.23003678, -0.23590765, -0.24175657,
       -0.24758282, -0.25338567, -0.25916441, -0.26491833, -0.27064671,
       -0.27634887, -0.28202411, -0.28767174, -0.29329109, -0.29888149,
       -0.30444227, -0.30997277, -0.43329212, -0.43876014, -0.44419597,
       -0.44959898, -0.45496855, -0.46030407, -0.46560495, -0.4708706 ,
       -0.47610043, -0.48129387, -0.48645036, -0.49156934, -0.49665027,
       -0.50169261, -0.50669584, -0.51165943, -0.5165829 , -0.52146573,
       -0.52630744, -0.53110757, -0.53586563, -0.54058119, -0.5452538 ,
       -0.54988301, -0.55446842, -0.55900962, -0.56350619, -0.56795776,
       -0.57236394, -0.57672438, -0.58103871, -0.5853066 , -0.58952772,
       -0.59370173, -0.59782835, -0.60190727, -0.60593821, -0.60992091,
       -0.61385509, -0.61774051, -0.62157694, -0.70641934, -0.71015714,
       -0.71384532, -0.71748369, -0.72107209, -0.72461035, -0.72809833,
       -0.73153589, -0.73492291, -0.73825929, -0.74154492, -0.74477972,
       -0.74796362, -0.75109657, -0.75417851, -0.75720942, -0.76018928])
    q2Ver = np.array([ -0.00495843, -0.00991687, -0.0148753 , -0.01983373,
       -0.02479216, -0.0297506 , -0.03470903, -0.03966746, -0.0446259 ,
       -0.04958433, -0.05454276, -0.0595012 , -0.06445963, -0.06941806,
       -0.07437649, -0.07933493, -0.08429336, -0.08925179, -0.09421023,
       -0.09916866, -0.10412709, -0.10908553, -0.11404396, -0.11900239,
       -0.12396082, -0.12891926, -0.13387769, -0.13883612, -0.14379456,
       -0.14875299, -0.15371142, -0.15866986, -0.16362829, -0.16858672,
       -0.17354515, -0.17850359, -0.18346202, -0.18842045, -0.19337889,
       -0.19833732, -0.20329575, -0.20825419, -0.21321262, -0.21817105,
       -0.22312948, -0.22808792, -0.23304635, -0.23800478, -0.24296322,
       -0.24792165, -0.25288008, -0.25783851, -0.26279695, -0.26775538,
       -0.27271381, -0.27767225, -0.28263068, -0.28758911, -0.29254755,
       -0.29750598, -0.30246441, -0.30742284, -0.31238128, -0.31733971,
       -0.32229814, -0.32725658, -0.33221501, -0.33717344, -0.34213188,
       -0.34709031, -0.35204874, -0.35700717, -0.36196561, -0.36692404,
       -0.37188247, -0.37684091, -0.38179934, -0.38675777, -0.39171621,
       -0.39667464, -0.40163307, -0.4065915 , -0.41154994, -0.41650837,
       -0.4214668 , -0.42642524, -0.43138367, -0.4363421 , -0.44130054,
       -0.44625897, -0.4512174 , -0.45617583, -0.46113427, -0.4660927 ,
       -0.47105113, -0.47600957, -0.480968  , -0.48592643, -0.49088487])
    q3Ver = np.array([-1.25107415e+00, -1.25841826e+00, -1.26544351e+00, -1.27214989e+00,
       -1.27853741e+00, -1.28460607e+00, -1.29035586e+00, -1.29578678e+00,
       -1.30089885e+00, -1.30569204e+00, -1.34725909e+00, -1.35141455e+00,
       -1.35525116e+00, -1.35876890e+00, -1.36196778e+00, -1.36484779e+00,
       -1.36740893e+00, -1.36965122e+00, -1.37157464e+00, -1.37317919e+00,
       -1.37446488e+00, -1.37543170e+00, -1.37607967e+00, -1.37640876e+00,
       -1.37641899e+00, -1.37611036e+00, -1.36898314e+00, -1.36803678e+00,
       -1.36677155e+00, -1.36518746e+00, -1.36328451e+00, -1.36106269e+00,
       -1.35852201e+00, -1.35566246e+00, -1.35248405e+00, -1.34898678e+00,
       -1.34517064e+00, -1.34103563e+00, -1.33658176e+00, -1.33180903e+00,
       -1.32671743e+00, -1.32130697e+00, -1.31557764e+00, -1.25943729e+00,
       -1.25307023e+00, -1.24638432e+00, -1.23937953e+00, -1.23205588e+00,
       -1.22441337e+00, -1.21645200e+00, -1.20817175e+00, -1.19957265e+00,
       -1.19065468e+00, -1.18141785e+00, -1.17186215e+00, -1.16198758e+00,
       -1.15179416e+00, -1.14128187e+00, -1.13045071e+00, -1.02561609e+00,
       -1.01414720e+00, -1.00235946e+00, -9.90252842e-01, -9.77827365e-01,
       -9.65083023e-01, -9.52019817e-01, -9.38637746e-01, -9.24936811e-01,
       -9.10917012e-01, -8.96578348e-01, -8.81920820e-01, -8.66944428e-01,
       -8.51649171e-01, -8.36035050e-01, -8.20102064e-01, -6.66573177e-01,
       -6.50002463e-01, -6.33112884e-01, -6.15904441e-01, -5.98377134e-01,
       -5.80530962e-01, -5.62365926e-01, -5.43882026e-01, -5.25079261e-01,
       -5.05957632e-01, -4.86517139e-01, -4.66757781e-01, -4.46679559e-01,
       -4.26282473e-01, -4.05566522e-01, -3.84531706e-01, -1.82308553e-01,
       -1.60636009e-01, -1.38644601e-01, -1.16334329e-01, -9.37051915e-02,
       -7.07571902e-02, -4.74903245e-02, -2.39045944e-02,  5.55111512e-17])
    q4Ver = np.array([5.55111512e-17, 2.39045944e-02, 4.74903245e-02, 7.07571902e-02,
       9.37051915e-02, 1.16334329e-01, 1.38644601e-01, 1.60636009e-01,
       1.82308553e-01, 3.84531706e-01, 4.05566522e-01, 4.26282473e-01,
       4.46679559e-01, 4.66757781e-01, 4.86517139e-01, 5.05957632e-01,
       5.25079261e-01, 5.43882026e-01, 5.62365926e-01, 5.80530962e-01,
       5.98377134e-01, 6.15904441e-01, 6.33112884e-01, 6.50002463e-01,
       6.66573177e-01, 8.20102064e-01, 8.36035050e-01, 8.51649171e-01,
       8.66944428e-01, 8.81920820e-01, 8.96578348e-01, 9.10917012e-01,
       9.24936811e-01, 9.38637746e-01, 9.52019817e-01, 9.65083023e-01,
       9.77827365e-01, 9.90252842e-01, 1.00235946e+00, 1.01414720e+00,
       1.02561609e+00, 1.13045071e+00, 1.14128187e+00, 1.15179416e+00,
       1.16198758e+00, 1.17186215e+00, 1.18141785e+00, 1.19065468e+00,
       1.19957265e+00, 1.20817175e+00, 1.21645200e+00, 1.22441337e+00,
       1.23205588e+00, 1.23937953e+00, 1.24638432e+00, 1.25307023e+00,
       1.25943729e+00, 1.31557764e+00, 1.32130697e+00, 1.32671743e+00,
       1.33180903e+00, 1.33658176e+00, 1.34103563e+00, 1.34517064e+00,
       1.34898678e+00, 1.35248405e+00, 1.35566246e+00, 1.35852201e+00,
       1.36106269e+00, 1.36328451e+00, 1.36518746e+00, 1.36677155e+00,
       1.36803678e+00, 1.36898314e+00, 1.37611036e+00, 1.37641899e+00,
       1.37640876e+00, 1.37607967e+00, 1.37543170e+00, 1.37446488e+00,
       1.37317919e+00, 1.37157464e+00, 1.36965122e+00, 1.36740893e+00,
       1.36484779e+00, 1.36196778e+00, 1.35876890e+00, 1.35525116e+00,
       1.35141455e+00, 1.34725909e+00, 1.30569204e+00, 1.30089885e+00,
       1.29578678e+00, 1.29035586e+00, 1.28460607e+00, 1.27853741e+00,
       1.27214989e+00, 1.26544351e+00, 1.25841826e+00, 1.25107415e+00])
    q5Ver = np.array([ 0       , -0.00495843, -0.00991687, -0.0148753 , -0.01983373,
       -0.02479216, -0.0297506 , -0.03470903, -0.03966746, -0.0446259 ,
       -0.04958433, -0.05454276, -0.0595012 , -0.06445963, -0.06941806,
       -0.07437649, -0.07933493, -0.08429336, -0.08925179, -0.09421023,
       -0.09916866, -0.10412709, -0.10908553, -0.11404396, -0.11900239,
       -0.12396082, -0.12891926, -0.13387769, -0.13883612, -0.14379456,
       -0.14875299, -0.15371142, -0.15866986, -0.16362829, -0.16858672,
       -0.17354515, -0.17850359, -0.18346202, -0.18842045, -0.19337889,
       -0.19833732, -0.20329575, -0.20825419, -0.21321262, -0.21817105,
       -0.22312948, -0.22808792, -0.23304635, -0.23800478, -0.24296322,
       -0.24792165, -0.25288008, -0.25783851, -0.26279695, -0.26775538,
       -0.27271381, -0.27767225, -0.28263068, -0.28758911, -0.29254755,
       -0.29750598, -0.30246441, -0.30742284, -0.31238128, -0.31733971,
       -0.32229814, -0.32725658, -0.33221501, -0.33717344, -0.34213188,
       -0.34709031, -0.35204874, -0.35700717, -0.36196561, -0.36692404,
       -0.37188247, -0.37684091, -0.38179934, -0.38675777, -0.39171621,
       -0.39667464, -0.40163307, -0.4065915 , -0.41154994, -0.41650837,
       -0.4214668 , -0.42642524, -0.43138367, -0.4363421 , -0.44130054,
       -0.44625897, -0.4512174 , -0.45617583, -0.46113427, -0.4660927 ,
       -0.47105113, -0.47600957, -0.480968  , -0.48592643, -0.49088487])
    q6Ver = np.array([0.76018928, 0.75720942, 0.75417851, 0.75109657, 0.74796362,
       0.74477972, 0.74154492, 0.73825929, 0.73492291, 0.73153589,
       0.72809833, 0.72461035, 0.72107209, 0.71748369, 0.71384532,
       0.71015714, 0.70641934, 0.62157694, 0.61774051, 0.61385509,
       0.60992091, 0.60593821, 0.60190727, 0.59782835, 0.59370173,
       0.58952772, 0.5853066 , 0.58103871, 0.57672438, 0.57236394,
       0.56795776, 0.56350619, 0.55900962, 0.55446842, 0.54988301,
       0.5452538 , 0.54058119, 0.53586563, 0.53110757, 0.52630744,
       0.52146573, 0.5165829 , 0.51165943, 0.50669584, 0.50169261,
       0.49665027, 0.49156934, 0.48645036, 0.48129387, 0.47610043,
       0.4708706 , 0.46560495, 0.46030407, 0.45496855, 0.44959898,
       0.44419597, 0.43876014, 0.43329212, 0.30997277, 0.30444227,
       0.29888149, 0.29329109, 0.28767174, 0.28202411, 0.27634887,
       0.27064671, 0.26491833, 0.25916441, 0.25338567, 0.24758282,
       0.24175657, 0.23590765, 0.23003678, 0.2241447 , 0.21823215,
       0.21229988, 0.20634862, 0.20037915, 0.19439221, 0.18838857,
       0.182369  , 0.17633427, 0.17028516, 0.16422244, 0.15814689,
       0.15205931, 0.14596048, 0.13985119, 0.13373224, 0.12760442,
       0.12146854, 0.11532538, 0.10917577, 0.1030205 , 0.09686037,
       0.09069619, 0.08452878, 0.07835894, 0.07218748, 0.0660152 ])

    Sz = 1
    Sy = 0
    discret = Discretization(99,198,99,198)
    q1,q2,q3,q4,sVec1,sVec2,sVec3,sVec4 = calcShFlow(craft.ha,craft.ca,craft.tsk,craft.tsp,craft.tst,craft.hst,craft.wst,craft.nst,Sz,Sy, discret.n1,discret.n2,discret.n3,discret.n4)

    q4VerM = np.concatenate((-q1Ver[::-1],-q6Ver[::-1]))
    q1VerM = (-q3Ver[::-1])
    q3VerM = (-q4Ver[::-1])
    q2VerM = np.concatenate((-q2Ver[::-1],q5Ver))

    relError1 = np.divide((q1-q1VerM), q1VerM, out=np.zeros_like((q1-q1VerM)), where=q1VerM != 0) * 100
    relError2 = np.divide((q2 - q2VerM), q2VerM, out=np.zeros_like((q2 - q2VerM)), where=q2VerM != 0) * 100
    relError3 = np.divide((q3 - q3VerM), q3VerM, out=np.zeros_like((q3 - q3VerM)), where=q3VerM != 0) * 100
    relError4 = np.divide((q4 - q4VerM), q4VerM, out=np.zeros_like((q4 - q4VerM)), where=q4VerM != 0) * 100
    #relError1 = ((q1-q1VerM)/q1VerM) * 100
    #relError2 = ((q2-q2VerM)/q2VerM) * 100
    #relError3 = ((q3-q3VerM)/q3VerM) * 100
    #relError4 = ((q4-q4VerM)/q4VerM) * 100
    print(q1[1])
    print(q1VerM[1])
    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Relative error VS distance in s-direction (Sy=0,Sz=1)", fontsize=16)

    axs[0, 0].plot(sVec1, relError1)
    axs[0, 0].set_title('Region 1')
    axs[0,0].grid()
    axs[0,0].set_ylim(-1,1)

    axs[0, 1].plot(sVec2, relError2)
    axs[0, 1].set_title('Region 2')
    axs[0,1].grid()
    axs[0, 1].set_ylim(-5, 5)

    axs[1, 0].plot(sVec3, relError3)
    axs[1, 0].set_title('Region 3')
    axs[1,0].grid()
    axs[1, 0].set_ylim(-1, 1)

    axs[1, 1].plot(sVec4, relError4)
    axs[1, 1].set_title('Region 4')
    axs[1,1].grid()
    axs[0,0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,0].set_xlabel('s [m]',fontsize = 12.0)
    axs[1, 0].set_ylabel('Relative error [%]',fontsize = 12.0)
    axs[1,1].set_xlabel('s [m]',fontsize = 12.0)

    axs[0,0].tick_params(axis='both', which='major', labelsize=12)
    axs[0,0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.show()


    fig, axs = plt.subplots(2, 2)
    plt.suptitle("Shear flow VS distance in s-direction (Sy=0,Sz=1)", fontsize=16)
    plt.rcParams["font.size"] = "12"
    plt.rcParams["axes.labelsize"] = "12"
    axs[0,0].plot(sVec1,q1,label="Numerical")
    axs[0,0].plot(sVec1,q1VerM,label="Verification")
    axs[0,0].set_title('Region 1')
    axs[0, 0].grid()

    axs[0, 1].plot(sVec2, q2, label="Numerical")
    axs[0, 1].plot(sVec2, q2VerM, label="Verification")
    axs[0, 1].set_title('Region 2')
    axs[0, 1].grid()

    axs[1, 0].plot(sVec3, q3, label="Numerical")
    axs[1, 0].plot(sVec3, q3VerM, label="Verification")
    axs[1, 0].set_title('Region 3')
    axs[1, 0].grid()

    axs[1, 1].plot(sVec4, q4, label="Numerical")
    axs[1, 1].plot(sVec4, q4VerM, label="Verification")
    axs[1, 1].set_title('Region 4')
    axs[1, 1].grid()

    axs[0, 0].set_ylabel('Relative error [%]', fontsize=12.0)
    axs[1, 0].set_xlabel('s [m]', fontsize=12.0)
    axs[1, 0].set_ylabel('Relative error [%]', fontsize=12.0)
    axs[1, 1].set_xlabel('s [m]', fontsize=12.0)

    axs[0, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[0, 1].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 0].tick_params(axis='both', which='minor', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='major', labelsize=12)
    axs[1, 1].tick_params(axis='both', which='minor', labelsize=12)
    plt.legend()
    plt.show()
    return

