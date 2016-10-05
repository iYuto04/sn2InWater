import scipy as sp
import numpy as np
import copy
import math
import subprocess

def cal_ele(ab, bc,zi):
    r = ab*ab - bc*bc
    print("reaction cordinate r is", r)
    alpha = np.zeros(3)
    beta = np.zeros(3)
    gamma = np.zeros(3)

    alpha  = [-4.3688, 0.0, 4.3688]
    beta   = [1.1255, -2.2510, 1.1255]
    gamma = [-11.3645, 4.5020, -11.3645]
    k1 = 0.75
    k2 = 0.25
    print("==================================")
    print("Cl:",  0.05492*(-1.0*math.pi*alpha[0]/2.0 + gamma[0]))
    print("CH3:", 0.05492*(-1.0*math.pi*alpha[1]/2.0 + gamma[1]))
    #print("kakunin", gamma[1])
    print("Cl-:", 0.05492*(-1.0*math.pi*alpha[2]/2.0 + gamma[2]))
    print("==================================")

    for i in range(3):
        #if i == 0:
            #print("alpha[i]*arctan = ",alpha[i]*math.atan(k1*r))
            #print("beta[i]/(k2*r*r +1) =", (alpha[i]*math.atan(k1*r)+beta[i]/(k2*r*r +1.0) + gammma[i])*0.06829)
        zi[i] = (alpha[i]*math.atan(k1*r) + beta[i]/(k2*r*r +1.0) + gamma[i])*0.05492
        print(zi[i])
    return zi


def write1_ross(r_ab,r_bc,r_ac,zi):
    f.write("SN2 REACTION IN TIP3P WATER\n")
    f.write("sn2.out\n")
    f.write("../../exvv/spc/uv.dat\n")
    f.write("restart.dat\n")
    f.write("gr.mat\n")
    f.write("../../exvv/sn2_tip3p/dtuv.dat\n")
    f.write("../../exvv/sn2_tip3p/dduv.dat\n")
    f.write("resdt.dat\n")
    f.write("resdd.dat\n")
    f.write("dthuv.mat\n")
    f.write("ddhuv.mat\n")
    f.write("pmf.mat\n")
    f.write("3, 3\n")
    f.write("\'meth\', \'cl1  \', \'cl2 \'\n")
    f.write("\'O   \', \'H   \', \'H   \'\n")
    f.write(str(r_ab))
    f.write(", ")
    f.write(str(r_bc))
    f.write(", ")
    f.write(str(r_ac))
    f.write("\n")
    f.write("0.9572, 0.9572, 1.5139\n")
    f.write("1.908, 2.2890, 2.2890\n")
    f.write("1.575, 0.2, 0.2\n")
    f.write("0.75006d-15, 0.12689d-14, 0.12689d-14\n")
    f.write("1.0605d-14, 3.2093d-15, 3.2093d-15\n")
    #f.write("0.24, -0.25, -0.9\n")
    # for i in range(3):
    #     print(zi[i])
    f.write(str(zi[1]))
    f.write(",")
    f.write(str(zi[0]))
    f.write(",")
    f.write(str(zi[2]))
    f.write("\n")
    f.write("-0.834, 0.417, 0.417\n")
    f.write("298.15, 0.03334\n")
    f.write("1.0\n")

def write1_hira(r_ab,r_bc,r_ac,zi):
    f.write("SN2 REACTION IN TIP3P WATER\n")
    f.write("tip3p.out\n")
    f.write("../../exvv/sn2_tip3p/uv.dat\n")
    f.write("restart.dat\n")
    f.write("gr.mat\n")
    f.write("../../exvv/sn2_tip3p/dtuv.dat\n")
    f.write("../../exvv/sn2_tip3p/dduv.dat\n")
    f.write("resdt.dat\n")
    f.write("resdd.dat\n")
    f.write("dthuv.mat\n")
    f.write("ddhuv.mat\n")
    f.write("pmf.mat\n")
    f.write("3, 3\n")
    f.write("\'meth\', \'cl1  \', \'cl2 \'\n")
    f.write("\'O   \', \'H   \', \'H   \'\n")
    f.write(str(r_ab))
    f.write(", ")
    f.write(str(r_bc))
    f.write(", ")
    f.write(str(r_ac))
    f.write("\n")
    f.write("0.9572, 0.9572, 1.5139\n")
    f.write("1.908, 2.2890, 2.2890\n")
    f.write("1.575, 0.5, 0.5\n")
    f.write("0.75006d-15, 0.12689d-14, 0.12689d-14\n")
    f.write("1.0605d-14, 3.8058d-15, 3.8058d-15\n")
    #f.write("0.24, -0.25, -0.9\n")
    # for i in range(3):
    #     print(zi[i])
    f.write(str(zi[1]))
    f.write(",")
    f.write(str(zi[0]))
    f.write(",")
    f.write(str(zi[2]))
    f.write("\n")
    f.write("-0.834, 0.417, 0.417\n")
    f.write("298.15, 0.03334\n")
    f.write("1.0\n")


def write_for_0():
    f.write("0, 1, 1, 1, 0\n")
    f.write("1,1\n")
    f.write("1, 0, 0,0,0,0,0,0,0\n")
    f.write("-90.0\n")
    f.write("1.d-6\n")
    f.write("0.0\n")
    f.write("257.21d-6\n")

def write_others(i):
    f.write("1, 1, 1, 1, 0\n")
    f.write("1,1\n")
    f.write("1, 0, 0,0,0,0,0,0,0\n")
    f.write("-90.0\n")
    f.write("1.d-6\n")
    f.write(str(i*0.1)+"\n")
    f.write("257.21d-6\n")

def write_last():
    f.write("\n")
    f.write("c***  read molecular parameters\n")
    f.write("      read(1,*) nu, nv\n")
    f.write("      read(1,*) (siteu(i),i=1,nu)\n")
    f.write("      read(1,*) (sitev(i),i=1,nv)\n")
    f.write("      if (nu .ne. 1) then\n")
    f.write("        read(1,*) ((lu(i,j),j=i+1,nu),i=1,nu-1)\n")
    f.write("      endif\n")
    f.write("      if (nv .ne. 1) then\n")
    f.write("        read(1,*) ((lv(i,j),j=i+1,nv),i=1,nv-1)\n")
    f.write("      endif\n")
    f.write("      read(1,*) (ruu(i),i=1,nu)\n")
    f.write("      read(1,*) (rvv(i),i=1,nv)\n")
    f.write("      read(1,*) (epsu(i),i=1,nu)\n")
    f.write("      read(1,*) (epsv(i),i=1,nv)\n")
    f.write("      read(1,*) (zu(i),i=1,nu)\n")
    f.write("      read(1,*) (zv(i),i=1,nv)\n")
    f.write("      read(1,*) temp, dens\n")
    f.write("      read(1,*) d\n")
    f.write("\n")
    f.write("c***  read controll parameters\n")
    f.write("      read(1,*) istart, iclosure, intrau, intrav, irism\n")
    f.write("      read(1,*) ((ilj(i,j),j=1,nv),i=1,nu)\n")
    f.write("      read(1,*) potcut\n")
    f.write("      read(1,*) crite\n")
    f.write("      read(1,*) scale\n")
    f.write("      read(1,*) alpha\n")

def get_rho(t):
    rho = (999.83952 + 16.945176*t - 7.9870401*10**(-3)*t**2 \
            - 46.170461*10**(-6)*t**3 + 105.56302*10**(-9)*t**4\
            - 280.54253*10**(-12)*t**5)/(1 + 16.87985*10**(-3)*t)
    print("rho is ",rho)
    rho = rho*10**(-4)*6.02/18.0
    print("tanikannsanngo rho", rho)

def makeInputFile(r_ab,r_bc):
    switch = 1 #Rossky or hirata  
    #global zi
    zi = np.zeros(3)

#r_ab = 2.3742
#r_bc = 2.3699
# r_ab = 2.0178
# r_bc = 2.6809

    #r_ab = 2.3742 
    #r_bc = 2.3699 
    r_ac = r_ab + r_bc
    print("reaction coordinate is ", r_ab - r_bc)
#温度T
    T = 25
    get_rho(T)
    cal_ele(r_ab, r_bc,zi)
    for i in [0,1,2,3,4,6,8,10]:
        filename = "tip3p_" + str(i) +".dat"
        #print(filename)
        global f
        f = open(filename,"w")
        # if i == 0:
        #     cal_ele(r_ab, r_bc)
        #     print("in i =0zi is","\n")
        #     for j in range(3):
        #         print(zi[i])
        if switch == 0:
            write1_ross(r_ab,r_bc,r_ac,zi)
        else:
            write1_hira(r_ab,r_bc,r_ac,zi)

        if i == 0:
            write_for_0()
        else:
            write_others(i)
        #write_last()
        f.close()

#makeInputFile(2.3742,2.3699)
