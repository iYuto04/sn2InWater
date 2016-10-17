import variables as var
import numpy as np
import makeInputFile as makeFile
import subprocess
import math
import copy
import matplotlib.pyplot as plt


N_str = var.N_str
h = 0.001
global theta_val

def get_norm(x1,y1,z1,x2,y2,z2):
    norm = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    norm = norm**0.50
    return norm

def get_potential(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    D = [[234.524674, 234.524674, 64.925971],
         [220.244820, 220.244820, 284.999867]]

    beta = [[0.929968, 0.929968, 0.432955],
            [4.822681, 4.822681, 1.016811]]

    r_0 = [[1.776382, 1.776382, 2.094857],
           [1.785014, 1.785014, 2.186060]]
    Q = np.zeros(3)
    J = np.zeros(3)

    r1 = get_norm(x1,y1,z1,x2,y2,z2)
    r2 = get_norm(x2,y2,z2,x3,y3,z3)
    r3 = get_norm(x1,y2,z3,x3,y3,z3)
    r_val = [r1,r2,r3]
    E = np.zeros((2,3))
    for i in range(2): #singlet,triplet
        for j in range(3): #r1,r2,r3
            E[i][j] = D[i][j]*(1.0 - (-1)**(i)* math.exp(-beta[i][j]*(r_val[j] - r_0[i][j])))**2.0 - D[i][j]

    for i in range(3):
        Q[i] = (E[0][i] + E[1][i])/2.0
        J[i] = (E[0][i] - E[1][i])/2.0

    THETA = Q[0] + Q[1] + Q[2] - (J[0]*J[0]+ J[1]*J[1] + J[2]*J[2] - J[0]*J[1] - J[1]*J[2] - J[2]*J[0])**0.50 + 234.524671
    return THETA

def saddle():
    naiseki = 0.0
    phi_s_old = np.zeros(3)
    tau = np.zeros(3)
    test_nabla_theta = np.zeros(3)
    l_s = 0.0
    dt = 0.001
    for i in range(N_str + 1):
        W_phi[i] = get_potential(phi[i][0], 0, 0, phi[i][1], 0, 0, phi[i][2], 0, 0)

        makeFile.makeInputFile(phi[i][1] - phi[i][0], phi[i][2] - phi[i][1])
        subprocess.call(["zsh", "start_program.sh"])
        fileFreeEnergy = open("freeEnergy.dat", "r")
        freeEnergy = float(fileFreeEnergy.read()) + var.freeEnergyZeropoint
        fileFreeEnergy.close()
        W_phi[i] += freeEnergy
        if i == 0:
            endW = W_phi[i]

    phi_s_W = max(W_phi)

    for i in range(N_str + 1):
        if W_phi[i] == phi_s_W:
            s = i + 1  # ポテンシャルが最大となるインデックス
    l = s - 1  # sのひとつ左の点
    r = s + 1  # sのひとつ右の点
    # print("s,r,l")
    # print(s, r, l)
    for i in range(3):
        l_s += (phi[r][i] - phi[l][i]) * (phi[r][i] - phi[l][i])
    l_s = l_s ** 0.50
    # print("l_s is ", l_s)
    for i in range(3):
        tau[i] = (phi[r][i] - phi[l][i]) / l_s
    # print("in saddle function")
    # print(tau)
    n_step = int(math.log(N_str ** 4))
    # print("n_step is")
    # print(n_step)
    # n_step = 12
    for i in range(3):
        phi_s[i] = phi[s][i]
    for step in range(n_step):
        # get_nabla_THETA(phi_s[0], phi_s[1], phi_s[2], phi_s[3], phi_s[4], phi_s[5], phi_s[6], phi_s[7], phi_s[8])

        for j in range(3): #Cl, CH3,CL'
            new_plus = copy.deepcopy(phi_s)
            new_minus = copy.deepcopy(phi_s)

            for switch in [-1,1]:
                if switch == -1:
                    new_minus[j] += switch * h
                    # print(new_minus)
                    theta_minus = get_potential(new_minus[0], 0, 0, \
                                                new_minus[1], 0, 0, \
                                                new_minus[2], 0, 0)

                    makeFile.makeInputFile(new_minus[1] - new_minus[0],new_minus[2] - new_minus[1])
                    subprocess.call(["zsh", "start_program.sh"])
                    fileFreeEnergy = open("freeEnergy.dat", "r")
                    freeEnergy = float(fileFreeEnergy.read()) + var.freeEnergyZeropoint
                    fileFreeEnergy.close()
                    theta_minus += freeEnergy



                elif switch == 1:
                    new_plus[j] += switch * h

                    theta_plus = get_potential(new_plus[0], 0, 0, \
                                               new_plus[1], 0, 0, \
                                               new_plus[2], 0, 0)

                    makeFile.makeInputFile(phi_s[1] - phi_s[0], phi_s[2] - phi_s[1])
                    subprocess.call(["zsh", "start_program.sh"])
                    fileFreeEnergy = open("freeEnergy.dat", "r")
                    freeEnergy = float(fileFreeEnergy.read()) + var.freeEnergyZeropoint
                    fileFreeEnergy.close()
                    theta_plus += freeEnergy

            test_nabla_theta[j] = (theta_plus - theta_minus) / 2.0 / h

        # #while True:
        for i in range(3):  # cl,ch3,cl'
            # for j in range(3):  # x,y,z
                #naiseki += nabla_theta[i][j] * tau[3 * i + j]
            naiseki += test_nabla_theta[i] * tau[i]
        for i in range(3):
            # for j in range(3):
                # phi_s[3 * i + j] = phi_s[3 * i + j] + dt * (-nabla_theta[i][j] + 2.0 * naiseki * tau[3 * i + j])
            phi_s[i] = phi_s[i] + dt * (-test_nabla_theta[i] + 2.0 * naiseki * tau[i])
        naiseki = 0.0
    theta_val = get_potential(phi_s[0], 0, 0, phi_s[1], 0, 0, phi_s[2], 0, 0)

    makeFile.makeInputFile(phi_s[1] - phi_s[0], phi_s[2] - phi_s[1])
    subprocess.call(["zsh", "start_program.sh"])
    fileFreeEnergy = open("freeEnergy.dat", "r")
    freeEnergy = float(fileFreeEnergy.read()) + var.freeEnergyZeropoint
    fileFreeEnergy.close()
    theta_val += freeEnergy

    print("the potential at saddle point")
    print(theta_val)
    print(phi_s)
    print("barrier height is")
    print(theta_val - endW)
    return theta_val

def plot_pote():
    plt_x = np.zeros(N_str+1)
    plt_y = np.zeros(N_str+1)
    plt_s = np.zeros(2)
    all_plot_x = np.zeros(N_str+2)
    all_plot_y = np.zeros(N_str+2)

    for i in range(N_str+1):
        plot_pote_r1 = get_norm(phi[i][0],0,0,phi[i][1],0,0)
        plot_pote_r2 = get_norm(phi[i][1],0,0,phi[i][2],0,0)
        plt_x[i] = plot_pote_r1 - plot_pote_r2
        plt_y[i] = W_phi[i]
        #plt_y[i] = get_potential(phi[i][0],0,0,phi[i][1],0,0,phi[i][2],0,0)

    plot_s_r1 = get_norm(phi_s[0],0,0,phi_s[1],0,0)
    plot_s_r2 = get_norm(phi_s[1],0,0,phi_s[2],0,0)
    plt_s[0] = plot_s_r1 - plot_s_r2
    plt_s[1] = theta_val

    for i in range(N_str+1):
        all_plot_x[i] = plt_x[i]
        all_plot_y[i] = plt_y[i]
    all_plot_x[N_str+1] = plt_s[0]
    all_plot_y[N_str+1] = plt_s[1]
    print("-----------all_plot_x------------")
    print(all_plot_x)


    for i in range(N_str+1):
        for j in range(i+1,N_str+2):
            if all_plot_x[i] > all_plot_x[j]:
                change_x = all_plot_x[i]
                change_y = all_plot_y[i]
                all_plot_x[i] = all_plot_x[j]
                all_plot_y[i] = all_plot_y[j]
                all_plot_x[j] = change_x
                all_plot_y[j] = change_y
    print("-----------after sort all_plot_x------------")
    print(all_plot_x)
    plt.plot(plt_x,plt_y,"o",markersize = 13,color = "green")
    plt.scatter(plt_s[0],plt_s[1],marker = "*",color = "black",s = 50)
    plt.plot(all_plot_x,all_plot_y,linewidth = 2,color = "m")
    plt.legend()
    plt.savefig("plot_pote.png")
    plt.show()

def readFile():
    f = open("logLast25.dat","r")
    while True:
        readLine = f.readline()
        if readLine == "":
            break
        else:
            a,b,c = map(float,readLine.split())
            phi.append([a,b,c])
phi = []
W_phi = np.zeros(N_str + 1)
phi_s = np.zeros(3)
readFile()
theta_val = saddle()
plot_pote()
