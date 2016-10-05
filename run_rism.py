#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import sys
import os

erg   = 1.0                # unit
cm    = 1.0                # unit
statC = 1.0                # unit
kcal  = 4.184e+10          * erg
angs  = 1.e-8              * cm
mole  = 6.0221412927e+23
qe    = 4.803e-10          * statC
e_unitc   = qe / (kcal * angs / mole)**0.5 # 18.2218395 (kcal A / mol)^1/2
kcal_mole = kcal/mole                      # 6.94769484248e-14 erg

#---------------------------------------------------------------
def water_model(name) :
    nv           = 3
    site_name_v  = ['O   ','H   ','H   ']
    iljv         = [1, 0, 0]
    if   name == 'SPC' :
        title        = 'LIQUID WATER [SPC MODEL]'
        lv           = [1.0,     1.0,     1.633]             # rOH, rOH, rHH
        rvv          = [1.583, 0.5, 0.5]                     # sigma[O,H,H]/2
        epsv         = [1.08e-14, 3.79e-15, 3.79e-15]        # epsv[nv]
        zv           = [-0.82, 0.41, 0.41]                   # zv[nv]
    elif name == 'SPC2' :
        title        = 'LIQUID WATER [SPC2 MODEL]'
        lv           = [1.0,     1.0,     1.633]             # rOH, rOH, rHH
        rvv          = [1.583, 0.2, 0.2]                     # sigma[O,H,H]/2
        epsv         = [1.08e-14, 3.79e-15, 3.79e-15]        # epsv[nv]
        zv           = [-0.82, 0.41, 0.41]                   # zv[nv]
    elif name == 'PR-SPC/E' :
        title        = 'LIQUID WATER [PR-SPC/E MODEL]'
        lv           = [1.0,     1.0,     1.633]             # rOH, rOH, rHH
        rvv          = [1.5829, 0.4/2, 0.4/2]                # sigma[O,H,H]/2
        epsv         = [1.079e-14, 3.1959e-15, 3.1959e-15]   # epsv[nv]
        zv           = [-0.8476, 0.4238, 0.4238]             # zv[nv]
    elif name == 'PR-TIP3P' :
        title        = 'LIQUID WATER [PR-TOP3P MODEL]'
        lv           = [0.9572, 0.9572,  1.5631]             # rOH, rOH, rHH
        rvv          = [1.57535, 0.4/2, 0.4/2]               # sigma[O,H,H]/2
        epsv         = [1.056e-14, 3.1959e-15, 3.1959e-15]   # epsv[nv]
        zv           = [-0.8340, 0.8340, 0.4238]             # zv[nv]
    elif name =="tip3p" :
        title        = 'LIQUID WATER [tip3p MODEL]'
        lv           = [0.9572, 0.9572,  1.5139]             # rOH, rOH, rHH
        rvv          = [1.575, 0.4/2, 0.4/2]               # sigma[O,H,H]/2
        epsv         = [1.056e-14, 3.1975e-15, 3.1975e-15]   # epsv[nv]
        zv           = [-0.8340, 0.4170, 0.4170]             # zv[nv]

    else :
        title        = None
        lv           = None
        rvv          = None
        epsv         = None
        zv           = None
    return title, nv, site_name_v, lv, rvv, epsv, zv, iljv

#---------------------------------------------------------------
def cl_atom() :
    nu           = 1
    site_name_u  = ['Cl  ']

    lu           = []                                        # lu[i,j] i<j<nu
    ruu          = [2.205]                                   # sigma[nu]/2
    epsu         = [0.8202298e-14]                           # epsu[nu]
    zu           = [-1.0]                                    # zu[nu]
    ilju         = [1]                                       # ilj[nu]

    return nu, site_name_u, lu, ruu, epsu, zu, ilju

def sn2_mol() :
    nu = 3
    site_name_u = ["meth", "cl1 ", "cl2 "]
    lu   = [2.0178, 2.6809, 4.6987]
    ruu  = [1.908, 2.2890, 2.2890]
    epsu = [0.75006e-15, 0.12689e-14, 0.12689e-14]
    zu   = [0.21117, -0.326232, -0.8859]
    ilju = [1]
    return nu, site_name_u, lu, ruu, epsu, zu, ilju

#---------------------------------------------------------------
def make_exuv_config_file(title, water_model_name, water_dir, counter, cfg_file):
    dscale       = 0.1
    c_str        = "%02d"%(counter,)

    outfile      = 'output.dat'               # output
    uv_data_file = water_dir+'uv.dat'         # input (from solvent calculstion)
    restart_file = 'restart.dat'              # restart file
    gr_file      = 'gr.mat'                   # radial distribution function 
    dtuv_file    = water_dir+'dtuv.dat'       # input (d/dT)
    dduv_file    = water_dir+'dduv.dat'       # input (d/d(density))
    resdt_file   = 'resdt.dat'                # restart (d/dT)
    resdd_file   = 'resdd.dat'                # restart (d/d(density))
    dthuv_file   = 'dthuv.mat'                # d/dT[pair correlation function]
    ddhuv_file   = 'ddhuv.mat'                # d/dd[pair correlation function]
    pmf_file     = 'pmf.mat'           # interatomic potential of mean force

    water_title, nv, site_name_v, lv, rvv, epsv, zv, iljv \
        = water_model(water_model_name)
    if water_title == None : 
        print 'Water model name invalid'
        sys.exit(-1) 
    #nu, site_name_u, lu, ruu, epsu, zu, ilju = cl_atom()
    nu, site_name_u, lu, ruu, epsu, zu, ilju = sn2_mol()

    #temp, dense  = 298.15, 0.03336			     # temp, dense
    temp, dense  = 273.15, 0.03345			     # temp, dense
    d            = 1.0					     # d

    if counter == 0 : 
        istart, iclosure, intrau, intrav, irism = 0, 1, 0, 1, 0
    else :
        istart, iclosure, intrau, intrav, irism = 1, 1, 0, 1, 0
    idt, idd     = 1,1                                       # idt, idd
    ilj          = [ilju[i]*iljv[j] for j in range(nv) \
                        for i in range(nu)]                  # ilj[nu,nv]
    potcut       = -90.0                                     # potcut
    crite        = 1.e-6                                     # crite
    scale        = dscale*counter                            # scale
    alpha        = 257.21e-6                                 # alpha

    #---------------------------------------------------------------
    f = open(cfg_file, "w")
    # title
    f.write(str(title)+"\n")
    # filles
    for str1 in [outfile, uv_data_file, restart_file, gr_file,\
                     dtuv_file, dduv_file, resdt_file, resdd_file,\
                     dthuv_file, ddhuv_file, pmf_file] : \
                     f.write(str(str1)+"\n")
    #
    f.write(str((nu, nv))[1:-1]+"\n")
    f.write(str(site_name_u)[1:-1]+"\n")
    f.write(str(site_name_v)[1:-1]+"\n")
    if nu > 1 : f.write(str(lu)[1:-1]+"\n")
    if nv > 1 : f.write(str(lv)[1:-1]+"\n")
    for list1 in [ruu, rvv, epsu, epsv, zu, zv] : \
            f.write(str(list1)[1:-1]+"\n")
    f.write(str((temp, dense))[1:-1]+"\n")
    f.write(str(d)+"\n")
    #
    f.write(str((istart, iclosure, intrau, intrav, irism))[1:-1]+"\n")
    f.write(str((idt, idd))[1:-1]+"\n")
    f.write(str(ilj)[1:-1]+"\n")
    for val in [potcut, crite, scale, alpha] : f.write(str(val)+"\n")
    f.closed
    return
                
def read_config_file(cfg_file):
    f = open(cfg_file, "r")
    for str1 in f.readlines() : print str1[:-1]
    f.closed
    return

def main() :
    title        = 'sn2 IN tip3p WATER'
    water_model_name   = 'tip3p'
    water_dir    = '../../exvv/sn2_tip3p/'
    program = '../../src/exuv/exuv'
    cfg_file     = "tip3p.dat"

    print 'calculate "%s" using %s\n'%(title, program)
    for counter in [0,1,2,3,4,6,8,10] :
        make_exuv_config_file(title, water_model_name, water_dir, counter, cfg_file)
        #os.system(program + " %s | tee out%02d.dat" % (cfg_file, counter)) 
        os.system(program + " %s > out%02d.dat" % (cfg_file, counter)) 

#---------------------------------------------------------------
if __name__ == "__main__":
    main()
#---------------------------------------------------------------
