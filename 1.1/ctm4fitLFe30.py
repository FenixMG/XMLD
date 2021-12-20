# Debugged vertion finished with foor loops
import numpy as np
from scipy.io import loadmat
from math import factorial
from AtomicLib_LEdgeFe import AtomicLib_LEdgeFe
import Sticks2Band as S2B
from scipy.interpolate import PchipInterpolator as pchip
from scipy.linalg import eigh
import os
from scipy.sparse import csc_matrix

def ctm4fitLFe30(x, Nelec, Red_F2, Red_F2pd,tenDq, Ds, Dt, Bz, Temp, Sym, MinL, MaxL, Gauss, SplitE , ShiftE, Container, XASInt, DichrInt, XMLD_50, XMLD_90):
    Cont = np.array([0, 0, 0], dtype=object)
    mm = 0

    current_direction = os.path.dirname(os.path.realpath(__file__))

    Strm = np.array([Nelec, Red_F2, Red_F2pd, tenDq, Ds, Dt, Bz, Temp, Sym])

    try:
        str4cont = current_direction + "/temp_sim_Fe.npy"
        Cont = np.load(str4cont, allow_pickle=True)
        pass
    except:
        pass

    try:
        Id = Cont[Container - 1]['Id']

        if Strm[0] == Id[0] and Strm[1] == Id[1] and Strm[2] == Id[2] and Strm[3] == Id[3] and Strm[4] == Id[4] and \
                Strm[5] == Id[5] and Strm[6] == Id[6] and Strm[7] == Id[7] and Strm[8] == Id[8]:
            mm = 1
        pass
    except:
        pass

    if mm == 0:
        Id = Strm
        AP = AtomicLib_LEdgeFe(Nelec)

        Eshift = AP['Shift']

        OPLE = loadmat(current_direction + "/OperadoresFesp.mat")

        #Fe3 - 2p6 3d5
        OppF0dd_d5 = OPLE['OppF0dd_d5']
        OppF2dd_d5 = OPLE['OppF2dd_d5']
        OppF4dd_d5 = OPLE['OppF4dd_d5']
        OppLS_d_d5 = OPLE['OppLS_d_d5']

        OppSx_d5 = OPLE['OppSx_d5']
        OppSy_d5 = OPLE['OppSy_d5']
        OppSz_d5 = OPLE['OppSz_d5']

        OppC40_d5 = OPLE['OppC40_d5']
        OppC4_1_d5 = OPLE['OppC4_1_d5']
        OppC4_2_d5 = OPLE['OppC4_2_d5']
        OppC4_3_d5 = OPLE['OppC4_3_d5']
        OppC4_4_d5 = OPLE['OppC4_4_d5']
        OppC41_d5 = OPLE['OppC41_d5']
        OppC42_d5 = OPLE['OppC42_d5']
        OppC43_d5 = OPLE['OppC43_d5']
        OppC44_d5 = OPLE['OppC44_d5']
        OppC20_d5 = OPLE['OppC20_d5']
        OppC2_1_d5 = OPLE['OppC2_1_d5']
        OppC2_2_d5 = OPLE['OppC2_2_d5']
        OppC21_d5 = OPLE['OppC21_d5']
        OppC22_d5 = OPLE['OppC22_d5']

        # Fe3 - 2p5 3d6
        OppF0dd_p5d6 = OPLE['OppF0dd_p5d6']
        OppF2dd_p5d6 = OPLE['OppF2dd_p5d6']
        OppF4dd_p5d6 = OPLE['OppF4dd_p5d6']

        OppF0pd_p5d6 = OPLE['OppF0pd_p5d6']
        OppF2pd_p5d6 = OPLE['OppF2pd_p5d6']
        OppG1pd_p5d6 = OPLE['OppG1pd_p5d6']
        OppG3pd_p5d6 = OPLE['OppG3pd_p5d6']

        OppLS_d_p5d6 = OPLE['OppLS_d_p5d6']
        OppLS_p_p5d6 = OPLE['OppLS_p_p5d6']

        OppC40_p5d6 = OPLE['OppC40_p5d6']
        OppC4_1_p5d6 = OPLE['OppC4_1_p5d6']
        OppC4_2_p5d6 = OPLE['OppC4_2_p5d6']
        OppC4_3_p5d6 = OPLE['OppC4_3_p5d6']
        OppC4_4_p5d6 = OPLE['OppC4_4_p5d6']
        OppC41_p5d6 = OPLE['OppC41_p5d6']
        OppC42_p5d6 = OPLE['OppC42_p5d6']
        OppC43_p5d6 = OPLE['OppC43_p5d6']
        OppC44_p5d6 = OPLE['OppC44_p5d6']
        OppC20_p5d6 = OPLE['OppC20_p5d6']
        OppC2_1_p5d6 = OPLE['OppC2_1_p5d6']
        OppC2_2_p5d6 = OPLE['OppC2_2_p5d6']
        OppC21_p5d6 = OPLE['OppC21_p5d6']
        OppC22_p5d6 = OPLE['OppC22_p5d6']

        OppSx_p5d6 = OPLE['OppSx_p5d6']
        OppSy_p5d6 = OPLE['OppSy_p5d6']
        OppSz_p5d6 = OPLE['OppSz_p5d6']

        # Fe3 Transitions 2p6 3d5 - 2p5 3d6(dipoloelectrico)
        OppT1_1_d5 = OPLE['OppT1_1_d5'] # componente C1_1
        OppT11_d5 = OPLE['OppT11_d5'] # componente C11

        # Fe2 - 2p6 3d6
        OppF0dd_d6 = OPLE['OppF0dd_d6']
        OppF2dd_d6 = OPLE['OppF2dd_d6']
        OppF4dd_d6 = OPLE['OppF4dd_d6']
        OppLS_d_d6 = OPLE['OppLS_d_d6']

        OppSx_d6 = OPLE['OppSx_d6']
        OppSy_d6 = OPLE['OppSy_d6']
        OppSz_d6 = OPLE['OppSz_d6']

        OppC40_d6 = OPLE['OppC40_d6']
        OppC4_1_d6 = OPLE['OppC4_1_d6']
        OppC4_2_d6 = OPLE['OppC4_2_d6']
        OppC4_3_d6 = OPLE['OppC4_3_d6']
        OppC4_4_d6 = OPLE['OppC4_4_d6']
        OppC41_d6 = OPLE['OppC41_d6']
        OppC42_d6 = OPLE['OppC42_d6']
        OppC43_d6 = OPLE['OppC43_d6']
        OppC44_d6 = OPLE['OppC44_d6']
        OppC20_d6 = OPLE['OppC20_d6']
        OppC2_1_d6 = OPLE['OppC2_1_d6']
        OppC2_2_d6 = OPLE['OppC2_2_d6']
        OppC21_d6 = OPLE['OppC21_d6']
        OppC22_d6 = OPLE['OppC22_d6']

        # Fe2 - 2p5 3d7

        OppF0dd_p5d7 = OPLE['OppF0dd_p5d7']
        OppF2dd_p5d7 = OPLE['OppF2dd_p5d7']
        OppF4dd_p5d7 = OPLE['OppF4dd_p5d7']

        OppF0pd_p5d7 = OPLE['OppF0pd_p5d7']
        OppF2pd_p5d7 = OPLE['OppF2pd_p5d7']

        OppG1pd_p5d7 = OPLE['OppG1pd_p5d7']
        OppG3pd_p5d7 = OPLE['OppG3pd_p5d7']

        OppLS_d_p5d7 = OPLE['OppLS_d_p5d7']
        OppLS_p_p5d7 = OPLE['OppLS_p_p5d7']

        OppSx_p5d7 = OPLE['OppSx_p5d7']
        OppSy_p5d7 = OPLE['OppSy_p5d7']
        OppSz_p5d7 = OPLE['OppSz_p5d7']

        OppC40_p5d7 = OPLE['OppC40_p5d7']
        OppC4_1_p5d7 = OPLE['OppC4_1_p5d7']
        OppC4_2_p5d7 = OPLE['OppC4_2_p5d7']
        OppC4_3_p5d7 = OPLE['OppC4_3_p5d7']
        OppC4_4_p5d7 = OPLE['OppC4_4_p5d7']
        OppC41_p5d7 = OPLE['OppC41_p5d7']
        OppC42_p5d7 = OPLE['OppC42_p5d7']
        OppC43_p5d7 = OPLE['OppC43_p5d7']
        OppC44_p5d7 = OPLE['OppC44_p5d7']
        OppC20_p5d7 = OPLE['OppC20_p5d7']
        OppC2_1_p5d7 = OPLE['OppC2_1_p5d7']
        OppC2_2_p5d7 = OPLE['OppC2_2_p5d7']
        OppC21_p5d7 = OPLE['OppC21_p5d7']
        OppC22_p5d7 = OPLE['OppC22_p5d7']

        # Fe2 Transitions 2p6 3d6 - 2p5 3d7(dipolo electrico)
        OppT1_1_d6 = OPLE['OppT1_1_d6'] # componente C1_1
        OppT11_d6 = OPLE['OppT11_d6'] # componente C11


        # Final state
        # H^ee
        Udd = 0
        F2 = AP['Ex']['Fdd2'] * Red_F2 / 100 * 0.8
        F4 = AP['Ex']['Fdd4'] * Red_F2 / 100 * 0.8
        F0 = Udd + 2 / 63 * F2 + 2 / 63 * F4

        F2pd = AP['Ex']['Fpd'] * Red_F2pd / 100 * 0.8
        G1pd = AP['Ex']['Gpd1'] * Red_F2pd / 100 * 0.8
        G3pd = AP['Ex']['Gpd3'] * Red_F2pd / 100 * 0.8
        F0pd = 1 / 15 * G1pd + 3 / 70 * G3pd

        # H^SOC
        Zeta_d = AP['Ex']['SOC3d']
        Zeta_p = AP['Ex']['SOC2p']

        #Campo magnético

        V0 = np.array([0.9962, 0, 0.0867]) * Bz
        V5 = np.array([0.9924, 0.0868, 0.0867]) * Bz
        V50 = np.array([0.6403, 0.7631, 0.0867]) * Bz
        V90 = np.array([0, 0.9962, 0.0867]) * Bz

        #Crate Common Hamiltonian
        if Nelec == 5:
            Hamiltonian = OppF0dd_p5d6 * F0 + OppF2dd_p5d6 * F2 + OppF4dd_p5d6 * F4
            Hamiltonian = Hamiltonian + OppLS_d_p5d6 * Zeta_d + OppLS_p_p5d6 * Zeta_p
            Hamiltonian = Hamiltonian + OppG1pd_p5d6 * G1pd + OppG3pd_p5d6 * G3pd + OppF2pd_p5d6 * F2pd + OppF0pd_p5d6 * F0pd
        else:
            Hamiltonian = OppF0dd_p5d7 * F0 + OppF2dd_p5d7 * F2 + OppF4dd_p5d7 * F4
            Hamiltonian = Hamiltonian + OppLS_d_p5d7 * Zeta_d + OppLS_p_p5d7 * Zeta_p
            Hamiltonian = Hamiltonian + OppG1pd_p5d7 * G1pd + OppG3pd_p5d7 * G3pd + OppF2pd_p5d7 * F2pd + OppF0pd_p5d7 * F0pd

        if Sym == 1:  # D4h/ Oh - z del sitio local coincide con z del cristal

            if Nelec == 5:
                #Crystal Field D4h
                OpptenDq_p5d6 = 21 / 10 * (OppC40_p5d6 + np.sqrt(5 / 14) * (OppC4_4_p5d6 + OppC44_p5d6))
                OppDs_p5d6 = -7 * OppC20_p5d6
                OppDt_p5d6 = -21 * OppC40_p5d6

                Hamiltonian = Hamiltonian + OpptenDq_p5d6 * tenDq + OppDs_p5d6 * Ds + OppDt_p5d6 * Dt

                Hamiltonian0 = Hamiltonian + OppSz_p5d6 * V0[2] + 1j * OppSy_p5d6 * V0[1] + OppSx_p5d6 * V0[0]
                Hamiltonian5 = Hamiltonian + OppSz_p5d6*V5[2] + 1j*OppSy_p5d6*V5[1] + OppSx_p5d6*V5[0]
                Hamiltonian50 = Hamiltonian + OppSz_p5d6 * V50[2] + 1j * OppSy_p5d6 * V50[1] + OppSx_p5d6 * V50[0]
                Hamiltonian90 = Hamiltonian + OppSz_p5d6 * V90[2] + 1j * OppSy_p5d6 * V90[1] + OppSx_p5d6 * V90[0]
            else:
                #Cristal Field D4h
                OpptenDq_p5d7 = 21 / 10 * (OppC40_p5d7 + np.sqrt(5 / 14) * (OppC4_4_p5d7 + OppC44_p5d7))
                OppDs_p5d7 = -7 * OppC20_p5d7
                OppDt_p5d7 = -21 * OppC40_p5d7

                Hamiltonian = Hamiltonian + OpptenDq_p5d7 * tenDq + OppDs_p5d7 * Ds + OppDt_p5d7 * Dt

                Hamiltonian0 = Hamiltonian + OppSz_p5d7 * V0[2] + 1j * OppSy_p5d7 * V0[1] + OppSx_p5d7 * V0[0]
                Hamiltonian5 = Hamiltonian + OppSz_p5d7*V5[2] + 1j*OppSy_p5d7*V5[1] + OppSx_p5d7*V5[0]
                Hamiltonian50 = Hamiltonian + OppSz_p5d7 * V50[2] + 1j * OppSy_p5d7 * V50[1] + OppSx_p5d7 * V50[0]
                Hamiltonian90 = Hamiltonian + OppSz_p5d7 * V90[2] + 1j * OppSy_p5d7 * V90[1] + OppSx_p5d7 * V90[0]

            # D4h
            #Sparse to full Matrix
            Hamiltonian0 = Hamiltonian0.todense()
            Hamiltonian5 = Hamiltonian5.todense()
            Hamiltonian50 = Hamiltonian50.todense()
            Hamiltonian90 = Hamiltonian90.todense()
            #diagonalization is faster with Full matrix (eigs for sparse uses different approach).
            Ef0, CCf0 = eigh(Hamiltonian0)
            Ef5,CCf5 = eigh(Hamiltonian5)
            Ef50, CCf50 = eigh(Hamiltonian50)
            Ef90, CCf90 = eigh(Hamiltonian90)

            Eff0 = Ef0 + Eshift
            Eff5 = Ef5 + Eshift
            Eff50 = Ef50 + Eshift
            Eff90 = Ef90 + Eshift

        else: # D3d z of local site is different than in crystal
            #Fe2+ only (Fe3+ does not distort)

            #site 1
            OpptenDq_p5d7_site1 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_p5d7 + 3 / 10 * np.sqrt(35 / 2) * OppC44_p5d7
            OpptenDq_p5d7_site1 = OpptenDq_p5d7_site1 + 21 / 10 * OppC40_p5d7

            OppDs_p5d7_site1 = 7j / np.sqrt(6) * OppC2_2_p5d7 - 7j / np.sqrt(6) * OppC22_p5d7
            OppDs_p5d7_site1 = OppDs_p5d7_site1 + (7 - 7j) / np.sqrt(6) * OppC2_1_p5d7 - (7 + 7j) / np.sqrt(6) * OppC21_p5d7

            OppDt_p5d7_site1 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_p5d7 + (-7 / 6 - 7j / 6) * np.sqrt(35) * OppC4_3_p5d7
            OppDt_p5d7_site1 = OppDt_p5d7_site1 + 7j / 3 * np.sqrt(10) * OppC4_2_p5d7 + (-7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_p5d7
            OppDt_p5d7_site1 = OppDt_p5d7_site1 + 49 / 6 * OppC40_p5d7
            OppDt_p5d7_site1 = OppDt_p5d7_site1 - 7j / 3 * np.sqrt(10) * OppC42_p5d7 + (7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_p5d7
            OppDt_p5d7_site1 = OppDt_p5d7_site1 + 7 / 6 * np.sqrt(35 / 2) * OppC44_p5d7 + (7 / 6 - 7j / 6) * np.sqrt(35) * OppC43_p5d7

            Hamiltonian_site1 = Hamiltonian + OpptenDq_p5d7_site1 * tenDq + OppDs_p5d7_site1 * Ds + OppDt_p5d7_site1 * Dt

            Hamiltonian01 = Hamiltonian_site1 + OppSx_p5d7 * V0[0] + 1j * OppSy_p5d7 * V0[1] + OppSz_p5d7 * V0[2]
            Hamiltonian51 = Hamiltonian_site1 + OppSx_p5d7 * V5[0] + 1j * OppSy_p5d7 * V5[1] + OppSz_p5d7 * V5[2]
            Hamiltonian501 = Hamiltonian_site1 + OppSx_p5d7 * V50[0] + 1j * OppSy_p5d7 * V50[1] + OppSz_p5d7 * V50[2]
            Hamiltonian901 = Hamiltonian_site1 + OppSx_p5d7 * V90[0] + 1j * OppSy_p5d7 * V90[1] + OppSz_p5d7 * V90[2]

            # site 2

            OpptenDq_p5d7_site2 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_p5d7 + 3 / 10 * np.sqrt(35 / 2) * OppC44_p5d7
            OpptenDq_p5d7_site2 = OpptenDq_p5d7_site2 + 21 / 10 * OppC40_p5d7

            OppDs_p5d7_site2 = -7j / np.sqrt(6) * OppC2_2_p5d7 + 7j / np.sqrt(6) * OppC22_p5d7
            OppDs_p5d7_site2 = OppDs_p5d7_site2 - (7 + 7j) / np.sqrt(6) * OppC2_1_p5d7 + (7 - 7j) / np.sqrt(6) * OppC21_p5d7

            OppDt_p5d7_site2 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_p5d7 + (7 / 6 - 7j / 6) * np.sqrt(35) * OppC4_3_p5d7
            OppDt_p5d7_site2 = OppDt_p5d7_site2 - 7j / 3 * np.sqrt(10) * OppC4_2_p5d7 + (7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_p5d7
            OppDt_p5d7_site2 = OppDt_p5d7_site2 + 49 / 6 * OppC40_p5d7
            OppDt_p5d7_site2 = OppDt_p5d7_site2 + 7j / 3 * np.sqrt(10) * OppC42_p5d7 + (-7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_p5d7
            OppDt_p5d7_site2 = OppDt_p5d7_site2 + 7 / 6 * np.sqrt(35 / 2) * OppC44_p5d7 - (7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_p5d7

            Hamiltonian_site2 = Hamiltonian + OpptenDq_p5d7_site2 * tenDq + OppDs_p5d7_site2 * Ds + OppDt_p5d7_site2 * Dt

            Hamiltonian02 = Hamiltonian_site2 + OppSx_p5d7 * V0[0] + 1j * OppSy_p5d7 * V0[1] + OppSz_p5d7 * V0[2]
            Hamiltonian52 = Hamiltonian_site2 + OppSx_p5d7*V5[0] + 1j*OppSy_p5d7*V5[1] + OppSz_p5d7*V5[2]
            Hamiltonian502 = Hamiltonian_site2 + OppSx_p5d7 * V50[0] + 1j * OppSy_p5d7 * V50[1] + OppSz_p5d7 * V50[2]
            Hamiltonian902 = Hamiltonian_site2 + OppSx_p5d7 * V90[0] + 1j * OppSy_p5d7 * V90[1] + OppSz_p5d7 * V90[2]

            # site 3

            OpptenDq_p5d7_site3 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_p5d7 + 3 / 10 * np.sqrt(35 / 2) * OppC44_p5d7
            OpptenDq_p5d7_site3 = OpptenDq_p5d7_site3 + 21 / 10 * OppC40_p5d7

            OppDs_p5d7_site3 = 7j / np.sqrt(6) * OppC2_2_p5d7 - 7j / np.sqrt(6) * OppC22_p5d7
            OppDs_p5d7_site3 = OppDs_p5d7_site3 - (7 - 7j) / np.sqrt(6) * OppC2_1_p5d7 + (7 + 7j) / np.sqrt(6) * OppC21_p5d7

            OppDt_p5d7_site3 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_p5d7 + (7 / 6 + 7j / 6) * np.sqrt(35) * OppC4_3_p5d7
            OppDt_p5d7_site3 = OppDt_p5d7_site3 + 7j / 3 * np.sqrt(10) * OppC4_2_p5d7 + (7 / 6 - 7j / 6) * np.sqrt(5) * OppC4_1_p5d7
            OppDt_p5d7_site3 = OppDt_p5d7_site3 + 49 / 6 * OppC40_p5d7
            OppDt_p5d7_site3 = OppDt_p5d7_site3 - 7j / 3 * np.sqrt(10) * OppC42_p5d7 - (7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_p5d7
            OppDt_p5d7_site3 = OppDt_p5d7_site3 + 7 / 6 * np.sqrt(35 / 2) * OppC44_p5d7 + (-7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_p5d7

            Hamiltonian_site3 = Hamiltonian + OpptenDq_p5d7_site3 * tenDq + OppDs_p5d7_site3 * Ds + OppDt_p5d7_site3 * Dt

            Hamiltonian03 = Hamiltonian_site3 + OppSx_p5d7 * V0[0] + 1j * OppSy_p5d7 * V0[1] + OppSz_p5d7 * V0[2]
            Hamiltonian53 = Hamiltonian_site3 + OppSx_p5d7*V5[0] + 1j*OppSy_p5d7*V5[1] + OppSz_p5d7*V5[2]
            Hamiltonian503 = Hamiltonian_site3 + OppSx_p5d7 * V50[0] + 1j * OppSy_p5d7 * V50[1] + OppSz_p5d7 * V50[2]
            Hamiltonian903 = Hamiltonian_site3 + OppSx_p5d7 * V90[0] + 1j * OppSy_p5d7 * V90[1] + OppSz_p5d7 * V90[2]

            # site 4

            OpptenDq_p5d7_site4 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_p5d7 + 3 / 10 * np.sqrt(35 / 2) * OppC44_p5d7
            OpptenDq_p5d7_site4 = OpptenDq_p5d7_site4 + 21 / 10 * OppC40_p5d7

            OppDs_p5d7_site4 = -7j / np.sqrt(6) * OppC2_2_p5d7 + 7j / np.sqrt(6) * OppC22_p5d7
            OppDs_p5d7_site4 = OppDs_p5d7_site4 + (7 + 7j) / np.sqrt(6) * OppC2_1_p5d7 - (7 - 7j) / np.sqrt(6) * OppC21_p5d7

            OppDt_p5d7_site4 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_p5d7 + (-7 / 6 + 7j / 6) * np.sqrt(35) * OppC4_3_p5d7
            OppDt_p5d7_site4 = OppDt_p5d7_site4 - 7j / 3 * np.sqrt(10) * OppC4_2_p5d7 - (7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_p5d7
            OppDt_p5d7_site4 = OppDt_p5d7_site4 + 49 / 6 * OppC40_p5d7
            OppDt_p5d7_site4 = OppDt_p5d7_site4 + 7j / 3 * np.sqrt(10) * OppC42_p5d7 + (7 / 6 - 7j / 6) * np.sqrt(5) * OppC41_p5d7
            OppDt_p5d7_site4 = OppDt_p5d7_site4 + 7 / 6 * np.sqrt(35 / 2) * OppC44_p5d7 + (7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_p5d7

            Hamiltonian_site4 = Hamiltonian + OpptenDq_p5d7_site4 * tenDq + OppDs_p5d7_site4 * Ds + OppDt_p5d7_site4 * Dt

            Hamiltonian04 = Hamiltonian_site4 + OppSx_p5d7 * V0[0] + 1j * OppSy_p5d7 * V0[1] + OppSz_p5d7 * V0[2]
            Hamiltonian54 = Hamiltonian_site4 + OppSx_p5d7*V5[0] + 1j*OppSy_p5d7*V5[1] + OppSz_p5d7*V5[2]
            Hamiltonian504 = Hamiltonian_site4 + OppSx_p5d7 * V50[0] + 1j * OppSy_p5d7 * V50[1] + OppSz_p5d7 * V50[2]
            Hamiltonian904 = Hamiltonian_site4 + OppSx_p5d7 * V90[0] + 1j * OppSy_p5d7 * V90[1] + OppSz_p5d7 * V90[2]

            #sparse to full
            Hamiltonian01=Hamiltonian01.todense()
            Hamiltonian51=Hamiltonian51.todense()
            Hamiltonian501=Hamiltonian501.todense()
            Hamiltonian901=Hamiltonian901.todense()

            Hamiltonian02=Hamiltonian02.todense()
            Hamiltonian52 = Hamiltonian52.todense()
            Hamiltonian502 = Hamiltonian502.todense()
            Hamiltonian902=Hamiltonian902.todense()

            Hamiltonian03=Hamiltonian03.todense()
            Hamiltonian53 = Hamiltonian53.todense()
            Hamiltonian503 = Hamiltonian503.todense()
            Hamiltonian903=Hamiltonian903.todense()

            Hamiltonian04=Hamiltonian04.todense()
            Hamiltonian54 = Hamiltonian54.todense()
            Hamiltonian504 = Hamiltonian504.todense()
            Hamiltonian904=Hamiltonian904.todense()

            Ef01, CCf01 = eigh(Hamiltonian01)
            Ef51, CCf51 = eigh(Hamiltonian51)
            Ef501, CCf501 = eigh(Hamiltonian501)
            Ef901, CCf901 = eigh(Hamiltonian901)

            Ef02, CCf02 = eigh(Hamiltonian02)
            Ef52, CCf52 = eigh(Hamiltonian52)
            Ef502, CCf502 = eigh(Hamiltonian502)
            Ef902, CCf902 = eigh(Hamiltonian902)

            Ef03, CCf03 = eigh(Hamiltonian03)
            Ef53, CCf53 = eigh(Hamiltonian53)
            Ef503, CCf503 = eigh(Hamiltonian503)
            Ef903, CCf903 = eigh(Hamiltonian903)

            Ef04, CCf04 = eigh(Hamiltonian04)
            Ef54, CCf54 = eigh(Hamiltonian54)
            Ef504, CCf504 = eigh(Hamiltonian504)
            Ef904, CCf904 = eigh(Hamiltonian904)

            Eff01 = Ef01 + Eshift
            Eff51 = Ef51 + Eshift
            Eff501 = Ef501 + Eshift
            Eff901 = Ef901 + Eshift

            Eff02 = Ef02 + Eshift
            Eff52 = Ef52 + Eshift
            Eff502 = Ef502 + Eshift
            Eff902 = Ef902 + Eshift

            Eff03 = Ef03 + Eshift
            Eff53 = Ef53 + Eshift
            Eff503 = Ef503 + Eshift
            Eff903 = Ef903 + Eshift

            Eff04 = Ef04 + Eshift
            Eff54 = Ef54 + Eshift
            Eff504 = Ef504 + Eshift
            Eff904 = Ef904 + Eshift

        # Estado inicial

        # H^ee
        Udd = 0
        F2 = AP['Gd']['Fdd2'] * Red_F2 / 100 * 0.8
        F4 = AP['Gd']['Fdd4'] * Red_F2 / 100 * 0.8
        F0 = Udd + 2 / 63 * F2 + 2 / 63 * F4

        # H^SOC
        Zeta_d = AP['Gd']['SOC3d']

        #Create common Hamiltonian
        if Nelec==5:
            Hamiltonian = OppF0dd_d5 * F0 + OppF2dd_d5 * F2 + OppF4dd_d5 * F4
            Hamiltonian = Hamiltonian + OppLS_d_d5 * Zeta_d
        else:
            Hamiltonian = OppF0dd_d6 * F0 + OppF2dd_d6 * F2 + OppF4dd_d6 * F4
            Hamiltonian = Hamiltonian + OppLS_d_d6 * Zeta_d

        # Crando el Hamiltoneano
        if Sym == 1:  # D4h

            if Nelec == 5:
                OpptenDq = 21 / 10 * (OppC40_d5 + np.sqrt(5 / 14) * (OppC4_4_d5 + OppC44_d5))
                OppDs = -7 * OppC20_d5
                OppDt = -21 * OppC40_d5

                Hamiltonian = Hamiltonian + OpptenDq * tenDq + OppDs * Ds + OppDt * Dt
                Hamiltonian0 = Hamiltonian + OppSz_d5 * V0[2] + 1j * OppSy_d5 * V0[1] + OppSx_d5 * V0[0]
                Hamiltonian5 = Hamiltonian + OppSz_d5*V5[2] + 1j*OppSy_d5*V5[1] + OppSx_d5*V5[0]
                Hamiltonian50 = Hamiltonian + OppSz_d5 * V50[2] + 1j * OppSy_d5 * V50[1] + OppSx_d5 * V50[0]
                Hamiltonian90 = Hamiltonian + OppSz_d5 * V90[2] + 1j * OppSy_d5 * V90[1] + OppSx_d5 * V90[0]
            else:
                OpptenDq = 21 / 10 * (OppC40_d6 + np.sqrt(5 / 14) * (OppC4_4_d6 + OppC44_d6))
                OppDs = -7 * OppC20_d6
                OppDt = -21 * OppC40_d6

                Hamiltonian = Hamiltonian + OpptenDq * tenDq + OppDs * Ds + OppDt * Dt

                Hamiltonian0 = Hamiltonian + OppSz_d6 * V0[2] + 1j * OppSy_d6 * V0[1] + OppSx_d6 * V0[0]
                Hamiltonian5 = Hamiltonian + OppSz_d6*V5[2] + 1j*OppSy_d6*V5[1] + OppSx_d6*V5[0]
                Hamiltonian50 = Hamiltonian + OppSz_d6 * V50[2] + 1j * OppSy_d6 * V50[1] + OppSx_d6 * V50[0]
                Hamiltonian90 = Hamiltonian + OppSz_d6 * V90[2] + 1j * OppSy_d6 * V90[1] + OppSx_d6 * V90[0]

            # Sparse to full Matrix
            Hamiltonian0 = Hamiltonian0.todense()
            Hamiltonian5 = Hamiltonian5.todense()
            Hamiltonian50 = Hamiltonian50.todense()
            Hamiltonian90 = Hamiltonian90.todense()

            Ei0, CCi0 = eigh(Hamiltonian0)
            Ei5, CCi5 = eigh(Hamiltonian5)
            Ei50, CCi50 = eigh(Hamiltonian50)
            Ei90, CCi90 = eigh(Hamiltonian90)

            fB0 = np.zeros(Ei0.size, dtype='complex_')
            fB5 = np.zeros(Ei5.size, dtype='complex_')
            fB50 = np.zeros(Ei50.size, dtype='complex_')
            fB90 = np.zeros(Ei90.size, dtype='complex_')

            l = 0  # Variable auxialiar para recorrer fB y Ei

            for i in Ei0:
                fB0[l] = np.exp(-(i - Ei0[0]) / (8.61733035e-5 * Temp))
                fB5[l] = np.exp(-(Ei5[l] - Ei5[0]) / (8.61733035e-5 * Temp))
                fB50[l] = np.exp(-(Ei50[l] - Ei50[0]) / (8.61733035e-5 * Temp))
                fB90[l] = np.exp(-(Ei90[l] - Ei90[0]) / (8.61733035e-5 * Temp))
                l += 1

        else: #D3d

            # Crystal Field - Site 1

            OpptenDq_site1 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_d6 + 3 / 10 * np.sqrt(35 / 2) * OppC44_d6
            OpptenDq_site1 = OpptenDq_site1 + 21 / 10 * OppC40_d6

            OppDs_site1 = 7j / np.sqrt(6) * OppC2_2_d6 - 7j / np.sqrt(6) * OppC22_d6
            OppDs_site1 = OppDs_site1 + (7 - 7j) / np.sqrt(6) * OppC2_1_d6 - (7 + 7j) / np.sqrt(6) * OppC21_d6

            OppDt_site1 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_d6 + (-7 / 6 - 7j / 6) * np.sqrt(35) * OppC4_3_d6
            OppDt_site1 = OppDt_site1 + 7j / 3 * np.sqrt(10) * OppC4_2_d6 + (-7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_d6
            OppDt_site1 = OppDt_site1 + 49 / 6 * OppC40_d6
            OppDt_site1 = OppDt_site1 - 7j / 3 * np.sqrt(10) * OppC42_d6 + (7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_d6
            OppDt_site1 = OppDt_site1 + 7 / 6 * np.sqrt(35 / 2) * OppC44_d6 + (7 / 6 - 7j / 6) * np.sqrt(35) * OppC43_d6

            Hamiltonian_site1 = Hamiltonian + OpptenDq_site1 * tenDq + OppDs_site1 * Ds + OppDt_site1 * Dt

            Hamiltonian01 = Hamiltonian_site1 + OppSx_d6 * V0[0] + 1j * OppSy_d6 * V0[1] + OppSz_d6 * V0[2]
            Hamiltonian51 = Hamiltonian_site1 + OppSx_d6*V5[0] + 1j*OppSy_d6*V5[1] + OppSz_d6*V5[2]
            Hamiltonian501 = Hamiltonian_site1 + OppSx_d6 * V50[0] + 1j * OppSy_d6 * V50[1] + OppSz_d6 * V50[2]
            Hamiltonian901 = Hamiltonian_site1 + OppSx_d6 * V90[0] + 1j * OppSy_d6 * V90[1] + OppSz_d6 * V90[2]

            # site 2

            OpptenDq_site2 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_d6 + 3 / 10 * np.sqrt(35 / 2) * OppC44_d6
            OpptenDq_site2 = OpptenDq_site2 + 21 / 10 * OppC40_d6

            OppDs_site2 = -7j / np.sqrt(6) * OppC2_2_d6 + 7j / np.sqrt(6) * OppC22_d6
            OppDs_site2 = OppDs_site2 - (7 + 7j) / np.sqrt(6) * OppC2_1_d6 + (7 - 7j) / np.sqrt(6) * OppC21_d6

            OppDt_site2 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_d6 + (7 / 6 - 7j / 6) * np.sqrt(35) * OppC4_3_d6
            OppDt_site2 = OppDt_site2 - 7j / 3 * np.sqrt(10) * OppC4_2_d6 + (7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_d6
            OppDt_site2 = OppDt_site2 + 49 / 6 * OppC40_d6
            OppDt_site2 = OppDt_site2 + 7j / 3 * np.sqrt(10) * OppC42_d6 + (-7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_d6
            OppDt_site2 = OppDt_site2 + 7 / 6 * np.sqrt(35 / 2) * OppC44_d6 - (7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_d6

            Hamiltonian_site2 = Hamiltonian + OpptenDq_site2 * tenDq + OppDs_site2 * Ds + OppDt_site2 * Dt

            Hamiltonian02 = Hamiltonian_site2 + OppSx_d6 * V0[0] + 1j * OppSy_d6 * V0[1] + OppSz_d6 * V0[2]
            Hamiltonian52 = Hamiltonian_site2 + OppSx_d6*V5[0] + 1j*OppSy_d6*V5[1] + OppSz_d6*V5[2]
            Hamiltonian502 = Hamiltonian_site2 + OppSx_d6 * V50[0] + 1j * OppSy_d6 * V50[1] + OppSz_d6 * V50[2]
            Hamiltonian902 = Hamiltonian_site2 + OppSx_d6 * V90[0] + 1j * OppSy_d6 * V90[1] + OppSz_d6 * V90[2]

            # site 3

            OpptenDq_site3 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_d6 + 3 / 10 * np.sqrt(35 / 2) * OppC44_d6
            OpptenDq_site3 = OpptenDq_site3 + 21 / 10 * OppC40_d6

            OppDs_site3 = 7j / np.sqrt(6) * OppC2_2_d6 - 7j / np.sqrt(6) * OppC22_d6
            OppDs_site3 = OppDs_site3 - (7 - 7j) / np.sqrt(6) * OppC2_1_d6 + (7 + 7j) / np.sqrt(6) * OppC21_d6

            OppDt_site3 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_d6 + (7 / 6 + 7j / 6) * np.sqrt(35) * OppC4_3_d6
            OppDt_site3 = OppDt_site3 + 7j / 3 * np.sqrt(10) * OppC4_2_d6 + (7 / 6 - 7j / 6) * np.sqrt(5) * OppC4_1_d6
            OppDt_site3 = OppDt_site3 + 49 / 6 * OppC40_d6
            OppDt_site3 = OppDt_site3 - 7j / 3 * np.sqrt(10) * OppC42_d6 - (7 / 6 + 7j / 6) * np.sqrt(5) * OppC41_d6
            OppDt_site3 = OppDt_site3 + 7 / 6 * np.sqrt(35 / 2) * OppC44_d6 + (-7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_d6

            Hamiltonian_site3 = Hamiltonian + OpptenDq_site3 * tenDq + OppDs_site3 * Ds + OppDt_site3 * Dt

            Hamiltonian03 = Hamiltonian_site3 + OppSx_d6 * V0[0] + 1j * OppSy_d6 * V0[1] + OppSz_d6 * V0[2]
            Hamiltonian53 = Hamiltonian_site3 + OppSx_d6*V5[0] + 1j*OppSy_d6*V5[1] + OppSz_d6*V5[2]
            Hamiltonian503 = Hamiltonian_site3 + OppSx_d6 * V50[0] + 1j * OppSy_d6 * V50[1] + OppSz_d6 * V50[2]
            Hamiltonian903 = Hamiltonian_site3 + OppSx_d6 * V90[0] + 1j * OppSy_d6 * V90[1] + OppSz_d6 * V90[2]

            # site 4

            OpptenDq_site4 = 3 / 10 * np.sqrt(35 / 2) * OppC4_4_d6 + 3 / 10 * np.sqrt(35 / 2) * OppC44_d6
            OpptenDq_site4 = OpptenDq_site4 + 21 / 10 * OppC40_d6

            OppDs_site4 = -7j / np.sqrt(6) * OppC2_2_d6 + 7j / np.sqrt(6) * OppC22_d6
            OppDs_site4 = OppDs_site4 + (7 + 7j) / np.sqrt(6) * OppC2_1_d6 - (7 - 7j) / np.sqrt(6) * OppC21_d6

            OppDt_site4 = 7 / 6 * np.sqrt(35 / 2) * OppC4_4_d6 + (-7 / 6 + 7j / 6) * np.sqrt(35) * OppC4_3_d6
            OppDt_site4 = OppDt_site4 - 7j / 3 * np.sqrt(10) * OppC4_2_d6 - (7 / 6 + 7j / 6) * np.sqrt(5) * OppC4_1_d6
            OppDt_site4 = OppDt_site4 + 49 / 6 * OppC40_d6
            OppDt_site4 = OppDt_site4 + 7j / 3 * np.sqrt(10) * OppC42_d6 + (7 / 6 - 7j / 6) * np.sqrt(5) * OppC41_d6
            OppDt_site4 = OppDt_site4 + 7 / 6 * np.sqrt(35 / 2) * OppC44_d6 + (7 / 6 + 7j / 6) * np.sqrt(35) * OppC43_d6

            Hamiltonian_site4 = Hamiltonian + OpptenDq_site4 * tenDq + OppDs_site4 * Ds + OppDt_site4 * Dt

            Hamiltonian04 = Hamiltonian_site4 + OppSx_d6 * V0[0] + 1j * OppSy_d6 * V0[1] + OppSz_d6 * V0[2]
            Hamiltonian54 = Hamiltonian_site4 + OppSx_d6*V5[0] + 1j*OppSy_d6*V5[1] + OppSz_d6*V5[2]
            Hamiltonian504 = Hamiltonian_site4 + OppSx_d6 * V50[0] + 1j * OppSy_d6 * V50[1] + OppSz_d6 * V50[2]
            Hamiltonian904 = Hamiltonian_site4 + OppSx_d6 * V90[0] + 1j * OppSy_d6 * V90[1] + OppSz_d6 * V90[2]

            # sparse to full
            Hamiltonian01 = Hamiltonian01.todense()
            Hamiltonian51 = Hamiltonian51.todense()
            Hamiltonian501 = Hamiltonian501.todense()
            Hamiltonian901 = Hamiltonian901.todense()

            Hamiltonian02 = Hamiltonian02.todense()
            Hamiltonian52 = Hamiltonian52.todense()
            Hamiltonian502 = Hamiltonian502.todense()
            Hamiltonian902 = Hamiltonian902.todense()

            Hamiltonian03 = Hamiltonian03.todense()
            Hamiltonian53 = Hamiltonian53.todense()
            Hamiltonian503 = Hamiltonian503.todense()
            Hamiltonian903 = Hamiltonian903.todense()

            Hamiltonian04 = Hamiltonian04.todense()
            Hamiltonian54 = Hamiltonian54.todense()
            Hamiltonian504 = Hamiltonian504.todense()
            Hamiltonian904 = Hamiltonian904.todense()

            Ei01, CCi01 = eigh(Hamiltonian01)
            Ei51, CCi51 = eigh(Hamiltonian51)
            Ei501, CCi501 = eigh(Hamiltonian501)
            Ei901, CCi901 = eigh(Hamiltonian901)

            Ei02, CCi02 = eigh(Hamiltonian02)
            Ei52, CCi52 = eigh(Hamiltonian52)
            Ei502, CCi502 = eigh(Hamiltonian502)
            Ei902, CCi902 = eigh(Hamiltonian902)

            Ei03, CCi03 = eigh(Hamiltonian03)
            Ei53, CCi53 = eigh(Hamiltonian53)
            Ei503, CCi503 = eigh(Hamiltonian503)
            Ei903, CCi903 = eigh(Hamiltonian903)

            Ei04, CCi04 = eigh(Hamiltonian04)
            Ei54, CCi54 = eigh(Hamiltonian54)
            Ei504, CCi504 = eigh(Hamiltonian504)
            Ei904, CCi904 = eigh(Hamiltonian904)

            fB01 = np.zeros(Ei01.size, dtype='complex_')
            fB51 = np.zeros(Ei51.size,dtype='complex_')
            fB501 = np.zeros(Ei501.size, dtype='complex_')
            fB901 = np.zeros(Ei901.size, dtype='complex_')

            fB02 = np.zeros(Ei02.size, dtype='complex_')
            fB52 = np.zeros(Ei52.size, dtype='complex_')
            fB502 = np.zeros(Ei502.size, dtype='complex_')
            fB902 = np.zeros(Ei902.size, dtype='complex_')

            fB03 = np.zeros(Ei03.size, dtype='complex_')
            fB53 = np.zeros(Ei53.size, dtype='complex_')
            fB503 = np.zeros(Ei503.size, dtype='complex_')
            fB903 = np.zeros(Ei903.size, dtype='complex_')

            fB04 = np.zeros(Ei04.size, dtype='complex_')
            fB54 = np.zeros(Ei54.size, dtype='complex_')
            fB504 = np.zeros(Ei504.size, dtype='complex_')
            fB904 = np.zeros(Ei904.size, dtype='complex_')

            #Boltzman
            l = 0  # Variable auxialiar para recorrer fB y Ei

            for i in Ei01:
                fB01[l] = np.exp(-(i - Ei01[0]) / (8.61733035e-5 * Temp))
                fB51[l] = np.exp(-(Ei51[l] - Ei51[0]) / (8.61733035e-5 * Temp))
                fB501[l] = np.exp(-(Ei501[l] - Ei501[0]) / (8.61733035e-5 * Temp))
                fB901[l] = np.exp(-(Ei901[l] - Ei901[0]) / (8.61733035e-5 * Temp))

                fB02[l] = np.exp(-(Ei02[l] - Ei02[0]) / (8.61733035e-5 * Temp))
                fB52[l] = np.exp(-(Ei52[l] - Ei52[0]) / (8.61733035e-5 * Temp))
                fB502[l] = np.exp(-(Ei502[l] - Ei502[0]) / (8.61733035e-5 * Temp))
                fB902[l] = np.exp(-(Ei902[l] - Ei902[0]) / (8.61733035e-5 * Temp))

                fB03[l] = np.exp(-(Ei03[l] - Ei03[0]) / (8.61733035e-5 * Temp))
                fB53[l] = np.exp(-(Ei53[l] - Ei53[0]) / (8.61733035e-5 * Temp))
                fB503[l] = np.exp(-(Ei503[l] - Ei503[0]) / (8.61733035e-5 * Temp))
                fB903[l] = np.exp(-(Ei903[l] - Ei903[0]) / (8.61733035e-5 * Temp))

                fB04[l] = np.exp(-(Ei04[l] - Ei04[0]) / (8.61733035e-5 * Temp))
                fB54[l] = np.exp(-(Ei54[l] - Ei54[0]) / (8.61733035e-5 * Temp))
                fB504[l] = np.exp(-(Ei504[l] - Ei504[0]) / (8.61733035e-5 * Temp))
                fB904[l] = np.exp(-(Ei904[l] - Ei904[0]) / (8.61733035e-5 * Temp))

                l += 1

        #Transitions -Fermi Golden Rule
        #Polarizations
        if Nelec==5:
            opp_aux=(np.sqrt(3) + 1j) / np.sqrt(8) * OppT1_1_d5 + (1j-np.sqrt(3)) / np.sqrt(8) * OppT11_d5
            OppTH = np.array(opp_aux.todense()) #L30
            OppTr = np.array(OppT1_1_d5.todense()) #Circular a la derecha
            OppTl = np.array(OppT11_d5.todense()) #Circular a la izquierda
        else:
            opp_aux = (np.sqrt(3)+1j)/np.sqrt(8)*OppT1_1_d6 + (1j-np.sqrt(3))/np.sqrt(8)*OppT11_d6
            OppTH = np.array(opp_aux.todense())  # L30
            OppTr = np.array(OppT1_1_d6.todense())  # Circular a la derecha
            OppTl = np.array(OppT11_d6.todense())  # Circular a la izquierda

        if Sym == 1: #D4h

            fH0_ = np.array([])
            EnH0_ = np.array([])

            fH5_ = np.array([])
            EnH5_ = np.array([])

            fH50_ = np.array([])
            EnH50_ = np.array([])

            fH90_ = np.array([])
            EnH90_ = np.array([])

            fr_ = np.array([])
            Enr_ = np.array([])

            fl_ = np.array([])
            Enl_ = np.array([])

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei0:
                if fB0[k] < 0.01:
                    break

                Smed_H0 = abs(CCi0[:, k].conj().transpose() @ OppTH @ CCf0) ** 2 * fB0[k]
                Smed_r = abs(CCi0[:, k].conj().transpose() @ OppTr @ CCf0) ** 2 * fB0[k]
                Smed_l = abs(CCi0[:, k].conj().transpose() @ OppTl @ CCf0) ** 2 * fB0[k]

                fr_ = np.append(fr_, Smed_r[Smed_r > 1e-5].conj().transpose())
                Enr_ = np.append(Enr_, np.around((Eff0[Smed_r > 1e-5] - i) * 100) / 100)

                fl_ = np.append(fl_, Smed_l[Smed_l > 1e-5].conj().transpose())
                Enl_ = np.append(Enl_, np.around((Eff0[Smed_l > 1e-5] - i) * 100) / 100)

                fH0_ = np.append(fH0_, Smed_H0[Smed_H0 > 1e-5].conj().transpose())
                EnH0_ = np.append(EnH0_, np.around((Eff0[Smed_H0 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH0_, return_index=True)
            EnH0 = EnH0_[np.sort(b)]
            fH0 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH0:
                for j in EnH0_:
                    if i == j:
                        aux = aux + fH0_[aux1]
                    aux1 += 1
                fH0 = np.append(fH0, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enr_, return_index=True)
            Enr = Enr_[np.sort(b)]
            fr = np.array([])
            aux = 0
            aux1 = 0

            for i in Enr:
                for j in Enr_:
                    if i == j:
                        aux = aux + fr_[aux1]
                    aux1 += 1
                fr = np.append(fr, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enl_, return_index=True)
            Enl = Enl_[np.sort(b)]
            fl = np.array([])
            aux = 0
            aux1 = 0

            for i in Enl:
                for j in Enl_:
                    if i == j:
                        aux = aux + fl_[aux1]
                    aux1 += 1
                fl = np.append(fl, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei5:
                if fB5[k] < 0.01:
                    break
                Smed_H5 = abs(CCi5[:, k].conj().transpose() @ OppTH @ CCf5) ** 2 * fB5[k]
                fH5_ = np.append(fH5_, Smed_H5[Smed_H5 > 1e-5].conj().transpose())
                EnH5_ = np.append(EnH5_, np.around((Eff5[Smed_H5 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH5_, return_index=True)
            EnH5 = EnH5_[np.sort(b)]
            fH5 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH5:
                for j in EnH5_:
                    if i == j:
                        aux = aux + fH5_[aux1]
                    aux1 += 1
                fH5 = np.append(fH5, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei50:
                if fB50[k] < 0.01:
                    break

                Smed_H50 = abs(CCi50[:, k].conj().transpose() @ OppTH @ CCf50) ** 2 * fB50[k]
                fH50_ = np.append(fH50_, Smed_H50[Smed_H50 > 1e-5].conj().transpose())
                EnH50_ = np.append(EnH50_, np.around((Eff50[Smed_H50 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH50_, return_index=True)
            EnH50 = EnH50_[np.sort(b)]
            fH50 = np.array([])
            aux = 0
            aux1 = 0

            for i in EnH50:
                for j in EnH50_:
                    if i == j:
                        aux = aux + fH50_[aux1]
                    aux1 += 1
                fH50 = np.append(fH50, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei90:
                if fB90[k] < 0.01:
                    break

                Smed_H90 = abs(CCi90[:, k].conj().transpose() @ OppTH @ CCf90) ** 2 * fB90[k]
                fH90_ = np.append(fH90_, Smed_H90[Smed_H90 > 1e-5].conj().transpose())
                EnH90_ = np.append(EnH90_, np.around((Eff90[Smed_H90 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH90_, return_index=True)
            EnH90 = EnH90_[np.sort(b)]
            fH90 = np.array([])
            aux = 0
            aux1 = 0

            for i in EnH90:
                for j in EnH90_:
                    if i == j:
                        aux = aux + fH90_[aux1]
                    aux1 += 1
                fH90 = np.append(fH90, aux)
                aux = 0
                aux1 = 0
        else:
            #D3d
            fH01_ = np.array([])
            EnH01_ = np.array([])

            fH02_ = np.array([])
            EnH02_ = np.array([])

            fH03_ = np.array([])
            EnH03_ = np.array([])

            fH04_ = np.array([])
            EnH04_ = np.array([])

            fH51_ = np.array([])
            EnH51_ = np.array([])

            fH52_ = np.array([])
            EnH52_ = np.array([])

            fH53_ = np.array([])
            EnH53_ = np.array([])

            fH54_ = np.array([])
            EnH54_ = np.array([])

            fH501_ = np.array([])
            EnH501_ = np.array([])

            fH502_ = np.array([])
            EnH502_ = np.array([])

            fH503_ = np.array([])
            EnH503_ = np.array([])

            fH504_ = np.array([])
            EnH504_ = np.array([])

            fH901_ = np.array([])
            EnH901_ = np.array([])

            fH902_ = np.array([])
            EnH902_ = np.array([])

            fH903_ = np.array([])
            EnH903_ = np.array([])

            fH904_ = np.array([])
            EnH904_ = np.array([])

            fr1_ = np.array([])
            Enr1_ = np.array([])

            fr2_ = np.array([])
            Enr2_ = np.array([])

            fr3_ = np.array([])
            Enr3_ = np.array([])

            fr4_ = np.array([])
            Enr4_ = np.array([])

            fl1_ = np.array([])
            Enl1_ = np.array([])

            fl2_ = np.array([])
            Enl2_ = np.array([])

            fl3_ = np.array([])
            Enl3_ = np.array([])

            fl4_ = np.array([])
            Enl4_ = np.array([])

            #Sitio 1

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei01:
                if fB01[k] < 0.01:
                    break

                Smed_H01 = abs(CCi01[:, k].conj().transpose() @ OppTH @ CCf01) ** 2 * fB01[k]
                Smed_r1 = abs(CCi01[:, k].conj().transpose() @ OppTr @ CCf01) ** 2 * fB01[k]
                Smed_l1 = abs(CCi01[:, k].conj().transpose() @ OppTl @ CCf01) ** 2 * fB01[k]

                fH01_ = np.append(fH01_, Smed_H01[Smed_H01 > 1e-5].conj().transpose())
                EnH01_ = np.append(EnH01_, np.around((Eff01[Smed_H01 > 1e-5] - i) * 100)/100)

                fr1_ = np.append(fr1_, Smed_r1[Smed_r1 > 1e-5].conj().transpose())
                Enr1_ = np.append(Enr1_, np.around((Eff01[Smed_r1 > 1e-5] - i) * 100) / 100)

                fl1_ = np.append(fl1_, Smed_l1[Smed_l1 > 1e-5].conj().transpose())
                Enl1_ = np.append(Enl1_, np.around((Eff01[Smed_l1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH01_, return_index=True)
            EnH01 = EnH01_[np.sort(b)]
            fH01 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH01:
                for j in EnH01_:
                    if i == j:
                        aux = aux + fH01_[aux1]
                    aux1 += 1
                fH01 = np.append(fH01, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enr1_, return_index=True)
            Enr1 = Enr1_[np.sort(b)]
            fr1 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enr1:
                for j in Enr1_:
                    if i == j:
                        aux = aux + fr1_[aux1]
                    aux1 += 1
                fr1 = np.append(fr1, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enl1_, return_index=True)
            Enl1 = Enl1_[np.sort(b)]
            fl1 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enl1:
                for j in Enl1_:
                    if i == j:
                        aux = aux + fl1_[aux1]
                    aux1 += 1
                fl1 = np.append(fl1, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei51:
                if fB51[k] < 0.01:
                    break

                Smed_H51 = abs(CCi51[:, k].conj().transpose() @ OppTH @ CCf51) ** 2 * fB51[k]
                fH51_ = np.append(fH51_, Smed_H51[Smed_H51 > 1e-5].conj().transpose())
                EnH51_ = np.append(EnH51_, np.around((Eff51[Smed_H51 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH51_, return_index=True)
            EnH51 = EnH51_[np.sort(b)]
            fH51 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH51:
                for j in EnH51_:
                    if i == j:
                        aux = aux + fH51_[aux1]
                    aux1 += 1
                fH51 = np.append(fH51, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei501:
                if fB501[k] < 0.01:
                    break

                Smed_H501 = abs(CCi501[:, k].conj().transpose() @ OppTH @ CCf501) ** 2 * fB501[k]
                fH501_ = np.append(fH501_, Smed_H501[Smed_H501 > 1e-5].conj().transpose())
                EnH501_ = np.append(EnH501_, np.around((Eff501[Smed_H501 > 1e-5] - i) * 100)/100)

                k += 1

            a, b = np.unique(EnH501_, return_index=True)
            EnH501 = EnH501_[np.sort(b)]
            fH501 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH501:
                for j in EnH501_:
                    if i == j:
                        aux = aux + fH501_[aux1]
                    aux1 += 1
                fH501 = np.append(fH501, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei901:
                if fB901[k] < 0.01:
                    break

                Smed_H901 = abs(CCi901[:, k].conj().transpose() @ OppTH @ CCf901) ** 2 * fB901[k]
                fH901_ = np.append(fH901_, Smed_H901[Smed_H901 > 1e-5].conj().transpose())
                EnH901_ = np.append(EnH901_, np.around((Eff901[Smed_H901 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH901_, return_index=True)
            EnH901 = EnH901_[np.sort(b)]
            fH901 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH901:
                for j in EnH901_:
                    if i == j:
                        aux = aux + fH901_[aux1]
                    aux1 += 1
                fH901 = np.append(fH901, aux)
                aux = 0
                aux1 = 0

            # Sitio 2

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei02:
                if fB02[k] < 0.01:
                    break

                Smed_H02 = abs(CCi02[:, k].conj().transpose() @ OppTH @ CCf02) ** 2 * fB02[k]
                Smed_r2 = abs(CCi02[:, k].conj().transpose() @ OppTr @ CCf02) ** 2 * fB02[k]
                Smed_l2 = abs(CCi02[:, k].conj().transpose() @ OppTl @ CCf02) ** 2 * fB02[k]

                fH02_ = np.append(fH02_, Smed_H02[Smed_H02 > 1e-5].conj().transpose())
                EnH02_ = np.append(EnH02_, np.around((Eff02[Smed_H02 > 1e-5] - i) * 100) / 100)

                fr2_ = np.append(fr2_, Smed_r2[Smed_r2 > 1e-5].conj().transpose())
                Enr2_ = np.append(Enr2_, np.around((Eff02[Smed_r2 > 1e-5] - i) * 100) / 100)

                fl2_ = np.append(fl2_, Smed_l2[Smed_l2 > 1e-5].conj().transpose())
                Enl2_ = np.append(Enl2_, np.around((Eff02[Smed_l2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH02_, return_index=True)
            EnH02 = EnH02_[np.sort(b)]
            fH02 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH02:
                for j in EnH02_:
                    if i == j:
                        aux = aux + fH02_[aux1]
                    aux1 += 1
                fH02 = np.append(fH02, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enr2_, return_index=True)
            Enr2 = Enr2_[np.sort(b)]
            fr2 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enr2:
                for j in Enr2_:
                    if i == j:
                        aux = aux + fr2_[aux1]
                    aux1 += 1
                fr2 = np.append(fr2, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enl2_, return_index=True)
            Enl2 = Enl2_[np.sort(b)]
            fl2 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enl2:
                for j in Enl2_:
                    if i == j:
                        aux = aux + fl2_[aux1]
                    aux1 += 1
                fl2 = np.append(fl2, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei52:
                if fB52[k] < 0.01:
                    break

                Smed_H52 = abs(CCi52[:, k].conj().transpose() @ OppTH @ CCf52) ** 2 * fB52[k]
                fH52_ = np.append(fH52_, Smed_H52[Smed_H52 > 1e-5].conj().transpose())
                EnH52_ = np.append(EnH52_, np.around((Eff52[Smed_H52 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH52_, return_index=True)
            EnH52 = EnH52_[np.sort(b)]
            fH52 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH52:
                for j in EnH52_:
                    if i == j:
                        aux = aux + fH52_[aux1]
                    aux1 += 1
                fH52 = np.append(fH52, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei502:
                if fB502[k] < 0.01:
                    break

                Smed_H502 = abs(CCi502[:, k].conj().transpose() @ OppTH @ CCf502) ** 2 * fB502[k]
                fH502_ = np.append(fH502_, Smed_H502[Smed_H502 > 1e-5].conj().transpose())
                EnH502_ = np.append(EnH502_, np.around((Eff502[Smed_H502 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH502_, return_index=True)
            EnH502 = EnH502_[np.sort(b)]
            fH502 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH502:
                for j in EnH502_:
                    if i == j:
                        aux = aux + fH502_[aux1]
                    aux1 += 1
                fH502 = np.append(fH502, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei902:
                if fB902[k] < 0.01:
                    break

                Smed_H902 = abs(CCi902[:, k].conj().transpose() @ OppTH @ CCf902) ** 2 * fB902[k]
                fH902_ = np.append(fH902_, Smed_H902[Smed_H902 > 1e-5].conj().transpose())
                EnH902_ = np.append(EnH902_, np.around((Eff902[Smed_H902 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH902_, return_index=True)
            EnH902 = EnH902_[np.sort(b)]
            fH902 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH902:
                for j in EnH902_:
                    if i == j:
                        aux = aux + fH902_[aux1]
                    aux1 += 1
                fH902 = np.append(fH902, aux)
                aux = 0
                aux1 = 0

            # Sitio 3

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei03:
                if fB03[k] < 0.01:
                    break

                Smed_H03 = abs(CCi03[:, k].conj().transpose() @ OppTH @ CCf03) ** 2 * fB03[k]
                Smed_r3 = abs(CCi03[:, k].conj().transpose() @ OppTr @ CCf03) ** 2 * fB03[k]
                Smed_l3 = abs(CCi03[:, k].conj().transpose() @ OppTl @ CCf03) ** 2 * fB03[k]

                fH03_ = np.append(fH03_, Smed_H03[Smed_H03 > 1e-5].conj().transpose())
                EnH03_ = np.append(EnH03_, np.around((Eff03[Smed_H03 > 1e-5] - i) * 100) / 100)

                fr3_ = np.append(fr3_, Smed_r3[Smed_r3 > 1e-5].conj().transpose())
                Enr3_ = np.append(Enr3_, np.around((Eff03[Smed_r3 > 1e-5] - i) * 100) / 100)

                fl3_ = np.append(fl3_, Smed_l3[Smed_l3 > 1e-5].conj().transpose())
                Enl3_ = np.append(Enl3_, np.around((Eff03[Smed_l3 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH03_, return_index=True)
            EnH03 = EnH03_[np.sort(b)]
            fH03 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH03:
                for j in EnH03_:
                    if i == j:
                        aux = aux + fH03_[aux1]
                    aux1 += 1
                fH03 = np.append(fH03, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enr3_, return_index=True)
            Enr3 = Enr3_[np.sort(b)]
            fr3 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enr3:
                for j in Enr3_:
                    if i == j:
                        aux = aux + fr3_[aux1]
                    aux1 += 1
                fr3 = np.append(fr3, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enl3_, return_index=True)
            Enl3 = Enl3_[np.sort(b)]
            fl3 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enl3:
                for j in Enl3_:
                    if i == j:
                        aux = aux + fl3_[aux1]
                    aux1 += 1
                fl3 = np.append(fl3, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei53:
                if fB53[k] < 0.01:
                    break

                Smed_H53 = abs(CCi53[:, k].conj().transpose() @ OppTH @ CCf53) ** 2 * fB53[k]
                fH53_ = np.append(fH53_, Smed_H53[Smed_H53 > 1e-5].conj().transpose())
                EnH53_ = np.append(EnH53_, np.around((Eff53[Smed_H53 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH53_, return_index=True)
            EnH53 = EnH53_[np.sort(b)]
            fH53 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH53:
                for j in EnH53_:
                    if i == j:
                        aux = aux + fH53_[aux1]
                    aux1 += 1
                fH53 = np.append(fH53, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei503:
                if fB503[k] < 0.01:
                    break

                Smed_H503 = abs(CCi503[:, k].conj().transpose() @ OppTH @ CCf503) ** 2 * fB503[k]
                fH503_ = np.append(fH503_, Smed_H503[Smed_H503 > 1e-5].conj().transpose())
                EnH503_ = np.append(EnH503_, np.around((Eff503[Smed_H503 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH503_, return_index=True)
            EnH503 = EnH503_[np.sort(b)]
            fH503 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH503:
                for j in EnH503_:
                    if i == j:
                        aux = aux + fH503_[aux1]
                    aux1 += 1
                fH503 = np.append(fH503, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei903:
                if fB903[k] < 0.01:
                    break

                Smed_H903 = abs(CCi903[:, k].conj().transpose() @ OppTH @ CCf903) ** 2 * fB903[k]
                fH903_ = np.append(fH903_, Smed_H903[Smed_H903 > 1e-5].conj().transpose())
                EnH903_ = np.append(EnH903_, np.around((Eff903[Smed_H903 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH903_, return_index=True)
            EnH903 = EnH903_[np.sort(b)]
            fH903 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH903:
                for j in EnH903_:
                    if i == j:
                        aux = aux + fH903_[aux1]
                    aux1 += 1
                fH903 = np.append(fH903, aux)
                aux = 0
                aux1 = 0

            # Sitio 4

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei04:
                if fB04[k] < 0.01:
                    break

                Smed_H04 = abs(CCi04[:, k].conj().transpose() @ OppTH @ CCf04) ** 2 * fB04[k]
                Smed_r4 = abs(CCi04[:, k].conj().transpose() @ OppTr @ CCf04) ** 2 * fB04[k]
                Smed_l4 = abs(CCi04[:, k].conj().transpose() @ OppTl @ CCf04) ** 2 * fB04[k]

                fH04_ = np.append(fH04_, Smed_H04[Smed_H04 > 1e-5].conj().transpose())
                EnH04_ = np.append(EnH04_, np.around((Eff04[Smed_H04 > 1e-5] - i) * 100) / 100)

                fr4_ = np.append(fr4_, Smed_r4[Smed_r4 > 1e-5].conj().transpose())
                Enr4_ = np.append(Enr4_, np.around((Eff04[Smed_r4 > 1e-5] - i) * 100) / 100)

                fl4_ = np.append(fl4_, Smed_l4[Smed_l4 > 1e-5].conj().transpose())
                Enl4_ = np.append(Enl4_, np.around((Eff04[Smed_l4 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH04_, return_index=True)
            EnH04 = EnH04_[np.sort(b)]
            fH04 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH04:
                for j in EnH04_:
                    if i == j:
                        aux = aux + fH04_[aux1]
                    aux1 += 1
                fH04 = np.append(fH04, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enr4_, return_index=True)
            Enr4 = Enr4_[np.sort(b)]
            fr4 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enr4:
                for j in Enr4_:
                    if i == j:
                        aux = aux + fr4_[aux1]
                    aux1 += 1
                fr4 = np.append(fr4, aux)
                aux = 0
                aux1 = 0

            a, b = np.unique(Enl4_, return_index=True)
            Enl4 = Enl4_[np.sort(b)]
            fl4 = np.array([])
            aux = 0
            aux1 = 0

            for i in Enl4:
                for j in Enl4_:
                    if i == j:
                        aux = aux + fl4_[aux1]
                    aux1 += 1
                fl4 = np.append(fl4, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei54:
                if fB54[k] < 0.01:
                    break

                Smed_H54 = abs(CCi54[:, k].conj().transpose() @ OppTH @ CCf54) ** 2 * fB54[k]
                fH54_ = np.append(fH54_, Smed_H54[Smed_H54 > 1e-5].conj().transpose())
                EnH54_ = np.append(EnH54_, np.around((Eff54[Smed_H54 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH54_, return_index=True)
            EnH54 = EnH54_[np.sort(b)]
            fH54 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH54:
                for j in EnH54_:
                    if i == j:
                        aux = aux + fH54_[aux1]
                    aux1 += 1
                fH54 = np.append(fH54, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei504:
                if fB504[k] < 0.01:
                    break

                Smed_H504 = abs(CCi504[:, k].conj().transpose() @ OppTH @ CCf504) ** 2 * fB504[k]
                fH504_ = np.append(fH504_, Smed_H504[Smed_H504 > 1e-5].conj().transpose())
                EnH504_ = np.append(EnH504_, np.around((Eff504[Smed_H504 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH504_, return_index=True)
            EnH504 = EnH504_[np.sort(b)]
            fH504 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH504:
                for j in EnH504_:
                    if i == j:
                        aux = aux + fH504_[aux1]
                    aux1 += 1
                fH504 = np.append(fH504, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei904:
                if fB904[k] < 0.01:
                    break

                Smed_H904 = abs(CCi904[:, k].conj().transpose() @ OppTH @ CCf904) ** 2 * fB904[k]
                fH904_ = np.append(fH904_, Smed_H904[Smed_H904 > 1e-5].conj().transpose())
                EnH904_ = np.append(EnH904_, np.around((Eff904[Smed_H904 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH904_, return_index=True)
            EnH904 = EnH904_[np.sort(b)]
            fH904 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH904:
                for j in EnH904_:
                    if i == j:
                        aux = aux + fH904_[aux1]
                    aux1 += 1
                fH904 = np.append(fH904, aux)
                aux = 0
                aux1 = 0

        #crear diccionarios
        if Sym == 1:
            # D4h contenedores
            Cont[Container - 1] = {
                "Id": Id,
                "fH0": fH0,
                "EnH0": EnH0,
                "fH5": fH5,
                "EnH5": EnH5,
                "fH50": fH50,
                "EnH50": EnH50,
                "fH90": fH90,
                "EnH90": EnH90,
                "fr": fr,
                "Enr": Enr,
                "fl": fl,
                "Enl": Enl
                }
        else:
            #D3d
            Cont[Container - 1] = {
                "Id": Id,

                "fH01": fH01,
                "EnH01": EnH01,
                "fH51": fH51,
                "EnH51": EnH51,
                "fH501": fH501,
                "EnH501": EnH501,
                "fH901": fH901,
                "EnH901": EnH901,
                "fr1": fr1,
                "Enr1": Enr1,
                "fl1": fl1,
                "Enl1": Enl1,

                "fH02": fH02,
                "EnH02": EnH02,
                "fH52": fH52,
                "EnH52": EnH52,
                "fH502": fH502,
                "EnH502": EnH502,
                "fH902": fH902,
                "EnH902": EnH902,
                "fr2": fr2,
                "Enr2": Enr2,
                "fl2": fl2,
                "Enl2": Enl2,

                "fH03": fH03,
                "EnH03": EnH03,
                "fH53": fH53,
                "EnH53": EnH53,
                "fH503": fH503,
                "EnH503": EnH503,
                "fH903": fH903,
                "EnH903": EnH903,
                "fr3": fr3,
                "Enr3": Enr3,
                "fl3": fl3,
                "Enl3": Enl3,

                "fH04": fH04,
                "EnH04": EnH04,
                "fH54": fH54,
                "EnH54": EnH54,
                "fH504": fH504,
                "EnH504": EnH504,
                "fH904": fH904,
                "EnH904": EnH904,
                "fr4": fr4,
                "Enr4": Enr4,
                "fl4": fl4,
                "Enl4": Enl4
            }

        # salvar
        np.save(current_direction + '/temp_sim_Fe.npy', Cont)

    Cont = np.load(current_direction + '/temp_sim_Fe.npy', allow_pickle=True)

    if Sym == 1:
        # D4h
        fH0 = Cont[Container - 1]['fH0']
        EnH0 = Cont[Container - 1]['EnH0'] + ShiftE
        fH5 = Cont[Container - 1]['fH5']
        EnH5 = Cont[Container - 1]['EnH5'] + ShiftE
        fH50 = Cont[Container - 1]['fH50']
        EnH50 = Cont[Container - 1]['EnH50'] + ShiftE
        fH90 = Cont[Container - 1]['fH90']
        EnH90 = Cont[Container - 1]['EnH90'] + ShiftE
        fr = Cont[Container - 1]['fr']
        Enr = Cont[Container - 1]['Enr'] + ShiftE
        fl = Cont[Container - 1]['fl']
        Enl = Cont[Container - 1]['Enl'] + ShiftE
    else:
        #D3d
        fH01 = Cont[Container - 1]['fH01']
        EnH01 = Cont[Container - 1]['EnH01'] + ShiftE
        fH51 = Cont[Container - 1]['fH51']
        EnH51 = Cont[Container - 1]['EnH51'] + ShiftE
        fH501 = Cont[Container - 1]['fH501']
        EnH501 = Cont[Container - 1]['EnH501'] + ShiftE
        fH901 = Cont[Container - 1]['fH901']
        EnH901 = Cont[Container - 1]['EnH901'] + ShiftE
        fr1 = Cont[Container - 1]['fr1']
        Enr1 = Cont[Container - 1]['Enr1'] + ShiftE
        fl1 = Cont[Container - 1]['fl1']
        Enl1 = Cont[Container - 1]['Enl1'] + ShiftE

        fH02 = Cont[Container - 1]['fH02']
        EnH02 = Cont[Container - 1]['EnH02'] + ShiftE
        fH52 = Cont[Container - 1]['fH52']
        EnH52 = Cont[Container - 1]['EnH52'] + ShiftE
        fH502 = Cont[Container - 1]['fH502']
        EnH502 = Cont[Container - 1]['EnH502'] + ShiftE
        fH902 = Cont[Container - 1]['fH902']
        EnH902 = Cont[Container - 1]['EnH902'] + ShiftE
        fr2 = Cont[Container - 1]['fr2']
        Enr2 = Cont[Container - 1]['Enr2'] + ShiftE
        fl2 = Cont[Container - 1]['fl2']
        Enl2 = Cont[Container - 1]['Enl2'] + ShiftE

        fH03 = Cont[Container - 1]['fH03']
        EnH03 = Cont[Container - 1]['EnH03'] + ShiftE
        fH53 = Cont[Container - 1]['fH53']
        EnH53 = Cont[Container - 1]['EnH53'] + ShiftE
        fH503 = Cont[Container - 1]['fH503']
        EnH503 = Cont[Container - 1]['EnH503'] + ShiftE
        fH903 = Cont[Container - 1]['fH903']
        EnH903 = Cont[Container - 1]['EnH903'] + ShiftE
        fr3 = Cont[Container - 1]['fr3']
        Enr3 = Cont[Container - 1]['Enr3'] + ShiftE
        fl3 = Cont[Container - 1]['fl3']
        Enl3 = Cont[Container - 1]['Enl3'] + ShiftE

        fH04 = Cont[Container - 1]['fH04']
        EnH04 = Cont[Container - 1]['EnH04'] + ShiftE
        fH54 = Cont[Container - 1]['fH54']
        EnH54 = Cont[Container - 1]['EnH54'] + ShiftE
        fH504 = Cont[Container - 1]['fH504']
        EnH504 = Cont[Container - 1]['EnH504'] + ShiftE
        fH904 = Cont[Container - 1]['fH904']
        EnH904 = Cont[Container - 1]['EnH904'] + ShiftE
        fr4 = Cont[Container - 1]['fr4']
        Enr4 = Cont[Container - 1]['Enr4'] + ShiftE
        fl4 = Cont[Container - 1]['fl4']
        Enl4 = Cont[Container - 1]['Enl4'] + ShiftE

    #H
    x = np.reshape(x, x.size)
    if Sym==1:
        FFH0 = S2B.Sticks2Band(x, EnH0, fH0, MinL, MaxL, SplitE, 0, Gauss)
        FFH5 = S2B.Sticks2Band(x, EnH5, fH5, MinL, MaxL, SplitE, 0, Gauss)
        FFH50 = S2B.Sticks2Band(x, EnH50, fH50, MinL, MaxL, SplitE, 0, Gauss)
        FFH90 = S2B.Sticks2Band(x, EnH90, fH90, MinL, MaxL, SplitE, 0, Gauss)
        FFr = S2B.Sticks2Band(x, Enr, fr, MinL, MaxL, SplitE, 0, Gauss)
        FFl = S2B.Sticks2Band(x, Enl, fl, MinL, MaxL, SplitE, 0, Gauss)
    else:
        #Sitio 1
        FFH01 = S2B.Sticks2Band(x, EnH01, fH01, MinL, MaxL, SplitE, 0, Gauss)
        FFH51 = S2B.Sticks2Band(x, EnH51, fH51, MinL, MaxL, SplitE, 0, Gauss)
        FFH501 = S2B.Sticks2Band(x, EnH501, fH501, MinL, MaxL, SplitE, 0, Gauss)
        FFH901 = S2B.Sticks2Band(x, EnH901, fH901, MinL, MaxL, SplitE, 0, Gauss)
        #Right
        FFr1 = S2B.Sticks2Band(x, Enr1, fr1, MinL, MaxL, SplitE, 0, Gauss)
        #Left
        FFl1 = S2B.Sticks2Band(x, Enl1, fl1, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 2
        FFH02 = S2B.Sticks2Band(x, EnH02, fH02, MinL, MaxL, SplitE, 0, Gauss)
        FFH52 = S2B.Sticks2Band(x, EnH52, fH52, MinL, MaxL, SplitE, 0, Gauss)
        FFH502 = S2B.Sticks2Band(x, EnH502, fH502, MinL, MaxL, SplitE, 0, Gauss)
        FFH902 = S2B.Sticks2Band(x, EnH902, fH902, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr2 = S2B.Sticks2Band(x, Enr2, fr2, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl2 = S2B.Sticks2Band(x, Enl2, fl2, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 3
        FFH03 = S2B.Sticks2Band(x, EnH03, fH03, MinL, MaxL, SplitE, 0, Gauss)
        FFH53 = S2B.Sticks2Band(x, EnH53, fH53, MinL, MaxL, SplitE, 0, Gauss)
        FFH503 = S2B.Sticks2Band(x, EnH503, fH503, MinL, MaxL, SplitE, 0, Gauss)
        FFH903 = S2B.Sticks2Band(x, EnH903, fH903, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr3 = S2B.Sticks2Band(x, Enr3, fr3, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl3 = S2B.Sticks2Band(x, Enl3, fl3, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 4
        FFH04 = S2B.Sticks2Band(x, EnH04, fH04, MinL, MaxL, SplitE, 0, Gauss)
        FFH54 = S2B.Sticks2Band(x, EnH54, fH54, MinL, MaxL, SplitE, 0, Gauss)
        FFH504 = S2B.Sticks2Band(x, EnH504, fH504, MinL, MaxL, SplitE, 0, Gauss)
        FFH904 = S2B.Sticks2Band(x, EnH904, fH904, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr4 = S2B.Sticks2Band(x, Enr4, fr4, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl4 = S2B.Sticks2Band(x, Enl4, fl4, MinL, MaxL, SplitE, 0, Gauss)

    if Sym==1:
        # D4h
        FFXAS = (FFr + FFl)/2
        A = np.trapz(FFXAS, x)
        FFXAS = FFXAS / A * 3 / 7 * (10 - Nelec)
        FFXAS = FFXAS * XASInt

        FFXMCD = FFl - FFr
        FFXMCD = FFXMCD / A * 3 / 7 * (10 - Nelec)
        FFXMCD = FFXMCD * DichrInt

        FFXMLD50 = FFH50 - FFH5
        FFXMLD50 = FFXMLD50 / A * 3 / 7 * (10 - Nelec)

        FFXMLD50 = FFXMLD50 * XMLD_50

        FFXMLD90 = FFH90 - FFH0
        FFXMLD90 = FFXMLD90 / A * 3 / 7 * (10 - Nelec)

        FFXMLD90 = FFXMLD90 * XMLD_90

        mask = np.logical_and(x > 1000, x < 2000)

        xd = x[mask] - 1000
        L = pchip(x, FFXMLD50)
        FF3d = L(xd)

        mask = np.logical_and(x>2000,x<3000)
        xd = x[mask] - 2000
        K = pchip(x, FFXMLD90)
        FF4d = K(xd)

        xd = x[x>3000]-3000
        M = pchip(x,FFXAS)
        FF5d=M(xd)

        FF = FFXMCD[x < 1000]
        FF = np.append(FF, FF3d)
        FF = np.append(FF, FF4d)
        FF = np.append(FF,FF5d)

        return FF.real

    else:
        # D3d
        FFXAS = (FFr1 + FFl1 + FFr2 + FFl2 + FFr3 + FFl3 + FFr4 + FFl4)/8
        A = np.trapz(FFXAS, x)
        FFXAS = FFXAS / A * 3 / 7 * (10 - Nelec)
        FFXAS = FFXAS * XASInt

        FFXMCD = -(FFr1 - FFl1 + FFr2 - FFl2 + FFr3 - FFl3 + FFr4 - FFl4)/4
        FFXMCD = FFXMCD / A * 3 / 7 * (10 - Nelec)
        FFXMCD = FFXMCD * DichrInt

        FFXMLD50 = (FFH501 - FFH51+FFH502 - FFH52+FFH503 - FFH53+FFH504 - FFH54)/4
        FFXMLD50 = FFXMLD50 / A * 3 / 7 * (10 - Nelec)

        FFXMLD50 = FFXMLD50 * XMLD_50

        FFXMLD90 = (FFH901 - FFH01+FFH902 - FFH02+FFH903 - FFH03+FFH904 - FFH04)/4
        FFXMLD90 = FFXMLD90 / A * 3 / 7 * (10 - Nelec)

        FFXMLD90 = FFXMLD90 * XMLD_90

        mask = np.logical_and(x > 1000, x < 2000)

        xd = x[mask] - 1000
        L = pchip(x, FFXMLD50)
        FF3d = L(xd)

        mask = np.logical_and(x>2000,x<3000)
        xd = x[mask] - 2000
        K = pchip(x, FFXMLD90)
        FF4d = K(xd)

        xd = x[x>3000]-3000
        M=pchip(x,FFXAS)
        FF5d = M(xd)

        FF = FFXMCD[x < 1000]
        FF = np.append(FF, FF3d)
        FF = np.append(FF, FF4d)
        FF = np.append(FF,FF5d)

        return FF.real
