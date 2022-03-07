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

def ctm4fitLFePol(x, Nelec, Red_F2, Red_F2pd,tenDq, Ds, Dt, Bz, Temp, Sym, MinL, MaxL, Gauss, SplitE , ShiftE, Container, XASInt, DichrInt, XMLD1, XMLD2, Phi_1, Phi_2, Phi_Ref1, Phi_Ref2, Pol):
    Cont = np.array([0, 0, 0], dtype=object)
    mm = 0

    current_direction = os.path.dirname(os.path.realpath(__file__))

    Strm = np.array([Nelec, Red_F2, Red_F2pd, tenDq, Ds, Dt, Bz, Temp, Sym,Phi_1, Phi_2, Phi_Ref1, Phi_Ref2, Pol])

    try:
        str4cont = current_direction + "/temp_sim_Fe.npy"
        Cont = np.load(str4cont, allow_pickle=True)
        pass
    except:
        pass

    try:
        Id = Cont[Container - 1]['Id']

        if Strm[0] == Id[0] and Strm[1] == Id[1] and Strm[2] == Id[2] and Strm[3] == Id[3] and Strm[4] == Id[4] and \
                Strm[5] == Id[5] and Strm[6] == Id[6] and Strm[7] == Id[7] and Strm[8] == Id[8] and \
                Strm[9] == Id[9] and Strm[10] == Id[10] and Strm[11] == Id[11] and Strm[12] == Id[12] and \
                Strm[13] == Id[13]:
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


        # =========================Final state
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

        VPhi1 = np.array([0.9962*np.cos(np.radians(Phi_1)), 0.9962*np.sin(np.radians(Phi_1)), 0.0867]) * Bz
        VPhi2 = np.array([0.9962*np.cos(np.radians(Phi_2)), 0.9962*np.sin(np.radians(Phi_2)), 0.0867]) * Bz
        VPhiRef1 = np.array([0.9962*np.cos(np.radians(Phi_Ref1)), 0.9962*np.sin(np.radians(Phi_Ref1)), 0.0867]) * Bz
        VPhiRef2 = np.array([0.9962*np.cos(np.radians(Phi_Ref2)), 0.9962*np.sin(np.radians(Phi_Ref2)), 0.0867]) * Bz

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

                Hamiltonian_Phi1 = Hamiltonian + OppSz_p5d6 * VPhi1[2] + 1j * OppSy_p5d6 * VPhi1[1] + OppSx_p5d6 * VPhi1[0]
                Hamiltonian_Phi2 = Hamiltonian + OppSz_p5d6*VPhi2[2] + 1j*OppSy_p5d6*VPhi2[1] + OppSx_p5d6*VPhi2[0]
                Hamiltonian_PhiRef1 = Hamiltonian + OppSz_p5d6 * VPhiRef1[2] + 1j * OppSy_p5d6 * VPhiRef1[1] + OppSx_p5d6 * VPhiRef1[0]
                Hamiltonian_PhiRef2 = Hamiltonian + OppSz_p5d6 * VPhiRef2[2] + 1j * OppSy_p5d6 * VPhiRef2[1] + OppSx_p5d6 * VPhiRef2[0]
            else:
                #Cristal Field D4h
                OpptenDq_p5d7 = 21 / 10 * (OppC40_p5d7 + np.sqrt(5 / 14) * (OppC4_4_p5d7 + OppC44_p5d7))
                OppDs_p5d7 = -7 * OppC20_p5d7
                OppDt_p5d7 = -21 * OppC40_p5d7

                Hamiltonian = Hamiltonian + OpptenDq_p5d7 * tenDq + OppDs_p5d7 * Ds + OppDt_p5d7 * Dt

                Hamiltonian_Phi1 = Hamiltonian + OppSz_p5d7 * VPhi1[2] + 1j * OppSy_p5d7 * VPhi1[1] + OppSx_p5d7 * VPhi1[0]
                Hamiltonian_Phi2 = Hamiltonian + OppSz_p5d7*VPhi2[2] + 1j*OppSy_p5d7*VPhi2[1] + OppSx_p5d7*VPhi2[0]
                Hamiltonian_PhiRef1 = Hamiltonian + OppSz_p5d7 * VPhiRef1[2] + 1j * OppSy_p5d7 * VPhiRef1[1] + OppSx_p5d7 * VPhiRef1[0]
                Hamiltonian_PhiRef2 = Hamiltonian + OppSz_p5d7 * VPhiRef2[2] + 1j * OppSy_p5d7 * VPhiRef2[1] + OppSx_p5d7 * VPhiRef2[0]

            # D4h
            #Sparse to full Matrix
            Hamiltonian_Phi1 = Hamiltonian_Phi1.todense()
            Hamiltonian_Phi2 = Hamiltonian_Phi2.todense()
            Hamiltonian_PhiRef1 = Hamiltonian_PhiRef1.todense()
            Hamiltonian_PhiRef2 = Hamiltonian_PhiRef2.todense()
            #diagonalization is faster with Full matrix (eigs for sparse uses different approach).
            Ef_Phi1, CCf_Phi1 = eigh(Hamiltonian_Phi1)
            Ef_Phi2,CCf_Phi2 = eigh(Hamiltonian_Phi2)
            Ef_PhiRef1, CCf_PhiRef1 = eigh(Hamiltonian_PhiRef1)
            Ef_PhiRef2, CCf_PhiRef2 = eigh(Hamiltonian_PhiRef2)

            Eff_Phi1 = Ef_Phi1 + Eshift
            Eff_Phi2 = Ef_Phi2 + Eshift
            Eff_PhiRef1 = Ef_PhiRef1 + Eshift
            Eff_PhiRef2 = Ef_PhiRef2 + Eshift

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

            Hamiltonian_Phi1_1 = Hamiltonian_site1 + OppSx_p5d7 * VPhi1[0] + 1j * OppSy_p5d7 * VPhi1[1] + OppSz_p5d7 * VPhi1[2]
            Hamiltonian_Phi2_1 = Hamiltonian_site1 + OppSx_p5d7 * VPhi2[0] + 1j * OppSy_p5d7 * VPhi2[1] + OppSz_p5d7 * VPhi2[2]
            Hamiltonian_PhiRef1_1 = Hamiltonian_site1 + OppSx_p5d7 * VPhiRef1[0] + 1j * OppSy_p5d7 * VPhiRef1[1] + OppSz_p5d7 * VPhiRef1[2]
            Hamiltonian_PhiRef2_1 = Hamiltonian_site1 + OppSx_p5d7 * VPhiRef2[0] + 1j * OppSy_p5d7 * VPhiRef2[1] + OppSz_p5d7 * VPhiRef2[2]

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

            Hamiltonian_Phi1_2 = Hamiltonian_site2 + OppSx_p5d7 * VPhi1[0] + 1j * OppSy_p5d7 * VPhi1[1] + OppSz_p5d7 * VPhi1[2]
            Hamiltonian_Phi2_2 = Hamiltonian_site2 + OppSx_p5d7*VPhi2[0] + 1j*OppSy_p5d7*VPhi2[1] + OppSz_p5d7*VPhi2[2]
            Hamiltonian_PhiRef1_2 = Hamiltonian_site2 + OppSx_p5d7 * VPhiRef1[0] + 1j * OppSy_p5d7 * VPhiRef1[1] + OppSz_p5d7 * VPhiRef1[2]
            Hamiltonian_PhiRef2_2 = Hamiltonian_site2 + OppSx_p5d7 * VPhiRef2[0] + 1j * OppSy_p5d7 * VPhiRef2[1] + OppSz_p5d7 * VPhiRef2[2]

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

            Hamiltonian_Phi1_3 = Hamiltonian_site3 + OppSx_p5d7 * VPhi1[0] + 1j * OppSy_p5d7 * VPhi1[1] + OppSz_p5d7 * VPhi1[2]
            Hamiltonian_Phi2_3 = Hamiltonian_site3 + OppSx_p5d7*VPhi2[0] + 1j*OppSy_p5d7*VPhi2[1] + OppSz_p5d7*VPhi2[2]
            Hamiltonian_PhiRef1_3 = Hamiltonian_site3 + OppSx_p5d7 * VPhiRef1[0] + 1j * OppSy_p5d7 * VPhiRef1[1] + OppSz_p5d7 * VPhiRef1[2]
            Hamiltonian_PhiRef2_3 = Hamiltonian_site3 + OppSx_p5d7 * VPhiRef2[0] + 1j * OppSy_p5d7 * VPhiRef2[1] + OppSz_p5d7 * VPhiRef2[2]

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

            Hamiltonian_Phi1_4 = Hamiltonian_site4 + OppSx_p5d7 * VPhi1[0] + 1j * OppSy_p5d7 * VPhi1[1] + OppSz_p5d7 * VPhi1[2]
            Hamiltonian_Phi2_4 = Hamiltonian_site4 + OppSx_p5d7*VPhi2[0] + 1j*OppSy_p5d7*VPhi2[1] + OppSz_p5d7*VPhi2[2]
            Hamiltonian_PhiRef1_4 = Hamiltonian_site4 + OppSx_p5d7 * VPhiRef1[0] + 1j * OppSy_p5d7 * VPhiRef1[1] + OppSz_p5d7 * VPhiRef1[2]
            Hamiltonian_PhiRef2_4 = Hamiltonian_site4 + OppSx_p5d7 * VPhiRef2[0] + 1j * OppSy_p5d7 * VPhiRef2[1] + OppSz_p5d7 * VPhiRef2[2]

            #sparse to full
            Hamiltonian_Phi1_1=Hamiltonian_Phi1_1.todense()
            Hamiltonian_Phi2_1=Hamiltonian_Phi2_1.todense()
            Hamiltonian_PhiRef1_1=Hamiltonian_PhiRef1_1.todense()
            Hamiltonian_PhiRef2_1=Hamiltonian_PhiRef2_1.todense()

            Hamiltonian_Phi1_2=Hamiltonian_Phi1_2.todense()
            Hamiltonian_Phi2_2 = Hamiltonian_Phi2_2.todense()
            Hamiltonian_PhiRef1_2 = Hamiltonian_PhiRef1_2.todense()
            Hamiltonian_PhiRef2_2=Hamiltonian_PhiRef2_2.todense()

            Hamiltonian_Phi1_3=Hamiltonian_Phi1_3.todense()
            Hamiltonian_Phi2_3 = Hamiltonian_Phi2_3.todense()
            Hamiltonian_PhiRef1_3 = Hamiltonian_PhiRef1_3.todense()
            Hamiltonian_PhiRef2_3=Hamiltonian_PhiRef2_3.todense()

            Hamiltonian_Phi1_4=Hamiltonian_Phi1_4.todense()
            Hamiltonian_Phi2_4 = Hamiltonian_Phi2_4.todense()
            Hamiltonian_PhiRef1_4 = Hamiltonian_PhiRef1_4.todense()
            Hamiltonian_PhiRef2_4=Hamiltonian_PhiRef2_4.todense()

            Ef_Phi1_1, CCf_Phi1_1 = eigh(Hamiltonian_Phi1_1)
            Ef_Phi2_1, CCf_Phi2_1 = eigh(Hamiltonian_Phi2_1)
            Ef_PhiRef1_1, CCf_PhiRef1_1 = eigh(Hamiltonian_PhiRef1_1)
            Ef_PhiRef2_1, CCf_PhiRef2_1 = eigh(Hamiltonian_PhiRef2_1)

            Ef_Phi1_2, CCf_Phi1_2 = eigh(Hamiltonian_Phi1_2)
            Ef_Phi2_2, CCf_Phi2_2 = eigh(Hamiltonian_Phi2_2)
            Ef_PhiRef1_2, CCf_PhiRef1_2 = eigh(Hamiltonian_PhiRef1_2)
            Ef_PhiRef2_2, CCf_PhiRef2_2 = eigh(Hamiltonian_PhiRef2_2)

            Ef_Phi1_3, CCf_Phi1_3 = eigh(Hamiltonian_Phi1_3)
            Ef_Phi2_3, CCf_Phi2_3 = eigh(Hamiltonian_Phi2_3)
            Ef_PhiRef1_3, CCf_PhiRef1_3 = eigh(Hamiltonian_PhiRef1_3)
            Ef_PhiRef2_3, CCf_PhiRef2_3 = eigh(Hamiltonian_PhiRef2_3)

            Ef_Phi1_4, CCf_Phi1_4 = eigh(Hamiltonian_Phi1_4)
            Ef_Phi2_4, CCf_Phi2_4 = eigh(Hamiltonian_Phi2_4)
            Ef_PhiRef1_4, CCf_PhiRef1_4 = eigh(Hamiltonian_PhiRef1_4)
            Ef_PhiRef2_4, CCf_PhiRef2_4 = eigh(Hamiltonian_PhiRef2_4)

            Eff_Phi1_1 = Ef_Phi1_1 + Eshift
            Eff_Phi2_1 = Ef_Phi2_1 + Eshift
            Eff_PhiRef1_1 = Ef_PhiRef1_1 + Eshift
            Eff_PhiRef2_1 = Ef_PhiRef2_1 + Eshift

            Eff_Phi1_2 = Ef_Phi1_2 + Eshift
            Eff_Phi2_2 = Ef_Phi2_2 + Eshift
            Eff_PhiRef1_2 = Ef_PhiRef1_2 + Eshift
            Eff_PhiRef2_2 = Ef_PhiRef2_2 + Eshift

            Eff_Phi1_3 = Ef_Phi1_3 + Eshift
            Eff_Phi2_3 = Ef_Phi2_3 + Eshift
            Eff_PhiRef1_3 = Ef_PhiRef1_3 + Eshift
            Eff_PhiRef2_3 = Ef_PhiRef2_3 + Eshift

            Eff_Phi1_4 = Ef_Phi1_4 + Eshift
            Eff_Phi2_4 = Ef_Phi2_4 + Eshift
            Eff_PhiRef1_4 = Ef_PhiRef1_4 + Eshift
            Eff_PhiRef2_4 = Ef_PhiRef2_4 + Eshift

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
                Hamiltonian_Phi1 = Hamiltonian + OppSz_d5 * VPhi1[2] + 1j * OppSy_d5 * VPhi1[1] + OppSx_d5 * VPhi1[0]
                Hamiltonian_Phi2 = Hamiltonian + OppSz_d5*VPhi2[2] + 1j*OppSy_d5*VPhi2[1] + OppSx_d5*VPhi2[0]
                Hamiltonian_PhiRef1 = Hamiltonian + OppSz_d5 * VPhiRef1[2] + 1j * OppSy_d5 * VPhiRef1[1] + OppSx_d5 * VPhiRef1[0]
                Hamiltonian_PhiRef2 = Hamiltonian + OppSz_d5 * VPhiRef2[2] + 1j * OppSy_d5 * VPhiRef2[1] + OppSx_d5 * VPhiRef2[0]
            else:
                OpptenDq = 21 / 10 * (OppC40_d6 + np.sqrt(5 / 14) * (OppC4_4_d6 + OppC44_d6))
                OppDs = -7 * OppC20_d6
                OppDt = -21 * OppC40_d6

                Hamiltonian = Hamiltonian + OpptenDq * tenDq + OppDs * Ds + OppDt * Dt

                Hamiltonian_Phi1 = Hamiltonian + OppSz_d6 * VPhi1[2] + 1j * OppSy_d6 * VPhi1[1] + OppSx_d6 * VPhi1[0]
                Hamiltonian_Phi2 = Hamiltonian + OppSz_d6*VPhi2[2] + 1j*OppSy_d6*VPhi2[1] + OppSx_d6*VPhi2[0]
                Hamiltonian_PhiRef1 = Hamiltonian + OppSz_d6 * VPhiRef1[2] + 1j * OppSy_d6 * VPhiRef1[1] + OppSx_d6 * VPhiRef1[0]
                Hamiltonian_PhiRef2 = Hamiltonian + OppSz_d6 * VPhiRef2[2] + 1j * OppSy_d6 * VPhiRef2[1] + OppSx_d6 * VPhiRef2[0]

            # Sparse to full Matrix
            Hamiltonian_Phi1 = Hamiltonian_Phi1.todense()
            Hamiltonian_Phi2 = Hamiltonian_Phi2.todense()
            Hamiltonian_PhiRef1 = Hamiltonian_PhiRef1.todense()
            Hamiltonian_PhiRef2 = Hamiltonian_PhiRef2.todense()

            Ei_Phi1, CCi_Phi1 = eigh(Hamiltonian_Phi1)
            Ei_Phi2, CCi_Phi2 = eigh(Hamiltonian_Phi2)
            Ei_PhiRef1, CCi_PhiRef1 = eigh(Hamiltonian_PhiRef1)
            Ei_PhiRef2, CCi_PhiRef2 = eigh(Hamiltonian_PhiRef2)

            fB_Phi1 = np.zeros(Ei_Phi1.size, dtype='complex_')
            fB_Phi2 = np.zeros(Ei_Phi2.size, dtype='complex_')
            fB_PhiRef1 = np.zeros(Ei_PhiRef1.size, dtype='complex_')
            fB_PhiRef2 = np.zeros(Ei_PhiRef2.size, dtype='complex_')

            l = 0  # Variable auxialiar para recorrer fB y Ei

            for i in Ei_Phi1:
                fB_Phi1[l] = np.exp(-(i - Ei_Phi1[0]) / (8.61733035e-5 * Temp))
                fB_Phi2[l] = np.exp(-(Ei_Phi2[l] - Ei_Phi2[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef1[l] = np.exp(-(Ei_PhiRef1[l] - Ei_PhiRef1[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef2[l] = np.exp(-(Ei_PhiRef2[l] - Ei_PhiRef2[0]) / (8.61733035e-5 * Temp))
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

            Hamiltonian_Phi1_1 = Hamiltonian_site1 + OppSx_d6 * VPhi1[0] + 1j * OppSy_d6 * VPhi1[1] + OppSz_d6 * VPhi1[2]
            Hamiltonian_Phi2_1 = Hamiltonian_site1 + OppSx_d6*VPhi2[0] + 1j*OppSy_d6*VPhi2[1] + OppSz_d6*VPhi2[2]
            Hamiltonian_PhiRef1_1 = Hamiltonian_site1 + OppSx_d6 * VPhiRef1[0] + 1j * OppSy_d6 * VPhiRef1[1] + OppSz_d6 * VPhiRef1[2]
            Hamiltonian_PhiRef2_1 = Hamiltonian_site1 + OppSx_d6 * VPhiRef2[0] + 1j * OppSy_d6 * VPhiRef2[1] + OppSz_d6 * VPhiRef2[2]

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

            Hamiltonian_Phi1_2 = Hamiltonian_site2 + OppSx_d6 * VPhi1[0] + 1j * OppSy_d6 * VPhi1[1] + OppSz_d6 * VPhi1[2]
            Hamiltonian_Phi2_2 = Hamiltonian_site2 + OppSx_d6*VPhi2[0] + 1j*OppSy_d6*VPhi2[1] + OppSz_d6*VPhi2[2]
            Hamiltonian_PhiRef1_2 = Hamiltonian_site2 + OppSx_d6 * VPhiRef1[0] + 1j * OppSy_d6 * VPhiRef1[1] + OppSz_d6 * VPhiRef1[2]
            Hamiltonian_PhiRef2_2 = Hamiltonian_site2 + OppSx_d6 * VPhiRef2[0] + 1j * OppSy_d6 * VPhiRef2[1] + OppSz_d6 * VPhiRef2[2]

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

            Hamiltonian_Phi1_3 = Hamiltonian_site3 + OppSx_d6 * VPhi1[0] + 1j * OppSy_d6 * VPhi1[1] + OppSz_d6 * VPhi1[2]
            Hamiltonian_Phi2_3 = Hamiltonian_site3 + OppSx_d6*VPhi2[0] + 1j*OppSy_d6*VPhi2[1] + OppSz_d6*VPhi2[2]
            Hamiltonian_PhiRef1_3 = Hamiltonian_site3 + OppSx_d6 * VPhiRef1[0] + 1j * OppSy_d6 * VPhiRef1[1] + OppSz_d6 * VPhiRef1[2]
            Hamiltonian_PhiRef2_3 = Hamiltonian_site3 + OppSx_d6 * VPhiRef2[0] + 1j * OppSy_d6 * VPhiRef2[1] + OppSz_d6 * VPhiRef2[2]

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

            Hamiltonian_Phi1_4 = Hamiltonian_site4 + OppSx_d6 * VPhi1[0] + 1j * OppSy_d6 * VPhi1[1] + OppSz_d6 * VPhi1[2]
            Hamiltonian_Phi2_4 = Hamiltonian_site4 + OppSx_d6*VPhi2[0] + 1j*OppSy_d6*VPhi2[1] + OppSz_d6*VPhi2[2]
            Hamiltonian_PhiRef1_4 = Hamiltonian_site4 + OppSx_d6 * VPhiRef1[0] + 1j * OppSy_d6 * VPhiRef1[1] + OppSz_d6 * VPhiRef1[2]
            Hamiltonian_PhiRef2_4 = Hamiltonian_site4 + OppSx_d6 * VPhiRef2[0] + 1j * OppSy_d6 * VPhiRef2[1] + OppSz_d6 * VPhiRef2[2]

            # sparse to full
            Hamiltonian_Phi1_1 = Hamiltonian_Phi1_1.todense()
            Hamiltonian_Phi2_1 = Hamiltonian_Phi2_1.todense()
            Hamiltonian_PhiRef1_1 = Hamiltonian_PhiRef1_1.todense()
            Hamiltonian_PhiRef2_1 = Hamiltonian_PhiRef2_1.todense()

            Hamiltonian_Phi1_2 = Hamiltonian_Phi1_2.todense()
            Hamiltonian_Phi2_2 = Hamiltonian_Phi2_2.todense()
            Hamiltonian_PhiRef1_2 = Hamiltonian_PhiRef1_2.todense()
            Hamiltonian_PhiRef2_2 = Hamiltonian_PhiRef2_2.todense()

            Hamiltonian_Phi1_3 = Hamiltonian_Phi1_3.todense()
            Hamiltonian_Phi2_3 = Hamiltonian_Phi2_3.todense()
            Hamiltonian_PhiRef1_3 = Hamiltonian_PhiRef1_3.todense()
            Hamiltonian_PhiRef2_3 = Hamiltonian_PhiRef2_3.todense()

            Hamiltonian_Phi1_4 = Hamiltonian_Phi1_4.todense()
            Hamiltonian_Phi2_4 = Hamiltonian_Phi2_4.todense()
            Hamiltonian_PhiRef1_4 = Hamiltonian_PhiRef1_4.todense()
            Hamiltonian_PhiRef2_4 = Hamiltonian_PhiRef2_4.todense()

            Ei_Phi1_1, CCi_Phi1_1 = eigh(Hamiltonian_Phi1_1)
            Ei_Phi2_1, CCi_Phi2_1 = eigh(Hamiltonian_Phi2_1)
            Ei_PhiRef1_1, CCi_PhiRef1_1 = eigh(Hamiltonian_PhiRef1_1)
            Ei_PhiRef2_1, CCi_PhiRef2_1 = eigh(Hamiltonian_PhiRef2_1)

            Ei_Phi1_2, CCi_Phi1_2 = eigh(Hamiltonian_Phi1_2)
            Ei_Phi2_2, CCi_Phi2_2 = eigh(Hamiltonian_Phi2_2)
            Ei_PhiRef1_2, CCi_PhiRef1_2 = eigh(Hamiltonian_PhiRef1_2)
            Ei_PhiRef2_2, CCi_PhiRef2_2 = eigh(Hamiltonian_PhiRef2_2)

            Ei_Phi1_3, CCi_Phi1_3 = eigh(Hamiltonian_Phi1_3)
            Ei_Phi2_3, CCi_Phi2_3 = eigh(Hamiltonian_Phi2_3)
            Ei_PhiRef1_3, CCi_PhiRef1_3 = eigh(Hamiltonian_PhiRef1_3)
            Ei_PhiRef2_3, CCi_PhiRef2_3 = eigh(Hamiltonian_PhiRef2_3)

            Ei_Phi1_4, CCi_Phi1_4 = eigh(Hamiltonian_Phi1_4)
            Ei_Phi2_4, CCi_Phi2_4 = eigh(Hamiltonian_Phi2_4)
            Ei_PhiRef1_4, CCi_PhiRef1_4 = eigh(Hamiltonian_PhiRef1_4)
            Ei_PhiRef2_4, CCi_PhiRef2_4 = eigh(Hamiltonian_PhiRef2_4)

            fB_Phi1_1 = np.zeros(Ei_Phi1_1.size, dtype='complex_')
            fB_Phi2_1 = np.zeros(Ei_Phi2_1.size,dtype='complex_')
            fB_PhiRef1_1 = np.zeros(Ei_PhiRef1_1.size, dtype='complex_')
            fB_PhiRef2_1 = np.zeros(Ei_PhiRef2_1.size, dtype='complex_')

            fB_Phi1_2 = np.zeros(Ei_Phi1_2.size, dtype='complex_')
            fB_Phi2_2 = np.zeros(Ei_Phi2_2.size, dtype='complex_')
            fB_PhiRef1_2 = np.zeros(Ei_PhiRef1_2.size, dtype='complex_')
            fB_PhiRef2_2 = np.zeros(Ei_PhiRef2_2.size, dtype='complex_')

            fB_Phi1_3 = np.zeros(Ei_Phi1_3.size, dtype='complex_')
            fB_Phi2_3 = np.zeros(Ei_Phi2_3.size, dtype='complex_')
            fB_PhiRef1_3 = np.zeros(Ei_PhiRef1_3.size, dtype='complex_')
            fB_PhiRef2_3 = np.zeros(Ei_PhiRef2_3.size, dtype='complex_')

            fB_Phi1_4 = np.zeros(Ei_Phi1_4.size, dtype='complex_')
            fB_Phi2_4 = np.zeros(Ei_Phi2_4.size, dtype='complex_')
            fB_PhiRef1_4 = np.zeros(Ei_PhiRef1_4.size, dtype='complex_')
            fB_PhiRef2_4 = np.zeros(Ei_PhiRef2_4.size, dtype='complex_')

            #Boltzman
            l = 0  # Variable auxialiar para recorrer fB y Ei

            for i in Ei_Phi1_1:
                fB_Phi1_1[l] = np.exp(-(i - Ei_Phi1_1[0]) / (8.61733035e-5 * Temp))
                fB_Phi2_1[l] = np.exp(-(Ei_Phi2_1[l] - Ei_Phi2_1[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef1_1[l] = np.exp(-(Ei_PhiRef1_1[l] - Ei_PhiRef1_1[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef2_1[l] = np.exp(-(Ei_PhiRef2_1[l] - Ei_PhiRef2_1[0]) / (8.61733035e-5 * Temp))

                fB_Phi1_2[l] = np.exp(-(Ei_Phi1_2[l] - Ei_Phi1_2[0]) / (8.61733035e-5 * Temp))
                fB_Phi2_2[l] = np.exp(-(Ei_Phi2_2[l] - Ei_Phi2_2[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef1_2[l] = np.exp(-(Ei_PhiRef1_2[l] - Ei_PhiRef1_2[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef2_2[l] = np.exp(-(Ei_PhiRef2_2[l] - Ei_PhiRef2_2[0]) / (8.61733035e-5 * Temp))

                fB_Phi1_3[l] = np.exp(-(Ei_Phi1_3[l] - Ei_Phi1_3[0]) / (8.61733035e-5 * Temp))
                fB_Phi2_3[l] = np.exp(-(Ei_Phi2_3[l] - Ei_Phi2_3[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef1_3[l] = np.exp(-(Ei_PhiRef1_3[l] - Ei_PhiRef1_3[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef2_3[l] = np.exp(-(Ei_PhiRef2_3[l] - Ei_PhiRef2_3[0]) / (8.61733035e-5 * Temp))

                fB_Phi1_4[l] = np.exp(-(Ei_Phi1_4[l] - Ei_Phi1_4[0]) / (8.61733035e-5 * Temp))
                fB_Phi2_4[l] = np.exp(-(Ei_Phi2_4[l] - Ei_Phi2_4[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef1_4[l] = np.exp(-(Ei_PhiRef1_4[l] - Ei_PhiRef1_4[0]) / (8.61733035e-5 * Temp))
                fB_PhiRef2_4[l] = np.exp(-(Ei_PhiRef2_4[l] - Ei_PhiRef2_4[0]) / (8.61733035e-5 * Temp))

                l += 1

        #Transitions -Fermi Golden Rule
        #Polarizations
        if Nelec==5:
            opp_aux=(np.cos(np.radians(Pol))/np.sqrt(2)) * (OppT1_1_d5-OppT11_d5) + (np.sin(np.radians(Pol))*1j/np.sqrt(2)) * (OppT11_d5+OppT1_1_d5)
            OppTH = np.array(opp_aux.todense()) #L30
            OppTr = np.array(OppT1_1_d5.todense()) #Circular a la derecha
            OppTl = np.array(OppT11_d5.todense()) #Circular a la izquierda
        else:
            opp_aux = (np.cos(np.radians(Pol))/np.sqrt(2))*(OppT1_1_d6-OppT11_d6) + (np.sin(np.radians(Pol))*1j/np.sqrt(2))*(OppT11_d6+OppT1_1_d6)
            OppTH = np.array(opp_aux.todense())  # L30
            OppTr = np.array(OppT1_1_d6.todense())  # Circular a la derecha
            OppTl = np.array(OppT11_d6.todense())  # Circular a la izquierda

        if Sym == 1: #D4h
            fHPhi1_ = np.array([])
            EnHPhi1_ = np.array([])

            fHPhi2_ = np.array([])
            EnHPhi2_ = np.array([])

            fHPhiRef1_ = np.array([])
            EnHPhiRef1_ = np.array([])

            fHPhiRef2_ = np.array([])
            EnHPhiRef2_ = np.array([])

            fr_ = np.array([])
            Enr_ = np.array([])

            fl_ = np.array([])
            Enl_ = np.array([])

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei_Phi1:
                if fB_Phi1[k] < 0.01:
                    break

                Smed_HPhi1 = abs(CCi_Phi1[:, k].conj().transpose() @ OppTH @ CCf_Phi1) ** 2 * fB_Phi1[k]
                Smed_r = abs(CCi_Phi1[:, k].conj().transpose() @ OppTr @ CCf_Phi1) ** 2 * fB_Phi1[k]
                Smed_l = abs(CCi_Phi1[:, k].conj().transpose() @ OppTl @ CCf_Phi1) ** 2 * fB_Phi1[k]

                fr_ = np.append(fr_, Smed_r[Smed_r > 1e-5].conj().transpose())
                Enr_ = np.append(Enr_, np.around((Eff_Phi1[Smed_r > 1e-5] - i) * 100) / 100)

                fl_ = np.append(fl_, Smed_l[Smed_l > 1e-5].conj().transpose())
                Enl_ = np.append(Enl_, np.around((Eff_Phi1[Smed_l > 1e-5] - i) * 100) / 100)

                fHPhi1_ = np.append(fHPhi1_, Smed_HPhi1[Smed_HPhi1 > 1e-5].conj().transpose())
                EnHPhi1_ = np.append(EnHPhi1_, np.around((Eff_Phi1[Smed_HPhi1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhi1_, return_index=True)
            EnHPhi1 = EnHPhi1_[np.sort(b)]
            fHPhi1 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnHPhi1:
                for j in EnHPhi1_:
                    if i == j:
                        aux = aux + fHPhi1_[aux1]
                    aux1 += 1
                fHPhi1 = np.append(fHPhi1, aux)
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

            for i in Ei_Phi2:
                if fB_Phi2[k] < 0.01:
                    break
                Smed_HPhi2 = abs(CCi_Phi2[:, k].conj().transpose() @ OppTH @ CCf_Phi2) ** 2 * fB_Phi2[k]
                fHPhi2_ = np.append(fHPhi2_, Smed_HPhi2[Smed_HPhi2 > 1e-5].conj().transpose())
                EnHPhi2_ = np.append(EnHPhi2_, np.around((Eff_Phi2[Smed_HPhi2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhi2_, return_index=True)
            EnHPhi2 = EnHPhi2_[np.sort(b)]
            fHPhi2 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnHPhi2:
                for j in EnHPhi2_:
                    if i == j:
                        aux = aux + fHPhi2_[aux1]
                    aux1 += 1
                fHPhi2 = np.append(fHPhi2, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei_PhiRef1:
                if fB_PhiRef1[k] < 0.01:
                    break

                Smed_HPhiRef1 = abs(CCi_PhiRef1[:, k].conj().transpose() @ OppTH @ CCf_PhiRef1) ** 2 * fB_PhiRef1[k]
                fHPhiRef1_ = np.append(fHPhiRef1_, Smed_HPhiRef1[Smed_HPhiRef1 > 1e-5].conj().transpose())
                EnHPhiRef1_ = np.append(EnHPhiRef1_, np.around((Eff_PhiRef1[Smed_HPhiRef1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhiRef1_, return_index=True)
            EnHPhiRef1 = EnHPhiRef1_[np.sort(b)]
            fHPhiRef1 = np.array([])
            aux = 0
            aux1 = 0

            for i in EnHPhiRef1:
                for j in EnHPhiRef1_:
                    if i == j:
                        aux = aux + fHPhiRef1_[aux1]
                    aux1 += 1
                fHPhiRef1 = np.append(fHPhiRef1, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop

            for i in Ei_PhiRef2:
                if fB_PhiRef2[k] < 0.01:
                    break

                Smed_HPhiRef2 = abs(CCi_PhiRef2[:, k].conj().transpose() @ OppTH @ CCf_PhiRef2) ** 2 * fB_PhiRef2[k]
                fHPhiRef2_ = np.append(fHPhiRef2_, Smed_HPhiRef2[Smed_HPhiRef2 > 1e-5].conj().transpose())
                EnHPhiRef2_ = np.append(EnHPhiRef2_, np.around((Eff_PhiRef2[Smed_HPhiRef2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhiRef2_, return_index=True)
            EnHPhiRef2 = EnHPhiRef2_[np.sort(b)]
            fHPhiRef2 = np.array([])
            aux = 0
            aux1 = 0

            for i in EnHPhiRef2:
                for j in EnHPhiRef2_:
                    if i == j:
                        aux = aux + fHPhiRef2_[aux1]
                    aux1 += 1
                fHPhiRef2 = np.append(fHPhiRef2, aux)
                aux = 0
                aux1 = 0
        else:
            #D3d
            fHPhi1_1_ = np.array([])
            EnHPhi1_1_ = np.array([])

            fH_Phi1_2_ = np.array([])
            EnH_Phi1_2_ = np.array([])

            fHPhi1_3_ = np.array([])
            EnHPhi1_3_ = np.array([])

            fH_Phi1_4_ = np.array([])
            EnH_Phi1_4_ = np.array([])

            fH_Phi2_1_ = np.array([])
            EnH_Phi2_1_ = np.array([])

            fH_Phi2_2_ = np.array([])
            EnH_Phi2_2_ = np.array([])

            fH_Phi2_3_ = np.array([])
            EnH_Phi2_3_ = np.array([])

            fH_Phi2_4_ = np.array([])
            EnH_Phi2_4_ = np.array([])

            fH_PhiRef1_1_ = np.array([])
            EnH_PhiRef1_1_ = np.array([])

            fH_PhiRef1_2_ = np.array([])
            EnH_PhiRef1_2_ = np.array([])

            fH_PhiRef1_3_ = np.array([])
            EnH_PhiRef1_3_ = np.array([])

            fH_PhiRef1_4_ = np.array([])
            EnH_PhiRef1_4_ = np.array([])

            fH_PhiRef2_1_ = np.array([])
            EnH_PhiRef2_1_ = np.array([])

            fH_PhiRef2_2_ = np.array([])
            EnH_PhiRef2_2_ = np.array([])

            fH_PhiRef2_3_ = np.array([])
            EnH_PhiRef2_3_ = np.array([])

            fH_PhiRef2_4_ = np.array([])
            EnH_PhiRef2_4_ = np.array([])

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
            for i in Ei_Phi1_1:
                if fB_Phi1_1[k] < 0.01:
                    break

                Smed_HPhi1_1 = abs(CCi_Phi1_1[:, k].conj().transpose() @ OppTH @ CCf_Phi1_1) ** 2 * fB_Phi1_1[k]
                Smed_r1 = abs(CCi_Phi1_1[:, k].conj().transpose() @ OppTr @ CCf_Phi1_1) ** 2 * fB_Phi1_1[k]
                Smed_l1 = abs(CCi_Phi1_1[:, k].conj().transpose() @ OppTl @ CCf_Phi1_1) ** 2 * fB_Phi1_1[k]

                fHPhi1_1_ = np.append(fHPhi1_1_, Smed_HPhi1_1[Smed_HPhi1_1 > 1e-5].conj().transpose())
                EnHPhi1_1_ = np.append(EnHPhi1_1_, np.around((Eff_Phi1_1[Smed_HPhi1_1 > 1e-5] - i) * 100)/100)

                fr1_ = np.append(fr1_, Smed_r1[Smed_r1 > 1e-5].conj().transpose())
                Enr1_ = np.append(Enr1_, np.around((Eff_Phi1_1[Smed_r1 > 1e-5] - i) * 100) / 100)

                fl1_ = np.append(fl1_, Smed_l1[Smed_l1 > 1e-5].conj().transpose())
                Enl1_ = np.append(Enl1_, np.around((Eff_Phi1_1[Smed_l1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhi1_1_, return_index=True)
            EnHPhi1_1 = EnHPhi1_1_[np.sort(b)]
            fHPhi1_1 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnHPhi1_1:
                for j in EnHPhi1_1_:
                    if i == j:
                        aux = aux + fHPhi1_1_[aux1]
                    aux1 += 1
                fHPhi1_1 = np.append(fHPhi1_1, aux)
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
            for i in Ei_Phi2_1:
                if fB_Phi2_1[k] < 0.01:
                    break

                Smed_H_Phi2_1 = abs(CCi_Phi2_1[:, k].conj().transpose() @ OppTH @ CCf_Phi2_1) ** 2 * fB_Phi2_1[k]
                fH_Phi2_1_ = np.append(fH_Phi2_1_, Smed_H_Phi2_1[Smed_H_Phi2_1 > 1e-5].conj().transpose())
                EnH_Phi2_1_ = np.append(EnH_Phi2_1_, np.around((Eff_Phi2_1[Smed_H_Phi2_1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi2_1_, return_index=True)
            EnH_Phi2_1 = EnH_Phi2_1_[np.sort(b)]
            fH_Phi2_1 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi2_1:
                for j in EnH_Phi2_1_:
                    if i == j:
                        aux = aux + fH_Phi2_1_[aux1]
                    aux1 += 1
                fH_Phi2_1 = np.append(fH_Phi2_1, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef1_1:
                if fB_PhiRef1_1[k] < 0.01:
                    break

                Smed_H_PhiRef1_1 = abs(CCi_PhiRef1_1[:, k].conj().transpose() @ OppTH @ CCf_PhiRef1_1) ** 2 * fB_PhiRef1_1[k]
                fH_PhiRef1_1_ = np.append(fH_PhiRef1_1_, Smed_H_PhiRef1_1[Smed_H_PhiRef1_1 > 1e-5].conj().transpose())
                EnH_PhiRef1_1_ = np.append(EnH_PhiRef1_1_, np.around((Eff_PhiRef1_1[Smed_H_PhiRef1_1 > 1e-5] - i) * 100)/100)

                k += 1

            a, b = np.unique(EnH_PhiRef1_1_, return_index=True)
            EnH_PhiRef1_1 = EnH_PhiRef1_1_[np.sort(b)]
            fH_PhiRef1_1 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef1_1:
                for j in EnH_PhiRef1_1_:
                    if i == j:
                        aux = aux + fH_PhiRef1_1_[aux1]
                    aux1 += 1
                fH_PhiRef1_1 = np.append(fH_PhiRef1_1, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef2_1:
                if fB_PhiRef2_1[k] < 0.01:
                    break

                Smed_H_PhiRef2_1 = abs(CCi_PhiRef2_1[:, k].conj().transpose() @ OppTH @ CCf_PhiRef2_1) ** 2 * fB_PhiRef2_1[k]
                fH_PhiRef2_1_ = np.append(fH_PhiRef2_1_, Smed_H_PhiRef2_1[Smed_H_PhiRef2_1 > 1e-5].conj().transpose())
                EnH_PhiRef2_1_ = np.append(EnH_PhiRef2_1_, np.around((Eff_PhiRef2_1[Smed_H_PhiRef2_1 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef2_1_, return_index=True)
            EnH_PhiRef2_1 = EnH_PhiRef2_1_[np.sort(b)]
            fH_PhiRef2_1 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef2_1:
                for j in EnH_PhiRef2_1_:
                    if i == j:
                        aux = aux + fH_PhiRef2_1_[aux1]
                    aux1 += 1
                fH_PhiRef2_1 = np.append(fH_PhiRef2_1, aux)
                aux = 0
                aux1 = 0

            # Sitio 2

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_Phi1_2:
                if fB_Phi1_2[k] < 0.01:
                    break

                Smed_H_Phi1_2 = abs(CCi_Phi1_2[:, k].conj().transpose() @ OppTH @ CCf_Phi1_2) ** 2 * fB_Phi1_2[k]
                Smed_r2 = abs(CCi_Phi1_2[:, k].conj().transpose() @ OppTr @ CCf_Phi1_2) ** 2 * fB_Phi1_2[k]
                Smed_l2 = abs(CCi_Phi1_2[:, k].conj().transpose() @ OppTl @ CCf_Phi1_2) ** 2 * fB_Phi1_2[k]

                fH_Phi1_2_ = np.append(fH_Phi1_2_, Smed_H_Phi1_2[Smed_H_Phi1_2 > 1e-5].conj().transpose())
                EnH_Phi1_2_ = np.append(EnH_Phi1_2_, np.around((Eff_Phi1_2[Smed_H_Phi1_2 > 1e-5] - i) * 100) / 100)

                fr2_ = np.append(fr2_, Smed_r2[Smed_r2 > 1e-5].conj().transpose())
                Enr2_ = np.append(Enr2_, np.around((Eff_Phi1_2[Smed_r2 > 1e-5] - i) * 100) / 100)

                fl2_ = np.append(fl2_, Smed_l2[Smed_l2 > 1e-5].conj().transpose())
                Enl2_ = np.append(Enl2_, np.around((Eff_Phi1_2[Smed_l2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi1_2_, return_index=True)
            EnH_Phi1_2 = EnH_Phi1_2_[np.sort(b)]
            fH_Phi1_2 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi1_2:
                for j in EnH_Phi1_2_:
                    if i == j:
                        aux = aux + fH_Phi1_2_[aux1]
                    aux1 += 1
                fH_Phi1_2 = np.append(fH_Phi1_2, aux)
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
            for i in Ei_Phi2_2:
                if fB_Phi2_2[k] < 0.01:
                    break

                Smed_H_Phi2_2 = abs(CCi_Phi2_2[:, k].conj().transpose() @ OppTH @ CCf_Phi2_2) ** 2 * fB_Phi2_2[k]
                fH_Phi2_2_ = np.append(fH_Phi2_2_, Smed_H_Phi2_2[Smed_H_Phi2_2 > 1e-5].conj().transpose())
                EnH_Phi2_2_ = np.append(EnH_Phi2_2_, np.around((Eff_Phi2_2[Smed_H_Phi2_2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi2_2_, return_index=True)
            EnH_Phi2_2 = EnH_Phi2_2_[np.sort(b)]
            fH_Phi2_2 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi2_2:
                for j in EnH_Phi2_2_:
                    if i == j:
                        aux = aux + fH_Phi2_2_[aux1]
                    aux1 += 1
                fH_Phi2_2 = np.append(fH_Phi2_2, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef1_2:
                if fB_PhiRef1_2[k] < 0.01:
                    break

                Smed_H_PhiRef1_2 = abs(CCi_PhiRef1_2[:, k].conj().transpose() @ OppTH @ CCf_PhiRef1_2) ** 2 * fB_PhiRef1_2[k]
                fH_PhiRef1_2_ = np.append(fH_PhiRef1_2_, Smed_H_PhiRef1_2[Smed_H_PhiRef1_2 > 1e-5].conj().transpose())
                EnH_PhiRef1_2_ = np.append(EnH_PhiRef1_2_, np.around((Eff_PhiRef1_2[Smed_H_PhiRef1_2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef1_2_, return_index=True)
            EnH_PhiRef1_2 = EnH_PhiRef1_2_[np.sort(b)]
            fH_PhiRef1_2 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef1_2:
                for j in EnH_PhiRef1_2_:
                    if i == j:
                        aux = aux + fH_PhiRef1_2_[aux1]
                    aux1 += 1
                fH_PhiRef1_2 = np.append(fH_PhiRef1_2, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef2_2:
                if fB_PhiRef2_2[k] < 0.01:
                    break

                Smed_H_PhiRef2_2 = abs(CCi_PhiRef2_2[:, k].conj().transpose() @ OppTH @ CCf_PhiRef2_2) ** 2 * fB_PhiRef2_2[k]
                fH_PhiRef2_2_ = np.append(fH_PhiRef2_2_, Smed_H_PhiRef2_2[Smed_H_PhiRef2_2 > 1e-5].conj().transpose())
                EnH_PhiRef2_2_ = np.append(EnH_PhiRef2_2_, np.around((Eff_PhiRef2_2[Smed_H_PhiRef2_2 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef2_2_, return_index=True)
            EnH_PhiRef2_2 = EnH_PhiRef2_2_[np.sort(b)]
            fH_PhiRef2_2 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef2_2:
                for j in EnH_PhiRef2_2_:
                    if i == j:
                        aux = aux + fH_PhiRef2_2_[aux1]
                    aux1 += 1
                fH_PhiRef2_2 = np.append(fH_PhiRef2_2, aux)
                aux = 0
                aux1 = 0

            # Sitio 3

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_Phi1_3:
                if fB_Phi1_3[k] < 0.01:
                    break

                Smed_Phi1_3 = abs(CCi_Phi1_3[:, k].conj().transpose() @ OppTH @ CCf_Phi1_3) ** 2 * fB_Phi1_3[k]
                Smed_r3 = abs(CCi_Phi1_3[:, k].conj().transpose() @ OppTr @ CCf_Phi1_3) ** 2 * fB_Phi1_3[k]
                Smed_l3 = abs(CCi_Phi1_3[:, k].conj().transpose() @ OppTl @ CCf_Phi1_3) ** 2 * fB_Phi1_3[k]

                fHPhi1_3_ = np.append(fHPhi1_3_, Smed_Phi1_3[Smed_Phi1_3 > 1e-5].conj().transpose())
                EnHPhi1_3_ = np.append(EnHPhi1_3_, np.around((Eff_Phi1_3[Smed_Phi1_3 > 1e-5] - i) * 100) / 100)

                fr3_ = np.append(fr3_, Smed_r3[Smed_r3 > 1e-5].conj().transpose())
                Enr3_ = np.append(Enr3_, np.around((Eff_Phi1_3[Smed_r3 > 1e-5] - i) * 100) / 100)

                fl3_ = np.append(fl3_, Smed_l3[Smed_l3 > 1e-5].conj().transpose())
                Enl3_ = np.append(Enl3_, np.around((Eff_Phi1_3[Smed_l3 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnHPhi1_3_, return_index=True)
            EnHPhi1_3 = EnHPhi1_3_[np.sort(b)]
            fHPhi1_3 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnHPhi1_3:
                for j in EnHPhi1_3_:
                    if i == j:
                        aux = aux + fHPhi1_3_[aux1]
                    aux1 += 1
                fHPhi1_3 = np.append(fHPhi1_3, aux)
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
            for i in Ei_Phi2_3:
                if fB_Phi2_3[k] < 0.01:
                    break

                Smed_H_Phi2_3 = abs(CCi_Phi2_3[:, k].conj().transpose() @ OppTH @ CCf_Phi2_3) ** 2 * fB_Phi2_3[k]
                fH_Phi2_3_ = np.append(fH_Phi2_3_, Smed_H_Phi2_3[Smed_H_Phi2_3 > 1e-5].conj().transpose())
                EnH_Phi2_3_ = np.append(EnH_Phi2_3_, np.around((Eff_Phi2_3[Smed_H_Phi2_3 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi2_3_, return_index=True)
            EnH_Phi2_3 = EnH_Phi2_3_[np.sort(b)]
            fH_Phi2_3 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi2_3:
                for j in EnH_Phi2_3_:
                    if i == j:
                        aux = aux + fH_Phi2_3_[aux1]
                    aux1 += 1
                fH_Phi2_3 = np.append(fH_Phi2_3, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef1_3:
                if fB_PhiRef1_3[k] < 0.01:
                    break

                Smed_H_PhiRef1_3 = abs(CCi_PhiRef1_3[:, k].conj().transpose() @ OppTH @ CCf_PhiRef1_3) ** 2 * fB_PhiRef1_3[k]
                fH_PhiRef1_3_ = np.append(fH_PhiRef1_3_, Smed_H_PhiRef1_3[Smed_H_PhiRef1_3 > 1e-5].conj().transpose())
                EnH_PhiRef1_3_ = np.append(EnH_PhiRef1_3_, np.around((Eff_PhiRef1_3[Smed_H_PhiRef1_3 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef1_3_, return_index=True)
            EnH_PhiRef1_3 = EnH_PhiRef1_3_[np.sort(b)]
            fH_PhiRef1_3 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef1_3:
                for j in EnH_PhiRef1_3_:
                    if i == j:
                        aux = aux + fH_PhiRef1_3_[aux1]
                    aux1 += 1
                fH_PhiRef1_3 = np.append(fH_PhiRef1_3, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef2_3:
                if fB_PhiRef2_3[k] < 0.01:
                    break

                Smed_H_PhiRef2_3 = abs(CCi_PhiRef2_3[:, k].conj().transpose() @ OppTH @ CCf_PhiRef2_3) ** 2 * fB_PhiRef2_3[k]
                fH_PhiRef2_3_ = np.append(fH_PhiRef2_3_, Smed_H_PhiRef2_3[Smed_H_PhiRef2_3 > 1e-5].conj().transpose())
                EnH_PhiRef2_3_ = np.append(EnH_PhiRef2_3_, np.around((Eff_PhiRef2_3[Smed_H_PhiRef2_3 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef2_3_, return_index=True)
            EnH_PhiRef2_3 = EnH_PhiRef2_3_[np.sort(b)]
            fH_PhiRef2_3 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef2_3:
                for j in EnH_PhiRef2_3_:
                    if i == j:
                        aux = aux + fH_PhiRef2_3_[aux1]
                    aux1 += 1
                fH_PhiRef2_3 = np.append(fH_PhiRef2_3, aux)
                aux = 0
                aux1 = 0

            # Sitio 4

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_Phi1_4:
                if fB_Phi1_4[k] < 0.01:
                    break

                Smed_H_Phi1_4 = abs(CCi_Phi1_4[:, k].conj().transpose() @ OppTH @ CCf_Phi1_4) ** 2 * fB_Phi1_4[k]
                Smed_r4 = abs(CCi_Phi1_4[:, k].conj().transpose() @ OppTr @ CCf_Phi1_4) ** 2 * fB_Phi1_4[k]
                Smed_l4 = abs(CCi_Phi1_4[:, k].conj().transpose() @ OppTl @ CCf_Phi1_4) ** 2 * fB_Phi1_4[k]

                fH_Phi1_4_ = np.append(fH_Phi1_4_, Smed_H_Phi1_4[Smed_H_Phi1_4 > 1e-5].conj().transpose())
                EnH_Phi1_4_ = np.append(EnH_Phi1_4_, np.around((Eff_Phi1_4[Smed_H_Phi1_4 > 1e-5] - i) * 100) / 100)

                fr4_ = np.append(fr4_, Smed_r4[Smed_r4 > 1e-5].conj().transpose())
                Enr4_ = np.append(Enr4_, np.around((Eff_Phi1_4[Smed_r4 > 1e-5] - i) * 100) / 100)

                fl4_ = np.append(fl4_, Smed_l4[Smed_l4 > 1e-5].conj().transpose())
                Enl4_ = np.append(Enl4_, np.around((Eff_Phi1_4[Smed_l4 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi1_4_, return_index=True)
            EnH_Phi1_4 = EnH_Phi1_4_[np.sort(b)]
            fH_Phi1_4 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi1_4:
                for j in EnH_Phi1_4_:
                    if i == j:
                        aux = aux + fH_Phi1_4_[aux1]
                    aux1 += 1
                fH_Phi1_4 = np.append(fH_Phi1_4, aux)
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
            for i in Ei_Phi2_4:
                if fB_Phi2_4[k] < 0.01:
                    break

                Smed_H_Phi2_4 = abs(CCi_Phi2_4[:, k].conj().transpose() @ OppTH @ CCf_Phi2_4) ** 2 * fB_Phi2_4[k]
                fH_Phi2_4_ = np.append(fH_Phi2_4_, Smed_H_Phi2_4[Smed_H_Phi2_4 > 1e-5].conj().transpose())
                EnH_Phi2_4_ = np.append(EnH_Phi2_4_, np.around((Eff_Phi2_4[Smed_H_Phi2_4 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_Phi2_4_, return_index=True)
            EnH_Phi2_4 = EnH_Phi2_4_[np.sort(b)]
            fH_Phi2_4 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_Phi2_4:
                for j in EnH_Phi2_4_:
                    if i == j:
                        aux = aux + fH_Phi2_4_[aux1]
                    aux1 += 1
                fH_Phi2_4 = np.append(fH_Phi2_4, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef1_4:
                if fB_PhiRef1_4[k] < 0.01:
                    break

                Smed_H_PhiRef1_4 = abs(CCi_PhiRef1_4[:, k].conj().transpose() @ OppTH @ CCf_PhiRef1_4) ** 2 * fB_PhiRef1_4[k]
                fH_PhiRef1_4_ = np.append(fH_PhiRef1_4_, Smed_H_PhiRef1_4[Smed_H_PhiRef1_4 > 1e-5].conj().transpose())
                EnH_PhiRef1_4_ = np.append(EnH_PhiRef1_4_, np.around((Eff_PhiRef1_4[Smed_H_PhiRef1_4 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef1_4_, return_index=True)
            EnH_PhiRef1_4 = EnH_PhiRef1_4_[np.sort(b)]
            fH_PhiRef1_4 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef1_4:
                for j in EnH_PhiRef1_4_:
                    if i == j:
                        aux = aux + fH_PhiRef1_4_[aux1]
                    aux1 += 1
                fH_PhiRef1_4 = np.append(fH_PhiRef1_4, aux)
                aux = 0
                aux1 = 0

            k = 0  # Variable auxiliar para recorrer en for lop
            for i in Ei_PhiRef2_4:
                if fB_PhiRef2_4[k] < 0.01:
                    break

                Smed_H_PhiRef2_4 = abs(CCi_PhiRef2_4[:, k].conj().transpose() @ OppTH @ CCf_PhiRef2_4) ** 2 * fB_PhiRef2_4[k]
                fH_PhiRef2_4_ = np.append(fH_PhiRef2_4_, Smed_H_PhiRef2_4[Smed_H_PhiRef2_4 > 1e-5].conj().transpose())
                EnH_PhiRef2_4_ = np.append(EnH_PhiRef2_4_, np.around((Eff_PhiRef2_4[Smed_H_PhiRef2_4 > 1e-5] - i) * 100) / 100)

                k += 1

            a, b = np.unique(EnH_PhiRef2_4_, return_index=True)
            EnH_PhiRef2_4 = EnH_PhiRef2_4_[np.sort(b)]
            fH_PhiRef2_4 = np.array([])

            aux = 0
            aux1 = 0

            for i in EnH_PhiRef2_4:
                for j in EnH_PhiRef2_4_:
                    if i == j:
                        aux = aux + fH_PhiRef2_4_[aux1]
                    aux1 += 1
                fH_PhiRef2_4 = np.append(fH_PhiRef2_4, aux)
                aux = 0
                aux1 = 0

        #crear diccionarios
        if Sym == 1:
            # D4h contenedores
            Cont[Container - 1] = {
                "Id": Id,
                "fHPhi1": fHPhi1,
                "EnHPhi1": EnHPhi1,
                "fHPhi2": fHPhi2,
                "EnHPhi2": EnHPhi2,
                "fHPhiRef1": fHPhiRef1,
                "EnHPhiRef1": EnHPhiRef1,
                "fHPhiRef2": fHPhiRef2,
                "EnHPhiRef2": EnHPhiRef2,
                "fr": fr,
                "Enr": Enr,
                "fl": fl,
                "Enl": Enl
                }
        else:
            #D3d
            Cont[Container - 1] = {
                "Id": Id,

                "fHPhi1_1": fHPhi1_1,
                "EnHPhi1_1": EnHPhi1_1,
                "fH_Phi2_1": fH_Phi2_1,
                "EnH_Phi2_1": EnH_Phi2_1,
                "fH_PhiRef1_1": fH_PhiRef1_1,
                "EnH_PhiRef1_1": EnH_PhiRef1_1,
                "fH_PhiRef2_1": fH_PhiRef2_1,
                "EnH_PhiRef2_1": EnH_PhiRef2_1,
                "fr1": fr1,
                "Enr1": Enr1,
                "fl1": fl1,
                "Enl1": Enl1,

                "fH_Phi1_2": fH_Phi1_2,
                "EnH_Phi1_2": EnH_Phi1_2,
                "fH_Phi2_2": fH_Phi2_2,
                "EnH_Phi2_2": EnH_Phi2_2,
                "fH_PhiRef1_2": fH_PhiRef1_2,
                "EnH_PhiRef1_2": EnH_PhiRef1_2,
                "fH_PhiRef2_2": fH_PhiRef2_2,
                "EnH_PhiRef2_2": EnH_PhiRef2_2,
                "fr2": fr2,
                "Enr2": Enr2,
                "fl2": fl2,
                "Enl2": Enl2,

                "fHPhi1_3": fHPhi1_3,
                "EnHPhi1_3": EnHPhi1_3,
                "fH_Phi2_3": fH_Phi2_3,
                "EnH_Phi2_3": EnH_Phi2_3,
                "fH_PhiRef1_3": fH_PhiRef1_3,
                "EnH_PhiRef1_3": EnH_PhiRef1_3,
                "fH_PhiRef2_3": fH_PhiRef2_3,
                "EnH_PhiRef2_3": EnH_PhiRef2_3,
                "fr3": fr3,
                "Enr3": Enr3,
                "fl3": fl3,
                "Enl3": Enl3,

                "fH_Phi1_4": fH_Phi1_4,
                "EnH_Phi1_4": EnH_Phi1_4,
                "fH_Phi2_4": fH_Phi2_4,
                "EnH_Phi2_4": EnH_Phi2_4,
                "fH_PhiRef1_4": fH_PhiRef1_4,
                "EnH_PhiRef1_4": EnH_PhiRef1_4,
                "fH_PhiRef2_4": fH_PhiRef2_4,
                "EnH_PhiRef2_4": EnH_PhiRef2_4,
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
        fHPhi1 = Cont[Container - 1]['fHPhi1']
        EnHPhi1 = Cont[Container - 1]['EnHPhi1'] + ShiftE
        fHPhi2 = Cont[Container - 1]['fHPhi2']
        EnHPhi2 = Cont[Container - 1]['EnHPhi2'] + ShiftE
        fHPhiRef1 = Cont[Container - 1]['fHPhiRef1']
        EnHPhiRef1 = Cont[Container - 1]['EnHPhiRef1'] + ShiftE
        fHPhiRef2 = Cont[Container - 1]['fHPhiRef2']
        EnHPhiRef2 = Cont[Container - 1]['EnHPhiRef2'] + ShiftE
        fr = Cont[Container - 1]['fr']
        Enr = Cont[Container - 1]['Enr'] + ShiftE
        fl = Cont[Container - 1]['fl']
        Enl = Cont[Container - 1]['Enl'] + ShiftE
    else:
        #D3d
        fHPhi1_1 = Cont[Container - 1]['fHPhi1_1']
        EnHPhi1_1 = Cont[Container - 1]['EnHPhi1_1'] + ShiftE
        fH_Phi2_1 = Cont[Container - 1]['fH_Phi2_1']
        EnH_Phi2_1 = Cont[Container - 1]['EnH_Phi2_1'] + ShiftE
        fH_PhiRef1_1 = Cont[Container - 1]['fH_PhiRef1_1']
        EnH_PhiRef1_1 = Cont[Container - 1]['EnH_PhiRef1_1'] + ShiftE
        fH_PhiRef2_1 = Cont[Container - 1]['fH_PhiRef2_1']
        EnH_PhiRef2_1 = Cont[Container - 1]['EnH_PhiRef2_1'] + ShiftE
        fr1 = Cont[Container - 1]['fr1']
        Enr1 = Cont[Container - 1]['Enr1'] + ShiftE
        fl1 = Cont[Container - 1]['fl1']
        Enl1 = Cont[Container - 1]['Enl1'] + ShiftE

        fH_Phi1_2 = Cont[Container - 1]['fH_Phi1_2']
        EnH_Phi1_2 = Cont[Container - 1]['EnH_Phi1_2'] + ShiftE
        fH_Phi2_2 = Cont[Container - 1]['fH_Phi2_2']
        EnH_Phi2_2 = Cont[Container - 1]['EnH_Phi2_2'] + ShiftE
        fH_PhiRef1_2 = Cont[Container - 1]['fH_PhiRef1_2']
        EnH_PhiRef1_2 = Cont[Container - 1]['EnH_PhiRef1_2'] + ShiftE
        fH_PhiRef2_2 = Cont[Container - 1]['fH_PhiRef2_2']
        EnH_PhiRef2_2 = Cont[Container - 1]['EnH_PhiRef2_2'] + ShiftE
        fr2 = Cont[Container - 1]['fr2']
        Enr2 = Cont[Container - 1]['Enr2'] + ShiftE
        fl2 = Cont[Container - 1]['fl2']
        Enl2 = Cont[Container - 1]['Enl2'] + ShiftE

        fHPhi1_3 = Cont[Container - 1]['fHPhi1_3']
        EnHPhi1_3 = Cont[Container - 1]['EnHPhi1_3'] + ShiftE
        fH_Phi2_3 = Cont[Container - 1]['fH_Phi2_3']
        EnH_Phi2_3 = Cont[Container - 1]['EnH_Phi2_3'] + ShiftE
        fH_PhiRef1_3 = Cont[Container - 1]['fH_PhiRef1_3']
        EnH_PhiRef1_3 = Cont[Container - 1]['EnH_PhiRef1_3'] + ShiftE
        fH_PhiRef2_3 = Cont[Container - 1]['fH_PhiRef2_3']
        EnH_PhiRef2_3 = Cont[Container - 1]['EnH_PhiRef2_3'] + ShiftE
        fr3 = Cont[Container - 1]['fr3']
        Enr3 = Cont[Container - 1]['Enr3'] + ShiftE
        fl3 = Cont[Container - 1]['fl3']
        Enl3 = Cont[Container - 1]['Enl3'] + ShiftE

        fH_Phi1_4 = Cont[Container - 1]['fH_Phi1_4']
        EnH_Phi1_4 = Cont[Container - 1]['EnH_Phi1_4'] + ShiftE
        fH_Phi2_4 = Cont[Container - 1]['fH_Phi2_4']
        EnH_Phi2_4 = Cont[Container - 1]['EnH_Phi2_4'] + ShiftE
        fH_PhiRef1_4 = Cont[Container - 1]['fH_PhiRef1_4']
        EnH_PhiRef1_4 = Cont[Container - 1]['EnH_PhiRef1_4'] + ShiftE
        fH_PhiRef2_4 = Cont[Container - 1]['fH_PhiRef2_4']
        EnH_PhiRef2_4 = Cont[Container - 1]['EnH_PhiRef2_4'] + ShiftE
        fr4 = Cont[Container - 1]['fr4']
        Enr4 = Cont[Container - 1]['Enr4'] + ShiftE
        fl4 = Cont[Container - 1]['fl4']
        Enl4 = Cont[Container - 1]['Enl4'] + ShiftE

    #H
    x = np.reshape(x, x.size)
    if Sym==1:
        FFHPhi1 = S2B.Sticks2Band(x, EnHPhi1, fHPhi1, MinL, MaxL, SplitE, 0, Gauss)
        FFHPhi2 = S2B.Sticks2Band(x, EnHPhi2, fHPhi2, MinL, MaxL, SplitE, 0, Gauss)
        FFHPhiRef1 = S2B.Sticks2Band(x, EnHPhiRef1, fHPhiRef1, MinL, MaxL, SplitE, 0, Gauss)
        FFHPhiRef2 = S2B.Sticks2Band(x, EnHPhiRef2, fHPhiRef2, MinL, MaxL, SplitE, 0, Gauss)
        FFr = S2B.Sticks2Band(x, Enr, fr, MinL, MaxL, SplitE, 0, Gauss)
        FFl = S2B.Sticks2Band(x, Enl, fl, MinL, MaxL, SplitE, 0, Gauss)
    else:
        #Sitio 1
        FFH_Phi1_1 = S2B.Sticks2Band(x, EnHPhi1_1, fHPhi1_1, MinL, MaxL, SplitE, 0, Gauss)
        FFH_Phi2_1 = S2B.Sticks2Band(x, EnH_Phi2_1, fH_Phi2_1, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef1_1 = S2B.Sticks2Band(x, EnH_PhiRef1_1, fH_PhiRef1_1, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef2_1 = S2B.Sticks2Band(x, EnH_PhiRef2_1, fH_PhiRef2_1, MinL, MaxL, SplitE, 0, Gauss)
        #Right
        FFr1 = S2B.Sticks2Band(x, Enr1, fr1, MinL, MaxL, SplitE, 0, Gauss)
        #Left
        FFl1 = S2B.Sticks2Band(x, Enl1, fl1, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 2
        FFH_Phi1_2 = S2B.Sticks2Band(x, EnH_Phi1_2, fH_Phi1_2, MinL, MaxL, SplitE, 0, Gauss)
        FFH_Phi2_2 = S2B.Sticks2Band(x, EnH_Phi2_2, fH_Phi2_2, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef1_2 = S2B.Sticks2Band(x, EnH_PhiRef1_2, fH_PhiRef1_2, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef2_2 = S2B.Sticks2Band(x, EnH_PhiRef2_2, fH_PhiRef2_2, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr2 = S2B.Sticks2Band(x, Enr2, fr2, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl2 = S2B.Sticks2Band(x, Enl2, fl2, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 3
        FFH_Phi1_3 = S2B.Sticks2Band(x, EnHPhi1_3, fHPhi1_3, MinL, MaxL, SplitE, 0, Gauss)
        FFH_Phi2_3 = S2B.Sticks2Band(x, EnH_Phi2_3, fH_Phi2_3, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef1_3 = S2B.Sticks2Band(x, EnH_PhiRef1_3, fH_PhiRef1_3, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef2_3 = S2B.Sticks2Band(x, EnH_PhiRef2_3, fH_PhiRef2_3, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr3 = S2B.Sticks2Band(x, Enr3, fr3, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl3 = S2B.Sticks2Band(x, Enl3, fl3, MinL, MaxL, SplitE, 0, Gauss)

        # Sitio 4
        FFH_Phi1_4 = S2B.Sticks2Band(x, EnH_Phi1_4, fH_Phi1_4, MinL, MaxL, SplitE, 0, Gauss)
        FFH_Phi2_4 = S2B.Sticks2Band(x, EnH_Phi2_4, fH_Phi2_4, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef1_4 = S2B.Sticks2Band(x, EnH_PhiRef1_4, fH_PhiRef1_4, MinL, MaxL, SplitE, 0, Gauss)
        FFH_PhiRef2_4 = S2B.Sticks2Band(x, EnH_PhiRef2_4, fH_PhiRef2_4, MinL, MaxL, SplitE, 0, Gauss)
        # Right
        FFr4 = S2B.Sticks2Band(x, Enr4, fr4, MinL, MaxL, SplitE, 0, Gauss)
        # Left
        FFl4 = S2B.Sticks2Band(x, Enl4, fl4, MinL, MaxL, SplitE, 0, Gauss)

    if Sym==1:
        # D4h
        FFXAS = (FFr + FFl)/2
        #A = np.trapz(FFXAS, x)
        A = (np.sum(fr)+np.sum(fl))/2
        FFXAS = FFXAS / A * 3 / 7 * (10 - Nelec)
        FFXAS = FFXAS * XASInt

        FFXMCD = FFl - FFr
        FFXMCD = FFXMCD / A * 3 / 7 * (10 - Nelec)
        FFXMCD = FFXMCD * DichrInt

        FFXMLDPhi1 = FFHPhi1 - FFHPhiRef1
        FFXMLDPhi1 = FFXMLDPhi1 / A * 3 / 7 * (10 - Nelec)

        FFXMLDPhi1 = FFXMLDPhi1 * XMLD1

        FFXMLDPhi2 = FFHPhi2 - FFHPhiRef2
        FFXMLDPhi2 = FFXMLDPhi2 / A * 3 / 7 * (10 - Nelec)

        FFXMLDPhi2 = FFXMLDPhi2 * XMLD2

        mask = np.logical_and(x > 1000, x < 2000)

        xd = x[mask] - 1000
        L = pchip(x, FFXMLDPhi1)
        FF3d = L(xd)

        mask = np.logical_and(x>2000,x<3000)
        xd = x[mask] - 2000
        K = pchip(x, FFXMLDPhi2)
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
        #A = np.trapz(FFXAS, x)
        A = (np.sum(fr1) + np.sum(fl1) + np.sum(fr2) + np.sum(fl2) +np.sum(fr3) + np.sum(fl3) + np.sum(fr4) + np.sum(fl4))/8
        FFXAS = FFXAS / A * 3 / 7 * (10 - Nelec)
        FFXAS = FFXAS * XASInt

        FFXMCD = -(FFr1 - FFl1 + FFr2 - FFl2 + FFr3 - FFl3 + FFr4 - FFl4)/4
        FFXMCD = FFXMCD / A * 3 / 7 * (10 - Nelec)
        FFXMCD = FFXMCD * DichrInt

        FFXMLDPhi1 = (FFH_Phi1_1 - FFH_PhiRef1_1+FFH_Phi1_2 - FFH_PhiRef1_2+FFH_Phi1_3 - FFH_PhiRef1_3+FFH_Phi1_4 - FFH_PhiRef1_4)/4
        FFXMLDPhi1 = FFXMLDPhi1 / A * 3 / 7 * (10 - Nelec)

        FFXMLDPhi1 = FFXMLDPhi1 * XMLD1

        FFXMLDPhi2 = (FFH_Phi2_1 - FFH_PhiRef2_1+FFH_Phi2_2 - FFH_PhiRef2_2+FFH_Phi2_3 - FFH_PhiRef2_3+FFH_Phi2_4 - FFH_PhiRef2_4)/4
        FFXMLDPhi2 = FFXMLDPhi2 / A * 3 / 7 * (10 - Nelec)

        FFXMLDPhi2 = FFXMLDPhi2 * XMLD2

        mask = np.logical_and(x > 1000, x < 2000)

        xd = x[mask] - 1000
        L = pchip(x, FFXMLDPhi1)
        FF3d = L(xd)

        mask = np.logical_and(x>2000,x<3000)
        xd = x[mask] - 2000
        K = pchip(x, FFXMLDPhi2)
        FF4d = K(xd)

        xd = x[x>3000]-3000
        M=pchip(x,FFXAS)
        FF5d = M(xd)

        FF = FFXMCD[x < 1000]
        FF = np.append(FF, FF3d)
        FF = np.append(FF, FF4d)
        FF = np.append(FF,FF5d)

        return FF.real
