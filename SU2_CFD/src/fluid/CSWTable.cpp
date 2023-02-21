/*!
 * \file CPengRobinsoPerso.cpp
 * \brief Source of the Peng-Robinson model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/fluid/CSWTable.hpp"

// ATTENTION VERIFIER CE QUI EST INITIQLISE DANS IDEALGAS
CSWTable::CSWTable()
    : CIdealGas(std::nan(""), std::nan("")) {
        // ATTENTION TOUT CA CEST DEFINI POUR LE CO2. ON S'EN SERT JUSTE COMME GUESS VALUES POUR LA TABLE
        double Tstar = 304.128;
        double Pstar = 7377330;
        double w = 0.228;
        double R = 188.9241;
        
        a = 0.45724 * R * R * Tstar * Tstar / Pstar;
        b = 0.0778 * R * Tstar / Pstar;
        TstarCrit = Tstar;
        Zed = 1.0;

        if (w <= 0.49)
          k = 0.37464 + 1.54226 * w - 0.26992 * w * w;
        else
          k = 0.379642 + 1.48503 * w - 0.164423 * w * w + 0.016666 * w * w * w;
}

void CSWTable::SetTDState_rhoe(su2double rho, su2double e) {
    
  //su2double Cv;
    
  Density = rho;
  StaticEnergy = e;

  AD::StartPreacc();
  AD::SetPreaccIn(rho);
  AD::SetPreaccIn(e);
    
    double Eref_SW=506779.92063833564;
    double Sref_SW=2739.05;
    su2double vv=1.0/ rho;
    su2double ee= e-Eref_SW;
    double R = 188.9241;
    
    double pp, TT, cc, ss, x_out, a_out, dummy, cv_out, vp_out, dp_du_v, dp_dv_u, dedr_T, res2,res3, res4, cp_out;
    int flag;

    if(e!=e){
        cout << "NAN in SetTDState_rhoe : e" << endl;
    }
    if(e<-400e3){
        cout << "Too low e in SetTDState_rhoe : e = " << e << endl;
    }
    if(rho!=rho){
        cout << "NAN in SetTDState_rhoe : rho" << endl;
    }
    if(rho<=0.0){
        cout << "Negative rho in SetTDState_rhoe : rho = " << rho << endl;
    }
    
    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,&ee,&vv,&flag);
    Temperature = (su2double) TT;
    Pressure = (su2double) pp;
    Quality = x_out;
    Flag_loca = flag;
    __transprop_MOD_entropyco2(&ss, &vv, &vp_out, &x_out, &TT, &pp, &flag);
    Entropy=ss+Sref_SW;
    __transprop_MOD_cvco2(&cv_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
    Cv=(su2double) cv_out;
    SoundSpeed2 = cc*cc;
    
    __derivees_MOD_co2der(&dp_dv_u, &dp_du_v , &ee, &vv, &TT, &pp, &res2, &res3, &res4);
    __transprop_MOD_dedrco2(&dedr_T, &vv,&vp_out,&x_out, &TT, &pp, &flag);
    
    dTde_rho = 1 / Cv; // ok par definition
    dTdrho_e = -dTde_rho*dedr_T;
    dPdrho_e = -vv * vv* dp_dv_u;
    dPde_rho = dp_du_v;
    
    __transprop_MOD_cpco2(&cp_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
    //Cp = std::min(cp_out,1e4);
    Cp=cp_out;
    
    //Utilise juste pour les guess
    Zed = Pressure / (R * Temperature * rho);

  AD::SetPreaccOut(Temperature);
  AD::SetPreaccOut(SoundSpeed2);
  AD::SetPreaccOut(dPde_rho);
  AD::SetPreaccOut(dPdrho_e);
  AD::SetPreaccOut(dTde_rho);
  AD::SetPreaccOut(dTdrho_e);
  AD::SetPreaccOut(Zed);
  AD::SetPreaccOut(Cp);
  AD::SetPreaccOut(Cv);
  AD::SetPreaccOut(Quality);
  AD::SetPreaccOut(Pressure);
  AD::SetPreaccOut(Entropy);
  AD::EndPreacc();
}

void CSWTable::SetTDState_PT(su2double P, su2double T) {

    su2double toll = 1e-6;
    su2double A, B, Z, DZ = 1.0, F, F1;
    su2double rho, e, alpha2T;
    unsigned short nmax = 20, count = 0;
    double R = 188.9241;
    double Eref_SW=506779.92063833564;

    AD::StartPreacc();
    AD::SetPreaccIn(P);
    AD::SetPreaccIn(T);

    // Guess on rho using Peng-Robinson //
    
    alpha2T = (1 + k * (1 - sqrt(T / TstarCrit))) * (1 + k * (1 - sqrt(T / TstarCrit)));
    A = a * alpha2T * P / (T * R) / (T * R);
    B = b * P / (T * R);

    if (Zed > 0.1)
      Z = min(Zed, 0.99);
    else
      Z = 0.99;

    do {
      F = Z * Z * Z + Z * Z * (B - 1.0) + Z * (A - 2 * B - 3 * B * B) + (B * B * B + B * B - A * B);
      F1 = 3 * Z * Z + 2 * Z * (B - 1.0) + (A - 2 * B - 3 * B * B);
      DZ = F / F1;
      Z -= DZ;
    } while (abs(DZ) > toll && count < nmax);

    if (count == nmax) {
      cout << "Warning Newton-Raphson exceed number of max iteration in PT" << endl;
      cout << "Compressibility factor  " << Z << " would be substituted with " << Zed << endl;
    }
    // check if the solution is physical otherwise uses previous point  solution
    if (Z <= 1.0001 && Z >= 0.05 && count < nmax) Zed = Z;

    //Guess
    rho = P / (Zed * R * T);
    
    int MODE=3;
    su2double energy, v_out, out_2, p_in, T_in, vguess;
    double resnorm,out3;
    int Niter, exitflag;
    
    p_in = P;
    T_in = T;
    //vguess = 1.0/400;
    vguess = 1.0/rho;
    
    // we find v as a function of p and T with a guess on v
    // attention normalement c'est bon parce que c'est utilise just a l'inlet qui est supercritique mais ici on
    // tient pas compte du 2 phase donc si jamais ce sera faux
    __non_linear_solvers_MOD_new_rap1d(&MODE, &v_out, &out_2, &resnorm, &Niter, &exitflag, &p_in, &vguess, &T_in, &out3);
    
    if (Niter>=500 || v_out!=v_out|| v_out<=0.0){
        cout << "Wrong v in SetTDState_PT : v = " << v_out << endl;
    }
    
    // pareil ici
    __properties_MOD_inter_energy(&T_in, &v_out, &energy);
    
    if (energy!=energy){
        cout << "NAN for e in SetTDState_PT" << endl;
    }
    
    rho = 1.0/v_out;
    e = energy+Eref_SW;
    
  AD::SetPreaccOut(rho);
  AD::SetPreaccOut(e);
  AD::EndPreacc();

  SetTDState_rhoe(rho, e);
}

// passe par ici ok
void CSWTable::SetTDState_Prho(su2double P, su2double rho) {
    
  SetEnergy_Prho(P, rho);

  SetTDState_rhoe(rho, StaticEnergy);
}

void CSWTable::SetTDState_hs(su2double h, su2double s) {
    
    // GUESS values with Peng-Robinson: a tester mais si ca marche comme ca on peut gagner du temps
    // ...
    //
    
    int MODE=5;
    su2double T_out, v_out, energy, hh, ss;
    double resnorm, guess_1, guess_2;
    int Niter, exitflag;
    double Eref_SW=506779.92063833564;
    double Sref_SW=2739.05;
    hh=h-Eref_SW;
    ss=s-Sref_SW;
    
    //guess_1=313.0;
    //guess_2=1.0/200.0;
    
    /* OP1 8 MPa - 308 K */
    guess_1=308.0;
    guess_2=1.0/436.24;
    /* OP2 8 MPa - 318 K */
    //guess_1=318.0;
    //guess_2=1.0/241.86;
    /* OP3 10 MPa - 308 K */
    //guess_1=308.0;
    //guess_2=1.0/714.84;
    /* OP4 10 MPa - 318 K */
    //guess_1=318.0;
    //guess_2=1.0/502.56;
    /* OP6 9 MPa - 310 K */
    //guess_1=310.0;
    //guess_2=1.0/614.87;
    /* case 4 Lettieri 8 MPa - 311 K */
    //guess_1=311.0;
    //guess_2=1.0/306.7;
    
    __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&hh, &ss, &guess_1, &guess_2);
    if (Niter>=500 || T_out!=T_out || v_out!=v_out || v_out<=0.0 || T_out<=100.0){
        for(int i=1;i<20; i++){
            guess_2=guess_2/1.1;
            __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&hh, &ss, &guess_1, &guess_2);
            if (Niter<500 & T_out==T_out & v_out==v_out & v_out>0.0 & T_out>100.0){
                break;
            }
        }
    }
    if (Niter>=500 ){
        cout << "Max iteration reached in SetTDState_hs" << endl;
    }
    if (T_out!=T_out){
        cout << "NAN T in SetTDState_hs" << endl;
    }
    if (v_out!=v_out){
        cout << "NAN v in SetTDState_hs" << endl;
    }
    if (v_out<=0.0){
        cout << "Negative v in SetTDState_hs : v = " << v_out << endl;
    }
    if (T_out<=100.0){
        cout << "Too low T in SetTDState_hs : T = " << T_out << endl;
    }
    __properties_MOD_inter_energy(&T_out, &v_out, &energy);
    
    if (energy!=energy){
        cout << "NAN e in SetTDState_hs" << endl;
    }
    
    su2double rho=(su2double) 1.0/v_out;
    su2double e=(su2double)energy+Eref_SW;
    
    SetTDState_rhoe(rho, e);

}

void CSWTable::SetEnergy_Prho(su2double P, su2double rho) {
    
    AD::StartPreacc();
    AD::SetPreaccIn(P);
    AD::SetPreaccIn(rho);

    
    su2double A, B, C, T, vb1, vb2;
    double R = 188.9241;

    vb1 = (1 / rho - b);
    vb2 = (1 / rho / rho + 2 * b / rho - b * b);

    A = R / vb1 - a * k * k / TstarCrit / vb2;

    B = 2 * a * k * (k + 1) / sqrt(TstarCrit) / vb2;

    C = -P - a * (1 + k) * (1 + k) / vb2;

    T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
    T *= T;
    
    int MODE=4;
    int MODE2=3;
    su2double T_out, v_in, out_2, energy, x_out,v_v,v_l,Tsat,u_v,u_l, p_in, Tguess;
    double resnorm, out3;
    double Pcrit=7377330;
    int Niter, exitflag;
    double Eref_SW=506779.92063833564;
    
    p_in = P;
    v_in = 1.0/rho;
    Tguess=T;
    
    // on calcule T en fonction de P et v avec un guess sur T
    // ici il faut tenir compte du 2 phase
    //checker si c'est 2 phase et si jamais Tsat(p)
    
    __properties_MOD_satprop(&MODE2, &p_in, &Tsat, &v_v, &v_l, &u_v, &u_l);
    
    if(v_in<v_v & v_in>v_l& p_in<Pcrit){
        __properties_MOD_inter_energy(&Tsat,&v_v,&u_v);
        __properties_MOD_inter_energy(&Tsat,&v_l,&u_l);
        x_out   = (v_in - v_l) / (v_v - v_l);
        energy = x_out * u_v + (1.0 - x_out) * u_l;
        T_out=Tsat;
    }
    else{
        __non_linear_solvers_MOD_new_rap1d(&MODE, &T_out, &out_2, &resnorm, &Niter, &exitflag, &p_in, &Tguess, &v_in, &out3);
        
        if (Niter>=500 || T_out!=T_out || T_out<=100.0){
            Tguess=270;
            for(int i=1;i<20; i++){
                Tguess=Tguess/1.01;
                __non_linear_solvers_MOD_new_rap1d(&MODE, &T_out, &out_2, &resnorm, &Niter, &exitflag, &p_in, &Tguess, &v_in, &out3);
                if (Niter<500 & T_out==T_out & T_out>100.0){
                    break;
                }
            }
        }
        __properties_MOD_inter_energy(&T_out, &v_in, &energy);
    }
    
    if(T_out!=T_out){
        cout << "NAN in T(p, v) in SetEnergy_Prho" << endl;
    }

    if(energy!=energy){
        cout << "NAN in e(T, v) in SetEnergy_Prho" << endl;
        T_out=300;
        __properties_MOD_inter_energy(&T_out, &v_in, &energy);
    }
    
    if(T_out<100.0){
        cout << "T_out too small in T(p, v) in SetEnergy_Prho" << endl;
    }

    StaticEnergy = (su2double) energy+Eref_SW;

    AD::SetPreaccOut(StaticEnergy);
    AD::EndPreacc();
    
}

void CSWTable::SetTDState_rhoT(su2double rho, su2double T) {
    
    su2double energy;
    su2double T_in = T;
    su2double v_in =1.0/rho;
    __properties_MOD_inter_energy(&T_in, &v_in, &energy);
  SetTDState_rhoe(rho, energy);
}

void CSWTable::SetTDState_Ps(su2double P, su2double s) {
    
    int MODE=6;
    su2double T_out, v_out, energy;
    double resnorm, guess_1, guess_2;
    int Niter, exitflag;
    
    //guess_1=313.0;
    //guess_2=1.0/250.0;
    
    /* OP1 8 MPa - 308 K */
    guess_1=308.0;
    guess_2=1.0/436.24;
    /* OP2 8 MPa - 318 K */
    //guess_1=318.0;
    //guess_2=1.0/241.86;
    /* OP3 10 MPa - 308 K */
    //guess_1=308.0;
    //guess_2=1.0/714.84;
    /* OP4 10 MPa - 318 K */
    //guess_1=318.0;
    //guess_2=1.0/502.56;
    /* OP6 9 MPa - 310 K */
    //guess_1=310.0;
    //guess_2=1.0/614.87;
    /* case 4 Lettieri 8 MPa - 311 K */
    //guess_1=311.0;
    //guess_2=1.0/306.7;
    
    __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&P, &s, &guess_1, &guess_2);
    if (Niter>=500 || T_out!=T_out || v_out!=v_out || v_out<=0.0 || T_out<=100.0){
        for(int i=1;i<10; i++){
            guess_2=guess_2/1.4;
            __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&P, &s, &guess_1, &guess_2);
            if (Niter<500 & T_out==T_out & v_out==v_out & v_out>0.0 & T_out>100.0){
                break;
            }
        }
    }
    if (Niter>=500 ){
        cout << "Max iteration reached in SetTDState_Ps" << endl;
    }
    if (T_out!=T_out){
        cout << "NAN T in SetTDState_Ps" << endl;
    }
    if (v_out!=v_out){
        cout << "NAN v in SetTDState_Ps" << endl;
    }
    if (v_out<=0.0){
        cout << "Negative v in SetTDState_Ps : v = " << v_out << endl;
    }
    if (T_out<=100.0){
        cout << "Too low T in SetTDState_Ps : T = " << T_out << endl;
    }
    __properties_MOD_inter_energy(&T_out, &v_out, &energy);
    
    if (energy!=energy){
        cout << "NAN for e in SetTDState_Ps" << endl;
    }
    
    su2double rho = (su2double) 1.0/v_out;
    su2double e = (su2double)energy;
    
    SetTDState_rhoe(rho, e);
}

void CSWTable::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  su2double dPdT_rho, dPdrho_T, dPds_rho;

  SetTDState_Prho(P, rho);
    
  su2double temp = Temperature;
  su2double vol = 1.0/rho;
  double dpdv_T;
  double dpdT_v;
    
  __properties_MOD_dpdv_t(&temp, &vol, &dpdv_T);
  dPdrho_T = -1.0/rho/rho*dpdv_T;
  __properties_MOD_dpdt_v(&temp, &vol, &dpdT_v);
  dPdT_rho = dpdT_v;

  dhdrho_P = -dPdrho_e / dPde_rho - P / rho / rho;
  dhdP_rho = 1.0 / dPde_rho + 1.0 / rho;
  dPds_rho = rho * rho * (SoundSpeed2 - dPdrho_T) / dPdT_rho; // ???
  dsdP_rho = 1.0 / dPds_rho;
  dsdrho_P = -SoundSpeed2 / dPds_rho;
}

su2double CSWTable::GetLaminarViscosity() {
    double vis;
    su2double vv = 1.0/Density;
    su2double TT = Temperature;
    su2double x_in = Quality;
    su2double pp = Pressure;
    int flag = Flag_loca;
    __transprop_MOD_co2visco(&vis, &vv, &TT, &x_in, &pp, &flag);
    Mu = vis;
    dmudrho_T = 0.0;
    dmudT_rho = 0.0;
    
    if(vis<0.0){
        cout << "Negative mu in GetLaminarViscosity : mu = " << vis << endl;
    }
    if(vis!=vis){
        cout << "NaN for mu in GetLaminarViscosity" << endl;
    }
    return Mu;
  }

su2double CSWTable::GetThermalConductivity() {
    double lambda_out,vp_out;
    su2double vv = 1.0/Density;
    su2double TT = Temperature;
    su2double x_in = Quality;
    su2double pp = Pressure;
    int flag = Flag_loca;
  __transprop_MOD_co2conduc2phase(&lambda_out, &vv, &vp_out, &x_in, &TT, &pp, &flag);
    Kt=lambda_out;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;
    if(lambda_out<0.0){
        cout << "Negative Kt in GetLaminarViscosity : Kt = " << lambda_out << endl;
    }
    if(lambda_out!=lambda_out){
        cout << "NaN for Kt in GetLaminarViscosity" << endl;
    }
  return Kt;
}


