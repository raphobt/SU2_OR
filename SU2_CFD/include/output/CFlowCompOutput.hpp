/*!
 * \file CFlowCompOutput.hpp
 * \brief  Headers of the compressible flow output.
 * \author R. Sanchez, T. Albring.
 * \version 7.5.0 "Blackbird"
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

#pragma once

#include "CFlowOutput.hpp"

#include "../../include/solvers/CSolver.hpp"
#include "../../include/fluid/CSWTable.hpp"

class CVariable;

/*! \class CFlowCompOutput
 *  \brief Output class for compressible flow problems.
 *  \author R. Sanchez, T. Albring.
 *  \date May 30, 2018.
 */
class CFlowCompOutput final: public CFlowOutput {
private:
  TURB_MODEL turb_model; //!< Kind of turbulence model

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CFlowCompOutput(const CConfig *config, unsigned short nDim);

  /*!
   * \brief Load the history output field values
   * \param[in] config - Definition of the particular problem.
   */
  void LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) override;

  /*!
   * \brief Set the available volume output fields
   * \param[in] config - Definition of the particular problem.
   */
  void SetVolumeOutputFields(CConfig *config) override;

  /*!
   * \brief Set the values of the volume output fields for a point.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - The container holding all solution data.
   * \param[in] iPoint - Index of the point.
   */
  void LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint) override;

/**********************/

  /*!
   * \brief Compute entropy generation rate by direct dissipation SDD
   * \return Value of SDD at the node
   */
  template<class T>
  su2double Get_SDD_compressible(const T& VelocityGradient, CSolver **solver, unsigned long iPoint) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();

    // Make a 3D copy of the gradient so we do not have worry about nDim

    su2double Grad_Vel[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

    // Definition of the modulus of the mean rate of strain tensor Smean
    su2double tauScaDoubleGradUTr = 2 * ( pow(Grad_Vel[0][0], 2) + pow(Grad_Vel[1][1], 2) + pow(Grad_Vel[2][2], 2) ) \
    + pow(Grad_Vel[0][1] + Grad_Vel[1][0], 2) \
    + pow(Grad_Vel[0][2] + Grad_Vel[2][0], 2) \
    + pow(Grad_Vel[1][2] + Grad_Vel[2][1], 2) \
    - 2/3 * pow(Grad_Vel[0][0] + Grad_Vel[1][1] + Grad_Vel[2][2], 2);

    // Entropy production rate by direct dissipation SDD (Kock and Herwig, 2004-2006)
    // SDD = molecular_visco / T * ||S||
    
    su2double SDD_compressible = Node_Flow->GetLaminarViscosity(iPoint) / Node_Flow->GetTemperature(iPoint) * tauScaDoubleGradUTr;
    
    return SDD_compressible;
  }

  /*!
   * \brief Compute entropy generation rate by direct dissipation SDD
   * \return Value of SDD at the node
   */
  template<class T>
  su2double Get_SDD(const T& VelocityGradient, CSolver **solver, unsigned long iPoint) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();

    // Make a 3D copy of the gradient so we do not have worry about nDim

    su2double Grad_Vel[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

    // Definition of the modulus of the mean rate of strain tensor Smean
    su2double SmeanModulus = pow(2 * ( pow(Grad_Vel[0][0], 2) + pow(Grad_Vel[1][1], 2) + pow(Grad_Vel[2][2], 2) ) \
    + pow(Grad_Vel[0][1] + Grad_Vel[1][0], 2) \
    + pow(Grad_Vel[0][2] + Grad_Vel[2][0], 2) \
    + pow(Grad_Vel[1][2] + Grad_Vel[2][1], 2), 0.5);

    // Entropy production rate by direct dissipation SDD (Kock and Herwig, 2004-2006)
    // SDD = molecular_visco / T * ||S||
    
    su2double SDD = Node_Flow->GetLaminarViscosity(iPoint) / Node_Flow->GetTemperature(iPoint) * SmeanModulus;
    
    return SDD;
  }

  /*!
   * \brief Compute entropy generation rate by indirect dissipation SID
   * \return Value of SID at the node
   */
  //template<class T>
  //su2double Get_SID(const T& VelocityGradient, CSolver **solver, unsigned long iPoint, CConfig *config) const {
  su2double Get_SID(CSolver **solver, unsigned long iPoint, CConfig *config) const {

    const auto* turb_solver = solver[TURB_SOL];
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    const auto* Node_Turb = (config->GetKind_Turb_Model() != TURB_MODEL::NONE) ? turb_solver->GetNodes() : nullptr;
    
    su2double betaStar = 0.09; // Impove: retrieve the betaStar value from the turbulence model (SST k-om)
    su2double omega = Node_Turb->GetSolution(iPoint, 1);
    su2double k = Node_Turb->GetSolution(iPoint, 0);
    su2double epsilon = betaStar * omega * k;

    // Entropy production rate by indirect (turbulent flow) dissipation SID (Kock and Herwig, 2004-2006)
    // SID = rho * epsilon / T

    su2double SID = Node_Flow->GetSolution(iPoint, 0) * epsilon / Node_Flow->GetTemperature(iPoint);
    
    return SID;
  }

  /*!
   * \brief Compute entropy generation rate by the mean flow temperature gradients SDT
   * \return Value of SDT at the node
   */
  template<class T>
  su2double Get_SDT(const T& TemperatureGradient, CSolver **solver, unsigned long iPoint) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    
    // Make a 3D copy of the gradient so we do not have worry about nDim
    su2double Grad_Temp[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Temp[iDim][jDim] = TemperatureGradient[iDim][jDim];

    // Entropy production rate by SDT (Kock and Herwig, 2004-2006)
    // SDT =  ( thermal_cond / T*T ) * || grad(T) || ** 2 
    
    su2double SDT = ( Node_Flow->GetThermalConductivity(iPoint) / Node_Flow->GetTemperature(iPoint) / Node_Flow->GetTemperature(iPoint) ) * ( Grad_Temp[0][0] * Grad_Temp[0][0] + Grad_Temp[0][1] * Grad_Temp[0][1] + Grad_Temp[0][2] * Grad_Temp[0][2] );
    
    return SDT;
  }

  /*!
   * \brief Compute entropy generation rate by the instantaneous flow temperature gradients SIT
   * \return Value of SIT at the node. Here GetSpecificHeatCp() is defined in CSWTable.cpp as for turb visc and thermal cond from SW table
   */
  template<class T>
  su2double Get_SIT(const T& TemperatureGradient, CSolver **solver, unsigned long iPoint, CConfig *config, CGeometry *geometry) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    //nPointDomain = geometry->GetnPointDomain();
    //vector<su2double> lambda_turb[nPointDomain];
    su2double lambda_turb;
    
    // Make a 3D copy of the gradient so we do not have worry about nDim
    su2double Grad_Temp[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Temp[iDim][jDim] = TemperatureGradient[iDim][jDim];
   
    if( Node_Flow->GetSpecificHeatCp(iPoint) > 0.0 ){ // && Node_Flow->GetSpecificHeatCp(iPoint) < 3.0e4 ){
      //lambda_turb = Node_Flow->GetEddyViscosity(iPoint) * Node_Flow->GetSpecificHeatCp(iPoint) / config->GetPrandtl_Turb();
      lambda_turb = Node_Flow->GetEddyViscosity(iPoint) * Node_Flow->GetSpecificHeatCp(iPoint) / 0.85;
    }
    //else{
      //cout << "cp negative for computing turbulent thermal conductivity: cp = " << Node_Flow->GetSpecificHeatCp(iPoint) << " J/kg/K"  << endl;
      //cout << "Turbulent thermal conductivity: lambda_turb = " << lambda_turb << " W/m/K"  << endl;
      //cout << "Values close to 0 for SIT may not be correct..."  << endl;
    //}

    // Entropy production rate by SIT (Kock and Herwig, 2004-2006)
    // SIT =  ( turb_thermal_cond / T*T ) * || grad(T) || ** 2 
    
    su2double SIT = ( lambda_turb / Node_Flow->GetTemperature(iPoint) / Node_Flow->GetTemperature(iPoint) ) * ( Grad_Temp[0][0] * Grad_Temp[0][0] + Grad_Temp[0][1] * Grad_Temp[0][1] + Grad_Temp[0][2] * Grad_Temp[0][2] );
    
    return SIT;
  }

  /*!
   * \brief Compute total enthalpy, total pressure and total temperature of the flow 
   * \return Value of h0 at the node.
   */
  su2double Get_h0(CSolver **solver, unsigned long iPoint) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    
    int flag;
    double pp, TT, cc, x_out, a_out, dummy;
    double Eref_SW=506779.92063833564;
    su2double vv = 1.0/Node_Flow->GetSolution(iPoint, 0);
    su2double ee;
    if (nDim == 3){
    ee = Node_Flow->GetSolution(iPoint, 4)*vv - Eref_SW - Node_Flow->GetVelocity2(iPoint)/2.0; // internal energy
    } else {
    ee = Node_Flow->GetSolution(iPoint, 3)*vv - Eref_SW - Node_Flow->GetVelocity2(iPoint)/2.0;
    }
    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,&ee,&vv,&flag);

    // specific total enthalpy (based on total energy) 
    su2double h0=Node_Flow->GetSolution(iPoint, 4)*vv+pp*vv;

    return h0;
  }

  /*!
   * \brief Compute total pressure or total temperature of the flow depending on choice value
   * \return Value of p0 or T0 at the node.
   */
  su2double Get_p0_T0(CSolver **solver, unsigned long iPoint, int choice) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    
    int Niter, flag, exitflag, MODE=5;
    double pp, TT, ss, cc, x_out, a_out, dummy, resnorm, guess_1, guess_2;
    double Eref_SW=506779.92063833564;
    
    su2double T_out, v_out, energy, hh;
    su2double vv = 1.0/Node_Flow->GetSolution(iPoint, 0);
    su2double ee;
    if (nDim == 3){ ee = Node_Flow->GetSolution(iPoint, 4)*vv - Eref_SW - Node_Flow->GetVelocity2(iPoint)/2.0; // internal energy
    } 
    else { ee = Node_Flow->GetSolution(iPoint, 3)*vv - Eref_SW - Node_Flow->GetVelocity2(iPoint)/2.0;
    }
    
    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,&ee,&vv,&flag);

    // specific entropy
    __transprop_MOD_entropyco2(&ss, &vv, &dummy, &x_out, &TT, &pp, &flag);

    hh=Get_h0(solver, iPoint)-Eref_SW;
  
    guess_1 = Node_Flow->GetTemperature(iPoint);//+Node_Flow->GetVelocity2(iPoint)/2.0/19000.0;
    guess_2 = 1/Node_Flow->GetSolution(iPoint, 0);

    __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&hh, &ss, &guess_1, &guess_2);
    if (Niter>=500 || T_out!=T_out || v_out!=v_out || v_out<=0.0 || T_out<=100.0){
    for(int i=1;i<20; i++){
      guess_2=guess_2*1.001;
      guess_1=guess_1*1.001;
      __non_linear_solvers_MOD_new_rap2d(&MODE, &T_out, &v_out, &resnorm, &Niter, &exitflag,&hh, &ss, &guess_1, &guess_2);
      if (Niter<500 & T_out==T_out & v_out==v_out & v_out>0.0 & T_out>100.0){
          break;
      }
    }
    }
    if (Niter>=500 ){ cout << "Max iteration reached in SetTDState_hs in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
    if (T_out!=T_out){ cout << "NAN T in SetTDState_hs in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
    if (v_out!=v_out){ cout << "NAN v in SetTDState_hs in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
    if (v_out<=0.0){ cout << "Negative v in SetTDState_hs : v = " << v_out << " in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
    if (T_out<=100.0){ cout << "Too low T in SetTDState_hs : T = " << T_out << " in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
  
    __properties_MOD_inter_energy(&T_out, &v_out, &energy);
    if (energy!=energy){ cout << "NAN e in SetTDState_hs in CFlowCompOutput.hpp for p0, T0 field values" << endl;
    }
    
    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,&energy,&v_out,&flag);

    su2double p0 = pp;
    su2double T0 = TT;

    if (choice == 0){ return p0;
    }
    else if (choice == 1){ return T0;
    }
    else { 
      cout << "Check the choice int value in CFlowCompOutput.cpp to output p0, T0 field values..."<< endl;
      return 0.0;
    }

  }

  /*!
   * \brief Output the mean rate of strain tensor and its components
   * \return Value of SmeanModulus and components at the node
   */
  template<class T>
  su2double Get_SmeanModulus(const T& VelocityGradient, CSolver **solver, int choiceS) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();

    // Make a 3D copy of the gradient so we do not have worry about nDim

    su2double Grad_Vel[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = VelocityGradient[iDim][jDim];

    // Definition of the modulus of the mean rate of strain tensor Smean
  /*  su2double SmeanModulus = pow(2 * ( pow(Grad_Vel[0][0], 2) + pow(Grad_Vel[1][1], 2) + pow(Grad_Vel[2][2], 2) ) \
    + pow(Grad_Vel[0][1] + Grad_Vel[1][0], 2) \
    + pow(Grad_Vel[0][2] + Grad_Vel[2][0], 2) \
    + pow(Grad_Vel[1][2] + Grad_Vel[2][1], 2), 0.5);*/

    su2double duxdx2 = Grad_Vel[0][0] * Grad_Vel[0][0];
    su2double duydy2 = Grad_Vel[1][1] * Grad_Vel[1][1];
    su2double duzdz2 = Grad_Vel[2][2] * Grad_Vel[2][2];
    su2double duxdy_duydx2 = (Grad_Vel[0][1] + Grad_Vel[1][0]) * (Grad_Vel[0][1] + Grad_Vel[1][0]); 
    su2double duxdz_duzdx2 = (Grad_Vel[0][2] + Grad_Vel[2][0]) * (Grad_Vel[0][2] + Grad_Vel[2][0]); 
    su2double duydz_duzdy2 = (Grad_Vel[1][2] + Grad_Vel[2][1]) * (Grad_Vel[1][2] + Grad_Vel[2][1]); 

    su2double SmeanModulus = pow(2 * (duxdx2 + duydy2 + duzdz2) + duxdy_duydx2 + duxdz_duzdx2 + duydz_duzdy2, 0.5);    

    if (choiceS == 0){ return SmeanModulus;
    }
    else if (choiceS == 1){ return duxdx2; 
    }
    else if (choiceS == 2){ return duydy2; 
    }
    else if (choiceS == 3){ return duzdz2; 
    }
    else if (choiceS == 4){ return duxdy_duydx2; 
    }
    else if (choiceS == 5){ return duxdz_duzdx2; 
    }
    else if (choiceS == 6){ return duydz_duzdy2; 
    }
    else { 
      cout << "Check the choiceS int value in CFlowCompOutput.cpp to output SmeanModulus and component values..."<< endl;
      return 0.0;
    }

  }

  /*!
   * \brief Compute temperature gradient modulus and its components
   * \return Value of temperature gradient modulus and components at the node.
   */
  template<class T>
  su2double Get_gradTmodulus(const T& TemperatureGradient, CSolver **solver, int choiceT) const {
    
    const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
    
    // Make a 3D copy of the gradient so we do not have worry about nDim
    su2double Grad_Temp[3][3] = {{0.0}};

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Temp[iDim][jDim] = TemperatureGradient[iDim][jDim];
 
    su2double dTdx2 = Grad_Temp[0][0]*Grad_Temp[0][0];
    su2double dTdy2 = Grad_Temp[0][1]*Grad_Temp[0][1];
    su2double dTdz2 = Grad_Temp[0][2]*Grad_Temp[0][2];
 
    su2double gradTmodulus = pow(dTdx2 + dTdy2 + dTdz2, 0.5);
    
    if (choiceT == 0){ return gradTmodulus;
    }
    else if (choiceT == 1){ return dTdx2; 
    }
    else if (choiceT == 2){ return dTdy2; 
    }
    else if (choiceT == 3){ return dTdz2; 
    }
    else { 
      cout << "Check the choiceT int value in CFlowCompOutput.cpp to output gradTmodulus and component values..."<< endl;
      return 0.0;
    }


  }
/**********************/

  /*!
   * \brief Set the available history output fields
   * \param[in] config - Definition of the particular problem.
   */
  void SetHistoryOutputFields(CConfig *config) override;

  /*!
   * \brief Check whether the base values for relative residuals should be initialized
   * \param[in] config - Definition of the particular problem.
   * \return <TRUE> if the residuals should be initialized.
   */
  bool SetInit_Residuals(const CConfig *config) override;

  /*!
   * \brief Write any additional output defined for the current solver.
   * \param[in] config - Definition of the particular problem per zone.
   */
  void SetAdditionalScreenOutput(const CConfig *config) override;

  /*!
   * \brief Determines if the history file output.
   * \param[in] config - Definition of the particular problem.
   */
  bool WriteHistoryFile_Output(const CConfig *config) override;
};
