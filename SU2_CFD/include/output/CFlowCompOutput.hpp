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
    
    su2double SDT = ( Node_Flow->GetThermalConductivity(iPoint) / Node_Flow->GetTemperature(iPoint) / Node_Flow->GetTemperature(iPoint) ) * ( Grad_Temp[0][0] * Grad_Temp[0][0] + Grad_Temp[1][0] * Grad_Temp[1][0] + Grad_Temp[2][0] * Grad_Temp[2][0] );
    
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
   
    if( Node_Flow->GetSpecificHeatCp(iPoint) > 0.0 && Node_Flow->GetSpecificHeatCp(iPoint) < 3.0e4 ){
      lambda_turb = Node_Flow->GetEddyViscosity(iPoint) * Node_Flow->GetSpecificHeatCp(iPoint) / config->GetPrandtl_Turb();
    }
    else{
      cout << "cp out of range for computing turbulent thermal conductivity: cp = " << Node_Flow->GetSpecificHeatCp(iPoint) << " J/kg/K"  << endl;
      cout << "Turbulent thermal conductivity: lambda_turb = " << lambda_turb << " W/m/K"  << endl;
    }

    // Entropy production rate by SIT (Kock and Herwig, 2004-2006)
    // SIT =  ( turb_thermal_cond / T*T ) * || grad(T) || ** 2 
    
    su2double SIT = ( lambda_turb / Node_Flow->GetTemperature(iPoint) / Node_Flow->GetTemperature(iPoint) ) * ( Grad_Temp[0][0] * Grad_Temp[0][0] + Grad_Temp[1][0] * Grad_Temp[1][0] + Grad_Temp[2][0] * Grad_Temp[2][0] );
    
    return SIT;
  }

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
