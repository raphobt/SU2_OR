extern "C" 
{
//void __properties_MOD_pressure(double* ii,double* jj, double* kk);
void __properties_MOD_dpdr_u(double* TT, su2double* vv, double* deriv);
void __properties_MOD_dpdu_v(double* T_in, su2double* v_in, double* deriv);
void __properties_MOD_dpdv_t(su2double* temp,su2double* vol, double* deriv);
void __properties_MOD_dpdt_v(su2double* temp,su2double* vol, double* deriv);
void __saturation_MOD_sat_curve();
void __saturation_MOD_saturation_curve();
void __properties_MOD_satprop(int* mode, su2double* psat, su2double* Tsat, su2double* vvsat, su2double* vlsat, su2double* uvsat, su2double* ulsat);
void __properties_MOD_inter_energy(su2double* temp,su2double* vol, su2double* energy);
void __properties_MOD_entropy(double* T_in, su2double* v_in,double* s_out);
void __non_linear_solvers_MOD_new_rap1d(int* MODE, su2double* out_1, su2double* out_2, double* resnorm,
                                        int* Niter, int* exitflag, su2double* GGG, su2double* X0, su2double* in_1, double* out3);
void __non_linear_solvers_MOD_new_rap2d(int* MODE, su2double* out1_NewRap, su2double* out2_NewRap, double* resnorm, int* Niter,
					  int* exitflag,su2double* in_1, su2double* in_2, double* guess_1, double* guess_2);
void __transprop_MOD_entropyco2(double* s_out, su2double* v_in, double* vp_out, double* x_in, double* T_in, double* p_in, int* flag_loca);
void __grid_MOD_make_grid();
void __interp_table_MOD_co2bllt_equi(double* p_out, double* T_out,double* c_out,
                                     double* x_out, double* a_out,double* dummy,
                                     su2double* u_in,  su2double* v_in, int*    flag);
void __solver_eos_MOD_brentroots2(int* MODE, double* out_1,double* out_2,double* residue,
                                  int* Niter,double* GGG,double* lower,double* upper,double* in_1,double* in_2);

void __solver_eos_MOD_eos_1d(int* MODE,double* es,double* out_2,double* residue,int* Niter,int* flag,double* Qp5,double* eguess,double*vbc,double* in_2);
void __derivees_MOD_co2der(double* dp_dv_u,double* dp_du_v ,su2double* e0, su2double* v0,double* T0,double* p0,double* res2,double* res3,double* res4);
void __transprop_MOD_co2visco(double* vis,su2double* vv,su2double* TT,su2double* x_out,su2double* pp,int* flag);
void __transprop_MOD_cpco2(double* cp_out,su2double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_cvco2(double* cv_out,su2double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_dedrco2(double* dedr_out,su2double* v_in,double* vp_out,double* x_in,double* T_in,double* p_in,int* flag_loca);
void __transprop_MOD_co2conduc2phase(double* lambda_out,su2double* v_in,double* vp_out,su2double* x_in,su2double* T_in,su2double* p_in,int* flag_loca);
//void __properties_MOD_heat_cap_p(double* temp,double* vol, double* cp)
//void __properties_MOD_heat_cap_v(double* temp,double* vol, double* cv)
//void grid_construction_left_low_();
//void grid_construction_left_high_();
//void grid_construction_high_temperature_();
}
//extern double* __def_variables_MOD_ccc_ll;


