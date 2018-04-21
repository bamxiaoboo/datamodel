/***************************************************************
  This is a source file of GAMIL, which registers all model 
  subroutines of exchanging surface data into C-Coupler library. 
  This file was initially finished by Dr. Li Liu. If you have any 
  problem, please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
 ***************************************************************/


#include "coupling_interface.h"


extern "C" void initialize_online_lsm_ccpl_interface_();
extern "C" void calculate_ocn_albedo_ccpl_interface_();
extern "C" void calculate_sice_albedo_ccpl_interface_();
extern "C" void initialize_data_ocn_ccpl_interface_();
extern "C" void initialize_data_sice_ccpl_interface_();
extern "C" void merge_ts_ccpl_interface_();
extern "C" void get_initial_surface_data_from_coupler_ccpl_interface_();
extern "C" void get_forcing_data_sst_ccpl_interface_();
extern "C" void get_forcing_data_sice_ccpl_interface_();
extern "C" void srfflx_state_reset_ccpl_interface_();
extern "C" void update_srf_fractions_ccpl_interface_();
extern "C" void process_surface_data_with_online_lsm_ccpl_interface_();
extern "C" void process_surface_data_with_forcing_data_sst_ccpl_interface_();
extern "C" void process_surface_data_with_forcing_data_sice_ccpl_interface_();
extern "C" void process_surface_data_with_coupler_ccpl_interface_();  
extern "C" void adjust_land_ocn_sice_fraction_ccpl_interface_();


extern "C" void register_gamil_surface_subroutines_()
{
   C_Coupler_interface_register_model_algorithm("initialize_online_lsm", initialize_online_lsm_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("calculate_ocn_albedo", calculate_ocn_albedo_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("calculate_sice_albedo", calculate_sice_albedo_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("initialize_data_ocn", initialize_data_ocn_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("initialize_data_sice", initialize_data_sice_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("merge_ts", merge_ts_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("get_initial_surface_data_from_coupler", get_initial_surface_data_from_coupler_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("get_forcing_data_sst", get_forcing_data_sst_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("get_forcing_data_sice", get_forcing_data_sice_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("srfflx_state_reset", srfflx_state_reset_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("update_srf_fractions", update_srf_fractions_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("process_surface_data_with_online_lsm", process_surface_data_with_online_lsm_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("process_surface_data_with_forcing_data_sst", process_surface_data_with_forcing_data_sst_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("process_surface_data_with_forcing_data_sice", process_surface_data_with_forcing_data_sice_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("process_surface_data_with_coupler", process_surface_data_with_coupler_ccpl_interface_);
   C_Coupler_interface_register_model_algorithm("adjust_land_ocn_sice_fraction", adjust_land_ocn_sice_fraction_ccpl_interface_);
}
