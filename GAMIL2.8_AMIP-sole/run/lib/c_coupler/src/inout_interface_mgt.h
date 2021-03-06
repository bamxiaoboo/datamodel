/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef INOUT_INTERFACE_MGT_H
#define INOUT_INTERFACE_MGT_H


#include "restart_mgt.h"
#include <vector>
#include "common_utils.h"
#include "memory_mgt.h"
#include "timer_mgt.h"
#include "runtime_trans_algorithm.h"
#include "runtime_datatype_transformer.h"
#include "runtime_cumulate_average_algorithm.h"
#include "runtime_remap_algorithm.h"
#include "Runtime_Algorithm_Basis.h"


#define INTERFACE_TYPE_REGISTER             ((int)0)
#define INTERFACE_TYPE_IO_OUTPUT            ((int)1)
#define INTERFACE_TYPE_IO_WRITE             ((int)2)


class Inout_interface;
class Coupling_connection;


class Connection_field_time_info
{
	public:
		int current_year;
		int current_month;
		int current_day;
		int current_second;
		int current_num_elapsed_days;
		int last_timer_num_elapsed_days;
		int last_timer_second;
		int next_timer_num_elapsed_days;
		int next_timer_second;
		int time_step_in_second;
		int inst_or_aver;
		int lag_seconds;
		Coupling_timer *timer;
		Inout_interface *inout_interface;

		Connection_field_time_info(Inout_interface*, Coupling_timer*, int, int);
		void get_time_of_next_timer_on(bool);
		void reset_last_timer_info() { last_timer_num_elapsed_days = -1; last_timer_second = -1; }
		void write_into_array_for_restart(char**, long&, long&);		
		void import_restart_data(const char *, long &, const char *);
};


class Connection_coupling_procedure
{
	private:
		std::vector<Field_mem_info *> fields_mem_registered;
		std::vector<Field_mem_info *> fields_mem_inner_step_averaged;
		std::vector<Field_mem_info *> fields_mem_inter_step_averaged;
		std::vector<Field_mem_info *> fields_mem_remapped;
		std::vector<Field_mem_info *> fields_mem_datatype_transformed;
		std::vector<Field_mem_info *> fields_mem_unit_transformed;
		std::vector<Field_mem_info *> fields_mem_transfer;
		std::vector<int> field_interface_local_index;
		Connection_field_time_info * fields_time_info_src;
		Connection_field_time_info * fields_time_info_dst;
		long last_remote_fields_time;
		long current_remote_fields_time;
		Coupling_connection *coupling_connection;
		Inout_interface *inout_interface;
		std::vector<Runtime_cumulate_average_algorithm*> runtime_inner_averaging_algorithm;
		std::vector<Runtime_cumulate_average_algorithm*> runtime_inter_averaging_algorithm;
		std::vector<Runtime_remap_algorithm*> runtime_remap_algorithms;
		std::vector<Runtime_algorithm_basis*> runtime_unit_transform_algorithms;
		std::vector<Runtime_datatype_transformer*> runtime_datatype_transform_algorithms;
		Runtime_trans_algorithm *runtime_data_transfer_algorithm;
		bool finish_status;
		bool transfer_data;
		
	public:
		Connection_coupling_procedure(Inout_interface*, Coupling_connection*);
		~Connection_coupling_procedure();
		void add_data_transfer_algorithm(Runtime_trans_algorithm * runtime_algorithm) { runtime_data_transfer_algorithm = runtime_algorithm; }
		void execute(bool, int*, const char*);
		void send_fields(bool);
		Field_mem_info *get_data_transfer_field_instance(int); 
		int get_num_runtime_remap_algorithms() { return runtime_remap_algorithms.size(); }
		Runtime_remap_algorithm *get_runtime_remap_algorithm(int i) { return runtime_remap_algorithms[i]; }
		Runtime_remapping_weights *get_runtime_remapping_weights(int i);
		bool get_finish_status() { return finish_status; }
		void write_into_array_for_restart(char**, long&, long&);
		void import_restart_data(const char *, long &, const char *);
};


class Inout_interface
{
	private:
		char interface_name[NAME_STR_SIZE];
		char comp_full_name[NAME_STR_SIZE];
		int interface_id;
		int interface_type;
		int comp_id;
		int import_or_export_or_remap;     // 0: import; 1: export; 2: remap
		Time_mgt *time_mgr;
		Coupling_timer *timer;
		int inst_or_aver;
		std::vector<Field_mem_info *> fields_mem_registered;
		std::vector<const char*> fields_name;
		std::vector<Connection_coupling_procedure*> coupling_procedures;
		std::vector<Inout_interface *> children_interfaces;           // only for remap interface 
		int execution_checking_status;
		long last_execution_time;
		char fixed_remote_comp_full_name[NAME_STR_SIZE];
		char fixed_remote_interface_name[NAME_STR_SIZE];
		char *inversed_dst_fraction;
		long bypass_counter;

	public:
		Inout_interface(const char*, long&);
		Inout_interface(const char *, int, int, int *, int *, int, int, int, int, const char*, const char *);
		Inout_interface(const char*, int, int, int, int*, int, int, int, const char *, const char*, const char*, int, int);
		~Inout_interface();
		void initialize_data(const char *, int, int, int, int, int *, int, const char *);	
		void common_checking_for_interface_registration(int, int *, int, int, int, int, const char *, int, int, const char *, const char *);
		const char *get_interface_name() { return interface_name; }
		const char *get_comp_full_name() { return comp_full_name; }
		int get_comp_id() { return comp_id; }
		int get_interface_id() { return interface_id; }
		int get_interface_type() { return interface_type; }
		int get_import_or_export_or_remap() { return import_or_export_or_remap; }
		void report_common_field_instances(const Inout_interface*);
		void get_fields_name(std::vector<const char*>*);
		const char *get_field_name(int);
		int get_num_dst_fields();
		void transform_interface_into_array(char**, long&, long&);
		void write_into_array_for_restart(char**, long&, long&);
		void import_restart_data(const char *, long &, const char *);
		Field_mem_info *search_registered_field_instance(const char*, int &);
		Coupling_timer *get_timer() { return timer; }
		void add_coupling_procedure(Connection_coupling_procedure*);
		int get_inst_or_aver() { return inst_or_aver; } 
		void execute(bool, int*, int, const char*);
		Inout_interface *get_child_interface(int i);
		const char *get_fixed_remote_comp_full_name() { return fixed_remote_comp_full_name; }
		const char *get_fixed_remote_interface_name() { return fixed_remote_interface_name; }
		int get_num_coupling_procedures() { return coupling_procedures.size(); }
		void add_remappling_fraction_processing(void *, void *, int, int, const char *, const char *, const char *);		
		void preprocessing_for_frac_based_remapping();
		void postprocessing_for_frac_based_remapping(bool);
		long get_bypass_counter() { return bypass_counter; } 
};


class Inout_interface_mgt
{
	private:
		std::vector<Inout_interface*> interfaces;
		std::vector<Runtime_trans_algorithm*> all_runtime_receive_algorithms;
		char *temp_array_buffer;
		long buffer_max_size;
		long buffer_content_size;
		long buffer_content_iter;
		std::vector<MPI_Win> all_MPI_wins;

	public:
		Inout_interface_mgt(const char*, long);
		Inout_interface_mgt();
		~Inout_interface_mgt();
		int register_inout_interface(const char*, int, int, int*, int, int, int, const char*, const char*, int);
		void generate_remapping_interface_connection(Inout_interface *, int, int *, bool);
		int register_normal_remap_interface(const char *, int, int *, int *, int, int, int, int, const char *, const char *);
		int register_frac_based_remap_interface(const char *, int, int *, int *, int, int, int, int, void *, void *, int, int, const char *, const char *, const char *);
		int get_next_interface_id();
		bool is_interface_id_legal(int);
		Inout_interface *get_interface(int);
		Inout_interface *get_interface(const char*, const char*);
		Inout_interface *get_interface(int, const char*);
		void get_all_import_interfaces_of_a_component(std::vector<Inout_interface*>&, int);	
		void get_all_unconnected_fixed_interfaces(std::vector<Inout_interface*> &, int, int, const char*);
		void merge_unconnected_inout_interface_fields_info(int);
		void write_all_interfaces_fields_info();
		const char *get_temp_array_buffer() { return temp_array_buffer; } 
		int get_buffer_content_size()  { return buffer_content_size; }
		void execute_interface(int, bool, int*, int, int*, const char*);
		void execute_interface(int, const char*, bool, int *, int, int*, const char*);
		void add_runtime_receive_algorithm(Runtime_trans_algorithm *new_algorithm) { all_runtime_receive_algorithms.push_back(new_algorithm); }
		void erase_runtime_receive_algorithm(Runtime_trans_algorithm *);
		void runtime_receive_algorithms_receive_data();
		void add_MPI_win(MPI_Win mpi_win) { all_MPI_wins.push_back(mpi_win); }
		void free_all_MPI_wins(); 
		void write_into_restart_buffers(std::vector<Restart_buffer_container*> *, int);
		void import_restart_data(Restart_mgt *, const char *, int, const char *);
};

#endif
