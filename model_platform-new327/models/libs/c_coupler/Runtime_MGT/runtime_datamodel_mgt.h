/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef RUNTIME_DATAMODEL_MGT
#define RUNTIME_DATAMODEL_MGT

#include "global_data.h"
#include "memory_mgt.h"
#include "timer_mgt.h"
#include "common_utils.h"
#include <string.h>

struct Datamodel_field_info {
	//char cfg_info_remain_line[NAME_STR_SIZE];
	char field_name_in_coupling[NAME_STR_SIZE];
    char field_name_in_IO_file[NAME_STR_SIZE];
    char field_datatype_in_IO_file[NAME_STR_SIZE];
    Field_mem_info *field_data_mem;
    //double add_offset;//????
    char field_unit[NAME_STR_SIZE];
    char field_unit_in_IO_file[NAME_STR_SIZE];
    //bool have_scale_factor;
    //double scale_factor;//????
};

struct Time_location_in_files//???
{
	int date;
	int num_elapsed_day;
	int second_in_day;
	char file_name[NAME_STR_SIZE];
	int time_dim_num;
};


/*struct Time_field
{
	char time_field_name[16][NAME_STR_SIZE];//xml, time_field
	char time_field_format[16][NAME_STR_SIZE];//xml, time_field
	int time_field_format_id[16];//xml, time_field[16]
};*/


struct Horizontal_grid
{
	char grid_name[NAME_STR_SIZE];
	bool specification;//0 for file name, 1 for file field
	char grid_file_name[NAME_STR_SIZE];
	char grid_file_dir[NAME_STR_SIZE];
	char edge_type[NAME_STR_SIZE];
	char coord_unit_name[NAME_STR_SIZE];
	bool cyclic_or_acyclic;//0 for cyclic, 1 for acyclic
	int dim_size1;
	int dim_size2;
	char dim_size1_name[NAME_STR_SIZE];
	char dim_size2_name[NAME_STR_SIZE];
	char min_lon_field_name[NAME_STR_SIZE];
	char max_lon_field_name[NAME_STR_SIZE];
	char min_lat_field_name[NAME_STR_SIZE];
	char max_lat_field_name[NAME_STR_SIZE];
	char center_lon_field_name[NAME_STR_SIZE];
	char center_lat_field_name[NAME_STR_SIZE];
	char mask_field_name[NAME_STR_SIZE];
	char area_field_name[NAME_STR_SIZE];
	char vertex_lon_field_name[NAME_STR_SIZE];
	char vertex_lat_field_name[NAME_STR_SIZE];
	char horizontal_grid_annotation[256];
};

struct z_grid
{
	char z_coord_values_field_name[NAME_STR_SIZE];
	char z_grid_annotation[256];
};

struct sigma_grid
{
	char s_coord_values_field_name[NAME_STR_SIZE];
	char sigma_grid_annotation[256];
	char top_value_field_name[NAME_STR_SIZE];
	char sigma_values_field_name[NAME_STR_SIZE];
};

struct hybrid_grid
{
	char hybrid_coord_values_field_name[NAME_STR_SIZE];
	char top_value_field_name[NAME_STR_SIZE];
	char coef_a_field_name[NAME_STR_SIZE];
	char coef_b_field_name[NAME_STR_SIZE];
	char hybrid_grid_annotation[256];
};

struct Vertical_grid
{
	char grid_name[NAME_STR_SIZE];
	char grid_type[16];
	char coord_unit_name[NAME_STR_SIZE];
	union vertical_type_grid
	{
		struct sigma_grid;
		struct hybrid_grid;
		struct z_grid;
	};
};

struct 3D_grid
{
	char horizontal_sub_grid_name[NAME_STR_SIZE];
	char vertical_sub_grid_name[NAME_STR_SIZE];
	char mid_point_grid_name[NAME_STR_SIZE];
	char surface_field_type[NAME_STR_SIZE];
	char surface_field_variable[NAME_STR_SIZE];
	bool generate_mid_point_grid;//0 for not to generate, 1 for to generate
};

class Datamodel_field_read_handler
{
	struct Field_read_info {
		int last_pos_time_interval_left;
		int last_pos_time_interval_right;
		long last_read_full_time;
		int reference_count;
		Field_mem_info *temp_field_before_time_remap_left;
		Field_mem_info *temp_field_before_time_remap_right;
		Field_mem_info *temp_field_readin;//
		Field_mem_info *output_field_instance;//
		Runtime_remap_algorithm *remap_algorithm_for_readin_left;
		Runtime_remap_algorithm *remap_algorithm_for_readin_right;
	};
	
	private:
		/* Information for file reading */
		char handler_name[NAME_STR_SIZE];

		int id_time_format_in_file_name;
		char file_dir[NAME_STR_SIZE];//xml, file name
		char file_name_pattern[NAME_STR_SIZE];//xml, file name
		char time_format_type_in_file_name[NAME_STR_SIZE];//xml, file name

		std::vector<Time_location_in_files> time_filename_map; //xml, time_field

		int num_time_field;//xml, time_field
		char time_field_name[16][NAME_STR_SIZE];//
		char time_field_format[16][NAME_STR_SIZE];//
		int time_field_format_id[16];//xml, time_field


		int pos_time_interval_left;
		int pos_time_interval_right;
		int time_point_type; //1: start, 2: middle, 3: end. Default is 1

		void initialize_time_filename_map();//file_dir???
		bool get_time_from_string(const char*, const char*, int, int, bool, int&, int&);
		int check_time_format(const char*, const char*);
		
	public:	
		Datamodel_field_read_handler(const char*);
		long search_current_time_interval(int);
		const char *get_handler_name() { return handler_name; }
		void update_one_field(int);
		int register_read_info_for_a_field(Field_mem_info*, Remap_weight_of_strategy_class*);
};

class Datamodel_field_read_handler_mgt
{
    private:
		std::vector<Datamodel_field_read_handler*> datamodel_field_read_handlers;

	public:
		Datamodel_field_read_handler_mgt() {}
		~Datamodel_field_read_handler_mgt();
		Datamodel_field_read_handler *get_a_handler(const char*);
};


class Runtime_datamodel
{
    private:
        std::vector<Datamodel_field_info*> datamodel_fields;
		std::vector<int> field_read_info_indexes;
		
        char datamodel_name[NAME_STR_SIZE];
        char IO_file_name[NAME_STR_SIZE];
        IO_netcdf *netcdf_file_object;
        //Coupling_timer *io_timer;//needed
		Datamodel_field_read_handler *field_read_handler;
		long last_read_full_time;
		std::vector<Horizontal_grid> horizontal_grid_info;
		std::vector<Vertical_grid> vertical_grid_info;
		std::vector<3D_grid> 3d_grid_info;
		
        void datamodel_read(void);
        void datamodel_check(void);
		void generate_algorithm_info_from_cfg_file();
		void allocate_one_field(int);

    public:
        Runtime_datamodel_algorithm(const char *);
        ~Runtime_datamodel_algorithm();
        void change_IO_file_name_for_restart(const char *);//??
        virtual void run(bool);//???
        void allocate_src_dst_fields(bool);
		const char *get_algorithm_cfg_name() { return algorithm_cfg_name; }
};

class Datamodel_mgt
{
    private:
        std::vector<Runtime_datamodel*> Runtime_datamodels;
    public:
        Datamodel_mgt() {};
        ~Datamodel_mgt();
};

class Datamodel_instance
{
private:
    char offset_unit[NAME_STR_SIZE];
    int offset_count;

	int period;// 0: acyclic; 1: year; 2: month; 3: day
    int period_start_time;
    char period_unit[NAME_STR_SIZE];
    int period_count;

	char datamodel_instance_name[NAME_STR_SIZE];
    char datamodel_name[NAME_STR_SIZE];
	//Runtime_datamodel_algorithm* Runtime_datamodel;
public:
    Datamodel_instance() {};
    ~Datamodel_instance();
};

class Datamodel_instances_mgt
{
private:
	std::vector<Datamodel_instance*> Runtime_datamodel_instances;
public:
	Datamodel_instances_mgt() {};
	~Datamodel_instances_mgt();
	Runtime_datamodel_instances *get_a_datamodel_instance(int, const char*, TiXmlNode*, &producers_info);//whether to return?the type
    void *load_datamodel_instatnces_configuration(int, const char*);
};

#endif
