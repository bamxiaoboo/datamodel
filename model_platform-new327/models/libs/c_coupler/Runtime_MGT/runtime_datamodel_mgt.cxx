/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <dirent.h>
#include <stdlib.h>
#include "runtime_datamodel_mgt.h"
#include "memory_mgt.h"
#include "field_info_mgt.h"

void *load_datamodel_instatnces_configuration(int comp_id, const char *comp_full_name, std::vector<std::pair<const char*, const char*> > &producers_info)
{
    char XML_file_name[NAME_STR_SIZE];
    int line_number;

    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_full_name);
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, MPI_COMM_NULL, false);
    if (XML_file == NULL)
        return;

    TiXmlElement *root_XML_element, *Data_inst_XML_element;
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file), *Data_inst_XML_element_node = NULL;
    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {
        if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        root_XML_element = root_XML_element_node->ToElement();
        if (words_are_the_same(root_XML_element->Value(),"datamodel_instances"))
            break;

    }
    if (root_XML_element_node == NULL)
        return;

    for (XML_element_node = root_XML_element->FirstChilde(); XML_element_node != NULL; XML_element_node = XML_element_node->NextSibling()) {
        data_inst_XML_element = XML_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_inst_XML_element->Value(),"datamodel_instance"), "The XML element for specifying the configuration information of a Datamodel instance in the XML configuration file \"%s\" should be named \"datamodel_instance\". Please verify the XML file around the line number %d.", XML_file_name, data_inst_XML_element->Row());
        const char *datamodel_instance_name = get_XML_attribute(comp_id, 80, data_inst_XML_element, "name", XML_file_name, line_number, "the \"name\" of an Datamodel instance", "Datamodel instance configuration file", true);//should be false?
        if (!is_XML_setting_on(comp_id, data_inst_XML_element, XML_file_name, "the \"status\" of the Datamodel instance for an component model", "Datamodel instance configuration file"))
            continue;
        check_and_verify_name_format_of_string_for_XML(-1, datamodel_instance_name, "datamodel instance", XML_file_name, line_number);//-1?
        Datamodel_instances_mgr->config_a_datamodel_instance(comp_id, comp_full_name, XML_file_name, *datamodel_instance_name, *XML_element_node, producers_info);
    }
    delete XML_file;
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish loading the configuration of datamodel instance from the XML file %s", XML_file_name);
}


Datamodel_instances_mgt::Datamodel_instances_mgt(int comp_id, char *datamodel_instance_name, TiXmlNode *data_inst_XML_element, producers_info){
    int line_number;
    TiXmlNode *datamodel_name_element_node = data_inst_XML_element->FirstChild();
    if (inst_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
        return;//??
    for (TiXmlNode *time_mapping_element_node = data_inst_XML_element->NextSibling(); time_mapping_element_node != NULL; time_mapping_element_node = time_mapping_element_node->NextSilbling()) {
        if (time_mapping_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        TiXmlElement *time_mapping_element = time_mapping_element_node->ToElement();
        EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(time_mapping_element->Value(), "time_mapping_configuration"), "When setting the datamodel instance configuration of the component \"%s\" in the XML file \"%s\", the XML element for specifying the datamodel instance configuration should be named \"time_mapping_configuration\". Pleas verify the XML file arround the line number %d.", comp_full_name, XML_file_name, time_mapping_element->Row());
        if (!is_XML_setting_on(comp_id, time_mapping_element, XML_file_name, "the status of some time mapping configurations for a datamodel instance", "datamodel instance configuration file"))
            continue;
        const char *specification_name = get_XML_attribute(comp_id, 80, time_mapping_element, "specification", XML_file_name, line_number, "the way of \"specification\" of a datamodel instance", "datamodel instance configuration file", true);
        if (specification_name == "period"){
            specification = 0;
            for (TiXmlNode *period_setting_element_node = time_mapping_element->FirstChild(); period_setting_element_node != NULL; period_setting_element_Node = period_setting_element_node->NextSibling()) {
                TiXmlElement *period_setting_element = period_setting_element_node->ToElement();

                if (!is_XML_setting_on(comp_id, period_setting_element, XML_file_name, "the status of period setting configuration for a datamodel instance", "datamodel instance configuration file"))
                    continue;
                else
                    EXECUTION_REPORT(ERROR, comp_id, !(is_XML_setting_on(comp_id, period_setting_element, XML_file_name, "the status of period setting configuration for a datamodel instance", "datamodel instance configuration file") && is_period_set), "When setting the datamodel instance configuration of the component \"%s\" in the XML file \"%s\", the XML element for period setting should only be set once. Pleas verify the XML file arround the line number %d.", comp_full_name, XML_file_name, period_setting_element->Row());
                is_period_set = 1;
                char *period_unit_config = get_XML_attribue(comp_id, 80, period_setting_element, "period_unit", XML_file_name, period_setting_element->Row(), "period unit of a datamodel instance", "datamodel instance configuration file", true);
                Period_setting.period_unit = check_unit_format(*period_unit_config);
                char *period_start_time_config = get_XML_attribute(comp_id, 80, period_setting_element, "period_start_time", XML_file_name, period_setting_element->Row(), "period start time of a datamodel instance", "datamodel instance configuration file", true);
                Period_setting.period_start_time = check_time_format(period_start_time_config, Period_setting.period_unit);
                char *period_count_config = get_XML_attribute(comp_id, 80, period_setting_element, "period_count", XML_file_name, period_setting_element->Row(), "period count of a datamodel instance", "datamodel instance configuration file", true);
                Period_setting.period_count = atoi(period_count_config);
            }
        }
        else if (specification_name == "offset") {
            specification = 1;
            //for ()
        }
        else if (specification_name == "default")
            specification = 2;
        else
            EXECUTION_REPORT(REPORT_ERROR, comp_id, false,"The specification for time mapping configuration of component \"%s\" in the XML file \"%s\" should only be \"default\", \"period\" or \"offset\". Please verify the XML file arround the line number %d.", comp_full_name, XML_file_name, time_mapping_element->Row());
    }
}


void *Datamodel_instances_mgt::config_a_datamodel_instance(int comp_id, char *comp_full_name, char *XML_file_name, char *Datamodel_instance_name, TiXmlNode *data_inst_XML_element, producers_info){
    for (int i=0; i < datamodel_instances_mgr.size(); i++)
        if (words_are_the_same(Datamodel_instance_name, datamodel_instances_mgr[i]->get_datamodel_instance_name()))
            return;
    datamodel_instances_mgr.push_back(new Datamodel_instance(comp_id, comp_full_name, XML_file_name, Datamodel_instance_name, data_inst_XML_element, producers_info));//parameter comp_id may not be needed,instance belongs to the component model that uses it,so it uses the comp id of the component which it belongs to
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish configuring datamodel instance %s",Datamodel_instance_name);

    return; 
}

int Datamodel_instance::check_time_format(const char *Period_start_time_config, int Period_start_time) {
}


int Datamodel_instance::check_unit_format(const char *Period_unit_config) {
}
