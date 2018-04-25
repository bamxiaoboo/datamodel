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

void *Datamodel_instances_mgt::load_datamodel_instatnces_configuration(int comp_id, const char *comp_full_name, std::vector<std::pair<const char*, const char*> > &producers_info)
{
    char XML_file_name[NAME_STR_SIZE];
    int line_number;

    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_full_name);
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, MPI_COMM_NULL, false);
    if (XML_file == NULL)
        return;

    TiXmlElement *root_XML_element, *Data_inst_XML_element, *detailed_XML_element;
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file), *Data_inst_XML_element_node = NULL, *detailed_XML_element_node = NULL;
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
        Datamodel_instances_mgr->get_a_datamodel_instance(comp_id, *datamodel_instance_name, *data_inst_XML_element, producers_info);//??
    }
    delete XML_file;
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish loading the configuration of datamodel instance from the XML file %s", XML_file_name);
}


Datamodel_instances_mgt::Datamodel_instances_mgt(){
}


void *Datamodel_instances_mgt::get_a_datamodel_instance(int comp_id, char *Datamodel_instance_name, TiXmlNode *data_inst_XML_element, producers_info){
    for (int i=0; i < datamodel_instances_mgr.size(); i++)
        if (words_are_the_same(Datamodel_instance_name, datamodel_instances_mgr->get_datamodel_instance_name()))
            return datamodel_instances_mgr[i];
    datamodel_instances_mgr.push_back(new Datamodel_instance());

    return datamodel_instances_mgr[datamodel_instances_mgr.size()-1];
}
