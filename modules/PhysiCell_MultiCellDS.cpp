/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
 
#include "PhysiCell_MultiCellDS.h"
// #include "../core/PhysiCell_cell.h"    //rwh
#ifdef ADDON_PHYSIBOSS
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"	
#endif
namespace PhysiCell{

void add_PhysiCell_cell_to_open_xml_pugi(  pugi::xml_document& xml_dom, Cell& C ); // not implemented -- future edition 

void add_PhysiCell_cells_to_open_xml_pugi( pugi::xml_document& xml_dom, std::string filename_base, Microenvironment& M  )
{
	std::cout << "Warning: " << __FUNCTION__ << " is deprecated and has been removed." << std::endl; 
		
	return; 
}

void add_PhysiCell_to_open_xml_pugi( pugi::xml_document& xml_dom , std::string filename_base, double current_simulation_time , Microenvironment& M );

void save_PhysiCell_to_MultiCellDS_xml_pugi( std::string filename_base , Microenvironment& M , double current_simulation_time)
{
	std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; 

	// start with a standard BioFVM save
	
	add_BioFVM_to_open_xml_pugi( BioFVM::biofvm_doc , filename_base , current_simulation_time , M ); 
	
	// now, add the PhysiCell data 

	add_PhysiCell_cells_to_open_xml_pugi( BioFVM::biofvm_doc , filename_base , M  ); 
	// add_PhysiCell_cells_to_open_xml_pugi_v2( BioFVM::biofvm_doc , filename_base , M  ); 
		
	// Lastly, save to the indicated filename 

	char filename[1024]; 
	sprintf( filename , "%s.xml" , filename_base.c_str() ); 
	BioFVM::biofvm_doc.save_file( filename );

	return; 
}


void save_PhysiCell_to_MultiCellDS_v2( std::string filename_base , Microenvironment& M , double current_simulation_time)
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this one July 2024

	// set some metadata

	BioFVM::MultiCellDS_version_string = "2"; 
	BioFVM::BioFVM_metadata.program.program_name = "PhysiCell"; 
	BioFVM::BioFVM_metadata.program.program_version = PhysiCell_Version; 
	BioFVM::BioFVM_metadata.program.program_URL = "http://physicell.org"; 

	BioFVM::BioFVM_metadata.program.creator.type = "creator"; 	
	BioFVM::BioFVM_metadata.program.creator.surname = "Macklin"; 	
	BioFVM::BioFVM_metadata.program.creator.given_names = "Paul"; 	
	BioFVM::BioFVM_metadata.program.creator.email = "macklinp@iu.edu"; 	
	BioFVM::BioFVM_metadata.program.creator.URL = "http://MathCancer.org"; 	
	BioFVM::BioFVM_metadata.program.creator.organization = "Indiana University & PhysiCell Project"; 	
	BioFVM::BioFVM_metadata.program.creator.department = "Intelligent Systems Engineering"; 	
	BioFVM::BioFVM_metadata.program.creator.ORCID = "0000-0002-9925-0151"; 	
	
	BioFVM::BioFVM_metadata.program.citation.DOI = "10.1371/journal.pcbi.1005991"; 
	BioFVM::BioFVM_metadata.program.citation.PMID = "29474446"; 
	BioFVM::BioFVM_metadata.program.citation.PMCID = "PMC5841829"; 
	BioFVM::BioFVM_metadata.program.citation.text = "A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: 10.1371/journal.pcbi.1005991"; 
	BioFVM::BioFVM_metadata.program.citation.notes = ""; 
	BioFVM::BioFVM_metadata.program.citation.URL = "https://dx.doi.org/PMC5841829"; 

	// start with a standard BioFVM save
		// overall XML structure 
	add_MultiCellDS_main_structure_to_open_xml_pugi( BioFVM::biofvm_doc ); 
		// save metadata 
	BioFVM_metadata.add_to_open_xml_pugi( current_simulation_time , BioFVM::biofvm_doc ); 
		// save diffusing substrates 
	add_BioFVM_substrates_to_open_xml_pugi( BioFVM::biofvm_doc , filename_base, M  ); 

		// add_BioFVM_agents_to_open_xml_pugi( xml_dom , filename_base, M); 
	
	// now, add the PhysiCell data 

	// add_PhysiCell_cells_to_open_xml_pugi( BioFVM::biofvm_doc , filename_base , M  ); 
	add_PhysiCell_cells_to_open_xml_pugi_v2( BioFVM::biofvm_doc , filename_base , M  ); 

	// Lastly, save to the indicated filename 

	char filename[1024]; 
	sprintf( filename , "%s.xml" , filename_base.c_str() ); 
	BioFVM::biofvm_doc.save_file( filename );

	return; 
}

/* look here */

int total_data_size( std::vector<int>& data_sizes )
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this one July 2024

	// current index: 
	int total = 0; 
	for( int i = 0 ; i < data_sizes.size(); i++ )
	{ total += data_sizes[i]; }
	return total; 
}

void add_variable_to_labels( std::vector<std::string>& data_names , 
	std::vector<std::string>& data_units, 
	std::vector<int>& data_start_indices, 
	std::vector<int>& data_sizes, 
	std::string var_name, std::string var_units, int var_size )
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this one July 2024

	// current index: 
	int index = 0; 
	for( int i = 0 ; i < data_sizes.size(); i++ )
	{ index += data_sizes[i]; }

	data_names.push_back( var_name ); 
	data_units.push_back( var_units ); 
	data_sizes.push_back( var_size ); 
	data_start_indices.push_back( index ); 

	return; 
}
void dump_cell_mat_vars(int number_of_data_entries )
{
    Cell* pCell; 
    double dTemp;

    double *dp; 
    std::cout << "\n---------- START dump_cell_mat_vars:\n";
	for( int i=0; i < number_of_data_entries ; i++ )
	{
		pCell = (*all_cells)[i]; 

		int writes = 0; 

		// compatibilty : first 17 entries 
		// ID 					<label index="0" size="1">ID</label>
		// double ID_temp = (double) (*all_cells)[i]->ID;
		// fwrite( (char*) ID_temp << std::endl; 

		// name = "ID"; 
		dTemp = (double) pCell->ID;
		std::cout << "ID= " << dTemp << std::endl; 
		// name = "position";    NOTE very different syntax for writing vectors!
        // std::cout <<  "pos= " << pCell->position.data()  << std::endl;
        dp = pCell->position.data();
		std::cout <<  "pos= " << dp[0]<<", "<<dp[1]<<", "<<dp[2]  << std::endl; 
		// name = "total_volume"; 
		std::cout <<  "vol= " << pCell->phenotype.volume.total << std::endl; 
		// name = "cell_type"; 
		dTemp = (double) pCell->type;
		std::cout <<  "type= " << dTemp << std::endl; 
		// name = "cycle_model"; 
		dTemp = (double) pCell->phenotype.cycle.model().code; 
		std::cout <<  "cycle model= " << dTemp << std::endl; // cycle model 
		// name = "current_phase"; 
		dTemp = (double) pCell->phenotype.cycle.current_phase().code; 
		std::cout <<  "cycle phase= " << dTemp << std::endl; // cycle model 
		// name = "elapsed_time_in_phase"; 
		std::cout <<  "elapsed time= " << pCell->phenotype.cycle.data.elapsed_time_in_phase << std::endl; 
		// name = "nuclear_volume"; 
		std::cout <<  "volume amts: " << pCell->phenotype.volume.nuclear << std::endl;   
		// name = "cytoplasmic_volume"; 
		std::cout <<  pCell->phenotype.volume.cytoplasmic << std::endl;
		// name = "fluid_fraction"; 
		std::cout <<  pCell->phenotype.volume.fluid_fraction << std::endl;
		// name = "calcified_fraction"; 
		std::cout <<  pCell->phenotype.volume.calcified_fraction << std::endl; 
		// name = "orientation"; 
		// std::cout <<  "orientation= " << pCell->state.orientation.data()  << std::endl; 
		double *dp = pCell->state.orientation.data(); 
		std::cout <<  "orientation= " << dp[0]<<", "<<dp[1]<<", "<<dp[2]  << std::endl; 
		// name = "polarity"; 
		std::cout <<  "polarity= " << pCell->phenotype.geometry.polarity << std::endl; 

 /* state variables to save */ 
// state
        std::cout << "------ state:\n";
		// name = "velocity"; 
		// std::cout <<  pCell->velocity.data()  << std::endl; 
		dp = pCell->velocity.data() ; 
		std::cout << " velocity= " << dp[0]<<", "<<dp[1]<<", "<<dp[2]  << std::endl; 
		// name = "pressure"; 
		std::cout <<  pCell->state.simple_pressure << std::endl; 
		// name = "number_of_nuclei"; 
		dTemp = (double) pCell->state.number_of_nuclei; 
		std::cout <<  dTemp << std::endl; 
		// // name = "damage"; 
		// std::cout <<  pCell->phenotype.integrity.damage << std::endl; 
		// name = "total_attack_time"; 
		std::cout <<  pCell->state.total_attack_time << std::endl; 
		// name = "contact_with_basement_membrane"; 
		dTemp = (double) pCell->state.contact_with_basement_membrane; 
		std::cout <<  dTemp << std::endl; 

/* now go through phenotype and state */ 
// cycle 
  // current exit rate // 1 
        std::cout << "------ cycle:\n";
		// name = "current_cycle_phase_exit_rate"; 
		int phase_index = pCell->phenotype.cycle.data.current_phase_index; 
		std::cout <<  "current phase exit rate= " <<pCell->phenotype.cycle.data.exit_rate(phase_index) << std::endl; 
		// name = "elapsed_time_in_phase"; 
		std::cout << "elapsed time= " << pCell->phenotype.cycle.data.elapsed_time_in_phase << std::endl; 

// death 
  // live or dead state // 1 
        std::cout << "------ death:\n";
		// name = "dead"; 
		dTemp = (double) pCell->phenotype.death.dead; 
		std::cout << " dead=" << dTemp << std::endl; 
		// name = "current_death_model"; // 
		dTemp = (double) pCell->phenotype.death.current_death_model_index; 
		std::cout <<  " current_death_model_index=" << dTemp << std::endl; 
		// name = "death_rates"; 
		// std::cout <<  pCell->phenotype.death.rates.data() << std::endl; //, sizeof(double) , nd , fp ); 
		dp = pCell->phenotype.death.rates.data(); 
        int nd = pCell->phenotype.death.rates.size(); 
		std::cout << " phenotype.death.rates:"<< std::endl; //, sizeof(double) , nd , fp ); 
        for (int idx=0; idx<nd; idx++)
		    std::cout <<  "  "<< dp[idx] << std::endl;
		
	// volume ()
        std::cout << "------ volume:\n";
		// name = "cytoplasmic_biomass_change_rate"; 
		std::cout <<  pCell->phenotype.volume.cytoplasmic_biomass_change_rate << std::endl; 
		// name = "nuclear_biomass_change_rate"; 
		std::cout <<  pCell->phenotype.volume.nuclear_biomass_change_rate << std::endl; 
		// name = "fluid_change_rate"; 
		std::cout <<  pCell->phenotype.volume.fluid_change_rate << std::endl; 
		// name = "calcification_rate"; 
		std::cout <<  pCell->phenotype.volume.calcification_rate << std::endl; 
		// name = "target_solid_cytoplasmic"; 
		std::cout <<  pCell->phenotype.volume.target_solid_cytoplasmic << std::endl; 
		// name = "target_solid_nuclear"; 
		std::cout <<  pCell->phenotype.volume.target_solid_nuclear << std::endl; 
		// name = "target_fluid_fraction"; 
		std::cout <<  pCell->phenotype.volume.target_fluid_fraction << std::endl; 

  // geometry 
        std::cout << "------ geometry:\n";
     // radius //1 
		// name = "radius"; 
		std::cout <<  pCell->phenotype.geometry.radius << std::endl; 
		// name = "nuclear_radius"; 
		std::cout <<  pCell->phenotype.geometry.nuclear_radius << std::endl; 
		// name = "surface_area"; 
		std::cout <<  pCell->phenotype.geometry.surface_area << std::endl; 

  // mechanics 
        std::cout << "------ mechanics:\n";
	// cell_cell_adhesion_strength; // 1
		// name = "cell_cell_adhesion_strength"; 
		std::cout <<  pCell->phenotype.mechanics.cell_cell_adhesion_strength << std::endl; 
		// name = "cell_BM_adhesion_strength"; 
		std::cout <<  pCell->phenotype.mechanics.cell_BM_adhesion_strength << std::endl; 
		// name = "cell_cell_repulsion_strength"; 
		std::cout <<  pCell->phenotype.mechanics.cell_cell_repulsion_strength << std::endl; 
		// name = "cell_BM_repulsion_strength"; 
		std::cout <<  pCell->phenotype.mechanics.cell_BM_repulsion_strength << std::endl; 
		// name = "cell_adhesion_affinities"; 
		std::cout <<"cell_adhesion_affinities\n"; 
		dp = pCell->phenotype.mechanics.cell_adhesion_affinities.data();
		int ncelltypes = pCell->phenotype.mechanics.cell_adhesion_affinities.size();
        for (int idx=0; idx<ncelltypes; idx++)
		    std::cout <<  "   " <<dp[idx] << std::endl;  // , sizeof(double) , n , fp ); 
		// name = "relative_maximum_adhesion_distance"; 
		std::cout <<  pCell->phenotype.mechanics.relative_maximum_adhesion_distance << std::endl; 
		// name = "maximum_number_of_attachments"; 
		dTemp = (double) pCell->phenotype.mechanics.maximum_number_of_attachments; 
		std::cout <<  dTemp << std::endl; 
		// name = "attachment_elastic_constant"; 
		std::cout <<  pCell->phenotype.mechanics.attachment_elastic_constant << std::endl; 
		// name = "attachment_rate"; 
		std::cout <<  pCell->phenotype.mechanics.attachment_rate << std::endl; 
 		// name = "detachment_rate"; 
		std::cout <<  pCell->phenotype.mechanics.detachment_rate << std::endl; 

 // Motility
        std::cout << "------ motility:\n";
 		// name = "is_motile"; 
		dTemp = (double) pCell->phenotype.motility.is_motile; 
		std::cout <<  dTemp << std::endl; 
 		// name = "persistence_time"; 
		std::cout <<  pCell->phenotype.motility.persistence_time << std::endl; 
 		// name = "migration_speed"; 
		std::cout <<  pCell->phenotype.motility.migration_speed << std::endl; 
 		// name = "migration_bias_direction"; 
 		std::cout << "migration_bias_direction\n"; 
		dp = pCell->phenotype.motility.migration_bias_direction.data(); 
		std::cout <<  dp[0]<< ", "<< dp[1]<< ", "<<dp[2]  << std::endl; 
 		// name = "migration_bias"; 
		std::cout <<  pCell->phenotype.motility.migration_bias << std::endl; 
 		// name = "motility_vector"; 
 		std::cout << "motility_vector\n"; 
		dp = pCell->phenotype.motility.motility_vector.data(); 
		std::cout <<  dp[0]<< ", "<< dp[1]<< ", "<<dp[2]  << std::endl; 
 		// name = "chemotaxis_index"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_index; 
		std::cout <<  dTemp << std::endl; 
 		// name = "chemotaxis_direction"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_direction; 
		std::cout <<  dTemp << std::endl; 
 		// name = "chemotactic_sensitivities"; 
		std::cout <<  pCell->phenotype.motility.chemotactic_sensitivities.data() << std::endl; // , sizeof(double) , m , fp ); 

// secretion 
        std::cout << "------ secretion:\n";
 		// name = "secretion_rates"; 
		std::cout <<  pCell->phenotype.secretion.secretion_rates.data() << std::endl; // , sizeof(double) , m , fp ); 
	 	// name = "uptake_rates"; 
		std::cout <<  pCell->phenotype.secretion.uptake_rates.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "saturation_densities"; 
		std::cout <<  pCell->phenotype.secretion.saturation_densities.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "net_export_rates"; 
		std::cout <<  pCell->phenotype.secretion.net_export_rates.data() << std::endl; // , sizeof(double) , m , fp ); 

// molecular 
        std::cout << "------ molecular:\n";
 		// name = "internalized_total_substrates"; 
		std::cout <<  pCell->phenotype.molecular.internalized_total_substrates.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "fraction_released_at_death"; 
		std::cout <<  pCell->phenotype.molecular.fraction_released_at_death.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "fraction_transferred_when_ingested"; 
		std::cout <<  pCell->phenotype.molecular.fraction_transferred_when_ingested.data() << std::endl; // , sizeof(double) , m , fp ); 

// interactions 
        std::cout << "------ interactions:\n";
	/*
 		// name = "dead_phagocytosis_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.dead_phagocytosis_rate << std::endl; 
	*/
 		// name = "apoptotic_phagocytosis_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.apoptotic_phagocytosis_rate << std::endl; 
 		// name = "necrotic_phagocytosis_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.necrotic_phagocytosis_rate << std::endl; 
 		// name = "other_dead_phagocytosis_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.other_dead_phagocytosis_rate << std::endl; 
 		// name = "live_phagocytosis_rates"; 
		std::cout <<  pCell->phenotype.cell_interactions.live_phagocytosis_rates.data() << std::endl; // , sizeof(double) , n , fp ); 

 		// name = "attack_rates"; 
		std::cout <<  pCell->phenotype.cell_interactions.attack_rates.data() << std::endl; // , sizeof(double) , n , fp ); 
 		// name = "immunogenicities"; 
		std::cout <<  pCell->phenotype.cell_interactions.immunogenicities.data() << std::endl; // , sizeof(double) , n , fp ); 
 		// name = "attack_target"; 
		Cell* pTarget = pCell->phenotype.cell_interactions.pAttackTarget; 
		int AttackID = -1; 
		if( pTarget )
		{ AttackID = pTarget->ID; }
		dTemp = (double) AttackID; 
		std::cout <<  "AttackID= " << dTemp << std::endl;
 		// name = "attack_damage_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.attack_damage_rate << std::endl; 
 		// name = "attack_duration"; 
		std::cout <<  pCell->phenotype.cell_interactions.attack_duration << std::endl; 
 		// name = "total_damage_delivered"; 
		std::cout <<  pCell->phenotype.cell_interactions.total_damage_delivered << std::endl; 

 		// name = "fusion_rates"; 
		std::cout <<  pCell->phenotype.cell_interactions.fusion_rates.data() << std::endl; // , sizeof(double) , n , fp ); 

// transformations 
        std::cout << "------ transformation_rates:\n";
  		// name = "transformation_rates"; 
		std::cout <<  pCell->phenotype.cell_transformations.transformation_rates.data() << std::endl; // , sizeof(double) , n , fp ); 

// asymmetric division
        std::cout << "------ :\n";
		// name = "asymmetric_division_rate"; 
		std::cout <<  "asymmetric_division_probabilities= " << pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.data() << std::endl; // , sizeof(double) , n, fp );

	// cell integrity 
        std::cout << "------ cell integrity:\n";
 		// name = "damage"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage << std::endl; 
 		// name = "damage_rate"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage_rate << std::endl; 
 		// name = "damage_repair_rate"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage_repair_rate << std::endl; 

// custom 
		// custom scalar variables 
        std::cout << "------ custom data: scalar:\n";
		for( int j=0 ; j < (*all_cells)[0]->custom_data.variables.size(); j++ )
		{ std::cout <<  pCell->custom_data.variables[j].value << std::endl; }

		// custom vector variables 
        std::cout << "------ custom data: vectors:\n";
		for( int j=0 ; j < (*all_cells)[0]->custom_data.vector_variables.size(); j++ )
		{
			int size_temp = pCell->custom_data.vector_variables[j].value.size(); 
			std::cout <<  pCell->custom_data.vector_variables[j].value.data() << std::endl; // , sizeof(double) , size_temp , fp );
		}
	}
    std::cout << "---------- END dump_cell_mat_vars:\n";
}


void add_PhysiCell_cells_to_open_xml_pugi_v2( pugi::xml_document& xml_dom, std::string filename_base, Microenvironment& M  ) 
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this one July 2024

	// get number of substrates 
	static int m =  microenvironment.number_of_densities(); // number_of_substrates  
	// get number of cell types
	static int n = cell_definition_indices_by_name.size(); // number_of_cell_types
	// get number of death models 
	static int nd = (*all_cells)[0]->phenotype.death.rates.size(); // 
	// get number of custom data 
	static int nc = 0; // 
	static int nc_scalar = 0; 
	static int nc_vector = 0; 

	static int cell_data_size = 0; 

	static bool legend_done = false; 
	static std::vector<std::string> data_names; 
	static std::vector<std::string> data_units; 
	static std::vector<int> data_start_indices; 
	static std::vector<int> data_sizes; 

	// set up the cell types labels 

	static bool cell_types_legend_done = false; 
	static std::vector<std::string> cell_type_names; 
	static std::vector<int> cell_type_indices; 
	static std::vector<int> cell_type_IDs; 

	if( cell_types_legend_done == false )
	{
		cell_type_names.clear(); 
		cell_type_IDs.clear(); 
		cell_type_indices.clear(); 

		for( int j=0; j < cell_definitions_by_index.size() ; j++ )
		{
			Cell_Definition* pCD = cell_definitions_by_index[j]; 
			int index = j; 
			int type = pCD->type; 
			std::string name = pCD->name; 
			cell_type_names.push_back( name ); 
			cell_type_IDs.push_back( type ); 
			cell_type_indices.push_back( index); 
		}

		cell_types_legend_done = true; 
	}

	// set up the labels 
	if( legend_done == false )
	{
		data_names.clear(); 
		data_sizes.clear(); 
		data_start_indices.clear(); 
		data_units.clear(); 

		std::string name; 
		std::string units; 
		int size; 
		int index = 0; 

// compatibilty : first 17 entries 
		// ID 					<label index="0" size="1">ID</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"ID" , "none" , 1 ) ; 

		//					<label index="1" size="3">position</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"position" , "microns" , 3 ); 

		//					<label index="4" size="1">total_volume</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"total_volume" , "cubic microns" , 1 ); 
		
		//					<label index="5" size="1">cell_type</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_type" , "none" , 1 ); 

		//					<label index="6" size="1">cycle_model</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cycle_model" , "none" , 1 ); 

		//					<label index="7" size="1">current_phase</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"current_phase" , "none" , 1 ); 

		//					<label index="8" size="1">elapsed_time_in_phase</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"elapsed_time_in_phase" , "min" , 1 ); 

		//					<label index="9" size="1">nuclear_volume</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"nuclear_volume" , "cubic microns" , 1 ); 
		//					<label index="10" size="1">cytoplasmic_volume</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cytoplasmic_volume" , "cubic microns" , 1 ); 

		//					<label index="11" size="1">fluid_fraction</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"fluid_fraction" , "none" , 1 ); 

		//					<label index="12" size="1">calcified_fraction</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"calcified_fraction" , "none" , 1 ); 

		//					<label index="13" size="3">orientation</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"orientation" , "none" , 3 ); 
			
		//					<label index="16" size="1">polarity</label>
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"polarity" , "none" , 1 ); 

	/* state variables to save */ 
	// state 
		// velocity // 3 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"velocity" , "micron/min" , 3 ); 

		// pressure // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"pressure" , "none" , 1 ); 

		// number of nuclei // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"number_of_nuclei" , "none" , 1 ); 

		// damage // 1 // this is in cell_integrity now 
		// add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
		//	"damage" , "none" , 1 ); 


		// total attack time // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"total_attack_time" , "min" , 1 ); 

		// contact_with_basement_membrane // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"contact_with_basement_membrane" , "none" , 1 ); 

	/* now go through phenotype and state */ 
		// cycle 
		// cycle model // already above 
		// current phase // already above 
		// current exit rate // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"current_cycle_phase_exit_rate" , "1/min" , 1 ); 

	  // elapsed time in phase // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"elapsed_time_in_phase" , "min" , 1 ); 

		// death 
		// live or dead state // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"dead" , "none" , 1 ); 

	    // current death model // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"current_death_model" , "none" , 1 ); 

	    // death rates // nd 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"death_rates" , "1/min" , nd ); 
		// 

		// volume ()
  		// cytoplasmic_biomass_change_rate // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cytoplasmic_biomass_change_rate" , "1/min" , 1 ); 

  		//  nuclear_biomass_change_rate;  // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"nuclear_biomass_change_rate" , "1/min" , 1 ); 

	    //  fluid_change_rate; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"fluid_change_rate" , "1/min" , 1 ); 

		// calcification_rate; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"calcification_rate" , "1/min" , 1 ); 

		// target_solid_cytoplasmic; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"target_solid_cytoplasmic" , "cubic microns" , 1 ); 

		// target_solid_nuclear; // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"target_solid_nuclear" , "cubic microns" , 1 ); 

		// target_fluid_fraction; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"target_fluid_fraction" , "none" , 1 ); 

	// geometry 
		// radius //1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"radius" , "microns" , 1 ); 

		// nuclear_radius // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"nuclear_radius" , "microns" , 1 ); 

	 	// surface_area //1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"surface_area" , "square microns" , 1 ); 

		// polarity // arleady done 

	// mechanics 
		// cell_cell_adhesion_strength; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_cell_adhesion_strength" , "micron/min" , 1 ); 

		// cell_BM_adhesion_strength; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_BM_adhesion_strength" , "micron/min" , 1 ); 
  	  		
		// cell_cell_repulsion_strength; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_cell_repulsion_strength" , "micron/min" , 1 ); 

		// cell_BM_repulsion_strength; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_BM_repulsion_strength" , "micron/min" , 1 ); 
 
		// std::vector<double> cell_adhesion_affinities; // n 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"cell_adhesion_affinities" , "none" , n ); 

		// relative_maximum_adhesion_distance; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"relative_maximum_adhesion_distance" , "none" , 1 ); 

		// maximum_number_of_attachments; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"maximum_number_of_attachments" , "none" , 1 ); 

		// attachment_elastic_constant; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attachment_elastic_constant" , "1/min" , 1 ); 

		// attachment_rate; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attachment_rate" , "1/min" , 1 ); 

		// detachment_rate; // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"detachment_rate" , "1/min" , 1 ); 

	 // Motility
		// is_motile // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"is_motile" , "none" , 1 ); 

		//  persistence_time; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"persistence_time" , "min" , 1 ); 

		// migration_speed; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"migration_speed" , "micron/min" , 1 ); 

		// std::vector<double> migration_bias_direction; // 3 // motility_bias_direction originally 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"migration_bias_direction" , "none" , 3 ); 

		// migration_bias; //1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"migration_bias" , "none" , 1 ); 

		// std::vector<double> motility_vector; // 3
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"motility_vector" , "micron/min" , 3 ); 

		// chemotaxis_index; // 1
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"chemotaxis_index" , "none" , 1 ); 

		// chemotaxis_direction; // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"chemotaxis_direction" , "none" , 1 ); 

		// advanced chemotaxis 
		// std::vector<double> chemotactic_sensitivities;  // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"chemotactic_sensitivities" , "none" , m ); 

		// secretion 
		// std::vector<double> secretion_rates; // m
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"secretion_rates" , "1/min" , m ); 
	
		// std::vector<double> uptake_rates; // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"uptake_rates" , "1/min" , m ); 

		// std::vector<double> saturation_densities; // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"saturation_densities" , "stuff/cubic micron" , m ); 

		// std::vector<double> net_export_rates; // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"net_export_rates" , "stuff/min" , m ); 

	// molecular 
		// internalized_total_substrates // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"internalized_total_substrates" , "stuff" , m ); 

		// 	fraction_released_at_death // m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"fraction_released_at_death" , "none" , m ); 

		// fraction_transferred_when_ingested //m 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"fraction_transferred_when_ingested" , "none" , m ); 

	// interactions 
/*	
		// dead_phagocytosis_rate; // 1 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"dead_phagocytosis_rate" , "1/min" , 1 ); 
*/

		// apoptotic phagocytosis_rate // new 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"apoptotic_phagocytosis_rate" , "1/min" , 1 ); 

		// necrotic phagocytosis rate // new 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"necrotic_phagocytosis_rate" , "1/min" , 1 ); 

		// other dead phagocytosis rate // new 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"other_dead_phagocytosis_rate" , "1/min" , 1 ); 

		// std::vector<double> live_phagocytosis_rates; // n 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"live_phagocytosis_rates" , "1/min" , n ); 

		// std::vector<double> attack_rates; // n 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attack_rates" , "1/min" , n ); 

		// std::vector<double> immunogenicities;  n // was missing  
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"immunogenicities" , "none" , n ); 

		// pAttackTarget;  1 // new -- use the cell ID 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attack_target" , "none" , 1 ); 

		// double (attack_)damage_rate;  1 // changed from damage_rate 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attack_damage_rate" , "1/min" , 1 ); 

		// double attack_duration;  1 // new 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attack_duration" , "min" , 1 ); 

		// double total_damage_delivered;  1 // new 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"attack_total_damage_delivered" , "none" , 1 ); 

		// std::vector<double> fusion_rates; // n 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"fusion_rates" , "1/min" , n ); 

	// transformations 
		// std::vector<double> transformation_rates; // n 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"transformation_rates" , "1/min" , n ); 

	// asymmetric division
		// std::vector<double> asymmetric_division_probabilities; // n
		add_variable_to_labels( data_names, data_units, data_start_indices, data_sizes, 
			"asymmetric_division_probabilities" , "none" , n );

	// cell integrity 

		// double damage; 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"damage" , "none" , 1 ); 

		// double damage_rate; // new use of old name! now the rate of undergoing damage (not by attack)
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"damage_rate" , "1/min" , 1 ); 

		// double damage_repair_rate; 
		add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
			"damage_repair_rate" , "1/min" , 1 ); 

// custom 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.variables.size(); j++ )
		{		
			name = (*all_cells)[0]->custom_data.variables[j].name; 
			units = (*all_cells)[0]->custom_data.variables[j].units; 
			add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
				name,units,1 ); 
		}
		
		// custom vector variables 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.vector_variables.size(); j++ )
		{
			name = (*all_cells)[0]->custom_data.vector_variables[j].name; 
			units = (*all_cells)[0]->custom_data.vector_variables[j].units; 
			size = (*all_cells)[0]->custom_data.vector_variables[j].value.size(); 
			add_variable_to_labels( data_names,data_units,data_start_indices,data_sizes, 
				name,units,size ); 
		}

		cell_data_size = total_data_size( data_sizes ); 
		legend_done = true; 
	}

	// get ready for XML navigation 
	// pugi::xml_document& xml_dom = BioFVM::biofvm_doc; 

	pugi::xml_node root = xml_dom.child("MultiCellDS") ; 
	pugi::xml_node node = root.child( "cellular_information" ); // root = cellular_information
	root = node; 	

	// Let's reduce memory allocations and sprintf calls. 
	// This reduces execution time. 
	static char* temp; 
	static bool initialized = false; 
	
	static char rate_chars [1024]; 
	static char volume_chars [1024]; 
	static char diffusion_chars [1024]; 
	if( !initialized )
	{ 
		temp = new char [1024]; 
		initialized = true; 
		
		sprintf( rate_chars, "1/%s" , M.time_units.c_str() ); 
		sprintf( volume_chars, "%s^3" , M.spatial_units.c_str() ); 
		sprintf( diffusion_chars , "%s^2/%s", M.spatial_units.c_str() , M.time_units.c_str() ); 
	}

	node = node.child( "cell_populations" ); 
	if( !node )
	{
		node = root.append_child( "cell_populations" ); 
	}
	root = node; // root = cellular_information.cell_populations 

	node = root.child( "cell_population" ); 
	if( !node )
	{
		node = root.append_child( "cell_population" ); 
		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "individual" ); 
	}
	root = node; // root = cellular_information.cell_populations.cell_population  

	node = root.child( "custom" ); 
	if( !node )
	{
		node = root.append_child( "custom" ); 
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom 

	node = root.child( "simplified_data" ); 
	if( !node )
	{
		node = root.append_child( "simplified_data" ); 
		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "matlab" ); 

		attrib = node.append_attribute( "source" ); 
		attrib.set_value( "PhysiCell" ); 

		attrib = node.append_attribute( "data_version" ); 
		attrib.set_value( "2" ); 
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.simplified_data 

	// write cell definiton labels // new in July 2024
	node = root.child( "cell_types" ); 
	if( !node )
	{
		node = root.append_child( "cell_types" ); 
		root = node; // root = cellular_information.cell_populations.cell_population.custom.simplified_data.cell_types   

		for( int i=0; i < cell_type_names.size(); i++ )
		{
			node = root.append_child( "type" ); 	

			pugi::xml_attribute attrib = node.append_attribute( "ID" ); 
			attrib.set_value( cell_type_indices[i] ); 	

			attrib = node.append_attribute( "type" ); 
			attrib.set_value( cell_type_IDs[i] ); 				

			node.append_child( pugi::node_pcdata ).set_value( cell_type_names[i].c_str() ); 
		}
		root = root.parent(); // root = cellular_information.cell_populations.cell_population.custom.simplified_data  

		// <type ID="" type="">name</type> 
	} 

	// write legend 

	node = root.child( "labels" ); 
	if( !node )
	{
		node = root.append_child( "labels" ); 
		root = node; // root = cellular_information.cell_populations.cell_population.custom.simplified_data.labels  

		for( int i=0; i < data_names.size(); i++ )
		{
			node = root.append_child( "label" ); 	
			pugi::xml_attribute attrib = node.append_attribute( "index" ); 
			attrib.set_value( data_start_indices[i] ); 	

			attrib = node.append_attribute( "size" ); 
			attrib.set_value( data_sizes[i] ); 				

			attrib = node.append_attribute( "units" ); 
			attrib.set_value( data_units[i].c_str() ); 	

			node.append_child( pugi::node_pcdata ).set_value( data_names[i].c_str() ); 
		}
		root = root.parent(); // root = cellular_information.cell_populations.cell_population.custom.simplified_data  
	}

	// write data 
	node = root.child( "filename" ); 
	if( !node )
	{
		node = root.append_child( "filename" ); 
	}

	// write the filename 

	// next, filename 
	char filename [1024]; 
	sprintf( filename , "%s_cells.mat" , filename_base.c_str() ); 
	
	/* store filename without the relative pathing (if any) */ 
	char filename_without_pathing [1024];
	char* filename_start = strrchr( filename , '/' ); 
	if( filename_start == NULL )
	{ filename_start = filename; }
	else	
	{ filename_start++; } 
	strcpy( filename_without_pathing , filename_start );  
	
	if( !node.first_child() )
	{
		node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
	}
	else
	{
		node.first_child().set_value( filename_without_pathing ); // filename ); 
	}

	// now write the actual data 

	int size_of_each_datum = cell_data_size;
	int number_of_data_entries = (*all_cells).size();  

	FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "cells" );  
	if( fp == NULL )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for MAT writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 

	Cell* pCell; 
	
	double dTemp; 
	// storing data as cols (each column is a cell)
	for( int i=0; i < number_of_data_entries ; i++ )
	{
		pCell = (*all_cells)[i]; 

		int writes = 0; 

		// compatibilty : first 17 entries 
		// ID 					<label index="0" size="1">ID</label>
		// double ID_temp = (double) (*all_cells)[i]->ID;
		// fwrite( (char*) &( ID_temp ) , sizeof(double) , 1 , fp ); 

		// name = "ID"; 
		dTemp = (double) pCell->ID;
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// name = "position";    NOTE very different syntax for writing vectors!
        std::fwrite( pCell->position.data() , sizeof(double) , 3 , fp );
		// name = "total_volume"; 
		std::fwrite( &( pCell->phenotype.volume.total ) , sizeof(double) , 1 , fp ); 
		// name = "cell_type"; 
		dTemp = (double) pCell->type;
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// name = "cycle_model"; 
		dTemp = (double) pCell->phenotype.cycle.model().code; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); // cycle model 
		// name = "current_phase"; 
		dTemp = (double) pCell->phenotype.cycle.current_phase().code; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); // cycle model 
		// name = "elapsed_time_in_phase"; 
		std::fwrite( &( pCell->phenotype.cycle.data.elapsed_time_in_phase ) , sizeof(double) , 1 , fp ); 
		// name = "nuclear_volume"; 
		std::fwrite( &( pCell->phenotype.volume.nuclear ) , sizeof(double) , 1 , fp );   
		// name = "cytoplasmic_volume"; 
		std::fwrite( &( pCell->phenotype.volume.cytoplasmic ) , sizeof(double) , 1 , fp );
		// name = "fluid_fraction"; 
		std::fwrite( &( pCell->phenotype.volume.fluid_fraction ) , sizeof(double) , 1 , fp );
		// name = "calcified_fraction"; 
		std::fwrite( &( pCell->phenotype.volume.calcified_fraction ) , sizeof(double) , 1 , fp ); 
		// name = "orientation"; 
		std::fwrite( pCell->state.orientation.data() , sizeof(double) , 3 , fp ); 
		// name = "polarity"; 
		std::fwrite( &( pCell->phenotype.geometry.polarity ) , sizeof(double) , 1 , fp ); 

 /* state variables to save */ 
// state
		// name = "velocity"; 
		std::fwrite( pCell->velocity.data() , sizeof(double) , 3 , fp ); 
		// name = "pressure"; 
		std::fwrite( &( pCell->state.simple_pressure ) , sizeof(double) , 1 , fp ); 
		// name = "number_of_nuclei"; 
		dTemp = (double) pCell->state.number_of_nuclei; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// // name = "damage"; 
		// std::fwrite( &( pCell->phenotype.integrity.damage ) , sizeof(double) , 1 , fp ); 
		// name = "total_attack_time"; 
		std::fwrite( &( pCell->state.total_attack_time ) , sizeof(double) , 1 , fp ); 
		// name = "contact_with_basement_membrane"; 
		dTemp = (double) pCell->state.contact_with_basement_membrane; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 

/* now go through phenotype and state */ 
// cycle 
  // current exit rate // 1 
		// name = "current_cycle_phase_exit_rate"; 
		int phase_index = pCell->phenotype.cycle.data.current_phase_index; 
		std::fwrite( &( pCell->phenotype.cycle.data.exit_rate(phase_index) ) , sizeof(double) , 1 , fp ); 
		// name = "elapsed_time_in_phase"; 
		std::fwrite( &( pCell->phenotype.cycle.data.elapsed_time_in_phase ) , sizeof(double) , 1 , fp ); 

// death 
  // live or dead state // 1 
		// name = "dead"; 
		dTemp = (double) pCell->phenotype.death.dead; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// name = "current_death_model"; // 
		dTemp = (double) pCell->phenotype.death.current_death_model_index; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// name = "death_rates"; 
		std::fwrite( pCell->phenotype.death.rates.data() , sizeof(double) , nd , fp ); 
		
	// volume ()
		// name = "cytoplasmic_biomass_change_rate"; 
		std::fwrite( &( pCell->phenotype.volume.cytoplasmic_biomass_change_rate ) , sizeof(double) , 1 , fp ); 
		// name = "nuclear_biomass_change_rate"; 
		std::fwrite( &( pCell->phenotype.volume.nuclear_biomass_change_rate ) , sizeof(double) , 1 , fp ); 
		// name = "fluid_change_rate"; 
		std::fwrite( &( pCell->phenotype.volume.fluid_change_rate ) , sizeof(double) , 1 , fp ); 
		// name = "calcification_rate"; 
		std::fwrite( &( pCell->phenotype.volume.calcification_rate ) , sizeof(double) , 1 , fp ); 
		// name = "target_solid_cytoplasmic"; 
		std::fwrite( &( pCell->phenotype.volume.target_solid_cytoplasmic ) , sizeof(double) , 1 , fp ); 
		// name = "target_solid_nuclear"; 
		std::fwrite( &( pCell->phenotype.volume.target_solid_nuclear ) , sizeof(double) , 1 , fp ); 
		// name = "target_fluid_fraction"; 
		std::fwrite( &( pCell->phenotype.volume.target_fluid_fraction ) , sizeof(double) , 1 , fp ); 

  // geometry 
     // radius //1 
		// name = "radius"; 
		std::fwrite( &( pCell->phenotype.geometry.radius ) , sizeof(double) , 1 , fp ); 
		// name = "nuclear_radius"; 
		std::fwrite( &( pCell->phenotype.geometry.nuclear_radius ) , sizeof(double) , 1 , fp ); 
		// name = "surface_area"; 
		std::fwrite( &( pCell->phenotype.geometry.surface_area ) , sizeof(double) , 1 , fp ); 

  // mechanics 
	// cell_cell_adhesion_strength; // 1
		// name = "cell_cell_adhesion_strength"; 
		std::fwrite( &( pCell->phenotype.mechanics.cell_cell_adhesion_strength ) , sizeof(double) , 1 , fp ); 
		// name = "cell_BM_adhesion_strength"; 
		std::fwrite( &( pCell->phenotype.mechanics.cell_BM_adhesion_strength ) , sizeof(double) , 1 , fp ); 
		// name = "cell_cell_repulsion_strength"; 
		std::fwrite( &( pCell->phenotype.mechanics.cell_cell_repulsion_strength ) , sizeof(double) , 1 , fp ); 
		// name = "cell_BM_repulsion_strength"; 
		std::fwrite( &( pCell->phenotype.mechanics.cell_BM_repulsion_strength ) , sizeof(double) , 1 , fp ); 
		// name = "cell_adhesion_affinities"; 
		std::fwrite( pCell->phenotype.mechanics.cell_adhesion_affinities.data() , sizeof(double) , n , fp ); 
		// name = "relative_maximum_adhesion_distance"; 
		std::fwrite( &( pCell->phenotype.mechanics.relative_maximum_adhesion_distance ) , sizeof(double) , 1 , fp ); 
		// name = "maximum_number_of_attachments"; 
		dTemp = (double) pCell->phenotype.mechanics.maximum_number_of_attachments; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
		// name = "attachment_elastic_constant"; 
		std::fwrite( &( pCell->phenotype.mechanics.attachment_elastic_constant ) , sizeof(double) , 1 , fp ); 
		// name = "attachment_rate"; 
		std::fwrite( &( pCell->phenotype.mechanics.attachment_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "detachment_rate"; 
		std::fwrite( &( pCell->phenotype.mechanics.detachment_rate ) , sizeof(double) , 1 , fp ); 

 // Motility
 		// name = "is_motile"; 
		dTemp = (double) pCell->phenotype.motility.is_motile; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
 		// name = "persistence_time"; 
		std::fwrite( &( pCell->phenotype.motility.persistence_time ) , sizeof(double) , 1 , fp ); 
 		// name = "migration_speed"; 
		std::fwrite( &( pCell->phenotype.motility.migration_speed ) , sizeof(double) , 1 , fp ); 
 		// name = "migration_bias_direction"; 
		std::fwrite( pCell->phenotype.motility.migration_bias_direction.data() , sizeof(double) , 3 , fp ); 
 		// name = "migration_bias"; 
		std::fwrite( &( pCell->phenotype.motility.migration_bias ) , sizeof(double) , 1 , fp ); 
 		// name = "motility_vector"; 
		std::fwrite( pCell->phenotype.motility.motility_vector.data() , sizeof(double) , 3 , fp ); 
 		// name = "chemotaxis_index"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_index; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
 		// name = "chemotaxis_direction"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_direction; 
		std::fwrite( &( dTemp ) , sizeof(double) , 1 , fp ); 
 		// name = "chemotactic_sensitivities"; 
		std::fwrite( pCell->phenotype.motility.chemotactic_sensitivities.data() , sizeof(double) , m , fp ); 

// secretion 
 		// name = "secretion_rates"; 
		std::fwrite( pCell->phenotype.secretion.secretion_rates.data() , sizeof(double) , m , fp ); 
	 	// name = "uptake_rates"; 
		std::fwrite( pCell->phenotype.secretion.uptake_rates.data() , sizeof(double) , m , fp ); 
 		// name = "saturation_densities"; 
		std::fwrite( pCell->phenotype.secretion.saturation_densities.data() , sizeof(double) , m , fp ); 
 		// name = "net_export_rates"; 
		std::fwrite( pCell->phenotype.secretion.net_export_rates.data() , sizeof(double) , m , fp ); 

// molecular 
 		// name = "internalized_total_substrates"; 
		std::fwrite( pCell->phenotype.molecular.internalized_total_substrates.data() , sizeof(double) , m , fp ); 
 		// name = "fraction_released_at_death"; 
		std::fwrite( pCell->phenotype.molecular.fraction_released_at_death.data() , sizeof(double) , m , fp ); 
 		// name = "fraction_transferred_when_ingested"; 
		std::fwrite( pCell->phenotype.molecular.fraction_transferred_when_ingested.data() , sizeof(double) , m , fp ); 

// interactions 
	/*
 		// name = "dead_phagocytosis_rate"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.dead_phagocytosis_rate ) , sizeof(double) , 1 , fp ); 
	*/
 		// name = "apoptotic_phagocytosis_rate"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.apoptotic_phagocytosis_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "necrotic_phagocytosis_rate"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.necrotic_phagocytosis_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "other_dead_phagocytosis_rate"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.other_dead_phagocytosis_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "live_phagocytosis_rates"; 
		std::fwrite( pCell->phenotype.cell_interactions.live_phagocytosis_rates.data() , sizeof(double) , n , fp ); 

 		// name = "attack_rates"; 
		std::fwrite( pCell->phenotype.cell_interactions.attack_rates.data() , sizeof(double) , n , fp ); 
 		// name = "immunogenicities"; 
		std::fwrite( pCell->phenotype.cell_interactions.immunogenicities.data() , sizeof(double) , n , fp ); 
 		// name = "attack_target"; 
		Cell* pTarget = pCell->phenotype.cell_interactions.pAttackTarget; 
		int AttackID = -1; 
		if( pTarget )
		{ AttackID = pTarget->ID; }
		dTemp = (double) AttackID; 
		std::fwrite( &(dTemp) , sizeof(double) , 1 , fp ); 
 		// name = "attack_damage_rate"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.attack_damage_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "attack_duration"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.attack_duration ) , sizeof(double) , 1 , fp ); 
 		// name = "total_damage_delivered"; 
		std::fwrite( &( pCell->phenotype.cell_interactions.total_damage_delivered ) , sizeof(double) , 1 , fp ); 

 		// name = "fusion_rates"; 
		std::fwrite( pCell->phenotype.cell_interactions.fusion_rates.data() , sizeof(double) , n , fp ); 

// transformations 
  		// name = "transformation_rates"; 
		std::fwrite( pCell->phenotype.cell_transformations.transformation_rates.data() , sizeof(double) , n , fp ); 

// asymmetric division
		// name = "asymmetric_division_rate"; 
		std::fwrite( pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.data() , sizeof(double) , n, fp );

	// cell integrity 
 		// name = "damage"; 
		std::fwrite( &( pCell->phenotype.cell_integrity.damage ) , sizeof(double) , 1 , fp ); 
 		// name = "damage_rate"; 
		std::fwrite( &( pCell->phenotype.cell_integrity.damage_rate ) , sizeof(double) , 1 , fp ); 
 		// name = "damage_repair_rate"; 
		std::fwrite( &( pCell->phenotype.cell_integrity.damage_repair_rate ) , sizeof(double) , 1 , fp ); 

// custom 
		// custom scalar variables 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.variables.size(); j++ )
		{ std::fwrite( &( pCell->custom_data.variables[j].value ) , sizeof(double) , 1 , fp ); }

		// custom vector variables 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.vector_variables.size(); j++ )
		{
			int size_temp = pCell->custom_data.vector_variables[j].value.size(); 
			std::fwrite( pCell->custom_data.vector_variables[j].value.data() , sizeof(double) , size_temp , fp );
		}
	}

	fclose( fp ); 
    dump_cell_mat_vars(number_of_data_entries);

#ifdef ADDON_PHYSIBOSS

	// PhysiBoSS Intracellular Data
	node = node.parent().parent();  // custom 

	root = node; 
	node = node.child( "boolean_intracellular_data" );  
	if( !node )
	{
		node = root.append_child( "boolean_intracellular_data" ); 

		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "text" ); 		

		attrib = node.append_attribute( "source" ); 
		attrib.set_value( "PhysiBoSS" ); 		

		attrib = node.append_attribute( "data_version" ); 
		attrib.set_value( "2" ); 	
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.intracellular_data
	node = root.child( "filename"); 
	if( !node )
	{
		node = root.append_child( "filename" ); 

	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.intracellular_data.filename


	// next, filename 
	sprintf( filename , "%s_boolean_intracellular.csv" , filename_base.c_str() ); 
		
	/* store filename without the relative pathing (if any) */ 
	filename_start = strrchr( filename , '/' ); 
	if( filename_start == NULL )
	{ filename_start = filename; }
	else	
	{ filename_start++; } 
	strcpy( filename_without_pathing , filename_start );  
	
	if( !node.first_child() )
	{
		node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
	}
	else
	{
		node.first_child().set_value( filename_without_pathing ); // filename ); 
	}	

	MaBoSSIntracellular::save( filename );

#endif

	// neighbor graph 
	node = node.parent().parent();  // custom 

	root = node; 
	node = node.child( "neighbor_graph" );  
	if( !node )
	{
		node = root.append_child( "neighbor_graph" ); 

		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "text" ); 		

		attrib = node.append_attribute( "source" ); 
		attrib.set_value( "PhysiCell" ); 		

		attrib = node.append_attribute( "data_version" ); 
		attrib.set_value( "2" ); 	
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.neighbor_graph
	node = root.child( "filename"); 
	if( !node )
	{
		node = root.append_child( "filename" ); 

	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.neighbor_graph.filename 


	// next, filename 
	sprintf( filename , "%s_cell_neighbor_graph.txt" , filename_base.c_str() ); 
		
	/* store filename without the relative pathing (if any) */ 
	filename_start = strrchr( filename , '/' ); 
	if( filename_start == NULL )
	{ filename_start = filename; }
	else	
	{ filename_start++; } 
	strcpy( filename_without_pathing , filename_start );  
	
	if( !node.first_child() )
	{
		node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
	}
	else
	{
		node.first_child().set_value( filename_without_pathing ); // filename ); 
	}	

	write_neighbor_graph( filename ); 


	// attached cell graph 
	node = root; 
	node = node.parent().parent(); // root = cellular_information.cell_populations.cell_population.custom

	root = node; 
	node = node.child( "attached_cells_graph" );  
	if( !node )
	{
		node = root.append_child( "attached_cells_graph" ); 

		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "text" ); 		

		attrib = node.append_attribute( "source" ); 
		attrib.set_value( "PhysiCell" ); 		

		attrib = node.append_attribute( "data_version" ); 
		attrib.set_value( "2" ); 	
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.attached_cells_graph
	node = root.child( "filename"); 
	if( !node )
	{ node = root.append_child( "filename" ); }
	root = node; // root = cellular_information.cell_populations.cell_population.custom.attached_cells_graph.filename 


	// next, filename 
	sprintf( filename , "%s_attached_cells_graph.txt" , filename_base.c_str() ); 
		
	/* store filename without the relative pathing (if any) */ 
	filename_start = strrchr( filename , '/' ); 
	if( filename_start == NULL )
	{ filename_start = filename; }
	else	
	{ filename_start++; } 
	strcpy( filename_without_pathing , filename_start );  
	
	if( !node.first_child() )
	{
		node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
	}
	else
	{
		node.first_child().set_value( filename_without_pathing ); // filename ); 
	}	

	write_attached_cells_graph( filename );	

	// spring attached cell graph 
	node = root; 
	node = node.parent().parent(); // root = cellular_information.cell_populations.cell_population.custom

	root = node; 
	node = node.child( "spring_attached_cells_graph" );  
	if( !node )
	{
		node = root.append_child( "spring_attached_cells_graph" ); 

		pugi::xml_attribute attrib = node.append_attribute( "type" ); 
		attrib.set_value( "text" ); 		

		attrib = node.append_attribute( "source" ); 
		attrib.set_value( "PhysiCell" ); 		

		attrib = node.append_attribute( "data_version" ); 
		attrib.set_value( "2" ); 	
	}
	root = node; // root = cellular_information.cell_populations.cell_population.custom.spring_attached_cells_graph
	node = root.child( "filename"); 
	if( !node )
	{ node = root.append_child( "filename" ); }
	root = node; // root = cellular_information.cell_populations.cell_population.custom.spring_attached_cells_graph.filename 


	// next, filename 
	sprintf( filename , "%s_spring_attached_cells_graph.txt" , filename_base.c_str() ); 
		
	/* store filename without the relative pathing (if any) */ 
	filename_start = strrchr( filename , '/' ); 
	if( filename_start == NULL )
	{ filename_start = filename; }
	else	
	{ filename_start++; } 
	strcpy( filename_without_pathing , filename_start );  
	
	if( !node.first_child() )
	{
		node.append_child( pugi::node_pcdata ).set_value( filename_without_pathing ); // filename ); 
	}
	else
	{
		node.first_child().set_value( filename_without_pathing ); // filename ); 
	}	

	write_spring_attached_cells_graph( filename ); 

	return; 
}


/* end of new stuff July 2024*/


void write_neighbor_graph( std::string filename )
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // We use this July 2024

	std::ofstream of( filename , std::ios::out ); 
	std::stringstream buffer; 

	for( int i=0 ; i < (*all_cells).size(); i++ )
	{
		buffer << (*all_cells)[i]->ID << ": " ; 
		int size = (*all_cells)[i]->state.neighbors.size(); 
		for( int j=0 ; j < size; j++ )
		{
			buffer << (*all_cells)[i]->state.neighbors[j]->ID; 
			if( j != size-1 )
			{ buffer << ","; }
		}
		if( i != (*all_cells).size()-1 )
		{ buffer << std::endl; }
		of << buffer.rdbuf(); 
	}
	of.close(); 

	return; 
}


void write_attached_cells_graph( std::string filename ) 
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this July 2024

	std::ofstream of( filename , std::ios::out ); 
	std::stringstream buffer; 

	for( int i=0 ; i < (*all_cells).size(); i++ )
	{
		buffer << (*all_cells)[i]->ID << ": " ; 
		int size = (*all_cells)[i]->state.attached_cells.size(); 
		for( int j=0 ; j < size; j++ )
		{
			buffer << (*all_cells)[i]->state.attached_cells[j]->ID; 
			if( j != size-1 )
			{ buffer << ","; }
		}
		if( i != (*all_cells).size()-1 )
		{ buffer << std::endl; }
		of << buffer.rdbuf(); 
	}
	of.close(); 

	return; 
}
 
void write_spring_attached_cells_graph( std::string filename ) 
{
	// std::cout << __LINE__ << " " << __FUNCTION__ << std::endl; // we use this July 2024 

	std::ofstream of( filename , std::ios::out ); 
	std::stringstream buffer; 

	for( int i=0 ; i < (*all_cells).size(); i++ )
	{
		buffer << (*all_cells)[i]->ID << ": " ; 
		int size = (*all_cells)[i]->state.spring_attachments.size();
		for( int j=0 ; j < size; j++ )
		{
			buffer << (*all_cells)[i]->state.spring_attachments[j]->ID; 
			if( j != size-1 )
			{ buffer << ","; }
		}
		if( i != (*all_cells).size()-1 )
		{ buffer << std::endl; }
		of << buffer.rdbuf(); 
	}
	of.close(); 

	return; 
}

void read_PhysiCell_cells_from_matlab_v2(std::string filename, Microenvironment& M)
{
    std::cout << "------- read_PhysiCell_cells_from_matlab_v2 ----------\n";

    // Get number of substrates, cell types, death models
    static int m = microenvironment.number_of_densities();
    static int n = cell_definition_indices_by_name.size();
    std::cout << "------- number_of_densities= " << m << std::endl;
    std::cout << "------- number of cell types= " << n << std::endl;
    static int nd = 0; // will be set from first cell
    
    // Open the MAT file for reading
    FILE* fp = fopen(filename.c_str(), "rb");
    if (fp == NULL)
    {
        std::cout << std::endl << "Error: Failed to open " << filename << " for reading." << std::endl << std::endl;
        return;
    }
    
    // Read MATLAB v4 header (5 integers)
    int32_t header[5];
    fread(header, sizeof(int32_t), 5, fp);
    
    int32_t type = header[0];        // Data type (should be 0 for double)
    int32_t mrows = header[1];       // Number of rows
    int32_t ncols = header[2];       // Number of columns
    int32_t imagf = header[3];       // Imaginary flag (should be 0)
    int32_t namelen = header[4];     // Length of name + 1
    
    // Read variable name
    char* var_name = new char[namelen];
    fread(var_name, sizeof(char), namelen, fp);
    
    int size_of_each_datum = mrows;
    int number_of_data_entries = ncols;
    
    std::cout << "Reading " << number_of_data_entries << " cells with " 
              << size_of_each_datum << " data points each..." << std::endl;
    
    delete[] var_name;
    
    // Clear existing cells
    // delete_all_cells();  //rwh
    
    double dTemp;
    double position[3];
    int cell_ID;
    double cell_vol;
    double orientation[3];
    double velocity[3];
    double migration_bias_direction[3];
    double motility_vector[3];
    
    // Store attack target IDs for later resolution
    std::vector<int> attack_target_ids;
    
    // Cell_Definition* pCD; 
    // Read each cell
    std::cout << "creating cells..." << std::endl;
    for (int i = 0; i < number_of_data_entries; i++)
    {
        std::cout << " ---  creating cell # " << i << std::endl;
        // Cell* pCell = create_cell();
        
        // ID
        fread(&dTemp, sizeof(double), 1, fp);
        // pCell->ID = (int)dTemp;
        cell_ID = (int)dTemp;
        std::cout << "   cell_ID = " << cell_ID  << std::endl;
        
        // position
        fread(position, sizeof(double), 3, fp);
        std::cout << "   position = " << position[0]<<", "<< position[1] << std::endl;
        // pCell->position[0] = position[0];
        // pCell->position[1] = position[1];
        // pCell->position[2] = position[2];
         // rwh - do these 2 things
        // pCell->assign_position(position[0], position[1], position[2]);
        // pCell->update_voxel_index();
        
        // total_volume
        // fread(&(pCell->phenotype.volume.total), sizeof(double), 1, fp);
        fread(&cell_vol, sizeof(double), 1, fp);
        std::cout << "   cell_vol = " << cell_vol  << std::endl;
        
        // cell_type
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "   cell_type(dTemp) = " << dTemp  << std::endl;
        // pCell->type = (int)dTemp;
        int i_cell_type = (int)dTemp;
        std::cout << "   i_cell_type = " << i_cell_type  << std::endl;

        Cell_Definition* pCD = cell_definitions_by_type[i_cell_type]; 
        // pCD = cell_definitions_by_type[i_cell_type]; 
        Cell* pCell = create_cell( *pCD );

        pCell->ID = cell_ID;

        pCell->assign_position(position[0], position[1], position[2]);
        pCell->update_voxel_index();
        
        // cycle_model
        fread(&dTemp, sizeof(double), 1, fp);
        int cycle_model_code = (int)dTemp;
        
        // current_phase
        fread(&dTemp, sizeof(double), 1, fp);
        int current_phase_code = (int)dTemp;
        
        // elapsed_time_in_phase
        fread(&(pCell->phenotype.cycle.data.elapsed_time_in_phase), sizeof(double), 1, fp);
        
        // nuclear_volume
        fread(&(pCell->phenotype.volume.nuclear), sizeof(double), 1, fp);
        
        // cytoplasmic_volume
        fread(&(pCell->phenotype.volume.cytoplasmic), sizeof(double), 1, fp);
        
        // fluid_fraction
        fread(&(pCell->phenotype.volume.fluid_fraction), sizeof(double), 1, fp);
        
        // calcified_fraction
        fread(&(pCell->phenotype.volume.calcified_fraction), sizeof(double), 1, fp);
        
        // orientation
        fread(orientation, sizeof(double), 3, fp);
        pCell->state.orientation[0] = orientation[0];
        pCell->state.orientation[1] = orientation[1];
        pCell->state.orientation[2] = orientation[2];
        
        // polarity
        fread(&(pCell->phenotype.geometry.polarity), sizeof(double), 1, fp);
        
        // velocity
        fread(velocity, sizeof(double), 3, fp);
        pCell->velocity[0] = velocity[0];
        pCell->velocity[1] = velocity[1];
        pCell->velocity[2] = velocity[2];
        
        // pressure
        fread(&(pCell->state.simple_pressure), sizeof(double), 1, fp);
        
        // number_of_nuclei
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->state.number_of_nuclei = (int)dTemp;
        
        // total_attack_time
        fread(&(pCell->state.total_attack_time), sizeof(double), 1, fp);
        
        // contact_with_basement_membrane
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->state.contact_with_basement_membrane = (bool)dTemp;
        
        // current_cycle_phase_exit_rate
        int phase_index = pCell->phenotype.cycle.data.current_phase_index;
        fread(&(pCell->phenotype.cycle.data.exit_rate(phase_index)), sizeof(double), 1, fp);
        
        // elapsed_time_in_phase (duplicate)
        fread(&(pCell->phenotype.cycle.data.elapsed_time_in_phase), sizeof(double), 1, fp);
        
        // dead
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.death.dead = (bool)dTemp;
        
        // current_death_model
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.death.current_death_model_index = (int)dTemp;
        
        // Set nd from first cell if not set
        if (nd == 0)
        {
            nd = pCell->phenotype.death.rates.size();
        }
        
        // death_rates
        fread(pCell->phenotype.death.rates.data(), sizeof(double), nd, fp);
        
        // Volume rates
        fread(&(pCell->phenotype.volume.cytoplasmic_biomass_change_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.nuclear_biomass_change_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.fluid_change_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.calcification_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.target_solid_cytoplasmic), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.target_solid_nuclear), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.volume.target_fluid_fraction), sizeof(double), 1, fp);
        
        // Geometry
        fread(&(pCell->phenotype.geometry.radius), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.geometry.nuclear_radius), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.geometry.surface_area), sizeof(double), 1, fp);
        
        // Mechanics
        fread(&(pCell->phenotype.mechanics.cell_cell_adhesion_strength), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.cell_BM_adhesion_strength), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.cell_cell_repulsion_strength), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.cell_BM_repulsion_strength), sizeof(double), 1, fp);
        fread(pCell->phenotype.mechanics.cell_adhesion_affinities.data(), sizeof(double), n, fp);
        fread(&(pCell->phenotype.mechanics.relative_maximum_adhesion_distance), sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.mechanics.maximum_number_of_attachments = (int)dTemp;
        fread(&(pCell->phenotype.mechanics.attachment_elastic_constant), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.attachment_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.detachment_rate), sizeof(double), 1, fp);
        
        // Motility
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.motility.is_motile = (bool)dTemp;
        fread(&(pCell->phenotype.motility.persistence_time), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.motility.migration_speed), sizeof(double), 1, fp);
        fread(migration_bias_direction, sizeof(double), 3, fp);
        pCell->phenotype.motility.migration_bias_direction[0] = migration_bias_direction[0];
        pCell->phenotype.motility.migration_bias_direction[1] = migration_bias_direction[1];
        pCell->phenotype.motility.migration_bias_direction[2] = migration_bias_direction[2];
        fread(&(pCell->phenotype.motility.migration_bias), sizeof(double), 1, fp);
        fread(motility_vector, sizeof(double), 3, fp);
        pCell->phenotype.motility.motility_vector[0] = motility_vector[0];
        pCell->phenotype.motility.motility_vector[1] = motility_vector[1];
        pCell->phenotype.motility.motility_vector[2] = motility_vector[2];
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.motility.chemotaxis_index = (int)dTemp;
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.motility.chemotaxis_direction = (int)dTemp;
        fread(pCell->phenotype.motility.chemotactic_sensitivities.data(), sizeof(double), m, fp);
        
        // Secretion
        fread(pCell->phenotype.secretion.secretion_rates.data(), sizeof(double), m, fp);
        fread(pCell->phenotype.secretion.uptake_rates.data(), sizeof(double), m, fp);
        fread(pCell->phenotype.secretion.saturation_densities.data(), sizeof(double), m, fp);
        fread(pCell->phenotype.secretion.net_export_rates.data(), sizeof(double), m, fp);
        
        // Molecular
        fread(pCell->phenotype.molecular.internalized_total_substrates.data(), sizeof(double), m, fp);
        fread(pCell->phenotype.molecular.fraction_released_at_death.data(), sizeof(double), m, fp);
        fread(pCell->phenotype.molecular.fraction_transferred_when_ingested.data(), sizeof(double), m, fp);
        
        // Interactions
        fread(&(pCell->phenotype.cell_interactions.apoptotic_phagocytosis_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_interactions.necrotic_phagocytosis_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_interactions.other_dead_phagocytosis_rate), sizeof(double), 1, fp);
        fread(pCell->phenotype.cell_interactions.live_phagocytosis_rates.data(), sizeof(double), n, fp);
        fread(pCell->phenotype.cell_interactions.attack_rates.data(), sizeof(double), n, fp);
        fread(pCell->phenotype.cell_interactions.immunogenicities.data(), sizeof(double), n, fp);
        
        // attack_target (stored as ID, need to resolve later)
        fread(&dTemp, sizeof(double), 1, fp);
        int attack_target_id = (int)dTemp;
        attack_target_ids.push_back(attack_target_id);
        
        fread(&(pCell->phenotype.cell_interactions.attack_damage_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_interactions.attack_duration), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_interactions.total_damage_delivered), sizeof(double), 1, fp);
        fread(pCell->phenotype.cell_interactions.fusion_rates.data(), sizeof(double), n, fp);
        
        // Transformations
        fread(pCell->phenotype.cell_transformations.transformation_rates.data(), sizeof(double), n, fp);
        
        // Asymmetric division
        fread(pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.data(), sizeof(double), n, fp);
        
        // Cell integrity
        fread(&(pCell->phenotype.cell_integrity.damage), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_integrity.damage_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.cell_integrity.damage_repair_rate), sizeof(double), 1, fp);
        
        // Custom scalar variables
        for (int j = 0; j < pCell->custom_data.variables.size(); j++)
        {
            fread(&(pCell->custom_data.variables[j].value), sizeof(double), 1, fp);
        }
        
        // Custom vector variables
        for (int j = 0; j < pCell->custom_data.vector_variables.size(); j++)
        {
            int size_temp = pCell->custom_data.vector_variables[j].value.size();
            fread(pCell->custom_data.vector_variables[j].value.data(), sizeof(double), size_temp, fp);
        }
        
        // Add cell to the simulation
        (*all_cells).push_back(pCell);
    }
    
    fclose(fp);
    
    // Resolve attack target pointers
    for (int i = 0; i < (*all_cells).size(); i++)
    {
        if (attack_target_ids[i] >= 0)
        {
            // Find cell with matching ID
            Cell* pTarget = NULL;
            for (int j = 0; j < (*all_cells).size(); j++)
            {
                if ((*all_cells)[j]->ID == attack_target_ids[i])
                {
                    pTarget = (*all_cells)[j];
                    break;
                }
            }
            (*all_cells)[i]->phenotype.cell_interactions.pAttackTarget = pTarget;
        }
        else
        {
            (*all_cells)[i]->phenotype.cell_interactions.pAttackTarget = NULL;
        }
    }
    
    std::cout << "Loaded " << number_of_data_entries << " cells from " << filename << std::endl;
}

};