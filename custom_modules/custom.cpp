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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

//------------------------------------------------------------
void dump_cells_mat(std::string filename, Microenvironment& M, bool create_cells)
{
    std::cout << "\n------- read START dump_cells_mat (custom.cpp) ----------\n";

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
    int cell_ID, cell_type;
    double cell_vol;
    double orientation[3];
    double velocity[3];
    double migration_bias_direction[3];
    double motility_vector[3];
    double params[10];
    double *dptr;
    
    // Store attack target IDs for later resolution
    std::vector<int> attack_target_ids;
    
    Cell_Definition* pCD; 
    Cell* pCell;
    // Read each cell
    // std::cout << "creating cells..." << std::endl;
    std::cout << "reading cell data..." << std::endl;
    for (int i = 0; i < number_of_data_entries; i++)
    {
        std::cout << "\n ------  reading cell # " << i << std::endl;
        
        fread(&dTemp, sizeof(double), 1, fp);
        cell_ID = (int)dTemp;
        std::cout << "ID= " << cell_ID << std::endl; 
        
        fread(position, sizeof(double), 3, fp);
        std::cout << "pos= " << position[0]<<", "<< position[1] << std::endl;
        
        fread(&cell_vol, sizeof(double), 1, fp);
        std::cout << "total vol= " << cell_vol  << std::endl;
        
        fread(&dTemp, sizeof(double), 1, fp);
        cell_type = (int)dTemp;
        std::cout << "type= " << cell_type  << std::endl;

        if (create_cells)
        {
            // Cell_Definition* pCD = cell_definitions_by_type[cell_type]; 
            pCD = cell_definitions_by_type[cell_type];    // rwh: better?
            pCell = create_cell( *pCD );

            pCell->ID = cell_ID;
            pCell->type = cell_type;
            pCell->phenotype.volume.total = cell_vol;

            pCell->assign_position(position[0], position[1], position[2]);
            pCell->update_voxel_index();
        }
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.cycle.model().code= " << int(dTemp)  << std::endl;
        if (create_cells)
            pCell->phenotype.cycle.model().code = int(dTemp);
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.cycle.current_phase().code = current_phase(?) = " << int(dTemp)  << std::endl;
        if (create_cells)
            pCell->phenotype.cycle.current_phase().code = int(dTemp);
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.cycle.data.elapsed_time_in_phase= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.cycle.data.elapsed_time_in_phase = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.volume.nuclear= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.volume.nuclear = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.volume.cytoplasmic= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.volume.cytoplasmic = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.volume.fluid_fraction= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.volume.fluid_fraction = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.volume.calcified_fraction= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.volume.calcified_fraction = dTemp;
        
        fread(orientation, sizeof(double), 3, fp);
        std::cout << "orientation= " << orientation[0]<<", "<<orientation[1]<<", " << orientation[2]  << std::endl;
        if (create_cells)
        {
            pCell->state.orientation[0] = orientation[0];
            pCell->state.orientation[1] = orientation[1];
            pCell->state.orientation[2] = orientation[2];
        }
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.geometry.polarity= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.geometry.polarity = dTemp;
        

        std::cout << "------ state:\n";
        fread(velocity, sizeof(double), 3, fp);
        std::cout << "velocity= " << velocity[0]<<", "<<velocity[1]<<", " << velocity[2]  << std::endl;
        if (create_cells)
        {
            pCell->velocity[0] = velocity[0];
            pCell->velocity[1] = velocity[1];
            pCell->velocity[2] = velocity[2];
        }
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "state.simple_pressure= " << dTemp  << std::endl;
        if (create_cells)
            pCell->state.simple_pressure = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "state.number_of_nuclei= " << (int)dTemp  << std::endl;
        if (create_cells)
            pCell->state.number_of_nuclei = (int)dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "state.total_attack_time= " << dTemp  << std::endl;
        if (create_cells)
            pCell->state.total_attack_time = dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "state.contact_with_basement_membrane= " << bool(dTemp)  << std::endl;
        if (create_cells)
            pCell->state.contact_with_basement_membrane = (bool)dTemp;
        

        std::cout << "------ cycle:\n";
        // current_cycle_phase_exit_rate
        // rwh - what do I use??
        // int phase_index = pCell->phenotype.cycle.data.current_phase_index;
        // int phase_index = pCell->phenotype.cycle.current_phase().code;
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.cycle.data.exit_rate(phase_index)= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.cycle.data.exit_rate(phase_index) = dTemp;
        
        // elapsed_time_in_phase (duplicate: rf. https://github.com/MathCancer/PhysiCell/pull/380)
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << "phenotype.cycle.data.elapsed_time_in_phase= " << dTemp  << std::endl;
        if (create_cells)
            pCell->phenotype.cycle.data.elapsed_time_in_phase = dTemp;
        

        std::cout << "------ death:\n";
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << " phenotype.death.dead= " << bool(dTemp)  << std::endl;
        if (create_cells)
            pCell->phenotype.death.dead = (bool)dTemp;
        
        fread(&dTemp, sizeof(double), 1, fp);
        std::cout << " phenotype.death.current_death_model_index= " << int(dTemp)  << std::endl;
        if (create_cells)
            pCell->phenotype.death.current_death_model_index = (int)dTemp;
        
        // Set nd from first cell if not set
        // if (nd == 0)
        // {
        //     nd = pCell->phenotype.death.rates.size();
        // }
        int ndeath = pCell->phenotype.death.rates.size();
        // std::cout << " phenotype.death.rates.size()= " <<  pCell->phenotype.death.rates.size()  << std::endl;
        std::cout << " phenotype.death.rates.size()= " <<  ndeath  << std::endl;
        fread(dptr, sizeof(double), ndeath, fp);
        for (idx=0; idx < ndeath); idx++)
        {
            fread(&dTemp, sizeof(double), 1, fp);
            std::cout << " phenotype.death.rates[" << idx << "] = " << dTemp  << std::endl;
            if (create_cells)
                pCell->phenotype.death.rates[idx] = dTemp;
        }
        
        // death_rates
        //rwh
        // fread(pCell->phenotype.death.rates.data(), sizeof(double), nd, fp);
        
        if (create_cells)
        {
            // pCell->phenotype.death.rates.data()
            for (int idx=0; idx<ndeath; idx++)
                pCell->phenotype.death.rates[idx] = dptr[idx]
        }

        // Volume rates
        fread(params, sizeof(double), 7, fp);
        std::cout << "phenotype.volume.cytoplasmic_biomass_change_rate= " << params[0] << std::endl;
        std::cout << "phenotype.volume.nuclear_biomass_change_rate= " << params[1] << std::endl;
        std::cout << "phenotype.volume.fluid_change_rate= " << params[2] << std::endl;
        std::cout << "phenotype.volume.calcification_rate= " << params[3] << std::endl;
        std::cout << "phenotype.volume.target_solid_cytoplasmic= " << params[4] << std::endl;
        std::cout << "phenotype.volume.target_solid_nuclear= " << params[5] << std::endl;
        std::cout << "phenotype.volume.target_fluid_fraction= " << params[6] << std::endl;
        if (create_cells)
        {
            pCell->phenotype.volume.cytoplasmic_biomass_change_rate = params[0];
            pCell->phenotype.volume.nuclear_biomass_change_rate = params[1];
            pCell->phenotype.volume.fluid_change_rate = params[2];
            pCell->phenotype.volume.calcification_rate = params[3];
            pCell->phenotype.volume.target_solid_cytoplasmic = params[4];
            pCell->phenotype.volume.target_solid_nuclear = params[5];
            pCell->phenotype.volume.target_fluid_fraction = params[6];
        }

        
        // Geometry
        fread(params, sizeof(double), 3, fp);
        std::cout << "phenotype.geometry.radius= " << params[0] << std::endl;
        std::cout << "phenotype.geometry.nuclear_radius= " << params[1] << std::endl;
        std::cout << "phenotype.geometry.surface_area= " << params[2] << std::endl;
        if (create_cells)
        {
            pCell->phenotype.geometry.radius= params[0];
            pCell->phenotype.geometry.nuclear_radius= params[1];
            pCell->pCell->phenotype.geometry.surface_area = params[2];
        }
        
        // Mechanics
        fread(params, sizeof(double), 4, fp);
        std::cout << "phenotype.mechanics.cell_cell_adhesion_strength= " << params[0] << std::endl;
        std::cout << "phenotype.mechanics.cell_BM_adhesion_strength= " << params[1] << std::endl;
        std::cout << "phenotype.mechanics.cell_cell_repulsion_strength= " << params[2] << std::endl;
        std::cout << "phenotype.mechanics.cell_BM_repulsion_strength= " << params[3] << std::endl;
        if (create_cells)
        {
            pCell->phenotype.mechanics.cell_cell_adhesion_strength= params[0];
            pCell->phenotype.mechanics.cell_BM_adhesion_strength= params[1];
            pCell->phenotype.mechanics.cell_cell_repulsion_strength = params[2];
            pCell->phenotype.mechanics.cell_BM_repulsion_strength = params[3];
        }

        fread(pCell->phenotype.mechanics.cell_adhesion_affinities.data(), sizeof(double), n, fp);  // NOTE

        fread(&(pCell->phenotype.mechanics.relative_maximum_adhesion_distance), sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->phenotype.mechanics.maximum_number_of_attachments = (int)dTemp;
        fread(&(pCell->phenotype.mechanics.attachment_elastic_constant), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.attachment_rate), sizeof(double), 1, fp);
        fread(&(pCell->phenotype.mechanics.detachment_rate), sizeof(double), 1, fp);
        
        // Motility
        fread(&dTemp, sizeof(double), 1, fp);
        // pCell->phenotype.motility.is_motile = (bool)dTemp;
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        fread(migration_bias_direction, sizeof(double), 3, fp);
        // pCell->phenotype.motility.migration_bias_direction[0] = migration_bias_direction[0];
        // pCell->phenotype.motility.migration_bias_direction[1] = migration_bias_direction[1];
        // pCell->phenotype.motility.migration_bias_direction[2] = migration_bias_direction[2];
        fread(&dTemp, sizeof(double), 1, fp);
        fread(motility_vector, sizeof(double), 3, fp);
        // pCell->phenotype.motility.motility_vector[0] = motility_vector[0];
        // pCell->phenotype.motility.motility_vector[1] = motility_vector[1];
        // pCell->phenotype.motility.motility_vector[2] = motility_vector[2];
        fread(&dTemp, sizeof(double), 1, fp);
        // pCell->phenotype.motility.chemotaxis_index = (int)dTemp;
        fread(&dTemp, sizeof(double), 1, fp);
        // pCell->phenotype.motility.chemotaxis_direction = (int)dTemp;

        // rwh - fix!
        // fread(pCell->phenotype.motility.chemotactic_sensitivities.data(), sizeof(double), m, fp);
        
        // // Secretion
        // fread(pCell->phenotype.secretion.secretion_rates.data(), sizeof(double), m, fp);
        // fread(pCell->phenotype.secretion.uptake_rates.data(), sizeof(double), m, fp);
        // fread(pCell->phenotype.secretion.saturation_densities.data(), sizeof(double), m, fp);
        // fread(pCell->phenotype.secretion.net_export_rates.data(), sizeof(double), m, fp);
        
        // // Molecular
        // fread(pCell->phenotype.molecular.internalized_total_substrates.data(), sizeof(double), m, fp);
        // fread(pCell->phenotype.molecular.fraction_released_at_death.data(), sizeof(double), m, fp);
        // fread(pCell->phenotype.molecular.fraction_transferred_when_ingested.data(), sizeof(double), m, fp);
        
        // Interactions
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        // fread(pCell->phenotype.cell_interactions.live_phagocytosis_rates.data(), sizeof(double), n, fp);
        // fread(pCell->phenotype.cell_interactions.attack_rates.data(), sizeof(double), n, fp);
        // fread(pCell->phenotype.cell_interactions.immunogenicities.data(), sizeof(double), n, fp);
        
        // attack_target (stored as ID, need to resolve later)
        fread(&dTemp, sizeof(double), 1, fp);
        int attack_target_id = (int)dTemp;
        attack_target_ids.push_back(attack_target_id);
        
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        // fread(pCell->phenotype.cell_interactions.fusion_rates.data(), sizeof(double), n, fp);
        
        // Transformations
        // fread(pCell->phenotype.cell_transformations.transformation_rates.data(), sizeof(double), n, fp);
        
        // Asymmetric division
        // fread(pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.data(), sizeof(double), n, fp);
        
        // Cell integrity
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        fread(&dTemp, sizeof(double), 1, fp);
        
        // Custom scalar variables
        // for (int j = 0; j < pCell->custom_data.variables.size(); j++)
        // {
        //     fread(&dTemp, sizeof(double), 1, fp);
        // }
        
        // Custom vector variables
        // for (int j = 0; j < pCell->custom_data.vector_variables.size(); j++)
        // {
        //     int size_temp = pCell->custom_data.vector_variables[j].value.size();
        //     fread(pCell->custom_data.vector_variables[j].value.data(), sizeof(double), size_temp, fp);
        // }
        
        // Add cell to the simulation
        // (*all_cells).push_back(pCell);
    }
    
    fclose(fp);
    
    // Resolve attack target pointers
    // for (int i = 0; i < (*all_cells).size(); i++)
    // {
    //     if (attack_target_ids[i] >= 0)
    //     {
    //         // Find cell with matching ID
    //         Cell* pTarget = NULL;
    //         for (int j = 0; j < (*all_cells).size(); j++)
    //         {
    //             if ((*all_cells)[j]->ID == attack_target_ids[i])
    //             {
    //                 pTarget = (*all_cells)[j];
    //                 break;
    //             }
    //         }
    //         (*all_cells)[i]->phenotype.cell_interactions.pAttackTarget = pTarget;
    //     }
    //     else
    //     {
    //         (*all_cells)[i]->phenotype.cell_interactions.pAttackTarget = NULL;
    //     }
    // }
    
    std::cout << "read " << number_of_data_entries << " cells from " << filename << std::endl;
    std::cout << "------- END dump_cells_mat ----------\n";
}