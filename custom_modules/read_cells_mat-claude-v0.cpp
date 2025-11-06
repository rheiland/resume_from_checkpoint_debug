void read_PhysiCell_cells_from_matlab_v2(std::string filename, Microenvironment& M)
{
    // Get number of substrates, cell types, death models
    static int m = microenvironment.number_of_densities();
    static int n = cell_definition_indices_by_name.size();
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
    delete_all_cells();
    
    double dTemp;
    double position[3];
    double orientation[3];
    double velocity[3];
    double migration_bias_direction[3];
    double motility_vector[3];
    
    // Store attack target IDs for later resolution
    std::vector<int> attack_target_ids;
    
    // Read each cell
    for (int i = 0; i < number_of_data_entries; i++)
    {
        Cell* pCell = create_cell();
        
        // ID
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->ID = (int)dTemp;
        
        // position
        fread(position, sizeof(double), 3, fp);
        pCell->position[0] = position[0];
        pCell->position[1] = position[1];
        pCell->position[2] = position[2];
        
        // total_volume
        fread(&(pCell->phenotype.volume.total), sizeof(double), 1, fp);
        
        // cell_type
        fread(&dTemp, sizeof(double), 1, fp);
        pCell->type = (int)dTemp;
        
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
