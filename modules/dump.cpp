void dump_cell_mat_vars(int number_of_data_entries )
{
    double dTemp;

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
		std::cout << dTemp << std::endl; 
		// name = "position";    NOTE very different syntax for writing vectors!
        std::cout <<  pCell->position.data()  << std::endl;
		// name = "total_volume"; 
		std::cout <<  pCell->phenotype.volume.total << std::endl; 
		// name = "cell_type"; 
		dTemp = (double) pCell->type;
		std::cout <<  dTemp << std::endl; 
		// name = "cycle_model"; 
		dTemp = (double) pCell->phenotype.cycle.model().code; 
		std::cout <<  dTemp << std::endl; // cycle model 
		// name = "current_phase"; 
		dTemp = (double) pCell->phenotype.cycle.current_phase().code; 
		std::cout <<  dTemp << std::endl; // cycle model 
		// name = "elapsed_time_in_phase"; 
		std::cout <<  pCell->phenotype.cycle.data.elapsed_time_in_phase << std::endl; 
		// name = "nuclear_volume"; 
		std::cout <<  pCell->phenotype.volume.nuclear << std::endl;   
		// name = "cytoplasmic_volume"; 
		std::cout <<  pCell->phenotype.volume.cytoplasmic << std::endl;
		// name = "fluid_fraction"; 
		std::cout <<  pCell->phenotype.volume.fluid_fraction << std::endl;
		// name = "calcified_fraction"; 
		std::cout <<  pCell->phenotype.volume.calcified_fraction << std::endl; 
		// name = "orientation"; 
		std::cout <<  pCell->state.orientation.data()  << std::endl; 
		// name = "polarity"; 
		std::cout <<  pCell->phenotype.geometry.polarity << std::endl; 

 /* state variables to save */ 
// state
		// name = "velocity"; 
		std::cout <<  pCell->velocity.data()  << std::endl; 
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
		// name = "current_cycle_phase_exit_rate"; 
		int phase_index = pCell->phenotype.cycle.data.current_phase_index; 
		std::cout <<  pCell->phenotype.cycle.data.exit_rate(phase_index) << std::endl; 
		// name = "elapsed_time_in_phase"; 
		std::cout <<  pCell->phenotype.cycle.data.elapsed_time_in_phase << std::endl; 

// death 
  // live or dead state // 1 
		// name = "dead"; 
		dTemp = (double) pCell->phenotype.death.dead; 
		std::cout <<  dTemp << std::endl; 
		// name = "current_death_model"; // 
		dTemp = (double) pCell->phenotype.death.current_death_model_index; 
		std::cout <<  dTemp << std::endl; 
		// name = "death_rates"; 
		std::cout <<  pCell->phenotype.death.rates.data() , sizeof(double) , nd , fp ); 
		
	// volume ()
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
     // radius //1 
		// name = "radius"; 
		std::cout <<  pCell->phenotype.geometry.radius << std::endl; 
		// name = "nuclear_radius"; 
		std::cout <<  pCell->phenotype.geometry.nuclear_radius << std::endl; 
		// name = "surface_area"; 
		std::cout <<  pCell->phenotype.geometry.surface_area << std::endl; 

  // mechanics 
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
		std::cout <<  pCell->phenotype.mechanics.cell_adhesion_affinities.data() << std::endl;  // , sizeof(double) , n , fp ); 
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
 		// name = "is_motile"; 
		dTemp = (double) pCell->phenotype.motility.is_motile; 
		std::cout <<  dTemp << std::endl; 
 		// name = "persistence_time"; 
		std::cout <<  pCell->phenotype.motility.persistence_time << std::endl; 
 		// name = "migration_speed"; 
		std::cout <<  pCell->phenotype.motility.migration_speed << std::endl; 
 		// name = "migration_bias_direction"; 
		std::cout <<  pCell->phenotype.motility.migration_bias_direction.data()  << std::endl; 
 		// name = "migration_bias"; 
		std::cout <<  pCell->phenotype.motility.migration_bias << std::endl; 
 		// name = "motility_vector"; 
		std::cout <<  pCell->phenotype.motility.motility_vector.data()  << std::endl; 
 		// name = "chemotaxis_index"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_index; 
		std::cout <<  dTemp << std::endl; 
 		// name = "chemotaxis_direction"; 
		dTemp = (double) pCell->phenotype.motility.chemotaxis_direction; 
		std::cout <<  dTemp << std::endl; 
 		// name = "chemotactic_sensitivities"; 
		std::cout <<  pCell->phenotype.motility.chemotactic_sensitivities.data() << std::endl; // , sizeof(double) , m , fp ); 

// secretion 
 		// name = "secretion_rates"; 
		std::cout <<  pCell->phenotype.secretion.secretion_rates.data() << std::endl; // , sizeof(double) , m , fp ); 
	 	// name = "uptake_rates"; 
		std::cout <<  pCell->phenotype.secretion.uptake_rates.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "saturation_densities"; 
		std::cout <<  pCell->phenotype.secretion.saturation_densities.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "net_export_rates"; 
		std::cout <<  pCell->phenotype.secretion.net_export_rates.data() << std::endl; // , sizeof(double) , m , fp ); 

// molecular 
 		// name = "internalized_total_substrates"; 
		std::cout <<  pCell->phenotype.molecular.internalized_total_substrates.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "fraction_released_at_death"; 
		std::cout <<  pCell->phenotype.molecular.fraction_released_at_death.data() << std::endl; // , sizeof(double) , m , fp ); 
 		// name = "fraction_transferred_when_ingested"; 
		std::cout <<  pCell->phenotype.molecular.fraction_transferred_when_ingested.data() << std::endl; // , sizeof(double) , m , fp ); 

// interactions 
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
		std::cout <<  &(dTemp) , sizeof(double) , 1 , fp ); 
 		// name = "attack_damage_rate"; 
		std::cout <<  pCell->phenotype.cell_interactions.attack_damage_rate << std::endl; 
 		// name = "attack_duration"; 
		std::cout <<  pCell->phenotype.cell_interactions.attack_duration << std::endl; 
 		// name = "total_damage_delivered"; 
		std::cout <<  pCell->phenotype.cell_interactions.total_damage_delivered << std::endl; 

 		// name = "fusion_rates"; 
		std::cout <<  pCell->phenotype.cell_interactions.fusion_rates.data() << std::endl; // , sizeof(double) , n , fp ); 

// transformations 
  		// name = "transformation_rates"; 
		std::cout <<  pCell->phenotype.cell_transformations.transformation_rates.data() << std::endl; // , sizeof(double) , n , fp ); 

// asymmetric division
		// name = "asymmetric_division_rate"; 
		std::cout <<  pCell->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.data() << std::endl; // , sizeof(double) , n, fp );

	// cell integrity 
 		// name = "damage"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage << std::endl; 
 		// name = "damage_rate"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage_rate << std::endl; 
 		// name = "damage_repair_rate"; 
		std::cout <<  pCell->phenotype.cell_integrity.damage_repair_rate << std::endl; 

// custom 
		// custom scalar variables 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.variables.size(); j++ )
		{ std::cout <<  pCell->custom_data.variables[j].value << std::endl; }

		// custom vector variables 
		for( int j=0 ; j < (*all_cells)[0]->custom_data.vector_variables.size(); j++ )
		{
			int size_temp = pCell->custom_data.vector_variables[j].value.size(); 
			std::cout <<  pCell->custom_data.vector_variables[j].value.data() << std::endl; // , sizeof(double) , size_temp , fp );
		}
	}
    std::cout << "---------- END dump_cell_mat_vars:\n";
}
