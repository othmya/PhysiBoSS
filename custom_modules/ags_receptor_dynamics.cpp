/*
	
	ags_receptor_dynamics.cpp

	each drug that enters does so with its own mechanism
*/

#include "./ags_receptor_dynamics.h" 
#include "./ags_boolean_model_interface.h" 

using namespace PhysiCell; 

Submodel_Information ags_receptor_info;

void ags_receptor_model_setup()
{
    ags_receptor_info.name = "AGS model"; 
	ags_receptor_info.version = "0.0.1";
	
    ags_receptor_info.main_function = ags_receptor_model; 

	// what custom data do I need?
	ags_receptor_info.cell_variables.push_back( "activation_threshold" );

	// dA (Densities & kinetics)

	ags_receptor_info.cell_variables.push_back( "Initial_I_dA" );
	ags_receptor_info.cell_variables.push_back( "Initial_E_dA" );
	ags_receptor_info.cell_variables.push_back( "k_dA" );
	ags_receptor_info.cell_variables.push_back( "Rc_dA" );
	ags_receptor_info.cell_variables.push_back( "Rcb_dA" );
	ags_receptor_info.cell_variables.push_back( "dA_Km" );
	ags_receptor_info.cell_variables.push_back( "dA_vmax" );

	// dB (Densities & kinetics)

	ags_receptor_info.cell_variables.push_back( "Initial_I_dB" );
	ags_receptor_info.cell_variables.push_back( "Initial_E_dB" );
	ags_receptor_info.cell_variables.push_back( "k_dB" );
	ags_receptor_info.cell_variables.push_back( "Rc_dB" );
	ags_receptor_info.cell_variables.push_back( "Rcb_dB" );
	ags_receptor_info.cell_variables.push_back( "dB_Km" );
	ags_receptor_info.cell_variables.push_back( "dB_vmax" );


	ags_receptor_info.register_model();

	return;
}

// Function that actually does the computation of kinetics of the movement of drugs from or to the agents
// Follows the same code I employ for the transport models

void ags_receptor_model( Cell* pCell, Phenotype& phenotype, double dt )
{

	if( phenotype.death.dead == true )
	{ return; } 

	// Fetch microenvironment variables
	static int dA_ix = microenvironment.find_density_index( "dA" );
	static int dB_ix = microenvironment.find_density_index( "dB" );


	// Internal, Extenral densities
	static int E_dA_ix = pCell->custom_data.find_variable_index( "E_dA" );
	static int E_dA_near_ix = pCell->custom_data.find_variable_index( "E_dA_near" );
	static int I_dA_ix = pCell->custom_data.find_variable_index( "I_dA" );
	
	static int E_dB_ix = pCell->custom_data.find_variable_index( "E_dB" );
	static int E_dB_near_ix = pCell->custom_data.find_variable_index( "E_dB_near" );
	static int I_dB_ix = pCell->custom_data.find_variable_index( "I_dB" );

	static int total_I_dA_ix = pCell->custom_data.find_variable_index( "total_I_dA" );
	static int total_E_dA_ix = pCell->custom_data.find_variable_index( "total_E_dA" );
	static int total_dA_ix = pCell->custom_data.find_variable_index( "total_dA" );
	static int total_I_dB_ix = pCell->custom_data.find_variable_index( "total_I_dB" );
	static int total_E_dB_ix = pCell->custom_data.find_variable_index( "total_E_dB" );
	static int total_dB_ix = pCell->custom_data.find_variable_index( "total_dB" );

	// // Receptor densities
    // static int nRc_dA_ix = pCell->custom_data.find_variable_index( "Rc_dA" ); 
	// static int nRcb_dA_ix = pCell->custom_data.find_variable_index( "Rcb_dA" );
	// static int nRc_total_dA_ix = pCell->custom_data.find_variable_index( "total_Rc_dA" );
	
    // static int nRc_dB_ix = pCell->custom_data.find_variable_index( "Rc_dB" ); 
	// static int nRcb_dB_ix = pCell->custom_data.find_variable_index( "Rcb_dB" );
	// static int nRc_total_dB_ix = pCell->custom_data.find_variable_index( "total_Rc_dB" );

	// // Kinetics

    // static int dA_bind_ix = pCell->custom_data.find_variable_index( "dA_k1" );
    // static int dA_recycle_ix = pCell->custom_data.find_variable_index( "dA_kn1" );
    // static int dA_endocyte_ix = pCell->custom_data.find_variable_index( "dA_k2" );
    // static int dA_internal_bind_ix = pCell->custom_data.find_variable_index( "dA_kn2" );

	// static int dA_Km_ix = pCell->custom_data.find_variable_index( "dA_Km" );
	// static int dA_vmax_ix = pCell->custom_data.find_variable_index( "dA_vmax" );

	static int k_dA_ix = pCell->custom_data.find_variable_index( "k_dA" );
	static int k_dB_ix = pCell->custom_data.find_variable_index( "k_dB" );

	// static int dB_bind_ix = pCell->custom_data.find_variable_index( "dB_k1" );
	// static int dB_recycle_ix = pCell->custom_data.find_variable_index( "dB_kn1" );
	// static int dB_endocyte_ix = pCell->custom_data.find_variable_index( "dB_k2" );
	// static int dB_internal_bind_ix = pCell->custom_data.find_variable_index( "dB_kn2" );

	// static int dB_Km_ix = pCell->custom_data.find_variable_index( "dB_Km" );
	// static int dB_vmax_ix = pCell->custom_data.find_variable_index( "dB_vmax" );

	// drug fluxes

	static int dA_flux_ix = pCell->custom_data.find_variable_index( "dA_flux" );
	static int dB_flux_ix = pCell->custom_data.find_variable_index( "dB_flux" );

	// Concentration gradient

	static int D_dA_ix = pCell->custom_data.find_variable_index( "D_dA" );
	static int D_dB_ix = pCell->custom_data.find_variable_index( "D_dB" );

	// static double R_binding_rate = pCell->custom_data[nR_bind]; Why does he include this here?

	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;


	// Constants 

	double DC_dA = microenvironment.diffusion_coefficients[dA_ix];
	double V_cell = pCell->phenotype.volume.total;
	double R_cell = pCell->phenotype.geometry.radius;
	double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Receptor concentrations
	// double Rc_dA = pCell->custom_data.variables[nRc_dA_ix].value;	  // [E] Receptor alone (mM)
	// double Rcb_dA = pCell->custom_data.variables[nRcb_dA_ix].value; // [ES] Bound receptor (mM)
	// double Rc_dB = pCell->custom_data.variables[nRc_dB_ix].value;	 
	// double Rcb_dB = pCell->custom_data.variables[nRcb_dB_ix].value; 

	// Substrate concentrations
	double E_dA_near = pCell->nearest_density_vector()[dA_ix];	 // A density (mM)
	double E_dB_near = pCell->nearest_density_vector()[dB_ix];				
	double I_dA = pCell->phenotype.molecular.internalized_total_substrates[dA_ix] / V_cell; // Convert to density (mM)
	double I_dB = pCell->phenotype.molecular.internalized_total_substrates[dB_ix] / V_cell;

	double D_dA = I_dA - E_dA_near;
	double D_dB = I_dB - E_dB_near;
	/*

		Transport model -- FDC

	*/

	// Obtaining the total net internal and external amounts through iterating over the cells

	// double total_E_dA;
	// double total_E_dB;
	// double total_I_dA;
	// double total_I_dB;

	// std::vector<double> I_dA_vector(microenvironment.number_of_voxels());
	// std::vector<double> I_dB_vector(microenvironment.number_of_voxels());

	// Counting net internal amount
	// #pragma omp parallel for
	// for (int i = 0; i < (*all_cells).size(); i++)
	// {
	// 	Cell *pC = (*all_cells)[i];
	// 	I_dA_vector[i] = pC->phenotype.molecular.internalized_total_substrates[dA_ix];
	// 	total_I_dA += pC->phenotype.molecular.internalized_total_substrates[dA_ix];

	// 	I_dB_vector[i] = pC->phenotype.molecular.internalized_total_substrates[dB_ix];
	// 	total_I_dB += pC->phenotype.molecular.internalized_total_substrates[dB_ix];
	// }

	// double I_dA_sum = std::accumulate(I_dA_vector.begin(), I_dA_vector.end(), 0.0);
	// double I_dB_sum = std::accumulate(I_dB_vector.begin(), I_dB_vector.end(), 0.0);


	// std::vector<double> E_dA_vector(microenvironment.number_of_voxels());
	// std::vector<double> E_dB_vector(microenvironment.number_of_voxels());
	
	// #pragma omp parallel for
	// for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	// {
	// 	E_dA_vector[n] = (microenvironment.density_vector(n)[dA_ix]);
	// 	E_dB_vector[n] = (microenvironment.density_vector(n)[dB_ix]);

	// 	total_E_dA += (microenvironment.density_vector(n)[dA_ix]) * V_voxel; // A net amount of dA
	// 	total_E_dB += (microenvironment.density_vector(n)[dB_ix]) * V_voxel; // A net amount of dB
	// }

	// double E_dA_sum = std::accumulate(E_dA_vector.begin(), E_dA_vector.end(), 0.0);
	// double E_dA_mean = E_dA_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	// double E_dB_sum = std::accumulate(E_dB_vector.begin(), E_dB_vector.end(), 0.0);
	// double E_dB_mean = E_dB_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels



	// replace hard-coded initial values at each time-step with the new values

	// // Substrate concentrations
	
	// pCell->custom_data.variables[E_dA_near_ix].value = E_dA_near; // Mean of near_density_vector
	// // pCell->custom_data.variables[E_dA_ix].value = E_dA_mean;		 // Mean of all voxels
	// pCell->custom_data.variables[I_dA_ix].value = I_dA;			 // Density

	// pCell->custom_data.variables[E_dB_near_ix].value = E_dB_near;
	// // pCell->custom_data.variables[E_dB_ix].value = E_dB_mean;
	// pCell->custom_data.variables[I_dB_ix].value = I_dB;

	// pCell->custom_data.variables[D_dA_ix].value = D_dA;	 // Concentration gradient
	// pCell->custom_data.variables[D_dB_ix].value = D_dB;


	
	// fdc reduced model -- can be obtained from kinetic constants or from Km / vmax	
	// double dA_vmax_f = parameters.doubles("dA_vmax_f");
	// double dA_vmax_r = parameters.doubles("dA_vmax_r");
	// double dA_Km1 = parameters.doubles("dA_Km1");
	// double dA_Km2 = parameters.doubles("dA_Km2");

	// double dA_flux = -((dA_vmax_f * E_dA_near) / dA_Km1 - (dA_vmax_r * I_dA) / dA_Km2) / (1 + (I_dA / dA_Km2) + (E_dA_near / dA_Km1));
	// dA_flux *= V_cell;
	
	// double dB_vmax_f = parameters.doubles("dB_vmax_f");
	// double dB_vmax_r = parameters.doubles("dB_vmax_r");
	// double dB_Km1 = parameters.doubles("dB_Km1");
	// double dB_Km2 = parameters.doubles("dB_Km2");

	// double dB_flux = -((dB_vmax_f * E_dB_near) / dB_Km1 - (dB_vmax_r * I_dB) / dB_Km2) / (1 + (I_dB / dB_Km2) + (E_dB_near / dB_Km1));
	// dB_flux *= V_cell;

	// double dA_monod = vf_max * (I_dA / (Km_1 + I_dA)); // Using the "forward" (entry) kinetic values
	// dA_monod *= diffusion_dt;

	// double dB_monod = vb_max * (I_dB / (Km_2 + I_dB)); // Using the "forward" (entry) kinetic values
	// dB_monod *= diffusion_dt;

	// Receptor concentrations
	// pCell->custom_data.variables[nRc_total_dA_ix].value = Rcb_dA + Rc_dA;
	// pCell->custom_data.variables[nRc_total_dB_ix].value = Rcb_dB + Rc_dB;

	// pCell->custom_data.variables[nRc_dA_ix].value = Rc_dA;
	// pCell->custom_data.variables[nRcb_dA_ix].value = Rcb_dA;
	// pCell->custom_data.variables[nRc_dB_ix].value = Rc_dB;
	// pCell->custom_data.variables[nRcb_dB_ix].value = Rcb_dB;


	/*

		Transport model -- Simple Diffusion

	*/

	double k_dA = parameters.doubles("k_dA");
	double dA_flux = k_dA * D_dA * A_cell;
	double adjusted_flux_dA = dA_flux;
	if ( dA_flux < 0.0 ) { // Uptake
		adjusted_flux_dA = - std::min( std::abs(dA_flux), std::abs(E_dA_near * V_voxel) );
	}
	else if ( dA_flux > 0.0 ) { // Secretion
		adjusted_flux_dA = std::min(dA_flux, std::abs(I_dA * V_cell ) );
	}

	double k_dB = parameters.doubles("k_dB");
	double dB_flux = k_dB * D_dB * A_cell;
	double adjusted_flux_dB = dB_flux;
	if ( dB_flux < 0.0 ) { // Uptake
		adjusted_flux_dB = - std::min( std::abs(dB_flux), std::abs(E_dB_near * V_voxel) );
	}
	else if ( dB_flux > 0.0 ) { // Secretion
		adjusted_flux_dB = std::min(dB_flux, std::abs(I_dB * V_cell ) );
	}




	// Fluxes
	// pCell->custom_data.variables[dA_flux_ix].value = adjusted_flux_dA;
	// pCell->custom_data.variables[dB_flux_ix].value = adjusted_flux_dB;


	// // pCell->custom_data.variables[dB_flux_ix].value = dB_flux;
	// // pCell->custom_data.variables[dA_flux_explicit_index].value = v_formation_explicit / diffusion_dt;
	
	// // Net total amounts
	// // pCell->custom_data.variables[total_E_dA_ix].value = E_dA_sum;			 
	// // pCell->custom_data.variables[total_I_dA_ix].value = I_dA_sum;			  
	// // pCell->custom_data.variables[total_dA_ix].value = total_E_dA + total_I_dA; 

	// // pCell->custom_data.variables[total_E_dB_ix].value = E_dB_sum;
	// // pCell->custom_data.variables[total_I_dB_ix].value = I_dB_sum;
	// // pCell->custom_data.variables[total_dB_ix].value = total_E_dB + total_I_dB;
	



	// // Mapping to net export rate of agent


	pCell->phenotype.secretion.net_export_rates[dA_ix] = adjusted_flux_dA;
	pCell->phenotype.secretion.net_export_rates[dB_ix] = adjusted_flux_dB;
	// pCell->phenotype.secretion.net_export_rates[dA_ix] = dA_flux; // If FDC model is employed

	// pCell->phenotype.secretion.net_export_rates[dB_ix] = dB_flux;	

	return;
}

void ags_receptor_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		if( pC->phenotype.death.dead == false )
		{ 
			ags_receptor_model( pC, pC->phenotype , dt );
		}
	}
	
	return;
}


// ETC code

	// double binding_component;
	// double recycling_component;
	// double endocytosis_component;
	// double internal_binding_component;

	// binding_component = diffusion_dt * dA_binding_rate * Rc_dA * E_dA_near;
	// if (binding_component > Rc_dA || binding_component > E_dA_near )
	// { binding_component = std::min(Rc_dA, E_dA_near); } // Get the limiting step 

	// Rcb_dA += binding_component;
	// Rc_dA -= binding_component;
	// // E_dA_near -= binding_component;
	// if(Rc_dA < 0.0) { Rc_dA = 0.0; }
	// if(E_dA_near < 0.0) { E_dA_near = 0.0; }

	// recycling_component = diffusion_dt * dA_recycling_rate * Rcb_dA;
	// if (recycling_component > Rcb_dA){ recycling_component = Rcb_dA; } 
	// Rc_dA += recycling_component;
	// Rcb_dA -= recycling_component;
	// // E_dA_near += recycling_component;
	// if(Rcb_dA < 0.0) { Rcb_dA = 0.0; }

	// endocytosis_component = diffusion_dt * dA_endocytosis_rate * Rcb_dA;
	// if (endocytosis_component > Rcb_dA){ endocytosis_component = Rcb_dA; } 
	// I_dA += endocytosis_component;
	// Rcb_dA -= endocytosis_component;
	// Rc_dA += endocytosis_component;
	// if(Rcb_dA < 0.0) { Rcb_dA = 0.0; }
	// // if(I_dA > E_dA_near) { endocytosis_component = 0.0; } // hard-code gradient

	// internal_binding_component = diffusion_dt * dA_internal_binding_rate * Rc_dA * I_dA;
	// if (internal_binding_component > Rcb_dA || internal_binding_component > I_dA )
	// { internal_binding_component = std::min(Rcb_dA, I_dA); } // Get the limiting step 

	// I_dA -= internal_binding_component;
	// Rcb_dA += internal_binding_component;
	// Rc_dA -= internal_binding_component;

	// if(Rc_dA < 0.0) { Rc_dA = 0.0; }
	// if(I_dA < 0.0) { I_dA = 0.0; }

	// double vf_max = dA_endocytosis_rate * (Rc_dA + Rcb_dA); // both are constant given that total receptor is constant
	// double vb_max = dA_recycling_rate * (Rc_dA + Rcb_dA);
	// double Km_1 = (dA_recycling_rate + dA_endocytosis_rate) / dA_binding_rate;
	// double Km_2 = (dA_recycling_rate + dA_endocytosis_rate) / dA_internal_binding_rate;



	// double v_formation_explicit = - endocytosis_component + internal_binding_component;
	// v_formation_explicit *= V_cell;