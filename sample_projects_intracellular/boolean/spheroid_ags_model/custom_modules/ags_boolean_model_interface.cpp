/*
 * ags_boolean_model_interface.cpp
 *
 */


#include "./ags_boolean_model_interface.h"
#include <math.h>

using namespace PhysiCell; 

Submodel_Information ags_bm_interface_info;

void ags_boolean_model_interface_setup()
{
    ags_bm_interface_info.name = "AGS Boolean model interface"; 
	ags_bm_interface_info.version = "0.0.1";
	
    ags_bm_interface_info.main_function= ags_bm_interface_main; 

	// These are just auxiliary variables to keep track of some BN nodes

    ags_bm_interface_info.cell_variables.push_back( "mek_node" );
    ags_bm_interface_info.cell_variables.push_back( "pi3k_node" );
    ags_bm_interface_info.cell_variables.push_back( "tak1_node" );
    ags_bm_interface_info.cell_variables.push_back( "akt_node" );

    ags_bm_interface_info.cell_variables.push_back( "anti_mek_node" );
    ags_bm_interface_info.cell_variables.push_back( "anti_pi3k_node" );
    ags_bm_interface_info.cell_variables.push_back( "anti_tak1_node" );
    ags_bm_interface_info.cell_variables.push_back( "anti_akt_node" );


	ags_bm_interface_info.register_model();
}


void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{

    if( pCell->phenotype.death.dead == true )
	{ return; }

    // Fetching variables from custom data
    
    static int dA_ix = microenvironment.find_density_index( "dA" );
	static int dB_ix = microenvironment.find_density_index( "dB" );

    static int mek_threshold_ix = pCell->custom_data.find_variable_index( "MEK_activation_threshold" );
    static int pi3k_threshold_ix = pCell->custom_data.find_variable_index( "PI3K_activation_threshold" );
    static int tak1_threshold_ix = pCell->custom_data.find_variable_index( "TAK1_activation_threshold" );
    static int akt_threshold_ix = pCell->custom_data.find_variable_index( "AKT_activation_threshold" );

    static int activation_threshold_ix = pCell->custom_data.find_variable_index( "activation_threshold" );
    
    // Connecting the Internal quantity of drug to the boolean model 

	double V_cell = pCell->phenotype.volume.total;

    double activation_threshold = pCell->custom_data[activation_threshold_ix];

    double I_dA = pCell->phenotype.molecular.internalized_total_substrates[dA_ix] / V_cell; // Convert to density (mM)
	double I_dB = pCell->phenotype.molecular.internalized_total_substrates[dB_ix] / V_cell;

    // std::cout << I_dA << std::endl;

    double GI50_dA = parameters.doubles("GI50_dA");
    double GI50_dB = parameters.doubles("GI50_dB");
    double dA_anti_activation_scaling = parameters.doubles("dA_antinode_activ_scaling");
    double dB_anti_activation_scaling = parameters.doubles("dB_antinode_activ_scaling");

    double dA_deactiv_prob = boolean_node_deactivation_prob(I_dA, dA_anti_activation_scaling, GI50_dA);
    double dB_deactiv_prob = boolean_node_deactivation_prob(I_dB, dB_anti_activation_scaling, GI50_dB);

    if (parameters.strings("node_affected_dA") != "none" && PhysiCell_globals.current_time >= parameters.doubles("dA_pulse_period") )
        {
            pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dA"), dA_deactiv_prob);
            // std::cout << "dA -- Drug density: " << I_dA << " Prob: " << dA_deactiv_prob << std::endl; 
        }

    if (parameters.strings("node_affected_dB") != "none"  && PhysiCell_globals.current_time >= parameters.doubles("dB_pulse_period") )
        {
            pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dB"), dB_deactiv_prob);
            // std::cout << "For dB -- Drug density: " << I_dB << " Prob: " << dB_deactiv_prob << std::endl; 
        }

    return;

}


void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{	
    static int dA_ix = microenvironment.find_density_index( "dA" );
	static int dB_ix = microenvironment.find_density_index( "dB" );

    static int dA_export_rate = pCell->custom_data.find_variable_index( "dA_net_production_rate" );

    static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
    
    // Getting the state of the boolean model readouts (Readout can be in the XML)

    bool casp37_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b1" );
    bool casp37_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b2" );
    bool FOXO = pCell->phenotype.intracellular->get_boolean_variable_value( "FOXO" );

    bool antisurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    bool antisurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    bool antisurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    double anti_w1 = 0.1;
    double anti_w2 = 0.2;
    double anti_w3 = 0.7;
    double S_anti = (anti_w1*antisurvival_b1) + (anti_w2 * antisurvival_b2) + (anti_w3 * antisurvival_b3);
    double S_anti_real = (anti_w1*casp37_b1) + (anti_w2 * casp37_b2) + (anti_w3 * FOXO);

    bool prosurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    bool prosurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    bool prosurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    bool cMYC = pCell->phenotype.intracellular->get_boolean_variable_value( "cMYC" );
    bool CCND_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b1" );
    bool CCND_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b2" );

    double pro_w1 = 0.2;
    double pro_w2 = 0.3;
    double pro_w3 = 0.5;
    double S_pro = (pro_w1*prosurvival_b1) + (pro_w2 * prosurvival_b2) + (pro_w3 * prosurvival_b3);
    double S_pro_real = (pro_w1*cMYC) + (pro_w2 * CCND_b1) + (pro_w3 * CCND_b2);

    // Connect output from model to actual cell variables

    double apoptosis_rate_basal = parameters.doubles("apoptosis_rate_basal");
    double maximum_apoptosis_rate = parameters.doubles("max_apoptosis_rate");
    double hill_coeff_apoptosis = parameters.doubles("hill_coeff_apoptosis");
    double K_half_apoptosis = parameters.doubles("K_half_apoptosis");


    double growth_rate_basal = parameters.doubles("basal_growth_rate");
    double p = pCell->state.simple_pressure; 
    double hill_coeff_pressure = parameters.doubles("hill_coeff_pressure");
    double pressure_half = parameters.doubles("pressure_half");
    double pressure_effect = pressure_effect_growth_rate(p, hill_coeff_pressure, pressure_half );    

    double hill_coeff_growth = parameters.doubles("hill_coeff_growth");
    double K_half_growth = parameters.doubles("K_half_growth");

    // std::cout << "Pressure is: " << p << std::endl;


    // Actual mapping

    pCell->phenotype.death.rates[apoptosis_model_index] = apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, hill_coeff_apoptosis, K_half_apoptosis, S_anti);
    // std::cout << "Apoptosis rate is: " << apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, hill_coeff_apoptosis, K_half_apoptosis, S_anti_real) << std::endl;
    // std::cout << "s_anti_real is: " << S_anti_real << std::endl;

    double growth_rate = phenotype.cycle.data.transition_rate(0,0);
    growth_rate = growth_mapping_logistic(growth_rate_basal, hill_coeff_growth, K_half_growth, S_pro);
    // std::cout << "s_pro_real is: " << S_pro_real << std::endl;

    // growth_rate *= 1 - pressure_effect;
    if (growth_rate < 0.0) growth_rate = 0.0; // Sanity check

    // std::cout << "Pressure effect is: " << pressure_effect << std::endl;
    std::cout << "Growth rate: " << growth_rate << std::endl;
    
    return;
}


void update_monitor_variables(Cell* pCell )
{
	static int mek_node_ix = pCell->custom_data.find_variable_index("mek_node");
	static int akt_node_ix = pCell->custom_data.find_variable_index("akt_node");
	static int pi3k_node_ix = pCell->custom_data.find_variable_index("pi3k_node");
	static int tak1_node_ix = pCell->custom_data.find_variable_index("tak1_node");

    static int anti_mek_node_ix = pCell->custom_data.find_variable_index("anti_mek_node");
    static int anti_akt_node_ix = pCell->custom_data.find_variable_index("anti_akt_node");
    static int anti_pi3k_node_ix = pCell->custom_data.find_variable_index("anti_pi3k_node");
    static int anti_tak1_node_ix = pCell->custom_data.find_variable_index("anti_tak1_node");

    static int antisurvival_b1_ix = pCell->custom_data.find_variable_index("antisurvival_b1_node");
    static int antisurvival_b2_ix = pCell->custom_data.find_variable_index("antisurvival_b2_mode");
    static int antisurvival_b3_ix = pCell->custom_data.find_variable_index("antisurvival_b3_node");

    static int prosurvival_b1_ix = pCell->custom_data.find_variable_index("prosurvival_b1_node");
    static int prosurvival_b2_ix = pCell->custom_data.find_variable_index("prosurvival_b2_node");
    static int prosurvival_b3_ix = pCell->custom_data.find_variable_index("prosurvival_b3_node");

	pCell->custom_data[mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "MEK" );
    pCell->custom_data[akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "AKT" );
    pCell->custom_data[pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "PI3K" );
    pCell->custom_data[tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "TAK1" );

    pCell->custom_data[anti_mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_MEK" );
    pCell->custom_data[anti_akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_AKT" );
    pCell->custom_data[anti_pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_PI3K" );
    pCell->custom_data[anti_tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_TAK1" );

    pCell->custom_data[antisurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    pCell->custom_data[antisurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    pCell->custom_data[antisurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    pCell->custom_data[prosurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    pCell->custom_data[prosurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    pCell->custom_data[prosurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    return;
}


void ags_bm_interface_main(Cell* pCell, Phenotype& phenotype, double dt)
{
    // std::cout << "Updating Boolean Model" << std::endl;

    // pCell->phenotype.intracellular->print_current_nodes();

    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

    if ( pCell->phenotype.intracellular->need_update() )
    {

        // std::cout << "This is the time: " << PhysiCell_globals.current_time << std::endl;

        // if (PhysiCell_globals.current_time >= parameters.ints("time_add_dA") - 10.0 && (PhysiCell_globals.current_time <= parameters.ints("time_add_dA") + 10.0)) 
        // {
        //     std::cout << "Restarting BM " << std::endl;
        //     pCell->phenotype.intracellular->start();
        // }

        // First we update the Boolean Model inputs
        update_boolean_model_inputs(pCell, phenotype, dt );
    
        // Run maboss to update the boolean state of the cell
        pCell->phenotype.intracellular->update();

        // update the cell fate based on the boolean outputs
        update_cell_from_boolean_model(pCell, phenotype, dt);

        // Get track of some boolean node values for debugging
        update_monitor_variables(pCell);
    }

    return;
}


bool boolean_node_deactivation_prob(double drug_density, double scaling, double GI50 ){

    // Logistic function

    double node_probability = 1 / (1 + std::exp(- log10(drug_density / GI50 ) * scaling));
    
    std::default_random_engine generator;
    std::bernoulli_distribution distro(node_probability);

    if (distro(generator)){ node_probability = 1;} else { node_probability = 0;}
    // if (node_probability <= 0.5){ node_probability = 0;} else { node_probability = 1;}

    return node_probability;

}

double growth_mapping_logistic(double doubling_time, double hill_coeff, double K_half, double S_value){

    // double growth_logistic_function = doubling_time / (1 + std::exp(- log10(readout_value) * scaling));

    double growth_logistic_function;
    growth_logistic_function = (doubling_time * std::pow(S_value, hill_coeff ) ) / (K_half + std::pow(S_value, hill_coeff) ) ;



    // if (readout_value == 1){ 

    // } else if (readout_value == 2){
    //     S_value = w1*1 + w2*2;
    //     growth_logistic_function = (doubling_time * std::pow(hill_coeff, S_value) ) / (K_half + std::pow(hill_coeff, S_value) ) ;

    // } else if (readout_value == 3){
    //     S_value = w1*1 + w2*2 + w3*3;
    //     growth_logistic_function = (doubling_time * std::pow(hill_coeff, S_value) ) / (K_half + std::pow(hill_coeff, S_value) ) ;
    // }

    return growth_logistic_function;

}


double apoptosis_mapping_logistic(double basal_apoptosis_rate, double maximum_apoptosis_rate, double hill_coeff, double K_half, double S_value){

    // double apoptosis_rate_logistic_function = (maximum_apoptosis_rate / (1 + std::exp(- log10(readout_value) * scaling)));
    // if (readout_value == 3) return 0.1;

    // return apoptosis_rate_logistic_function;

    double apoptosis_mapping_function;
    apoptosis_mapping_function = (maximum_apoptosis_rate * std::pow(S_value, hill_coeff) ) / (K_half + std::pow(S_value, hill_coeff) ) ;

    return apoptosis_mapping_function +  basal_apoptosis_rate;

    // if (readout_value == 1){ 
    //     S_value = w1*readout_value;

    // } else if (readout_value == 2){
    //     S_value = w1*1 + w2*2;
    //     apoptosis_mapping_function = (maximum_apoptosis_rate * std::pow(hill_coeff, S_value) ) / (K_half + std::pow(hill_coeff, S_value) ) ;

    // } else if (readout_value == 3){
    //     S_value = w1*1 + w2*2 + w3*3;
    //     apoptosis_mapping_function = (maximum_apoptosis_rate * std::pow(hill_coeff, S_value) ) / (K_half + std::pow(hill_coeff, S_value) ) ;
    // }
    
}

double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half){

    // double pressure_exponential_function = std::pow(6e-03, pressure);
    double pressure_exponential_function =  std::pow(pressure, hill_coeff) / (pressure_half + std::pow(pressure, hill_coeff));
    // if (pressure_exponential_function > 1) pressure_exponential_function = 1.0;
    return pressure_exponential_function;
}



// if(casp37_b1){ std::cout << "node caspase37 b1 is ON " << std::endl; }
// if(casp37_b2){ std::cout << "node caspase37 b2 is ON " << std::endl; }
// if(FOXO){ std::cout << "node FOXO is ON " << std::endl; }

// if(cMYC){ std::cout << "node cMYC is ON " << std::endl; }
// if(CCND_b1){ std::cout << "node CCND_b1 is ON " << std::endl; }
// if(CCND_b2){ std::cout << "node CCND_b2 is ON " << std::endl; }

// bool casp37_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b1" );
// bool casp37_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b2" );
// bool FOXO = pCell->phenotype.intracellular->get_boolean_variable_value( "FOXO" );

// bool cMYC = pCell->phenotype.intracellular->get_boolean_variable_value( "cMYC" );
// bool CCND_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND_b1" );
// bool CCND_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND_b2" );


// Drug A
// if ( I_dA >= GI50_dA - 1e-05 && PhysiCell_globals.current_time >= parameters.ints("time_add_dA")  && parameters.strings("node_affected_dA") != "none" ) 
// {
//     // std::cout << distro(generator) << std::endl;

//     if (distro(generator)){
//         pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dA"), 1);
//     } else {
//         pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dA"), 0);
//     }
    
// } 


// // else { 
// //     // pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dA"), 0);
// // }

// // Drug B interface
// if ( I_dA >= GI50_dB - 1e-05 && PhysiCell_globals.current_time >= parameters.ints("time_add_dB")  && parameters.strings("node_affected_dB") != "none" ) // || pCell->custom_data[I_dB] > pCell->custom_data[activation_threshold_ix]
// {
//     if (distro(generator)){
//         pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dB"), 1);
//     } else {
//         pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dB"), 0);
//     }
// } 
// else { 
//     pCell->phenotype.intracellular->set_boolean_variable_value(parameters.strings("node_affected_dB"), 0);
// }


    // if (antisurvival_b1){
    //     // std::cout << "node anti b1 ON" << std::endl;
    //     // pCell->phenotype.death.rates[apoptosis_model_index] = 0.005;
    //     // // growth_rate = 0.023;
    //     pCell->phenotype.death.rates[apoptosis_model_index] = apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, 1, hill_coeff_apoptosis, K_half_apoptosis, S_anti);

    // } 
    
    // if ( antisurvival_b2){
    //     // std::cout << "node anti b2 ON" << std::endl;
    //     // pCell->phenotype.death.rates[apoptosis_model_index] = 0.009;
    //     // growth_rate = 0.013;
    //     pCell->phenotype.death.rates[apoptosis_model_index] = apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, 2, hill_coeff_apoptosis, K_half_apoptosis, S_anti);


    // }  
    
    // if (antisurvival_b3){
    //     // std::cout << "node anti b3 ON" << std::endl;
    //     // pCell->phenotype.death.rates[apoptosis_model_index] = 0.1;
    //     // growth_rate = 0.001;
    //     pCell->phenotype.death.rates[apoptosis_model_index] = apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, 3, hill_coeff_apoptosis, K_half_apoptosis);

    // }







    // if (prosurvival_b1){
    //     // std::cout << "node pro b1 ON" << std::endl;
    //     // growth_rate = 0.09;
    //     growth_rate = growth_mapping_logistic(growth_rate_basal, hill_coeff_growth, K_half_growth, S_pro, 1);
    // }
    
    // if ( prosurvival_b2){
    //     // std::cout << "node pro b2 ON" << std::endl;
    //     // growth_rate = 0.05;
    //     growth_rate = growth_mapping_logistic(growth_rate_basal, hill_coeff_growth, K_half_growth, 2);
    // } 
    
    // if (prosurvival_b3){
    //     // std::cout << "node pro b3 ON" << std::endl;
    //     // growth_rate = 0.046;
    //     growth_rate = growth_mapping_logistic(growth_rate_basal, hill_coeff_growth, K_half_growth, 3);
    // }