#include "maboss_network.h"

/* Default constructor */
void MaBoSSNetwork::init_maboss( std::string networkFile, std::string configFile)
{
	if (this->network != NULL) {
		delete this->network;
	}
	
	if (this->config != NULL) {
		delete this->config;
	}
	
	if (this->engine != NULL) {
		delete this->engine;
	}
	// Initialize MaBoSS Objects for a model
	this->network = new Network();
	this->network->parse(networkFile.c_str());

	this->config = new RunConfig();
	this->config->parse(this->network, configFile.c_str());

	IStateGroup::checkAndComplete(this->network);

	engine = new StochasticSimulationEngine(this->network, this->config, PhysiCell::UniformInt());

	this->update_time_step = this->config->getMaxTime();
	
	// Building map of nodes for fast later access 
	for (auto node : this->network->getNodes()) {
		this->nodesByName[node->getLabel()] = node;
	}
	
	// Building map of parameters for fast later access
	for (auto parameter : this->network->getSymbolTable()->getSymbolsNames()) {
		if (parameter[0] == '$')
			this->parametersByName[parameter] = this->network->getSymbolTable()->getSymbol(parameter);
	}
	
	for (auto node : network->getNodes())
      if (!node->isInternal()) 
        output_mask.setNodeState(node, true);

}

void MaBoSSNetwork::mutate(std::map<std::string, double> mutations) 
{
	for (auto mutation : mutations) {
		nodesByName[mutation.first]->mutate(mutation.second);
	}
}

void MaBoSSNetwork::set_parameters(std::map<std::string, double> parameters) 
{	
	for (auto parameter: parameters) {
		set_parameter_value(parameter.first, parameter.second);
	}
}

double MaBoSSNetwork::get_parameter_value(std::string name) 
{
	return network->getSymbolTable()->getSymbolValue(parametersByName[name]);
}


void MaBoSSNetwork::set_parameter_value(std::string name, double value) 
{
	network->getSymbolTable()->setSymbolValue(parametersByName[name], value);
	network->getSymbolTable()->unsetSymbolExpressions();
}

/* Reset a vector of bools to the init state of the network */
void MaBoSSNetwork::restart_node_values()
{
	// NetworkState network_state;
	this->network->initStates(state, engine->random_generator);
	
	for (auto initial_value : initial_values) {
		state.setNodeState(nodesByName[initial_value.first], PhysiCell::UniformRandom() < initial_value.second);
	}
	
	this->set_time_to_update();
}

/* Run a MaBoSS simulation with the input values*/
void MaBoSSNetwork::run_simulation()
{	
	engine->setMaxTime(time_to_update/scaling);
	state = engine->run(state, NULL);
	this->set_time_to_update();

}

bool MaBoSSNetwork::has_node( std::string name ) {
	return nodesByName.find(name) != nodesByName.end();
}

void MaBoSSNetwork::set_node_value(std::string name, bool value) {
	state.setNodeState(nodesByName[name], value);
}

bool MaBoSSNetwork::get_node_value(std::string name) {
	return state.getNodeState(nodesByName[name]);
}

std::string MaBoSSNetwork::get_state() {
	return NetworkState(state.getState() & output_mask.getState()).getName(network);
}

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_nodes()
{
	const std::string filename = "/home/oth/BHS/TFM/TFM_modelling/PhysiBoSS/output/BM_readouts.txt";
	std::ofstream readouts_file(const std::string filename);

	int i = 0;
	std::vector<Node*> nodes = this->network->getNodes();
	for ( auto node: nodes )
	{
		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		// readouts_file << std::to_string(node->getLabel()) << "=" << std::to_string(state.getNodeState(node)) << "; " << std::endl;
		i++;
	}
	std::cout << std::endl;

	// readouts_file.close();

}


// added

/* Print current state of all the nodes of the network */
void MaBoSSNetwork::print_survival_nodes()
{

	// const std::string filename = "/home/oth/BHS/TFM/TFM_modelling/PhysiBoSS/output/BM_readouts.txt";
	
	// std::ofstream readouts_file(const std::string filename);

	// readouts_file << "Node = Status;" << std::endl;

	// int i = 0;
	// std::vector<Node*> nodes = this->network->getNodes();

	// const std::string pro_b1 = "Prosurvival_b1";
	// const std::string pro_b2 = "Prosurvival_b2";
	// const std::string pro_b3 = "Prosurvival_b3";

	// const std::string anti_b1 = "Antisurvival_b1";
	// const std::string anti_b2 = "Antisurvival_b2";
	// const std::string anti_b3 = "Antisurvival_b3";

	// for ( auto node: nodes )
	// {

	// 	// PRO nodes

	// 	if(node->getLabel() == pro_b1)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		// readouts_file << static_cast<MaBoSSNetwork*>(this->network->node)->getLabel() << std::endl;
	// 		// readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		
	// 	if(node->getLabel() == pro_b2)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		// readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	
	// 	if(node->getLabel() == pro_b3)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		// readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		
	// 	// ANTI nodes

	// 	if(node->getLabel() == anti_b1)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		
	// 	if(node->getLabel() == anti_b2)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		
	// 	if(node->getLabel() == anti_b3)
	// 		std::cout << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
	// 		readouts_file << node->getLabel() << "=" << state.getNodeState(node) << "; " << std::endl;
		
		
	// 	i++;
	// }

	
	// std::cout << std::endl;

	// readouts_file.close();
}