// much of this adapted with little modification from SimOptions.h in Fang Da's implementation of "Surface-Only Liquids"

#include "simoptions.h"

std::map<std::string, SimOptions::SimOption> SimOptions::sim_options;

void SimOptions::addStringOption(const std::string & key, const std::string & default_value){
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = STRING;
    o.str_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addIntegerOption(const std::string & key, int default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = INTEGER;
    o.int_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addDoubleOption(const std::string & key, double default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = DOUBLE;
    o.double_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addBooleanOption(const std::string & key, bool default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = BOOLEAN;
    o.bool_value = default_value;
    
    sim_options[key] = o;
}

const std::string & SimOptions::strValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == STRING);          // verify this option has the correct type
    return sim_options[key].str_value;
}

int SimOptions::intValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == INTEGER);         // verify this option has the correct type
    return sim_options[key].int_value;
}

double SimOptions::doubleValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == DOUBLE);          // verify this option has the correct type
    return sim_options[key].double_value;
}

bool SimOptions::boolValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == BOOLEAN);         // verify this option has the correct type
    return sim_options[key].bool_value;
}