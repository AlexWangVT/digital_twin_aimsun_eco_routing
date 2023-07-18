print("##### This is the start of script_Modify_CAV_percentage_all_scenario #####")

random_seed_list = [59940404]
demand_percentage = 100
prediction_horizon = 0
vehicle_types = ["ICE", "BEV", "PHEV", "HFCV"]
vehicle_ratios = [94, 3, 2, 1]
CAV_penetration_list = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
eco_routing_with_travel_time = False

if target == None:
	scenario = model.getCatalog().find(394099)
else:
	scenario = model.getCatalog().find(target.getId())

scenario_name = scenario.getName()

if scenario_name.find("_P") == -1:
	print("No valid scenario found! Need to have a name end with _P{x}, where {x} is the prediction horizon. Reset to default scenario.")
	scenario = model.getCatalog().find(394099)
else:
	prediction_horizon = int(scenario_name[scenario_name.find("_P")+2:])
	print("The predict horizon for current scenario {} is {}.".format(scenario.getName(), prediction_horizon))

if scenario.getName().find("Travel_Time") != -1:
	eco_routing_with_travel_time = True
	print("Will do eco routing with travel time for scenario {}".format(scenario.getName()))

print("Will rename and change the demand percentage for scenario {}'s all experiments:".format(scenario.getName()))

# change demand percentage for the scenario
scenario.setValueForVariable("$demand_percentage", str(demand_percentage))


# get attributes in the experiments
attribute_demand = model.getColumn("GKExperiment::demand_percentage")
attribute_prediction = model.getColumn("GKExperiment::prediction_horizon")
attribute_cav_penetration = model.getColumn("GKExperiment::cav_penetration")
attribute_eco_routing_with_travel_time = model.getColumn("GKExperiment::eco_routing_with_travel_time")

for idx, exp in enumerate(scenario.getExperiments()):
	CAV_penetration = CAV_penetration_list[idx]
	exp.setName("Experiment_P{}_CAV{}".format(prediction_horizon, CAV_penetration))
	
	# change ratio variable for each vehicle type in the experiments
	cav_ratio = CAV_penetration / 100.0
	noncav_ratio = (100 - CAV_penetration) / 100.0
	for idx in range(len(vehicle_types)):
		exp.setValueForVariable("$"+vehicle_types[idx]+"_percentage", str(cav_ratio * vehicle_ratios[idx]))
		exp.setValueForVariable("$"+vehicle_types[idx]+"_NONCAV_percentage", str(noncav_ratio * vehicle_ratios[idx]))
	# set attributes in the experiments for correct logging in the API
	exp.setDataValueDouble(attribute_demand, demand_percentage)
	exp.setDataValueDouble(attribute_prediction, prediction_horizon)
	exp.setDataValueDouble(attribute_cav_penetration, CAV_penetration)
	exp.setDataValue(attribute_eco_routing_with_travel_time, eco_routing_with_travel_time)
	
	for re_idx, repli in enumerate(exp.getReplications()):
		repli.setName("Replication_P{}_CAV{}_{}".format(prediction_horizon, CAV_penetration, re_idx))
		repli.setRandomSeed(random_seed_list[re_idx % len(random_seed_list)])


print("##### End of script_Modify_CAV_percentage_all_scenario #####")

model.getCommander().addCommand( None )