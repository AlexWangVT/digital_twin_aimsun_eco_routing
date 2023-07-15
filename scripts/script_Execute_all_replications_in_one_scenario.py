import random
import time

scenario = model.getCatalog().find(394325)
N = 11
CAV_range = [0, 100, 10]
random.seed(time.time())
seed_list = []
all_replications = []
all_averages = []
repli_map = {}
start_time = time.time()

def get_padding(per):
	return str(per).zfill(3)

def get_shortname(fullname):
	for penetration in range(CAV_range[0], CAV_range[1]+CAV_range[2], CAV_range[2]):
		if get_padding(penetration) in fullname:
			return get_padding(penetration)

# init replication map and seed list
for penetration in range(CAV_range[0], CAV_range[1]+CAV_range[2], CAV_range[2]):
	repli_map[get_padding(penetration)] = []
	seed_list.append(random.randrange(100000000))

# collect replications into corresponding dictionaries
for exp in scenario.getExperiments():
	for repli in exp.getReplications():
		if not repli.isA("GKExperimentResult"):
			all_replications.append(repli)
			repli_map[get_shortname(exp.getName())].append(repli)
		else:
			all_averages.append(repli)

# use the same random seeds in all CAV penetrations
#for shortname in repli_map:
#	for idx, repli in enumerate(repli_map[shortname]):
#		repli_map[shortname][idx].setRandomSeed(seed_list[idx])

print("Will execute all replications in these experiments:")
for exp in scenario.getExperiments():
	print(exp.getName())
print("Total number of replications are {}.".format(len(all_replications)))

# this will execute all replications, no easy way to stop 
#GKSystem.getSystem().executeAction( "execute", all_replications[0], all_replications, "" )

 
#GKSystem.getSystem().executeAction( "execute", repli_map["010"][0], repli_map["010"], "" )
#for penetration in range(10, 110, 10):
#	GKSystem.getSystem().executeAction( "execute", repli_map[get_padding(penetration)][0], repli_map[get_padding(penetration)], "")
#	GKSystem.getSystem().executeAction( "unload", repli_map[get_padding(penetration)][0], repli_map[get_padding(penetration)], "")

#GKSystem.getSystem().executeAction( "play", replication, [], "" )
#GKSystem.getSystem().executeAction( "retrieve", all_replications[0], all_replications, "" )

 
#GKSystem.getSystem().executeAction( "execute", all_averages[0], all_averages, "" )
#GKSystem.getSystem().executeAction( "unload", all_averages[0], all_averages, "" )
#GKSystem.getSystem().executeAction( "retrieve", all_averages[0], all_averages, "" )

end_time = time.time()
print("Total time elapsed: {} seconds".format(end_time - start_time))