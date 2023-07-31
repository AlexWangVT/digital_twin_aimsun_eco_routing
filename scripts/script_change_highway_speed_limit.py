target_speed = 120
cnt = 0
sectionType = model.getType( "GKSection" )
for types in model.getCatalog().getUsedSubTypesFromType( sectionType ):
	for s in types.values():
		secid=s.getId();
		road_type = s.getRoadType().getName()
		if road_type == "Motorway":
			cnt+=1
			s.setSpeed(target_speed)
print("Total Motorway: {}".format(cnt))
print("Speed limit set to: {}".format(target_speed))
