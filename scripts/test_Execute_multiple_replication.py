#for s in GKSection.values():
#	print("ID: %i Name: %s" % (s.getId(), s.getName()))
#print(target.getId())
#for object in selection:
#	print("ID: %i Name: %s" % (object.getId(), object.getName()))
replications = [model.getCatalog().find(393921)]
replications.append(model.getCatalog().find(394033))

for re in replications:
	print(re.getName(), re.getId())
GKSystem.getSystem().executeAction( "execute", replications[0], replications, "" )
#GKSystem.getSystem().executeAction( "play", replication, [], "" )
#GKSystem.getSystem().executeAction( "retrieve", replication, [], "" )
