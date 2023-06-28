updated_dynamic_function = model.getCatalog().find(394219)

turningType = model.getType( "GKTurning" )
for types in model.getCatalog().getUsedSubTypesFromType( turningType ):
    for turning in types.values():
        print("original dynamic function is : " + str(turning.getDynamicFunction().getId()) )
        turning.setDynamicFunction(updated_dynamic_function)
        print("after modinication, dynamic function is : " + str(turning.getDynamicFunction().getId()) )
    print("Total number of turnings is: " + str(len(types)))

#section0 = model.getCatalog().find(3306)
#node0 = model.getCatalog().find(234228)
#print(len(section0.getDestTurnings()))
#print(section0.getDestTurnings()[0].getDynamicFunction().getId())
model.getCommander().addCommand( None )