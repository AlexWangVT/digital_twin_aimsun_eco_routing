sectionType = model.getType( "GKSection" )
for types in model.getCatalog().getUsedSubTypesFromType( sectionType ):
    for s in types.values():
        print(("Section ID: %i Name: %s") % (s.getId(), s.getName()))
        print(s.__class__.__name__)
    #print(s.getId(), s.getName())