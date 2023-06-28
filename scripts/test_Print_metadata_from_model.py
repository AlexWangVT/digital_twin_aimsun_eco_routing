section_test = model.getCatalog().find(3306)
print(section_test.getType().getName())
print(section_test.getName())

section_type = section_test.getType()
cols=section_type.getColumns(GKType.eSearchOnlyThisType)
for col in cols:
    if col.getColumnType() != GKColumn._GKTimeSerie and section_test.getDataValue(col)[0] != None:
        print(col.getName()+": "+str(section_test.getDataValue(col)[0]))
