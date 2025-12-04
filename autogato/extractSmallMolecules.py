
import sys, os
import libsbml

path = "./XML-Files/"
speciesList = os.listdir(path)
speciesXMLDict = {}
metabolites = set()
print(speciesList)
for species in speciesList:    
    print(species)
    if species == ".DS_Store":
        continue
    if species == "HelicobacterPylori":
        continue
    if species == "EColiCore":
        continue
    if species == "arabidopsisThaliana":
        continue
    speciesPath = path+species+"/"
    xmlFileList =  os.listdir(speciesPath)
    for xmlFile in xmlFileList:
    # 2. Read Model
        xmlFilePath = speciesPath + xmlFile
        reader = libsbml.SBMLReader()
        document = reader.readSBML(xmlFilePath)
        model = document.getModel()
        for metaboliteReference in model.getListOfSpecies():
            metaboliteID = metaboliteReference.getId()
            metabliteName = metaboliteReference.getName()
            metabolites.add(str(metaboliteID + "\t \t \t \t \t " + metabliteName))
metabolites = sorted(list(metabolites))
with open ("allMetabolitesThermococcus.txt", "w") as file:
    for m in metabolites:
        file.write(m)
        file.write("\n")
    file.close()

