classifyStructures <-
function(structures,classes,classValues,classToExclude,cancerIdentifier,classOther){
	message("classify structures")
	allStructures=unique(structures)
#take only real structures
	allStructures=allStructures[allStructures != 0]
	for (structure in allStructures){
		cellsInStructure=classes[structures==structure]
		#cell types in structure
		cellTypes=unique(cellsInStructure)
		#i need to fix that
		for (cellT in cellTypes){
			if(cellT == "a"){ #this is hard coded, has to be changed!!!
				cellsInStructure[cellsInStructure==cellT]=cancerIdentifier
			}
		}
		classes[structures==structure]=cellsInStructure
	}
	classes
}
