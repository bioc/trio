.onAttach <- function(libname, pkgname){
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui"){
		winMenuAdd("Vignettes")
		winMenuAdd("Vignettes/trio")
		winMenuAddItem("Vignettes/trio", "Vignette", "vignette('trio', package='trio')")
		winMenuAddItem("Vignettes/trio", "R Code", "edit(vignette('trio', package='trio'))")
	}
}


