.onAttach <- function(libname, pkgname){
	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui"){
		winMenuAdd("Vignettes")
		winMenuAdd("Vignettes/trio9")
		winMenuAddItem("Vignettes/trio9", "Vignette", "vignette('trio', package='trio9')")
		winMenuAddItem("Vignettes/trio9", "R Code", "edit(vignette('trio', package='trio9'))")
	}
}
