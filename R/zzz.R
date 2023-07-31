.onLoad <- function(CAT, pkgname) {
  logo_info <- c("
                      ______             ___           ____________
                     /      |           /   \\\         |____    ____|
                    |  ,----'          /  ^  \\\             |  |
                    |  |              /  /_\\\  \\\            |  |
                    |  `----.        /  _____  \\\           |  |
                     \\______|       /__/     \\__\\          |__|
")

  message(logo_info)
  message("Author: Tao Zhu")
  message(sprintf("CAT : Version %s",packageVersion("CAT")))
  #message("For more information see our website : www.ArchRProject.com")
  message("If you encounter a bug please report : https://github.com/tzhu-bio/CAT/issues")
}
