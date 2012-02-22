.onLoad <-function(lib,pkg){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),
  "Version")
  packageStartupMessage("blockcluster version ", as.character(ver), " loaded\n\nCondition of use\n----------------\nCopyright (C)  <MODAL team @INRIA,Lille & U.M.R. C.N.R.S. 6599 Heudiasyc, UTC>\nblockcluster is publicly available under the GPL license (see www.gnu.org/copyleft/gpl.html)\nYou can redistribute it and/or modify it under the terms of the GPL-3 license.\nPlease be informed that there may still be bugs and errors. Use it at your own risk.\nPlease post questions and bugs at: <https://gforge.inria.fr/forum/forum.php?forum_id=11190&group_id=3679>\n")
}
