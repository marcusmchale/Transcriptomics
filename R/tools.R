# Title     : Tools to support general procedures
# Objective : A place to gather general tools that don't fit into the main project functions
# Created by: marcus
# Created on: 10/02/2020

# make the pipe operator available without loading dplyr
# working out how to avoiding library so easier to package code later
`%>%` <- magrittr::`%>%`

Pathing <- R6::R6Class("Pathing", list(
  base_path = NULL,
  current_path = NULL,
  initialize = function(analysis_root='analysis', extra_path=as.Date(Sys.Date(), format = '%Y%m%d')) {
	self$base_path <- file.path(analysis_root, extra_path)
	dir.create(self$base_path, recursive = TRUE, showWarnings = FALSE)
	self$current_path <- self$base_path
	invisible(self)
  },
  append_to_path = function(folder) {
	self$current_path <- file.path(self$current_path, folder)
	if (!dir.exists(self$current_path)) {
	  dir.create(self$current_path, recursive=TRUE, showWarnings = FALSE)
	}
	invisible(self) # for a side-effect methods like this try to return self object to allow chaining
  },
  return_to_base = function() {
	self$current_path = self$base_path
	invisible(self)
  },
  # function to recursively split path to get list of directories starting with base and ending with the root
  split_path = function(path) {
	if (dirname(path) %in% c("..", path)) return(basename(path))
	return(c(basename(path), split_path(dirname(path))))
  }
))
# R6 class usage e.g.
# pathing <- Pathing$new()
# pathing$split_path('asdf/asdf/asdf')
