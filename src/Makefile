##
## redirect 'make' requests to base directory
##
#
# This file is part of Velour and is distributed under the University of
# Illinois Open Source License. See LICENSE.txt for details.
#

BASEDIR = ..

.PHONY: default

default:
	@$(MAKE) -C $(BASEDIR) PWD="$(PWD)/.."

%:
	@$(MAKE) -C $(BASEDIR) $@ PWD="$(PWD)/.."

