# svn $Id: Module.mk 975 2009-05-05 22:51:13Z kate $
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2009 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

local_sub  := ROMS/Fish

local_lib  := libFish.a
local_src  := $(wildcard $(local_sub)/*.F)
local_c_lib  := libFish_c.a
local_c_src  := $(wildcard $(local_sub)/*.cc)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(call make-c-library,$(local_c_lib),$(local_c_src)))

$(eval $(compile-rules))

#$(eval $(c-compile-rules))
Build/c_tree.o: $(local_sub)/c_tree.cc
	cd $(SCRATCH_DIR); $(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $(notdir $@) ../$<

