# (C) Copyright 1989- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# 
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction

ecbuild_info("[surf]")

ecbuild_list_add_pattern(LIST surf.${PREC}_src GLOB
  surf/module/* 
  surf/external/*
)


ecbuild_add_library(TARGET surf.${PREC}
  DEFINITIONS ${IFS_DEFINITIONS}
  #SOURCES_GLOB surf/module/* surf/external/*
  SOURCES ${surf.${PREC}_src} 
  PUBLIC_INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/surf/interface>
  PRIVATE_INCLUDES surf/function 
  PUBLIC_LIBS ${IFSAUX_LIBRARIES})
