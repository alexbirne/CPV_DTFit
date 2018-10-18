# Library
root_generate_dictionaries_local(decratecoeff_bd_dict DecRateCoeff_Bd.h)
add_library(decratecoeff_bd SHARED DecRateCoeff_Bd.h DecRateCoeff_Bd.cxx DecRateCoeff_Bd_dict.h DecRateCoeff_Bd_dict.cc)
target_link_libraries(decratecoeff_bd ${ALL_LIBRARIES})

# Executables

# Aliases for module libs

# Linker information
add_dependencies(decratecoeff_bd decratecoeff_bd_dict)
