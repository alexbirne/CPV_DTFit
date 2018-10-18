# Library
add_library(CPTimeFit CPTimeFit.h
                      CPTimeFit.cpp
                      BuildPDF/build_PDF.cpp
                      BuildPDF/build_acceptance.cpp
                      BuildPDF/build_resolution.cpp
                      BuildPDF/build_mistaghistogram.cpp
                      BuildPDF/calculateSweightCorrection.cpp
                      Fit/Fit.cpp
                      Fit/buildConstraint.cpp
                      Fit/build4DConstraint.cpp
                      SetData/setData.cpp
                      SetData/doBootstrapping.cpp
                      SetData/GenerateToy.cpp
                      SetData/TrueToyTagger.cpp
                      Plot/Plot.cpp)

# Executables

# Aliases for module libs

# Linker information
target_link_libraries(CPTimeFit decratecoeff_bd Espresso ${ALL_LIBRARIES})
# target_link_libraries(CPTimeFit decratecoeff_bd ${ALL_LIBRARIES})

add_subdirectory(DecRateCoeff)
