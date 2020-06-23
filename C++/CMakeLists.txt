Project(ThighMuscleandFatQuantification)

cmake_minimum_required(VERSION 2.6)

SET (LIBRARY_OUTPUT_PATH ${ThighMuscleandFatQuantification_BINARY_DIR}/bin/ 
        CACHE PATH 
        "Single output directory for building all libraries."
        FORCE)
SET (EXECUTABLE_OUTPUT_PATH ${ThighMuscleandFatQuantification_BINARY_DIR}/bin 
        CACHE PATH 
        "Single output directory for building all executables."
        FORCE)
SET (BUILD_SHARED_LIBS FALSE CACHE BOOL "build shared libs" )

#Find ITK
FIND_PACKAGE( ITK REQUIRED)
IF ( ITK_FOUND )
  INCLUDE( ${ITK_USE_FILE} )
ELSE (ITK_FOUND)
  MESSAGE(ERROR "Cannot build without a ITK build tree.")
ENDIF ( ITK_FOUND )

SET( SEGMENTATION_ITK_LIBS  
     ITKAlgorithms
     ITKBasicFilters
     ITKStatistics
     ITKNumerics
     ITKCommon
     ITKIO
)

ADD_EXECUTABLE(ThighMuscleandFatQuantification ThighMuscleandFatQuantification.cxx)

TARGET_LINK_LIBRARIES(ThighMuscleandFatQuantification ${SEGMENTATION_ITK_LIBS} )



