/*===========================================================================

  Program:   Quantification of muscle and intramusclurar fat from MRI.
  Module:    $RCSfile: ThighMuscleQuantification.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/16 11:12:32 $
  Version:   $Revision: 0.1 $

  3T MRI Facility National Institute on Aging/National Institutes of Health.

=============================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageFileReader.h" 
#include "itkImageFileWriter.h" 

#include "itkOrientedImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"

#include "itkSubtractImageFilter.h"
//#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include <itkLaplacianSharpeningImageFilter.h>
#include "itkMinimumMaximumImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

// Doxygen comments.
/*! \mainpage Muscle Quantification Index Page
 *
 */

const unsigned int Dimension=3;
typedef float FloatPixelType;
typedef unsigned char MaskPixelType;
typedef unsigned int LabelPixelType;
typedef itk::Image<FloatPixelType, Dimension> FloatImageType;
typedef itk::Image<MaskPixelType, Dimension> MaskImageType;
typedef itk::Image<LabelPixelType, Dimension> LabelImageType;


// Typical dilation followed by erosion.

MaskImageType::Pointer BinaryMorphologicalClosing( MaskImageType::Pointer InputImage,
						   float structureElementRadius )
{

  // First, generate a structuring element for binary morphology.
  typedef itk::BinaryBallStructuringElement<MaskPixelType, Dimension> 
    StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius( structureElementRadius );
  structuringElement.CreateStructuringElement();
  std::cout << "Using a circle str. element with radius "
	    << structureElementRadius
	    << std::endl;

  
  // Morphological binary dilation.
  typedef itk::BinaryDilateImageFilter<MaskImageType,
    MaskImageType,
    StructuringElementType>
    BinaryDilateFilterType;
  BinaryDilateFilterType::Pointer binaryDilateFilter =
    BinaryDilateFilterType::New();
  binaryDilateFilter->SetKernel( structuringElement );
  // binaryDilateFilter->SetInput( binaryErodeFilter->GetOutput() );
  binaryDilateFilter->SetInput( InputImage );
  binaryDilateFilter->SetDilateValue( 1 );
  binaryDilateFilter->Update();
  std::cout << "Binary Dilation, done" << std::endl;


  // Morphological binary erosion filter to disconnect regions which 
  // are bridged by partial voluming effect.
  typedef itk::BinaryErodeImageFilter<MaskImageType,
    MaskImageType,
    StructuringElementType>
    BinaryErodeFilterType;
  BinaryErodeFilterType::Pointer binaryErodeFilter =
    BinaryErodeFilterType::New();
  binaryErodeFilter->SetKernel( structuringElement );
  // binaryErodeFilter->SetInput( scalarImageKmeansImageFilter->GetOutput() );
  binaryErodeFilter->SetInput( binaryDilateFilter->GetOutput() );
  binaryErodeFilter->SetErodeValue( 1 );
  binaryErodeFilter->Update();
  std::cout << "Binary Erosion, done" << std::endl;

  return binaryErodeFilter->GetOutput();

}


// Typical erosion followed by dilation.

MaskImageType::Pointer BinaryMorphologicalOpening( MaskImageType::Pointer InputImage,
						   float structureElementRadius )
{

  // First, generate a structuring element for binary morphology.
  typedef itk::BinaryBallStructuringElement<MaskPixelType, Dimension> 
    StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius( structureElementRadius );
  structuringElement.CreateStructuringElement();
  std::cout << "Using a circle str. element with radius "
	    << structureElementRadius
	    << std::endl;

  
  // Morphological binary erosion filter to disconnect regions which 
  // are bridged by partial voluming effect.
  typedef itk::BinaryErodeImageFilter<MaskImageType,
    MaskImageType,
    StructuringElementType>
    BinaryErodeFilterType;
  BinaryErodeFilterType::Pointer binaryErodeFilter =
    BinaryErodeFilterType::New();
  binaryErodeFilter->SetKernel( structuringElement );
  binaryErodeFilter->SetInput( InputImage );
  binaryErodeFilter->SetErodeValue( 1 );
  binaryErodeFilter->Update();
  std::cout << "Binary Erosion, done" << std::endl;

  // Morphological binary dilation.
  typedef itk::BinaryDilateImageFilter<MaskImageType,
    MaskImageType,
    StructuringElementType>
    BinaryDilateFilterType;
  BinaryDilateFilterType::Pointer binaryDilateFilter =
    BinaryDilateFilterType::New();
  binaryDilateFilter->SetKernel( structuringElement );
  binaryDilateFilter->SetInput( binaryErodeFilter->GetOutput() );
  binaryDilateFilter->SetDilateValue( 1 );
  binaryDilateFilter->Update();
  std::cout << "Binary Dilation, done" << std::endl;

  return binaryDilateFilter->GetOutput();

}


// Main function.

int main(int argc, char *argv[])
{

  // Inputs will be: Fat suppressed image, non fat suppressed image.

  std::string nonfatsuppressedVolumeFilename;
  std::string fatsuppressedVolumeFilename;
  std::string outputFileName = "MuscleRegion.mhd";

  FloatImageType::Pointer nonfatsuppressedVolume;
  FloatImageType::Pointer fatsuppressedVolume;
  float structureElementRadius = 2.0;
  
  if ( argc < 4 )
    {
    std::cout << "Wrong syntax" << std::endl;
    std::cout << "Usage: " << argv[0] 
	      << " <non-fat suppressed volume> <structuring element radius>" 
	      << std::endl;
    return EXIT_FAILURE;
    } 
  else 
    { 
      nonfatsuppressedVolumeFilename = argv[1]; 
      fatsuppressedVolumeFilename = argv[2];
      structureElementRadius = atof( argv[3] );
    }

 
  // Read-in the volumes.

  typedef itk::ImageFileReader<FloatImageType> FloatReaderType;

  FloatReaderType::Pointer reader = FloatReaderType::New();

  reader->SetFileName( nonfatsuppressedVolumeFilename.c_str() ); 
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  nonfatsuppressedVolume = reader->GetOutput();
  reader = 0;

  std::cout << nonfatsuppressedVolumeFilename 
	    << " Non fat-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties:" << std::endl
	    << nonfatsuppressedVolume->GetBufferedRegion() 
	    << std::endl;
  std::cout << "Spacing:" << std::endl
	    << nonfatsuppressedVolume->GetSpacing() 
	    << std::endl;

  reader = FloatReaderType::New();
  reader->SetFileName( fatsuppressedVolumeFilename.c_str() ); 
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject &ex)
    {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  fatsuppressedVolume = reader->GetOutput();
  reader = 0;


  std::cout << fatsuppressedVolumeFilename 
	    << " Fat-suppressed image read in." 
	    << std::endl;
  std::cout << "Matrix properties:" << std::endl
	    << fatsuppressedVolume->GetBufferedRegion() 
	    << std::endl;
  std::cout << "Spacing:" << std::endl
	    << fatsuppressedVolume->GetSpacing() 
	    << std::endl;

  // First subtract the two volumes.

  typedef itk::SubtractImageFilter<FloatImageType,
    FloatImageType,
    FloatImageType> SubtractImageFilterType;

  SubtractImageFilterType::Pointer subtractImageFilter = 
    SubtractImageFilterType::New();

  subtractImageFilter->SetInput1( nonfatsuppressedVolume );
  subtractImageFilter->SetInput2( fatsuppressedVolume );

  subtractImageFilter->Update();

  std::cout << "Image subtraction, done" << std::endl;
  
  itk::ImageFileWriter<FloatImageType>::Pointer floatWriter = 
    itk::ImageFileWriter<FloatImageType>::New();
  floatWriter->SetInput( subtractImageFilter->GetOutput() ); 
  floatWriter->SetFileName( "difference_image.mhd" );
  floatWriter->Update();
  floatWriter = 0;


  // // Then:
  // // Apply contrast enhancement. (Use itkAdaptiveHistogramEqualizationImageFilter)
  // typedef itk::AdaptiveHistogramEqualizationImageFilter<FloatImageType> 
  //   ImageEnhancementFilterType;
  // ImageEnhancementFilterType::Pointer ImageEnhancementFilter =
  //   ImageEnhancementFilterType::New();

  // ImageEnhancementFilter->SetInput(subtractImageFilter->GetOutput());
  // // Set alpha.
  // float alpha = 0;
  // ImageEnhancementFilter->SetAlpha(alpha);

  // // Set beta.
  // float beta = 0;
  // ImageEnhancementFilter->SetBeta(beta);

  // // Set window size;
  // ImageEnhancementFilterType::ImageSizeType radius;
  // radius.Fill( 3.0 );
  // ImageEnhancementFilter->SetRadius(radius);
  // ImageEnhancementFilter->Update();

  // std::cout << "Adaptive histogram equalization, done" << std::endl;

  // typedef itk::LaplacianSharpeningImageFilter< FloatImageType, FloatImageType >
  //   ImageEnhancementFilterType;
  // ImageEnhancementFilterType::Pointer ImageEnhancementFilter =
  //    ImageEnhancementFilterType::New();
  // ImageEnhancementFilter->SetInput(subtractImageFilter->GetOutput());
  // ImageEnhancementFilter->SetUseImageSpacingOn();
  // ImageEnhancementFilter->Update();
  // std::cout << "Laplacian image sharpening, done" << std ::endl;

  // floatWriter = 
  //   itk::ImageFileWriter<FloatImageType>::New();
  // floatWriter->SetInput( ImageEnhancementFilter->GetOutput() ); 
  // floatWriter->SetFileName( "enhanced_difference_image.mhd" );
  // floatWriter->Update();
  // floatWriter = 0;


  // Clustering (fcm or k-means).
  // Use 3 classes.

  typedef itk::ScalarImageKmeansImageFilter<FloatImageType> 
    ScalarImageKmeansImageFilterType;
  ScalarImageKmeansImageFilterType::Pointer scalarImageKmeansImageFilter = 
    ScalarImageKmeansImageFilterType::New();
  scalarImageKmeansImageFilter->SetInput( subtractImageFilter->GetOutput() );
  //scalarImageKmeansImageFilter->SetInput( nonfatsuppressedVolume );
  scalarImageKmeansImageFilter->SetDebug( true );
  scalarImageKmeansImageFilter->AddClassWithInitialMean( (float) 0.0 );
  scalarImageKmeansImageFilter->AddClassWithInitialMean( (float) 0.0 );
  //scalarImageKmeansImageFilter->AddClassWithInitialMean( (float) 0.0 );
  scalarImageKmeansImageFilter->Update();

  std::cout << "K-means clustering, done" << std::endl;
  std::cout << "Final means are: " << std::endl;
  std::cout << scalarImageKmeansImageFilter->GetFinalMeans() << std::endl;

  itk::ImageFileWriter<MaskImageType>::Pointer maskWriter = 
    itk::ImageFileWriter<MaskImageType>::New();
  maskWriter->SetInput( scalarImageKmeansImageFilter->GetOutput() ); 
  maskWriter->SetFileName( "k_means_output.mhd" );
  maskWriter->Update();
  maskWriter = 0;
  

  // Morphological operations.

  // MaskImageType::Pointer morphologyClosedVolume =
  //   BinaryMorphologicalClosing( scalarImageKmeansImageFilter->GetOutput(),
  // 				structureElementRadius );

  // MaskImageType::Pointer morphologyOpenedVolume = 
  //   BinaryMorphologicalOpening( scalarImageKmeansImageFilter->GetOutput(),
  // 				structureElementRadius );


  // Connected component labeling.

  typedef itk::ConnectedComponentImageFilter<MaskImageType, 
    LabelImageType, 
    MaskImageType> 
    ConnectedComponentLabelFilterType;
  ConnectedComponentLabelFilterType::Pointer labelMaskFilter = 
    ConnectedComponentLabelFilterType::New();
  labelMaskFilter->SetInput( scalarImageKmeansImageFilter->GetOutput() );
  labelMaskFilter->SetMaskImage( scalarImageKmeansImageFilter->GetOutput() );
  labelMaskFilter->SetFullyConnected( true );
  labelMaskFilter->Update();


  // Rank components wrt to size and relabel.

  typedef itk::RelabelComponentImageFilter<LabelImageType, 
    LabelImageType> 
    RelabelFilterType;
  RelabelFilterType::Pointer  sortLabelsImageFilter = 
    RelabelFilterType::New();
  sortLabelsImageFilter->SetInput( labelMaskFilter->GetOutput() );
  sortLabelsImageFilter->Update();

  itk::ImageFileWriter<LabelImageType>::Pointer labelWriter = 
    itk::ImageFileWriter<LabelImageType>::New();
  labelWriter->SetInput( sortLabelsImageFilter->GetOutput() ); 
  labelWriter->SetFileName( "all_components.mhd" );
  labelWriter->Update();
  labelWriter = 0;


  // Get and display the outer and inner fat sizes and volumes.
  // These steps are based on the assumption that the
  // peripheral fat corresponds to the 2nd largest connected component
  // after the background.

  // Largest non-background component.

  int nComponents = 
    sortLabelsImageFilter->GetNumberOfObjects();
  
  // Muscle region.

  int muscleRegionLabel = 1;
  int muscleSizeinPixels = 
    sortLabelsImageFilter->GetSizeOfObjectsInPixels()[muscleRegionLabel];
  float muscleSizeinPhysicalUnits = 
    sortLabelsImageFilter->GetSizeOfObjectsInPhysicalUnits()[muscleRegionLabel];


  // Display results.

  std::cout << "=======Muscle Quantification results======" 
	    << std::endl;
  std::cout << "Number of spatially connected regions: " 
	    << nComponents
	    << std::endl;
  std::cout << "Muscle component (in pixels) " 
	    << muscleSizeinPixels
	    << ", in Physical Units (mm^3) "
	    << muscleSizeinPhysicalUnits
	    << std::endl;
  std::cout << "=======================================" 
	    << std::endl;


  // Select second largest region (?) and
  // apply threshold to mask out other connected components.
  typedef itk::BinaryThresholdImageFilter<LabelImageType, 
    LabelImageType> 
    Select2ndLargestScoreRegionFilterType;
  Select2ndLargestScoreRegionFilterType::Pointer muscleRegionThresholdFilter = 
    Select2ndLargestScoreRegionFilterType::New();
  muscleRegionThresholdFilter->SetInput( sortLabelsImageFilter->GetOutput() );
  muscleRegionThresholdFilter->SetInsideValue( 255 );
  muscleRegionThresholdFilter->SetOutsideValue( 0 );
  muscleRegionThresholdFilter->SetLowerThreshold( muscleRegionLabel );
  muscleRegionThresholdFilter->SetUpperThreshold( muscleRegionLabel );  // number of regions we want detected.
  muscleRegionThresholdFilter->Update();

  // Save output to a volume.
  labelWriter = 
    itk::ImageFileWriter<LabelImageType>::New();
  labelWriter->SetInput( muscleRegionThresholdFilter->GetOutput() ); 
  labelWriter->SetFileName( outputFileName.c_str() );
  labelWriter->Update();
  labelWriter = 0;

  muscleRegionThresholdFilter = 0;

  return 0;
}

