/*
 * This Java code is designed to perform image alignment and registration tasks using the Metadata of dcm files. 
 * It begins by defining arrays of required additional translation values for alignment along the X, Y or Z axes. (Mainly in Z) 
 * 
 * TODO: finish this
 */



package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.FolderOpener;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.TransformUtils;

public class ScriptAlignmentWorkingOn22WithNoRE {
    static double[]requiredAdditionalTranslationForWeirdOneZ=new double[]{
        -21.8        /* 318*/,
        -1.2         /* 322*/,
        0            /* 323*/,
        0            /* 330*/,
        -3.5         /* 335*/,
        -5.4         /* 764B*/,
        0            /* 1181*/,
        -17.4        /* 1193*/,
    };

    static double[]requiredAdditionalTranslationForWeirdOneX=new double[]{
        0           /* 318*/,
        -3.6        /* 322*/,
        0           /* 323*/,
        0           /* 330*/,
        0           /* 335*/,
        0           /* 764B*/,
        0           /* 1181*/,
        0           /* 1193*/,
    };


    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        //int[] indexSpec = {0,1,2,3,4,5,6,7};
        //int indexSpec =1;
        String year="2022";
        String rawDir = "/home/phukon/Desktop/MRI/02_ceps/"+year+"-03_CEPS_suivi_clinique/0_RAW_data_MRI_";
        //String[] specimen=GargeeUtils.getSpecimenList();
        String specimen = "1193";
        //////////Access the data and metadata
        // T1
        //for(int i = 0; i<specimen.length; i++){
            
            System.out.println("Specimen = "+specimen);
            //System.out.println("Index for specimen "+specimen[i]+" is ="+indexSpec[i]);
            Object[]objs = TransformUsingMetadata.getNamesOfSlicesFilesOrderedByDepth( rawDir+year+"/"+specimen+"/T1/"); //2022
            //objs = TransformUsingMetadata.getNamesOfSlicesFilesOrderedByDepthYear2024(rawDir+year+"/"+specimen+"/T1/"); //2024
            String strFirstT1=rawDir+year+"/"+specimen+"/T1/"+(String)objs[1];//2022
            //String strFirstT1=rawDir+year+"/"+specimen+"/T1/"+(String)objs[0];//2024
            ImagePlus imgT1 = FolderOpener.open(rawDir+year+"/"+specimen+"/T1/", "");

            // T2
            objs = TransformUsingMetadata.getNamesOfSlicesFilesOrderedByDepth(rawDir+year+"/"+specimen+"/T2/");//2022
            //objs = TransformUsingMetadata.getNamesOfSlicesFilesOrderedByDepthYear2024(rawDir+year+"/"+specimen+"/T2/");//2024
            String strFirstT2=rawDir+year+"/"+specimen+"/T2/"+(String)objs[1];//2022
            //String strFirstT2=rawDir+year+"/"+specimen+"/T2/"+(String)objs[0];//2024
            ImagePlus imgT2 = FolderOpener.open(rawDir+year+"/"+specimen+"/T2/", "");

            //Gathering metadata (vx, vy and pos)
            String[]strTabFirstT1=getPosAndOrientationSliceThicknessAndPixSizeFromDicomAsString(strFirstT1);
            String[]strTabFirstT2=getPosAndOrientationSliceThicknessAndPixSizeFromDicomAsString(strFirstT2);

            /////////// Process the metadata
            //Compute and set the voxel size
            double voxSZT1=Double.parseDouble(strTabFirstT1[2]);
            double voxSXT1=Double.parseDouble(strTabFirstT1[3].split("\\\\")[0]);
            VitimageUtils.adjustVoxelSize(imgT1,new double[]{voxSXT1,voxSXT1,voxSZT1});

            double voxSZT2=Double.parseDouble(strTabFirstT2[2]);
            double voxSXT2=Double.parseDouble(strTabFirstT2[3].split("\\\\")[0]);
            VitimageUtils.adjustVoxelSize(imgT2,new double[]{voxSXT2,voxSXT2,voxSZT2});

            //Read the position vectors
            String[]strTabPosFirst=strTabFirstT1[0].split("\\\\");
            double[]tabPosFirstT1=new double[]{Double.parseDouble(strTabPosFirst[0]),Double.parseDouble(strTabPosFirst[1]),Double.parseDouble(strTabPosFirst[2])};

            strTabPosFirst=strTabFirstT2[0].split("\\\\");
            double[]tabPosFirstT2=new double[]{Double.parseDouble(strTabPosFirst[0]),Double.parseDouble(strTabPosFirst[1]),Double.parseDouble(strTabPosFirst[2])};

            //Read the orientation vectors and compute the third one
            String[]strTabOrientFirst=strTabFirstT1[1].split("\\\\");
            double[]vectXT1=new double[]{Double.parseDouble(strTabOrientFirst[0]),Double.parseDouble(strTabOrientFirst[1]),Double.parseDouble(strTabOrientFirst[2])};
            double[]vectYT1=new double[]{Double.parseDouble(strTabOrientFirst[3]),Double.parseDouble(strTabOrientFirst[4]),Double.parseDouble(strTabOrientFirst[5])};
            double[]vectZT1=TransformUtils.vectorialProduct(vectXT1, vectYT1);

            strTabOrientFirst=strTabFirstT2[1].split("\\\\");
            double[]vectXT2=new double[]{Double.parseDouble(strTabOrientFirst[0]),Double.parseDouble(strTabOrientFirst[1]),Double.parseDouble(strTabOrientFirst[2])};
            double[]vectYT2=new double[]{Double.parseDouble(strTabOrientFirst[3]),Double.parseDouble(strTabOrientFirst[4]),Double.parseDouble(strTabOrientFirst[5])};
            double[]vectZT2=TransformUtils.vectorialProduct(vectXT2, vectYT2);

            
            //Assemble the matrices with respect to the inversion of y z
            double[]array16elementsT1=new double[]{
                vectXT1[0],    vectYT1[0],  vectZT1[0],  tabPosFirstT1[0],
                -vectXT1[1],  -vectYT1[1], -vectZT1[1], -tabPosFirstT1[1],
                -vectXT1[2],  -vectYT1[2], -vectZT1[2], -tabPosFirstT1[2]
            };
            ItkTransform trT1=ItkTransform.array16ElementsToItkTransform(array16elementsT1);

            double[]array16elementsT2=new double[]{
                vectXT2[0],    vectYT2[0],  vectZT2[0],  tabPosFirstT2[0],
                -vectXT2[1],  -vectYT2[1], -vectZT2[1], -tabPosFirstT2[1],
                -vectXT2[2],  -vectYT1[2], -vectZT2[2], -tabPosFirstT2[2]
            };
            ItkTransform trT2=ItkTransform.array16ElementsToItkTransform(array16elementsT2);

            ItkTransform trFromT2toT1=new ItkTransform(trT2.getInverse());
            trFromT2toT1.addTransform(trT1);

            // ItkTransform trAddAlignementZ=ItkTransform.array16ElementsToItkTransform(new double[]{1,0,0,requiredAdditionalTranslationForWeirdOneX[indexSpec[i]],0,1,0,0,0,0,1,requiredAdditionalTranslationForWeirdOneZ[indexSpec[i]]});
            // trFromT2toT1.addTransform(trAddAlignementZ);

            ImagePlus resultT2OnT1 = trFromT2toT1.transformImage(imgT1, imgT2);// Here we try that T1inv o T2 applied to  T1 to register T2 onto T1
            ImagePlus imgFuse=VitimageUtils.compositeNoAdjustOf(imgT1, resultT2OnT1);
            imgFuse.setTitle(specimen);
            imgFuse.show();
        
            //Prepare a BlockmMatchingRegistration bm between imgRef and imgInit, with trFromT2toT1 as init tr,

            ItkTransform trTranslation = RegistrationAutomate.autoLinearRegistration(imgT1,imgT2,trFromT2toT1,specimen);
            trTranslation.writeMatrixTransformToFile("/home/phukon/Desktop/MRItest/regTranslation/"+specimen+"_transformation.txt");
            ImagePlus finalResult=trTranslation.transformImage(imgT1,imgT2);
            IJ.saveAsTiff(finalResult,"/home/phukon/Desktop/MRItest/regTranslation/"+specimen+"_translated.tif");
            ImagePlus imgComposite=VitimageUtils.compositeNoAdjustOf(imgT1, finalResult);
            IJ.saveAsTiff(imgComposite, "/home/phukon/Desktop/MRItest/regTranslation/"+specimen+"_translatedComposite.tif");
            imgComposite.show();
   
        }
    
    
 
   

    public static String[]getPosAndOrientationSliceThicknessAndPixSizeFromDicomAsString(String pathToReferenceDicom){
        ImagePlus imgDcm = IJ.openImage(pathToReferenceDicom);
        String[] paramsWanted = new String[]{"0020,0032","0020,0037","0018,0050","0028,0030"};
        return TransformUsingMetadata.getParamsFromDCM(imgDcm, paramsWanted);
    }

}



