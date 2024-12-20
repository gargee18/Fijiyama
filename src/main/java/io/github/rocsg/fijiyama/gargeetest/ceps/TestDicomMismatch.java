package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class TestDicomMismatch {
    
    // Test to display dicom tags for all specimens and also calculate the voxel size 
    public static void displayAllRelevantInformationForInvestigation(){
        //String[] paramsWanted = new String[]{"0020,0032","0020,0037","0018,0050","0028,0030","0018,1030", "0020,1041", "0018,1250", "0018,1251", "0018,1314", "0051,100C"};
        String[] paramsTranslation = new String[]{"Position", "Orientation", "Slice Thickness", "Pixel Spacing", "Protocol Name", "Slice Location", "Transmitting Coil", "Receiving Coil", "Flip angle", "FoV" };
        String year="2022";
        String rawDir = "/home/phukon/Desktop/MRI/02_ceps/"+year+"-03_CEPS_suivi_clinique/0_RAW_data_MRI_";
        Object[]objs=null;
        int nbRowsCsv=1+2*GeneralUtils.specimenListMRI.length;
        int nbColumnsCsv=1+GeneralUtils.paramsWanted.length;
        String[][]strTab=new String[nbRowsCsv][nbColumnsCsv];
        for(int i=0;i<GeneralUtils.specimenListMRI.length;i++){strTab[1+2*i][0]=GeneralUtils.specimenListMRI[i]+"_T1";strTab[1+2*i+1][0]=GeneralUtils.specimenListMRI[i]+"_T2";};
        for(int j=0;j<GeneralUtils.paramsWanted.length;j++)strTab[0][1+j]=paramsTranslation[j];


        for(int i=0;i<GeneralUtils.specimenListMRI.length;i++){
            String specimen=GeneralUtils.specimenListMRI[i];
            for(String mod : GeneralUtils.modalities){
                int rowIndex=1+ i*2+(  (mod.contains("T1")) ? 0 : 1);
                objs = TransformUsingMetadata.getNamesOfSlicesFilesOrderedByDepth( rawDir+year+"/"+specimen+"/"+mod+"/");
                String str=rawDir+year+"/"+specimen+"/"+mod+"/"+(String)objs[1];
                ImagePlus imgDcm = IJ.openImage(str);
                String[]strMetaData = TransformUsingMetadata.getParamsFromDCM(imgDcm, GeneralUtils.paramsWanted);
                System.out.print("Metadata of specimen "+specimen +": \n");
                for (int j = 0; j < strMetaData.length; j++) {
                    strTab[rowIndex][j+1]=strMetaData[j];
                    System.out.println(paramsTranslation[j]+": "+strMetaData[j]);
                }
                double vox = TransformUsingMetadata.getSlicesToPerformVoxelCalculation(mod);
                System.out.print("Voxel size "+VitimageUtils.dou(vox));
                System.out.println("\n\n");
            }
        }
        String pathToTheCsv="/home/phukon/Desktop/Temp/test.csv";
        VitimageUtils.writeStringTabInCsv2(strTab,pathToTheCsv);
    }

    public static void main(String[] args) {
        displayAllRelevantInformationForInvestigation();
    }


    public static String[]getPosAndOrientationSliceThicknessAndPixSizeFromDicomAsString(String pathToReferenceDicom){
        ImagePlus imgDcm = IJ.openImage(pathToReferenceDicom);
        String[] paramsWanted = new String[]{"0020,0032","0020,0037","0018,0050","0028,0030"};
        return TransformUsingMetadata.getParamsFromDCM(imgDcm, paramsWanted);
    }

}

