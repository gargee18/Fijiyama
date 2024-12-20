/*
 * 
 * 
 * 
 */
package io.github.rocsg.fijiyama.gargeetest.ceps;

import java.io.File;
import java.util.Arrays;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.TransformUtils;

public class TransformUsingMetadata {
    
    static String mainDir="/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/0_RAW_data_MRI_2022/";
    static String tiffDir="/home/phukon/Desktop/MRI/02_ceps/2022-03_CEPS_suivi_clinique/1_STACK_16-bits/";
    static String specimen = "318";
    static String[] modality = GeneralUtils.modalities;   
    static String year = GeneralUtils.years[2];
    
    // To get T1(first slice and last slice) and T2(first slice and last slice) metadata,  voxel size, tiff images, and the transformation matrices
    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        System.out.println("Starting test...\n");
        t1t2FinalMatrices();

    }

    public static Object[] t1t2FinalMatrices(){

        ImagePlus imgT1 = IJ.openImage(getPathToT1stack());
        ImagePlus imgT2 = IJ.openImage(getPathToT2stack());
        for(int i = 0; i<modality.length;i++) {
            double vox = getSlicesToPerformVoxelCalculation(modality[i]);
            System.out.println("\nVoxel Size for Modality " + modality[i] + ": " + vox+"\n");

            // Open image and set voxel size
            if (modality[i].equals("T1")) {
                if (imgT1 != null) {
                    double[]voxSizes=VitimageUtils.getVoxelSizes(imgT1);
                    VitimageUtils.printImageResume(imgT1);
                    voxSizes[2]=vox;
                    VitimageUtils.adjustVoxelSize(imgT1, voxSizes);
                }
            }
            if (modality[i].equals("T2")) {
                if (imgT2 != null) {
                    double[]voxSizes=VitimageUtils.getVoxelSizes(imgT2);
                    VitimageUtils.printImageResume(imgT2);
                    voxSizes[2]=vox;
                    VitimageUtils.adjustVoxelSize(imgT2, voxSizes);
                }
            }
                
        }

        double[][] allVectors = getTrMatrices(getPathToT1dcm(), getPathToT2dcm());
        //            0      1      2     3     4     5      6     7    8    9   10    11
        // returns {posT1,oriT1X,oriT1Y,posT2,oriT2X,oriT2Y,t1VX,t1VY,t1VZ,t2VX,t2VY,t2VZ}


        //Take the orthonormal basis (ex,ey,ez) (same for img1 and img2)
        double[] eX = allVectors[1];
        double[] eY = allVectors[2];
        double[] eZ = TransformUtils.vectorialProduct(allVectors[1],allVectors[2]);
        TransformUtils.printVector(eZ, "eZ");

        //Take the pos x1,y1,z1 (guessed expressed in this basis)
        double[] pos1 = allVectors[0];
        //Take the pos x2,y2,z2 (guessed expressed in this basis)
        double[] pos2 = allVectors[3];
        //Compute the delta (x1-x2, y1-y2, z1-z2)
        double[] delta = new double[]{pos1[0]- pos2[0], pos1[1]- pos2[1], pos1[2]- pos2[2]};
        TransformUtils.printVector(delta, "Delta");

        //Convert this delta in real world ()
        //deltaX * ex   (3[]vector)
        //deltaY * ey   (3[]vector)
        //deltaZ * ez   (3[]vector)
        //and sum them pointwise
        double[] deltaRCS = new double[3];
        for (int i = 0; i < 3; i++) {
            deltaRCS[i] = delta[0] * eX[i] + delta[1] * eY[i] + delta[2] * eZ[i];
        }
        TransformUtils.printVector(deltaRCS, "DeltaRCS");


        ItkTransform trVectorTransformT1 = ItkTransform.itkTransformFromDICOMVectors(allVectors[6],allVectors[7],allVectors[8],allVectors[0]);
        ItkTransform trVectorTransformT2 = ItkTransform.itkTransformFromDICOMVectors(allVectors[9],allVectors[10],allVectors[11],allVectors[3]);


        //System.out.println(trVectorTransformT1);
        ItkTransform trT1Inv=new ItkTransform(trVectorTransformT1.getInverse());
        trT1Inv.addTransform(trVectorTransformT2);
        trT1Inv=trT1Inv.simplify();

        trT1Inv.writeMatrixTransformToFile("/home/phukon/Desktop/MRItest/"+specimen+"_T1T2/TransformationMatrixAfterAlignment1.txt");
        ImagePlus resultT2OnT1 = trT1Inv.transformImage(imgT1, imgT2);// Here we try that T1inv o T2 applied to  T1 to register T2 onto T1
        ImagePlus imgFuse=VitimageUtils.compositeNoAdjustOf(imgT1, resultT2OnT1);
        imgFuse.setTitle("trT1Inv");
        imgFuse.show();

        return new Object[]{imgT1, imgT2, trT1Inv};
    }
    


    public static ItkTransform gettingTheItkTransformFromAdicomAsAtrueBelieverInThePowerOfChampagne(String stringPos,String stringOrient){
        String[]coordsPos=stringPos.split("\\\\");
        double[] translation=new double[] {Double.parseDouble(coordsPos[0]),Double.parseDouble(coordsPos[1]),Double.parseDouble(coordsPos[2])};


		//Matrice d orientation 
		String[]coordsOrient=stringOrient.split("\\\\");		
		double[]vectX=new double[] {Double.parseDouble(coordsOrient[0]),Double.parseDouble(coordsOrient[1]),Double.parseDouble(coordsOrient[2])};
		double[]vectY=new double[] {Double.parseDouble(coordsOrient[3]),Double.parseDouble(coordsOrient[4]),Double.parseDouble(coordsOrient[5])};
		double[]vectZ=TransformUtils.vectorialProduct(vectX,vectY);
		ItkTransform tr=ItkTransform.itkTransformFromDICOMVectorsNOAXIAL(vectX,vectY,vectZ,translation);
		System.out.println("In gettingTheItkTra, Matrix = "+tr.drawableString());
		return tr;
        
    }


    
    public static double getSlicesToPerformVoxelCalculation(String mod)
    {
        String inputDir=mainDir+specimen+"/"+mod+"/";
        Object[]objs=getNamesOfSlicesFilesOrderedByDepth(inputDir);
        int nbSlices=(Integer)(objs[2])-(Integer)(objs[0]);
        String pathToLastSlice=(String)objs[3]; //gets name of the image
        String pathToFirstSlice=(String)objs[1];
        Object[] vecs =getTrMatrices(inputDir+pathToFirstSlice, inputDir+pathToLastSlice);
        double[]tabPosFirst=(double[])(vecs[0]);//Getmetadatafrom file in pathToFirstSlice
        double[]tabPosLast=(double[])(vecs[3]);//Getmetadatafrom file in pathToLastSlice
        double[]vectDelta=TransformUtils.vectorialSubstraction(tabPosLast, tabPosFirst);
        double[]vectDeltaPerSlice=divideVector(vectDelta, nbSlices);  
        double voxelSizeZ=TransformUtils.
        norm(vectDeltaPerSlice);
        return voxelSizeZ;
    }

    // To get the T1 DICOM image
    public static String getPathToT1dcm(){  
        //return mainDir+specimen+"/T1/322.MR.RACHIS_LOMBAIRE.6001.1.2022.03.29.21.32.02.571.30590813.dcm";
        return mainDir+specimen+"/T1/318.MR.RACHIS_LOMBAIRE.5001.1.2022.03.29.22.47.09.452.11742060.dcm";
        //return mainDir+specimen+"/T1/"+specimen+".MR.RACHIS_LOMBAIRE.6001.1.2022.03.29.21.28.45.38.16959216.dcm";
    }

    // To get the T2 DICOM image
    public static String getPathToT2dcm(){
        //return mainDir+specimen+"/T2/322.MR.RACHIS_LOMBAIRE.7001.1.2022.03.29.21.31.24.633.13150705.dcm";
        return mainDir+specimen+"/T2/318.MR.RACHIS_LOMBAIRE.6001.1.2022.03.29.22.46.39.609.74920114.dcm";
        //return mainDir+specimen+"/T2/323.MR.RACHIS_LOMBAIRE.5001.1.2022.03.29.21.29.17.66.28442283.dcm";
        
    }
    
    // To get the T1 stack
    public static String getPathToT1stack(){
        return tiffDir+"CEP_"+specimen+"_"+year+"_MRI_T1.tif";
    }
    
    // To get the T2 stack
    public static String getPathToT2stack(){
        return tiffDir+"CEP_"+specimen+"_"+year+"_MRI_T2.tif";
    }
    
    // To get the first and last slices (for year 2022)
    public static Object[]getNamesOfSlicesFilesOrderedByDepth(String inputDir){
        //System.out.println(inputDir);
        File f=new File(inputDir);
        String[]tabNames=f.list();
        if(tabNames.length==1){IJ.showMessage("Some critical fail will happen");}//TODO
        int min=10000000;
        int max=-10000000;
        String strMin="";
        String strMax="";
        for(String name:tabNames){
            int n=getNumberFromWeirdNamingConventionOfDicom(name);
            if(n<min){
                min=n;
                strMin=name;
            }
            if(n>max){
                max=n;
                strMax=name;
            }
        }
        //min is the min number and strMin is the name of the slice with the minimal slice number
        //max is the max number and strMax is the name of the slice with the minimal slice number
        return new Object[]{min,strMin,max,strMax};
    }
    //for year 2024
    public static Object[]getNamesOfSlicesFilesOrderedByDepthYear2024(String inputDir){
        //System.out.println(inputDir);
        File f = new File(inputDir);
        String firstSlice = null;
        String lastSlice = null;
        String[] tabNames = f.list(); // Changed the type to String[]
        if (tabNames != null && tabNames.length > 0) {
        firstSlice = tabNames[0];
        lastSlice = tabNames[tabNames.length-1];
        }
        System.out.println("FS:"+firstSlice);
        System.out.println("LS:"+lastSlice);
        return new Object[]{firstSlice,lastSlice};

    }
    
    public static int getNumberFromWeirdNamingConventionOfDicom(String name){
        String[]tabStr=name.split("\\.");
        int ret=Integer.parseInt(tabStr[4]);
        return ret;
    }
       
    // To get the metadata
    public static String[] getParamsFromDCM(ImagePlus imgPar,String []paramsWanted){
        if(imgPar==null) System.out.println("Warning in TransformUsingMetadata: image used for parameters detection is detected as null image. Computation will fail");
        String paramGlob=imgPar.getInfoProperty();
        String[]paramLines=paramGlob.split("\n");
        String[]ret=new String[paramsWanted.length];
        for(int i=0;i<paramsWanted.length;i++) {
            ret [i]="NOT DETECTED";
            for(int j=0;j<paramLines.length ;j++) {
                if(paramLines[j].split(": ")[0].indexOf(paramsWanted[i])>=0)ret[i]=paramLines[j].split(": ")[1];
            }
        }
        return ret; 
	}

    // To make metadata readable
    public static double[][] metaDataConversion(String[] paramsMetadata){

        String metadata = Arrays.toString(paramsMetadata);
        
        // Remove square brackets
        metadata = metadata.replace("[", "").replace("]", "");

        // Split the string into two separate vectors
        String[] vectors = metadata.split(", ");

        // Extract coordinates for vector 1
        String[] vecPosition = vectors[0].split("\\\\");
        double[] vecPos = new double[3];
        vecPos[0]= Double.parseDouble(vecPosition[0]);
        vecPos[1] = Double.parseDouble(vecPosition[1]);
        vecPos[2]= Double.parseDouble(vecPosition[2]);

        // Extract coordinates for vector 2
        String[] vecOrientationX = vectors[1].split("\\\\");
        double[] vecOriX = new double[3];
        vecOriX[0] = Double.parseDouble(vecOrientationX[0]);
        vecOriX[1] = Double.parseDouble(vecOrientationX[1]);
        vecOriX[2] = Double.parseDouble(vecOrientationX[2]);
        double[] vecOriY = new double[3];
        vecOriY[0] = Double.parseDouble(vecOrientationX[3]);
        vecOriY[1] = Double.parseDouble(vecOrientationX[4]);
        vecOriY[2] = Double.parseDouble(vecOrientationX[5]);
     
        //System.out.println("In metaDataConversion Position Vector: (" + vecPos[0] + ", " + vecPos[1] + ", " + vecPos[2] + ")");
        //System.out.println("Orientation Vector (x1,y1,z1):(" + vecOriX[0] + ", " + vecOriX[1] + ", " + vecOriX[2] + ")");
        //System.out.println("Orientation Vector (x2,y2,z2): (" + vecOriY[0] + ", " + vecOriY[1] + ", " + vecOriY[2] + ")");

        return new double[][] {vecPos,vecOriX,vecOriY};
    }
        

    public static String[]getPosAndOrientationFromDicomAsString(String pathToReferenceDicom){
        ImagePlus imgDcm = IJ.openImage(pathToReferenceDicom);
        String[] paramsWanted = new String[]{"0020,0032","0020,0037"};
        return getParamsFromDCM(imgDcm, paramsWanted);
    }
   


    public static double[][]getTransformMatrixFromDicom(String pathToReferenceDicom){
        ImagePlus imgDcm = IJ.openImage(pathToReferenceDicom);
        String[] paramsWanted = new String[]{"0020,0032","0020,0037"};
        String[] paramsRef = getParamsFromDCM(imgDcm, paramsWanted);
        double[][] matrixT1 = metaDataConversion(paramsRef);
        double[] pos = matrixT1[0];
        double[] oriX = TransformUtils.normalize(matrixT1[1]);   // a 
        double[] oriY = TransformUtils.normalize(matrixT1[2]);   // b
        double[] vX = new double[]{oriX[0],oriY[0]};
        double[] vY = new double[]{oriX[1],oriY[1]};
        double[] vZ = new double[]{oriX[2],oriY[2]};
        return new double[][]{pos,oriX,oriY,vX,vY,vZ};
    }


    // To get the transformation matrices using the metadata of the given images 
    public static double[][] getTrMatrices(String pathtoT1, String pathtoT2){
        double[][]matT1=getTransformMatrixFromDicom(pathtoT1);
        double[][]matT2=getTransformMatrixFromDicom(pathtoT2);
        return new double[][]{
            matT1[0],matT1[1],matT1[2],      matT2[0],matT2[1],matT2[2] ,    matT1[3],matT1[4],matT1[5],      matT2[3],matT2[4],matT2[5]
        };  //posT1,oriT1X,oriT1Y,           posT2,oriT2X,oriT2Y,            t1VX,t1VY,t1VZ,                  t2VX,t2VY,t2VZ} ;
    }

    // To get cross products of vectors
    public static double[] getCrossProd(double[] vec1, double[] vec2){
        double[] crossProd= TransformUtils.vectorialProduct(vec1,vec2); 

        System.out.print("Cross Product: (");
        for (int i = 0; i < crossProd.length; i++) {
            System.out.print(crossProd[i]);
            if (i < crossProd.length - 1) {
                System.out.print(", ");
            }
        }
        System.out.println(")");

        return crossProd;
    }
    
    // To get division
    public static double[] divideVector(double[] v, double divisor) {
        double[] ret = new double[3];
        ret[0] = v[0] / divisor;
        ret[1] = v[1] / divisor;
        ret[2] = v[2] / divisor;
        return ret;
    }



}
