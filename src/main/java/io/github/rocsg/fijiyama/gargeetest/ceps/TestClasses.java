package io.github.rocsg.fijiyama.gargeetest.ceps;

import java.awt.Image;
import java.util.Arrays;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.io.ImportDialog;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import ij.plugin.filter.ImageMath;

public class TestClasses {
    
   
    public static String year = "2022"; 

    public static void main(String[] args) {
        ImageJ ij = new ImageJ();

        String[] specimen = GeneralUtils.specimenListXR;
        String mainDir = "/home/phukon/Desktop/";
        String outDir = mainDir + "cropdata"+year +"/";
        String rawData = mainDir+"XR_"+year+"_raw/";
        String cropCoords = mainDir + year+"_cropcoords.csv";
        String sizeOfCroppedImage = mainDir + "croppedImageDimensions" + year + ".csv";
        
        ImagePlus img = new ImagePlus("/home/phukon/Desktop/XR_2023_raw/CEP_1266A_2023_XR.tif");
        //getCoordinatesandCrop(rawData, outDir, specimen, cropCoords, sizeOfCroppedImage);
        Roi roi = new Roi(113, 109, 312, 304);
        ImagePlus imgc = CropImagesWithRoi.crop(img, roi, "stack");
        imgc.show();
        IJ.saveAsTiff(imgc,"/home/phukon/Desktop/cropdata2023/CEP_1266A_2023_XR_crop.tif");

        
       
    
    }

    public static ImagePlus getCoordinatesandCrop(String inputDir, String outputDir, String[] specimen, String pathToCSVFileforCropCoords, String pathtoCSVFileforCropSize){
        int numberOfSpecimens = specimen.length;
        int numberOfCoordinates = 3; //x,y,z
        ImagePlus imgCropped = null;
        int[][] trCoords = new int[numberOfSpecimens+1][numberOfCoordinates+1];
        int[][] cropImgdim = new int[numberOfSpecimens+1][numberOfCoordinates+1];
        
        String[][] cropCoordsData = VitimageUtils.readStringTabFromCsv(pathToCSVFileforCropCoords);
        String[][] sizeCropData = VitimageUtils.readStringTabFromCsv(pathtoCSVFileforCropSize);


        for (int i = 1; i < numberOfSpecimens+1; i++) {
                trCoords[i][1] = (int)Double.parseDouble(cropCoordsData[i][1]); // x coordinate
                trCoords[i][2] = (int)Double.parseDouble(cropCoordsData[i][2]); // y coordinate
                trCoords[i][3] = (int)Double.parseDouble(cropCoordsData[i][3]); // z coordinate


                cropImgdim[i][1] = (int)Double.parseDouble(sizeCropData[i][1]); // x coordinate
                cropImgdim[i][2] = (int)Double.parseDouble(sizeCropData[i][2]); // y coordinate
                cropImgdim[i][3] = (int)Double.parseDouble(sizeCropData[i][3]); // z coordinate

        }
        System.out.println("Crop coordinates:");
        for (int i = 1; i < numberOfSpecimens+1; i++) {
            System.out.println("Specimen " + specimen[i-1] + ": x = " + trCoords[i][1] + ", y = " + trCoords[i][2] + ", z = " + trCoords[i][3]);

        }


        System.out.println("Crop Dimensions:");
        for (int i = 1; i < numberOfSpecimens+1; i++) {
            System.out.println("Specimen " + specimen[i-1] + ": x = " + cropImgdim[i][1] + ", y = " + cropImgdim[i][2] + ", z = " + cropImgdim[i][3]);
        }

        // Cropping 
        for (int i = 1; i < specimen.length+1; i++){
            ImagePlus img = new ImagePlus(inputDir + "CEP_" + specimen[i-1] + "_" + year + "_XR_unsigned.tif");
            imgCropped = cropImg(img, trCoords[i][1], trCoords[i][2], trCoords[i][3], cropImgdim[i][1], cropImgdim[i][2], cropImgdim[i][3]);
            IJ.saveAsTiff(imgCropped, outputDir + "CEP_" + specimen[i-1] + "_" + year + "_XR_crop_us.tif");
        }

        return imgCropped;
    }
    
    public static ImagePlus cropImg(ImagePlus img, int x0, int y0 , int z0, int dimXX, int dimYY, int dimZZ){
        
        ImagePlus imgCrop =cropImageShort(img,x0,y0,z0,dimXX,dimYY,dimZZ); 
        return imgCrop;
    }

    public static void imgDimensionsforCroppedData(String inputDir, String[] specimen){

        int nbRowsCsv=1+GeneralUtils.specimenListXR.length;
        int nbColumnsCsv=4;
        String[][]strTab=new String[nbRowsCsv][nbColumnsCsv];
        strTab[0][0] = "CEP";
        strTab[0][1] = "X";
        strTab[0][2] = "Y";
        strTab[0][3] = "Z";
        for (int i = 0; i < specimen.length; i ++){
            String imgPath = inputDir+"cep_"+specimen[i]+"/raw/"+specimen[i]+"_"+year+"_crop.tif";
            System.out.println(imgPath);
            ImagePlus img= new ImagePlus(imgPath);
            int[] dim = VitimageUtils.getDimensions(img);
            System.out.println("Parameters of specimen "+specimen[i]+ ": " +Arrays.toString(dim));
            strTab[i + 1][0] = specimen[i];
            for(int j=0; j<dim.length;j++){
                strTab[i+1][j+1] = Double.toString(dim[j]);
            }
        }
        VitimageUtils.writeStringTabInCsv2(strTab,"/home/phukon/Desktop/croppedImageDimensions2024.csv");
    }

    public static void imgVoxSize(String inputDir, String[] specimen){

        int nbRowsCsv=1+GeneralUtils.specimenListXR.length;
        int nbColumnsCsv=4;
        String[][]strTab=new String[nbRowsCsv][nbColumnsCsv];
        strTab[0][0] = "CEP";
        strTab[0][1] = "X";
        strTab[0][2] = "Y";
        strTab[0][3] = "Z";
        for (int i = 0; i < specimen.length; i ++){
            String imgPath = inputDir + "CEP_"+specimen[i]+"_"+year+"_XR.tif";
            ImagePlus img= new ImagePlus(imgPath);
            double[]voxs=VitimageUtils.getVoxelSizes(img);
            System.out.println("Parameters of specimen "+specimen[i]+ ": " +Arrays.toString(voxs));
            strTab[i + 1][0] = specimen[i];
            for(int j=0; j<voxs.length;j++){
                strTab[i+1][j+1] = Double.toString(voxs[j]);
            }
        }
        VitimageUtils.writeStringTabInCsv2(strTab,inputDir+".csv");
    }

    public static ImagePlus createImageShort(ImagePlus imgCropped, String title, int width, int height, int slices, int options) {

        ImageStack stack = new ImageStack(width, height);

        for (int i = 0; i < slices; i++) {
            short[] pixels = new short[width * height]; // Each slice is an array of pixels
            ImageProcessor ip = new ShortProcessor(width, height, pixels, null);
            stack.addSlice("Slice " + (i + 1), ip);
        }

        ImagePlus imp = new ImagePlus(title, stack);
        return imp;
    }

    public static ImagePlus cropImageShort(ImagePlus img,int x0,int y0,int z0,int dimXX,int dimYY,int dimZZ) {
		if(img.getType()!=ImagePlus.GRAY16)return null;
		int[]dims=VitimageUtils.getDimensions(img);
		int X=dims[0];int Y=dims[1];int Z=dims[2];
		int xm=(int)Math.max(0, x0);
		int xM=(int)Math.min(xm+dimXX-1,X-1);
		int dimX=xM-xm+1;
		int ym=(int)Math.max(0, y0);
		int yM=(int)Math.min(ym+dimYY-1,Y-1);
		int dimY=yM-ym+1;
		int zm=(int)Math.max(0, z0);
		int zM=(int)Math.min(zm+dimZZ-1,Z-1);
		int dimZ=zM-zm+1;		
			
		ImagePlus out=createImageShort(img,"CroppedImage", dimX, dimY, dimZ, ij.gui.NewImage.FILL_BLACK);	

		//VitimageUtils.adjustImageCalibration(out, img);

		for(int z=zm;z<zm+dimZ;z++) {
			short[] valsImg=(short[])img.getStack().getProcessor(z+1).getPixels();
			short[] valsOut=(short[])out.getStack().getProcessor(z-zm+1).getPixels();
			for(int x=xm;x<xm+dimX;x++) {
				for(int y=ym;y<ym+dimY;y++){
					valsOut[dimX*(y-ym)+(x-xm)]=((short)(valsImg[X*y+x]));
				}			
			}
			out.getStack().setSliceLabel(img.getStack().getSliceLabel(z+1), z-zm+1);
		}
		VitimageUtils.transferProperties(img, out);
		return out;
	}

  
    public static ImagePlus getRegisteredImage(String inputDir, String[] specimen){
        
        ImagePlus imgRegistered = new ImagePlus();

        for (int i = 0; i < specimen.length; i ++){
            RegistrationHighRes.getHighResTransform(specimen[i]);
            // Get images year 1 & year 2
            ImagePlus imgOriginal = new ImagePlus(inputDir+"CEP_"+specimen[i]+"_"+year+"_XR_us.tif");
            ImagePlus imgCropped = new ImagePlus(inputDir+"CEP_"+specimen[i]+"_"+year+"_XR_crop_us.tif");

            // Get  transformation matrices for low res
            //double trMatYear1 = 

            // Get crop coordinates from csv


            // Multiply transformations


            //Save tr matrix
            String pathtosave = "/home/phukon/Desktop/cropdata2022/"+specimen[i]+"_cropped.tif";


            // save composite
            IJ.saveAsTiff(imgCropped,pathtosave);
        }

        return imgRegistered;
        
    }

    


}