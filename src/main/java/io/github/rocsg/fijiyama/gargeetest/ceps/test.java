package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.HyperStackConverter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.registration.ItkTransform;

public class test {


    public static void main(String[] args) {
        String[] specimen = GeneralUtils.specimenListXR;
        String year = "2024";
        ImageJ ij=new ImageJ();
        int targetValue=163;

        String dirdir = "/home/phukon/Desktop/Registered_2022_2023_2024_highres_UpdateJuly/";
    
        String[][]tabData=VitimageUtils.readStringTabFromCsv("/home/phukon/Downloads/ROI measurements (all ceps in 3 years) - Pixel intensities for different classes.csv");

        for (int i=1;i<2/*specimen.length*/;i++){
            ImagePlus img=IJ.openImage(dirdir+"CEP_"+specimen[i]+"_"+year+"_XR.tif");
            System.out.println(dirdir+"CEP_"+specimen[i]+"_"+year+"_XR.tif");
            int specimenNumber=GeneralUtils.getIndexFromXRSpecimenName(specimen[i]);
            int startingrowNumber=1+27*specimenNumber;
            int bgRow22=startingrowNumber+19;
            int fgRow22=startingrowNumber+1;
            int bgRow23=startingrowNumber+22;
            int fgRow23=startingrowNumber+4;
            int bgRow24=startingrowNumber+25;
            int fgRow24=startingrowNumber+7;

            int min22=(int) Double.parseDouble(tabData[bgRow22][8]);
            int max22=(int) Double.parseDouble(tabData[fgRow22][8]);

            int min23=(int) Double.parseDouble(tabData[bgRow23][8]);
            int max23=(int) Double.parseDouble(tabData[fgRow23][8]);

            int min24=(int) Double.parseDouble(tabData[bgRow24][8]);
            int max24=(int) Double.parseDouble(tabData[fgRow24][8]);


            //System.out.println("The min and max number I got for "+specimen[i]+" are : ["+min24+"],["+max24+"]");

            VitimageUtils.printImageResume(img);
            System.out.println("Val of min "+min24);
            System.out.println("Val of max "+max24);
            System.out.println("Val of target "+targetValue);

            ImagePlus imgByte = ShortToByteImage.way2ToConvert(img, min24, max24, targetValue);
            img = imgByte;
            int[]tabFirstNoFuck=new int[256];
            byte[] dataR;
            for(int z=0;z<img.getStackSize();z++) {
                dataR=(byte[])(img.getStack().getProcessor(z+1).getPixels());
            for(int x=0;x<img.getWidth();x++) {
                for(int y=0;y<img.getHeight();y++){
                    int val=(int)(dataR[x+img.getWidth()*y] & 0xff);
                    tabFirstNoFuck[val]++;
                    
                }
            }
        }
    int countNoWarn=0;
    int countWarn=0;
    for(int j=1;j<255;j++){
            double warningArea=0.1;
            double delta=(tabFirstNoFuck[j])*2.0/(tabFirstNoFuck[j-1]+tabFirstNoFuck[j+1]);
            if(Math.abs(delta-1.0)>warningArea){
                System.out.println("Warning : "+delta + " at index="+j);
                countWarn++;
            }
            else countNoWarn++;
//            System.out.println(delta);

            //ImagePlus imgByte = ShortToByteImage.convertShortToByteStack(img,min24, max24, targetValue);
            
            //convertShortToByteWithoutDynamicChanges(imgByte);
            IJ.saveAsTiff(imgByte, "/home/phukon/Desktop/nomalizedToByte/CEP_"+specimen[i]+"_"+year+"_XR.tif");
        }
        System.out.println("Warnings: "+countWarn+" \n no warnings : "+countNoWarn);
    }

    
    }
        
        






        public static ImagePlus makeCompositeandHyperstack(ImagePlus imgYear2022, ImagePlus imgYear2023, ImagePlus imgYear2024){

            

            ImagePlus hyperImage = new ImagePlus();
            hyperImage = VitimageUtils.compositeRGBDoubleJet(imgYear2022,imgYear2023,imgYear2024, "Composite 3 years", true,2);
            hyperImage.show();
            hyperImage = HyperStackConverter.toHyperStack(hyperImage, 0, 0, 0);

            return hyperImage;
        }

        public static void autoregister(String specimen, String year1, String year2, ImagePlus imgRef, ImagePlus imgMov, ItkTransform trInit){
            ItkTransform trFinal = RegistrationAutomate.autoLinearRegistration(imgRef,imgMov, trInit,specimen);
            trFinal.writeMatrixTransformToFile("/home/phukon/Desktop/highresTransformationMatricesAndRegisteredImages/registered/"+specimen+"_transform_global_img_moving_"+year1+"_"+year2+".txt");
            ImagePlus finalResult = trFinal.transformImage(imgRef,imgMov);
            IJ.saveAsTiff(finalResult, "/home/phukon/Desktop/highresTransformationMatricesAndRegisteredImages/registered/"+specimen+"_Registered_"+year1+"_"+year2+".tif");

        }
}


   
        // for(int i = 0; i<specimen.length;i++){
        //     ImagePlus img=IJ.openImage(dirdir+"CEP_"+specimen[i]+"_"+year+"_XR.tif");
        //     System.out.println(dir+"CEP_"+specimen[i]+"_"+year+"_XR.tif");
        //     System.out.println(img.getCalibration());
        //     img.getCalibration().setFunction(Calibration.STRAIGHT_LINE, new double[]{0,1}, "Gray Value");
        // //ItkTransform tr=ItkTransform.itkTransformFromCoefs(new double[]{1,0,0,0,0,1,0,0,0,0,1,0});

        
        // //ImagePlus imgRes1=tr.transformImage(img, img);

        //     //img.show();
        //     //IJ.saveAsTiff(img, dirdir+"CEP_"+specimen[i]+"_"+year+"_XR.tif");
        // }
