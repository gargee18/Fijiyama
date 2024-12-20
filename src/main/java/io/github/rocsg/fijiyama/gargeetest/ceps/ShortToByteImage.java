package io.github.rocsg.fijiyama.gargeetest.ceps;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Duplicator;
import ij.process.ByteProcessor;
import io.github.rocsg.fijiyama.common.VitimageUtils;


public class ShortToByteImage {

    public static void main2(String[] args) {
        String[] specimen = GeneralUtils.specimenListXR;
        String year = "2022";

        ImageJ ij=new ImageJ();

        String dirdir = "/home/phukon/Desktop/Registered_2022_2023_2024_highres_UpdateJuly/";
        ImagePlus img=IJ.openImage(dirdir+"CEP_313B_2022_XR.tif");
        int targetValue=160;

        
    }

    public static void main(String[] args) {
        
    }



    //We want that min becomes 0 and max becomes 255
    public static ImagePlus way1ToConvert(ImagePlus img,double min, double max,double targetValueOfMax){
        ImagePlus ret=VitimageUtils.imageCopy(img);
        IJ.run(ret,"32-bit","");
        ret.setDisplayRange(min, max*(255.0/targetValueOfMax));
        //VitimageUtils.convertToGray8(ret);
        IJ.run(ret,"8-bit","");
        return ret;
    }


    //We want that min becomes 0 and max becomes 255
    public static ImagePlus way2ToConvert(ImagePlus img,double min, double max,double targetValueOfMax){
        ImagePlus ret=VitimageUtils.imageCopy(img);
        ret=VitimageUtils.convertToFloat(ret);

        //add a normal random noise with an amplitude that is (max-min)/500;
        double stdNoise=(max-min)/1000.0;
        IJ.run(ret, "Add Specified Noise...", "stack standard="+stdNoise);

        ret=VitimageUtils.makeOperationOnOneImage(ret, 1, -min, false);
        //After this operation, min becomes 0, and max becomes (max-min)


        //Now we would like (max-min) becomes 255
        ret=VitimageUtils.makeOperationOnOneImage(ret,2, targetValueOfMax/(max-min), false);
        if(ret.getType()==ImagePlus.GRAY8)return ret;
        else if(ret.getType()==ImagePlus.GRAY16)return convertShortToByteWithoutDynamicChanges(ret);
        else if(ret.getType()==ImagePlus.GRAY32)return VitimageUtils.convertFloatToByteWithoutDynamicChanges(ret);
        else{IJ.showMessage("Are you tryingg to convert a RGB image with a min max ? This is not good !!!");}
        return ret;
    }



    public static ImagePlus way3ToConvert(ImagePlus ip,double min, double max, double targetValue) {
        if(targetValue<1)targetValue=1;
        if(targetValue>255)targetValue=255;
        int width = ip.getWidth();
        int height = ip.getHeight();
        int nSlices = ip.getStackSize();
        ImageStack stack16 = ip.getStack();
        ImageStack stack8 = new ImageStack(width, height);

        //Scale factor required to jump from a [min-max] range to [0 - 160] range
        double scale = targetValue*1.0 / (max - min);
        //System.out.println(scale);

        // Convert each slice in the stack
        for (int s = 1; s <= nSlices; s++) {
            short[] pixels16 = (short[]) stack16.getProcessor(s).getPixels();
            byte[] pixels8 = new byte[width * height];

            for (int i = 0; i < pixels16.length; i++) {
                int value = (int) (((int)(pixels16[i] & 0xffff)) - min);
                if (value < 0) value = 0;
                value = (int) Math.round(value * scale);
                if (value > 255) value = 255;
                pixels8[i] = (byte) value;
            }

            ByteProcessor byteProcessor = new ByteProcessor(width, height, pixels8, null);
            stack8.addSlice(stack16.getSliceLabel(s), byteProcessor);
        }

        return new ImagePlus("8-bit Image", stack8);
    }

    public static ImagePlus convertShortToByteWithoutDynamicChanges(ImagePlus imgIn) {
        ImagePlus ret=new Duplicator().run(imgIn);
        if(imgIn.getType()==ImagePlus.GRAY8)return ret;
        ret=IJ.createImage("", "8-bit", ret.getWidth(), ret.getHeight(), ret.getNChannels(), ret.getNSlices(), ret.getNFrames());
        short[][] in=new short[imgIn.getStackSize()][];
        byte[][] out=new byte[ret.getStackSize()][];
        int X=imgIn.getWidth();
        int res=0;
        for(int z=0;z<imgIn.getStackSize();z++) {
            in[z]=(short []) imgIn.getStack().getProcessor(z+1).getPixels();
            out[z]=(byte []) ret.getStack().getProcessor(z+1).getPixels();

            for(int x=0;x<imgIn.getWidth();x++) {
                for(int y=0;y<imgIn.getHeight();y++) {
                    res=((int)(Math.round(in[z][y*X+x])));
                    if(res<0)res=0;
                    if(res>255)res=255;
                    out[z][y*X+x]=(byte)res;
                }			
            }
            ret.getStack().setSliceLabel(imgIn.getStack().getSliceLabel(z+1), z+1);
        }
        VitimageUtils.transferProperties(imgIn, ret);
        
        return ret;
    }



    public static ImagePlus convertShortToByteStack(ImagePlus ip,double min, double max, double targetValue) {
        if(targetValue<1)targetValue=1;
        if(targetValue>255)targetValue=255;
        int width = ip.getWidth();
        int height = ip.getHeight();
        int nSlices = ip.getStackSize();
        ImageStack stack16 = ip.getStack();
        ImageStack stack8 = new ImageStack(width, height);

        //Scale factor required to jump from a [min-max] range to [0 - 160] range
        double scale = targetValue / (max - min);
        //System.out.println(scale);

        // Convert each slice in the stack
        for (int s = 1; s <= nSlices; s++) {
            short[] pixels16 = (short[]) stack16.getProcessor(s).getPixels();
            byte[] pixels8 = new byte[width * height];

            for (int i = 0; i < pixels16.length; i++) {
                int value = (int) (((int)(pixels16[i] & 0xffff)) - min);
                if (value < 0) value = 0;
                value = (int) (value * scale);
                if (value > 255) value = 255;
                pixels8[i] = (byte) value;
            }

            ByteProcessor byteProcessor = new ByteProcessor(width, height, pixels8, null);
            stack8.addSlice(stack16.getSliceLabel(s), byteProcessor);
        }

        return new ImagePlus("8-bit Image", stack8);
    }
    
}
    

        
 
    

