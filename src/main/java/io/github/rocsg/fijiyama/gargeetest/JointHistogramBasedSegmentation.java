package io.github.rocsg.fijiyama.gargeetest;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import io.github.rocsg.fijiyama.common.VitimageUtils;

public class JointHistogramBasedSegmentation {

    public static void main(String[] args) {
        ImageJ ij=new ImageJ();
        String imgPathR="/home/phukon/Desktop/test_quantification/313B_imgRef.tif";//Red channel
        String imgPathG="/home/phukon/Desktop/test_quantification/313B_TransformedImage.tif";//Green channel
        ImagePlus imgR=IJ.openImage(imgPathR);
        ImagePlus imgG=IJ.openImage(imgPathG);
        ImagePlus imgComposite = VitimageUtils.compositeNoAdjustOf(imgR, imgG);
        imgComposite.show();

        ImagePlus[]jointHistograms=computeAndDisplayJointHistograms(imgR,imgG);
        ImagePlus segmentation=makeSegmentationBasedOnRoiInputFromTheUser(imgR,imgG,jointHistograms[0],jointHistograms[1]);
        segmentation.setTitle("Segmentation");
        IJ.run(segmentation,"Fire","");
        segmentation.show();
        // ImagePlus imgSegComp = VitimageUtils.compositeNoAdjustOf(imgComposite,imgG);
        // imgSegComp.setTitle("Composite of Segmented and Registered");
        // imgSegComp.show();

        }



        public static ImagePlus makeSegmentationBasedOnRoiInputFromTheUser(ImagePlus imgR,ImagePlus imgG,ImagePlus histo,ImagePlus histoLog){
                
            //Ask for three rois
            RoiManager roiManager=RoiManager.getRoiManager();
            int nbAreas=3;
            boolean finished=false;

            do {
                VitimageUtils.waitFor(1000);
                if(roiManager.getCount()==nbAreas)finished=true;
                    System.out.println("Waiting "+nbAreas+". Current number="+roiManager.getCount());
                }while (!finished);

            Roi roi1=roiManager.getRoi(0); //index for first polygon
            Roi roi2=roiManager.getRoi(1); //index for second polygon
            Roi roi3=roiManager.getRoi(2); //index for thrid polygon
            

            // Get access to original images (imgR and imgG) as dataR[] and dataG[]
            int X = imgR.getWidth();
            int Y = imgR.getHeight();
            int Z = imgR.getStackSize();

            // Initialize arrays to store pixel data
            byte[] dataR;
            byte[] dataG;

            // Create a destination image (float), of the size of the original image
            ImagePlus segmentation=VitimageUtils.nullImage(imgG);
            segmentation=VitimageUtils.convertByteToFloatWithoutDynamicChanges(segmentation);

            //Take access to the pixels of segmentation as dataSegmentation[]
            for(int z=0; z<Z; z++ ){
                float[] dataSegmentation = (float[]) segmentation.getStack().getProcessor(z + 1).getPixels();
                //For each pixl (x,y,z)
                //Collect the value of imgR and imgG. It yields equivalent coordinates within the histogram (x,y)
                dataR = (byte[]) imgR.getStack().getProcessor(z + 1).getPixels();
                dataG = (byte[]) imgG.getStack().getProcessor(z + 1).getPixels();
                for (int y = 0; y < Y; y++) {
                    for (int x = 0; x < X; x++) {
                        // Collect the pixel values of imgR and imgG
                        int valR=((int)(dataR[x+X*y] & 0xff));
                        int valG=((int)(dataG[x+X*y] & 0xff));

                        int xHisto=valR;
                        int yHisto=255-valG;
                        
                        //For each roi, check if the point (x,y) is contained by the Roi
                        // Check if the pixel is contained by ROI 1
                        if (roi1.contains(xHisto, yHisto)) {
                            dataSegmentation[x+X*y] = 1;
                        }
                        // Check if the pixel is contained by ROI 2
                        else if (roi2.contains(xHisto, yHisto)) {
                            dataSegmentation[x+X*y] = 2;
                        }
                        // Check if the pixel is contained by ROI 3
                        else if (roi3.contains(xHisto, yHisto)) {
                            dataSegmentation[x+X*y] = 3;
                        }
                        // If the pixel is not contained by any ROI
                        else {
                            dataSegmentation[x+X*y] = 0;
                        }
                    }
                }
               
            }
            segmentation.setDisplayRange(0, 4);

            return segmentation;
        }
    

        public JointHistogramBasedSegmentation(){
        //Do nothing
        }

        public static ImagePlus[] computeAndDisplayJointHistograms(ImagePlus imgR,ImagePlus imgG){
            ImagePlus imgProba=ij.gui.NewImage.createImage("Proba",256,256,1,32,ij.gui.NewImage.FILL_BLACK);
            ImagePlus imgLogProba=ij.gui.NewImage.createImage("Log-Proba",256,256,1,32,ij.gui.NewImage.FILL_BLACK);
            ImagePlus imgMutualProba=ij.gui.NewImage.createImage("Mutual Probability", 256, 256, 1, 32, ij.gui.NewImage.FILL_BLACK);
            float[] dataProba=(float[])(imgProba.getStack().getProcessor(1).getPixels());
            float[] dataLogProba=(float[])(imgLogProba.getStack().getProcessor(1).getPixels());
            float[] dataMutualProba=(float[])(imgMutualProba.getStack().getProcessor(1).getPixels());


            //Normally we should verify that they have the same size
            int X=imgR.getWidth();
            int Y=imgR.getHeight();
            int Z=imgR.getStackSize();
            int Npix=X*Y*Z;
            double fuzzyValue=Math.log(1/Npix)-5;

            //We know that values are between 0 and 255
            int[][]jointHistogramSumOfPixels=new int[256][256];
            double[][]jointHistogramProbability=new double[256][256];
            double[][]jointHistogramLogProbability=new double[256][256];
            double[] marginalDistributionX = new double[256];
            double[] marginalDistributionY = new double[256];
            double[][]mutualProba = new double[256][256];

            //Gain access to pixels of image
            byte[] dataR;
            byte[] dataG;
            for(int z=0;z<Z;z++) {
                dataR=(byte[])(imgR.getStack().getProcessor(z+1).getPixels());
                dataG=(byte[])(imgG.getStack().getProcessor(z+1).getPixels());
                for(int x=0;x<X;x++) {
                    for(int y=0;y<Y;y++){
                        int valR=((int)(dataR[x+X*y] & 0xff));
                        int valG=((int)(dataG[x+X*y] & 0xff));
                        jointHistogramSumOfPixels[valR][valG]+=1;
                    }
                }
            }
            //Normalize for getting a probability
        
            for(int r=0;r<256;r++)for(int g=0;g<256;g++){
                jointHistogramProbability[r][g]=jointHistogramSumOfPixels[r][g]*1.0/Npix;

                //write things in the proba image
                dataProba[r+256*(255-g)]=(float)jointHistogramProbability[r][g];

                if(jointHistogramSumOfPixels[r][g]!=0){
                    jointHistogramLogProbability[r][g]=Math.log(jointHistogramSumOfPixels[r][g]*1.0/Npix);
                    // Calculate marginal distribution along X-axis
                    marginalDistributionX[r] += jointHistogramProbability[r][g];
                    
                    // Calculate marginal distribution along Y-axis
                    marginalDistributionY[g] += jointHistogramProbability[r][g];

                    // Calculate mutual information
                    if (jointHistogramSumOfPixels[r][g] != 0) {
                        double pxy = jointHistogramProbability[r][g];
                        double px = marginalDistributionX[r];
                        double py = marginalDistributionY[g];
                        
                        mutualProba[r][g] = pxy * Math.log(pxy / (px * py));
                    }
                }
                    else{
                    jointHistogramLogProbability[r][g]=fuzzyValue;
                    //mutualProba[r][g] =fuzzyValue;
                }

    
               

                //write things in the logproba image
                dataLogProba[r+256*(255-g)]=(float)jointHistogramLogProbability[r][g];
                

                //write things in mutual image
                // dataMutualProba[r+256*(255-g)]=(float)(jointHistogramProbability[r][g] *(jointHistogramLogProbability[r][g])/(marginalDistributionX[r] * marginalDistributionY[g]));
                dataMutualProba[r+256*(255-g)]=(float)(mutualProba[r][g]);

            }
                imgProba.show();
                imgLogProba.show();
                imgMutualProba.show();
                IJ.run(imgLogProba,"Fire","");
                IJ.run(imgMutualProba,"Fire","");
                imgLogProba.setDisplayRange(-18, 4);
                imgMutualProba.setDisplayRange(10e-14, 10e-5);
                return new ImagePlus[]{imgProba,imgLogProba,imgMutualProba};
            }
        
}
