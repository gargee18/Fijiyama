 //Input :  point 1 is top cutting, point 2 is bottom cutting, point 3 is inoculation point (gathered with ImageJ), given in RealWorldCoordinates (in mm for example)
 //Output : the ItkTransform that can be applied to an image in order to get it resampled, with the inoculation point at the given coordinates (in RWC)

package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.registration.ItkTransform;
import io.github.rocsg.fijiyama.registration.TransformUtils;
import ij.IJ;
import ij.ImagePlus;


public class Step_1_InocAlignment implements PipelineStep{

    public static void main(String[] args) throws Exception{
        Specimen spec= new Specimen("B_201");
        new Step_1_InocAlignment().execute(spec,true); 
    }

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }

    public void execute(Specimen specimen,boolean testing) throws Exception {
        String[] timestamps = Config.timestamps;
        int N = timestamps.length;

        for(int n=0; n< (testing ? 1 : N); n++){
            String filePath = Config.mainDir+"TargetCoordsForAlignment.csv";

            // Read the CSV file into a 2D array (String[][])
            String[][] data = VitimageUtils.readStringTabFromCsv(filePath);
            
            //Step 1 : gather the coordinates and the image
            double[]Pinoc=getInocCoord(data,specimen);
            double[]Ptop=getTopCoord(data,specimen);            
            double[]Pbot=getBotCoord(data,specimen);
            double[]PtargetPinoc=getTargetCoord();

            //Step 2 : call the computebasis function and gather the ItkTransform
            ItkTransform trAlignMat = computeBasisChangeAccordingToThreeGivenPoints(Ptop,Pbot,Pinoc,PtargetPinoc);
            trAlignMat.writeMatrixTransformToFile(Config.getPathToInoculationAlignmentTransformation(specimen,n));
            System.out.println(trAlignMat);

            //Step 3 : apply to the image and display it
            ImagePlus imgTest=IJ.openImage(Config.getPathToNormalizedImage(specimen,n));
            ImagePlus imgAligned=trAlignMat.transformImage(imgTest,imgTest);
            imgAligned.setTitle(specimen + "_aligned");
            IJ.saveAsTiff(imgAligned, Config.getPathToInocAlignedImage(specimen, n));
            if(testing)imgAligned.show();

        }
        
    }
    
   
    public static ItkTransform computeBasisChangeAccordingToThreeGivenPoints(double[]Ptop,double[]Pbot,double[]Pinoc,double[]PtargetPinoc){
        //Step 1 : compute vectors
        double[] Vz = TransformUtils.normalize(TransformUtils.vectorialSubstraction(Ptop, Pbot)); //u
        double[] Vv = TransformUtils.vectorialSubstraction(Pinoc, Pbot); //v
        double [] VBotToProj=TransformUtils.multiplyVector(Vz,TransformUtils.scalarProduct(Vv, Vz)); // projection of vector bottom to inoculation point
        double[] Pproj = TransformUtils.vectorialAddition(Pbot, VBotToProj);
        
        //Step 2 : compute the basis (the normalized vectors of the orientation of the object)
        double[] Vy = TransformUtils.normalize(TransformUtils.vectorialSubstraction(Pinoc, Pproj));
        double[] Vx = TransformUtils.normalize(TransformUtils.vectorialProduct(Vy, Vz));

        //Step 3 : Compute the translation
        double[] tr = computeTranslation(Pinoc, PtargetPinoc, Vx, Vy, Vz) ;

        //Step 4 : Create an ItkTransform from it
        ItkTransform trAlignMat=ItkTransform.array16ElementsToItkTransform(new double[]{Vx[0],Vy[0],Vz[0],tr[0],   Vx[1],Vy[1],Vz[1],tr[1],   Vx[2],Vy[2],Vz[2],tr[2],   0,0,0,1 });
        return trAlignMat;
    }

  

    // pTar : the target coordinates where the inoculation point is intended to be positioned
    public static double[] getTargetCoord(){
        double[] pTarPix = new double[]{256,384,512};
        double voxSize = 0.0351563;
        double[] pTar = new double[pTarPix.length];
        for (int i = 0; i < pTar.length; i++) {
            pTar[i] = pTarPix[i] * voxSize;
        }
        return pTar;
    } 

   // Tx Ty Tz
    public static double[] computeTranslation(double[] pInoc, double[] pTar, double[] vecX, double[] vecY, double[] vecZ){
        double X = pInoc[0], Y = pInoc[1], Z = pInoc[2];
        double Xp = pTar[0], Yp = pTar[1], Zp = pTar[2];
        double tx = X - (vecX[0] * Xp + vecY[0] * Yp + vecZ[0] * Zp);
        double ty = Y - (vecX[1] * Xp + vecY[1] * Yp + vecZ[1] * Zp);
        double tz = Z - (vecX[2] * Xp + vecY[2] * Yp + vecZ[2] * Zp);
        return new double[]{tx, ty, tz};
    }

    public static double[] extractCoordinates(String[][] data, String rowName) {
        for (String[] row : data) {
            if (row[1].equals(rowName)) {  
                double x = Double.parseDouble(row[17]); // X(mm)_sub is at index 17
                double y = Double.parseDouble(row[18]); // Y(mm)_sub is at index 18
                double z = Double.parseDouble(row[19]); // Z(mm)_sub is at index 19
                return new double[]{x, y, z}; 
            }
        }
        return null;
    }

    public static double[] getTopCoord(String[][] data, Specimen specimen) {
        return extractCoordinates(data, specimen.getName() + "_J001_sub222_top");
    }

    public static double[] getBotCoord(String[][] data, Specimen specimen) {
        return extractCoordinates(data, specimen.getName() + "_J001_sub222_bot");
    }

    public static double[] getInocCoord(String[][] data, Specimen specimen) {
            return extractCoordinates(data, specimen.getName() + "_J001_sub222_inoc");
    }
 

}





