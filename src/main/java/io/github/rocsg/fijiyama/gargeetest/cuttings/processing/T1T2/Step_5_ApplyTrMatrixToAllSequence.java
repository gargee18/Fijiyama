package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.T1T2;
//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.registration.ItkTransform;
public class Step_5_ApplyTrMatrixToAllSequence implements PipelineStep{

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
    public static void main(String[] args) throws Exception{
        ImageJ ij = new ImageJ();
        Specimen spec= new Specimen("B_239");
        new Step_5_ApplyTrMatrixToAllSequence ();
        
        Step_5_ApplyTrMatrixToAllSequence.execute(spec,true); 
    }
 

    /**
     * Execute the step of applying a transformation matrix to a sequence of images (all sequences at a given time).
     * @param specimen the specimen to process
     * @param testing whether this is a test run
     * @return the transformed image
     * @throws Exception
     */
    public static void execute(Specimen specimen, boolean testing) throws Exception{
       
        String[] timestamp = Config.timestamps;
        int time = timestamp.length;
        for(int t = 0; t<time; t++){
            // open the raw image
            ImagePlus openImgSeq = IJ.openImage(Config.getRawT1T2sequence(specimen, t));
            openImgSeq.show();
            // read the transformation matrix
            ItkTransform trt = ItkTransform.readTransformFromFile(Config.getPathToSerialRegistrationTRMatrix(specimen, t));
            
            // apply the transformation matrix to the image
            ImagePlus imgTransformed = openImgSeq;
            openImgSeq.show();
            imgTransformed = trt.transformImage(openImgSeq,imgTransformed);
            imgTransformed.show();
            
            // save the transformed image
            IJ.saveAsTiff(imgTransformed, Config.getPathToRegisteredT1T2sequence(specimen, t));
            
        }
    }

    
}
