package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.T1T2;
//specific libraries
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
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
    public static ImagePlus execute(Specimen specimen, boolean testing) throws Exception{
        // open the raw image
        ImagePlus openImgSeq = IJ.openImage("/home/phukon/Desktop/206_T1_T2/206_J141_hyperMap.tif");
        openImgSeq.show();
        
        // read the transformation matrix
        ItkTransform tr001 = ItkTransform.readTransformFromFile("/home/phukon/Desktop/206_T1_T2/test/output/Exported_data/transform_global_206_J141_TR2400_TE12.txt");
        
        // apply the transformation matrix to the image
        ImagePlus imgTransformed001 = openImgSeq;
        openImgSeq.show();
        imgTransformed001 = tr001.transformImage(openImgSeq,imgTransformed001);
        imgTransformed001.show();
        
        // save the transformed image
        IJ.saveAsTiff(imgTransformed001, "/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/206_T1_T2/raw_registered/206_J141_hyperMap_reg.tif");
        return imgTransformed001;
    }

    
}
