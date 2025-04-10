package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;

//specific libraries
import ij.IJ;
import ij.ImagePlus;
//my libraries
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
import io.github.rocsg.fijiyama.registration.ItkTransform;

public class Step_4_Hyperstack implements PipelineStep{

    public static void main(String[] args) throws Exception{
        Specimen spec= new Specimen("B_239");
        new Step_4_Hyperstack().execute(spec,true); 
    }
 

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }
    
    public void execute(Specimen specimen, boolean testing) throws Exception {
        String[] days=Config.timestamps;

        ImagePlus[] tabImgover4setsofdays = new ImagePlus[4];
        int N=days.length;

        for (int i = 0; i<N; i++){

            // raw
            ImagePlus img= IJ.openImage(Config.getPathToNormalizedImage(specimen,i));
            ImagePlus imgTransformed = img;

            // open tr and compose
            if(i==0){
                ItkTransform trInoc = ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen, i));
                imgTransformed =trInoc.transformImage(img,imgTransformed);

            }
           
            else {
                // Apply the inoculation alignment first
                ItkTransform trInocMov = ItkTransform.readTransformFromFile(Config.getPathToInoculationAlignmentTransformation(specimen, i));
                imgTransformed = trInocMov.transformImage(img, imgTransformed);
        
                // Apply all previous rigid transformations in order
                for (int j = 0; j < i; j++) {
                    ItkTransform trRigid = ItkTransform.readTransformFromFile(Config.getPathToAutoRigidRegistrationMatrix(specimen, j, j + 1));
                    imgTransformed = trRigid.transformImage(img, imgTransformed);
                }
            }



            // save to tab
            tabImgover4setsofdays[i] = imgTransformed;
            System.out.println(tabImgover4setsofdays[i]);
            }
        ImagePlus hyperStackXR = VitimageUtils.hyperStackingChannels(tabImgover4setsofdays);
        ImagePlus hyperFrame = VitimageUtils.hyperStackChannelToHyperStackFrame(hyperStackXR);
        hyperFrame.setTitle(specimen+"_Hyperstack");
        hyperFrame.show();
        hyperFrame.setSlice(512);
        IJ.saveAsTiff(hyperFrame, Config.getPathToHyperstack(specimen));
        VitimageUtils.waitFor(5000);
        hyperFrame.close();
        }
    

}


