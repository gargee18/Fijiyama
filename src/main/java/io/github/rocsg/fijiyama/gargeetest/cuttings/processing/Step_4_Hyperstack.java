package io.github.rocsg.fijiyama.gargeetest.cuttings.processing;
import java.io.File;

//specific libraries
import ij.IJ;
import ij.ImagePlus;
//my libraries
import io.github.rocsg.fijiyama.common.VitimageUtils;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_4_Hyperstack implements PipelineStep{
    public void execute(Specimen specimen, boolean testing) throws Exception {
        Specimen spec= new Specimen("B_201");
        new Step_4_Hyperstack().execute(spec,true); 
        ImagePlus hyperStack=null;
        String dir="";
        //Test if files are available to reach step 4
        if(false){
            //If so, compute hyperstack to level 4

        }

        else{
            //Test if files are available to reach step 3
            if(false){
                //If so, compute hyperstack to level 3

            }
    
            else{
                //Test if files are available to reach step 2
                if(new File(dir+"Transform_3raw_to_2raw.txt").exists()){
                    //If so, compute hyperstack to level 2

                }
        
                else{
                    //Test if files are available to reach step 1 ()
                    if(new File(dir+"Transform_3raw_to_2raw.txt").exists()){
                        //If so, compute hyperstack to level 1
    
                    }
            
                    else{
                        //compute level 0 hyperstack (all the raws)
                        createHyperStacks(spec);
                    }
                        
                }
                
            }
    
        }

       //createHyperStacks(specimen);
    }


    // public static void main(String[] args) {

    //     ImageJ ij = new ImageJ();
    //     computeAllHyperStacks();
    // }


    public static ImagePlus createHyperStacks(Specimen specimen) {
        ImagePlus imgB001 = IJ.openImage(Config.mainDir + specimen + "/raw/" + specimen + "_J001_aligned.tif");
        ImagePlus imgB029 = IJ.openImage(Config.mainDir + specimen + "/outputs/" + specimen + "_001_029_TransformedImage_Reg_Aligned.tif");
        ImagePlus imgB077 = IJ.openImage(Config.mainDir + specimen + "/outputs/" + specimen + "_001_077_TransformedImage_Reg_Aligned.tif");
        ImagePlus imgB141 = IJ.openImage(Config.mainDir + specimen + "/outputs/" + specimen + "_001_0141_TransformedImage_Reg_Aligned.tif");
        ImagePlus[] tabImgover4setsofdays = new ImagePlus[] { imgB001, imgB029, imgB077, imgB141 };
        ImagePlus hyperStackXR = VitimageUtils.hyperStackingChannels(tabImgover4setsofdays);
        ImagePlus hyperFrame = VitimageUtils.hyperStackChannelToHyperStackFrame(hyperStackXR);
        // hyperFrame.setTitle(specimen+"_Hyperstack");
        // hyperFrame.show();
        IJ.saveAsTiff(hyperFrame, Config.mainDir + specimen + "/hyperimage/" + specimen + "_Hyperstack.tif");
        return hyperFrame;
    }

    @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen,false);
    }

}

