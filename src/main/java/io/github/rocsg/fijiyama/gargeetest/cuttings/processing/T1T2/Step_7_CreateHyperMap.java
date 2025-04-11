package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.T1T2;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.Concatenator;
import ij.plugin.HyperStackConverter;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Config;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;

public class Step_7_CreateHyperMap implements PipelineStep{
    @Override
    public void execute(Specimen specimen) throws Exception {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'execute'");

    }
  
    public static void main(String[] args) throws Exception {
        ImageJ ij = new ImageJ();           
        Specimen spec= new Specimen("B_223");
        new Step_7_CreateHyperMap ();
        Step_7_CreateHyperMap.execute(spec,true); 
    }

    /**
     * Create a hypermap for a given specimen by concatinating registered and normalized hyperstacks and hyperframes.
     * 
     * @param specimen the specimen to process
     * @param testing whether this is a test run
     * @return the hypermap
     * @throws Exception
     */
    public static ImagePlus execute(Specimen specimen, boolean testing) throws Exception {
        ImageJ ij = new ImageJ();
        String[] timestamp = Config.timestamps;
        int N = timestamp.length;
        ImagePlus[] transformedImages = new ImagePlus[N];
        for(int t = 0; t<N; t++){
            // Load the image at the given timepoint
            ImagePlus img = IJ.openImage("/home/phukon/Desktop/206_T1_T2/206_"+timestamp[t]+"_hyperMap.tif");
            transformedImages[t] = img;
            System.out.println("Loaded and transformed: " + timestamp[t]);
        }
        // Create a hyperMap from the transformed images
        Concatenator con=new Concatenator();
        con.setIm5D(true);
        ImagePlus hypTemp=con.concatenate(transformedImages,true);
        ImagePlus hyperMap = HyperStackConverter.toHyperStack(hypTemp, 19, 40, 4, "Grayscale");
        // Save the hyperstack
        hyperMap.show();
        IJ.saveAsTiff(hyperMap,"/mnt/41d6c007-0c9e-41e2-b2eb-8d9c032e9e53/gargee/T1_T2_Specimen_Analysis/results/"+specimen.getName()+"/"+specimen.getName()+"_hypermap_registered.tif" );
        return hyperMap;
    }

   
    
}
