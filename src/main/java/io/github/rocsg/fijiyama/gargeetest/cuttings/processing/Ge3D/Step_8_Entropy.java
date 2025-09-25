package io.github.rocsg.fijiyama.gargeetest.cuttings.processing.Ge3D;
import ij.ImageJ;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.PipelineStep;
import io.github.rocsg.fijiyama.gargeetest.cuttings.core.Specimen;
public class Step_8_Entropy implements PipelineStep {
     @Override
    public void execute(Specimen specimen) throws Exception {
        execute(specimen, true);
    }

    public static void main (String[] args) throws Exception {
        Specimen spec = new Specimen("B_201");
        
        ImageJ ij = new ImageJ();
        new  Step_8_Entropy().execute(spec,true);
    }

    public void execute (Specimen specimen, boolean testing) throws Exception {
        

    }

   
}
