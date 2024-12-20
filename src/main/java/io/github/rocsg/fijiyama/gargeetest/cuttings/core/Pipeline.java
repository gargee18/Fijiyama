package io.github.rocsg.fijiyama.gargeetest.cuttings.core;

import java.util.List;

public class Pipeline {
    private List<PipelineStep> steps;

    public Pipeline(List<PipelineStep> steps) {
        this.steps = steps;
    }

    public void run(Specimen specimen) throws Exception {
        for (PipelineStep step : steps) {
            step.execute(specimen);
        }
    }
}
