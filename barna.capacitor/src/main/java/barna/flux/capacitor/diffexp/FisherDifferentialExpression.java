package barna.flux.capacitor.diffexp;

import barna.commons.Execute;
import barna.commons.log.Log;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

public class FisherDifferentialExpression {
    private QuantificationModel model;
    private Quantification source;
    private Quantification target;
    private File sourceFile;
    private File targetFile;
    private File modelFile;

    public FisherDifferentialExpression(File source, File target, File model) {
        this.sourceFile = source;
        this.targetFile = target;
        this.modelFile = model;
    }

    public List<DifferentialExpression> getTranscriptDifferences() throws Exception {
        if(source == null || target == null) init();

        double sourceReads = source.getTotalReads();
        double targetReads = target.getTotalReads();
        Log.info("Computing Differential Expression Transcripts");
        Log.progressStart("Computing");

        List<Future<DifferentialExpression>> jobs = new ArrayList<Future<DifferentialExpression>>();
        List<DifferentialExpression> exp = new ArrayList<DifferentialExpression>();

        // pump source -> target comparisons
        for (Transcript sourceTranscript : source.transcripts()) {
            Transcript targetTranscript = target.getTranscript(sourceTranscript.getId());
            TranscriptComparator job = new TranscriptComparator(sourceTranscript, targetTranscript, sourceReads, targetReads);
            jobs.add(Execute.getExecutor().submit(job));
        }
        // pump target->null jobs (in target but not in source
        for (Transcript targetTranscript : target.transcripts()) {
            Transcript sourceTranscript = source.getTranscript(targetTranscript.getId());
            if(sourceTranscript == null){
                TranscriptComparator job = new TranscriptComparator(sourceTranscript, targetTranscript, sourceReads, targetReads);
                jobs.add(Execute.getExecutor().submit(job));
            }
        }

        // get the results
        long count = 0;
        for (Future<DifferentialExpression> job : jobs) {
            exp.add(job.get());
            Log.progress(count++, jobs.size());
        }
        Log.progressFinish("Done", true);
        return exp;
    }

    /**
     * Initialize and read input
     */
    void init() throws IOException {
        Log.info("Reading input");
        if(modelFile != null){
            Log.progressStart("Reading quantification model");
            model = QuantificationModel.read(modelFile);
            Log.progressFinish("Done", true);
        }
        Log.progressStart("Reading source");
        source = Quantification.read(sourceFile);
        Log.progressFinish("Done", true);

        Log.progressStart("Reading target");
        target = Quantification.read(targetFile);
        Log.progressFinish("Done", true);

    }

    public static void main(String[] args) throws Exception {
        String model = "/Users/thasso/data/annotations/gencode_v12.gtf.gz";
        String file_1 = "/Users/thasso/data/quantifications/data_1.bam_gencode_v12.gtf.gtf";
        String file_2 = "/Users/thasso/data/quantifications/data_2.bam_gencode_v12.gtf.gtf";
        FisherDifferentialExpression ff = new FisherDifferentialExpression(new File(file_1), new File(file_2), new File(model));
        Execute.initialize(4);
        List<DifferentialExpression> expressions = ff.getTranscriptDifferences();
        System.out.println("Got " + expressions.size() + " values");
        Execute.shutdown();
    }



}
