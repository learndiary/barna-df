/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.diffexp;

import barna.commons.Execute;
import barna.commons.log.Log;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

/**
 * Do fisher exact test to get p-values for differential expression between
 * two quantifications created with the capacitor
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
public class FisherDifferentialExpression {
    private QuantificationModel model;
    private File modelFile;

    /**
     * Create differential expression without a model
     */
    public FisherDifferentialExpression() {
        this(null);
    }

    /**
     * Create a differential expression with a gene model. The model
     * can be used to annotate entries from the quantification.
     *
     * @param model the GTF model file
     */
    public FisherDifferentialExpression(File model) {
        this.modelFile = model;
    }

    /**
     * Compute the differential expression that compares source against target
     *
     * @param sourceFile the source GTF file
     * @param targetFile the target GTF file
     * @return diffExp list of differential expression values sorted by P-Value
     * @throws Exception in case an error occurs
     */
    List<DifferentialExpression> computeDifferentialExpression(File sourceFile, File targetFile) throws Exception {
        return computeDifferentialExpression(sourceFile, targetFile, 2);
    }
    /**
     * Compute the differential expression that compares source against target
     *
     * @param sourceFile the source GTF file
     * @param targetFile the target GTF file
     * @param samples number of samples used for the bonferroni correction
     * @return diffExp list of differential expression values sorted by P-Value
     * @throws Exception in case an error occurs
     */
    List<DifferentialExpression> computeDifferentialExpression(File sourceFile, File targetFile, int samples) throws Exception {
        if(sourceFile == null || targetFile == null) throw new NullPointerException();
        Log.info("Reading input");
        Log.progressStart("Reading source");
        Quantification source = Quantification.read(sourceFile);
        Log.progressFinish("Done", true);

        Log.progressStart("Reading target");
        Quantification target = Quantification.read(targetFile);
        Log.progressFinish("Done", true);



        double sourceReads = source.getTotalReads();
        double targetReads = target.getTotalReads();
        Log.info("Computing Differential Expression Transcripts");
        Log.progressStart("Computing");

        List<Future<DifferentialExpression>> jobs = new ArrayList<Future<DifferentialExpression>>();
        List<DifferentialExpression> exp = new ArrayList<DifferentialExpression>();

        // pump source -> target comparisons
        for (QuantificationEntry sourceTranscript : source.transcripts()) {
            QuantificationEntry targetTranscript = target.getTranscript(sourceTranscript.getKey());
            QuantificationEntryComparator job = new QuantificationEntryComparator(sourceTranscript, targetTranscript, sourceReads, targetReads);
            jobs.add(Execute.getExecutor().submit(job));
        }
        // pump target->null jobs (in target but not in source
        for (QuantificationEntry targetTranscript : target.transcripts()) {
            QuantificationEntry sourceTranscript = source.getTranscript(targetTranscript.getKey());
            if(sourceTranscript == null){
                QuantificationEntryComparator job = new QuantificationEntryComparator(sourceTranscript, targetTranscript, sourceReads, targetReads);
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

        Log.info("Correcting P-Values");
        Corrections.fdr(exp);
        Corrections.bonferroni(exp, samples);

        return exp;
    }

    /**
     * Write the expression to the writer
     *
     * @param writer the target writer
     * @param expressions the expressions
     */
    public void write(Writer writer, List<DifferentialExpression> expressions) throws IOException{
        if(writer == null) throw new NullPointerException("NULL writer not permitted");
        if(expressions == null) throw new NullPointerException("NULL expression not permitted");
        BufferedWriter b = new BufferedWriter(writer);
        for (DifferentialExpression expression : expressions) {
            StringBuilder s = new StringBuilder();
            s.append(expression.getSource().getId()).append("\t");
            s.append(expression.getSource().getName()).append("\t");
            s.append(expression.getSource().getReadCount()).append("\t");
            s.append(expression.getSource().getRpkm()).append("\t");
            s.append(expression.getTarget().getReadCount()).append("\t");
            s.append(expression.getTarget().getRpkm()).append("\t");
            s.append(expression.getDifference()).append("\t");
            double foldChange = expression.getFoldChange();
            s.append(foldChange == Double.NEGATIVE_INFINITY ? "-" : foldChange == Double.POSITIVE_INFINITY ? "+": foldChange).append("\t");
            s.append(expression.getP()).append("\t");
            s.append(expression.getFdrP()).append("\t");
            s.append(expression.getBonferroniP()).append("\t");
            s.append("\n");
            b.write(s.toString());
        }
    }

    /**
     * Get the gene model if the model file was specified.
     *
     * @return model the gene model
     */
    QuantificationModel getModel() {
        if(modelFile != null){
            if(model == null){
                Log.info("Reading model");
                Log.progressStart("Reading quantification model");
                try {
                    model = QuantificationModel.read(modelFile);
                } catch (IOException e) {
                    Log.progressFailed("Error while reading model from " + modelFile.getAbsolutePath() + ": " + e.getMessage());
                    throw new RuntimeException("Error while reading model from " + modelFile.getAbsolutePath() + ": " + e.getMessage(), e);
                }
                Log.progressFinish("Done", true);
            }
        }
        return null;
    }

    public static void main(String[] args) throws Exception {
        String model = "/Users/thasso/data/annotations/hg/gencode_v12.gtf";
//        String file_1 = "/Users/thasso/data/quantifications/data_1.bam_gencode_v12.gtf.gtf";
//        String file_2 = "/Users/thasso/data/quantifications/data_2.bam_gencode_v12.gtf.gtf";
        String file_1 = "/Users/thasso/data/quantifications/data_1.bam_hg19_ref_ucsc120203.gtf.gz.gtf";
        String file_2 = "/Users/thasso/data/quantifications/data_2.bam_hg19_ref_ucsc120203.gtf.gz.gtf";
        String output = "/Users/thasso/data/quantifications/diff.txt";
        FisherDifferentialExpression ff = new FisherDifferentialExpression(new File(model));
        Execute.initialize(4);

        List<DifferentialExpression> expressions = ff.computeDifferentialExpression(new File(file_1), new File(file_2));
        ff.write(new FileWriter(output), expressions);

        Execute.shutdown();
    }



}
