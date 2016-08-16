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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.diffexp;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Flux Tool implementation to expose differential expression analysis
 * to the command line
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class DiffExpTool implements Tool<Void> {
    private static final String longDescription = "Calculate differential expression between two quantificatons.\n" +
            "\n" +
            "Both quantifications must be done with the capacitor on the same\n" +
            "annotation.\n" +
            "A Fisher Exact Test is performed and the tools prints a tab separated\n" +
            "table with the at least the following entries in order\n" +
            "\n" +
            "id           : The transcript id\n" +
            "               (note that this might not be unique (RefSeq), but\n" +
            "                you can combine it with the location)\n" +
            "location     : Location as <chr>:<start>:<end>:<strand>\n" +
            "S1_reads     : #Reads Sample 1\n" +
            "S1_rpkm      : RPKM Sample 1\n" +
            "S2_reads     : #Reads Sample 2\n" +
            "S2_rpkm      : RPKM Sample 2\n" +
            "difference   : Absolute RPKM Difference (S2_rpkm - S1_rpkm)\n" +
            "fold_change  : log_2 fold change\n" +
            "p-value      : Fisher Exact P-Value\n" +
            "p-fdr        : P-Value corrected with FDR correction (Benjamini-Hochberg)\n" +
            "p-bonferroni : P-value corrected with Bonferroni ( 2 Samples by default, see options)\n" +
            "                The tool works in a pairwise fashion on 2 samples, but if you\n" +
            "                know already that you have X samples, you can specify --samples 12\n" +
            "                to compute Bonferroni correction for 12 samples.\n" +
            "\n" +
            "In addition, you can specify a model to annotate the entries with further information.\n" +
            "The model can be either a GTF file (i.e. Gencode annotation) or a tab separated\n" +
            "table where the first line contains the attribute names, first column is the ID\n" +
            "and second column is the type (currently only 'transcript' is supported). All additional\n" +
            "columns can contains arbitrary information that you can pick using the model-attributes\n" +
            "option.\n" +
            "\n" +
            "NOTE: The tool supports multiple threads, activate multi-threading with the --threads option\n";

    private File sample2;
    private File sample1;
    private File model;
    private String[] attributes;
    private int samples = 2;
    private Writer output;

    @Override
    public String getName() {
        return "diffexp";
    }

    @Override
    public String getDescription() {
        return "Differential Expression between two samples";
    }

    @Override
    public String getLongDescription(){
        return longDescription;
    }

    @Override
    public List<Parameter> getParameter() {
        List<Parameter> params = new ArrayList<Parameter>();
        params.add(JSAPParameters.flaggedParameter("sample-1", '1')
                .valueName("gtf")
                .help("Sample 1 Quantification")
                .required()
                .get());
        params.add(JSAPParameters.flaggedParameter("sample-2", '2')
                .valueName("gtf")
                .help("Sample 2 Quantification")
                .required()
                .get());
        params.add(JSAPParameters.flaggedParameter("output", 'o')
                .valueName("output")
                .help("Target file, defaults to stdout")
                .get());
        params.add(JSAPParameters.flaggedParameter("model", 'm')
                .valueName("model")
                .help("Model file (GTF or tab separated table)")
                .get());
        params.add(JSAPParameters.flaggedParameter("model-attributes", 'a')
                .valueName("<atr>[,<atr>]*")
                .help("Comma separated list of attributes that are taken from the model file")
                .get());
        params.add(JSAPParameters.flaggedParameter("samples", 's')
                .defaultValue("2")
                .valueName("#samples")
                .help("Number of overall samples used for Bonferroni Correction")
                .get());
        return params;
    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        if(!args.userSpecified("sample-1") || !args.userSpecified("sample-2")){
            Log.error("You have to specify two sample files !");
            return false;
        }
        if(args.userSpecified("model-attributes") && !args.userSpecified("model")){
            Log.error("You have to specify a model file if you are using additional attributes");
            return false;
        }

        sample1 = new File(args.getString("sample-1"));
        sample2 = new File(args.getString("sample-2"));
        if(!sample1.exists()){
            Log.error("Sample file " + sample1 + " not found!");
            return false;
        }
        if(!sample2.exists()){
            Log.error("Sample file " + sample2 + " not found!");
            return false;
        }

        if(args.userSpecified("model")){
            model = new File(args.getString("model"));
            if(!model.exists()){
                Log.error("Model file " + model + " not found!");
                return false;
            }
        }

        if(args.userSpecified("output")){
            try {
                output = new FileWriter(new File(args.getString("output")));
            } catch (IOException e) {
                Log.error("Error while preparing output file : " + e.getMessage());
                return false;
            }

        }else{
            // default to std out
            output = new PrintWriter(System.out);
        }

        if(args.userSpecified("model-attributes")){
            attributes = args.getString("model-attributes").split(",");
        }
        samples = 2;
        if(args.userSpecified("samples")){
            try {
                samples = Integer.parseInt(args.getString("samples"));
                if(samples < 2){
                    Log.error("Number of samples has to be >= 2!");
                    return false;
                }
            } catch (NumberFormatException e) {
                Log.error("Unable to parse number of samples : " + args.getString("samples"));
            }
        }
        return true;
    }

    @Override
    public Void call() throws Exception {
        FisherDifferentialExpression ff = new FisherDifferentialExpression();
        List<DifferentialExpression> differentialExpressions = ff.computeDifferentialExpression(sample1, sample2, samples);
        ff.write(output, differentialExpressions, model, attributes);
        return null;
    }
}
