package barna.astafunk.utils;

import barna.commons.parameters.*;
import java.io.File;

/**
 * @version 2
 * @autor vitorlc
 * @since 11/09/14
 */

/**
 * This class describe ASTAFUNK settings. FunkSettingsOld extends ParameterSchema to
 * register parameters in JSAP
 */
public class FunkSettings extends ParameterSchema {

    /**
     * Very high value that represents "infinity" compared to the scores used in the model.
     */
    public static final double INF = 9999999;

//    /**
//     * Print all hits founded on entire sequence.
//     */
//    public static final Parameter<Boolean> ALLHITS = Parameters.booleanParameter("ALLHITS",
//            "Complete output of ASTAFUNK. Include the whole sequence " +
//                    "and alternative splicing region hits (default: only AS region hits output)",
//            false, null)
//            .longOption("all").shortOption('a');

    /**
     * Number of threads to run ASTAFUNK.
     */
    public static final Parameter<Integer> CPU = Parameters.intParameter("CPU",
            "Number of threads to run (Default: 1)",
            1, null)
            .longOption("cpu");

    /**
     * Performs an exhaustive search against HMM database (without heuristic table).
     */
    public static final Parameter<Boolean> EXHAUSTIVE =
            Parameters.booleanParameter("EXHAUSTIVE",
            "Perform exhaustive search against HMM database (default: heuristic search)", false, null)
            .longOption("exh").shortOption('e');


    /**
     * Path to the directory with the genomic sequences,
     * i.e., one fasta file per chromosome/scaffold/contig
     * with a file name corresponding to the identifiers of
     * the first column in the GTF annotation.
     */
    public static final Parameter<File> GENOME = Parameters.fileParameter("GENOME",
            "Path to the directory with the genomic sequences,\n" +
                    "i.e., one fasta file per chromosome/scaffold/contig\n" +
                    "with a file name corresponding to the identifiers of\n" +
                    "the first column in the GTF annotation", null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                    File genomeFile = (File) schema.get(parameter);
                    if (genomeFile == null) {
                        throw new ParameterException("You have to specify a genome directory");
                    }
                }
            }, null)
            .longOption("genome");


    /**
     * Annotation file in Gene Transfer Format (GTF)
     */
    public static final Parameter<File> GTF = Parameters.fileParameter("GTF",
            "Path to the GTF reference annotation",
            null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                    File refFile = (File) schema.get(parameter);
                    if (refFile == null) {
                        throw new ParameterException("Reference annotation cannot be null!");
                    }
                    if (!refFile.exists()) {
                        throw new ParameterException("The reference annotation " + refFile.getAbsolutePath() + " could not be found!");
                    }
                }
            }, null)
            .longOption("gtf");

    /**
     * Path to the profile HMM file.
     */
    public static final Parameter<File> HMM_FILE = Parameters.fileParameter("HMM_FILE",
            "Path to the profile HMM file",
            null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter)
                        throws ParameterException {
                    File refFile = (File) schema.get(parameter);
                    if (refFile == null) {
                        throw new ParameterException("Profile HMM file cannot be null!");
                    }

                    if (!refFile.exists()) {
                        throw new ParameterException("The profile HMM file " +
                                refFile.getAbsolutePath() + " could not be found!");
                    }
                }
            })
            .longOption("hmm");

    /**
     * Local search.
     */
    public static final Parameter<Boolean> LOCAL = Parameters.booleanParameter("LOCAL",
            "Run local search (Default: glocal", false, null)
            .longOption("local").shortOption('l');

    /**
     * If true performs search on all transcripts of a gene with AS.
     */

    public static final Parameter<Boolean> NAIVE_SEARCH = Parameters.booleanParameter("NAIVE_SEARCH",
            "Run Näive search. Search domains against all genes with alternaive splicig. This search uses a reference" +
                    "domain file (Method to obtain results for the paper).", false, null)
            .longOption("naive");

    /**
     * Overlapping hit ratio allowed.
     */
    public static final Parameter<Double> OVERLAPPING = Parameters.doubleParameter("OVERLAPPING",
            "Hit overlapping allowed (default: 0)",
            0, 0, 1, null)
            .shortOption('o');

    /**
     * Path to the reference domain annotation.
     */
    public static final Parameter<File> REFERENCE_FILE = Parameters.fileParameter("REFERENCE_FILE",
            "Path to the reference domain file",
            null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter)
                        throws ParameterException {
                    File refFile = (File) schema.get(parameter);
                    if (refFile == null) {
                        if (!schema.get(EXHAUSTIVE)) {
                            throw new ParameterException("Reference domain file cannot be null " +
                                    "on heuristic search!");
                        }
                    } else if (!refFile.exists()) {
                        throw new ParameterException("The reference domain file "
                                + refFile.getAbsolutePath() + " could not be found!");
                    }
                }
            }, null)
            .longOption("reference").shortOption('r');

    /**
     * Additional tool: print reference protein on FASTA format.
     */
    public static final Parameter<Boolean> REFTRANSCRIPT = Parameters
            .booleanParameter("REFTRANSCRIPT",
                    "Output FASTA file with reference transcript of each gene. " +
                            "This parameter is only used with" +
                            "--gtf and --genome parameters", false, null)
            .longOption("tref");

    public static final Parameter<Boolean> CONSTITUTIVE = Parameters
            .booleanParameter("CONSTITUTIVE",
                    "Performs a domain search only on constituve regions of all genes " +
                            "(Method to obtain results for the paper)", false, null)
            .longOption("const");

    /**
     * If true performs search on all transcripts of a gene with AS.
     */

    public static final Parameter<Boolean> OUTPUT_HITS_PER_GENE = Parameters.booleanParameter("OUTPUT_HITS_PER_GENE",
            "Output best non-overlapped domain hits of the alternatively spliced (AS) gene." +
                    "(Default: output best non-overlapped domain hits of each variant).", false, null)
            .shortOption('g');

    public static final Parameter<Boolean> OUTPUT_ALL_TRANSCRIPT_HITS = Parameters.booleanParameter("OUTPUT_ALL_TRANSCRIPT_HITS",
            "Output the all domain hits of each alternative variant."+
                    "(Default: output best non-overlapped domain hits of each variant).", false, null)
            .longOption("all");

    /**
     * Verbose mode.
     */
    public static final Parameter<Boolean> VERBOSE = Parameters.booleanParameter("VERBOSE",
            "Verbose", false, null)
            .longOption("verbose");


}