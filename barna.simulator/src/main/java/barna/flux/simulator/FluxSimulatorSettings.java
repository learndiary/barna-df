/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
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

package barna.flux.simulator;

import barna.commons.parameters.*;
import barna.io.FileHelper;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;

/**
 * Flux Simulator settings
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 * @author Micha Sammeth (gmicha@gmail.com)
 */
public class FluxSimulatorSettings extends ParameterSchema {

    /**
     * Instance to check the validity of the maximum memory size.
     */
    static ParameterValidator motifHeapValidator= new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {

            // check whether memory is realistic
            File f = schema.get((Parameter<File>) parameter);
            if (f== null)
                return;   // valid
            long heapMaxSize = Runtime.getRuntime().maxMemory();
            if (heapMaxSize<= HEAP_SIZE_MOTIFS) {
                DecimalFormat df = new DecimalFormat("#.##");
                throw new ParameterException("Due to the parameter value " +
                        parameter.getName() +
                        " heap size of > " +
                        df.format(HEAP_SIZE_MOTIFS/ (double) 1000000000)  +
                        " GB is required, but currently only "+ df.format(heapMaxSize/ (double) 1000000000) +
                        " have been provided. "
                );
            }

        }
    };

    /**
     * Maximum heap size required when motifs are used
     */
    protected static long HEAP_SIZE_MOTIFS= 2500000000l;

    /**
     * Loci separator
     */
    public static final char SEP_LOC_TID = '@'; // todo refactor this to profiler and use the profilers getGlobalID where possible

    /**
     * Helper to parse relative filenames
     */
    protected static RelativePathParser relativePathParser = new RelativePathParser();


    /**
     * Possible substrates
     */
    public static enum Substrate {
        DNA, RNA
    }

    /**
     * Available fragmentation methods
     */
    public static enum FragmentationMethod {
        /**
         * Nebulization fragmentation method.
         */
        NB,
        /**
         * Uniformal random fragmentation method.
         */
        UR,
        /**
         * Enzymatic digestion as fragmentation method.
         */
        EZ,
        /**
         * Disable fragmentation
         */
        NONE

    }

    public static enum RtranscriptionMode {
        /**
         * PAR_RT_MODE_POLY_DT
         */
        PDT,
        /**
         * PAR_RT_MODE_RANDOM
         */
        RH,
    }

    public static enum SizeSamplingModes {
        /**
         * rejection
         */
        RJ,
        /**
         * Acception
         */
        AC,
        /**
         * Metropolis hastings
         */
        MH
    }


    /*
    File locations
     */


    /**
     * The parameter file
     */
    private File parameterFile;

    /**
     * Path to the GTF reference annotation, either absolute
     * or relative to the location of the parameter file.
     */
    public static final Parameter<File> REF_FILE = Parameters.fileParameter("REF_FILE_NAME",
            "path to the GTF reference annotation, either absolute\n" +
            "or relative to the location of the parameter file",
            null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                    File refFile = (File) schema.get(parameter);
                    if (refFile == null) {
                        throw new ParameterException("You have to specify a reference file");
                    }
                    if (!refFile.exists()) {
                        throw new ParameterException("The reference file " + refFile.getAbsolutePath() + " could not be found!");
                    }
                }
            }, relativePathParser);
    /**
     * Path to the profile of the run, either absolute
     * or relative to the location of the parameter file;
     * the default profile uses the name of the parameter
     * file with the extension .pro.
     */
    public static final Parameter<File> PRO_FILE = Parameters.fileParameter("PRO_FILE_NAME",
            "path to the profile of the run, either absolute\n" +
            "or relative to the location of the parameter file;\n" +
            "the default profile uses the name of the parameter\n" +
            "file with the extension .pro", null, new FileValidator("pro"), relativePathParser);

    /**
     * Path to the library file of the run, either absolute
     * or relative to the location of the parameter file;
     * the default profile uses the name of the parameter
     * file with the extension .lib.
     */
    public static final Parameter<File> LIB_FILE = Parameters.fileParameter("LIB_FILE_NAME",
            "path to the library file of the run, either absolute\n" +
            "or relative to the location of the parameter file;\n" +
            "the default profile uses the name of the parameter\n" +
            "file with the extension .lib", null, new FileValidator("lib"), relativePathParser);

    /**
     * Path to the sequencing file of the run, either absolute
     * or relative to the location of the parameter file;
     * the default profile uses the name of the parameter
     * file with the extension .bed.
     */
    public static final Parameter<File> SEQ_FILE = Parameters.fileParameter("SEQ_FILE_NAME",
            "path to the sequencing file of the run, either absolute\n" +
            "or relative to the location of the parameter file;\n" +
            "the default profile uses the name of the parameter\n" +
            "file with the extension .bed", null, new FileValidator("bed"), relativePathParser);

    /**
     * Path to the directory with the genomic sequences,
     * i.e., one fasta file per chromosome/scaffold/contig
     * with a file name corresponding to the identifiers of
     * the first column in the GTF annotation.
     */
    public static final Parameter<File> GEN_DIR = Parameters.fileParameter("GEN_DIR",
            "path to the directory with the genomic sequences,\n" +
            "i.e., one fasta file per chromosome/scaffold/contig\n" +
            "with a file name corresponding to the identifiers of\n" +
            "the first column in the GTF annotation", null, new ParameterValidator() {
                @Override
                public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                    File genomeFile = (File) schema.get(parameter);
                    boolean req= requiresGenomicSequence(schema);
                    if (req&& genomeFile == null) {
                        throw new ParameterException("You have to specify a genome directory");
                    }
                    if (req&& !genomeFile.exists()) {
                        throw new ParameterException("The genome directory " + genomeFile.getAbsolutePath() + " could not be found!");
                    }

                }
            }, relativePathParser);

    /**
     * Temporary directory
     */
    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR",
            "Temporary directory", new File(System.getProperty("java.io.tmpdir")), new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
             // ISSUE IO-11 is related to this, we have to make sure the TMP dir exists and is writable
            File tmp = parameterSchema.get(TMP_DIR);
            if(tmp == null){
                throw new ParameterException("No temp directory specified!");
            }
            if(!tmp.canWrite()){
                throw new ParameterException("The temp-directory " + tmp.getAbsolutePath() + " does not exist or is not writable!");
            }

        }
    });

    /**
     * Checks by the current parameter settings whether a genomic sequence is required for the run.
     * @param schema an instance with the current parameter settings
     * @return <code>true</code> if the genomic sequence is required for the current run, <code>false</code> otherwise.
     */
    public static final boolean requiresGenomicSequence(ParameterSchema schema) {
        boolean req= (schema.get(FRAG_EZ_MOTIF)!= null)
                || (schema.get(RT_MOTIF)!= null)
                || (schema.get(GC_MEAN)!= null&& !Double.isNaN(schema.get(GC_MEAN)))
                || (schema.get(FASTA)!= null&& schema.get(FASTA)== true);
        return req;
    }



    /*
    Expression parameters
     */

    /**
     * Coding messengers, i.e., transcripts that have an annotated CDS, are extracted from the cell.
     */
    public static final Parameter<Boolean> LOAD_CODING = Parameters.booleanParameter("LOAD_CODING",
            "coding messengers, i.e., transcripts\n" +
            "that have an annotated CDS, are extracted\n" +
            "from the cell", true, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            if (!schema.get(LOAD_CODING) && !schema.get(LOAD_NONCODING)) {
                throw new ParameterException("Sorry, but either LOAD_CODING or LOAD_NONCODING has to be enabled !");
            }
        }
    });
    /**
     * Non-coding RNAs, i.e., transcripts without an annotated ORF are extracted from the cell.
     */
    public static final Parameter<Boolean> LOAD_NONCODING = Parameters.booleanParameter("LOAD_NONCODING",
            "non-coding RNAs, i.e., transcripts\n" +
            "without an annotated ORF are extracted\n" +
            "from the cell", true, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            if (!schema.get(LOAD_CODING) && !schema.get(LOAD_NONCODING)) {
                throw new ParameterException("Sorry, but either LOAD_CODING or LOAD_NONCODING has to be enabled !");
            }
        }
    });
    /**
     * Number of RNA molecules initially in the experiment.
     */
    public static final Parameter<Long> NB_MOLECULES = Parameters.longParameter("NB_MOLECULES",
            "number of RNA molecules initially in the experiment", 5000000);
    /**
     * Exponent of power-law underlying the expression profile.
     */
    public static final Parameter<Double> EXPRESSION_K = Parameters.doubleParameter("EXPRESSION_K",
            "exponent of power-law underlying the expression profile", -0.6, -1, 0, null);
    /**
     * Linear parameter of the exponential decay.
     */
    public static final Parameter<Double> EXPRESSION_X0 = Parameters.doubleParameter("EXPRESSION_X0",
            "linear parameter of the exponential decay", 9500, 1, Double.MAX_VALUE, null);
    /**
     * linear parameter of the exponential decay
     */
    public static final Parameter<Double> EXPRESSION_X1 = Parameters.doubleParameter("EXPRESSION_X1",
            "quadratic parameter of the exponential decay", Math.pow(9500, 2), 1, Double.MAX_VALUE, null);

    /**
     * Average deviation from the annotated transcription start site (TSS),
     * set to 'NaN' to deactivate simulated transcription start variability.
     */
    public static final Parameter<Double> TSS_MEAN = Parameters.doubleParameter("TSS_MEAN",
            "average deviation from the annotated transcription start site (TSS),\n" +
            "set to 'NaN' to deactivate simulated transcription start variability",
            25d, new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            // FIX #55 and make sure that NaN is supported to disable TSS_MEAN
            if (!Double.isNaN(schema.get(TSS_MEAN)) && schema.get(TSS_MEAN) <= 0) {
                throw new ParameterException("TSS_MEAN must be > 0");
            }
        }
    });
    /**
     * Shape of the Weibull distribution describing poly-A tail sizes,
     * set to 'NaN' to disable simulated poly-A tails.
     */
    public static final Parameter<Double> POLYA_SHAPE = Parameters.doubleParameter("POLYA_SHAPE",
            "shape of the Weibull distribution describing poly-A tail sizes,\n" +
            "set to 'NaN' to disable simulated poly-A tails.",
            2d, 0.0, Double.MAX_VALUE, null);
    /**
     * Scale of the Weibull distribution, shifts the average length of poly-A tail sizes,
     * set to 'NaN' to disable simulated poly-A tails.
     */
    public static final Parameter<Double> POLYA_SCALE = Parameters.doubleParameter("POLYA_SCALE",
            "scale of the Weibull distribution, shifts the average length of poly-A tail sizes,\n" +
            "set to 'NaN' to disable simulated poly-A tails.",
            300d, 0.0, Double.MAX_VALUE, null);

    /*
    Fragementation
     */
    /**
     * Turn fragmentation on/off.
     */
    public static final Parameter<Boolean> FRAGMENTATION = Parameters.booleanParameter("FRAGMENTATION",
            "turn fragmentation on/off", true);

    /**
     * Fragmentation method employed:
     * [NB] Fragmentation by nebulization
     * [UR] Uniformal random fragmentation
     * [EZ] Fragmentation by enzymatic digestion
     */
    public static final Parameter<FragmentationMethod> FRAG_METHOD = Parameters.enumParameter("FRAG_METHOD", "" +
            "Fragmentation method employed:\n" +
            "[EZ] Fragmentation by enzymatic digestion\n" +
            "[NB] Fragmentation by nebulization\n" +
            "[UR] Uniformal random fragmentation",
            FragmentationMethod.UR, new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            if (schema.get(FRAG_METHOD) == FragmentationMethod.EZ) {
                if (schema.get(FRAG_SUBSTRATE) == Substrate.RNA) {
                    throw new ParameterException("Enzymatic digestion is not supported for RNA substrate!");
                }
                if (schema.get(FRAG_EZ_MOTIF) == null) {
                    throw new ParameterException("You have to specify FRAG_EZ_MOTIF in order ot use Enzymatic digestion");
                }
            }

            if (schema.get(FRAG_METHOD) == FragmentationMethod.NB) {
                if (schema.get(FRAG_SUBSTRATE) == Substrate.RNA) {
                    throw new ParameterException("Sorry, but nebulizing RNA is not supported!");
                }
            }

        }
    });

    /**
     * Substrate of fragmentation, determines the order
     * of fragmentation and reverse transcription (RT):
     * for substrate DNA, fragmentation is carried out
     * after RT, substrate RNA triggers fragmentation
     * before RT.
     */
    public static final Parameter<Substrate> FRAG_SUBSTRATE = Parameters.enumParameter("FRAG_SUBSTRATE",
            "substrate of fragmentation, determines the order\n" +
            "of fragmentation and reverse transcription (RT):\n" +
            "for substrate DNA, fragmentation is carried out\n" +
            "after RT, substrate RNA triggers fragmentation\n" +
            "before RT", Substrate.DNA, null);


    /* Fragmentation by Enzymatic Digestikon */

    /**
     * Sequence motif caused by selective restriction
     * with an enzyme, choose pre-defined NlaIII, DpnII,
     * or a file with a custom position weight matrix.
     */
    public static final Parameter<File> FRAG_EZ_MOTIF = Parameters.fileParameter("FRAG_EZ_MOTIF",
            "sequence motif caused by selective restriction\n" +
            "with an enzyme, choose pre-defined NlaIII, DpnII,\n" +
            "or a file with a custom position weight matrix",
            null,
            motifHeapValidator,
            relativePathParser);


    /* Nebulization */

    /**
     * Threshold on molecule length that cannot
     * be broken by the shearfield of nebulization.
     */
    public static final Parameter<Double> FRAG_NB_LAMBDA = Parameters.doubleParameter("FRAG_NB_LAMBDA",
                    "Threshold on molecule length that cannot " +
                    "be broken by the shearfield of nebulization.",
                    900);
    /**
     * Threshold on the fraction of the molecule population;
     * if less molecules break per time unit, convergence
     * to steady state is assumed.
     */
    public static final Parameter<Double> FRAG_NB_THOLD = Parameters.doubleParameter("FRAG_NB_THOLD",
                    "threshold on the fraction of the molecule population;\n" +
                    "if less molecules break per time unit, convergence \n" +
                    "to steady state is assumed", 0.1, 0, 1, null);

    /**
     * Strength of the nebulization shearfield (i.e., rotor speed).
     */
    public static final Parameter<Double> FRAG_NB_M = Parameters.doubleParameter("FRAG_NB_M", "Parameter " +
            "Strength of the nebulization shearfield (i.e., rotor speed)", 1.0);

    /**
     * Average expected framgent size after fragmentations,
     * i.e., number of breaks per unit length (exhautiveness of fragmentation);
     * NaN optimizes the fragmentation process w.r.t. the size filtering.
     */
    public static final Parameter<Double> FRAG_UR_ETA = Parameters.doubleParameter("FRAG_UR_ETA",
            "Average expected framgent size after fragmentations,\n" +
            "i.e., number of breaks per unit length (exhautiveness of fragmentation);\n" +
            "NaN optimizes the fragmentation process w.r.t. the size filtering",
            Double.NaN);

    /**
     * Geometry of molecules in the UR process:
     * NaN= depends logarithmically on molecule length,
     * 1= always linear, 2= surface-diameter, 3= volume-diameter,
     * etc.
     */
    public static final Parameter<Double> FRAG_UR_DELTA = Parameters.doubleParameter("FRAG_UR_DELTA",
            "Geometry of molecules in the UR process:\n" +
            "NaN= depends logarithmically on molecule length,\n" +
            "1= always linear, 2= surface-diameter, 3= volume-diameter, etc.",
            Double.NaN);

    /**
     * Minimum length of fragments produced by UR fragmentation.
     */
    public static final Parameter<Double> FRAG_UR_D0 = Parameters.doubleParameter("FRAG_UR_D0",
            "Minimum length of fragments produced by UR fragmentation",
            1.0, 1.0, Double.MAX_VALUE, null);


    /*
     Reverse Transcription
     */


    /**
     * Switch on/off Reverse Transcription.
     */
    public static final Parameter<Boolean> RTRANSCRIPTION = Parameters.booleanParameter("RTRANSCRIPTION",
            "Switch on/off Reverse Transcription", true);

    /**
     * Primers used for first strand synthesis:
     * [RH] for random hexamers or
     * [PDT] for poly-dT primers
     */
    public static final Parameter<RtranscriptionMode> RT_PRIMER = Parameters.enumParameter("RT_PRIMER",
            "Primers used for first strand synthesis:\n" +
            "[RH] for random hexamers or\n" +
            "[PDT] for poly-dT primers", RtranscriptionMode.RH, null);

    /**
     * Minimum fragment length observed after reverse
     * transcription of full-length transcripts.
     */
    public static final Parameter<Integer> RT_MIN = Parameters.intParameter("RT_MIN",
            "Minimum fragment length observed after reverse\n" +
            "transcription of full-length transcripts.",
            500, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int min = schema.get(RT_MIN);
            int max = schema.get(RT_MAX);
            if (min > max) {
                throw new ParameterException("RT_MIN must me <= RT_MAX");
            }
        }
    });
    /**
     * Maximum fragment length observed after reverse
     * transcription of full-length transcripts.
     */
    public static final Parameter<Integer> RT_MAX = Parameters.intParameter("RT_MAX",
            "Maximum fragment length observed after reverse\n" +
            "transcription of full-length transcripts.", 5500, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int min = schema.get(RT_MIN);
            int max = schema.get(RT_MAX);
            if (min > max) {
                throw new ParameterException("RT_MIN must me <= RT_MAX");
            }
        }
    });

    /* PCR Amplification */

    /**
     * Mean value of a gaussian distribution that reflects GC bias amplification probability,
     * set this to 'NaN' to disable GC biases.
     */
    public static final Parameter<Double> GC_MEAN = Parameters.doubleParameter("GC_MEAN",
            "Mean value of a gaussian distribution that reflects GC bias amplification chance,\n" +
            "set to 'NaN' to disable GC biases.", 0.5, 0.0, 1.0, null);
    /**
     * Standard deviation of a gaussian distribution that reflects GC bias amplification probability,
     * inactive if GC_MEAN is set to NaN.
     */
    public static final Parameter<Double> GC_SD = Parameters.doubleParameter("GC_SD",
            "Standard deviation of a gaussian distribution that reflects GC bias amplification chance,\n" +
            "inactive if GC_MEAN is set to NaN.", 0.1, 0.0, 1.0, null);
    /**
     * PCR duplication probability when GC filtering is disabled by setting GC_MEAN to NaN.
     */
    public static final Parameter<Double> PCR_PROBABILITY = Parameters.doubleParameter("PCR_PROBABILITY",
            "PCR duplication probability\n" +
            "when GC filtering is disabled by setting GC_MEAN to NaN", 0.7, 0.0, 1.0, null);
    /**
     * PCR distribution file, 'default' to use
     * a distribution with 15 rounds and 20 bins,
     * 'none' to disable amplification.
     */
    public static final Parameter<String> PCR_DISTRIBUTION = Parameters.stringParameter("PCR_DISTRIBUTION", "PCR distribution file, 'default' to use .\n" +
            "a distribution with 15 rounds and 20 bins,\n" +
            "'none' to disable amplification.", "default", new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            String dist = schema.get(FluxSimulatorSettings.PCR_DISTRIBUTION);
            if(dist != null){
                // FIX SIMULATOR-2
                if(!dist.equals("default") && !dist.equals("none")){
                    if(!new File(dist).exists()) {
                        throw new ParameterException("Distribution file " + dist  + " not found");
                    }
                }
            }
        }
    });

    /**
     * Flag to force every molecule to be reversely transcribed.
     */
    public static final Parameter<Boolean> RT_LOSSLESS = Parameters.booleanParameter("RT_LOSSLESS", "Flag to force every molecule to be reversely transcribed", true);

    /**
     * Position weight matrix (PWM) used during reverse/transcription and adapter ligation.
     * This is disabled by default (value 'null'), by the value 'default' a PWM derived from
     * the current Illumina protocol is used. Optionally, a file containing a custom matrix
     * may be provided.
     */
    public static final Parameter<File> RT_MOTIF = Parameters.fileParameter(
            "RT_MOTIF",
            "Position weight matrix (PWM) used during reverse/transcription and adapter ligation.\n" +
                    "This is disabled by default (value 'null'), by the value 'default' a PWM derived from \n" +
                    "the current Illumina protocol is used. Optionally, a file containing a custom matrix \n" +
                    "may be provided.",
            null,
            motifHeapValidator,
            relativePathParser);

    /*
    Size Selection
     */
    /**
     * Switches size selection on/off.
     */
    public static final Parameter<Boolean> FILTERING = Parameters.booleanParameter("FILTERING",
            "switches size selection on/off", false);

    /**
     * Size distribution of fragments after filtering,
     * either specified by the fully qualified path of a file with an empirical distribution
     * where each line represents the length of a read, no ordering required<br>
     * </br>
     * or attributes of a gaussian distribution (mean and standard deviation) in the form:<br>
     * <br>
     * N(mean, sd)<br>
     * <br>
     * for example: N(800, 200)<br>
     * If no size distribution is provided, an empirical Illumina fragment size distribution is employed.
     */
    public static final Parameter<String> SIZE_DISTRIBUTION = Parameters.stringParameter("SIZE_DISTRIBUTION",
            "Size distribution of fragments after filtering,\n" +
            "either specified by the fully qualified path of a file with an empirical distribution" +
            "where each line represents the length of a read, no ordering required,\n" +
            "\n" +
            "or attributes of a gaussian distribution (mean and standard deviation) in the form:\n" +
            "\n" +
            "N(mean, sd) \n" +
            "\n" +
            "for example: N(800, 200)\n" +
            "If no size distribution is provided, an empirical Illumina fragment size distribution is employed."
    );

    /**
     * Method for sub-sampling fragments according to the characteristics of (see SIZE_DISTRIBUTION):<br>
     *
     * MH the Metropolis-Hastings algorithm is used for filtering<br>
     * RJ rejection sampling, employing probability directly from the distribution<br>
     * AC (acceptance) transforms the probability distribution, s.t. the 'most likely'
     *    element in the distribution has a probability of 1.0.
     */
    public static final Parameter<SizeSamplingModes> SIZE_SAMPLING = Parameters.enumParameter("SIZE_SAMPLING",
            "Method for sub-sampling fragments according to the characteristics of (see SIZE_DISTRIBUTION) \n" +
                    "\n" +
                    "MH the Metropolis-Hastings algorithm is used for filtering\n" +
                    "RJ rejection sampling, employing probability directly from the distribution\n" +
                    "AC (acceptance) transforms the probability distribution, s.t. the 'most likely'\n" +
                    "   element in the distribution has a probability of 1.0.\n", SizeSamplingModes.AC, null);

    /*
      Sequencing
     */
    /**
     * Number of reads sequenced
     */
    public static final Parameter<Long> READ_NUMBER = Parameters.longParameter("READ_NUMBER",
            "Number of reads sequenced", 5000000);
    /**
     * Length of the reads.
     */
    public static final Parameter<Integer> READ_LENGTH = Parameters.intParameter("READ_LENGTH",
            "Length of the reads", 36);
    /**
     * Switch on/off paired-end reads.
     */
    public static final Parameter<Boolean> PAIRED_END = Parameters.booleanParameter("PAIRED_END",
            "Switch on/off paired-end reads", false);
    /**
     * Creates .fasta/.fastq output.
     * Requires the genome sequences in a folder specified by GEN_DIR.
     * If a quality model is provided by parameter ERR_FILE, a .fastq
     * file is produced. Otherwise read sequences are given as .fasta.
     */
    public static final Parameter<Boolean> FASTA = Parameters.booleanParameter("FASTA",
            "creates .fasta/.fastq output.\n" +
            "Requires the genome sequences in a folder specified by GEN_DIR.\n" +
            "If a quality model is provided by parameter ERR_FILE, a .fastq\n" +
            "file is produced. Otherwise read sequences are given as .fasta.", false);

    /**
     * Path to the file with the error model.<br>
     *
     * With the values '35' or '76', default error models are provided
     * for the corresponding read lengths, otherwise the path to a
     * custom error model file is expected.
     */
    public static final Parameter<String> ERR_FILE = Parameters.stringParameter("ERR_FILE",
            "path to the file with the error model\n" +
            "\n" +
            "for the values '35' or '76', default error models are provided for the corresponding read lengths,\n" +
            "otherwise the path to a custom error model file is expected\n", null, new ParameterValidator(){
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            String v = schema.get(FluxSimulatorSettings.ERR_FILE);

            if(v != null && v.length() > 0){
                if(!v.equals("35") && !v.equals("76")){
                    // check file
                    if((v.startsWith("/") && !new File(v).canRead())){
                        throw new ParameterException("Unable to read error model from " + v + "\n" + "Use either the defaults '36' od '76' or specify an error model file");
                    }
                }
            }
        }
    });



    /**
     * Load the setting from a file
     *
     * @param f the file
     * @return settings the loaded and validated settings
     * @throws Exception in case of a parser or validation error
     */
    public static FluxSimulatorSettings createSettings(File f) throws Exception {
        if (f == null) {
            throw new NullPointerException("Null parameter file not permitted!");
        }
        if (!f.exists()) {
            throw new IllegalArgumentException("Parameter file " + f.getAbsolutePath() + " can not be found!");
        }
        InputStream in = null;
        try {

            f = new File(f.getCanonicalPath());    // kill Win32ShellFolder instances, they fuck up relative path conversion
            relativePathParser.parentDir = f.getParentFile();
            FluxSimulatorSettings settings = new FluxSimulatorSettings();
            settings.parameterFile = f;
            in = new FileInputStream(f);
            settings.parse(in);
            settings.validate();
            return settings;
        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (IOException e) {
                }
            }
        }
    }


    public int getMaxThreads() {
        return 1;
    }


    // todo : can we refactor this somehow ?
    public void setReadLength(int readLength) {
        set(READ_LENGTH, readLength);
    }


    // todo check if we can refactor this
    public void setRefFile(File refFile) {
        set(REF_FILE, refFile);
    }


    // todo check if we can refactor this
    public void setFragUReta(double fragUReta) {
        set(FRAG_UR_ETA, fragUReta);
    }

    /**
     * Get the parameter file or null
     *
     * @return parameterFile the parameter file or null
     */
    public File getParameterFile() {
        return parameterFile;
    }

    /**
     * Helper class to allow us to parse file names relative to a parent directory
     */
    static class RelativePathParser implements FileNameParser {
        /**
         * the Parent directory
         */
        File parentDir = null;

        @Override
        public File parse(String string) {
            return FileHelper.fromRelative(string, parentDir);
        }
    }


    /**
     * Validator that sets a default file name using the flux settings parameter file and a custom extension
     */
    static class FileValidator implements ParameterValidator {
        private String suffix;
        private String message;


        FileValidator(String suffix) {
            this(suffix, "You have to specify the ." + suffix + " file ");
        }

        FileValidator(String suffix, String message) {
            this.suffix = suffix;
            this.message = message;
        }

        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            if (schema.get(parameter) == null) {
                // set default pro file from parameters file
                if (schema instanceof FluxSimulatorSettings) {
                    FluxSimulatorSettings s = (FluxSimulatorSettings) schema;
                    if (s.getParameterFile() != null) {
                        File file = new File(s.getParameterFile().getParent(), s.getParameterFile().getName());
                        if (!suffix.startsWith(".")) {
                            suffix = "." + suffix;
                        }
                        file = FileHelper.replaceSfx(file, suffix);
                        schema.set(parameter, file);
                        return;
                    }
                }
                throw new ParameterException(message + " (Parameter : " + parameter.getName() + ")");
            }

        }
    }

}