package fbi.genome.sequencing.rnaseq.simulation;


import fbi.commons.file.FileHelper;
import fbi.commons.parameters.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Flux Simulator settings
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 * @author Micha Sammeth (gmicha@gmail.com)
 */
public class FluxSimulatorSettings extends ParameterSchema {

    /**
     * Loci separator
     */
    public static final char SEP_LOC_TID = '@'; // todo refactor this to profiler and use the profilers getGlobalID where possible

    /**
     * Helper to parse relative filenames
     */
    private static RelativePathParser relativePathParser = new RelativePathParser();


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
    public static final Parameter<File> REF_FILE = Parameters.fileParameter("REF_FILE_NAME", "GTF reference file", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File refFile = (File) schema.get(parameter);
            if (refFile == null) throw new ParameterException("You have to specify a reference file");
            if (!refFile.exists())
                throw new ParameterException("The reference file " + refFile.getAbsolutePath() + " could not be found!");
        }
    }, relativePathParser);
    public static final Parameter<File> PRO_FILE = Parameters.fileParameter("PRO_FILE_NAME", "Target Profiler file", null, new FileValidator("pro"), relativePathParser);
    public static final Parameter<File> LIB_FILE = Parameters.fileParameter("LIB_FILE_NAME", "Target library file", null, new FileValidator("lib"), relativePathParser);
    public static final Parameter<File> SEQ_FILE = Parameters.fileParameter("SEQ_FILE_NAME", "Target sequences file", null, new FileValidator("bed"), relativePathParser);
    public static final Parameter<File> GEN_DIR = Parameters.fileParameter("GEN_DIR", "The Genome directory", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File genomeFile = (File) schema.get(parameter);
            if (genomeFile == null) throw new ParameterException("You have to specify a genome directory");
            if (!genomeFile.exists())
                throw new ParameterException("The genome directory " + genomeFile.getAbsolutePath() + " could not be found!");

        }
    }, relativePathParser);
    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR", "Temporary directory", new File(System.getProperty("java.io.tmpdir")));
    public static final Parameter<File> ERR_FILE = Parameters.fileParameter("ERR_FILE_NAME", "Error model file", relativePathParser);


    public static final Parameter<Boolean> FASTQ = Parameters.booleanParameter("FASTQ", "Create FastQ output", false);
    /*
    Expression parameters
     */
    public static final Parameter<Boolean> LOAD_CODING = Parameters.booleanParameter("LOAD_CODING", "", true, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            if(!schema.get(LOAD_CODING) && !schema.get(LOAD_NONCODING)) throw new ParameterException("Sorry, but either LOAD_CODING or LOAD_NONCODING has to be enabled !");
        }
    });
    public static final Parameter<Boolean> LOAD_NONCODING = Parameters.booleanParameter("LOAD_NONCODING", "", true, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            if(!schema.get(LOAD_CODING) && !schema.get(LOAD_NONCODING)) throw new ParameterException("Sorry, but either LOAD_CODING or LOAD_NONCODING has to be enabled !");
        }
    });
    public static final Parameter<Long> NB_MOLECULES = Parameters.longParameter("NB_MOLECULES", "", 5000000);
    public static final Parameter<Double> EXPRESSION_K = Parameters.doubleParameter("EXPRESSION_K", "", -0.6);
    public static final Parameter<Double> EXPRESSION_X0 = Parameters.doubleParameter("EXPRESSION_X0", "", 50000000);
    public static final Parameter<Double> EXPRESSION_X1 = Parameters.doubleParameter("EXPRESSION_X1", "", 9500);
    public static final Parameter<Double> TSS_MEAN = Parameters.doubleParameter("TSS_MEAN", "", 25d, new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            if(schema.get(TSS_MEAN) <= 0){
                throw new ParameterException("TSS_MEAN must be > 0");
            }
        }
    });
    public static final Parameter<Double> POLYA_SHAPE = Parameters.doubleParameter("POLYA_SHAPE", "", 2d, 0.0, Double.MAX_VALUE, null);
    public static final Parameter<Double> POLYA_SCALE = Parameters.doubleParameter("POLYA_SCALE", "", 300d, 0.0, Double.MAX_VALUE, null);

    /*
    Fragementation
     */
    public static final Parameter<Boolean> FRAGMENTATION = Parameters.booleanParameter("FRAGMENTATION", "", true);
    public static final Parameter<FragmentationMethod> FRAG_METHOD = Parameters.enumParameter("FRAG_METHOD", "" +
            "Method applied for Fragmentation\n" +
            "[NB] Nebulization fragmentation method.\n" +
            "[UR] Uniformal random fragmentation method.\n" +
            "[EZ] Enzymatic digestion as fragmentation method.", FragmentationMethod.NB, new ParameterValidator() {
        @Override
        public void validate(final ParameterSchema schema, final Parameter parameter) throws ParameterException {
            if(schema.get(FRAG_METHOD) == FragmentationMethod.EZ){
                if(schema.get(FRAG_SUBSTRATE) == Substrate.RNA){
                    throw new ParameterException("Enzymatic digestion is not supported for RNA substrate!");
                }
                if(schema.get(FRAG_EZ_MOTIF) == null || !schema.get(FRAG_EZ_MOTIF).canRead()){
                    throw new ParameterException("You have to specify FRAG_EZ_MOTIF in order ot use Enzymatic digestion");
                }
            }

            if(schema.get(FRAG_METHOD) == FragmentationMethod.NB){
                if(schema.get(FRAG_SUBSTRATE) == Substrate.RNA){
                    throw new ParameterException("Sorry, but nebulizing RNA is not supported!");
                }
            }

        }
    });
    public static final Parameter<Substrate> FRAG_SUBSTRATE = Parameters.enumParameter("FRAG_SUBSTRATE", " Parameter specifying the substrate of fragmentation.", Substrate.DNA, null);

    /*
    Enzymatic
     */
    public static final Parameter<File> FRAG_EZ_MOTIF = Parameters.fileParameter("FRAG_EZ_MOTIF", "The motif description for enzymatic digestion", relativePathParser);


    /*
    Nebulization
     */
    /**
     * Parameter for the threshold on molecule length that cannot
     * be broken by the shearfield during nebulization.
     */
    public static final Parameter<Double> FRAG_NB_LAMBDA = Parameters.doubleParameter("FRAG_NB_LAMBDA",
            "Parameter for the threshold on molecule length that cannot " +
                    "be broken by the shearfield during nebulization.",
            900);

    public static final Parameter<Double> FRAG_NB_THOLD = Parameters.doubleParameter("FRAG_NB_THOLD",
            "Parameter denoting the threshold on molecule population still " +
                    "breaking when determining convergence of iterative nebulizaiton.", 0.1);

    /**
     * Parameter for providing the fraction of molecule that marks
     * the standard deviation of the breakpoint in nebulization
     * (mean=center of the molecule).
     */
    public static final Parameter<Double> FRAG_NB_SIGMA = Parameters.doubleParameter("FRAG_NB_SIGMA",
            "Parameter for providing the fraction of molecule that marks " +
                    "the standard deviation of the breakpoint in nebulization " +
                    "(mean=center of the molecule).",
            0.05);

    public static final Parameter<Double> FRAG_NB_M = Parameters.doubleParameter("FRAG_NB_M", "Parameter " +
            "specifying the strength of the " +
            "nebulization shearfield.", Double.NaN); // todo : check default

    /*
     RNA Hydrolysis (Uniform-Random)
     */
    public static final Parameter<Double> FRAG_UR_ETA = Parameters.doubleParameter("FRAG_UR_ETA", "", Double.NaN);
    public static final Parameter<Double> FRAG_UR_DELTA = Parameters.doubleParameter("FRAG_UR_DELTA", "", Double.NaN);
    public static final Parameter<Double> FRAG_UR_D0 = Parameters.doubleParameter("FRAG_UR_D0", "", 1.0, 1.0, Double.MAX_VALUE, null);

    /*
     Reverse Transcription
     */
    public static final Parameter<Boolean> RTRANSCRIPTION = Parameters.booleanParameter("RTRANSCRIPTION", "Switch on/off Reverse Transcription", true);// todo : default ?
    public static final Parameter<RtranscriptionMode> RT_PRIMER = Parameters.enumParameter("RT_PRIMER", "", RtranscriptionMode.RH, null);
    public static final Parameter<Integer> RT_MIN = Parameters.intParameter("RT_MIN", "Minimum length observed after " +
            "reverse transcription of full-length transcripts.", 500, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int min = schema.get(RT_MIN);
            int max = schema.get(RT_MAX);
            if(min > max) throw new ParameterException("RT_MIN must me <= RT_MAX");
        }
    });
    public static final Parameter<Integer> RT_MAX = Parameters.intParameter("RT_MAX", "Maximum length observed after " +
            "reverse transcription of full-length transcripts.", 5500,  new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int min = schema.get(RT_MIN);
            int max = schema.get(RT_MAX);
            if(min > max) throw new ParameterException("RT_MIN must me <= RT_MAX");
        }
    });
    public static final Parameter<Double> RT_GC_LO = Parameters.doubleParameter("RT_GC_LO", "Minimum GC content for RNA " +
            "stretches to be successfully reversely transcribed.", 0.4, 0.0, 1.0, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            double lo = schema.get(RT_GC_LO);
            double hi = schema.get(RT_GC_HI);
            if(lo > hi) throw new ParameterException("RT_GC_HI must be >= RT_GC_LO");
        }
    });
    // todo : remove ?
    public static final Parameter<Double> RT_GC_HI = Parameters.doubleParameter("RT_GC_HI",
            "GC content where reverse transcription saturates", 0.7, 0.0, 1.0, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            double lo = schema.get(RT_GC_LO);
            double hi = schema.get(RT_GC_HI);
            if(lo > hi) throw new ParameterException("RT_GC_HI must be >= RT_GC_LO");
        }
    });

    public static final Parameter<Boolean> RT_LOSSLESS = Parameters.booleanParameter("RT_LOSSLESS", "Always force RT ", true);

    public static final Parameter<File> RT_MOTIF = Parameters.fileParameter("RT_MOTIF", "", relativePathParser);

    /*
    Size Selection
     */
    public static final Parameter<Boolean> FILTERING = Parameters.booleanParameter("FILTERING", "", false);
    public static final Parameter<File> SIZE_DISTRIBUTION = Parameters.fileParameter("SIZE_DISTRIBUTION", "Describes " +
            "the distribution of fragments after filtering", relativePathParser);
    public static final Parameter<SizeSamplingModes> SIZE_SAMPLING = Parameters.enumParameter("SIZE_SAMPLING",
            "Describes the method for subsampling fragments in order to meet the characteristics " +
                    "of the filter Distribution (see SIZE_DISTRIBUTION)", SizeSamplingModes.AC, null);

    /*
      Sequencing
     */
    public static final Parameter<Long> READ_NUMBER = Parameters.longParameter("READ_NUMBER", "Number of reads", 5000000);
    public static final Parameter<Integer> READ_LENGTH = Parameters.intParameter("READ_LENGTH", "The Read length", 36);
    public static final Parameter<Boolean> PAIRED_END = Parameters.booleanParameter("PAIRED_END", "Pair end reads", false);


    /**
     * The parameter file
     */
    private File parameterFile;

    /**
     * Load the setting from a file
     *
     * @param f the file
     * @return settings the loaded and validated settings
     * @throws Exception in case of a parser or validation error
     */
    public static FluxSimulatorSettings createSettings(File f) throws Exception {
        if (f == null) throw new NullPointerException("Null parameter file not permitted!");
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
            if (in != null) try {
                in.close();
            } catch (IOException e) {
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
                        if (!suffix.startsWith(".")) suffix = "." + suffix;
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