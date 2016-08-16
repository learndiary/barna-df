package barna.astalavista;

import barna.commons.parameters.*;
import barna.model.Transcript;
//import barna.model.gff.GTFschema;
import barna.model.constants.Constants;

import java.io.File;
import java.io.OutputStream;
import java.util.EnumSet;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 6/18/12
 * Time: 1:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class AStalavistaSettings extends ParameterSchema /*GTFschema*/ {

    // TODO {SplicingGraph.class, null, SJextractor.class, AttributeExtractor.class};	// null= LaVista.class

    /**
     * Print parameters and descriptions.
     */
//    public static final Parameter<Boolean> HELP = Parameters.booleanParameter("HELP",
//            "print parameters and descriptions", false, null).longOption("printParameters").shortOption('h');

    /**
     * Path to parameter file.
     */
    public static final Parameter<File> PAR_FILE = Parameters.fileParameter("PAR_FILE",
            "path to the parameter file",
            null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File f = (File) schema.get(parameter);
            if (f == null) {
                throw new ParameterException("Parameter file cannot be null!");
            }
            if (!f.exists()) {
                throw new ParameterException("The parameter file " + f.getAbsolutePath() + " could not be found!");
            }
        }
    }, null).longOption("par").shortOption('p');

    /**
     * Path to the reference annotation.
     */
    public static final Parameter<File> IN_FILE = Parameters.fileParameter("IN_FILE",
            "path to the GTF reference annotation",
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
    }, null).longOption("in").shortOption('i');


    /*
     * Parameter specifying a collection of 'source' tags that are read from the input,
     * possibly other flag for non-/coding transcripts to be considered
     *
    public static final Parameter<EnumSet<EventOptions>> IN_OPTIONS = Parameters.enumSetParameter(
            "IN_OPTIONS",
            "Toggle criteria for elements of the input",
            EnumSet.noneOf(EventOptions.class),
            EventOptions.class,
            null).longOption("io");
     */

    /**
     * Path to the directory with the genomic sequences,
     * i.e., one fasta file per chromosome/scaffold/contig
     * with a file name corresponding to the identifiers of
     * the first column in the GTF annotation.
     */
    public static final Parameter<File> CHR_SEQ = Parameters.fileParameter("CHR_SEQ",
                    "path to the directory with the genomic sequences,\n" +
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
    }, null).longOption("chr").shortOption('c');

    /**
     * The temporary directory.
     */
    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR", "The temporary directory", new File(System.getProperty("java.io.tmpdir")), new ParameterValidator() {
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
            if (file == null) {
                schema.set(parameter, new File(System.getProperty(Constants.PROPERTY_TMPDIR)));
            }
            if (!file.exists()) {
                throw new ParameterException("The temporary directory " + file.getAbsolutePath()
                        + " could not be found!");
            }
            if (!file.canWrite()) {
                throw new ParameterException("The temporary directory " + file.getAbsolutePath()
                        + " is not writable!");
            }

        }
    }).longOption("tmp");

    /**
     * Checks whether a folder with genomic sequences is necessary in order
     * to complete all tasks specified with <code>this</code> parameter set.
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    public boolean requiresGenomicSequence() {
        return false;
    }

}
