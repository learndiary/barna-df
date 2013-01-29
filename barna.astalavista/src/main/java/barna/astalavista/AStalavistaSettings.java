package barna.astalavista;

import barna.commons.log.Log;
import barna.commons.parameters.*;
import barna.model.ASEvent;
import barna.model.Graph;
import barna.model.Species;
import barna.model.Transcript;
//import barna.model.gff.GTFschema;
import barna.model.splicegraph.SplicingGraph;

import java.io.File;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.Iterator;

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
    public static final Parameter<Boolean> PRINT_PARAMETERS = Parameters.booleanParameter("PRINT_PARAMETERS",
            "print parameters and descriptions", false, null).longOption("printParameters").shortOption('h');

    /**
     * Path to the GTF reference annotation.
     */
    public static final Parameter<File> REF_FILE = Parameters.fileParameter("INPUT",
            "path to the GTF reference annotation",
            null, new ParameterValidator() {
/*      if (args[i].equals("-i")|| args[i].equals("--input")) {
            if (i+1>= args.length) {
                System.err.println("Hey, you forgot to specify the input file!");
                System.exit(-1);
            }
            file= new java.io.File(args[++i]);
            continue;
        }
*/
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
    }, null).longOption("input").shortOption('i');


    /**
     * File with GeneID parameters / splice site profiles.
     */
    public static final Parameter<File> GENEID_PARAM = Parameters.fileParameter("GENEID_PARAM",
            "name and path of a file with the GeneID models for splice sites",
            null, null, null).longOption("gparam");


    /**
     * File with variant information (vcf format).
     */
    public static final Parameter<File> VARIANTS = Parameters.fileParameter("VARIANTS",
            "name and path of a file with the variant information (vcf)",
            null, new ParameterValidator() {

            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File vcf = (File) schema.get(parameter);

                if (vcf == null || (!vcf.exists())) {
                    throw new ParameterException("VCF file not valid: "+ vcf== null? "null": vcf.getAbsolutePath());
                }
            }
    }).longOption("vcf");

    /**
     * Path to the GTF output annotation.
     */
/*   if (args[i].equals("-o")|| args[i].equals("--output")) {
        String s= args[++i];
        if (s.equalsIgnoreCase("stdout"))
            SplicingGraph.writeStdOut= true;
        else
            writerThread.outputFname= s+".gz";
        continue;
    } */
    public static final Parameter<File> OUT_FILE = Parameters.fileParameter("OUT_FILE",
            "keyword 'stdout' for standard output, or a path to the GTF output file",
            null, null, null).longOption("output").shortOption('o');


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

/*          if (args[i].equals("-g")|| args[i].equals("--genome")) {
                if (i+1>= args.length) {
                    System.err.println("You forgot to specify the genome!");
                    System.exit(-1);
                }
                MyFile checkFile= new MyFile(args[i+1]);
                if (checkFile.exists()&& checkFile.isDirectory())
                    Graph.overrideSequenceDirPath= checkFile.getAbsolutePath();
                else {
                    String[] s= args[++i].split("_");
                    if (s.length!= 2) {
                        System.err.println("Invalid genome directory: "+args[i+1]);
                        System.exit(-1);
                    }
                    Species spe= new Species(s[0]);
                    spe.setGenomeVersion(s[1]);
                    SplicingGraph.EventExtractor.setSpecies(spe);
                }
                acceptableIntrons= true;
                continue;
            } */

            File genomeFile = (File) schema.get(parameter);
            boolean req= requiresGenomicSequence(schema);
            if (!req)
                return;

            if (genomeFile == null) {
                throw new ParameterException("You have to specify a genome directory");
            }
            if (genomeFile.exists()) {
                Graph.overrideSequenceDirPath= genomeFile.getAbsolutePath();
            } else {
                Log.message("Trying to parse species_version pair");
                String[] s= genomeFile.getName().split("_");
                if (s.length!= 2) {
                    throw new ParameterException("The genome " + genomeFile.getAbsolutePath() + " could not be found!");
                }
                Species spe= new Species(s[0]);
                spe.setGenomeVersion(s[1]);
                AStalavista.setSpecies(spe);

            }

            // see AStalavista.validateParameter()
            // acceptableIntrons= true;

        }
    }, null).longOption("genome").shortOption('g');

    /**
     * Checks whether a folder with genomic sequences is necessary.
     * @param schema current parameter schema
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    protected static boolean requiresGenomicSequence(ParameterSchema schema) {
        if (schema.get(OUTPUT_SITESEQ))
            return true;
        if (ASEvent.isOutputFlankMode()&& Graph.overrideSequenceDirPath== null)
            return true;
        return false;
    }

    /**
     * Dimension of the AS events to be extracted,
     * retrieves 'complete' events <TM> for parameter
     * values < 2.
     */
    public static final Parameter<Integer> DIMENSION = Parameters.intParameter("DIMENSION",
            "Dimension of the AS events to be extracted, retrieves 'complete' events <TM>\n" +
            "for parameter values < 2",
            2, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {

/*          if (args[i].equals("-k")|| args[i].equals("--dimension")) {
                try {
                    int x= Integer.parseInt(args[++i]);
                    if (x< -1|| ((x> -1)&& (x< 2)))
                        System.err.println(args[i]+" is not a valid dimension, ignored");
                    SplicingGraph.EventExtractor.n= x;
                } catch (NumberFormatException e) {
                    System.err.println(args[i]+" is not a valid dimension, ignored"); // :)
                }
                continue;
            } */

            int v= (Integer) schema.get(parameter);
            if (v< 2)
                EventExtractor.n= (-1); // complete events
            else
                EventExtractor.n= 2;
        }
    }).longOption("dimension").shortOption('k');


    /**
     * Require 3'-complete transcripts.
     */
/*  if (args[i].equalsIgnoreCase("-3pc")) {
        ASEvent.check3Pcomplete= true;
        continue;
    } */
    public static final Parameter<Boolean> THREE_PRIME_COMPLETE = Parameters.booleanParameter("THREE_PRIME_COMPLETE",
            "require 3'-complete transcripts", false, null).longOption("3primeComplete");

    /**
     * Temporary directory
     */
    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR",
            "Temporary directory", new File(System.getProperty("java.io.tmpdir")), new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
/*          if (args[i].equals("-tmp")) {
                System.setProperty(Constants.PROPERTY_TMPDIR, args[++i]);
                continue;
            } */

            File tmp = parameterSchema.get(TMP_DIR);
            if(tmp == null){
                throw new ParameterException("No temp directory specified!");
            }
            if(!tmp.canWrite()){
                throw new ParameterException("The temp-directory " + tmp.getAbsolutePath() + " does not exist or is not writable!");
            }

        }
    }).longOption("tmp");

    /**
     * Switch on DS event retrieval.
     */
/*   if (args[i].equals("-ds")) {
        SplicingGraph.retrieveDSEvents= false;
     continue;
        } */
/*   if (args[i].equals("+ds")) {
        SplicingGraph.retrieveDSEvents= true;
        continue;
    } */
    public static final Parameter<Boolean> DS_EVENTS = Parameters.booleanParameter("DS_EVENTS",
            "do retrieve DS events", false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
            SplicingGraph.retrieveDSEvents= (Boolean) parameterSchema.get(parameter);
        }
    }).longOption("ds");

        /**
         * Switch on VS event retrieval.
         */
/*  if (args[i].equals("+vs")) {
        SplicingGraph.retrieveVSEvents= true;
        continue;
    }
    if (args[i].equals("-vs")) {
        SplicingGraph.retrieveVSEvents= false;
        continue; */
        public static final Parameter<Boolean> VS_EVENTS = Parameters.booleanParameter("VS_EVENTS",
                "do retrieve VS events", false, new ParameterValidator() {
            @Override
            public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
                SplicingGraph.retrieveVSEvents= (Boolean) parameterSchema.get(parameter);
            }
        }).longOption("vs");

    /**
     * Switch on external event retrieval.
     */
/*   if (args[i].equals("-ext")) {
        onlyInternal= true;	// net false
        continue;
    }
    if (args[i].equals("+ext")) {
        onlyInternal= false;
        continue;
    } */
    public static final Parameter<Boolean> EXT_EVENTS = Parameters.booleanParameter("EXT_EVENTS",
            "do retrieve external events", Boolean.FALSE, null).longOption("ext");


    /**
     * Different types of variation found in exon-intron structures
     * of transcripts:
     * <ul><li>AS= alternative splicing when comprising
     * at least one alternative splice site.
     * Types of alternative splicing can either be &quot;internal&quot;
     * and delimited by two common sites, or &quot;external&quot;
     * comprising at least one alternative splice site in addition
     * to alternative 5'- or 3' transcript structures.</li>
     * extending transcript structures by additional (splice) sites
     * to the 5'- or the 3'-end</li>
     * <li>VS= variable sites is any other form of sites that differ
     * between overlapping transcript structures</li>
     * @see barna.model.ASEvent#getType()
     * @see <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147">
     *     http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147</a><li>DS= additional splicing that are flanked by a common site and
     */
    public static enum EventTypes {ASExt,ASInt,DS,VS};


    /**
     * Parameter for counting reads that falls into specific elements
     */
    public static final Parameter<EnumSet<EventTypes>> EVENT_TYPES = Parameters.enumSetParameter(
            "EVENT_TYPES",
            " Type of events that is considered",
            EnumSet.of(EventTypes.ASInt),
            EventTypes.class,
            null).longOption("events").shortOption('e');

    /**
        * Flag to suppress AS event retrieval.
        */
/*   if (args[i].equals("-as")) {
        SplicingGraph.retrieveASEvents= false;
        continue;
    } */
/*  if (args[i].equals("+as")) {
        SplicingGraph.retrieveASEvents= true;
        continue;
    } */
    public static final Parameter<Boolean> NO_AS_EVENTS = Parameters.booleanParameter("NO_AS_EVENTS",
            "don't retrieve AS events", true, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
            SplicingGraph.retrieveASEvents= (Boolean) parameterSchema.get(parameter);
        }
    }).longOption("noas");


    /**
     * Check nonsense-mediated decay conditions.
     */
/*  if (args[i].equalsIgnoreCase("-nmd")) {
        ASEvent.checkNMD= true;
        continue;
    } */
    public static final Parameter<Boolean> NMD = Parameters.booleanParameter("NMD",
            "check nonsense-mediated decay conditions", false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
            ASEvent.checkNMD= (Boolean) parameterSchema.get(parameter);
        }
    }).longOption("nmd");

    /**
     * Confidence level of the intron to be retrieved.
     */
    public static final Parameter<Integer> INTRON_CONFIDENCE = Parameters.intParameter("INTRON_CONFIDENCE",
            "Confidence level of the intron to be retrieved",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {

        /*  if (args[i].equalsIgnoreCase("-ic")|| args[i].equalsIgnoreCase("--intronConfidence")) {
                if (i+1== args.length)
                    System.err.println("You did not provide an intron confidence.");
                try {
                    SplicingGraph.intronConfidenceLevel= Byte.parseByte(args[i+1]);
                    ++i;	// ignore if missing
                    acceptableIntrons= true;
                } catch (NumberFormatException e) {
                    System.err.println("Intron confidence must be an integer value, you gave me "+args[i+1]);
                }
                continue;
            } */
            int x= schema.get(INTRON_CONFIDENCE);
            if (x> Byte.MAX_VALUE|| x< 0)
                throw new IllegalArgumentException("Invalid confidence level "+ x);
            SplicingGraph.intronConfidenceLevel= (byte) x;
        }
    }).longOption("intronConfidence");

    /**
     * Confidence level of the intron to be retrieved.
     */
    public static final Parameter<Integer> EDGE_CONFIDENCE = Parameters.intParameter("EDGE_CONFIDENCE",
            "Confidence level of the intron to be retrieved",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {

/*  if (args[i].equalsIgnoreCase("-ec")|| args[i].equalsIgnoreCase("--edgeConfidence")) {
        if (i+1== args.length)
            System.err.println("You did not provide an edge confidence.");
        try {
            Transcript.setEdgeConfidenceLevel(Byte.parseByte(args[i + 1]));
            ++i;	// ignore if missing
        } catch (NumberFormatException e) {
            System.err.println("Exon confidence must be an integer value, you gave me "+args[i+1]);
        }
        continue;
    } */
            int x= schema.get(EDGE_CONFIDENCE);
            if (x> Byte.MAX_VALUE|| x< 0)
                throw new IllegalArgumentException("Invalid confidence level "+ x);
            Transcript.setEdgeConfidenceLevel((byte) x);
        }
    }).longOption("edgeConfidence");

    /**
     * Switch on flank type output.
     */
/*  if (args[i].equals("--flankType")) {
        ASEvent.setOutputFlankMode(true);
        continue;
    } */
    public static final Parameter<Boolean> FLANK_TYPE = Parameters.booleanParameter("FLANK_TYPE",
            "switch on flank type output", false, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema parameterSchema, Parameter parameter) throws ParameterException {
            ASEvent.setOutputFlankMode((Boolean) parameterSchema.get(parameter));
        }
    }).longOption("flankType");



/*
    for (int i = 0; i < args.length; i++) {

        // -c is now RESERVED for "command" multicaster
//				if (args[i].equals("-c")|| args[i].equals("--canonical")) {
//					canonicalSS= true;
//					continue;
//				}
        //			if (args[i].equals("-c")|| args[i].equals("--cluster")) {
        //				readAheadLimit= Integer.parseInt(args[++i]);
        //				continue;
        //			}


        }

        // reactivated 20100112



    }




    barna.io.gtf.GTFwrapper checkReader= new GTFwrapper(file.getAbsolutePath());
    //System.err.println("DEBUG -- temporarily deactivated file check");
    if (!checkReader.isApplicable()) {
        System.err.println("sorting input file, temporary directory "+System.getProperty(Constants.PROPERTY_TMPDIR));
        file= checkReader.sort();  // WAS: GFFreader.createSortedFile();
        System.err.println("Here is a sorted version of your file: "+file.getAbsolutePath());
    }


    return file;
*/

    @Override
    public void write(OutputStream out) {
        super.write(out);
/*        if (isPrintParameters()) {
            FluxSimulatorSettings settings = new FluxSimulatorSettings();
            settings.write(System.out);
            return false;

            System.err.println("Here is a list of the options I understand:\n");
            System.err.println("-i, --input <input file>");
            System.err.println("This is a bit important, I cannot work without an input annotation. I want a GTF file with " +
                    "transcript annotations (exon features, with a mandatory optional attribute named \'transcript_id\') " +
                    "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it has to be in " +
                    "all lines of the file in column #10. The rest of the file should comply with the standard as specified " +
                    "at http://mblab.wustl.edu/GTF2.html.\n"+
                    "There may also be CDS features, but they become only interesting when checking for additional things " +
                    "as NMD probability etc.."+
                    "\n");
            System.err.println("-o, --output <output file|\'stdout\'>");
            System.err.println("Optional, the name of the output file (fully qualified path) OR the keyword \'stdout\' for " +
                    "writing everything to the standard output stream. " +
                    "If nothing is specified, the output will be written to a file \'<input file>_astalavista.gtf.gz\'. " +
                    "\n");
            System.err.println("-g, --genome <path to directory>");
            System.err.println("Path to the directory containing sequence files corresponding to the <seqname> field " +
                    "in the input GTF. A genome directory is required if a intron confidence value is specified." +
                    "\n");
            System.err.println("-k, --dimension <int value>");
            System.err.println("Dimension >1 of the events to be extracted. Default is 2 (i.e., \'pairwise events\'). " +
                    "\n");
            System.err.println("-tmp");
            System.err.println("Set temporary directory" +
                    "\n");
            System.err.println("-ext, +ext");
            System.err.println("(De-)activate external events, i.e. events that include the transcript start or the " +
                    "poly-adenylation site" +
                    "\n");
            System.err.println("-ic, --intronConfidence [int value]");
            System.err.println("Level of intron confidence. The default is to trust all introns. Introns are assigned " +
                    "a confidency class:\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all introns of confidency level > intronConfidence will be checked for proper splice sites when extracting events." +
                    "\n");
            System.err.println("-as, +as");
            System.err.println("Deactivate (\'-as\') or activate (\'+as\') the retrieval of Alternative Splicing events. See documentation " +
                    "for the definition of events that suffice an alternative splicing event." +
                    "\n");
            System.err.println("-ds, +ds");
            System.err.println("Deactivate (\'-ds\') or activate (\'+ds\') the retrieval of aDditional splicing events. See documentation " +
                    "for the definition of events that suffice an additional splicing event." +
                    "\n");
            System.err.println("-s, --seqsite");
            System.err.println("Output splice site sequences with events. Requires a reference genome."+
                    "\n");
            System.err.println("--flankType");
            System.err.println("Output the type of the event flanks, i.e., \'constitutive\' or \'alternative\'."+
                    "\n");

            // reactivated on 20100112
            System.err.println("-ec, --edgeConfidence [int value]");
            System.err.println("Level of confidence for edges (i.e., annotated transcription starts/poly-adenylation sites). " +
                    "The default is to trust no annotated edge and to extend overlapping first/last exons of a transcript to " +
                    "their most extreme position. :\n" +
                    "\t 0 if 'RefSeq' appears in the source field of the annotation\n" +
                    "\t 1 if 'mRNA' appears in the source field of the annotation\n" +
                    "\t 2 if 'EST' appears in the source field of the annotation\n" +
                    "\t 3 if if none of the above applies\n" +
                    "all transcript edges of confidency level > edgeConfidence will be extended in case the annotation shows " +
                    "another exon with the same adjacent splice site and an earlier/later start/end." +
                    "\n");
            System.err.println("AStalavista.");
            System.exit(0);
        }


        return true;
    } */

    }

    /**
     * Output splice site sequences.
     * @deprecated replace by OutputEvent.SITESCORES
     */
/*  if (args[i].equals("-s")|| args[i].equals("--seqsite")) {
        SplicingGraph.outputSeq= true;
        continue;
    } */
    public static final Parameter<Boolean> OUTPUT_SITESEQ = Parameters.booleanParameter("OUTPUT_SITESEQ",
            "output splice site sequences", false, null).longOption("seqsite").shortOption('s');


    /**
     * Score splice site sequences.
     * @deprecated replace by OutputSite.SITESCORES
     */
    public static final Parameter<File> SCORE_SITES = Parameters.fileParameter("SCORE_SITES",
            "score splice site sequences", null, null, null).longOption("scoresites");


}
