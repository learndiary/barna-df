package barna.astafunk;

/**
 * @author vitorcoelho
 * @version 2
 * @since 09/11/15
 */

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.commons.parameters.ParameterException;
import barna.commons.parameters.ParameterSchema;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.parser.HMMParser;
import barna.astafunk.parser.HeuristicTableParser;
import barna.astafunk.utils.FunkSettings;
import barna.astafunk.utils.Utils;
import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.Graph;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Main class of ASTAFUNK. In this class, parameters are initialized; input files are read and parsed;
 * and search is performed.
 */
public class Astafunk implements Tool<Void> {

    public static FunkSettings funkSettings = null;
    /**
     * Array of genes from GTF annotation.
     */
    private static Gene[] genes = null;

    /**
     * Hashmap of the domain reference file.
     */
    private static HashMap<String, List<String>> heuristicTable;

    /**
     * Hashmap of the domain reference file
     */

    private static HashMap<String, ProfileHMM> hmmHash;

    /**
     * Wrapper for the GTF annotation.
     */
    private static GTFwrapper wrapper = null;

    /**
     * Main class.
     * @param args List of options for ASTAFUNK.
     */
    public static void main(String [] args) throws Exception {

        Astafunk myFunk= new Astafunk();

        // INIT
        // construct to register parameters in JSAP
        List<Parameter> parameter = JSAPParameters.getJSAPParameter(new FunkSettings());

        JSAP jsap = JSAPParameters.registerParameters(parameter);

        // parse
        try{
            JSAPResult toolParameter = jsap.parse(args);
            if (!myFunk.validateParameter(toolParameter)){
                funkSettings.write(System.out);
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Parameter error : " + e.getMessage(), e);
            e.printStackTrace();
            System.exit(-1);
        }
    }

    /**
     * Computes a result, or throws an exception if unable to do so.
     *
     * @return computed result
     * @throws Exception if unable to compute a result
     */
    @Override
    public Void call() throws Exception {

        Log.message("ASTAFUNK (v.1 beta) - " + getDescription() + "\nCoelho, VL and Sammeth, M. (vitorlimac2@gmail.com, micha@sammeth.net)");

        printParameters();

        long t0= System.currentTimeMillis();
        Log.message("# started\t" + new Date(t0));

        callBegin();

        ExecutorService executor = Executors.newFixedThreadPool(funkSettings.get(FunkSettings.CPU));

        if(!isTranscriptReferencePrint())
            Utils.printHeaderHitList();

        Log.message("Analysing genes...");

       // Utils.initBar();

        for(Gene gene: genes){
            Runnable search = new Tsearch(gene);
            executor.execute(search);
        }
        executor.shutdown();

        while (!executor.isTerminated()) {
     //       Utils.updateBar(Utils.finishedThread, genes.length);
        }

        callFinish();

        Log.progressFinish("done.", true);

        Log.message("took " + ((System.currentTimeMillis() - t0) / 1000) + " sec.");

        System.exit(0);

        return null;
    }

    public static void callBegin() throws Exception{

        /**
         * parse gtf and get genes
         */
        parseGTF();

        if(!isTranscriptReferencePrint()){

             /**
             * Parsing profile HMM(s) file.
             */
            setHmmHash(parseHMMs());

            if(!isExhaustive()){
                /**
                 * Parse heuristic table
                 */
                setHeuristicTable(parseHeuristicTable());
            }
        }
    }

    protected FunkSettings getSettings() {
        if (funkSettings== null) {
            return new FunkSettings();
        }
        return funkSettings;    // else
    }

    /**
     * @return a brief description of the tool's functionality
     */
    @Override
    public String getDescription() {
        return "Alternative Splicing Transcriptional Analyses with Functional Knowledge: " +
                "Search Pfam HMMs against alternatively spliced regions.";
    }

    /**
     * @return a verbose description of the tool
     */
    @Override
    public String getLongDescription() {
        return "Search HMM-profiles of protein families (Pfam) on alternatively spliced genes.";
    }

    /**
     * @return the unique name of the tool
     */
    @Override
    public String getName() {
        return "astafunk";
    }

    /**
     * List of parameters for the tool.
     * @return parameter list
     */
    @Override
    public List<Parameter> getParameter() {
        // converts parameter file parameters to CLI parameters

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        FunkSettings settings= getSettings();
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter parameter : pars) {

            Class c= parameter.getType();
            Parameter p;
            if (c.toString().toLowerCase().contains("enum"))
                System.currentTimeMillis();
            if (c.equals(Boolean.class)) {
                p= JSAPParameters.switchParameter(
                        parameter.getLongOption(),
                        parameter.getShortOption())
                        .defaultValue(parameter.getDefault().toString())
                        .type(c)
                        .help(parameter.getDescription())
                        .get();
            } else {
                p= JSAPParameters.flaggedParameter(
                        parameter.getLongOption(),
                        parameter.getShortOption())
                        .type(c)
                        .help(parameter.getDescription())
                        .valueName(parameter.getName())
                        .get();
            }
            // TODO required() not implemented
            if (parameter.getLongOption()!= null|| parameter.getShortOption()!= 0)
                parameters.add(p);
        }
        return parameters;
    }

    /**
     * Validate arguments from command line.
     * @param args result from parsing the command line
     * @return True if parameters are valid.
     */
    @Override
    public boolean validateParameter(JSAPResult args) {
        return validateParameter(new FunkSettings(), args);
    }

    /**
     * Checks CLI parameters with respect to their validity, fills the
     * settings instance provided with values from the parameter file
     * and from the command line.
     * @param schema a non-<code>null</code>
     * @param args result from parsing the command line
     * @return <code>true</code> if everything is ok with the parameters,
     * <code>false</code> otherwise
     */
    public boolean validateParameter(ParameterSchema schema, JSAPResult args) {

        // TODO think about pulling up to interface / abstract class in commons

        // non-null settings are to be assumed
        if (schema== null|| !(schema instanceof FunkSettings)) {
            Log.error("Must provide an instance of FunkSettings!");
            return false;
        }
        funkSettings = (FunkSettings) schema;

        try {
            funkSettings= (FunkSettings) ParameterSchema.create(funkSettings,
                    JSAPParameters.getParameterMap(funkSettings, args));
        } catch (ParameterException e) {
            Log.error(e.getMessage(), e);
            return false;
        }

        // CPU
        if(funkSettings.get(FunkSettings.CPU)< 0){
            Log.error("Option " + FunkSettings.CPU.getName()+ " must be a positive " +
                    "integer number.");
            return false;
        }

        //GENOME
        if (funkSettings.get(FunkSettings.GENOME)!= null) {
            File f= funkSettings.get(FunkSettings.GENOME);
            if (f.exists()&& f.isDirectory()) {
                Graph.overrideSequenceDirPath= f.getAbsolutePath();
            } else {
                Log.message("Trying to read genome FASTA files. Check directory or file names.");
                String[] s= f.getName().split("_");
                if (s.length!= 2) {
                    Log.error(f.getAbsolutePath() + " is not a valid species name");
                    return false;
                }
                Log.error("Systematic species names not supported!");
            }
        }

        //GTF
        if (funkSettings.get(FunkSettings.GTF)== null) {
            Log.error("You need GTF annotation file (--gtf <FILE>) with transcript annotations (exon features, with a \n" +
                            "mandatory optional attribute named \'transcript_id\'\n) " +
                            "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it\n" +
                            "has to be in all lines of the file in column #10. The rest of the file should comply with the\n" +
                            "standard as specified at http://mblab.wustl.edu/GTF2.html.\n"+
                            "There may also be CDS features, but they become only interesting when checking for additional things\n" +
                            "as NMD probability etc.. "
            );
            return false;
        }

        //HMM_FILE: if it is not reference transcript search, this parameter must exist
        if(!isTranscriptReferencePrint()) {
            if (funkSettings.get(FunkSettings.HMM_FILE) == null) {
                Log.error("You forgot a profile HMM file (--hmm <FILE>).");
                // print help
                return false;
            }
        }else{

            funkSettings.set("CPU", 1);
        }

        //OVERLAPPING
        if(funkSettings.get(FunkSettings.OVERLAPPING)<0
                || funkSettings.get(FunkSettings.OVERLAPPING)>1){
            Log.error("The overlapping hit ratio must be between 0.0 and 1.0.");
            return false;
        }

        //REFERENCE_FILE:
        if(!isExhaustive() && !isTranscriptReferencePrint()) {
            if (funkSettings.get(FunkSettings.REFERENCE_FILE) == null) {
                Log.error("You are running a heuristic search, " +
                        "so you forgot a domain reference file (--reference/-r <FILE>). " +
                        "If you do not have a reference domain file, " +
                        "you can run a exhaustive search (--exh/-e). See user's guide to create a " +
                        "domain reference file.");
                // print help
                return false;
            }
        }
        return true;
    }

    /**
     * Parse profile HMM file;
     * @return List of profile HMMs
     * @throws IOException
     */
    private static HashMap<String, ProfileHMM> parseHMMs() throws IOException {
        HMMParser hmmParser = new HMMParser(funkSettings.get(FunkSettings.HMM_FILE).getAbsolutePath());
        return hmmParser.parseAll2();
    }

    /**
     * Initialize ASta settings; parse GTF annotation file and retrieve the genes.
     */
    public static void parseGTF() {

        File f = funkSettings.get(FunkSettings.GTF);

        wrapper = new GTFwrapper(f);

        if (!wrapper.isApplicable()) {
            f = wrapper.sort();
            wrapper = new GTFwrapper(f);
        }

        //Read all genes
        wrapper.setReadAll(true);

        //wrapper.setChromosomeWise(false);
        //wrapper.setGeneWise(false);
        wrapper.setReadAheadTranscripts(-1);

        //Read the file
        wrapper.read();

        //Get genes
        genes = wrapper.getGenes();

        //Define path to fastas
        Graph.overrideSequenceDirPath =
                funkSettings.get(FunkSettings.GENOME).getAbsolutePath();
    }

    /**
     * Parse a file (heuristic table) with gene entries and the respective HMM ACCs.
     * @return A hashmap of gene ID and the respective profile HMM ACC list.
     * @throws IOException
     */
    public static HashMap<String, List<String>> parseHeuristicTable() throws IOException {
        if(isExhaustive())
            return null;
        HeuristicTableParser heuristicTableParser = new HeuristicTableParser(Astafunk.funkSettings
                .get(FunkSettings.REFERENCE_FILE).getAbsolutePath());
        return heuristicTableParser.parseAll();
    }

    /**
     * @return A hash map - heuristic table - where the key is a gene/transcript ID and the
     * value is a one or a list of Pfam accession numbers (ACCs).
     */
    public static HashMap<String, List<String>> getHeuristicTable() {
        return heuristicTable;
    }

    /**
     *
     * @return List of profile HMMs.
     */
    public static HashMap<String, ProfileHMM> getHmmHash() {
        return hmmHash;
    }

    /**
     * @param heuristicTable Hashmap of gene/transcript IDs and Pfam accession numbers (ACC) to set.
     */
    private static void setHeuristicTable(HashMap<String, List<String>> heuristicTable) {
        Astafunk.heuristicTable = heuristicTable;
    }

    /**
     *
     * @param hmmHash List of profile HMMs to set.
     */
    private static void setHmmHash(HashMap<String, ProfileHMM> hmmHash) {
        Astafunk.hmmHash = hmmHash;
    }

    /**
     * Clear hmm list and heuristic table; close reader.
     * @throws Exception
     */
    protected static void callFinish() throws Exception{
        wrapper.close();
    }

    /**
     *
     * @return true if verbose mode is activated.
     */
    public static boolean isVerbose(){
        return funkSettings.get(FunkSettings.VERBOSE);
    }

    /**
     *
     * @return true if local search is activated.
     */
    public static boolean isLocal(){ return funkSettings.get(FunkSettings.LOCAL);
    }

    /**
     *
     * @return True if transcript reference printer is selected.
     */
    public static boolean isTranscriptReferencePrint(){
        return funkSettings.get(FunkSettings.REFTRANSCRIPT);
    }

    /**
     *
     * @return True if exhaustive search is selected.
     */
    public static boolean isExhaustive(){
        return funkSettings.get(FunkSettings.EXHAUSTIVE);
    }

    /**
     *
     * @return Overlapping threshold.
     */
    public static double getOverlappingThreshold(){
        return funkSettings.get(FunkSettings.OVERLAPPING);
    }

    /**
     * Print ASTAFUNK options
     */

    public static void printParameters(){

        Log.message("Input:");
        Log.message("\tGenome files: " + funkSettings.get(FunkSettings.GENOME).getAbsolutePath());
        Log.message("\tAnnotation: " + funkSettings.get(FunkSettings.GTF).getAbsolutePath());

        if(isTranscriptReferencePrint()) {
            Log.message("\nOptions:");
            Log.message("\tMode: Reference transcript output.");
        }else {
            Log.message("\tHMM-Profiles: " + funkSettings.get(FunkSettings.HMM_FILE).getAbsolutePath());
            Log.message("\tReference domain file: " + funkSettings.get(FunkSettings.REFERENCE_FILE).getAbsolutePath());

            Log.message("\nOptions:");
            if (isLocal())
                Log.message("\tProfile search: Local");
            else
                Log.message("\tProfile search: Glocal");

            if (isExhaustive())
                Log.message("\tMode: Exhaustive");
            else if(isNaive())
                Log.message("\tMode: Näive");
            else if(isConstitutive())
                Log.message("\tMode: Search domains on constitutive genomic regions.");
            else
                Log.message("\tMode: Heuristic");

            Log.message("\tOverlapping hit threshold: " + funkSettings.get(FunkSettings.OVERLAPPING));

            if(isBestGeneHits()){
                Log.message("\tOutput: Non-overlapped hits of the AS genes.");
            }else if(isAllTranscriptHits()){
                Log.message("\tOutput: All hits of the AS transcripts.");
            }else{
                Log.message("\tOutput: Non-overlapped hits of each AS transcripts.");
            }

        }
        Log.message("\t# of threads : " + funkSettings.get(FunkSettings.CPU) + "\n");
    }

    /**
     *
     * @return True to perform search of domains on reference transcripts.
     */
    public static boolean isConstitutive(){
        return funkSettings.get(FunkSettings.CONSTITUTIVE);
    }

    /**
     *
     * @return True to perform search of domains on all transcripts of each gene with AS.
     */
    public static boolean isNaive(){
        return funkSettings.get(FunkSettings.NAIVE_SEARCH);
    }

    /**
     *
     * @return True to print non-overlapped hits per gene.
     */
    public static boolean isBestGeneHits(){
        return funkSettings.get(FunkSettings.OUTPUT_HITS_PER_GENE);
    }

    public static boolean isAllTranscriptHits(){
        return funkSettings.get(FunkSettings.OUTPUT_ALL_TRANSCRIPT_HITS);
    }

    public static Gene[] getGenes() {
        return genes;
    }
}