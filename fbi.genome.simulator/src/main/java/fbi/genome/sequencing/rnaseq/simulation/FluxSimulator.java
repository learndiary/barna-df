package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.tools.CommandLine;
import fbi.genome.io.SpliceGraphIO;
import fbi.genome.model.IntronModel;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.net.URL;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.jar.Attributes;
import java.util.jar.JarFile;
import java.util.jar.Manifest;

//import gphase.solexa.lp.GraphLPsolver5;
//import gphase.solexa.simulation.Nebulizer;

// TODO check whether selected fragments concord w pro file

public class FluxSimulator implements Callable<Void> {
    /**
     * Current Flux Simulator version
     */
    public static String FLUX_VERSION = "";
    /**
     * Current Flux Simulator revision
     */
    public static String FLUX_REVISION = "";

    /**
     * Default 5' flank
     */
    private static final int FLANK5_LEN= 50;
    /**
     * Default 3' flank
     */
    private static final int FLANK3_LEN= 50;

    /**
     * Expression mode
     */
	private boolean expression = false;
    /**
     * Library mode
     */
    private boolean library = false;
    /**
     * Sequence mode
     */
    private boolean sequence = false;
    /**
     * Extract splice junctions mode
     */
    private boolean extractSpliceJunctions = false;

    /**
     * The parameter file
     */
	private File file;
    /**
     * IModel file
     */
    private File modelFile;
    /**
     * Genome directory
     */
    private File genomeFile;
    /**
     * The flux settings
     */
	private FluxSimulatorSettings settings;
    /**
     * The Profiler
     */
	private Profiler profiler;
    /**
     * The Fragmenter
     */
	private Fragmenter fragmenter;
    /**
     * The sequencer
     */
	private Sequencer sequencer;
    /**
     * Store the flanks
     */
	private int[] eflanks= null;


    /**
     * Create a new simulator
     */
    public FluxSimulator() {
    }


    /**
     * Set the input file
     *
     * @param file the file
     */
    public void setFile(File file) {
        this.file= file;
    }

    /**
     * Set the model file
     *
     * @param modelFile absolute path to the model file
     */
    public void setImodel(File modelFile) {
        this.modelFile = modelFile;
    }

    /**
     * Activate splice junction extraction and use the given file as input
     *
     * @param extractSpliceJunctions activate/deactivate
     * @param file input file
     */
    public void setExtractSpliceJunctions(boolean extractSpliceJunctions, File file) {
        this.extractSpliceJunctions = extractSpliceJunctions;
        if(extractSpliceJunctions) setFile(file);
    }

    /**
     *
     * @return expressionMode true if expression mode is on
     */
    public boolean isExpression() {
        return expression;
    }

    /**
     * Activate expression mode
     *
     * @param expression activate/deactivate expression mode
     */
    public void setExpression(boolean expression) {
        this.expression = expression;
    }

    /**
     *
     * @return libraryMode true if library mode is active
     */
    public boolean isLibrary() {
        return library;
    }

    /**
     * Activate/Deactivate library mode
     *
     * @param library activate/deactivate library mode
     */
    public void setLibrary(boolean library) {
        this.library = library;
    }

    /**
     *
     * @return sequenceMode true if sequence mode is active
     */
    public boolean isSequence() {
        return sequence;
    }

    /**
     * Activate/deactivate sequence mode
     * @param sequence activate/deactivate sequence mode
     */
    public void setSequence(boolean sequence) {
        this.sequence = sequence;
    }

    /**
     * Get the current flanks. Returns a two element array, index 0 is 5' flank,  index 1 is the 3' flank
     *
     * @return flanks the flanks
     */
    public int[] getEFlanks() {
        if (eflanks == null) {
            eflanks = new int[2];
            eflanks[0]= -1;
            eflanks[1]= -1;
        }

        return eflanks;
    }

    /**
     * Set the genome directory
     *
     * @param genomeDir absolute path to the genome directory
     */
    public void setGenomeDir(File genomeDir) {
        if(genomeDir != null && !genomeDir.isDirectory()) throw new IllegalArgumentException(genomeDir + " is no directory!");
        this.genomeFile = genomeDir;
    }

    /**
     * Set the 5' flank
     *
     * @param length flank length  @code{> 0}
     */
    public void set5flank(int length) {
        getEFlanks()[0]= -1;
        if(length <= 0){
            throw new IllegalArgumentException("Not a valid length for 5' exon flank: " + length);
        }
        getEFlanks()[0]= length;
    }

    /**
     * Set the 3' flank
     *
     * @param length flank length @code{> 0}
     */
    public void set3flank(int length) {
        getEFlanks()[1]= -1;
        if(length <= 0){
            throw new IllegalArgumentException("Not a valid length for 3' exon flank: " + length);
        }
        getEFlanks()[1]= length;
    }


    /**
     * Executes the pipeline
     *
     * @return void always returns null
     * @throws Exception in case of execution or configuration errors
     */
    public Void call() throws Exception {
        Log.info("I am the Flux Simulator (v"+ FLUX_VERSION +" build"+FLUX_REVISION+"), nice to meet you!\n");
        if (file == null || !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }

        if (extractSpliceJunctions) {
            if (genomeFile == null|| !genomeFile.exists()) {
                throw new RuntimeException("[AIAIII] I have no directory with the genomic sequences!");
            } else{
                // todo : refactor this to not use a static variable
                fbi.genome.model.Graph.overrideSequenceDirPath= genomeFile.getAbsolutePath();
            }

            int[] flanks= getEFlanks();
            flanks[0]= flanks[0]< 0? FLANK5_LEN: flanks[0];
            flanks[1]= flanks[1]< 0? FLANK3_LEN: flanks[1];

            IntronModel iModel= new IntronModel();
            if (modelFile != null){
                Log.info("Reading Model file " + modelFile.getAbsolutePath());
                iModel.read(modelFile);
            }
            Log.info("Extracting splice junctions, 5'sequence "+flanks[0]
                + ", 3'sequence "+eflanks[1]+", intron model "+ iModel.getName());
            SpliceGraphIO.extractSpliceJunctions(flanks[0], flanks[1], iModel, file, null);
            // todo: exit here ?
        }


        // init
        settings= FluxSimulatorSettings.createSettings(file);
        if (settings== null)
            throw new RuntimeException("Unable to create Setting file!");

        Log.info("[INIT] I am collecting information on the run.");
        profiler= new Profiler(settings);
        settings.setProfiler(profiler);
        if (settings.getProFile()!= null&& settings.getProFile().exists()&& settings.getProFile().canRead()) {
            profiler.loadStats();
            if (expression && profiler.isFinishedExpression()) {
                if (!CommandLine.confirm("[CAUTION] I overwrite the expression values in file "+settings.getProFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
                    return null;
                else {
                    boolean b= settings.getProFile().delete();	// TODO maybe only remove rfreqs..
                    if (settings.getProfiler()!= null)
                        settings.getProfiler().status= -1;
                }
            }
            Log.info("");
        }
        fragmenter= new Fragmenter(settings);
        if (settings.getFrgFile()!= null&& settings.getFrgFile().exists()&& settings.getFrgFile().canRead()) {
            if (expression || library) {
                //System.err.println("RE-USING library");
                // see doLib()
//				if (!userCLIconfirm("[WARNING] I will overwrite the library file "+settings.getFrgFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
//					System.exit(0);
                settings.getFrgFile().delete();
            } else
                fragmenter.loadStats();
            Log.info("");
        }
        sequencer= new Sequencer(settings);
        if (settings.getErrFile()!= null&& !sequencer.loadErrors())
            throw new RuntimeException(""); // todo: describe the problem
        if (settings.getSeqFile()!= null&& settings.getSeqFile().exists()&& settings.getSeqFile().canRead()) {
            if (expression || library || sequence) {
                if (!CommandLine.confirm("[ATTENTION] I am going to delete the sequencing file "+settings.getSeqFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
                    return null;
                settings.getSeqFile().delete();
            } else
                sequencer.loadStats();
            Log.info("");
        }
        Log.info("");

        // do
        long t0= System.currentTimeMillis();


        if (expression)
            doExpr();
        else{
            Log.info("you did not ask for expression, I skip it.\n");
        }


        if (library)
            doLib();
        else{
            Log.info("you did not want me to construct the library, I skip it.");
        }
        Log.info("");


        if (sequence)
            doSeq();
        else {
            Log.info("sequencing has not been demanded, skipped.");
        }


        Log.info("\n[END] I finished, took me " + (System.currentTimeMillis() - t0) / 1000 + " sec.");
        return null;
    }

    void doExpr() {
//		if (!profiler.isFinishedReadAnnotation())
//			profiler.readAnnotation();
//
//		if (!profiler.isReady()) {
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[MISSING] I lack parameters for performing expression.");
//			System.exit(-1);
//		}
        profiler.run();
//		if (!profiler.profile()) {
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[FATAL] Problem during expression, I exit.");
//			System.exit(-1);
//		}
    }

    void doLib() {
        if (!fragmenter.isReady()) {
            throw new RuntimeException("[WHATSUP] I am missing parameters for performing fragmentation.");
        }
        // see run()
/*		if (fragmenter.isFinished()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			if (!userCLIconfirm("[CAREFULLY] There are files describing a constructed libary, do you want to overwrite?\n(Yes,No,Don't care)"))
				System.exit(-1);
*/
        fragmenter.run();
    }

    void doSeq() {
        if (!sequencer.isReady()) {
            throw new RuntimeException("[NONONO] I am missing parameters for sequencing.");
        }
        sequencer.run();
    }








	public static void exit(int code) {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			
		}
	}


    /**
     * Read properties like version and build revision from jar file
     *
     * @return valid returns true if properties are valid and everything is fine
     */
	private static boolean readProperties() {
        /*
        Find the manifest file and extract version revision adn jdk information
         */
        URL location = FluxSimulator.class.getResource("FluxSimulator.class");
        String fileString = location.toExternalForm();
        File jar = null;
        if(fileString.startsWith("jar")){
            fileString = fileString.substring(9);
            fileString = fileString.substring(0, fileString.lastIndexOf("!"));
            jar = new File(fileString);
        }
        String buildVersion = "";
        String buildRevision = "";
        String buildJDK = "";
        if(jar != null){
            try {
                JarFile jf = new JarFile(jar);
                Manifest manifest = jf.getManifest();
                Attributes mainAttributes = manifest.getMainAttributes();
                Object v = mainAttributes.getValue("Build-Version");
                Object r = mainAttributes.getValue("Build-Revision");
                Object j = mainAttributes.getValue("Flux-JDK");
                if(v != null){
                    buildVersion = v.toString();
                }
                if(r != null){
                    buildRevision = r.toString();
                }
                if(j != null){
                    buildJDK = j.toString();
                }
            } catch (IOException e) {}
        }


		if(!buildVersion.isEmpty()){
            FLUX_VERSION = buildVersion;
        }

        if(!buildRevision.isEmpty()){
            FLUX_REVISION = buildRevision;
        }

        if (!buildJDK.isEmpty()) {
			try {
				float v= Float.parseFloat(buildJDK);
				String ver= System.getProperty("java.version");
				int p= ver.indexOf('.', 0);
				p= ver.indexOf('.', p+1);
				float v2= Float.parseFloat(ver.substring(0, p));
				if (v2< v) {
					Log.error("Wrong java version, I need "+v+" but I found "+v2+".");
                    return false;
				}
			} catch (Exception e) {
				; // :)
			}
			
		}
        return true;
		
	}
	
	public static void main(String[] args) {
		
		if(!readProperties()){
            System.exit(-1);
        }
		// prepare options
        FluxOptions options = new FluxOptions();
        try {
            if(!options.parseParameters(args)){
                System.exit(-1);
            }
        } catch (Exception e) {
            Log.error("Error while parsing command line parameters: " + e.getMessage());
            System.exit(-1);
        }


        FluxSimulator sim= new FluxSimulator();
        options.apply(sim);
        try {
            sim.call();
        } catch (Exception e) {
            Log.error(e.getMessage());
        }
    }
	
	
	public static boolean invertTable(File invFile) {
		
		System.err.println("[INFO] inverting .pro file");
		File tmpFile= new File(invFile.getAbsolutePath()+"_inv");
		try {
			Vector<StringTokenizer> lineTokis= new Vector<StringTokenizer>();
			BufferedReader buffy= new BufferedReader(new FileReader(invFile));
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) 
				lineTokis.add(new StringTokenizer(s));
			buffy.close();
			
			int c= -1;
			for (int i = 0; i < lineTokis.size(); i++) {
				if (c< 0)
					c= lineTokis.elementAt(i).countTokens();
				else if (c!= lineTokis.elementAt(i).countTokens()) {
					System.err.println("\t[OHNO] inconsistent column count "+c+" <> "+lineTokis.elementAt(i).countTokens());
					return false;
				}
			}
			
			BufferedWriter writer= new BufferedWriter(new FileWriter(tmpFile));
			while (true) {
				for (int i = 0; i < lineTokis.size(); i++) {
					writer.write(lineTokis.elementAt(i).nextToken());
					if (i< lineTokis.size()-1)
						writer.write("\t");
					else
						writer.write("\n");
				}
				if (!lineTokis.elementAt(0).hasMoreTokens()) 
					break;					
			}			
			writer.flush();
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
			if (tmpFile.exists())
				tmpFile.delete();
			return false;
		}
		
		if (!invFile.delete()) 
			System.err.println("\t[OHNO] failed to remove "+invFile);
		if (!tmpFile.renameTo(invFile)) 
			System.err.println("\t[OHNO] failed to move "+tmpFile+" to "+invFile);
		return true;
	}


}
