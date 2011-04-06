package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.CommandLine;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "simulator", description = "Flux Simulation Pipeline")
public class SimulationPipeline implements FluxTool<Void> {

    /**
     * Expression mode
     */
	private boolean expression;
    /**
     * Library mode
     */
    private boolean library;
    /**
     * Sequence mode
     */
    private boolean sequence;

    /**
     * The parameter file
     */
	private File file;
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
     * Returns true if expression mode os active
     *
     * @return expression true if expression mode is active
     */
    public boolean isExpression() {
        return expression;
    }

    /**
     * Activate/Deactivate expression mode
     *
     * @param expression activate/deactivate expression mode
     */
    @Option(name = "x", longName = "express", description = "simulate expression")
    public void setExpression(boolean expression) {
        this.expression = expression;
    }
    /**
     * Returns true if library mode os active
     *
     * @return library true if expression mode is active
     */
    public boolean isLibrary() {
        return library;
    }
    /**
     * Activate/Deactivate library mode
     *
     * @param library activate/deactivate library mode
     */
    @Option(name = "l", longName = "library", description = "simulate library construction")
    public void setLibrary(boolean library) {
        this.library = library;
    }
    /**
     * Returns true if sequence mode os active
     *
     * @return library true if expression mode is active
     */
    public boolean isSequence() {
        return sequence;
    }
    /**
     * Activate/Deactivate sequence mode
     *
     * @param sequence activate/deactivate sequence mode
     */
    @Option(name = "s", longName = "sequence", description = "simulate sequencing")
    public void setSequence(boolean sequence) {
        this.sequence = sequence;
    }

    /**
     * The parameter file
     *
     * @return parameterFile the parameter file
     */
    public File getFile() {
        return file;
    }

    /**
     * Set the parameter file
     *
     * @param file parameter file
     */
    @Option(name="p", longName = "parameter", description = "specify parameter file (PAR file)", displayName = "file", required = true)
    public void setFile(File file) {
        this.file = file;
    }

    /**
     * Get the flux simulator settings
     *
     * @return settings the flux simulator settings
     */
    public FluxSimulatorSettings getSettings() {
        if(settings == null){
            // init
            settings= FluxSimulatorSettings.createSettings(file);
            if (settings== null)
                throw new RuntimeException("Unable to load settings from " + file);
        }
        return settings;
    }

    /**
     * Get the profiler
     *
     * @return profiler the profiler
     */
    public Profiler getProfiler() {
        if(profiler == null){
            profiler= new Profiler(getSettings());
            getSettings().setProfiler(profiler);
        }
        return profiler;
    }

    /**
     * Get the fragmenter
     *
     * @return fragmenter the fragmenter
     */
    public Fragmenter getFragmenter() {
        if(fragmenter == null){
            fragmenter = new Fragmenter(getSettings());
        }
        return fragmenter;
    }

    /**
     * Get the sequencer
     *
     * @return sequencer the sequencer
     */
    public Sequencer getSequencer() {
        if(sequencer == null){
            sequencer = new Sequencer(getSettings());
        }
        return sequencer;
    }


    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
        if(getFile() == null){
            Log.error("");
            Log.error("No parameter file specified !");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }
        if(!getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }

        if(!isExpression() && !isLibrary() && !isSequence()){
            Log.error("");
            Log.error("You must enable at least one mode (-x -l -s)");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }
        return true;
    }

    /**
     * Execute the configured pipeline
     *
     * @return void always returns null
     * @throws Exception in case of any configuration errors
     */
    public Void call() throws Exception {
        if (file == null || !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }

        Log.info("[INIT] I am collecting information on the run.");
        FluxSimulatorSettings settings = getSettings(); // initialize settings

        if (settings.getProFile()!= null&& settings.getProFile().exists()&& settings.getProFile().canRead()) {
            getProfiler().loadStats();
            if (isExpression() && getProfiler().isFinishedExpression()) {
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



        if (settings.getFrgFile()!= null&& settings.getFrgFile().exists()&& settings.getFrgFile().canRead()) {
            if (isExpression() || isLibrary()) {
                //System.err.println("RE-USING library");
                // see doLib()
//				if (!userCLIconfirm("[WARNING] I will overwrite the library file "+settings.getFrgFile().getName()+", please confirm:\n\t(Yes,No,Don't know)"))
//					System.exit(0);
                settings.getFrgFile().delete();
            } else
                getFragmenter().loadStats();
            Log.info("");
        }

        if (settings.getErrFile()!= null&& !getSequencer().loadErrors())
            throw new RuntimeException("The sequencer produced errors !"); // todo: describe the problem
        if (settings.getSeqFile()!= null&& settings.getSeqFile().exists()&& settings.getSeqFile().canRead()) {
            if (isExpression() || isLibrary() || isSequence()) {
                if (!CommandLine.confirm("[ATTENTION] I am going to delete the sequencing file " + settings.getSeqFile().getName() + ", please confirm:\n\t(Yes,No,Don't know)"))
                    return null;
                settings.getSeqFile().delete();
            } else
                getSequencer().loadStats();
            Log.info("");
        }
        Log.info("");

        // now start the pipeline
        long t0= System.currentTimeMillis();

        if (isExpression())
            doExpr();
        else{
            Log.info("you did not ask for expression, I skip it.\n");
        }


        if (isLibrary())
            doLib();
        else{
            Log.info("you did not want me to construct the library, I skip it.");
        }
        Log.info("");


        if (isSequence())
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
        getProfiler().run();
//		if (!profiler.profile()) {
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[FATAL] Problem during expression, I exit.");
//			System.exit(-1);
//		}
    }

    void doLib() {
        if (!getFragmenter().isReady()) {
            throw new RuntimeException("[WHATSUP] I am missing parameters for performing fragmentation.");
        }
        // see run()
//		if (fragmenter.isFinished()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//			if (!userCLIconfirm("[CAREFULLY] There are files describing a constructed libary, do you want to overwrite?\n(Yes,No,Don't care)"))
//				System.exit(-1);

        getFragmenter().run();
    }

    void doSeq() {
        String message = getSequencer().isReady();
        if (message != null) {
            throw new RuntimeException(message);
        }
        getSequencer().run();
    }

}
