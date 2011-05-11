package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.CommandLine;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * Flux Tool that implements the simulation pipeline
 *
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
            try {
                settings= FluxSimulatorSettings.createSettings(file);
            } catch (Exception e) {
                throw new RuntimeException("Unable to load settings from " + file +"\n\n " +e.getMessage(), e);
            }
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
            fragmenter = new Fragmenter(getSettings(), getProfiler());
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
            sequencer = new Sequencer(getSettings(), getProfiler());
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
            Log.info("No mode selected, executing the full pipeline (-x -l -s)");
            setExpression(true);
            setLibrary(true);
            setSequence(true);
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

        Log.info("I am collecting information on the run.");
        FluxSimulatorSettings settings = getSettings(); // initialize settings
        if(settings == null){
            Log.error("No settings available");
            return null;
        }

        if (settings.get(FluxSimulatorSettings.PRO_FILE).exists()) {
            getProfiler().loadStats();
            if (isExpression() && getProfiler().isFinishedExpression()) {
                if (!CommandLine.confirm("[CAUTION] I overwrite the expression values in file "+ settings.get(FluxSimulatorSettings.PRO_FILE).getName()+", please confirm:\n\t(Yes,No,Don't know)"))
                    return null;
                else {
                    boolean b= settings.get(FluxSimulatorSettings.PRO_FILE).delete();	// TODO maybe only remove rfreqs..
                    if (getProfiler()!= null)
                        getProfiler().status= -1;
                }
            }
            Log.message("");
        }


        if (settings.get(FluxSimulatorSettings.LIB_FILE) != null&& settings.get(FluxSimulatorSettings.LIB_FILE).exists()&& settings.get(FluxSimulatorSettings.LIB_FILE).canRead()) {
            if (isExpression() || isLibrary()) {
                Log.info("Removing existing fragmentation file " + settings.get(FluxSimulatorSettings.LIB_FILE).getAbsolutePath());
                settings.get(FluxSimulatorSettings.LIB_FILE).delete();
            } else{
                Log.info("Loading Fragmentation from " + settings.get(FluxSimulatorSettings.LIB_FILE));
                getFragmenter().loadStats();
            }
        }

        if (settings.get(FluxSimulatorSettings.ERR_FILE) != null&& !getSequencer().loadErrors())
            throw new RuntimeException("The sequencer produced errors !"); // todo: describe the problem

        if (settings.get(FluxSimulatorSettings.SEQ_FILE) != null&& settings.get(FluxSimulatorSettings.SEQ_FILE).exists()&& settings.get(FluxSimulatorSettings.SEQ_FILE).canRead()) {
            if (isSequence()) {
                if (!CommandLine.confirm("[ATTENTION] I am going to delete the sequencing file " + settings.get(FluxSimulatorSettings.SEQ_FILE).getName() + ", please confirm:\n\t(Yes,No,Don't know)"))
                    return null;
                settings.get(FluxSimulatorSettings.SEQ_FILE).delete();
            } else{
                Log.info("Loading sequencing file " + settings.get(FluxSimulatorSettings.SEQ_FILE).getAbsolutePath());
                getSequencer().loadStats();
            }
            Log.message("");
        }
        Log.message("");

        // now start the pipeline
        long t0= System.currentTimeMillis();

        if (isExpression())
            getProfiler().call();
        else{
            Log.info("you did not ask for expression, I skip it.\n");
        }


        if (isLibrary()){
            String message = getFragmenter().isReady();
            if (message != null) {
                throw new RuntimeException(message);
            }
            getFragmenter().run();
        }else{
            Log.info("you did not want me to construct the library, I skip it.\n");
        }
        Log.message("");


        if (isSequence()){
            String message = getSequencer().isReady();
            if (message != null) {
                throw new RuntimeException(message);
            }
            getSequencer().run();
        }else {
            Log.info("sequencing has not been demanded, skipped.\n");
        }


        Log.message("\n[END] I finished, took me " + (System.currentTimeMillis() - t0) / 1000 + " sec.");
        return null;
    }

}
