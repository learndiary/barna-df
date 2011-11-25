/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.commons.tools.CommandLine;
import fbi.genome.io.FileHelper;
import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.sequencing.rnaseq.simulation.fragmentation.Fragmenter;
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
     * Print default parameters
     */
    private boolean printParameters;

    /**
     * True if default parameters are printed
     *
     * @return print print parameters
     */
    public boolean isPrintParameters() {
        return printParameters;
    }

    /**
     * Enable default parameter printing
     *
     * @param printParameters enable disable
     */
    @Option(name = "o", longName = "printParameters", description = "Print default parameters", required = false)
    public void setPrintParameters(final boolean printParameters) {
        this.printParameters = printParameters;
    }

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
    @Option(name = "p", longName = "parameter", description = "specify parameter file (PAR file)", displayName = "file", required = true)
    public void setFile(File file) {
        this.file = file;
    }

    /**
     * Get the flux simulator settings
     *
     * @return settings the flux simulator settings
     */
    public FluxSimulatorSettings getSettings() {
        if (settings == null) {
            // init
            try {
                settings = FluxSimulatorSettings.createSettings(file);
            } catch (Exception e) {
                throw new RuntimeException("Unable to load settings from " + file + "\n\n " + e.getMessage(), e);
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
        if (profiler == null) {
            profiler = new Profiler(getSettings());
        }
        return profiler;
    }

    /**
     * Get the fragmenter
     *
     * @return fragmenter the fragmenter
     */
    public Fragmenter getFragmenter() {
        if (fragmenter == null) {
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
        if (sequencer == null) {
            sequencer = new Sequencer(getSettings(), getProfiler());
        }
        return sequencer;
    }


    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
        if(isPrintParameters()){
            FluxSimulatorSettings settings = new FluxSimulatorSettings();
            settings.write(System.out);
            return false;
        }

        if (getFile() == null) {
            Log.error("");
            Log.error("No parameter file specified !");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }
        if (!getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }

        if (!isExpression() && !isLibrary() && !isSequence()) {
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
        if (settings == null) {
            Log.error("No settings available");
            return null;
        }

        // fix issue #60 and transfer the temp file
        File settingsTmp = settings.get(FluxSimulatorSettings.TMP_DIR);
        if(settingsTmp != null){
            FileHelper.tempDirectory = settingsTmp;
        }


        File profilerFile = settings.get(FluxSimulatorSettings.PRO_FILE);
        if (profilerFile.exists()) {
            // initialize the profiler
            if(!getProfiler().initializeProfiler(profilerFile)){
                throw new RuntimeException("Error while initializing Profiler!");
            }
            if (isExpression() && getProfiler().isFinishedExpression()) {
                if (!CommandLine.confirm("[CAUTION] I overwrite the expression values in file " + profilerFile.getName() + ", please confirm:\n\t(Yes,No,Don't know)")) {
                    return null;
                } else {
                    getProfiler().resetProfile();
                }
            }
            Log.message("");
        }


        if (settings.get(FluxSimulatorSettings.LIB_FILE) != null && settings.get(FluxSimulatorSettings.LIB_FILE).exists() && settings.get(FluxSimulatorSettings.LIB_FILE).canRead()) {
            if (isExpression() || isLibrary()) {
                Log.info("Removing existing fragmentation file " + settings.get(FluxSimulatorSettings.LIB_FILE).getAbsolutePath());
                settings.get(FluxSimulatorSettings.LIB_FILE).delete();
            } else {
                Log.info("Loading Fragmentation from " + settings.get(FluxSimulatorSettings.LIB_FILE));
                getFragmenter().loadStats(settings.get(FluxSimulatorSettings.LIB_FILE));
            }
        }

        if (settings.get(FluxSimulatorSettings.ERR_FILE) != null && settings.get(FluxSimulatorSettings.ERR_FILE).length() > 0 && !getSequencer().loadErrors()) {
            throw new RuntimeException("Unable to load the error model! Specify a valid model or disable fasta/fastq output");
        }

        if (settings.get(FluxSimulatorSettings.SEQ_FILE) != null && settings.get(FluxSimulatorSettings.SEQ_FILE).exists() && settings.get(FluxSimulatorSettings.SEQ_FILE).canRead()) {
            if (isSequence()) {
                if (!CommandLine.confirm("[ATTENTION] I am going to delete the sequencing file " + settings.get(FluxSimulatorSettings.SEQ_FILE).getName() + ", please confirm:\n\t(Yes,No,Don't know)")) {
                    return null;
                }
                settings.get(FluxSimulatorSettings.SEQ_FILE).delete();
            }
            Log.message("");
        }
        Log.message("");

        // sort the GTF file
        sortGTFReference();

        // now start the pipeline
        long t0 = System.currentTimeMillis();

        if (isExpression()) {
            getProfiler().call();
        } else {
            Log.info("you did not ask for expression, I skip it.\n");
        }


        if (isLibrary()) {
            String message = getFragmenter().isReady();
            if (message != null) {
                throw new RuntimeException(message);
            }
            getFragmenter().call();
        } else {
            Log.info("you did not want me to construct the library, I skip it.\n");
        }
        Log.message("");


        if (isSequence()) {
            String message = getSequencer().isReady();
            if (message != null) {
                throw new RuntimeException(message);
            }
            getSequencer().call();
        } else {
            Log.info("sequencing has not been demanded, skipped.\n");
        }


        Log.message("\n[END] I finished, took me " + (System.currentTimeMillis() - t0) / 1000 + " sec.");
        return null;
    }


    /**
     * Ensure that we are working on a sorted GTF file. If the file is not sorted,
     * it will be and the file in the settings is replaced and the user is informed about
     * the change.
     * <p>
     * We also check if a sorted file exists, following the naming schema, {@code <originalName>_sorted.<extension>}
     * </p>
     */
    protected void sortGTFReference() {
        File refFile = settings.get(FluxSimulatorSettings.REF_FILE);
        // Fix Issue #58 and make sure the sorted file is used and differs in name
        // Fix Simulator-8 and make sure we use a path relative to the parameters file
        String sortedFileName = settings.getParameterFile().getParent() + File.separator + FileHelper.append(refFile.getName(), "_sorted");
        File sorted = new File(sortedFileName);
        // the sorted file will not be gzipped
        if(sortedFileName.toLowerCase().endsWith(".gz")){
            sorted = new File(sortedFileName.substring(0, sortedFileName.length()-3));
        }


        GTFwrapper gffReader = new GTFwrapper(refFile.getAbsolutePath());
        // make sure the gtf is valid and sorted
        Log.info("Checking GTF file");
        if (!gffReader.isApplicable()) {
            gffReader.close();

            // okey its not sorted, check if there is a sorted version

            if (sorted.exists()) {
                gffReader = new GTFwrapper(sortedFileName);
                if (gffReader.isApplicable()) {
                    // found a sorted file
                    // inform the user and switch
                    Log.warn("GTF FILE", "The GTF reference file given is not sorted, but we found a sorted version.");
                    Log.warn("GTF FILE", "The Simulator will use " + sortedFileName);
                    Log.warn("GTF FILE", "You might want to update your parameters file");
                    settings.setRefFile(sorted);
                }
                gffReader.close();
                return;
            }
            // sort the file
            Log.warn("GTF FILE", "The GTF reference file given is not sorted, sorting it right now...");
            gffReader = new GTFwrapper(refFile.getAbsolutePath());
            gffReader.setStars(true);
            gffReader.setSilent(false);
            gffReader.sort(sorted);
            gffReader.close();
            settings.setRefFile(sorted);

            Log.warn("GTF FILE", "The Simulator will use " + sortedFileName);
            Log.warn("GTF FILE", "You might want to update your parameters file");
        }
    }


}
