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

import barna.commons.Execute;
import barna.commons.RandomFactory;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.CommandLine;
import barna.commons.launcher.Flux;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.flux.simulator.fragmentation.Fragmenter;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

/**
 * Flux Tool that implements the simulation pipeline
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SimulationPipeline implements Tool<Void> {

	public static void main(String[] args) {
		Execute.initialize(2);
		
		try {

			final SimulationPipeline mySimulator= new SimulationPipeline();

            JSAP jsap= new JSAP();
            List<Parameter> parameter = mySimulator.getParameter();
            if(parameter != null){
                try{
                    for (Parameter p : parameter) {
                        jsap.registerParameter(p);
                    }
                    JSAPResult toolParameter = jsap.parse(args);
                    if (!mySimulator.validateParameter(toolParameter)){
                        Flux.printUsage(mySimulator, jsap, Flux.findTools(), null, false);
                    }
                } catch (Exception e) {
                    Log.error("Parameter error : " + e.getMessage(), e);
                    System.exit(-1);
                }

            }

		    // run
			// mySimulator.call();
			Future f= Execute.getExecutor().submit(mySimulator);
			f.get();
			
		} catch (Throwable t) {
			System.err.println(t.getMessage());
			if (t instanceof Exception)
				((Exception) t).printStackTrace();
			else if (t instanceof Error)
				((Error) t).printStackTrace();
			
		} finally {
			Execute.shutdown();
		}

	}
	
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


    @Override
    public String getName() {
        return "simulator";
    }

    @Override
    public String getDescription() {
        return "The Flux Simulator";
    }

    @Override
    public String getLongDescription() {
        return null;
    }


    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("parameter", 'p').type(File.class).help("specify parameter file (PAR file)").valueName("file").get());
        parameters.add(JSAPParameters.switchParameter("express", 'x').help("Simulate Expression").get());
        parameters.add(JSAPParameters.switchParameter("library", 'l').help("Simulate Library Construction").get());
        parameters.add(JSAPParameters.switchParameter("sequence", 's').help("Simulate Sequencing").get());
        parameters.add(JSAPParameters.switchParameter("printParameters", 'o').help("Print default parameters").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setFile(args.getFile("parameter"));
        setPrintParameters(args.userSpecified("printParameters"));
        setLibrary(args.userSpecified("library"));
        setExpression(args.userSpecified("express"));
        setSequence(args.userSpecified("sequence"));

        if(isPrintParameters()){
            FluxSimulatorSettings settings = new FluxSimulatorSettings();
            settings.write(System.out);
            return true;
        }

        if (getFile() == null) {
            Log.error("");
            Log.error("No parameter file specified !");
            Log.error(barna.commons.system.OSChecker.NEW_LINE);
            return false;
        }
        if (!getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error(barna.commons.system.OSChecker.NEW_LINE);
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
        if(isPrintParameters()){
            return null;
        }

        if (file == null || !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }

        Log.info("I am collecting information on the run.");
        FluxSimulatorSettings settings = getSettings(); // initialize settings
        if (settings == null) {
            Log.error("No settings available");
            return null;
        }
        // BARNA-306 initialize the global seed
        if(settings.get(FluxSimulatorSettings.SEED) != 0){
            RandomFactory.SEED = settings.get(FluxSimulatorSettings.SEED);
            Log.info("Random seed set to : " + RandomFactory.SEED);
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
        if(sortedFileName.toLowerCase().endsWith(".gz")){
            sortedFileName = sortedFileName.substring(0, sortedFileName.length()-3);
        }
        File sorted = new File(sortedFileName);

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
