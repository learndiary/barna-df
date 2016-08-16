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

package barna.astalavista;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.commons.parameters.ParameterException;
import barna.commons.parameters.ParameterSchema;
import barna.io.GeneAheadReaderThread;
import barna.io.gtf.GTFwrapper;
import barna.model.*;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.*;
import java.util.*;


public abstract class AStalavista implements Tool<Void> {

    public void setSettings(AStalavistaSettings settings) {
        this.settings = settings;
    }

    /**
     * Parameters for the AStalavista run
     */
    protected AStalavistaSettings settings= null;

    /**
     * Thread handling GTF parsing and building up gene models
     */
    protected GeneAheadReaderThread readerThread= null;

    /**
     * Tasks to be executed <b>before</b> iterating over loci
     * @throws Exception if something went wrong
     */
    protected void callBegin() throws Exception {
        // print parameters
        ArrayList<barna.commons.parameters.Parameter> pars= new ArrayList<barna.commons.parameters.Parameter>(
                settings.getParameters().values());
        Collections.sort(pars, new barna.commons.parameters.Parameter.ParameterByNameComparator());
        for (barna.commons.parameters.Parameter p : pars) {
            Log.message("# "+ p.getName()+ "\t"+ settings.get(p));
        }

        // init input / reader thread
        readerThread= getGeneAheadReaderThread();

    }

    /**
     * Tasks to carry out when solving one locus (gene).
     * @param g a gene
     * @throws Exception if something went wrong
     */
    protected void callLoop(Gene g) throws Exception {
        //g.setSpecies(species);
    }

    /**
     * Tasks to carry out before solving a batch.
     * @param g sample Gene from the batch
     * @throws Exception if something went wrong
     */
    protected void callBatchBegin(Gene g) throws Exception {
    }

    /**
     * Tasks to carry out after solving a batch.
     * @param g sample Gene from the batch
     * @throws Exception if something went wrong
     */
    protected void callBatchFinish(Gene g) throws Exception {
    }

    /**
     * Tasks to carry out <b>after</b> iteration over loci stopped.
     * @throws Exception if something went wrong
     */
    protected void callFinish() throws Exception {
        getReader().close();
    }


    @Override
    public Void call() throws Exception {


        // INIT
        long t0= System.currentTimeMillis();
        Log.message("# started\t" + new Date(t0));

        callBegin();


        // MAIN LOOP
        Log.progressStart("Iterating Annotation");
        long inBytes= settings.get(AStalavistaSettings.IN_FILE).length();
        Log.progress(0, inBytes);
        while (true) {
            Gene[] g= readerThread.getGenes();
            Log.progress(readerThread.getBytesRead(),inBytes);
            if (g== null)
                break;

            Gene ge= null;
            for (int i = 0; i < g.length; i++) {
                callLoop(g[i]);
                ge= g[i];
                g[i] = null;    // de-reference (if memory is an issue)
            }

            callBatchFinish(ge);
        }
       Log.progressFinish("done.", true);
       Log.info("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");

       callFinish();

       return null;
    }


    private GeneAheadReaderThread getGeneAheadReaderThread() {

        if (readerThread== null) {

            readerThread= new GeneAheadReaderThread(getReader());
            readerThread.setOutput(false);
            readerThread.setOutput2(true);
            readerThread.start();
        }

        return readerThread;
    }

    private GTFwrapper getReader() {
        GTFwrapper reader= new GTFwrapper(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        if (!reader.isApplicable()) {
            File f= reader.sort();
            settings.set(AStalavistaSettings.IN_FILE, f);
            reader= new GTFwrapper(settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        }
//        if (readAheadLimit> 0)
//            reader.setReadAheadLimit(readAheadLimit);
        reader.setNoIDs(null);
        return reader;
    }




    /**
     * Lazily create new instance, to be overwritten by sub-classes with
     * extended settings.
     * @return an instance of the settings associated, possibly an empty stub
     */
    protected AStalavistaSettings getSettings() {
        if (settings== null)
            return new AStalavistaSettings();
        return settings;
    }


    @Override
    public List<Parameter> getParameter() {

        // converts parameter file parameters to CLI parameters

        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        AStalavistaSettings settings= getSettings();
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter parameter : pars) {

            Class c= parameter.getType();
            Parameter p= null;
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

    @Override
    public boolean validateParameter(JSAPResult args) {
        return validateParameter(new AStalavistaSettings(), args);
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

        // output help\
        if (args != null && args.userSpecified("printParameters")) {
            // TODO
            // check why the following line does not work
            // because parameter is disabled in AStalavistaSettings
            //if (settings.get(AStalavistaSettings.HELP)) {
            settings.write(System.out);
            //}
            return false;
        }

        // non-null settings are to be assumed
        if (schema == null || !(schema instanceof AStalavistaSettings)) {
            Log.error("Must provide an instance of AStalavistaSettings!");
            return false;
        }
        settings = (AStalavistaSettings) schema;

        if (args != null)
            try {
                settings = (AStalavistaSettings) ParameterSchema.create(settings,
                        JSAPParameters.getParameterMap(settings, args));
            } catch (ParameterException e) {
                Log.error(e.getMessage(), e);
                return false;
            }


        return validateSettings(settings);
    }

    public boolean validateSettings(AStalavistaSettings settings) {

        // init values from parameter file, if any
        if (settings.get(AStalavistaSettings.PAR_FILE)!= null) {
            InputStream in = null;
            try {
                File f= settings.get(AStalavistaSettings.PAR_FILE);
                f = new File(f.getCanonicalPath());    // kill Win32ShellFolder instances, they fuck up relative path conversion
                //relativePathParser.setParentDir(f.getParentFile());
                in = new FileInputStream(f);
                settings.parse(in);
                settings.validate();
            } catch(Exception e) {
                Log.error(e.getMessage(), e);
                return false;
            } finally {
                if (in != null) {
                    try {
                        in.close();
                    } catch (IOException e) {
                    }
                }
            }
        }

        if (settings.get(AStalavistaSettings.IN_FILE)== null) {
            Log.error("Hey, you forgot to specify a valid input file!\n"+
                    "This is a bit important, I cannot work without an input annotation. I want a GTF file with\n" +
                    "transcript annotations (exon features, with a mandatory optional attribute named \'transcript_id\'\n) " +
                    "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it\n" +
                    "has to be in all lines of the file in column #10. The rest of the file should comply with the\n" +
                    "standard as specified at http://mblab.wustl.edu/GTF2.html.\n"+
                    "There may also be CDS features, but they become only interesting when checking for additional things\n" +
                    "as NMD probability etc.."
            );
            return false;
        }

        if (settings.get(AStalavistaSettings.CHR_SEQ)!= null) {

            File f= settings.get(AStalavistaSettings.CHR_SEQ);
            if (f.exists()&& f.isDirectory()) {
                Graph.overrideSequenceDirPath= f.getAbsolutePath();
            } else {
                Log.message("Trying to parse species_version pair");
                String[] s= f.getName().split("_");
                if (s.length!= 2) {
                    Log.error(f.getAbsolutePath() + " is not a valid species name");
                    return false;
                }
                Log.error("Systematic species names not supported!");
//                Species spe= new Species(s[0]);
//                spe.setGenomeVersion(s[1]);
//                AStalavista.setSpecies(spe);


            }
        }

        return true;
    }

    /**
     * Method to convert command line arguments to parameter file stub,
     * possibly overwriting values set via additional parameter file.
     *
     * @param settings pre-defined settings or <code>null</code>
     * @param args parsed command line arguments
     * @return newly created or extended settings
     * @throws ParameterException in case a parameter does not get what
     * it expects
     */
    AStalavistaSettings createSettings(AStalavistaSettings settings, JSAPResult args) throws ParameterException {

        // copy
        Collection<barna.commons.parameters.Parameter> pars=
                settings.getParameters().values();
        for (barna.commons.parameters.Parameter p : pars) {
            if (args.userSpecified(p.getLongOption())) {
                p.parse(args.getObject(p.getLongOption()).toString());
            }
        }

        return settings;
    }

    /**
     * Summarizes the settings before the run
     * @param settings the settings
     */
    @Deprecated
    protected void printSettings(AStalavistaSettings settings) {

        Log.message("# started\t" + new Date(System.currentTimeMillis()));
        Log.message("# INPUT #");
        Log.message("# annotation\t" + settings.get(AStalavistaSettings.IN_FILE).getAbsolutePath());
        Log.message("# chromosomes\t"+ barna.model.Graph.overrideSequenceDirPath);

    }

}
