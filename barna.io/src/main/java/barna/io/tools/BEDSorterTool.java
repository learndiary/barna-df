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

package barna.io.tools;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.io.bed.BEDwrapper;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Sort BED Files from command line
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BEDSorterTool implements FluxTool {
    /**
     * The source file
     */
    private File inFile;
    /**
     * The output file
     */
    private File outFile;

    /**
     * Get the input file
     *
     * @return input file
     */
    public File getInFile() {
        return inFile;
    }

    /**
     * Set the input file
     *
     * @param inFile the input file
     */
    public void setInFile(final File inFile) {
        this.inFile = inFile;
    }

    /**
     * Get the output file
     *
     * @return output file
     */
    public File getOutFile() {
        return outFile;
    }

    /**
     * Set the output file
     *
     * @param outFile the output file
     */

    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

    @Override
    public String getName() {
        return "sortBED";
    }

    @Override
    public String getDescription() {
        return "Sort a BED file. If no output file is specified, result is printed to standard out";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("BED input file").valueName("bed").required().get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("BED output file. Sorts to stdout if no file is given").valueName("bed").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setInFile(args.getFile("input"));
        if(args.userSpecified("output"))
            setOutFile(args.getFile("output"));

        if (getInFile() == null || !getInFile().canRead()) {
            Log.error("Unable to read input file " + getInFile().getAbsolutePath());
            return false;
        }
        return true;


    }

    @Override
    public Object call() throws Exception {
        BEDwrapper w= new BEDwrapper(inFile);
        if(getOutFile() != null){
            Log.info("SORT", "Sorting " + getInFile().getAbsolutePath() +" to " + getOutFile().getAbsolutePath());
            Log.progressStart("Sorting " + getInFile().getName());
            w.sort(outFile);
            Log.progressFinish(StringUtils.OK, true);
        }else{
            w.sort(System.out);
        }
        return null;
    }
}
