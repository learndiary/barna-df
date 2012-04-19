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
import barna.io.gtf.GTFwrapper;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Sort GTF Files from command line
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GTFSorterTool implements FluxTool {
    /**
     * The source file
     */
    private File inFile;
    /**
     * The output file
     */
    private File outFile;

    /**
     * Check if the file is sorted first
     */
    private boolean check;

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

    /**
     * Returns true if file is checked before sorting
     *
     * @return check check befor esorting
     */
    public boolean isCheck() {
        return check;
    }

    /**
     * Check if the file is already sorted before sorting
     *
     * @param check check if its sorted before sorting
     */
    public void setCheck(boolean check) {
        this.check = check;
    }

    @Override
    public String getName() {
        return "sortGTF";
    }

    @Override
    public String getDescription() {
        return "Sort a GTF file. If no output file is specified, result is printed to standard out";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("GTF input file").valueName("gtf").required().get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("GTF output file. Sorts to stdout if no file is given").valueName("output").get());
        parameters.add(JSAPParameters.switchParameter("check", 'c').help("Check if the file is sorted before sorting").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setInFile(args.getFile("input"));
        if(args.userSpecified("output"))
            setOutFile(args.getFile("output"));
        if(args.userSpecified("check")) setCheck(true);

        if (getInFile() == null || !getInFile().canRead()) {
            Log.error("Unable to read input file " + getInFile().getAbsolutePath());
            return false;
        }
        return true;
    }

    @Override
    public Object call() throws Exception {
        GTFwrapper w= new GTFwrapper(inFile);
        if(isCheck()){
            if(w.isApplicable()){                
                if(getOutFile() != null){
                    Log.info("File is already sorted, skip sorting");
                }
                return null;
            }
        }
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
