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
