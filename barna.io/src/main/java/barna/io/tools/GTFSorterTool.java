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

import barna.commons.launcher.FluxTool;
import barna.commons.launcher.HelpPrinter;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.io.gtf.GTFwrapper;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * Sort GTF Files from command line
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "sortGTF", description = "Sort a GTF file. If no output file is specified, result is printed to standard out")
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
    @Option(name = "i", longName = "input", description = "GTF input file", required = true)
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
    @Option(name = "o", longName = "output", description = "GTF output file. Sorts to stdout if no file is given", required = false)
    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

    @Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        if (getInFile() == null) {
            printer.out.println("Please specify an input file");
            return false;
        }
        return true;
    }

    @Override
    public Object call() throws Exception {
        GTFwrapper w= new GTFwrapper(inFile);
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
