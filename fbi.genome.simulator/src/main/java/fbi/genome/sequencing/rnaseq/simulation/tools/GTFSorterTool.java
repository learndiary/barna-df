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

package fbi.genome.sequencing.rnaseq.simulation.tools;

import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.genome.io.gtf.GTFwrapper;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * Sort GTF Files
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "sortGtf", description = "Sort a GTF file")
public class GTFSorterTool implements FluxTool {

    private File intFile;
    private File outFile;

    public File getIntFile() {
        return intFile;
    }

    @Option(name = "i", longName = "input", description = "GTF input file")
    public void setIntFile(final File intFile) {
        this.intFile = intFile;
    }

    public File getOutFile() {
        return outFile;
    }

    @Option(name = "o", longName = "output", description = "GTF output file")
    public void setOutFile(final File outFile) {
        this.outFile = outFile;
    }

    @Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        if (getIntFile() == null) {
            printer.out.println("Please specify an input file");
            return false;
        }

        if (getOutFile() == null) {
            printer.out.println("Please specify an output file");
            return false;
        }
        return true;
    }

    @Override
    public Object call() throws Exception {
        Log.progressStart("Sorting " + getIntFile().getName());
        GTFwrapper w= new GTFwrapper(intFile);
        w.sort(outFile);
        Log.progressFinish(StringUtils.OK, true);
        return null;
    }
}
