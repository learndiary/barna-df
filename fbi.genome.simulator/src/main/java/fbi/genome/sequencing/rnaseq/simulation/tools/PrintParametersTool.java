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
import fbi.genome.io.gff.GFFSorter;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.*;

/**
 * Print parameters to standard output.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name = "printParameters", description = "Print available parameters")
public class PrintParametersTool implements FluxTool {

    @Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        return true;
    }

    @Override
    public Object call() throws Exception {
        FluxSimulatorSettings settings = new FluxSimulatorSettings();
        settings.write(System.out);
        return null;
    }
}
