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

package barna.commons.launcher;

import org.cyclopsgroup.jcli.ArgumentProcessor;

import java.util.concurrent.Callable;

/**
 * Base interface for flux tools. Implement this interface to add new flux tool.
 * If the default main class is used, the tools are automatically registered. You can use the
 * {@code @Cli} annotation on the class to give it a name and a description and then annotate
 * setter methods with {@code @Option} to add them as command line parameters.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FluxTool<T> extends Callable<T> {

    /**
     * This method is called after the command line arguments are processed. Implementations
     * should validate the parameters. If a required parameter is not set, return false and use
     * the printer to print information
     *
     * @param printer       the printer to print the actual help message
     * @param toolArguments the argument processor
     * @return valid true if valid arguments
     */
    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments);

}
