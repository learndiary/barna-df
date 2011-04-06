package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.options.HelpPrinter;
import org.cyclopsgroup.jcli.ArgumentProcessor;

import java.util.concurrent.Callable;

/**
 * Base interface for flux simulator tools. Implement this interface to add new flux tools.
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
     *
     * @param printer the printer to print the actual help message
     * @param toolArguments the argument processor
     * @return valid true if valid arguments
     */
    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments);

}
