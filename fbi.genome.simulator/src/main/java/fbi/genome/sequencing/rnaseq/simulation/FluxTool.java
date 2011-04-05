package fbi.genome.sequencing.rnaseq.simulation;

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

}
