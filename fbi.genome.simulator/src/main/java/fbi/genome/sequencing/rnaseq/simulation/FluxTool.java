package fbi.genome.sequencing.rnaseq.simulation;

import java.util.concurrent.Callable;

/**
 * Base interface for flux simulator tools
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FluxTool<T> extends Callable<T> {

}
