package fbi.commons;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Helper class that provides a configured executor to run stuff in a background task.
 * Note that you MUST explicitly initializes AND shutdown to use this class ! If you
 * do not shutdown, your process will not terminate.
 *
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Execute {
    /**
     * The executor
     */
    private static ExecutorService executor;

    /**
     * Shutdown the executor
     */
    public static void shutdown(){
        if(executor != null){
            executor.shutdownNow();
            executor = null;
        }
    }

    /**
     * Initialize the executor
     *
     * @param maxThreads number of threads
     */
    public static void initialize(int maxThreads){
        if(maxThreads <= 0) throw new IllegalArgumentException("Max threads must be >= 1");
        if(executor == null){
            executor = Executors.newFixedThreadPool(maxThreads);
        }else{
            throw new RuntimeException("Executor is already initialized!");
        }
    }

    /**
     * Get the executor
     *
     * @return executor the executor
     */
    public static ExecutorService getExecutor(){
        if(executor == null){
            throw new RuntimeException("The executor was not initialized properly! Make sure you call initialize first, and remember to call shutdown at the end!");
        }
        return executor;
    }
}
