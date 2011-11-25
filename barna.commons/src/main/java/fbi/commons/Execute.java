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

package fbi.commons;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Helper class that provides a configured executor to run stuff in a background task.
 * <B>Note that you MUST explicitly initializes this class!</B> After initialization, this adds
 * a system shutdown hook to kill the executor on system shutdown, but you are
 * advised to explicitly call {@link #shutdown()}.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Execute {
    /**
     * The executor
     */
    private static ExecutorService executor;
    /**
     * The current shutdown thread
     */
    private static ExecutorShutdown shutdown;

    /**
     * Shutdown the executor
     */
    public static void shutdown(){
        if(executor != null){
            executor.shutdown();
            try {
                executor.awaitTermination(30, TimeUnit.MINUTES);
            } catch (InterruptedException ignore) {
                // ignore
            }
            executor.shutdownNow();
            executor = null;
            if(shutdown != null){
                try{
                    Runtime.getRuntime().removeShutdownHook(shutdown);
                }catch (IllegalStateException inShutdown){
                    if(!inShutdown.getMessage().equals("Shutdown in progress")){
                        throw inShutdown;
                    }
                }
            }
            shutdown = null;
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
            shutdown = new ExecutorShutdown();
            Runtime.getRuntime().addShutdownHook(shutdown);
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
            Log.warn("The executor was not initialized properly! Make sure you call initialize first, and remember to call shutdown at the end!\n" +
                    "The executor will be initialized with 2 threads now!");
            initialize(2);
            //throw new RuntimeException("The executor was not initialized properly! Make sure you call initialize first, and remember to call shutdown at the end!");
        }
        return executor;
    }

    /**
     * Shutdown thread to be added as shutdown hook
     */
    private static class ExecutorShutdown extends Thread{
        @Override
        public void run() {
            super.run();
            Execute.shutdown();
        }
    }
}
