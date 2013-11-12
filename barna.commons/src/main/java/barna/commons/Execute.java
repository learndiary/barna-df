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

package barna.commons;

import barna.commons.log.Log;
import barna.commons.system.OSChecker;

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
            // NOTE: we have to make sure we have at least 3 threads in the pool. This is
            // because for BARNA-309 we ensured everything goes through the fixed thread pool.
            // In sorter operations for example in the flux, this can lock with < 3 threads, because
            // there will threads dropped into the pool in this order:
            //
            //  - one thread for main (capacitor call())
            //  - one thread for the background sorter
            //  - n threads submitted from the background sorter (per merge chunk)
            //
            // With this, we need the ability to run at least 3 jobs in parallel without blocking.
            executor = Executors.newFixedThreadPool(Math.max(maxThreads, 3));
            shutdown = new ExecutorShutdown();
            Runtime.getRuntime().addShutdownHook(shutdown);
        }else{
            Log.debug("Executor is already initialized!");
        }
    }

    /**
     * Get the executor
     *
     * @return executor the executor
     */
    public static ExecutorService getExecutor(){
        if(executor == null){
            Log.warn("The executor was not initialized properly! Make sure you call initialize first, and remember to call shutdown at the end!"+ OSChecker.NEW_LINE +
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
