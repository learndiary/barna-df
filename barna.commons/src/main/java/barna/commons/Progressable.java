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

/**
 * Step based progress. Use {@link #steps()}  to find the number of steps provided
 * by the implementation. Call {@link #start(String)} to reset the current steps and start
 * a new progress. If a current progress is running, it is canceled without a message.
 *
 * @author Thasso Griebel (Thasso.griebel@googlemail.com)
 */
public interface Progressable {
    /**
     * Start a new progress run. This ends running progresses and resets the step counter.
     * Implementation must also store the start time to be able to print the start time
     * when the progress is finished
     *
     * @param message the message (null permitted)
     */
    public void start(String message);

    /**
     * Do one progress step if more steps are available
     */
    public void progress();

    /**
     * Finish the progress. This print the rest of the available steps and then ends the progress
     */
    public void finish();

    /**
     * Finish the progress and print the optional message
     *
     * @param msg  the message (null permitted)
     * @param time print the running time
     */
    public void finish(String msg, boolean time);

    /**
     * Progress process failed. This ends the progress and prints the (optional) message
     *
     * @param msg the message (null permitted)
     */
    public void failed(String msg);

    /**
     * Returns the number of steps provided by this progressable.
     *
     * @return steps the number of steps
     */
    public int steps();

    /**
     * Returns the current step
     *
     * @return step the current step
     */
    public int currentStep();
}
