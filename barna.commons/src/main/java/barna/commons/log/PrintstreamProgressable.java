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

package barna.commons.log;

import barna.commons.Progressable;

import java.io.PrintStream;


/**
 * Print progress to a stream, i.e. System.out or System.err
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class PrintstreamProgressable implements Progressable {
    /**
     * The progress character
     */
    private char progressChar;
    /**
     * The stream
     */
    private PrintStream stream;
    /**
     * The maximum value
     */
    private int steps = 10;
    /**
     * Last printed value
     */
    private int currentStep = 0;

    /**
     * The start time of the current progress
     */
    private long startTime;


    /**
     * Create a new progress printer. The default minimum value is 0 and the maximum is 9.
     *
     * @param stream the stream
     */
    public PrintstreamProgressable(PrintStream stream) {
        this('*', stream);
    }

    /**
     * Create a new instance using the given character and stream. The minimum value is 0 and the maximum
     * is 9.
     *
     * @param progressChar the character
     * @param stream       the stream
     */
    public PrintstreamProgressable(char progressChar, PrintStream stream) {
        this(progressChar, stream, 10);
    }

    /**
     * Create a new instance using the given character to print progress. The progress is printed to
     * the given stream and is limited by the given min and max values
     *
     * @param progressChar the progress character
     * @param stream       the stream to print to
     * @param steps        the maximum value
     */
    public PrintstreamProgressable(char progressChar, PrintStream stream, int steps) {
        if (stream == null) {
            throw new NullPointerException("Null stream not permitted");
        }
        this.progressChar = progressChar;
        this.stream = stream;
        this.steps = steps;
    }

    public void start(String message) {
        startTime = System.currentTimeMillis();
        if (message != null) {
            stream.print("\t" + message + " ");
        } else {
            stream.print("\t ");
        }
        stream.flush();
        currentStep = 0;
    }

    public void progress() {
        if (currentStep < steps) {
            stream.print(progressChar);
            stream.flush();
            ++currentStep;
        }
    }

    public void finish() {
        /*
         * Print the progress if progress already started
         */
        if (currentStep > 0) {
            for (int i = currentStep; i < steps; i++) {
                stream.print("*");
            }

            currentStep = 0;
        }
        /*
         * new line
         */
        stream.println();
        stream.flush();
    }


    public void finish(String msg, boolean printTime) {

        long time = 0;
        if (printTime) {
            time = System.currentTimeMillis() - startTime;
        }
        if (currentStep > 0) {
            for (int i = currentStep; i < steps; i++) {
                stream.print(progressChar);
            }

            currentStep = 0;
        }

        if (msg != null) {
            stream.print(" " + msg);
        }

        if (time >= 0) {
            int hh = (int) (time / 3600000);
            if (hh > 0) {
                time %= (hh * 3600000);
            }
            int mm = (int) (time / 60000);
            if (mm > 0) {
                time %= (mm * 60000);
            }
            int ss = (int) (time / 1000);
            if (ss > 0) {
                time %= (ss * 1000);
            }
            String s = hh < 10 ? "0" + Integer.toString(hh) : Integer.toString(hh);
            s += ":";
            s += mm < 10 ? "0" + Integer.toString(mm) : Integer.toString(mm);
            s += ":";
            s += ss < 10 ? "0" + Long.toString(ss) : Long.toString(ss);
            stream.print(" (" + s + ")");
        }
        stream.println();
        stream.flush();
    }

    public void failed(String msg) {
        if (msg != null) {
            stream.print(" " + msg);
        }
        currentStep = 0;
        stream.println();
        stream.flush();
    }

    public int steps() {
        return steps;
    }

    public int currentStep() {
        return currentStep;
    }
}
