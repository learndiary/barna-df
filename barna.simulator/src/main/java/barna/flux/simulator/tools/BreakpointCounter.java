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

package barna.flux.simulator.tools;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Count and print breakpoint distribution
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

public class BreakpointCounter implements FluxTool {

    private File bedFile;
    private int bins = 20;
    private boolean bounds;
    private int minLength = 0;
    private int maxLength = Integer.MAX_VALUE;



    public File getBedFile() {
        return bedFile;
    }

    public void setBedFile(final File bedFile) {
        this.bedFile = bedFile;
    }

    public boolean isBounds() {
        return bounds;
    }


    public void setBounds(final boolean bounds) {
        this.bounds = bounds;
    }

    public int getMinLength() {
        return minLength;
    }


    public void setMinLength(final int minLength) {
        this.minLength = minLength;
    }

    public int getMaxLength() {
        return maxLength;
    }


    public void setMaxLength(final int maxLength) {
        this.maxLength = maxLength;
    }

    @Override
    public String getName() {
        return "breakpoint";
    }

    @Override
    public String getDescription() {
        return "Plot the breakpoint distribution of a .bed file";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();

        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help(".bed input file").valueName("bed").required().get());
        parameters.add(JSAPParameters.switchParameter("bounds", 'b').type(File.class).help("Respect transcript bounds and ignore reads outside (negative breakpoints or breakpoints after transcript end)").get());
        parameters.add(JSAPParameters.flaggedParameter("min").help("Minimum transcript length").type(Integer.class).defaultValue("0").get());
        parameters.add(JSAPParameters.flaggedParameter("max").help("Maximum transcript length").type(Integer.class).defaultValue("100000000").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setBedFile(args.getFile("input"));
        setMaxLength(Integer.parseInt(args.getString("max")));
        setMinLength(Integer.parseInt(args.getString("min")));
        setBounds(args.userSpecified("bounds"));
        if (bedFile == null || !bedFile.exists()) {
            Log.error("Please specify a .bed file");
            return false;
        }
        return true;
    }

    @Override
    public Object call() throws Exception {
        BufferedReader reader = null;
        long[] sense = new long[bins];
        long[] antiSense = new long[bins];

        double b = 100.0 / bins;
        try {
            reader = new BufferedReader(new FileReader(getBedFile()));

            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }

                String[] e = line.split("\\s+");
                if (e.length < 4) {
                    throw new RuntimeException("Unable to find FRMD name for a read. Are you sure you passed a .bed file generated by the Flux Simulator ?");
                }
                String frmd = e[3];
                String[] elements = frmd.split("\\:");
                try {
                    int length = Integer.parseInt(elements[4]);
                    int breakPoint = Integer.parseInt(elements[7]);
                    char direction = Character.toUpperCase(elements[8].charAt(0));

                    if (length < minLength || length > maxLength) {
                        continue;
                    }
                    // negative starts
                    if (!bounds) {
                        if (breakPoint < 0) {
                            breakPoint = 0;
                        }
                        if (breakPoint > length) {
                            breakPoint = length;
                        }
                    } else {
                        if (breakPoint < 0 || breakPoint > length) {
                            continue;
                        }
                    }

                    // percent
                    double p = (breakPoint / (double) length) * 100.0;

                    double binD = Math.ceil(p / b) - 1.0;
                    if (binD < 0) {
                        binD = 0;
                    }
                    int bin = (int) binD;
                    if (direction == 'A') {
                        antiSense[bin]++;
                    } else if (direction == 'S') {
                        sense[bin]++;
                    } else {
                        throw new RuntimeException("Unknown direction " + direction);
                    }
                } catch (Exception error) {
                    throw new RuntimeException("Unable to parse FRMD name for a read.\n" +
                            " Are you sure you passed a .bed file generated by the Flux Simulator ?\n\n" +
                            "The line was: " + line + barna.commons.system.OSChecker.NEW_LINE +
                            "With FRMD: " + frmd, error);
                }
            }

        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        // print

        int x = (int) Math.floor(b);
        Log.println("Bin\tSense\tAntisense");
        for (int i = 0; i < bins; i++) {
            Log.println(((i + 1) * x) + "\t" + sense[i] + "\t" + antiSense[i]);
        }

        return null;
    }
}
