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

package barna.io;

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.IOHandler;
import barna.commons.io.IOHandlerFactory;
import barna.commons.utils.ArrayUtils;
import barna.commons.utils.LineComparator;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;

/**
 * Simplest sorter that first reads the complete file content before sorting. This is memory intensive!
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class MemoryStreamSorter implements StreamSorter {

    private int field;
    private boolean numeric;
    private String fieldSeparator;

    public MemoryStreamSorter(int field, boolean numeric, String fieldSeparator) {
        this.field = field;
        this.numeric = numeric;
        this.fieldSeparator = fieldSeparator;
    }

    public void sort(InputStream input, OutputStream output) throws IOException {
        ArrayList<String> lines = new ArrayList<String>();

        IOHandler io = IOHandlerFactory.createDefaultHandler();
        io.addStream(input);
        io.addStream(output);

        Comparator<String> comparator = new LineComparator(numeric, fieldSeparator, field);
        try {
            // read and sort
            ByteArrayCharSequence line;
            while ((line = io.readLine(input)) != null) {
                ArrayUtils.addSorted(lines, line.toString(), comparator);
            }
            // write
            for (String l : lines) {
                io.writeLine(l, output);
            }
        } catch (IOException e) {
            throw e;
        } catch (Exception error) {
            throw new RuntimeException("Error while sorting : " + error.getMessage(), error);
        } finally {
            io.close();
        }
    }

}
