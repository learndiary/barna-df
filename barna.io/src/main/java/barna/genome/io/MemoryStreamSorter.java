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

package barna.genome.io;

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.IOHandler;
import barna.commons.io.IOHandlerFactory;
import barna.commons.tools.ArrayUtils;
import barna.commons.tools.LineComparator;

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
