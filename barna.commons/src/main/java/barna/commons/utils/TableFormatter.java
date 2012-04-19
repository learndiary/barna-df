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

package barna.commons.utils;

import java.util.ArrayList;
import java.util.List;

/**
 * Format a table to simplify printing as string.
 *
 * @author Micha Sammeth (gmicha@googlemail.com)
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class TableFormatter {
    /**
     * Fill character
     */
    private static final char CHAR_FILL = ' ';
    /**
     * Space character
     */
    private static final char CHAR_SPACE = ' ';
    /**
     * Horizontal bar character
     */
    private static final char CHAR_HBAR = '=';

    /**
     * The table rows
     */
    private List<String[]> rows;
    /**
     * Number of columns
     */
    private final int col;
    /**
     * Store the maximum length of each column
     */
    private final int[] max;
    /**
     * Print first row as header
     */
    private final boolean header;

    /**
     * Create a new TableFormatter for {@code col} columns
     *
     * @param col the number of columns
     */
    public TableFormatter(int col) {
        this(col, false);
    }

    /**
     * Create a new TableFormatter for {@code col} columns
     *
     * @param col    the number of columns
     * @param header print first row as header
     */
    public TableFormatter(int col, boolean header) {
        if (col <= 0) {
            throw new IllegalArgumentException("#Columns must be > 0");
        }
        this.header = header;
        this.col = col;
        rows = new ArrayList<String[]>();
        max = new int[col];
        for (int i = 0; i < max.length; i++) {
            max[i] = 0;
        }
    }

    /**
     * Add a row to the table. The row mus provide an element for each column, otherwise an
     * exception will be thrown. You can pass null values.
     *
     * @param ss the row (mus consist of |columns| elements)
     */
    public void addRow(String... ss) {
        if (ss.length != col) {
            throw new RuntimeException("The passed row does not consist of " + col + " elements. You have to explicitly set a value for each column!");
        }
        for (int i = 0; i < ss.length; i++) {
            String s = ss[i];
            if (s != null) {
                max[i] = s.length() > max[i] ? s.length() : max[i];
            }
        }
        rows.add(ss);
    }

    @Override
    public String toString() {
        if (rows.size() == 0) {
            return "";
        }

        StringBuilder s = new StringBuilder();
        for (int i = 0; i < rows.size(); i++) {
            s.append("\t");
            for (int j = 0; j < rows.get(i).length; j++) {
                s.append(rows.get(i)[j]);
                for (int m = rows.get(i)[j].length(); m < max[j]; m++) {
                    s.append(CHAR_FILL);
                }
                s.append(CHAR_SPACE);
            }
            s.append("\n");
            if (i == 0 && header) {
                for (int j = 0; j < max.length; j++) {
                    for (int m = 0; m < max[j] + 1; m++)    // +1 for spacer
                    {
                        s.append(CHAR_HBAR);
                    }
                }
                s.append("\n");
            }
        }
        return s.toString();
    }

}
