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

package barna.commons.io;

import barna.commons.ByteArrayCharSequence;

import java.io.IOException;
import java.io.InputStream;

/**
 * Reads data from a given input stream and fills a {@link ByteArrayCharSequence}.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ByteArrayInputStream extends InputStream {

    /**
     * The sequence to fill
     */
    private ByteArrayCharSequence seq;

    /**
     * The delegate stream
     */
    private InputStream inputStream;

    public ByteArrayInputStream(InputStream inputStream) {
        this(new ByteArrayCharSequence(128), inputStream);
    }

    public ByteArrayInputStream(ByteArrayCharSequence seq, InputStream inputStream) {
        if (seq == null || inputStream == null) {
            throw new NullPointerException();
        }
        this.seq = seq;
        this.inputStream = inputStream;
    }

    @Override
    public int read() throws IOException {
        int r = this.inputStream.read();
        if (r >= 0) {
            seq.append((byte) r);
        }
        return r;
    }

    /**
     * Reads bytes until a newline character is found. Returns the number of bytes read or -1
     * if nothing could be read
     *
     * @return length number of bytes read or -1
     * @throws IOException in case of an error
     */
    public int readLine() throws IOException {

        int c = read();
        if (c < 0) {
            return -1;
        }
        int l = 0;
        while (c >= 0 && c != '\n') {
            l++;
            c = read();
        }
        if (c == '\n') {
            // last character was a newline
            // remove this from the sequence
            seq.end = seq.end - 1;
        }
        return l;
    }

    @Override
    public void close() throws IOException {
        super.close();
        inputStream.close();
    }

    @Override
    public long skip(long n) throws IOException {
        return inputStream.skip(n);
    }

    @Override
    public int available() throws IOException {
        return inputStream.available();
    }

    @Override
    public boolean markSupported() {
        return inputStream.markSupported();
    }

    @Override
    public void mark(int readlimit) {
        super.mark(readlimit);
        inputStream.mark(readlimit);
    }

    @Override
    public void reset() throws IOException {
        inputStream.reset();
    }

    /**
     * Get the sequence
     *
     * @return sequence the sequence
     */
    public ByteArrayCharSequence getSequence() {
        return seq;
    }

    /**
     * Set the sequence that is filled by the read methods
     *
     * @param sequence the sequence
     */
    public void setSequence(ByteArrayCharSequence sequence) {
        if (sequence == null) {
            throw new NullPointerException("Null sequence is not permitted!");
        }
        this.seq = sequence;
    }
}
