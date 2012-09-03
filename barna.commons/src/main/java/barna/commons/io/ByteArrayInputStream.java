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
            if(seq.length() > 0 && seq.charAt(seq.end-1) == '\r'){
                seq.end -= 1;
                l--;
            }
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
