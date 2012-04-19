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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;

/**
 * InputStream on a random access file. You can use {@link #position(long)} to jump to
 * an absolute position. The {@link #skip(long)} method, as defined by InputStream, skips
 * relative to the current position
 *
 * @author Thasso Griebel
 */
public class RandomAccessInputStream extends InputStream {

    /**
     * The random access file
     */
    private RandomAccessFile file;

    /**
     * When mark is called, save here the current position so we can go back
     * on reset.
     */
    private long mark = -1;

    /**
     * True if we are to close the underlying random access file when this
     * stream is closed.
     */
    private boolean closeRandomAccess;

    /**
     * Create a new input stream on the given random access file. The
     * random access file will not be closed if the stream is closed and
     * the initial offeset is 0.
     *
     * @param file the file
     * @throws IOException in case of an error
     */
    public RandomAccessInputStream(RandomAccessFile file) throws IOException {
        this(file, false, 0);
    }

    /**
     * Creates a new read-only random access file and
     * the file is closed when this stream is closed
     *
     * @param file the source file
     * @throws IOException in case of an error
     */
    public RandomAccessInputStream(final File file) throws IOException {
        this(new RandomAccessFile(file, "r"), true, 0);
    }

    /**
     * Creates a new read-only random access file and
     * the file is closed when this stream is closed.
     *
     * @param file the source file
     * @param offset the start offset
     * @throws IOException in case of an error
     */
    public RandomAccessInputStream(final File file, final long offset)
            throws IOException {
        this(new RandomAccessFile(file, "r"), true, offset);
    }

    /**
     * Create a new input stream on the given random access file
     * and specify weather the file should be closed when this stream is closed.
     *
     * @param file the source file
     * @param closeRandomAccess close random access file if stream is closed
     * @param offset initial offset
     *
     * @throws IOException in case of an error
     */
    public RandomAccessInputStream(final RandomAccessFile file,
                                   final boolean closeRandomAccess, final long offset)
            throws IOException {
        super();
        this.closeRandomAccess = closeRandomAccess;
        this.file = file;
        if (offset > 0) {
            this.file.seek(offset);
        }
    }

    public int read() throws IOException {
        return this.file.read();
    }

    public int read(byte[] b, int off, int len) throws IOException {
        return this.file.read(b, off, len);
    }

    public int read(byte[] b) throws IOException {
        return this.file.read(b);
    }

    public long skip(long n) throws IOException {
        this.file.seek(this.file.getFilePointer() + n);
        return n;
    }

    /**
     * Returns the current position
     *
     * @return long current position
     * @throws IOException in case of an error
     */
    public long position() throws IOException {
        return this.file.getFilePointer();
    }

    /**
     * Set the current absolute position. In contrast to skip, this
     * takes the absolute position
     *
     * @param position absolute position
     * @throws IOException in case of an error
     */
    public void position(long position) throws IOException {
        this.file.seek(position);
    }

    public int available() throws IOException {
        long amount = this.file.length() - this.position();
        return (amount >= Integer.MAX_VALUE) ? Integer.MAX_VALUE
                : (int) amount;
    }

    public boolean markSupported() {
        return true;
    }

    public synchronized void mark(int readlimit) {
        try {
            this.mark = position();
        } catch (IOException e) {
            this.mark = -1;
        }
    }

    public synchronized void reset() throws IOException {
        if (this.mark == -1) {
            throw new IOException("Mark has not been set.");
        }
        position(this.mark);
    }

    public void close() throws IOException {
        try {
            super.close();
        } finally {
            if (this.closeRandomAccess) {
                this.file.close();
            }
        }
    }
}