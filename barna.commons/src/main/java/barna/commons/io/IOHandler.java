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
import java.io.OutputStream;

/**
 * Manage disk access. Instances can manage {@link InputStream}s and {@link OutputStream}s
 * and then read or write data.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface IOHandler {
    /**
     * New line
     */
    static final byte BYTE_NL = '\n';
    /**
     * Windows new line
     */
    static byte BYTE_CR = '\r';
    /**
     * The default buffer size
     */
    static final int DEFAULT_BUFFER_SIZE = 10 * 1024;


    /**
     * Add an {@link InputStream} or {@link OutputStream} to the list of managed streams
     *
     * @param stream the stream
     */
    public void addStream(Object stream);

    /**
     * Add an {@link InputStream} or {@link OutputStream} to the list of managed streams using
     * the given buffer size to cache data
     *
     * @param stream     the stream
     * @param bufferSize the buffer size
     */
    public void addStream(Object stream, int bufferSize);

    /**
     * Remove a given {@link InputStream} or {@link OutputStream} from the list of managed streams
     *
     * @param stream
     */
    public void removeStream(Object stream);

    /**
     * Close the handler
     *
     * @return success returns true if the handler was successfully closed
     */
    public boolean close();

    /**
     * Write bytes to the given stream
     *
     * @param source   the bytes
     * @param position the start index in the byte[]
     * @param length   the length of the bytes
     * @param stream   the target stream
     * @throws IOException in case of a problem while writing to the stream
     */
    public void write(byte[] source, int position, int length, OutputStream stream) throws IOException;

    /**
     * Write the given character sequence as line
     *
     * @param cs  the line to write
     * @param out the line content
     * @throws IOException in case of an error
     */
    public void writeLine(ByteArrayCharSequence cs, OutputStream out) throws IOException;

    /**
     * Write the given character sequence as line
     *
     * @param object the line to write
     * @param out    the line content
     * @throws IOException in case of an error
     */
    public void writeLine(Object object, OutputStream out) throws IOException;


    /**
     * Read a line from the given stream and load it into the returned the sequence. The sequence
     * contains only the line. If the method returns null, nothing was read.
     *
     * @param stream the stream
     * @return charSequence returns the read line as a {@link ByteArrayCharSequence} or null if nothing could be read (EOF)
     * @throws IOException in case of any errors
     */
    public ByteArrayCharSequence readLine(InputStream stream) throws IOException;

    /**
     * Reads the line to the given sequence and returns the number of bytes read. If -1 is returned,
     * this stream end was reached and there is nothing to read.
     *
     * @param in the input stream
     * @param cs the sequence
     * @return index number of bytes read or -1 if nothing could be read (EOF)
     * @throws IOException
     */
    public int readLine(InputStream in, ByteArrayCharSequence cs) throws IOException;
}
