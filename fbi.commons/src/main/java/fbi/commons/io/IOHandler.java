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

package fbi.commons.io;

import fbi.commons.ByteArrayCharSequence;

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
