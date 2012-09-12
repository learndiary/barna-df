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
import barna.commons.system.OSChecker;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A direct IO handler that does not perform any caching
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class SimpleIOHandler implements IOHandler {
    /**
     * Force that only buffered streams are used
     */
    private static final boolean FORCE_BUFFERED_STREAMS = true;
    /**
     * The list of managed input streams
     */
    private List<InputStream> inputStreams;
    /**
     * The list of managed output streams
     */
    private Map<OutputStream, OutputStream> outputStreams;

    /**
     * Cached input streams that are used to quickly fill the byte arrays char seqs
     */
    private Map<InputStream, ByteArrayInputStream> cachedStreams;
    /**
     * Default buffer size to use
     */
    private int bufferSize;

    /**
     * The buffer
     */
    private ByteArrayCharSequence bufferSequence;

    /**
     * Create a new handler
     */
    SimpleIOHandler() {
        this(IOHandler.DEFAULT_BUFFER_SIZE);
    }

    /**
     * Create a new handler using the given size for the initial buffers
     *
     * @param bufferSize the buffers size
     */
    SimpleIOHandler(int bufferSize) {
        this.bufferSize = bufferSize;
        this.inputStreams = new ArrayList<InputStream>();
        this.outputStreams = new HashMap<OutputStream, OutputStream>();
        this.cachedStreams = new HashMap<InputStream, ByteArrayInputStream>();
    }

    public void addStream(Object stream) {
        addStream(stream, bufferSize);
    }

    public void addStream(Object stream, int bufferSize) {
        if (stream == null) {
            throw new NullPointerException();
        }
        if (bufferSize <= 0) {
            throw new IllegalArgumentException("Stream buffer size must be > 0");
        }
        if (stream instanceof InputStream) {
            if (!inputStreams.contains(stream)) {

                InputStream toAdd = (InputStream) stream;
                if (FORCE_BUFFERED_STREAMS && !(stream instanceof BufferedInputStream)) {
                    toAdd = new BufferedInputStream((InputStream) stream, IOHandler.DEFAULT_BUFFER_SIZE);
                }
                inputStreams.add((InputStream) stream);
                cachedStreams.put((InputStream) stream, new ByteArrayInputStream(new ByteArrayCharSequence(bufferSize), toAdd));
            }
        } else if (stream instanceof OutputStream) {
            if (!outputStreams.containsKey(stream)) {
                OutputStream toAdd = (OutputStream) stream;
                if (FORCE_BUFFERED_STREAMS && !(stream instanceof BufferedOutputStream)) {
                    toAdd = new BufferedOutputStream((OutputStream) stream);
                }
                outputStreams.put((OutputStream) stream, toAdd);
            }
        } else {
            throw new IllegalArgumentException("The given Object is neither an InputStream or an OutputStream");
        }
    }

    public void removeStream(Object stream) {
        if (stream == null) {
            return;
        }
        if (stream instanceof InputStream) {
            inputStreams.remove(stream);
            cachedStreams.remove(stream);
        } else if (stream instanceof OutputStream) {
            outputStreams.remove(stream);
        }
    }

    public boolean close() {
        // close all streams
        for (InputStream inputStream : inputStreams) {
            try {
                ByteArrayInputStream cc = cachedStreams.get(inputStream);
                if (cc != null) {
                    cc.close();
                } else {
                    inputStream.close();
                }
            } catch (IOException e) {
                // ignore this one
            }
        }
        for (OutputStream outputStream : outputStreams.values()) {
            try {
                outputStream.flush();
                outputStream.close();
            } catch (IOException e) {
                // ignore this one
            }
        }
        return true;
    }

    public void write(byte[] source, int position, int length, OutputStream stream) throws IOException {
        // get the stream
        OutputStream outputStream = outputStreams.get(stream);
        if (outputStream == null) {
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to write!");
        }
        outputStream.write(source, position, length);
    }


    public void writeLine(ByteArrayCharSequence cs, OutputStream out) throws IOException {
        OutputStream outputStream = outputStreams.get(out);
        if (outputStream == null) {
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to write!");
        }

        outputStream.write(cs.chars, 0, cs.length());
        outputStream.write(OSChecker.NEW_LINE.getBytes());
    }

    public void writeLine(Object object, OutputStream out) throws IOException {
        if (object == null) {
            throw new NullPointerException();
        }
        if (bufferSequence == null) {
            bufferSequence = new ByteArrayCharSequence(object.toString());
        } else {
            bufferSequence.clear();
            bufferSequence.append(object.toString());
        }
        writeLine(bufferSequence, out);
    }

    public ByteArrayCharSequence readLine(InputStream stream) throws IOException {
        ByteArrayInputStream cc = cachedStreams.get(stream);
        if (cc == null) {
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to read!");
        }

        // reset the sequence
        cc.getSequence().clear();
        int read = cc.readLine();
        if (read < 0) {
            return null;
        }
        return cc.getSequence();
    }

    public int readLine(InputStream in, ByteArrayCharSequence cs) throws IOException {
        ByteArrayInputStream cc = cachedStreams.get(in);
        if (cc == null) {
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to read!");
        }
        ByteArrayCharSequence old = cc.getSequence();
        cs.clear();
        cc.setSequence(cs);
        int read = cc.readLine();
        cc.setSequence(old);
        return read;
    }

}
