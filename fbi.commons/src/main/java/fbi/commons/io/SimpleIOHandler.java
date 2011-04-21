package fbi.commons.io;

import fbi.commons.ByteArrayCharSequence;

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
class SimpleIOHandler implements IOHandler{
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
    private Map<OutputStream,OutputStream> outputStreams;

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
        if(stream == null) throw new NullPointerException();
        if (bufferSize <=0) throw new IllegalArgumentException("Stream buffer size must be > 0");
        if(stream instanceof InputStream){
            if(!inputStreams.contains(stream)){

                InputStream toAdd = (InputStream) stream;
                if(FORCE_BUFFERED_STREAMS && ! (stream instanceof BufferedInputStream)){
                    toAdd = new BufferedInputStream((InputStream) stream, IOHandler.DEFAULT_BUFFER_SIZE);
                }
                inputStreams.add((InputStream) stream);
                cachedStreams.put((InputStream) stream, new ByteArrayInputStream(new ByteArrayCharSequence(bufferSize), toAdd));
            }
        }else if(stream instanceof OutputStream){
            if(!outputStreams.containsKey(stream)){
                OutputStream toAdd = (OutputStream) stream;
                if(FORCE_BUFFERED_STREAMS && ! (stream instanceof BufferedOutputStream)){
                    toAdd = new BufferedOutputStream((OutputStream) stream);
                }
                outputStreams.put((OutputStream) stream, toAdd);
            }
        }else{
            throw new IllegalArgumentException("The given Object is neither an InputStream or an OutputStream");
        }
    }

    public void removeStream(Object stream) {
        if(stream == null) return;
        if(stream instanceof InputStream){
            inputStreams.remove(stream);
            cachedStreams.remove(stream);
        }else if(stream instanceof OutputStream){
            outputStreams.remove(stream);
        }
    }

    public boolean close() {
        // close all streams
        for (InputStream inputStream : inputStreams) {
            try {
                ByteArrayInputStream cc = cachedStreams.get(inputStream);
                if(cc != null)cc.close();
                else inputStream.close();
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
        if (outputStream == null){
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to write!");
        }
        outputStream.write(source, position, length);
    }


    public void writeLine(ByteArrayCharSequence cs, OutputStream out) throws IOException{
        OutputStream outputStream = outputStreams.get(out);
        if (outputStream == null){
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to write!");
        }

        outputStream.write(cs.a, 0, cs.length());
        outputStream.write(BYTE_NL);
    }

    public void writeLine(Object object, OutputStream out) throws IOException {
        if(object == null) throw new NullPointerException();
        if(bufferSequence == null){
            bufferSequence = new ByteArrayCharSequence(object.toString());
        }else{
            bufferSequence.reset();
            bufferSequence.append(object.toString());
        }
        writeLine(bufferSequence, out);
    }

    public ByteArrayCharSequence readLine(InputStream stream) throws IOException {
        ByteArrayInputStream cc = cachedStreams.get(stream);
        if(cc== null){
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to read!");
        }

        // reset the sequence
        cc.getSequence().reset();
        int read = cc.readLine();
        if(read < 0) return null;
        return cc.getSequence();
    }

    public int readLine(InputStream in, ByteArrayCharSequence cs) throws IOException {
        ByteArrayInputStream cc = cachedStreams.get(in);
        if(cc== null){
            throw new IllegalArgumentException("The stream was not registered with this handler or it was removed. Pleas ensure you call addStream before you try to read!");
        }
        ByteArrayCharSequence old = cc.getSequence();
        cs.reset();
        cc.setSequence(cs);
        int read = cc.readLine();
        cc.setSequence(old);
        return read;
    }
}
