package barna.flux.capacitor.diffexp;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Iterator;

/**
 * Read from a GFF reader and provides iterator functions over the entries.
 * <p>
 * The reader is closed in case of an error and at the end of the
 * iteration. If you break the iteration manually, you have to close the reader.
 * </p>
 */
class GFFIterator implements Iterable<GFFEntry>, Iterator<GFFEntry>{
    /**
     * The input stream
     */
    private BufferedReader reader;
    /**
     * The next line to parse
     */
    private String line;
    /**
     *
     */
    private long lineCounter;

    /**
     * Create a new iterator
     *
     * @param reader the reader
     */
    public GFFIterator(Reader reader) {
        if(reader == null) throw new NullPointerException();
        this.reader = new BufferedReader(reader);
        this.lineCounter = 0;
    }

    @Override
    public Iterator<GFFEntry> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        if(reader == null) return false;
        try {
            while((line = reader.readLine()) != null){
                lineCounter+=1;
                line = line.trim();
                if(!line.isEmpty()){
                    return true;
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Error while reading from GFF stream!", e);
        }
        try {reader.close();} catch (IOException e1) {}
        return false;
    }

    @Override
    public GFFEntry next() {
        if(line != null){
            try {
                return GFFEntry.parse(line);
            } catch (Exception e) {
                try {reader.close();} catch (IOException e1) {}
                throw new RuntimeException("Error while parsing GFF line " + lineCounter + " : " + e.getMessage(),e);
            }
        }
        return null;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Close the underlying reader
     */
    public void close(){
        if(reader != null){
            try {reader.close();} catch (IOException e1) {}
        }
    }
}
