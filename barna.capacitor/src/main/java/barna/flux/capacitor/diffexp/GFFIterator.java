/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
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
