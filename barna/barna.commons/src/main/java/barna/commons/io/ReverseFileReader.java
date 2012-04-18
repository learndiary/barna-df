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

/*******************************************************************
 * Author:		Ryan D. Emerle
 * Date:			10.12.2004
 * Desc:			Reverse file reader.  Reads a file from the end to the
 *						beginning
 *
 * Known Issues:
 *						Does not support unicode!
 *******************************************************************/

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Read a file reverse, line by line. Note that this wil not work properly for non ASCII characters
 */
public class ReverseFileReader {
    /**
     * The file reader
     */
    private RandomAccessFile randomfile;
    /**
     * Current position
     */
    private long position;
    /**
     * The file
     */
    private File file;

    /**
     * Create the reader
     *
     * @param filename the filename
     */
    public ReverseFileReader(String filename) {
        this(new File(filename));
    }

    /**
     * The file
     *
     * @param file the file
     */
    public ReverseFileReader(File file) {
        if(file == null) throw new NullPointerException();
        this.file = file;
    }

    /**
     * Reads the next line
     *
     * @return line the next line
     * @throws Exception in case of any errors
     */
    public String readLine() throws IOException {
        if(randomfile==null){
            // initialize
            // Open up a random access file
            this.randomfile = new RandomAccessFile(file, "r");
            // Set our seek position to the end of the file
            this.position = this.randomfile.length();

            // Seek to the end of the file
            this.randomfile.seek(this.position);

            //Move our pointer to the first valid position at the end of the file.
            String thisLine = this.randomfile.readLine();
            while (thisLine == null || thisLine.equals("")) {
                this.position--;
                this.randomfile.seek(this.position);
                thisLine = this.randomfile.readLine();
                this.randomfile.seek(this.position);
            }
        }

        int thisCode;
        char thisChar;
        String finalLine = "";

        // If our position is less than zero already, we are at the beginning
        // with nothing to return.
        if (this.position < 0) {
            return null;
        }

        for (; ;) {
            // we've reached the beginning of the file
            if (this.position < 0) {
                break;
            }
            // Seek to the current position
            this.randomfile.seek(this.position);

            // Read the data at this position
            thisCode = this.randomfile.readByte();
            thisChar = (char) thisCode;


            // A little warning: the code posted there encodes bytes to chars by simply casting them.
            // This should work correctly if the encoding of the file is ISO-8859-1.
            // If it's something else, the non-ASCII chars will be mangled. As the author says, works fine for him.
            // If this is a line break or carrige return, stop looking
            if (thisCode == 13 || thisCode == 10) {
                // See if the previous character is also a line break character.
                // this accounts for crlf combinations
                this.randomfile.seek(this.position - 1);
                int nextCode = this.randomfile.readByte();
                if ((thisCode == 10 && nextCode == 13) || (thisCode == 13 && nextCode == 10)) {
                    // If we found another linebreak character, ignore it
                    this.position = this.position - 1;
                }
                // Move the pointer for the next readline
                this.position--;
                break;
            } else {
                // This is a valid character append to the string
                finalLine = thisChar + finalLine;
            }
            // Move to the next char
            this.position--;
        }
        // return the line
        return finalLine;
    }

    /**
     * Close the reader
     */
    public void close(){
        if(randomfile != null){
            try {
                randomfile.close();
            } catch (IOException ignored) {
                // ignore
            }
        }
    }
}