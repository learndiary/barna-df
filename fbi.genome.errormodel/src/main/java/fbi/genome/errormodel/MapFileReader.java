package fbi.genome.errormodel;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.commons.tools.Qualities;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * Read .map files
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class MapFileReader {
    /**
     * The input file
     */
    private File file;
    private Qualities.Technology qualityTechnology;
    /**
     * Current read
     */
    private Read read;
    /**
     * The input stream
     */
    private BufferedInputStream stream;
    /**
     * The io handler
     */
    private IOHandler io;

    public MapFileReader(File file, Qualities.Technology qualityTechnology) {
        this.file = file;
        this.qualityTechnology = qualityTechnology;
        this.read = new Read();
    }

    Read parseNext(boolean skip) throws IOException {
        if (stream == null){
            stream = new BufferedInputStream(new FileInputStream(file));
            io = IOHandlerFactory.createDefaultHandler();
            io.addStream(stream);
        }

        ByteArrayCharSequence line = io.readLine(stream);
        if(line == null) return null;
        if(skip) return read;


        // parse the line and
        String[] split = line.toString().split("\t");
        read.setName(split[0]);
        read.setSequence(split[1]);
        String quals = split[2];

        for (int i = 0; i < read.getLength(); i++) {
            char q = quals.charAt(i);
            read.getQualities()[i] = Qualities.quality(qualityTechnology, q);
        }

        return read;
    }


}
