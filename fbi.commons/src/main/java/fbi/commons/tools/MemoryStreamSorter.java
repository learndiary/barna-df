package fbi.commons.tools;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;

/**
 * Simplest sorter that first reads the complete file content before sorting. This is memory intensive!
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class MemoryStreamSorter implements StreamSorter{


    public void sort(InputStream input, OutputStream output, int field, boolean numeric, String fieldSeparator) throws IOException {
        ArrayList<String> lines = new ArrayList<String>();

        IOHandler io = IOHandlerFactory.createDefaultHandler();
        io.addStream(input);
        io.addStream(output);

        Comparator<String> comparator = new LineComparator(field, numeric, fieldSeparator);
        try{
            // read and sort
            ByteArrayCharSequence line;
            while( (line = io.readLine(input)) != null){
                ArrayUtils.addSorted(lines, line.toString(), comparator);
            }
            // write
            for (String l : lines) {
                io.writeLine(l, output);
            }
        }catch (IOException e){
            throw e;
        }catch(Exception error ){
            throw new RuntimeException("Error while sorting : " + error.getMessage(), error);
        }finally {
            io.close();
        }


    }

}
