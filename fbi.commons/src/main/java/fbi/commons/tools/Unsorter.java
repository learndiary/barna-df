package fbi.commons.tools;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Read a file and randomize the order of lines
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Unsorter {

    public void unsort(File source, File target) throws IOException {
        if (target == null){
            target = new File(source.getParentFile(), source.getName()+"_unsorted");
        }

        IOHandler handler = IOHandlerFactory.createDefaultHandler();

        FileInputStream in = new FileInputStream(source);
        FileOutputStream out = new FileOutputStream(target);
        handler.addStream(in);
        handler.addStream(out);

        ByteArrayCharSequence ss = new ByteArrayCharSequence(1024);
        int line = handler.readLine(in, ss);
    }
}
