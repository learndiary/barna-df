package fbi.commons.tools;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Sort the data from the input stream line by line using the specified field.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface StreamSorter {

    void sort(InputStream input, OutputStream output) throws IOException;

}
