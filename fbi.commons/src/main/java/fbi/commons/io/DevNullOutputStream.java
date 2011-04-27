package fbi.commons.io;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Dummy output that moves all data to /dev/null
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 * `
 */
public class DevNullOutputStream extends OutputStream {

    @Override
    public void write(byte[] b) throws IOException {}

    @Override
    public void write(byte[] b, int off, int len) throws IOException {}

    @Override
    public void write(int b) throws IOException {}
}
