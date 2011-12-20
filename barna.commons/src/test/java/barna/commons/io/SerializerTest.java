package barna.commons.io;

import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SerializerTest {

    @Test
    public void testArraySerialization(){
        int[] ii = {1,2,3};
        final StringWriter w = new StringWriter();
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                w.write(b);
            }
        };
        Serializer.save(ii, out);

        String string = "<int-array>\n" +
                "  <int>1\t2\t3</int>\n" +
                "</int-array>";
        assertEquals(string, w.toString());
    }

    @Test
    public void testArraySerializationRead(){
        int[] ii = {1,2,3};
        String string = "<int-array>\n" +
                "  <int>1\t2\t3</int>\n" +
                "</int-array>";
        final StringReader r = new StringReader(string);
        InputStream in  = new InputStream() {
            @Override
            public int read() throws IOException {
                return r.read();
            }
        };

        int[] loaded = (int[]) Serializer.load(in);
        for (int i = 0; i < loaded.length; i++) {
            assertEquals(ii[i], loaded[i]);
        }
    }

}
