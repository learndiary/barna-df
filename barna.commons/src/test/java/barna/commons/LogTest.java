package barna.commons;

import org.junit.Test;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.logging.LogManager;
import java.util.logging.Logger;

import static org.junit.Assert.assertEquals;

/**
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class LogTest {

    @Test
    public void testJavaUtilLogConfig(){
        // load manually
        try{
            LogManager.getLogManager().readConfiguration(LogTest.class.getResourceAsStream("/logging.properties"));
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        final StringWriter out = new StringWriter();
        Log.logStream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                out.write(b);
            }
        });

        Log.setLogLevel(Log.Level.INFO);
        // issue a log message
        Logger.getLogger("MYLOGGER").info("first message");
        assertEquals("[MYLOGGER] first message\n", out.toString());
    }
}
