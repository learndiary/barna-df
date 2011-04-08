package fbi.commons.file;

import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FileHelperTest {

    @Test
    public void testGuessFileSeparator(){
        File unixFile = new File(FileHelperTest.class.getResource("/TestFileSepUnix.txt").getFile());
        File nosepFile = new File(FileHelperTest.class.getResource("/TestFileSepNoSep.txt").getFile());
        File winFile = new File(FileHelperTest.class.getResource("/TestFileSepWindows.txt").getFile());

        String unix = FileHelper.guessFileSep(unixFile);
        String win = FileHelper.guessFileSep(winFile);
        String no = FileHelper.guessFileSep(nosepFile);

        assertEquals("", no);
        assertEquals("\r\n", win);
        assertEquals("\n", unix);
    }
}
