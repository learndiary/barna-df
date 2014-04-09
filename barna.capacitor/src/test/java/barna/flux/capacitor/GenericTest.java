package barna.flux.capacitor;

import barna.commons.Execute;
import barna.commons.system.OSChecker;
import barna.io.FileHelper;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;

import java.io.File;
import java.util.Locale;
import junit.framework.Assert;
import org.junit.*;

/**
 * Generic test environment setup.
 * User: micha
 */
public class GenericTest {

    @BeforeClass
    public static void initExecuter() {
        //Force en-US locale to use "." as the decimal separator in Windows OS
        if (OSChecker.isWindows()) {
            Locale.setDefault(new Locale("en", "US"));
        }
        Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

    File currentTestDirectory = null;

    @Before
    public void setUpTest() throws Exception {
        currentTestDirectory = FileHelper.createTempDir(getClass().getSimpleName(), "", null);
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory);
        }
    }

    /**
     * Dummy test for IDE to keep imports.
     */
    @Test
    public void testDummy() {
        Assert.assertTrue(true);
    }

}
