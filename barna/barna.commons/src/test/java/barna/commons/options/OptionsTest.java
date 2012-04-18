package barna.commons.options;

import barna.commons.launcher.Options;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * Test the Options class
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class OptionsTest {

    @Test
    public void testString(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addParameter("setValue", "Test String", 'c', "string");

        try {
            options.parse(new String[]{"--string", "test"});
            assertEquals("test", t.getValue());
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testStringShort(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addParameter("setValue", "Test String", 'c', "string");

        try {
            options.parse(new String[]{"-c", "test"});
            assertEquals("test", t.getValue());
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }
    @Test
    public void testStringLongShort(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addParameter("setValue", "Test String", 'c', "string");

        try {
            options.parse(new String[]{"-string", "test"});
            assertEquals("test", t.getValue());
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testBoolean(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addOption("setBooleanValue", "Test String", 'c', "string");

        try {
            options.parse(new String[]{"-string"});
            assertEquals(true, t.isBooleanValue());
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testBooleanParameterMethod(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addOption("sb", "Test String", 'c', "string");

        try {
            options.parse(new String[]{"-string"});
            assertEquals(true, t.isBooleanValue());
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }




    @Test
    public void testMultipleMappings(){
        OptionTarget t = new OptionTarget();
        Options options = new Options(t);
        options.addParameter("setValue", "Test String", 'c', "string");
        try{
            options.addParameter("setValue", "Test String", 'c', "string");
            fail();
        }catch (Exception e){

        }
    }


    public static class OptionTarget{
        String value;
        boolean booleanValue;

        public String getValue() {
            return value;
        }

        public void setValue(String value) {
            this.value = value;
        }

        public boolean isBooleanValue() {
            return booleanValue;
        }

        public void setBooleanValue() {
            this.booleanValue = true;
        }
        public void sb(boolean bb) {
            this.booleanValue = bb;
        }

    }
}
