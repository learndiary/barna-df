package barna.commons.parameters;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class IntegerParameterTest {

    @Test
    public void testParse(){
        IntegerParameter test = new IntegerParameter("TEST");

        try {
            test.parse("10");
            assertEquals(10, (int) test.get());
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("10.3");
            assertEquals(10, (int) test.get());
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("ERR");
            fail();
        } catch (ParameterException e) {
            assertEquals("Unable to parse parameter TEST with value ERR", e.getMessage());
        }
    }

    @Test
    public void testValidate(){
        IntegerParameter test = new IntegerParameter("TEST", "", 3, 3, 10);


        assertEquals(3, (int) test.get());
        try {
            test.parse("10");
            test.validate(null);
            assertEquals(10, (int) test.get());
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("10.3");
            test.validate(null);
            assertEquals(10, (int) test.get());
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("1");
            test.validate(null);
            fail();
        } catch (ParameterException e) {
            assertEquals("Parameter TEST value must be >= 3", e.getMessage());
        }

        try {
            test.parse("11");
            test.validate(null);
            fail();
        } catch (ParameterException e) {
            assertEquals("Parameter TEST value must be <= 10", e.getMessage());
        }


    }

}
