package fbi.commons.parameters;

import org.junit.Test;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class DoubleParameterTest {

    @Test
    public void testParse(){
        DoubleParameter test = new DoubleParameter("TEST");

        try {
            test.parse("10");
            assertEquals(10d, test.get(), 0.000001);
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("10.3");
            assertEquals(10.3, test.get(), 0.000001);
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("NaN");
            assertTrue(test.get().isNaN());
        } catch (ParameterException e) {
            e.printStackTrace();
            fail();
        }

        try {
            test.parse("inf");
            assertTrue(test.get().isInfinite());
        } catch (ParameterException e) {
            e.printStackTrace();
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
        DoubleParameter test = new DoubleParameter("TEST", "", 3d, 3d, 10d);


        assertEquals(3d, test.get(), 0.000001);
        try {
            test.parse("10");
            test.validate(null);
            assertEquals(10d,  test.get(), 0.000001);
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("9.3");
            test.validate(null);
            assertEquals(9.3, test.get(), 0.000001);
        } catch (ParameterException e) {
            fail();
        }

        try {
            test.parse("1");
            test.validate(null);
            fail();
        } catch (ParameterException e) {
            assertEquals("Parameter TEST value must be >= 3.0", e.getMessage());
        }

        try {
            test.parse("11");
            test.validate(null);
            fail();
        } catch (ParameterException e) {
            assertEquals("Parameter TEST value must be <= 10.0", e.getMessage());
        }


    }

}
