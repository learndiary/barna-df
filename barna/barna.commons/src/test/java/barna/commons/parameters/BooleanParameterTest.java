package barna.commons.parameters;

import org.junit.Test;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BooleanParameterTest {

    @Test
    public void testParse(){
        BooleanParameter test = new BooleanParameter("TEST");

        try{
            test.parse("YES");
            assertTrue(test.get());
            test.parse("YEs");
            assertTrue(test.get());
            test.parse("1");
            assertTrue(test.get());
            test.parse("true");
            assertTrue(test.get());

            test.parse("NO");
            assertFalse(test.get());
            test.parse("nO");
            assertFalse(test.get());
            test.parse("no");
            assertFalse(test.get());
            test.parse("0");
            assertFalse(test.get());
            test.parse("false");
            assertFalse(test.get());

        }catch (ParameterException e){
            fail();
        }

        try {
            test.parse("hmm");
            fail();
        } catch (ParameterException e) {
            assertEquals("Unable to parse parameter TEST with value 'hmm'. Possible values are: [YES,NO]", e.getMessage());
        }


    }
}
