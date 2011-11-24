package fbi.commons.parameters;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class EnumParameterTest {


    private static enum Test{A, B}

    @org.junit.Test
    public void testParse(){
        EnumParameter<Test> p = new EnumParameter<Test>("TEST", "", Test.A);
        assertEquals(Test.A, p.get());

        try {
            p.parse("B");
            assertEquals(Test.B, p.get());
        } catch (ParameterException e) {
            e.printStackTrace();
            fail();
        }
    }
}
