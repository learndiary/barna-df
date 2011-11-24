package fbi.commons.parameters;

import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.*;

/**
 * Test string parameter
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class StringParameterTest {

    @Test
    public void testParse(){
        StringParameter test = new StringParameter("TEST");
        try{
            test.parse("VALUE");
            test.validate(null);
            assertEquals("VALUE", test.get());
        }catch (Exception e){fail();}
        try{
            test.parse("");
            test.validate(null);
            assertEquals("",test.get());
        }catch (Exception e){fail();}
        try{
            test.parse(null);
            test.validate(null);
            assertNull(test.get());
        }catch (Exception e){fail();}

        test = new StringParameter("TEST", "", "DEFAULT", Arrays.asList("DEFAULT", "A", "B", "C"));
        assertEquals("DEFAULT",test.get());

        try{
            test.parse("A");
            test.validate(null);
            assertEquals("A",test.get());
        }catch (Exception e){fail();}

        try{
            test.parse("B");
            test.validate(null);
            assertEquals("B",test.get());
        }catch (Exception e){fail();}

        try{
            test.parse("C");
            test.validate(null);
            assertEquals("C",test.get());
        }catch (Exception e){fail();}


        try{
            test.parse("ELSE");
            test.validate(null);
            fail();
        }catch (Exception e){
            assertEquals("Parameter TEST must be one of [DEFAULT, A, B, C]", e.getMessage());
        }
        try{
            test.parse(null);
            test.validate(null);
            fail();
        }catch (Exception e){
            assertEquals("Parameter TEST must be one of [DEFAULT, A, B, C]", e.getMessage());

        }
    }
}
