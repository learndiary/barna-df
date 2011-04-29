package fbi.commons.parameters;

import org.junit.Test;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ParameterSchemaTest {

    @Test
    public void testParameterParsing(){
        SimpleParameterSchema schema = new SimpleParameterSchema();
        assertEquals(4, schema.size());
    }


    @Test
    public void testDefaultAccess(){
        SimpleParameterSchema schema = new SimpleParameterSchema();
        assertEquals("A", schema.get(SimpleParameterSchema.STRING1));
        assertFalse(schema.get(SimpleParameterSchema.BOOLEAN1));
        assertEquals(0, (int)schema.get(SimpleParameterSchema.INT1));
        assertEquals(0d, schema.get(SimpleParameterSchema.DOUBLE1), 0.00001);
    }

    @Test
    public void testUnknownParameter(){
        SimpleParameterSchema schema = new SimpleParameterSchema();

        try{
            schema.get(new StringParameter("Unknown"));
            fail();
        }catch (IllegalArgumentException e){
            assertEquals("Unknown parameter 'Unknown'", e.getMessage());
        }

    }



    public static class SimpleParameterSchema extends ParameterSchema{
        public static final Parameter<String> STRING1 = Parameters.stringParameter("STRING1", "description", "A", "A", "B", "C");
        public static final Parameter<Boolean> BOOLEAN1 = Parameters.booleanParameter("BOOLEAN");
        public static final Parameter<Integer> INT1 = Parameters.intParameter("INT");
        public static final Parameter<Double> DOUBLE1 = Parameters.doubleParameter("DOUBLE");
    }

}
