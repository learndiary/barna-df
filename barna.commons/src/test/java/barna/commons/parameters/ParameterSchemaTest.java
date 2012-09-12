package barna.commons.parameters;

import barna.commons.system.OSChecker;
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
    public void testUnknownParameter(){
        SimpleParameterSchema schema = new SimpleParameterSchema();

        try{
            schema.get(new StringParameter("Unknown"));
            fail();
        }catch (IllegalArgumentException e){
            assertEquals("Unknown parameter 'Unknown'", e.getMessage());
        }
    }

    @Test
    public void testParse(){
        SimpleParameterSchema schema = new SimpleParameterSchema();

        // parameters
        java.io.ByteArrayInputStream in = new java.io.ByteArrayInputStream(
                ("" +
                        "# Comment"+ OSChecker.NEW_LINE +
                        "STRING1 B"+ OSChecker.NEW_LINE+
                        "BOOLEAN\tyes"+ OSChecker.NEW_LINE+
                        "INT\t20"+ OSChecker.NEW_LINE+
                        "DOUBLE\t5.8"+ OSChecker.NEW_LINE+
                "").getBytes());

        try {
            schema.parse(in);

            assertEquals("B",schema.get(SimpleParameterSchema.STRING1));
            assertEquals(20, (int) schema.get(SimpleParameterSchema.INT1)  );
            assertEquals(5.8d, schema.get(SimpleParameterSchema.DOUBLE1), 0.000001 );
            assertTrue(schema.get(SimpleParameterSchema.BOOLEAN1));

        } catch (ParameterException e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testParseErrorValue(){
        SimpleParameterSchema schema = new SimpleParameterSchema();

        // parameters
        java.io.ByteArrayInputStream in = new java.io.ByteArrayInputStream(
                ("" +
                        "# Comment"+ OSChecker.NEW_LINE +
                        "STRING1 B"+ OSChecker.NEW_LINE+
                        "BOOLEAN\tyes"+ OSChecker.NEW_LINE+
                        "INT\t20"+ OSChecker.NEW_LINE+
                        "DOUBLE\t5.8sa"+ OSChecker.NEW_LINE+
                "").getBytes());

        try {
            schema.parse(in);
            fail();
        } catch (ParameterException e) {
            assertEquals("Error while parsing line 5. Unable to parse parameter DOUBLE with value 5.8sa", e.getMessage());
        }
    }
    @Test
    public void testDefaultAccess(){
        SimpleParameterSchema schema = new SimpleParameterSchema();
        assertEquals("A", schema.get(SimpleParameterSchema.STRING1));
        assertFalse(schema.get(SimpleParameterSchema.BOOLEAN1));
        assertEquals(0, (int)schema.get(SimpleParameterSchema.INT1));
        assertEquals(0d, schema.get(SimpleParameterSchema.DOUBLE1), 0.00001);
    }



    public static class SimpleParameterSchema extends ParameterSchema{
        public static final Parameter<String> STRING1 = Parameters.stringParameter("STRING1", "description", "A", "A", "B", "C");
        public static final Parameter<Boolean> BOOLEAN1 = Parameters.booleanParameter("BOOLEAN");
        public static final Parameter<Integer> INT1 = Parameters.intParameter("INT");
        public static final Parameter<Double> DOUBLE1 = Parameters.doubleParameter("DOUBLE");
    }

}
