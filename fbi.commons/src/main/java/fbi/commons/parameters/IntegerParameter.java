package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class IntegerParameter extends NumberParameter<Integer>{
    public IntegerParameter(String name) {
        this(name, "");
    }

    public IntegerParameter(String name, String description) {
        this(name, description, 0, null, null);
    }

    public IntegerParameter(String name, String description, Integer defaultValue) {
        this(name, description, defaultValue, null, null);
    }
    public IntegerParameter(String name, String description, Integer defaultValue, Integer minimumValue, Integer maximumValue) {
        super(name, description, defaultValue, minimumValue, maximumValue, Integer.class);
    }
}
