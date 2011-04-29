package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class DoubleParameter extends NumberParameter<Double>{
    public DoubleParameter(String name) {
        this(name, "");
    }

    public DoubleParameter(String name, String description) {
        this(name, description, 0d, null, null);
    }

    public DoubleParameter(String name, String description, Double defaultValue) {
        this(name, description, defaultValue, null, null);
    }
    public DoubleParameter(String name, String description, Double defaultValue, Double minimumValue, Double maximumValue) {
        super(name, description, defaultValue, minimumValue, maximumValue, Double.class);
    }
}
