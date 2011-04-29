package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class NumberParameter<T extends Number> extends Parameter<T> {
    T value;
    protected T minimumValue;
    protected T maximumValue;

    protected NumberParameter(String name, String description, T defaultValue, T minimumValue, T maximumValue, Class<T> type) {
        super(name, description, defaultValue, type);
        this.minimumValue = minimumValue;
        this.maximumValue = maximumValue;
        if(defaultValue == null) throw new NullPointerException();
    }

    @Override
    T get() {
        if(value == null) return getDefault();
        return value;
    }

    @Override
    void parse(String value) throws ParameterException {

    }
}
