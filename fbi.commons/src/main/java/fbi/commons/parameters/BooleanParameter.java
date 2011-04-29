package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class BooleanParameter extends Parameter<Boolean>{

    private boolean value;

    protected BooleanParameter(String name) {
        this(name, "");
    }

    protected BooleanParameter(String name, String description) {
        this(name, description, false);
    }

    protected BooleanParameter(String name, String description, boolean defaultValue) {
        super(name, description, defaultValue, Boolean.class);
    }

    @Override
    Boolean get() {
        return value;
    }

    @Override
    void parse(String value) throws ParameterException {
    }
}
