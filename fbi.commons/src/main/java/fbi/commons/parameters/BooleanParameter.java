package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class BooleanParameter extends Parameter<Boolean>{
    private static final String[] TRUE = {"yes", "1", "true"};
    private static final String[] FALSE = {"no", "0", "false"};
    private boolean value;

    protected BooleanParameter(String name) {
        this(name, "");
    }

    protected BooleanParameter(String name, String description) {
        this(name, description, false);
    }

    protected BooleanParameter(String name, String description, boolean defaultValue) {
        this(name, description, defaultValue, null);
    }
    protected BooleanParameter(String name, String description, boolean defaultValue, ParameterValidator validator) {
        super(name, description, defaultValue, Boolean.class, validator);
        value = defaultValue;
    }

    @Override
    void set(Boolean value) {
        this.value = value;
    }

    @Override
    Boolean get() {
        return value;
    }

    @Override
    void parse(String value) throws ParameterException {
        String l = value.toLowerCase();
        for (String s : TRUE) {
            if(l.equals(s)){
                this.value = true;
                return;
            }
        }

        for (String s : FALSE) {
            if(l.equals(s)){
                this.value = false;
                return;
            }
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value '"+value + "'. Possible values are: [YES,NO]");
    }

    @Override
    public String getValuesString() {
        return "true|false or yes|no";
    }
}
