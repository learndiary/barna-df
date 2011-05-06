package fbi.commons.parameters;

import java.math.BigDecimal;

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
        this(name, description, defaultValue, minimumValue, maximumValue, null);
    }
    public IntegerParameter(String name, String description, Integer defaultValue, Integer minimumValue, Integer maximumValue, ParameterValidator validator) {
        super(name, description, defaultValue, minimumValue, maximumValue, Integer.class, validator);
    }
    @Override
    void parse(String value) throws ParameterException {
        try{
            this.value = new BigDecimal(value).intValue();
            return;
        }catch (Exception e){
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value "+value);
    }

    @Override
    public Parameter copy() {
        IntegerParameter intparameter = new IntegerParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        intparameter.set(get());
        return intparameter;
    }

}
