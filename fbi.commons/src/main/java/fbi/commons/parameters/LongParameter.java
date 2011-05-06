package fbi.commons.parameters;

import java.math.BigDecimal;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class LongParameter extends NumberParameter<Long>{
    public LongParameter(String name) {
        this(name, "");
    }

    public LongParameter(String name, String description) {
        this(name, description, 0l, null, null);
    }

    public LongParameter(String name, String description, Long defaultValue) {
        this(name, description, defaultValue, null, null);
    }
    public LongParameter(String name, String description, Long defaultValue, Long minimumValue, Long maximumValue) {
        this(name, description, defaultValue, minimumValue, maximumValue, null);
    }
    public LongParameter(String name, String description, Long defaultValue, Long minimumValue, Long maximumValue, ParameterValidator validator) {
        super(name, description, defaultValue, minimumValue, maximumValue, Long.class, validator);
    }
    @Override
    void parse(String value) throws ParameterException {
        try{
            this.value = new BigDecimal(value).longValue();
            return;
        }catch (Exception e){
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value "+value);
    }

    @Override
    public Parameter copy() {
        LongParameter longParameter = new LongParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        longParameter.set(get());
        return longParameter;
    }


}
