package fbi.commons.parameters;

import java.math.BigDecimal;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
abstract class NumberParameter<T extends Number> extends Parameter<T> {
    T value;
    protected T minimumValue;
    protected T maximumValue;

    protected NumberParameter(String name, String description, T defaultValue, T minimumValue, T maximumValue, Class<T> type, ParameterValidator validator) {
        super(name, description, defaultValue, type, validator);
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
    void validate(ParameterSchema schema) throws ParameterException {
        if(getValidator() == null && minimumValue != null || maximumValue != null){
            if(minimumValue != null){
                if(new BigDecimal(value.toString()).compareTo(new BigDecimal(minimumValue.toString())) < 0){
                    throw new ParameterException(this, value.toString(), "Parameter " + this+ " value must be >= "+ minimumValue);
                }
            }

            if(maximumValue != null){
                if(new BigDecimal(value.toString()).compareTo(new BigDecimal(maximumValue.toString())) > 0){
                    throw new ParameterException(this, value.toString(), "Parameter " + this+ " value must be <= "+ maximumValue);
                }
            }
        }else{
            super.validate(schema);
        }
    }
}
