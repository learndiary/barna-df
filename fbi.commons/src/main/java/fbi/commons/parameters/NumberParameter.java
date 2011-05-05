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
    void set(T value) {
        this.value = value;
    }

    @Override
    T get() {
        if(value == null) return getDefault();
        return value;
    }

    @Override
    void validate(ParameterSchema schema) throws ParameterException {
        if(minimumValue != null || maximumValue != null){
            if(value == null){
                value = getDefault();
            }
            if(value != null && !value.equals(Double.NaN)){
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
            }
        }
        super.validate(schema);
    }

    @Override
    public String getValuesString() {
        if (minimumValue != null && maximumValue != null){
            return minimumValue +"<= number <= "+maximumValue;
        }else if(minimumValue != null){
            return minimumValue +"<= number";
        }else if(maximumValue != null){
            return "number <= "+maximumValue;
        }
        return "number";
    }
}
