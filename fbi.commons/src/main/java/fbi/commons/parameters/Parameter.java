package fbi.commons.parameters;

/**
 * Represents a parameter parsed from string values
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public abstract class Parameter<T> {
    private String name;
    private String description;
    private T defaultValue;
    private Class<T> type;
    private ParameterValidator validator;

    protected Parameter(String name, String description, T defaultValue, Class<T> type, ParameterValidator validator) {
        this.name = name;
        this.description = description;
        this.defaultValue = defaultValue;
        this.type = type;
        this.validator = validator;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public Class<T> getType() {
        return type;
    }

    public T getDefault() {
        return defaultValue;
    }

    ParameterValidator getValidator() {
        return validator;
    }

    abstract T get();
    abstract void set(T value);
    abstract void parse(String value) throws ParameterException;
    void validate(ParameterSchema schema) throws ParameterException {
        if(validator != null){
            validator.validate(schema, this);
        }
    }

    @Override
    public String toString() {
        return name;
    }

    public String getValuesString(){
        return null;
    }
}
