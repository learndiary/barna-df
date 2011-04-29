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

    protected Parameter(String name, String description, T defaultValue, Class<T> type) {
        this.name = name;
        this.description = description;
        this.defaultValue = defaultValue;
        this.type = type;
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

    abstract T get();
    abstract void parse(String value) throws ParameterException;
}
