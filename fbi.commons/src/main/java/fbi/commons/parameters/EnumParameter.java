package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class EnumParameter<E extends Enum<E>> extends Parameter<Enum<E>>{

    private E[] values;
    private E value;



    public EnumParameter(String name, String description, E defaultValue) {
        this(name, description, defaultValue, null);
    }

    public EnumParameter(String name, String description, E defaultValue, ParameterValidator validator) {
        this(name, description, defaultValue, defaultValue.getDeclaringClass(), validator);
    }
    public EnumParameter(String name, String description, E defaultValue, Class<E> values, ParameterValidator validator) {
        super(name, description, defaultValue, (Class<Enum<E>>) defaultValue.getClass(), validator);
        this.values = values.getEnumConstants();
    }


    Enum<E> get() {
        return value == null ? getDefault() : value;
    }

    void parse(String value) throws ParameterException{
        for (E e : values) {
            if(e.name().equalsIgnoreCase(value)){
                this.value = e;
                return;
            }
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value "+value);
    }
}
