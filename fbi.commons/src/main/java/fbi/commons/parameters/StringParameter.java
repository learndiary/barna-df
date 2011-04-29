package fbi.commons.parameters;

import java.util.Collections;
import java.util.List;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class StringParameter extends Parameter<String>{

    private List<String> values;
    private String value;

    public StringParameter(String name) {
        this(name, "",  null);
    }

    public StringParameter(String name, String description) {
        this(name, description, "", null);
    }

    public StringParameter(String name, String description, String defaultValue) {
        this(name, description, defaultValue, null);
    }

    public StringParameter(String name, String description, String defaultValue, List<String> values) {
        super(name, description, defaultValue, String.class);
        this.values = values;
    }

    public List<String> getValues() {
        return Collections.unmodifiableList(values);
    }

    String get() {
        return value == null ? getDefault() : value;
    }

    void parse(String value) throws ParameterException{
        // if values are available, check that this is valid
        if(values != null && values.size() > 0){
            if(!values.contains(value)) throw new ParameterException(this, value);
        }
        this.value = value;
    }
}
