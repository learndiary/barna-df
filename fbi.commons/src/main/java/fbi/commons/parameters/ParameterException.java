package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ParameterException extends Exception {
    private Parameter parameter;

    public ParameterException(String message) {
        super(message);
    }

    public ParameterException(Parameter parameter, String value) {
        super("Invalid Parameter " + parameter.getName() + " with value "+ value );
        this.parameter = parameter;
    }

    public ParameterException(Parameter parameter, String value, String message) {
        super(message );
        this.parameter = parameter;
    }

}
