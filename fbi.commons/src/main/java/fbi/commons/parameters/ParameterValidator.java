package fbi.commons.parameters;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface ParameterValidator {

    void validate(ParameterSchema schema, Parameter parameter) throws ParameterException;
}
