package fbi.commons.parameters;

import java.io.File;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class FileParameter extends Parameter<File>{

    private File value;

    public FileParameter(String name) {
        this(name, "",  null);
    }

    public FileParameter(String name, String description) {
        this(name, description, new File("."), null);
    }

    public FileParameter(String name, String description, File defaultValue) {
        this(name, description, defaultValue, null);
    }

    public FileParameter(String name, String description, File defaultValue, ParameterValidator validator) {
        super(name, description, defaultValue, File.class, validator);
    }


    File get() {
        return value == null ? getDefault() : value;
    }

    void parse(String value) throws ParameterException{
        this.value = new File(value);
    }
}
