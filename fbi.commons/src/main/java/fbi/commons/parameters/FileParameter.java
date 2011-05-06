package fbi.commons.parameters;

import java.io.File;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class FileParameter extends Parameter<File> {

    private File value;
    private FileNameParser nameParser;

    public FileParameter(String name) {
        this(name, "", null);
    }

    public FileParameter(String name, String description) {
        this(name, description, new File("."), null);
    }

    public FileParameter(String name, String description, File defaultValue) {
        this(name, description, defaultValue, null);
    }

    public FileParameter(String name, String description, File defaultValue, ParameterValidator validator) {
        this(name, description, defaultValue, validator, null);
    }

    public FileParameter(String name, String description, File defaultValue, ParameterValidator validator, FileNameParser nameParser) {
        super(name, description, defaultValue, File.class, validator);
        this.nameParser = nameParser;
    }


    @Override
    void set(File value) {
        this.value = value;
    }

    File get() {
        return value == null ? getDefault() : value;
    }

    void parse(String value) throws ParameterException {
        if (nameParser != null) {
            this.value = nameParser.parse(value);
        } else {
            this.value = new File(value);
        }
    }

    @Override
    public Parameter copy() {
        FileParameter fileParameter = new FileParameter(getName(), getDescription(), getDefault(), getValidator(), nameParser);
        fileParameter.set(get());
        return fileParameter;
    }
}
