package fbi.commons.parameters;

import java.io.File;
import java.util.Arrays;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Parameters {

    public static Parameter<String> stringParameter(String name){
        return stringParameter(name, "");
    }

    public static Parameter<String> stringParameter(String name, String description){
        return stringParameter(name, description, null, null, (String[])null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, ParameterValidator validator){
        return stringParameter(name, description, defaultValue, validator, (String[])null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, String...values){
        return stringParameter(name, description, defaultValue, null, values);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue,ParameterValidator validator, String...values){
        return new StringParameter(name, description, defaultValue, Arrays.asList(values), validator);
    }

    public static Parameter<Boolean> booleanParameter(String name){
        return booleanParameter(name, "");
    }

    public static Parameter<Boolean> booleanParameter(String name, String description){
        return booleanParameter(name, description, false, null);
    }

    public static Parameter<Boolean> booleanParameter(String name, String description, boolean defautlValue, ParameterValidator validator){
        return new BooleanParameter(name, description, defautlValue, validator);
    }

    public static Parameter<Double> doubleParameter(String name){
            return doubleParameter(name, "");
        }

    public static Parameter<Double> doubleParameter(String name, String description){
        return doubleParameter(name, description, 0d, null);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue, ParameterValidator validator){
        return new DoubleParameter(name, description, defaultValue, null, null, validator);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue, double min, double max, ParameterValidator validator){
        return new DoubleParameter(name, description, defaultValue, min, max, validator);
    }




    public static Parameter<Integer> intParameter(String name){
        return intParameter(name, "");
    }

    public static Parameter<Integer> intParameter(String name, String description){
        return intParameter(name, description, 0, null);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue,ParameterValidator validator){
        return new IntegerParameter(name, description, defaultValue, null, null, validator);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue, int min, int max, ParameterValidator validator){
        return new IntegerParameter(name, description, defaultValue, min, max, validator);
    }



    public static Parameter<File> fileParameter(String name){
        return fileParameter(name, "");
    }

    public static Parameter<File> fileParameter(String name, String description){
        return fileParameter(name, description, new File("."));
    }

    public static Parameter<File> fileParameter(String name, String description, File file){
        return fileParameter(name, description, file, null);
    }

    public static Parameter<File> fileParameter(String name, String description, File file, ParameterValidator validator){
        return new FileParameter(name, description, file, validator);
    }

}
