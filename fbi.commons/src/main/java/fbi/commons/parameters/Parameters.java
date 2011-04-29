package fbi.commons.parameters;

import java.util.Arrays;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Parameters {

    public static Parameter<String> stringParameter(String name){
        return stringParameter(name, "");
    }

    public static Parameter<String> stringParameter(String name, String description){
        return stringParameter(name, description, null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue){
        return stringParameter(name, description, defaultValue, null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, String...values){
        return new StringParameter(name, description, defaultValue, Arrays.asList(values));
    }

    public static Parameter<Boolean> booleanParameter(String name){
        return booleanParameter(name, "");
    }

    public static Parameter<Boolean> booleanParameter(String name, String description){
        return booleanParameter(name, description, false);
    }

    public static Parameter<Boolean> booleanParameter(String name, String description, boolean defautlValue){
        return new BooleanParameter(name, description, defautlValue);
    }

    public static Parameter<Double> doubleParameter(String name){
            return doubleParameter(name, "");
        }

    public static Parameter<Double> doubleParameter(String name, String description){
        return doubleParameter(name, description, 0d);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue){
        return new DoubleParameter(name, description, defaultValue, null, null);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue, double min, double max){
        return new DoubleParameter(name, description, defaultValue, min, max);
    }




    public static Parameter<Integer> intParameter(String name){
        return intParameter(name, "");
    }

    public static Parameter<Integer> intParameter(String name, String description){
        return intParameter(name, description, 0);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue){
        return new IntegerParameter(name, description, defaultValue, null, null);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue, int min, int max){
        return new IntegerParameter(name, description, defaultValue, min, max);
    }
}
