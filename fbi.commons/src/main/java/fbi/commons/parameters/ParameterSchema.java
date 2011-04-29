package fbi.commons.parameters;

import fbi.commons.Log;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

public abstract class ParameterSchema {
    private Map<String, Parameter> parameters;

    protected ParameterSchema() {
        this.parameters = new HashMap<String, Parameter>();

        /*
        do the reflection magic
         */

        Class<? extends ParameterSchema> clazz = getClass();
        Field[] fields = clazz.getFields();
        for (Field field : fields) {
            if(Parameter.class.isAssignableFrom(field.getType())){
                try {
                    if(!field.isAccessible()){
                        Log.warn("Parameter field " + field.getName() + " in " + clazz.getName() + " is not accessible!\n" +
                                "You can not access the parameter in a typesafe way outside of your parameter class.\nThis" +
                                " is probably not what you want! Use public static final for your parameter fields and make\n" +
                                "sure your Parameter class is accessible.");
                    }
                    field.setAccessible(true);
                    Parameter p = (Parameter) field.get(null);
                    if(p != null){
                        register(p);
                    }
                } catch (IllegalAccessException e) {
                    throw new RuntimeException("Unable to access parameter field " + field.getName() + " in " + clazz.getName());
                }
            }
        }
    }

    /**
     * Returns the number of parameters registered in this set
     *
     * @return size the number of registered parameters
     */
    public int size() {
        return parameters.size();
    }

    public void register(Parameter parameter){
        if (parameters.containsKey(parameter.getName())) throw new IllegalArgumentException("Paramter "+ parameter.getName() + " already exists !");
        this.parameters.put(parameter.getName(), parameter);
    }


    public <T> T get(Parameter<T> parameter){
        // find the parameter
        Parameter<T> local = parameters.get(parameter.getName());
        if(local == null){
            throw new IllegalArgumentException("Unknown parameter '" + parameter.getName()+"'");
        }
        return local.get();
    }
}
