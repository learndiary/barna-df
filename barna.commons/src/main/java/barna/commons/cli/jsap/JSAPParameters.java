package barna.commons.cli.jsap;


import com.martiansoftware.jsap.*;

import java.math.BigDecimal;

public class JSAPParameters {
    /**
     * Create a builder for a flagged parameter
     *
     * @param longOption the long option
     * @param shortOption the short option
     * @return builder the builder
     */
    public static FlaggedParameterBuilder flaggedParameter(String longOption, char shortOption){
        return new FlaggedParameterBuilder(longOption, shortOption);
    }
    /**
     * Create a builder for a flagged parameter
     *
     * @param longOption the long option
     * @return builder the builder
     */
    public static FlaggedParameterBuilder flaggedParameter(String longOption){
        return new FlaggedParameterBuilder(longOption, (char)0);
    }

    /**
     * Create a builder for a switch parameter
     *
     * @param longOption the long option
     * @param shortOption the short option
     * @return builder the builder
     */
    public static SwitchParameterBuilder switchParameter(String longOption, char shortOption){
        return new SwitchParameterBuilder(longOption, shortOption);
    }
    /**
     * Create a builder for a switch parameter
     *
     * @param longOption the long option

     * @return builder the builder
     */
    public static  SwitchParameterBuilder switchParameter(String longOption){
        return new SwitchParameterBuilder(longOption, (char) 0);
    }
    /**
     * Create a builder for an unflagged parameter
     *
     * @param id the long option
     * @return builder the builder
     */
    public static UnflaggedParameterBuilder unflaggedParameter(String id){
        return new UnflaggedParameterBuilder(id);
    }


    public static class UnflaggedParameterBuilder extends ParameterBuilder<UnflaggedParameterBuilder>{
        protected boolean greedy;

        public UnflaggedParameterBuilder(String optionName) {
            super(optionName);
        }

        public UnflaggedParameterBuilder(String longOption, char shortOption) {
            super(longOption, shortOption);
        }

        public Parameter get(){
            return new UnflaggedOption(longOption, getParser(type), defaultValue, required, greedy, help);
        }

        /**
         * Make the parameter greedy
         *
         * @return builder the builder
         */
        public UnflaggedParameterBuilder greedy(){
            this.greedy = true;
            return this;
        }

    }

    public static class FlaggedParameterBuilder extends ParameterBuilder<FlaggedParameterBuilder>{
        String valueDescription;

        public FlaggedParameterBuilder(String optionName) {
            super(optionName);
        }

        public FlaggedParameterBuilder(String longOption, char shortOption) {
            super(longOption, shortOption);
        }

        public FlaggedParameterBuilder valueName(String valueName){
            this.valueDescription = valueName;
            return this;
        }

        public Parameter get(){
            FlaggedOption flaggedOption = new FlaggedOption(longOption, getParser(type), defaultValue, required, shortOption, longOption, help);
            if(valueDescription != null){
                flaggedOption.setUsageName(valueDescription);
            }
            return flaggedOption;
        }
    }

    public static class SwitchParameterBuilder extends ParameterBuilder<SwitchParameterBuilder>{
        public SwitchParameterBuilder(String optionName) {
            super(optionName);
        }

        public SwitchParameterBuilder(String longOption, char shortOption) {
            super(longOption, shortOption);
        }

        public Parameter get(){
            return new Switch(longOption, shortOption, longOption, help);
        }
    }

    public static abstract class ParameterBuilder<T extends ParameterBuilder>{
        protected String longOption;
        protected char shortOption;
        protected Class type = String.class;
        protected String help = "No help available";
        protected boolean required = false;
        protected String defaultValue;


        public ParameterBuilder(String optionName) {
            this.longOption = optionName;
        }

        public ParameterBuilder(String longOption, char shortOption) {
            this.longOption = longOption;
            this.shortOption = shortOption;
        }

        /**
         * Set the value type
         *
         * @param type value type
         * @return builder the builder
         */
        public T type(Class type){
            this.type = type;
            return (T) this;
        }
        /**
         * Set the help message
         *
         * @param help help message
         * @return builder the builder
         */
        public T help(String help){
            this.help = help;
            return (T) this;
        }
        /**
         * Set the default value
         *
         * @param defaultValue the default value
         * @return builder the builder
         */
        public T defaultValue(String defaultValue){
            this.defaultValue = defaultValue;
            return (T) this;
        }

        /**
         * Make the parameter required
         *
         * @return builder the builder
         */
        public T required(){
            this.required = true;
            return (T) this;
        }

        public abstract Parameter get();

        protected StringParser getParser(Class type) {
            if(type == String.class) return JSAP.STRING_PARSER;
            if(type == Boolean.class) return JSAP.BOOLEAN_PARSER;
            if(type == BigDecimal.class) return JSAP.BIGDECIMAL_PARSER;
            if(type == Integer.class) return JSAP.INTEGER_PARSER;
            if(type == Double.class) return JSAP.DOUBLE_PARSER;
            if(type == Character.class) return JSAP.CHARACTER_PARSER;
            if(type == Byte.class) return JSAP.BYTE_PARSER;
            if(type == Short.class) return JSAP.SHORT_PARSER;
            return JSAP.STRING_PARSER;
        }

    }

}
