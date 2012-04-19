/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.commons.cli.jsap;


import com.martiansoftware.jsap.*;

import java.io.File;
import java.math.BigDecimal;

public class JSAPParameters {
    private static StringParser FILE_PARSER = new StringParser() {
        @Override
        public Object parse(String s) throws ParseException {
            return new File(s);
        }
    };

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
            if(type == File.class) {
                return FILE_PARSER;
            }
            return JSAP.STRING_PARSER;
        }

    }

}
