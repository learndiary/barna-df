/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.commons.io;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.converters.Converter;
import com.thoughtworks.xstream.converters.MarshallingContext;
import com.thoughtworks.xstream.converters.UnmarshallingContext;
import com.thoughtworks.xstream.io.HierarchicalStreamReader;
import com.thoughtworks.xstream.io.HierarchicalStreamWriter;

import java.io.InputStream;
import java.io.OutputStream;

/**
 * Simple XML serializer that can be used to serialize object to and from XML
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Serializer {

    private Serializer() {
        // hide the class
    }

    /**
     * Save the given object to the given stream
     *
     * @param object the object
     * @param out    the output stream
     */
    public static void save(Object object, OutputStream out) {
        XStream stream = createXStream();
        stream.toXML(object, out);
    }

    /**
     * Load an object from the given stream
     *
     * @param in the input stream
     * @return object the object loaded from the stream
     */
    public static Object load(InputStream in) {
        XStream stream = createXStream();
        return stream.fromXML(in);
    }

    private static XStream createXStream() {
        XStream stream = new XStream();
        stream.registerConverter(new ArrayConverter());
        return stream;
    }


    private static class ArrayConverter implements Converter {
        public void marshal(Object o, HierarchicalStreamWriter writer, MarshallingContext context) {
            Class type = o.getClass().getComponentType();
            writer.startNode(type.getName());
            StringBuffer b = new StringBuffer();

            if (type == int.class) {
                for (int o1 : (int[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == byte.class) {
                for (byte o1 : (byte[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == short.class) {
                for (short o1 : (short[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == long.class) {
                for (long o1 : (long[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == double.class) {
                for (double o1 : (double[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == float.class) {
                for (float o1 : (float[]) o) {
                    b.append(o1).append("\t");
                }
            } else if (type == char.class) {
                for (char o1 : (char[]) o) {
                    b.append(o1).append("\t");
                }
            }

            if (b.length() > 0) {
                b.deleteCharAt(b.length() - 1);
            }
            writer.setValue(b.toString());

            writer.endNode();
        }

        public Object unmarshal(HierarchicalStreamReader reader, UnmarshallingContext context) {
            reader.moveDown();
            String name = reader.getNodeName();
            String value = reader.getValue();
            String[] split = value.split("\t");
            if (name.equals("int")) {
                int[] ii = new int[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Integer.parseInt(split[i]);
                }
                return ii;
            } else if (name.equals("byte")) {
                byte[] ii = new byte[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Byte.parseByte(split[i]);
                }
                return ii;
            } else if (name.equals("short")) {
                short[] ii = new short[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Short.parseShort(split[i]);
                }
                return ii;
            } else if (name.equals("long")) {
                long[] ii = new long[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Long.parseLong(split[i]);
                }
                return ii;
            } else if (name.equals("double")) {
                double[] ii = new double[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Double.parseDouble(split[i]);
                }
                return ii;
            } else if (name.equals("float")) {
                float[] ii = new float[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = Float.parseFloat(split[i]);
                }
                return ii;
            } else if (name.equals("char")) {
                char[] ii = new char[split.length];
                for (int i = 0; i < split.length; i++) {
                    ii[i] = split[i].charAt(0);
                }
                return ii;
            }
            return null;
        }

        public boolean canConvert(Class aClass) {
            if (aClass.isArray()) {
                Class type = aClass.getComponentType();
                return type.isPrimitive();
            }
            return false;
        }
    }
}
