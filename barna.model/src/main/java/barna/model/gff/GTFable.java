package barna.model.gff;

import barna.commons.parameters.ParameterSchema;

/**
 * Interface for serializing GTF strings according to a <code>ParameterSchema</code>.
 *
 * User: micha
 * Date: 11/20/12
 * Time: 2:47 PM
 */
public interface GTFable {

    /**
     * Coordinates the composition of a GTF String according to the parameters in the schema.
     * @param schema parameter set with output options
     * @return a string literal representanting <code>this</code> object's attributes as requested
     * by the schema
     */
    public String toStringGTF(ParameterSchema schema);
}
