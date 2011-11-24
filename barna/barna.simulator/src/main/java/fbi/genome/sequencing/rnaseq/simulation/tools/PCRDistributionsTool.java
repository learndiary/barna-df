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

package fbi.genome.sequencing.rnaseq.simulation.tools;

import com.thoughtworks.xstream.XStream;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.genome.sequencing.rnaseq.simulation.distributions.GCPCRDistribution;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Tool to create PCR distributions
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
@Cli(name="pcrdistributions", description = "Create a set of pcr distributions with n generations")
public class PCRDistributionsTool implements FluxTool<GCPCRDistribution>{
    /**
     * The number of bins
     */
    private int bins = 20;
    /**
     * The number of generations
     */
    private int generations = 15;
    /**
     * The output file
     */
    private File outputFile;

    /**
     * Validate
     */
    private String validate;

    /**
     * The probability to validate
     */
    private double validateProbability = 0.9;

    /**
     * Get the number of bins
     * @return bins the number of bins
     */
    public int getBins() {
        return bins;
    }

    /**
     * Set the number of bins
     *
     * @param bins the number of bins
     */
    @Option(name = "b", longName = "bins", description = "number of bins", required = false)
    public void setBins(final int bins) {
        if(bins <= 0 ) throw new IllegalArgumentException("Number of bins must be > 0");
        this.bins = bins;
    }

    /**
     * Get the number of generations
     * @return generations the number of generations
     */
    public int getGenerations() {
        return generations;
    }

    /**
     * Set the number of generations
     * @param generations the number of generations
     */
    @Option(name = "g", longName = "generations", description = "number of generations (PCR rounds)", required = true)
    public void setGenerations(final int generations) {
        if(generations <= 0 ) throw new IllegalArgumentException("Number of generations must be > 0");
        this.generations = generations;
    }

    /**
     * Get the output file
     *
     * @return output the
     */
    public File getOutputFile() {
        return outputFile;
    }

    /**
     * Set the target file
     *
     * param outputFile target the target file
     */
    @Option(name = "o", longName = "out", description = "output file", required = true)
    public void setOutputFile(final File outputFile) {
        if(outputFile == null) throw  new NullPointerException("You have to specify an output file!");
        this.outputFile = outputFile;
    }

    /**
     * Get the file name to validate
     *
     * @return name filename or 'default'
     */
    public String getValidate() {
        return validate;
    }

    /**
     * set the filename or 'default'
     *
     * @param validate filename
     */
    @Option(name = "v", longName = "validate", description = "validate file or 'default'", required = false)
    public void setValidate(final String validate) {
        this.validate = validate;
    }

    /**
     * The validation probability
     *
     * @return prob the validation probability
     */
    public double getValidateProbability() {
        return validateProbability;
    }

    /**
     * Set the validation probability, the PCR distribution for this probability is printed
     *
     * @param validateProbability the probability
     */
    @Option(name = "p", longName = "probability", description = "print distributions for this probability", required = false)
    public void setValidateProbability(final double validateProbability) {
        this.validateProbability = validateProbability;
    }

    @Override
    public boolean validateParameters(final HelpPrinter printer, final ArgumentProcessor toolArguments) {
        if(getValidate() == null){
            if (getOutputFile() == null) {
                printer.out.println("Please specify an output file");
                return false;
            }
        }
        return true;
    }

    @Override
    public GCPCRDistribution call() throws Exception {

        if(getValidate() != null){
            Log.info("PCR", "Validating " + getValidate());
            InputStream in = null;
            if(getValidate().equals("default")){
                in = getClass().getResourceAsStream("/pcr_15_20.dat");
            }else{
                in = new FileInputStream(new File(getValidate()));
            }

            GCPCRDistribution dist = load(in);
            Log.info("PCR", "Loaded distribution successfully");

            Log.println("PCR Distribution for " + StringUtils.fprint(getValidateProbability(), 3));

            dist.printBin(getValidateProbability());
            return null;
        }
        Log.info("PCR", "Create PCR Distributions for " + getBins() + " bins and " + getGenerations() + " generations");
        GCPCRDistribution dist = GCPCRDistribution.create(getBins(), getGenerations());
        Log.info("PCR", "Distributions created. Writing output.");

        XStream s = new XStream();
        OutputStream out = new GZIPOutputStream(new FileOutputStream(getOutputFile()));
        s.toXML(dist, out);
        out.close();

        Log.info("PCR", "Distributions written to "+ getOutputFile().getAbsolutePath());
        return dist;
    }

    /**
     * Load a distribution from file
     *
     * @param input the input stream
     * @return dist the distribution
     * @throws Exception in case of any errors
     */
    public static GCPCRDistribution load(InputStream input) throws  Exception{
        InputStream inputStream = new GZIPInputStream(input);
        XStream ss = new XStream();
        Object o = ss.fromXML(inputStream);
        return (GCPCRDistribution) o;
    }
}
