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

package barna.flux.simulator.tools;

import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.flux.simulator.distributions.GCPCRDistribution;
import barna.flux.simulator.distributions.PCRDistribution;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.thoughtworks.xstream.XStream;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Tool to create PCR distributions
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

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
    public void setValidateProbability(final double validateProbability) {
        this.validateProbability = validateProbability;
    }


    @Override
    public String getName() {
        return "pcrdistributions";
    }

    @Override
    public String getDescription() {
        return "Create a set of pcr distributions with n generations";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("out", 'o').type(File.class).help("Output File").required().valueName("file").get());
        parameters.add(JSAPParameters.flaggedParameter("generations", 'g').help("number of generations (PCR rounds)").required().valueName("rounds").get());
        parameters.add(JSAPParameters.flaggedParameter("bins", 'b').type(Integer.class).defaultValue("20").help("number of bins").get());
        parameters.add(JSAPParameters.flaggedParameter("validate").help("file to validate or 'default'").get());
        parameters.add(JSAPParameters.flaggedParameter("probability").type(Float.class).help("print distributions for this probability").get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setOutputFile(args.getFile("out"));
        setGenerations(args.getInt("generations"));
        setBins(args.getInt("bins"));
        if(args.userSpecified("validate")) setValidate(args.getString("validate"));
        if(args.userSpecified("probability")) setValidateProbability(args.getFloat("probability"));

        if(getValidate() == null){
            if (getOutputFile() == null) {
                Log.error("Please specify an output file");
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

        XStream s = createXStream();
        OutputStream out = new GZIPOutputStream(new FileOutputStream(getOutputFile()));
        s.toXML(dist, out);
        out.close();

        Log.info("PCR", "Distributions written to " + getOutputFile().getAbsolutePath());
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
        Object o = createXStream().fromXML(inputStream);
        return (GCPCRDistribution) o;
    }

    private static XStream createXStream() {
        XStream ss = new XStream();
        // make sure we use refactor save aliases
        ss.alias("GCPCRDistribution", GCPCRDistribution.class);
        ss.alias("PCRDistribution", PCRDistribution.class);
        return ss;
    }
}
