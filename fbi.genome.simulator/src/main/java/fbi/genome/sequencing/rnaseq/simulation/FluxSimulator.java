package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;
import org.cyclopsgroup.jcli.spi.ParsingContext;
import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.jar.Attributes;
import java.util.jar.JarFile;
import java.util.jar.Manifest;

//import gphase.solexa.lp.GraphLPsolver5;
//import gphase.solexa.simulation.Nebulizer;

// TODO check whether selected fragments concord w pro file

/**
 * Flux Simulator starter class
 *
 */
@Cli(name="flux", restrict = false)
public class FluxSimulator {
    /**
     * Current Flux Simulator version
     */
    public static String FLUX_VERSION = "";
    /**
     * Current Flux Simulator revision
     */
    public static String FLUX_REVISION = "";

    /**
     * The tool to be executed
     */
    private String toolName;

    /**
     * Show help message
     */
    private boolean help;


    /**
     * Start the Flux simulator
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /*
        Read properties
         */
        if(!readProperties()){
            System.exit(-1);
        }

        // register tools
        List<FluxTool> tools = findTools();

        // prepare the simulator
        FluxSimulator simulator = new FluxSimulator();
        ArgumentProcessor fluxArguments = ArgumentProcessor.newInstance(FluxSimulator.class);
        fluxArguments.process(args, simulator);

        // prepare tools
        List<org.cyclopsgroup.jcli.spi.Cli> toolClis = new ArrayList<org.cyclopsgroup.jcli.spi.Cli>();
        for (FluxTool fluxTool : tools) {
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
            ParsingContext context = toolArguments.createParsingContext();
            toolClis.add(context.cli());
        }



        // start
        Log.info("I am the Flux Simulator (v"+ FLUX_VERSION +" build"+FLUX_REVISION+"), nice to meet you!\n");

        // find the tool to start
        FluxTool tool = null;
        if(simulator.getToolName() != null){
            int i = 0;
            for (org.cyclopsgroup.jcli.spi.Cli cli : toolClis) {
                if(cli.getName().equals(simulator.getToolName())){
                    tool = tools.get(i);
                    break;
                }
                i++;
            }
        }

        if(simulator.isHelp()){
            // show help message
            PrintWriter stdoutWriter = new PrintWriter(System.err);
            if(tool == null){
                // show general help
                try {
                    fluxArguments.printHelp(stdoutWriter);
                    stdoutWriter.flush();

                    System.err.println("\tTools available:");
                    for (FluxTool fluxTool : tools) {
                        ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
                        ParsingContext context = toolArguments.createParsingContext();
                        System.err.println("\t\t" + context.cli().getName() + " - " + context.cli().getDescription());
                    }
                    System.err.println("\n");
                    System.err.println("To get help for a specific tool try -t <tool> --help");
                } catch (IOException e) {
                    Log.error("This is embarrassing. Something went terribly wrong and I can not even print the help message ! ", e);
                }
            }else{
                // show tool help
                // create the argument parser
                ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
                try {
                    toolArguments.printHelp(stdoutWriter);
                    stdoutWriter.flush();
                } catch (IOException e) {
                    Log.error("This is embarrassing. Something went terribly wrong and I can not even print the help message ! ", e);
                }
            }
            // exit after printing help
            System.exit(-1);
        }

        if(tool != null){
            // execute the tool
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
            toolArguments.process(args, tool);

            try {
                tool.call();
            } catch (Exception e) {
                Log.error("Error while executing "+ tool.getClass(), e);
            }
        }
    }

    /**
     * Search for flux tools. Currently the search is limited to classes in the fbi package.
     *
     * @return tools all detected flux tools
     */
    static List<FluxTool> findTools() {

        // scan the classpath to find tools
        List<FluxTool> tools = new ArrayList<FluxTool>();
        List<URL> fbiUrls = new ArrayList<URL>(ClasspathHelper.getUrlsForPackagePrefix("fbi"));

        ConfigurationBuilder config = new ConfigurationBuilder();
        config.useParallelExecutor();
        config.setScanners(new SubTypesScanner());
        config.setUrls(fbiUrls);

        Reflections reflections = new Reflections(config);


        Set<Class<? extends FluxTool>> toolClasses = reflections.getSubTypesOf(FluxTool.class);
        for (Class<? extends FluxTool> toolClass : toolClasses) {
            try {
                FluxTool fluxTool = toolClass.newInstance();
                tools.add(fluxTool);
            } catch (Exception e) {
                Log.error("Error while creating tool instance for " + toolClass.getName(), e);
                Log.error("Make sure the class exists and has a default constructor.");
            }
        }
        return tools;
    }

    /**
     * Get the name of the tool that is used
     *
     * @return toolName the name of the tool
     */
    public String getToolName() {
        return toolName;
    }

    /**
     * Set the name of the tool
     *
     * @param toolName toolname
     */
    @Option(name = "t", longName = "tool", description = "select the tool")
    public void setToolName(String toolName) {
        this.toolName = toolName;
    }

    /**
     * Should the help be displayed
     *
     * @return help true if help message should be displayed
     */
    public boolean isHelp() {
        return help;
    }

    /**
     * Show the help message
     *
     * @param help the help message
     */
    @Option(name = "h", longName = "help", description = "show help message")
    public void setHelp(boolean help) {
        this.help = help;
    }

    /**
     * Read properties like version and build revision from jar file
     *
     * @return valid returns true if properties are valid and everything is fine
     */
	private static boolean readProperties() {
        /*
        Find the manifest file and extract version revision adn jdk information
         */
        URL location = FluxSimulator.class.getResource("FluxSimulator.class");
        String fileString = location.toExternalForm();
        File jar = null;
        if(fileString.startsWith("jar")){
            fileString = fileString.substring(9);
            fileString = fileString.substring(0, fileString.lastIndexOf("!"));
            jar = new File(fileString);
        }
        String buildVersion = "";
        String buildRevision = "";
        String buildJDK = "";
        if(jar != null){
            try {
                JarFile jf = new JarFile(jar);
                Manifest manifest = jf.getManifest();
                Attributes mainAttributes = manifest.getMainAttributes();
                Object v = mainAttributes.getValue("Build-Version");
                Object r = mainAttributes.getValue("Build-Revision");
                Object j = mainAttributes.getValue("Flux-JDK");
                if(v != null){
                    buildVersion = v.toString();
                }
                if(r != null){
                    buildRevision = r.toString();
                }
                if(j != null){
                    buildJDK = j.toString();
                }
            } catch (IOException e) {}
        }


		if(!buildVersion.isEmpty()){
            FLUX_VERSION = buildVersion;
        }

        if(!buildRevision.isEmpty()){
            FLUX_REVISION = buildRevision;
        }

        if (!buildJDK.isEmpty()) {
			try {
				float v= Float.parseFloat(buildJDK);
				String ver= System.getProperty("java.version");
				int p= ver.indexOf('.', 0);
				p= ver.indexOf('.', p+1);
				float v2= Float.parseFloat(ver.substring(0, p));
				if (v2< v) {
					Log.error("Wrong java version, I need "+v+" but I found "+v2+".");
                    return false;
				}
			} catch (Exception e) {
				; // :)
			}
			
		}
        return true;
		
	}

	public static boolean invertTable(File invFile) {
		
		System.err.println("[INFO] inverting .pro file");
		File tmpFile= new File(invFile.getAbsolutePath()+"_inv");
		try {
			Vector<StringTokenizer> lineTokis= new Vector<StringTokenizer>();
			BufferedReader buffy= new BufferedReader(new FileReader(invFile));
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) 
				lineTokis.add(new StringTokenizer(s));
			buffy.close();
			
			int c= -1;
			for (int i = 0; i < lineTokis.size(); i++) {
				if (c< 0)
					c= lineTokis.elementAt(i).countTokens();
				else if (c!= lineTokis.elementAt(i).countTokens()) {
					System.err.println("\t[OHNO] inconsistent column count "+c+" <> "+lineTokis.elementAt(i).countTokens());
					return false;
				}
			}
			
			BufferedWriter writer= new BufferedWriter(new FileWriter(tmpFile));
			while (true) {
				for (int i = 0; i < lineTokis.size(); i++) {
					writer.write(lineTokis.elementAt(i).nextToken());
					if (i< lineTokis.size()-1)
						writer.write("\t");
					else
						writer.write("\n");
				}
				if (!lineTokis.elementAt(0).hasMoreTokens()) 
					break;					
			}			
			writer.flush();
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
			if (tmpFile.exists())
				tmpFile.delete();
			return false;
		}
		
		if (!invFile.delete()) 
			System.err.println("\t[OHNO] failed to remove "+invFile);
		if (!tmpFile.renameTo(invFile)) 
			System.err.println("\t[OHNO] failed to move "+tmpFile+" to "+invFile);
		return true;
	}
}
