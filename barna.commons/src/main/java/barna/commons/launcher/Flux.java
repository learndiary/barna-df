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

package barna.commons.launcher;

import barna.commons.Execute;
import barna.commons.log.Log;
import org.cyclopsgroup.caff.ref.AccessFailureException;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;
import org.cyclopsgroup.jcli.spi.ParsingContext;
import org.reflections.Reflections;
import org.reflections.scanners.SubTypesScanner;
import org.reflections.util.ClasspathHelper;
import org.reflections.util.ConfigurationBuilder;

import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.security.AccessControlException;
import java.security.Permission;
import java.util.*;


/**
 * Flux Simulator starter class. Contains the main method and parses command line arguments. During startup,
 * this checks for any {@link FluxTool} implementations and adds them to the list of available utils.
 */
@Cli(name = "flux", restrict = false)
public class Flux {
    /**
     * Required java version
     */
    private static final float FLUX_JAVA_VERSION = 1.6f;
    /**
     * The tool to be executed
     */
    private String toolName;

    /**
     * Show help message
     */
    private boolean help;

    /**
     * List available tools
     */
    private boolean listTools;

    /**
     * Number of executor threads
     */
    private int threads = 2;



    /**
     * Start the Flux simulator
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
    	
    	// set a system manager
    	SecurityManager sm= getSecurityManager();
    	if (sm!= null)
    		System.setSecurityManager(sm);
    	
        // load java.util.logger configuration
        Log.initialize();
        /*
        Check java version
         */
        checkJavaVersion();

        // register tools
        List<FluxTool> tools = findTools();

        // prepare the fluxInstance
        Flux fluxInstance = new Flux();
        ArgumentProcessor fluxArguments = ArgumentProcessor.newInstance(Flux.class);
        try {
            fluxArguments.process(args, fluxInstance);
        } catch (AccessFailureException ae) {
            Log.error("Error while processing arguments !");
            if (ae.getCause() instanceof InvocationTargetException) {
                Log.error(((InvocationTargetException) (ae.getCause())).getTargetException().getMessage());
            } else {
                Log.error(ae.getCause().getMessage());
            }
            System.exit(-1);
        } catch (Exception e) {
            Log.error("Error while processing arguments !");
            Log.error(e.getMessage());
            System.exit(-1);
        }


        // prepare tools
        List<org.cyclopsgroup.jcli.spi.Cli> toolClis = new ArrayList<org.cyclopsgroup.jcli.spi.Cli>();
        for (FluxTool fluxTool : tools) {
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
            ParsingContext context = toolArguments.createParsingContext();
            toolClis.add(context.cli());
        }

        if(fluxInstance.getToolName() == null){
            // check for a default tool
            String tool = System.getProperty("flux.tool");
            if(tool != null){
                fluxInstance.setToolName(tool);
            }
        }

        // find the tool to start
        FluxTool tool = null;
        if (fluxInstance.getToolName() != null) {
            int i = 0;
            for (org.cyclopsgroup.jcli.spi.Cli cli : toolClis) {
                if (cli.getName().equals(fluxInstance.getToolName())) {
                    tool = tools.get(i);
                    break;
                }
                i++;
            }
        }

        // delegate help prints to std err
        PrintWriter out = new PrintWriter(System.err);
        HelpPrinter printer = new HelpPrinter(out);


        if (fluxInstance.isHelp()) {
            // show help message
            if (tool == null) {
                printFluxHelp(tools, fluxArguments, out, printer);
            } else {
                out.println("General Options");
                printer.print(fluxArguments);
                out.println();
                out.println("Tool Options");
                out.println();
                ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
                printer.print(toolArguments);
                out.flush();
            }
            // exit after printing help
            System.exit(-1);
        }

        /**
         * Print tools
         */
        if(fluxInstance.isListTools()){
            printTools(tools, out);
            System.exit(-1);
        }


        if (tool != null) {
            // execute the tool
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(tool.getClass());
            toolArguments.process(args, tool);
            if (!tool.validateParameters(printer, toolArguments)) {
                out.flush();
                System.exit(-1);
            }

            try {
                // configure the executor
                Execute.initialize(fluxInstance.getThreads());

                tool.call();
            }catch (OutOfMemoryError outOfMemoryError){
                int mb = 1024*1024;
                long maxMemoryBytes = Runtime.getRuntime().maxMemory();
                long maxMemoryMB = (maxMemoryBytes/mb);
                Log.error("The Flux " + fluxInstance.getToolName() + " tool run into memory problems ! " +
                        "Please use the FLUX_MEM environment variable to increase the memory. For example: export FLUX_MEM=\"6G\"; flux -t capacitor ... to use" +
                        "6 GB of memory.");
                Log.error("Current memory setting : " + maxMemoryMB + " MB");
                System.exit(-1);
            } catch (IOException ioError) {
                // check for some specific errors
                if (ioError.getMessage().equals("No space left on device")) {
                    Log.error("[DISK] There is no space left on the device!");
                } else {
                    Log.error("Error while executing " + tool.getClass(), ioError);
                }
                System.exit(-1);
            } catch (Exception e) {
                Log.error(e.getMessage());
                Log.debug("\n\n");
                Log.debug("Error while executing " + tool.getClass() + " : " + e.getMessage(), e);
                System.exit(-1);
            } finally {
                // shutdown the executor
                Execute.shutdown();
            }
        } else {
            if (fluxInstance.getToolName() == null || fluxInstance.getToolName().isEmpty()) {
                Log.error("");
                Log.error("No tool specified!");
                Log.error("\n");
            } else {
                Log.error("");
                Log.error("Unable to find tool : " + fluxInstance.getToolName());
                Log.error("\n");
            }

            // warn and show help
            printFluxHelp(tools, fluxArguments, out, printer);
            System.exit(-1);
        }
    }

    /**
     * Sets a <code>SecurityManager</code> when required by property
     * variables.
     * @return a security manager or <code>null</code> if none is needed
     */
    protected static SecurityManager getSecurityManager() {
    	
    	final Properties props= System.getProperties();
    	Iterator iter= props.keySet().iterator();
    	boolean createSM= false;
    	for(;iter.hasNext()&& !createSM;) 
    		if (iter.next().toString().startsWith("flux.io.deny"))
    			createSM= true;
    	if (!createSM)
    		return null;
    	
		SecurityManager sm= new SecurityManager() {
			File systemTemp= new File(System.getProperty("java.io.tmpdir"));
			@Override
			public void checkPermission(Permission perm) {
				checkPermission(perm, null);
			};
			@Override
			public void checkPermission(Permission perm, Object context) {
				if (perm instanceof FilePermission) {
					if (props.containsKey("flux.io.deny.tmpdir")
								&& (perm.getActions().contains("write"))) {
							File f= new File(perm.getName());
							if (!(f.exists()|| f.isDirectory()))
								f= f.getParentFile();
							if (f.equals(systemTemp))
								throw new AccessControlException("access denied", perm);
					}
				}
			}
		};

		return sm;
	}

	private static void printFluxHelp(List<FluxTool> tools, ArgumentProcessor fluxArguments, PrintWriter out, HelpPrinter printer) {
        // show general help
        //fluxArguments.printHelp(out);
        printer.print(fluxArguments);
        printTools(tools, out);
        out.println();
        out.println("To get help for a specific tool try -t <tool> --help");
        out.println();
        out.flush();
    }

    /**
     * Print available flux tools
     *
     * @param tools the tools available flux tools
     * @param out the output stream
     */
    private static void printTools(List<FluxTool> tools, PrintWriter out) {
        out.println("\tTools available:");
        for (FluxTool fluxTool : tools) {
            ArgumentProcessor toolArguments = ArgumentProcessor.newInstance(fluxTool.getClass());
            ParsingContext context = toolArguments.createParsingContext();
            out.println("\t\t" + context.cli().getName() + " - " + context.cli().getDescription());
        }
        out.println();
    }

    /**
     * Search for flux tools. Currently the search is limited to classes in the barna package.
     *
     * @return tools all detected flux tools
     */
    public static List<FluxTool> findTools() {

        // scan the classpath to find tools
        List<FluxTool> tools = new ArrayList<FluxTool>();
        List<URL> fbiUrls = null;
        try {
        	fbiUrls= new ArrayList<URL>(ClasspathHelper.getUrlsForPackagePrefix("barna"));
        } catch (Throwable t) {
        	System.err.println(t);
        	System.currentTimeMillis();
        }

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
    @Option(name = "t", longName = "tool", description = "select the tool", displayName = "tool", required = true)
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
     * Set the log level. The Sting must be one of {@code NONE, INFO, ERROR, DEBUG}
     *
     * @param level the level
     */
    @Option(name = "", longName = "log", description = "Log level (NONE|INFO|ERROR|DEBUG)", defaultValue = "INFO", required = false)
    public void setLogLevel(String level) {
        Log.setLogLevel(level);
    }

    /**
     * Get the current log level
     *
     * @return level the log level
     */
    public String getLogLevel() {
        return Log.getLogLevel().toString();
    }

    /**
     * Disable interactivity
     *
     * @param detached detache
     */
    @Option(name = "f", longName = "force", description = "Disable interactivity. No questions will be asked", required = false)
    public void setDetached(boolean detached) {
        Log.setInteractive(!detached);
    }

    /**
     * Returns true if interactivity is disabled
     *
     * @return detached true if interactivity is disabled
     */
    public boolean isDetached() {
        return !Log.isInteractive();
    }

    /**
     * Returns the number of available background threads
     *
     * @return threads the number of background threads
     */
    public int getThreads() {
        return threads;
    }

    /**
     * Set the number of background threads
     *
     * @param threads number of background threads
     */
    @Option(name = "x", longName = "threads", description = "Number of threads, default is 2", required = false, defaultValue = "2")
    public void setThreads(final int threads) {
        this.threads = threads;
    }

    /**
     * Returs true if tools should be listed
     *
     * @return listTools list available tools
     */
    public boolean isListTools() {
        return listTools;
    }

    /**
     * Set listing tools
     *
     * @param listTools list tools
     */
    @Option(name = "", longName = "list-tools", description = "List available tools", required = false)
    public void setListTools(boolean listTools) {
        this.listTools = listTools;
    }

    /**
     * Read properties like version and build revision from jar file
     *
     * @return valid returns true if properties are valid and everything is fine
     */
    private static void checkJavaVersion() {
        try {
            float v = FLUX_JAVA_VERSION;
            String ver = System.getProperty("java.version");
            int p = ver.indexOf('.', 0);
            p = ver.indexOf('.', p + 1);
            float v2 = Float.parseFloat(ver.substring(0, p));
            if (v2 < v) {
                Log.error("Wrong java version, I need " + v + " but I found " + v2 + ".");

            }
        } catch (Exception e) {
            ; // :)
        }
    }

    /**
     * Cover build and version information
     */
    public static class FluxVersionInfo{
        private String libVersion;
        private String appVersion;
        private String buildDate;
        private String buildVersion;
        private String buildBranch;

        /**
         * Createa a new version info from a properties file. File path is relative to the classpath, Flux
         * is used to find the resource.
         *
         * @param fileName the file name (class path reference)
         */
        public FluxVersionInfo(String fileName){
            InputStream buildProperties = getClass().getResourceAsStream(fileName);
            if(buildProperties != null){
                Properties properties = new Properties();
                try {
                    properties.load(buildProperties);
                    libVersion = properties.getProperty("flux.version", "Unknown");
                    appVersion = properties.getProperty("flux.appversion", "Unknown");
                    buildDate = properties.getProperty("build.date", "Unknown");
                    buildBranch = properties.getProperty("build.branch", "Unknown");
                    buildVersion = properties.getProperty("build.version", "Unknown");
                    buildVersion = properties.getProperty("build.version", "Unknown");
                } catch (IOException ignore) {
                    // ignore
                }
            }
        }

        /**
         * Get the flux library version
         * @return libVersion the library version
         */
        public String getLibVersion() {
            return libVersion;
        }

        /**
         * Get the application version
         * @return appVersion the application version
         */
        public String getAppVersion() {
            return appVersion;
        }

        /**
         * Get the build date
         * @return buildDate the build date
         */
        public String getBuildDate() {
            return buildDate;
        }

        /**
         * Get the build version
         * @return buildVersion the build version
         */
        public String getBuildVersion() {
            return buildVersion;
        }

        /**
         * Get the branch name
         * @return buildBranch the branch name
         */
        public String getBuildBranch() {
            return buildBranch;
        }

        /**
         * Long multi-line string representation of all the version information
         * @return info the build info
         */
        public String toString(){
            StringBuilder builder = new StringBuilder();
            builder.append("Version ").append(appVersion).append('\n');
            builder.append("Flux Library ").append(appVersion).append('\n');
            builder.append("-----------------------------------------------\n");
            builder.append("Build Date ").append(buildDate).append('\n');
            builder.append("Build Version ").append(buildVersion).append('\n');
            builder.append("Build Branch ").append(buildBranch).append('\n');
            return builder.toString();
        }
        /**
         * Shor version info
         * @return info short version info
         */
        public String toShortString(){
            StringBuilder builder = new StringBuilder();
            builder.append("v").append(appVersion).append(" (Flux Library: ").append(libVersion).append(")\n");
            return builder.toString();
        }
    }
}
